#!/usr/bin/python3 -u

import numpy as np
from utils_wgbs import MAX_PAT_LEN, validate_files_list, splitextgz, IllegalArgumentError, color_text, load_borders
from genomic_region import GenomicRegion
from view import ViewPat
import os.path as op
import sys


str2int = {'C': 2, 'T': 3, '.': 4, 'D': 5, 'X': 6}
int2str = {2: 'C', 3: 'T', 4: '.', 5: 'D', 1: ' ', 0: '', 6: 'X'}

num2color_dict = {
    'C': '01;31',  # red
    'T': '01;92',   # green
    'X': '01;33'
}


class PatVis:
    def __init__(self, args, file):
        self.gr = GenomicRegion(args)
        self.max_reps = args.max_reps if args.max_reps > 0 else sys.maxsize
        self.strict = args.strict
        self.strip = args.strip
        self.min_len = args.min_len
        self.start, self.end = self.gr.sites
        self.file = file
        self.no_color = args.no_color
        self.max_width = self.end - self.start + 2 * MAX_PAT_LEN  # maximal width of the output (in characters)
        self.blocks_path = args.blocks_path
        self.no_dense = args.no_dense
        self.uxm = args.uxm

        self.fullres = self.get_block()

    def insert_borders(self, markers):
        ctable = self.fullres['table']

        borders = load_borders(self.blocks_path, self.gr)# + self.start - self.fullres['start']
        if not borders.size:
            return self.fullres['text'], markers
        # pad right columns with space, if there are missing sites before the last border/s
        missing_width = borders[-1] - ctable.shape[1]
        if missing_width > 0:
            charar = np.chararray((ctable.shape[0], missing_width))
            charar[:] = ' '
            ctable = np.concatenate([ctable, charar], axis=1)

        # insert the borders:
        table = np.insert(ctable, borders, '|', axis=1)
        txt = '\n'.join(''.join(line) for line in table)
        rmark = ''
        j = 0
        for i in range(ctable.shape[1] + 1):
            if table[0][i] == '|':
                rmark += '|'
                continue
            elif j < len(markers):
                rmark += markers[j]
            else:
                rmark += ' '
            j += 1
        return txt, rmark

    def print_results(self):

        res = self.fullres
        if not res:
            return
        if res['score'] != 'NA':
            print('Methylation average: {}%'.format(res['score']))

        # Markers for sites of interest:
        markers = ' ' * (self.start - res['start']) + '+' * (self.end - self.start)
        txt = res['text']

        # print(self.fullres['table'])

        # Insert borders
        if self.blocks_path:
            txt, markers = self.insert_borders(markers)

        # Color text
        if not self.no_color:
            txt = color_text(txt, num2color_dict)

        if self.uxm:
            print("U, X, M: ({}, {}, {})".format(res['uxm'][0], res['uxm'][1], res['uxm'][2]))
        print(markers)
        print(txt)

    def get_block(self):
        df = ViewPat(self.file, None, self.gr, self.strict, min_len=self.min_len,
                strip = self.strip).perform_view(dump=False)
        #print(df)
        if not df.empty:
            return self.cyclic_print(df)

    def cyclic_print(self, df):
        table = np.zeros((df['count'].sum() + 1, self.max_width), dtype=np.int8)
        first_to_show = df.loc[0, 'start']
        row = -1
        u_count = 0
        m_count = 0
        x_count = 0

        for idx, read in df.iterrows():
            # _, read_start, patt, count = row
            read_start = int(read[1])
            patt = read[2]
            count = int(read[3])

            u_sites = 0
            m_sites = 0

            # perform multiple times for reads with count > 1, but no more than "max_reps" times:
            for c in range(min(self.max_reps, count)):

                # find the relative starting point of the current read
                col = read_start - first_to_show
                if col < 0:
                    raise IllegalArgumentError('Error: Input file must be sorted by CpG_ID!')

                # find the first available row to insert current read:
                if self.no_dense: # no_dense: present each read in a new line
                    row += 1
                else:
                    row = np.argmin(table[:, col])

                # make sure the slots are free
                assert(table[row, col:col + len(patt)].sum() == 0)
                assert(row < table.shape[0])

                # insert read and spaces:
                for j, l in enumerate(patt):
                    if self.uxm:
                        if l == 'C':
                            m_sites += 1
                        elif l == 'T':
                            u_sites += 1
                    else:
                        table[row, col + j] = str2int[l]
                table[row, :col][table[row, :col] == 0] = 1
                table[row, col + len(patt)] = 1
                if self.uxm:
                    total = u_sites + m_sites
                    u_prop = u_sites / total if total > 0 else 0
                    m_prop = m_sites / total if total > 0 else 0
                    is_u = u_prop >= self.uxm
                    is_m = m_prop >= self.uxm
                    if is_u:
                        u_count += 1
                    elif is_m:
                        m_count += 1
                    else:
                        x_count += 1
                    for j, l in enumerate(patt):
                        if is_u:
                            table[row, col + j] = str2int['T']
                        elif is_m:
                            table[row, col + j] = str2int['C']
                        else:
                            table[row, col + j] = str2int['X']


        nr_lines = int(np.argmin(table[:, 0]))
        width = np.max(np.argmin(table, axis=1))
        table = table[:nr_lines, :width]
        table[table == 0] = 1

        # if reqested range starts before the actual data (i.e first read),
        # concat empty columns at the beginning of the table
        if first_to_show > self.start:
            table = np.concatenate([np.ones((table.shape[0], first_to_show - self.start), dtype=np.uint8), table], axis=1)

        # Translate ints table to characters table
        table = table.astype(np.str)
        for key in int2str.keys():
            table = np.core.defchararray.replace(table, str(key), int2str[key])

        # Convert table to one long string
        res = '\n'.join([''.join(row) for row in table])

        fullres = {'start': first_to_show,
                   'chr': df.loc[0, 'chr'],
                   'text': res,
                   'width': width,
                   'table': table,
                   'score': calc_score(df),
                   'uxm': (u_count, x_count, m_count)}
        return fullres


def calc_score(df):
    # count C's and T's for the score:
    nm = (df['pat'].str.count('C') * df['count']).sum()
    ntotal = nm + (df['pat'].str.count('T') * df['count']).sum()
    score = int(100 * nm / ntotal) if ntotal else 'NA'
    return score

def main(args):
    validate_files_list(args.input_files, '.pat.gz')

    gr = GenomicRegion(args)
    print(gr)
    for pat_file in args.input_files:
        print(splitextgz(op.basename(pat_file))[0])     # print file name
        PatVis(args, pat_file).print_results()
