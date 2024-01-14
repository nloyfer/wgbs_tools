#!/usr/bin/python3 -u

import numpy as np
import os.path as op
import sys
import pandas as pd
import argparse
import re
from operator import itemgetter
from utils_wgbs import validate_file_list, splitextgz, IllegalArgumentError, \
                       color_text, load_borders, drop_dup_keep_order, read_shell
from genomic_region import GenomicRegion
from cview import view_gr


FULL_CIRCLE = '\u25CF'
DASH = '\u2014'
BORDER = '|'
str2int = {c: i for i, c in enumerate([''] + list(' .CTUXMctga'))}
int2str = {v: k for k, v in str2int.items()}

num2color_dict = {
    'C': '01;31',   # red
    'T': '01;32',   # green
    'X': '01;33',   # yellow
    'M': '01;31',   # red
    'U': '01;32',   # green
    'c': '01;106',  # blue
    't': '01;90',   # black
    'g': '01;91',   # red
    'a': '01;92',   # green
}

num2color_dict2 = {
    'T': '01;34',   # blue
    'C': '01;33',   # yellow
    'X': '01;33',   # yellow
    'M': '01;31',   # red
    'U': '01;32',   # green
}

def table2text(table):
    """ join table to text.
        add strikethrough to borders if they cross a read """
    lines = []
    for line in table:
        nline = ''
        for i, ch in enumerate(line):
            nline += ch
            # add strikethrough to the borders
            if (len(line) - 1 > i > 0) and (ch == BORDER) \
                and (set((line[i - 1], line[i + 1])) <= set('CT.')):
                nline += '\u0336'
        lines.append(nline)
    return '\n'.join(lines)


class PatVis:
    def __init__(self, args, pat_path):
        self.gr = GenomicRegion(args)
        self.args = args
        self.max_reps = args.max_reps if args.max_reps > 0 else sys.maxsize
        self.start, self.end = self.gr.sites
        self.pat_path = pat_path
        self.blocks_path = args.blocks_path
        self.uxm = args.uxm
        self.uxm_counts = {'U': 0, 'X': 0, 'M': 0}
        self.fullres = self.get_block()

    def insert_borders(self, markers):
        ctable = self.fullres['table']
        start = self.fullres['start']

        # load borders from file
        # build gr to span the whole table
        bsites = '{}-{}'.format(start, start + ctable.shape[1])
        table_gr = GenomicRegion(sites=bsites, genome_name=self.gr.genome_name)
        borders = load_borders(self.blocks_path, table_gr, self.args.genome)
        if not borders.size:
            return self.fullres['text'], markers
        # pad right columns with space, if there are missing sites before the last border/s
        missing_width = borders[-1] - ctable.shape[1]
        if missing_width > 0:
            charar = np.chararray((ctable.shape[0], missing_width))
            charar[:] = ' '
            ctable = np.concatenate([ctable, charar], axis=1)

        # insert the borders:
        txt = table2text(np.insert(ctable, borders, BORDER, axis=1))

        # insert the borders to the markers line:
        markers_arr = np.array(list(markers.ljust( ctable.shape[1])))[:, None]
        rmark = ''.join(np.insert(markers_arr, borders, BORDER))
        return txt, rmark

    def print_results(self):

        res = self.fullres
        if not res:
            return
        if res['score'] != 'NA':
            to_print = f'Methylation average: {res["score"]}%'
            if self.uxm:
                u_thresh = round(self.uxm * 100, 2)
                m_thresh = round((1 - self.uxm) * 100, 2)
                to_print += f'\nUXM {u_thresh:g}/{m_thresh:g} '
                arr = np.array(itemgetter(*list('UXM'))(res['uxm']))
                to_print += '[{}/{}/{}]'.format(*arr)
                to_print += ' [{:.1%}/{:.1%}/{:.1%}]'.format(*(arr/arr.sum()))
            print(to_print)

        # Markers for sites of interest:
        markers = ' ' * (self.start - res['start']) + '+' * (self.end - self.start)
        txt = res['text']

        # print(self.fullres['table'])

        # Insert borders
        if self.blocks_path:
            txt, markers = self.insert_borders(markers)

        # Color text
        if not self.args.no_color:
            nd = num2color_dict2 if self.args.yebl else num2color_dict
            txt = color_text(txt, nd)
        if not self.args.text:
            txt = re.sub('[CTUXM]', FULL_CIRCLE, txt)               # letters -> circles
            txt = re.sub('\.', DASH, txt)                           # dots -> dashes
            for x, X in zip('acgt', 'ACGT'):
                txt = re.sub(x, X, txt)
            if self.args.strike:
                txt = txt.replace(FULL_CIRCLE, FULL_CIRCLE + '\u0336')  # strikethrough
                # txt = txt.replace(FULL_CIRCLE, '\u0336' + FULL_CIRCLE)  # strikethrough
        print(markers)
        print(txt, flush=True)  # the flush is needed to handle BrokenPipeError

    def get_block(self):
        cmd = view_gr(self.pat_path, self.args, get_cmd=True)
        df = read_shell(cmd, index_col=None)
        if not df.empty:
            names = ['chr', 'start', 'pat', 'count']
            df.columns = names + list(df.columns)[len(names):]
            return self.cyclic_print(df)

    def read_uxm(self, patt, count):
        u_sites = patt.count('T')
        m_sites = patt.count('C')
        total = u_sites + m_sites
        assert(total > 0)
        if u_sites / total >= self.uxm:
            uxm_status = 'U'
        elif m_sites / total >= self.uxm:
            uxm_status = 'M'
        else:
            uxm_status = 'X'
        self.uxm_counts[uxm_status] += count
        return uxm_status * len(patt)

    def insert_read_to_table(self, read, table, shift):
        read_start = int(read.iloc[1])
        patt = read.iloc[2]
        count = int(read.iloc[3])

        # skip empty (all dots) reads:
        if not patt.strip('.'):
            return

        if self.uxm:
            patt = self.read_uxm(patt, count)
        patt_ints = [str2int[l] for l in patt]

        # perform multiple times for reads with count > 1, but no more than "max_reps" times:
        for c in range(min(self.max_reps, count)):
            # find the relative starting point of the current read
            col = read_start - shift
            if col < 0:
                raise IllegalArgumentError('Error: Pat is not sorted!')

            # find the first available row to insert current read:
            if self.args.no_dense: # no_dense: present each read in a new line
                row = np.argmin(table.sum(axis=1))
            else:
                row = np.argmin(table[:, col])

            # make sure the slots are free
            assert(table[row, col:col + len(patt)].sum() == 0)
            assert(row < table.shape[0])

            # insert read and spaces:
            table[row, col:col + len(patt)] = patt_ints
            table[row, :col][table[row, :col] == 0] = 1  # before read
            table[row, col + len(patt)] = 1              # after read

    def cyclic_print(self, df):
        longest_read = df.pat.str.len().max()
        max_width = self.end - self.start + 2 * longest_read  # maximal width of the output (in characters)
        table = np.zeros((df['count'].sum() + 1, max_width), dtype=np.int8)
        first_to_show = df.loc[0, 'start']

        for _, read in df.iterrows():
            self.insert_read_to_table(read[:4], table, first_to_show)

        nr_lines = int(np.argmin(table[:, 0]))
        width = np.max(np.argmin(table, axis=1))
        table = table[:nr_lines, :width]
        table[table == 0] = 1

        # if reqested range starts before the actual data (i.e first read),
        # concat empty columns at the beginning of the table
        if first_to_show > self.start:
            table = np.concatenate([np.ones((table.shape[0], first_to_show - self.start), dtype=np.uint8), table], axis=1)

        int_table = table.copy()
        # Translate ints table to characters table
        table = pd.DataFrame(data=table).replace(int2str).values

        # Convert table to one long string
        res = '\n'.join([''.join(row) for row in table])

        fullres = {'start': first_to_show,
                   'chr': df.loc[0, 'chr'],
                   'text': res,
                   'table': table,
                   'int_table': int_table,
                   'score': calc_score(df),
                   'uxm': self.uxm_counts}
        return fullres


def calc_score(df):
    # count C's and T's for the score:
    nm = (df['pat'].str.count('C') * df['count']).sum()
    ntotal = nm + (df['pat'].str.count('T') * df['count']).sum()
    score = int(100 * nm / ntotal) if ntotal else 'NA'
    return score


def main(args):
    validate_file_list(args.input_files, '.pat.gz')

    # drop duplicated files, while keeping original order
    input_files = drop_dup_keep_order(args.input_files)

    gr = GenomicRegion(args)
    print(gr)
    for pat_file in input_files:
        print(splitextgz(op.basename(pat_file))[0])     # print file name
        PatVis(args, pat_file).print_results()

