#!/usr/bin/python3 -u

import argparse
import os
import numpy as np
import os.path as op
import pandas as pd
from utils_wgbs import validate_files_list
from multiprocessing import Pool
from os.path import splitext, basename
import sys
from utils_wgbs import load_beta_data, trim_to_uint8, default_blocks_path, eprint, GenomeRefPaths, splitextgz, \
                        IllegalArgumentError


def get_bins(df, genome):
    end = GenomeRefPaths(genome).count_nr_sites() + 1   # 28217449 for hg19
    arr = np.unique(np.concatenate([[1], df['startCpG'], df['endCpG'], [end]]))
    arr.sort()
    isin = np.isin(arr, np.concatenate([df['startCpG'], [df['endCpG'][df.shape[0] - 1]]]))
    return arr - 1, isin


def apply_filter_wrapper(args, blocks_bins, finds, beta_path, df):
    try:
        # load beta file:
        data = load_beta_data(beta_path)

        # reduce to blocks:
        blocks_bins[-1] -= 1
        reduced_data = np.add.reduceat(data, blocks_bins)[finds][:-1]

        # dump to file
        suff = '.lbeta' if args.lbeta else '.bin'
        out_name = splitext(splitext(basename(args.blocks_file))[0])[0]
        out_name = splitext(basename(beta_path))[0] + suff
        out_name = op.join(args.out_dir, out_name)

        trim_to_uint8(reduced_data, args.lbeta).tofile(out_name)
        eprint(out_name)

        if args.bedGraph:
            with np.errstate(divide='ignore', invalid='ignore'):
                beta_vals = reduced_data[:, 0] / reduced_data[:, 1]
                eprint(beta_vals.shape, df.shape)
            # beta_vals[reduced_data[:, 1] == 0] = np.nan
            df['beta'] = beta_vals
            df.to_csv(out_name.replace(suff, '.bedGraph'), sep='\t',
                      index=None, header=None, na_rep=-1,
                      float_format='%.2f')

    except Exception as e:
        print('Failed with beta', beta_path)
        print('Exception:', e)


def blocks_cleanup(df, args):
    # remove duplicated blocks
    nr_orig = df.shape[0]
    df.drop_duplicates(inplace=True)
    nr_non_dup = df.shape[0]
    if nr_orig != nr_non_dup:
        eprint('removed {:,} duplicated blocks'.format(nr_orig - nr_non_dup))
    
    # remove blocks with no CpGs
    nr_removed = df[df.startCpG == df.endCpG].shape[0]
    if nr_removed:
        eprint('removed {:,} regions with no CpGs'.format(nr_removed))

    if args.debug:
        eprint(df[df.startCpG == df.endCpG])

    df = df[df.startCpG < df.endCpG]

    if df.shape[0] != nr_orig:
        fixed_blocks_path = splitextgz(args.blocks_file)[0] + '.clean.bed.gz'
        eprint('Dumping cleaned blocks file to', fixed_blocks_path)
        df.to_csv(fixed_blocks_path, header=None, index=None, sep='\t', compression='gzip')

    return df

def validate_block_file(df):
    # file has >=5 columns
    # no duplicated lines
    # sorted by 4'th column (startCpG)
    # no overlap between regions
    return

# def is_blocks_file_nice(df):
    # nice = True
    # # see if the blocks file is "nice"

    # # startCpG and endCpG is monotonically increasing
    # if !pd.Index(df['startCpG']).is_monotonic_increasing():
        # eprint('startCpG is not monotonically increasing')
        # return False
    # if !pd.Index(df['endCpG']).is_monotonic_increasing():
        # eprint('endCpG is not monotonically increasing')
        # return False

    # if not (np.diff(df['endCpG'] - df['startCpG']) > 0).all():
        # eprint('Some blocks don\'t cover any CpGs')

def load_blocks_file(blocks_path):
    if not op.isfile(blocks_path):
        eprint('Invalid blocks file:', blocks_path)
        raise IllegalArgumentError('No blocks file')

    names = ['chr', 'start', 'end', 'startCpG', 'endCpG']
    df = pd.read_csv(blocks_path, sep='\t', usecols=[0, 1, 2, 3, 4], header=None, names=names)
    return df
    # starts = df['startCpG']
    # ends = df['endCpG']

    # # blocks file is valid
    # if not (np.diff(df['endCpG'] - df['startCpG']) >= 0).all():
        # eprint('Invalid blocks file: endCpG>=startCpG is required')
        # raise IllegalArgumentError('Invalid blocks file')



def slow_but_robust(df, args):
    for beta_path in args.input_files:

        try:
            res = np.zeros((df.shape[0], 2), dtype=int)

            # load beta file:
            data = load_beta_data(beta_path)

            # reduce to blocks:
            i = -1
            for _, row in df.iterrows():
                i += 1
                startCpG = row[3]
                endCpG = row[4]
                line = np.sum(data[startCpG - 1:endCpG - 1, :], axis=0)
                res[i, :] = line

            # dump to file
            suff = '.lbeta' if args.lbeta else '.bin'
            out_name = splitext(splitext(basename(args.blocks_file))[0])[0]
            out_name = splitext(basename(beta_path))[0] + suff
            out_name = op.join(args.out_dir, out_name)

            trim_to_uint8(res, args.lbeta).tofile(out_name)
            eprint(out_name)

            if args.bedGraph:
                with np.errstate(divide='ignore', invalid='ignore'):
                    beta_vals = res[:, 0] / res[:, 1]
                    eprint(beta_vals.shape, df.shape)
                # beta_vals[res[:, 1] == 0] = np.nan
                df['beta'] = beta_vals
                df.to_csv(out_name.replace(suff, '.bedGraph'), sep='\t',
                          index=None, header=None, na_rep=-1,
                          float_format='%.2f')

        except Exception as e:
            print('Failed with beta', beta_path)
            print('Exception:', e)


    return

def main():
    """
    Collapse beta file to blocks binary file, of the same beta format
    """

    args = parse_args()
    files = args.input_files
    validate_files_list(files, '.beta')

    # load blocks:
    df = load_blocks_file(args.blocks_file)

    if args.slow:
        slow_but_robust(df, args)
        return

    df = blocks_cleanup(df, args)

    blocks_bins, filtered_indices = get_bins(df, args.genome)

    with Pool() as p:
        for beta_path in files:
            params = (args, blocks_bins,
                      filtered_indices, beta_path, df[['chr', 'start', 'end']])
            p.apply_async(apply_filter_wrapper, params)
        p.close()
        p.join()

    # for beta_path in files:
    #     reduced_data = apply_filter_wrapper(beta_path, args.blocks, args.cov_thresh)


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('input_files', nargs='+', help='one or more beta files')
    parser.add_argument('-b', '--blocks_file', help='blocks path', default=default_blocks_path)
    parser.add_argument('-o', '--out_dir', help='output directory. Default is "."', default='.')
    parser.add_argument('-l', '--lbeta', action='store_true', help='Use lbeta file (uint16) instead of bin (uint8)')
    parser.add_argument('--bedGraph', action='store_true', help='output a text file in addition to binary file')
    parser.add_argument('--debug', '-d', action='store_true')
    parser.add_argument('--slow', action='store_true', help='slower but more robust implementation. Use this flag when input blocks file contains overlapping blocks')
    parser.add_argument('--genome', help='Genome reference name. Default is hg19.', default='hg19')

    return parser.parse_args()


if __name__ == '__main__':
    main()
