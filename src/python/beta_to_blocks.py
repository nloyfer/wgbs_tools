#!/usr/bin/python3 -u

import argparse
import sys
import os.path as op
from multiprocessing import Pool
import pandas as pd
import numpy as np
from utils_wgbs import load_beta_data, trim_to_uint8, \
                        IllegalArgumentError, add_multi_thread_args, \
                        splitextgz, validate_file_list, validate_single_file, \
                        eprint, validate_out_dir, COORDS_COLS5

def b2b_log(*args, **kwargs):
    print('[ wt beta_to_blocks ]', *args, file=sys.stderr, **kwargs)

######################################################
#                                                    #
#  Loading and Parsing Blocks file                   #
#                                                    #
######################################################

def is_block_file_nice(df):

    msg = ''
    # no empty blocks (noCpGs):
    # no missing values (NAs)
    if df[['startCpG', 'endCpG']].isna().values.sum() > 0:
        msg = 'Some blocks are empty (NA)'
    # no (startCpG==endCpG)
    elif not (df['endCpG'] - df['startCpG'] > 0).all():
        msg = 'Some blocks are empty (startCpG==endCpG)'
    # blocks are sorted
    # startCpG and endCpG are monotonically increasing
    elif not np.all(np.diff(df['startCpG'].values) >= 0):
        msg = 'startCpG is not monotonically increasing'
    elif not np.all(np.diff(df['endCpG'].values) >= 0):
        msg = 'endCpG is not monotonically increasing'
    # no duplicated blocks
    elif df.shape[0] != df.drop_duplicates().shape[0]:
        msg = 'Some blocks are duplicated'
    # no overlaps between blocks
    elif not (df['startCpG'][1:].values - df['endCpG'][:df.shape[0] - 1].values  >= 0).all():
        msg = 'Some blocks overlap'
    if msg:
        return False, msg
    return True, ''


def load_blocks_file(blocks_path, anno=False, nrows=None):
    # validate blocks_path
    validate_single_file(blocks_path)

    try:
        # see if blocks_path has a header:
        peek_df = pd.read_csv(blocks_path, sep='\t', nrows=1, header=None, comment='#')
        header = None if str(peek_df.iloc[0, 1]).isdigit() else 0

        names = COORDS_COLS5.copy()
        if anno:
            names += ['anno', 'gene']
        if len(peek_df.columns) < len(COORDS_COLS5):
            msg = f'Invalid blocks file: {blocks_path}. less than {len(names)} columns.\n'
            msg += f'Run wgbstools convert -L {blocks_path} -o OUTPUT_REGION_FILE to add the CpG columns'
            raise IllegalArgumentError(msg)
        elif len(peek_df.columns) < len(names):  # no annotations columns
            names = COORDS_COLS5

        # load
        # dtypes = {'chr':str, 'start', 'end', 'startCpG', 'endCpG'}
        dtypes = {'startCpG':'Int64', 'endCpG':'Int64'}
        df = pd.read_csv(blocks_path, sep='\t', usecols=range(len(names)), dtype=dtypes,
                         header=header, names=names, nrows=nrows, comment='#')

        # blocks start before they end - invalid file
        dfnona = df.dropna()    # allow blocks with missing values
        if not ((dfnona['endCpG'] -  dfnona['startCpG']) >= 0).all():
            raise IllegalArgumentError(f'Invalid CpG columns in blocks file {blocks_path}')

        if dfnona.shape[0] == df.shape[0]:
            df['startCpG'] = df['startCpG'].astype(int)
            df['endCpG'] = df['endCpG'].astype(int)

    except pd.errors.ParserError as e:
        eprint(f'Invalid input file.\n{e}')
        return pd.DataFrame()
    except pd.errors.EmptyDataError as e:
        eprint(f'Empty blocks file.\n{e}')
        return pd.DataFrame()

    return df


######################################################
#                                                    #
#    Perform the reduction                           #
#                                                    #
######################################################


def fast_method(data, df):
    df -= (df.startCpG.iloc[0] - 1)
    block_bins = np.sort(np.unique(df.values.flatten()))[:-1]
    filtered_indices = np.isin(block_bins, df['startCpG'])
    return np.add.reduceat(data, block_bins - 1)[filtered_indices]


def slow_method(data, df):
    reduced_data = np.zeros((df.shape[0], 2), dtype=int)
    for i, row in df.iloc[:, 3:5].iterrows():
        startCpG, endCpG = row.astype('Int64')
        if pd.isna(startCpG):
            reduced_data[i, :] = [0, 0]
            continue
        reduced_data[i, :] = np.sum(data[startCpG - 1:endCpG - 1, :], axis=0)
    return reduced_data


def reduce_data(beta_path, df, is_nice):
    if is_nice:
        df = df[['startCpG', 'endCpG']].astype(int)
        start = df['startCpG'].values[0]
        end = df['endCpG'].values[df.shape[0] - 1]
        return fast_method(load_beta_data(beta_path, (start, end)), df)

    return slow_method(load_beta_data(beta_path), df)


def collapse_process(beta_path, df, is_nice, lbeta=False, out_dir=None, bedGraph=False):
    try:
        # load beta file:
        reduced_data = reduce_data(beta_path, df.reset_index(drop=True), is_nice)
        if out_dir is None:
            return reduced_data
        return dump(df, reduced_data, beta_path, lbeta, out_dir, bedGraph)

    except Exception as e:
        b2b_log('Failed with beta', beta_path)
        b2b_log('Exception:', e)


######################################################
#                                                    #
#    Dump results                                    #
#                                                    #
######################################################

def dump(df, reduced_data, beta_path, lbeta, out_dir, bedGraph):

    bin_table = trim_to_uint8(reduced_data, lbeta)
    # dump to binary file
    name = op.splitext(op.basename(beta_path))[0]
    suff = '.lbeta' if lbeta else '.bin'
    prefix = op.join(out_dir, name)
    bin_table.tofile(prefix + suff)
    b2b_log(prefix + suff)

    # dump to bed
    if bedGraph:
        with np.errstate(divide='ignore', invalid='ignore'):
            df['beta'] = reduced_data[:, 0] / reduced_data[:, 1]
        df['coverage'] = reduced_data[:, 1]
        df[['chr', 'start', 'end', 'beta', 'coverage']].to_csv(prefix + '.bedGraph', sep='\t',
                  index=None, header=None, na_rep=-1, float_format='%.2f')


def filter_existing_files(files, out_dir, lbeta):
    files_to_process = []
    suff = '.lbeta' if lbeta else '.bin'
    for beta in files:
        prefix = op.join(out_dir, splitextgz(op.basename(beta))[0])
        if not op.isfile(prefix + suff):
            files_to_process.append(beta)
        else:
            b2b_log(f'Skipping {beta}. Use -f flag to overwrite')
    return files_to_process


def main():
    """
    Collapse beta file to blocks binary file, of the same beta format
    """

    args = parse_args()
    files = args.input_files
    validate_file_list(files)
    validate_out_dir(args.out_dir)

    if not args.force:
        files = filter_existing_files(files, args.out_dir, args.lbeta)

    # load blocks:
    # b2b_log('load blocks...')
    df = load_blocks_file(args.blocks_file)
    is_nice, msg = is_block_file_nice(df)
    if not is_nice:
        b2b_log(msg)
    params = [(b, df, is_nice, args.lbeta, args.out_dir, args.bedGraph)
              for b in files]
    if args.debug:
        _ = [collapse_process(*k) for k in params]
    else:
        p = Pool(args.threads)
        p.starmap(collapse_process, params)
    p.close()
    p.join()


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('input_files', nargs='+', help='one or more beta files')
    parser.add_argument('-b', '--blocks_file', help='blocks path', required=True)
    parser.add_argument('-o', '--out_dir', help='output directory. Default is "."', default='.')
    parser.add_argument('-l', '--lbeta', action='store_true', help='Use lbeta file (uint16) instead of bin (uint8)')
    parser.add_argument('--bedGraph', action='store_true', help='output a text file in addition to binary file')
    parser.add_argument('--force', '-f', action='store_true', help='Overwrite existing files if existed')
    parser.add_argument('--debug', '-d', action='store_true')
    add_multi_thread_args(parser)

    return parser.parse_args()


if __name__ == '__main__':
    main()
