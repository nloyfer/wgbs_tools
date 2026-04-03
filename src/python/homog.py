#!/usr/bin/python3 -u

import argparse
import os.path as op
import sys
import subprocess
from io import StringIO
import numpy as np
import pandas as pd
from beta_to_blocks import load_blocks_file
from utils_wgbs import IllegalArgumentError, homog_tool, main_script, \
        validate_file_list, COORDS_COLS5, validate_local_exe, \
        mkdirp, delete_or_skip, pretty_name


def homog_log(*args, **kwargs):
    print('[ wt homog ]', *args, file=sys.stderr, **kwargs)


def _blocks_are_sorted(df):
    """Return True if blocks are sorted by startCpG (and endCpG as tiebreaker)."""
    starts = df['startCpG'].values
    ends = df['endCpG'].values
    for i in range(1, len(starts)):
        if starts[i] < starts[i - 1]:
            return False
        if starts[i] == starts[i - 1] and ends[i] < ends[i - 1]:
            return False
    return True


def _validate_blocks_ignore_sort(df):
    """Check block validity except sort order. Overlapping blocks are allowed. Returns (ok, msg)."""
    if df[['startCpG', 'endCpG']].isna().values.sum() > 0:
        return False, 'Some blocks are empty (NA)'
    if not (df['endCpG'] - df['startCpG'] > 0).all():
        return False, 'Some blocks are empty (startCpG==endCpG)'
    return True, ''


######################################################
#                                                    #
#    wrap the c++ tool                               #
#                                                    #
######################################################


def trim_uxm_to_uint8(data, nr_bits):
    data = data.copy()
    if nr_bits == 16:
        dtype = np.uint16
    else:
        dtype = np.uint8
    max_val = 2 ** nr_bits - 1
    big_inds = np.argwhere(data.max(axis=1) > max_val).flatten()
    data[big_inds, :] = data[big_inds, :] / data.max(axis=1)[big_inds][:, None] * max_val
    res = data.astype(dtype)
    return res


def ctool_wrap(pat, name, blocks_path, rates_cmd, view_full, inclusive,
               sort_blocks=False, verbose=False):
    if view_full:
        cmd = f'gunzip -c {pat}'
    else:
        cmd = f'{main_script} cview {pat} -L {blocks_path}'
    cmd += f' | {homog_tool} -b {blocks_path} -n {name} {rates_cmd}'
    if inclusive:
        cmd += ' --inclusive'
    if sort_blocks:
        cmd += ' --sort_blocks'
    se = None if verbose else subprocess.PIPE
    txt = subprocess.check_output(cmd, shell=True, stderr=se).decode()
    if verbose:
        homog_log(cmd)
    names = list('UXM')
    df = pd.read_csv(StringIO(txt), sep='\t', header=None, names=names)
    if df.values.sum() == 0:
        homog_log(f' [ {name} ] WARNING: all zeros!')
    return df


def homog_process(pat, blocks, args, outdir, prefix,
                   sort_blocks=False, nodump=False):
    # generate output path:
    name = pretty_name(pat)
    if prefix is None:
        prefix = op.join(outdir, name)
    opath = prefix + '.uxm'
    if not args.binary:
        opath += '.bed.gz'
    if not delete_or_skip(opath, args.force):
        homog_log(f'skipping {name}. Use -f to overwrite')
        return

    # generate rate_cmd:
    l = args.rlen
    rate_cmd = f' -l {l} -r '
    if args.thresholds:
        rate_cmd += f'0,{args.thresholds},1'
    else:
        th1 = round(1 - (l - 1) / l, 3) + 0.001
        th2 = round((l - 1) / l, 3)
        rate_cmd += f'0,{th1},{th2},1 '

    # for a long marker file (>5K blocks),
    # parse the whole pat file instead of running "cview -L BED"
    view_full = blocks.shape[0] > 5e3

    df = ctool_wrap(pat, name, args.blocks_file, rate_cmd, view_full,
                    args.inclusive, sort_blocks=sort_blocks, verbose=args.verbose)

    if sort_blocks:
        # C++ output is in sorted CpG order — restore to original block order
        # Build mapping: for each row in sorted output, find its original position
        sorted_starts = blocks['startCpG'].values.argsort(kind='stable')
        inv_order = np.argsort(sorted_starts, kind='stable')
        df = df.iloc[inv_order].reset_index(drop=True)

    df = blocks.merge(
        pd.concat([blocks[COORDS_COLS5].reset_index(drop=True), df], axis=1),
        how='left', on=COORDS_COLS5,
    )

    if nodump:
        return df
    if args.binary:
        trim_uxm_to_uint8(df[list('UXM')].values, args.nr_bits).tofile(opath)
    else:
        df.to_csv(opath, sep='\t', header=None, index=None)
    return df


def parse_outdir_prefix(args):
    outdir = args.out_dir
    prefix = args.prefix
    if prefix is not None:
        outdir = op.dirname(prefix)
    if not outdir:
        outdir = '.'
    mkdirp(outdir)
    return outdir, prefix


def main():
    """
    Generage homog files. Given a blocks file and pat[s],
    count the number of U,X,M reads for each block for each file
    """

    args = parse_args()
    if args.nr_bits not in (8 , 16):
        raise IllegalArgumentError('nr_bits must be in {8, 16}')
    if args.rlen < 2:
        raise IllegalArgumentError('rlen must be >= 2')
    if args.thresholds is not None:
        th = args.thresholds.split(',')
        if not len(th) == 2: # and th[0].is_number():
            raise IllegalArgumentError('Invalid thresholds')
        th = float(th[0]), float(th[1])
        if not 1 > th[1] > th[0] > 0:
            raise IllegalArgumentError('Invalid thresholds')
    elif args.rlen == 2:
        raise IllegalArgumentError('for rlen==2, --thresholds must be specified')
    # make sure homog tool is valid:
    validate_local_exe(homog_tool)

    pats = args.input_files
    validate_file_list(pats, '.pat.gz')

    outdir, prefix = parse_outdir_prefix(args)

    # load blocks:
    blocks_df = load_blocks_file(args.blocks_file)

    # validate blocks (allowing unsorted)
    ok, msg = _validate_blocks_ignore_sort(blocks_df)
    if not ok:
        homog_log(msg)
        raise IllegalArgumentError(f'Invalid blocks file: {args.blocks_file}')

    # check if blocks need sorting (C++ binary handles it via --sort_blocks)
    sort_blocks = not _blocks_are_sorted(blocks_df)
    if sort_blocks:
        homog_log(f'WARNING: blocks file is not sorted by startCpG. '
                  f'C++ binary will sort internally.')

    for pat in sorted(pats):
        homog_process(pat, blocks_df, args, outdir, prefix,
                      sort_blocks=sort_blocks)


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('input_files', nargs='+', help='one or more pat files')
    parser.add_argument('-b', '--blocks_file', help='blocks path', required=True)
    output_parser = parser.add_mutually_exclusive_group(required=False)
    output_parser.add_argument('-o', '--out_dir', help='output directory. Default is "."')
    output_parser.add_argument('-p', '--prefix', help='output prefix')
    parser.add_argument('--force', '-f', action='store_true', help='Overwrite files if exist')
    parser.add_argument('--inclusive', action='store_true', help='consider the whole read. Opposite of "strict"')
    parser.add_argument('--verbose', '-v', action='store_true')
    parser.add_argument('--binary', action='store_true', help='Output binary files (uint8)')
    parser.add_argument('--genome', help='Genome reference name.')
    parser.add_argument('--nr_bits', type=int, default=8,
            help='For binary output, specify number of bits for the output format - 8 or 16. ' \
                 '(e.g. 8 statnds for uint8, which means values are trimmed to [0, 255])')
    parser.add_argument('--thresholds', '-t',
            help='UXM thresholds, LOW,HIGH. E.g, "0.3334,0.666".\n')
    parser.add_argument('--rlen', '-l', type=int, default=3,
            help='Minimal read length (in CpGs) to consider. Default is 3')
    parser.add_argument('--debug', '-d', action='store_true')
    return parser.parse_args()


if __name__ == '__main__':
    main()
