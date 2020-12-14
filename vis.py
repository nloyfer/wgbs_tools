#!/usr/bin/python3 -u

import argparse
from utils_wgbs import splitextgz, add_GR_args, default_blocks_path
from beta_vis import main as beta_vis_main
from pat_vis import main as pat_vis_main
# from pat_vis2 import main as pat_vis_main_d

def parse_args():  # todo: seperate args parsing for beta and pat
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('input_files', nargs='+', help='A pat.gz file or one or more beta files')
    parser.add_argument('-d', '--dists', action='store_true', help='print results with distances (kind of log scale)')
    parser.add_argument('-t', '--title', help='A text to be printed before the results.')
    parser.add_argument('-o', '--output', help='beta vis: save plot to file')
    parser.add_argument('-b', '--blocks_path', help='Display blocks borders. If [-b] is specified with no '
                                                    'blocks path, default blocks are used.',
                        nargs='?', const=default_blocks_path, default=False)
    parser.add_argument("--no_color", action='store_true', help='Print without colors.')
    #parser.add_argument('--debug', action='store_true', help='debug')
    parser.add_argument('--strict', action='store_true', help='Truncate reads that start/end outside the given region. '
                                                              'Only relevant for pat files.')
    parser.add_argument('--strip', action='store_true',
                        help='pat: Remove trailing dots (from beginning/end of reads).')
    parser.add_argument('--max_reps', '-m', type=int, default=10,
                        help='Pat vis: Display a read at most "max_reps" times, '
                             'if it is repeating itself. [10]')
    parser.add_argument('--min_len', type=int, default=1,
                        help='Pat vis: Display only reads covering at least MIN_LEN CpG sites [1]')
    parser.add_argument('--color_scheme', '-cs', type=int, default=256,
                        help='beta vis: Color scheme. Possible values: 16 or 256 [256]')
    parser.add_argument('--plot', action='store_true', help='beta vis: plot results in a heatmap.')
    parser.add_argument('--no_dense', action='store_true',
                        help='pat: Do not squeeze multiple reads to every line.\n'
                             'Each read appears in a different line.')
    parser.add_argument('--uxm', type=float, default=None,
                        help='Pat vis: Float between 0 and 1 where reads with methylation proportion above'
                             ' this value will be displayed as fully methylated, reads with unmethylated CpG site proporiton below 1 - value will be displayed as fully unmethylated,'
                             ' or otherwise as X. ')
    add_GR_args(parser, required=True)
    return parser.parse_args(), parser


def main():
    """
    Visualize wgbs files
    Possible inputs:
        - a pat.gz file
        - One or more beta files
    """

    args, parser = parse_args()
    if args.uxm and not (args.uxm <= 1 and args.uxm >= 0):
        parser.error("uxm value must be between 0 and 1")
    file_type = splitextgz(args.input_files[0])[1]

    # print title
    if args.title:
        print('{}'.format(args.title))

    if file_type in ('.beta', '.bin'):
        beta_vis_main(args)
    elif file_type == '.pat.gz':
        # if args.tmp:
            # pat_vis_main_d(args)
            # return
        pat_vis_main(args)
    else:
        print('Unsupported file type:', file_type)


if __name__ == '__main__':
    main()
