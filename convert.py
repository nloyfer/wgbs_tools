#!/usr/bin/python3 -u

import argparse
from utils_wgbs import add_GR_args, eprint
from genomic_region import GenomicRegion


def parse_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    add_GR_args(parser)

    # todo: allow bed file input
    # parser.add_argument('-L', '--bed_file',
    #                     help='convert all regions in a bed file')
    args = parser.parse_args()
    return args


def main():
    """
    Convert genomic region to CpG index range and vise versa
    """
    args = parse_args()

    # if args.bed_file and (args.region or args.sites):
    #     eprint('-L, -s and -r are mutually exclusive')
    #     return

    # bed_wrapper = BedFileWrap(args.bed_file) if args.bed_file else None
    gr = GenomicRegion(args)
    print(gr)


if __name__ == '__main__':
    main()
