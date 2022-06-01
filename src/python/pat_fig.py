#!/usr/bin/python3 -u

import os
import sys
import os.path as op
import argparse
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pat_vis import PatVis, int2str, str2int
from vis import pat_args
from utils_wgbs import add_GR_args, eprint, validate_file_list, \
        drop_dup_keep_order


def parse_args():
    parser = argparse.ArgumentParser()
    pat_args(parser)
    add_GR_args(parser, required=True, no_anno=True)
    parser.add_argument('pats', nargs='+')
    parser.add_argument('--outpath', '-o', required=True,
            help='Output path (e.g., path to a PDF or PNG file)')
    parser.add_argument('--top', type=int, default=1000,
            help='Output at most TOP reads for each pat file [Default: 1000]')
    parser.add_argument('--col_wrap', type=int, default=5,
            help='Wrap the columns at this width [Default: 5]')
    parser.add_argument('--space_cols', type=int, default=1,
            help='Add space between columns [Default: 1]')
    parser.add_argument('--space_rows', type=int, default=4,
            help='Add space between rows [Default: 4]')
    parser.add_argument('--circle_size', type=float, default=1.0,
            help='Size of the circles, relative to baseline. E.g., '
                 'for 1.1 increases circle size by 10%% [Default: 1.0]')
    parser.add_argument('--line_width', type=float, default=1.0,
            help='Line width of circle borders and strikethroughs '
                 'relative to baseline. E.g, 1.1 increases line width '
                 ' by 10%% [Default 1.0]')
    parser.add_argument('--font_size', type=float, default=1.0,
            help='Font size, relative to baseline. E.g., '
                 'for 1.1 increases font size by 10%% [Default: 1.0]')
    parser.add_argument('--title',
            help='Title for the figure. Default is details about the region')
    parser.add_argument('--fig_height', type=int, default=20)
    parser.add_argument('--blocks_path')
    parser.add_argument('--red_green', action='store_true')
    return parser.parse_args()


def get_strikes_coords(kf):
    kf[kf < 2] = 0
    kf[kf > 1] = 1
    z = np.zeros((kf.shape[0], 1))
    dif = np.diff(np.hstack([z, kf, z]))
    return np.hstack([np.argwhere(dif==1), np.argwhere(dif==-1)])[:, [0,1,3]].T



def plot(tf, headers, gr, args):

    # setup figure
    height, width = tf.shape
    fig = plt.figure(figsize=(args.fig_height * (width / height), args.fig_height))
    ax = fig.add_subplot(111)
    ax.set_ylim((-1, height + 1 + 3))
    ax.set_xlim((-1, width + 1))

    # strikes
    hly, xmins, xmaxs = get_strikes_coords(tf.copy())
    lw = args.line_width
    ax.hlines(height - hly, xmin=xmins - .5, xmax=xmaxs - .5, lw=0, zorder=-2)

    # set marker size and line widths
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    msize = (bbox.width / width * 43) * args.circle_size
    lw = msize / 5 * args.line_width
    ax.hlines(height - hly, xmin=xmins - .5, xmax=xmaxs - .5, lw=lw, color='black', zorder=-1)

    # circles
    def plot_circles(simb, color):
        x, y = np.argwhere(tf == simb).T[::-1]
        ax.plot(x, height - y, 'o', markersize=msize,
                markeredgewidth=lw,
                markeredgecolor='black', c=color)

    # TODO: 1. allow UXM and maybe SNPs, not just C's and T's
    #       2. allow color customization
    if args.red_green:
        plot_circles(3, 'red')
        plot_circles(4, 'green')
    else:
        plot_circles(3, 'yellow')
        plot_circles(4, 'blue')

    # headers
    fsize = msize * 1.5 * args.font_size
    for trio in headers:
        ax.text(*trio, color='black', fontsize=fsize)

    # title
    title = args.title
    if not title:
        title = str(gr).replace('\t', ' ')
    plt.title(title, size=fsize*1.2)

    plt.axis('off')
    plt.savefig(args.outpath)


def validate_args(args):
    validate_file_list(args.pats, '.pat.gz')
    def enforce_positive(var, name):
        if var <= 0:
            eprint('[wt vis] Invalid {name} flag: must be positive')
            exit()
    enforce_positive(args.col_wrap, 'col_wrap')
    enforce_positive(args.space_rows, 'space_rows')
    enforce_positive(args.space_cols, 'space_cols')
    enforce_positive(args.circle_size, 'circle_size')
    enforce_positive(args.font_size, 'font_size')
    enforce_positive(args.line_width, 'line_width')
    fig_suffs = tuple(plt.gcf().canvas.get_supported_filetypes().keys())
    if not args.outpath.endswith(fig_suffs):
        eprint(f'[wt vis] Invalid output flag: {args.outpath}.')
        eprint(f'         Must end with {fig_suffs}.')
        exit()


def pad(table, height=None, width=None):
    if height is None:
        height = table.shape[0]
    if width is None:
        width = table.shape[1]
    if height < table.shape[0]:
        eprint(f'[ wt vis] Error: unable to pad table with shape {table.shape}, height {height}')
        exit()
    if width < table.shape[1]:
        eprint(f'[ wt vis] Error: unable to pad table with shape {table.shape}, width {width}')
        exit()
    padz = np.zeros((height, width), dtype=int)
    padz[:table.shape[0], :table.shape[1]] = table
    return padz


def main():
    args = parse_args()
    validate_args(args)

    # drop duplicated files, while keeping original order
    pats = drop_dup_keep_order(args.pats)

    tables = []

    # load pat secions
    for pat in pats:
        pv = PatVis(args, pat)
        fr = pv.fullres
        if fr is None:
            t = np.zeros((0, 0))
        else:
            t = fr['int_table'][:args.top,]
        width = max(pv.gr.nr_sites + 1, t.shape[1]) + args.space_cols
        tables.append(pad(t, args.top + args.space_rows, width))

    tmp = []
    N = len(pats)
    step = args.col_wrap if args.col_wrap < N else N
    for i in range(0, N, step):
        row = np.hstack(tables[i:i+step])
        # trim tail (last zero rows), leave 2 lines
        nr_lines = np.argmin(row.sum(axis=1)) + args.space_rows
        tmp.append(row[:nr_lines, :])

    max_width = np.max([t.shape[1] for t in tmp])
    table = np.vstack([pad(t, None, max_width) for t in tmp])

    # trace positions of headers
    tc = []
    shifty = 0
    shiftx = 0
    s = 0
    for i in range(N):
        name = op.basename(pats[i])[:-7]
        tc.append((shiftx, table.shape[0] - shifty + 2, name))
        shiftx += tables[i].shape[1]
        if (i+1) % step == 0 and i > 0:
            shifty += tmp[s].shape[0]
            shiftx = 0
            s += 1

    if table.sum() == 0:
        eprint(f'[wt vis] WARNING: empty table for region {args.region}')
        return
    plot(table, tc, pv.gr, args)


if __name__ == '__main__':
    main()

