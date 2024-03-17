#!/usr/bin/python3 -u

import re
import os.path as op
import numpy as np
from utils_wgbs import load_borders, load_beta_data, validate_file_list, color_text, \
        beta2vec, catch_BrokenPipeError, drop_dup_keep_order
from genomic_region import GenomicRegion

FULL_SQUARE = '\u25A0'
# FULL_SQUARE = '\u2588'
MISSING_VAL_SIGN = ' '
NR_CHARS_PER_FNAME = 50
MISSING_VAL = '.'


class BetaVis:
    def __init__(self, args):
        self.gr = GenomicRegion(args)
        self.start, self.end = self.gr.sites
        self.args = args

        # drop duplicated files, while keeping original order
        self.files = drop_dup_keep_order(args.input_files)

        # load raw data:
        self.dsets = self.load_data()

        # load borders:
        self.borders = load_borders(args.blocks_path, self.gr, args.genome)

        # Generate colors dictionary
        self.num2color_dict = generate_colors_dict(args.color_scheme)

        self.print_all()
        if self.args.plot:
            self.plot_all()

    def load_data(self):
        # raw table from *beta files:
        dsets = np.zeros((len(self.files), self.end - self.start, 2))
        for i, fpath in enumerate(self.files):
            dsets[i] = load_beta_data(fpath, self.gr.sites)
        return dsets

    def build_vals_line(self, data):

        # build a list of single character values.
        with np.errstate(divide='ignore', invalid='ignore'):
            vec = np.round((data[:, 0] / data[:, 1] * 10), 0).astype(int)  # normalize to range [0, 10)
        vec[vec == 10] = 9
        vec[data[:, 1] < self.args.min_cov] = -1
        vals = [MISSING_VAL if x == -1 else str(int(x)) for x in vec]

        # insert borders:
        if self.borders.size:
            vals = np.insert(vals, self.borders, '|')


        return self.color_vals(vals)

    def color_vals(self, vals):
        # join vals to a string line and color it:
        line = ''.join(vals)
        if not self.args.no_color:
            line = color_text(line, self.num2color_dict, scheme=self.args.color_scheme)
            if self.args.heatmap:
                line = re.sub('m[0-9]', 'm' + FULL_SQUARE * 1, line)
                line = re.sub('\.', MISSING_VAL_SIGN, line)
        return line

    def print_all(self):
        print(self.gr)

        # set the fixed number of characters for fpath names:
        fname_len = min(NR_CHARS_PER_FNAME, max([len(op.basename(op.splitext(f)[0])) for f in self.files])) + 1

        for dset, fpath in zip(self.dsets, self.files):
            line = self.build_vals_line(dset)
            adj_fname = op.splitext(op.basename(fpath))[0][:fname_len].ljust(fname_len)
            print(adj_fname + ': ' + line)

        if self.args.colorbar:
            digits = '0123456789'
            print('colorbar')
            print(self.color_vals(digits))
            if self.args.heatmap:
                print(digits)

    def plot_all(self):
        import matplotlib.pyplot as plt

        fname_len = min(NR_CHARS_PER_FNAME, max([len(op.basename(op.splitext(f)[0])) for f in self.files]))
        ticks = [op.splitext(op.basename(f))[0][:fname_len].ljust(fname_len) for f in self.files]

        r = np.concatenate([beta2vec(d).reshape((1, -1)) for d in self.dsets])

        plt.imshow(1 - r, cmap='RdYlGn')
        # insert borders:
        if self.borders.size:
            plt.vlines(self.borders - .5, -.5, len(self.files) - .5)

        plt.yticks(np.arange(len(self.files)), ticks)
        if self.args.title:
            plt.title(self.args.title)
        if self.args.output is not None:
            plt.savefig(self.args.output)
        plt.show()


def generate_colors_dict(scheme=16):
    if scheme == 16:
        colors = [
            "01;92",  # bold light green
            "92",  # light green
            "32",  # green
            "32",  # green
            "34",  # blue
            "34",  # blue
            "02;31",  # dark red
            "02;31",  # dark red
            "31",  # red
            "01;31"  # bold red
        ]
    else:
        colors = [10, 47, 70, 28, 3, 3, 202, 204, 197, 196]
    return {str(i): colors[i] for i in range(10)}


def main(args):
    validate_file_list(args.input_files) #, '.beta')
    try:
        BetaVis(args)
    except BrokenPipeError:
        catch_BrokenPipeError()
