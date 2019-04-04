#!/usr/bin/python3 -u

from utils_wgbs import load_borders, load_dists, load_beta_data, validate_files_list, color_text, GenomeRefPaths
from genomic_region import GenomicRegion
import os.path as op
import numpy as np


NR_CHARS_PER_FNAME = 50
MISSING_VAL = '.'
# DISTS_STEPS = (10, 50, 100, 500, 1000, 5000)
DISTS_STEPS = [10 ** i for i in range(6)]


class BetaVis:
    # def __init__(self, sites, input_files, with_dists, blocks_path, region, title):
    def __init__(self, input_files, gr, with_dists, blocks_path, no_color, color_shceme=16, plot=False):
        self.no_color = no_color
        self.plot = plot
        self.gr = gr
        self.start, self.end = gr.sites
        self.nr_sites = self.end - self.start

        # load distances
        self.distances = self.load_pairwise_dists() if with_dists else None

        # drop duplicated files, while keeping original order
        seen = set()
        self.files = [x for x in input_files if not (x in seen or seen.add(x))]

        # load raw data:
        self.dsets = self.load_data()

        # load borders:
        self.borders = None if blocks_path is None else load_borders(blocks_path, gr)

        # Generate colors dictionary
        self.color_scheme = color_shceme
        self.num2color_dict = generate_colors_dict(color_shceme)

        self.print_all()
        if self.plot:
            self.plot_all()


    def load_pairwise_dists(self):
        """load distances between consecutive sites:"""
        pairwise_dists = load_dists(self.start, self.nr_sites, self.gr.genome)
        return [np.searchsorted(np.array(DISTS_STEPS), bp_dist) for bp_dist in pairwise_dists]

    def load_data(self):
        # raw table from *beta files:
        dsets = np.zeros((len(self.files), self.nr_sites, 2))
        for i, file in enumerate(self.files):
            dsets[i] = load_beta_data(file, (self.start, self.end))
        return dsets

    def build_vals_line(self, data):

        # build a list of single character values.
        data[data[:, 1] == 0] = np.nan
        data = (data[:, 0] / data[:, 1] * 9)        # normalize to range [0, 10)
        vals = [MISSING_VAL if np.isnan(x) else str(int(x)) for x in data]

        # insert distances:
        if self.distances is not None:
            vals = [c + ' ' * d for d, c in zip(self.distances, vals)]

        # insert borders:
        if self.borders is not None:
            # print(len(vals), len(self.borders))
            vals = np.insert(vals, self.borders, '|')

        # join vals to a string line and color it:
        line = ''.join(vals)
        if not self.no_color:
            line = color_text(line, self.num2color_dict, scheme=self.color_scheme)

        return line

    def print_all(self):
        print(self.gr)

        # set the fixed number of characters for file names:
        fname_len = min(NR_CHARS_PER_FNAME, max([len(op.basename(op.splitext(f)[0])) for f in self.files]))

        for dset, file in zip(self.dsets, self.files):
            line = self.build_vals_line(dset)
            adj_fname = op.splitext(op.basename(file))[0][:fname_len].ljust(fname_len)
            print(adj_fname + ': ' + line)


    def plot_all(self):
        import pandas as pd
        import matplotlib.pyplot as plt

        fname_len = min(NR_CHARS_PER_FNAME, max([len(op.basename(op.splitext(f)[0])) for f in self.files]))
        r = np.empty((len(self.files), self.nr_sites))
        i = 0
        ticks = []
        for dset, file in zip(self.dsets, self.files):
            vec = np.divide(dset[:, 0], dset[:, 1], where=dset[:, 1] > 0)  # normalize to range [0, 10)
            vec[dset[:, 1] == 0] = np.nan
            r[i, :] = vec
            i += 1
            ticks.append(op.splitext(op.basename(file))[0][:fname_len].ljust(fname_len))

        # insert distances:
        # if self.distances is not None:
        #     vals = [c + ' ' * d for d, c in zip(self.distances, vals)]

        plt.imshow(1 - r, cmap='RdYlGn')
        # insert borders:
        if self.borders is not None:
            plt.vlines(self.borders -.5, -.5, len(self.files)-.5)

        plt.yticks(np.arange(len(self.files)), ticks)
        plt.show()



def generate_colors_dict(scheme=16):
    if scheme == 16:
        colors = [
            "01;92",  # bold light green
            "92",     # light green
            "32",     # green
            "32",     # green
            "34",     # blue
            "34",     # blue
            "02;31",  # dark red
            "02;31",  # dark red
            "31",     # red
            "01;31"   # bold red
        ]
    else:
        colors = [10, 47, 70, 28, 3, 3, 202, 204, 197, 196]
    return dict([(str(i), colors[i]) for i in range(10)])


def main(args):
    validate_files_list(args.input_files, '.beta')
    BetaVis(args.input_files, GenomicRegion(args), args.dists, args.blocks_path, args.no_color,
            color_shceme=args.color_scheme, plot=args.plot)


