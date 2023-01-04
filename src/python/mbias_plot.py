#!/usr/bin/python3 -u

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os.path as op
import matplotlib
import argparse
from utils_wgbs import eprint
matplotlib.rc('font', size=15)


def arange_table(in_df, rn):
    tf = in_df.copy()
    tf.columns = ['M', 'U']
    tf['N'] = tf.sum(axis=1)
    tf['meth'] = tf.M / tf.N
    tf['pos'] = tf.index + 1
    cov = np.nanmedian(tf.N[:50]) / 10
    tf.loc[tf.N < cov, ['meth', 'N']] = np.nan
    del tf['M'], tf['U']
    tf['read'] = 'read #' + str(rn)
    return tf

def load_and_arange(mp, by, PE):
    df = pd.read_csv(mp, sep='\t')
    dfl = df.iloc[:, [0, 1]]
    dfr = df.iloc[:, [2, 3]]
    dd = pd.concat([arange_table(dfl, 1), arange_table(dfr, 2)])
    r = dd.melt(id_vars=['pos', 'read'], value_vars=[by])
    if not PE:
        r = r[r['read'] != 'read #2']
        del r['read']
    return r


def plot_mbias(mtables, out_dir, PE=True):
    # load and arange tables
    assert len(mtables) == 2
    # make sure the first table is OB, and the second is OT
    if mtables[0].endswith('.OT.txt'):
        mtables.reverse()
    tables = [load_and_arange(m, l, PE) for l in ('meth', 'N') for m in mtables]

    # plot 4 figures
    fig, axes = plt.subplots(2, 2, figsize=(10,10))
    ((ax1, ax2), (ax3, ax4)) = axes
    hue = 'read' if PE else None
    for dd, ax in zip(tables, axes.flatten()):
        sns.lineplot(data=dd, x='pos', y='value', hue=hue, ax=ax)

    # set titles
    ax1.set_title('OT / CTOT' if PE else 'OB')
    ax2.set_title('OB / CTOB' if PE else 'OT')
    name = op.basename(mtables[0])[:-13]
    fig.suptitle(f'{name}.bam: Methylation Bias')

    # set axes labels
    ax1.set_ylabel('Average methylation')
    ax3.set_ylabel('Number of observations')
    ax2.set_ylabel('')
    ax1.set_xlabel('')
    ax2.set_xlabel('')
    ax4.set_ylabel('')

    # set shared axes 
    ax3.sharey(ax4)
    ax1.sharey(ax2)
    ax3.sharex(ax1)
    ax4.sharex(ax2)

    # set axes limits
    for ax in axes.flatten():
        ax.set_xlim(0, None)
    for ax in axes[0]:
        ax.set_ylim(0, 1)
    outpath = op.join(out_dir, name) + '.pdf'
    eprint(f'[wt bam2pat] [mbias] dumped figure to {outpath}')
    plt.savefig(outpath)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('mbias_tables', nargs=2)
    parser.add_argument('--out_dir', '-o', default='.')
    parser.add_argument('-PE', action='store_true')
    parser.add_argument('--debug', '-d', action='store_true')
    return parser.parse_args()


def main():
    args = parse_args()
    plot_mbias(args.mbias_tables, args.out_dir, args.PE)


if __name__ == '__main__':
    main()
