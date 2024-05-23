#!/usr/bin/python3 -u

import os.path as op
import numpy as np
import pandas as pd
from utils_wgbs import eprint, IllegalArgumentError

um_ind_dict = {'U': 0, 'M': 2}

def load_uxm(binpath, blocks_df, um, min_cov):
    data = np.fromfile(binpath, np.uint8).reshape((-1, 3))[blocks_df.index, :]
    covs = data.sum(axis=1)
    cond = covs > min_cov
    r = np.divide(data[:, um_ind_dict[um]], covs, where=cond)
    r[~cond] = np.nan
    return r.astype(float)

#######################################################
#                                                     #
#          Parse groups and bins                      #
#                                                     #
#######################################################

def load_gfile_helper(groups_file):
    # load and validate csv
    gf = pd.read_csv(groups_file, index_col=False, comment='#')
    if 'group' not in gf.columns:
        raise IllegalArgumentError('gropus file must have a column named "group"')
    # drop samples where include==False
    if 'include' in gf.columns:
        if gf['include'].dtype != bool:
            eprint('Invalid group file')
            raise IllegalArgumentError('Invalid group file. Include column must be boolean')
        gf = gf[gf['include']]
    # drop unnecessary columns
    gf = gf.rename(columns={gf.columns[0]: 'fname'})
    gf = gf[['fname', 'group']].dropna().reset_index(drop=True)
    return gf


def match_prefix_to_bin(prefixes, bins, suff=None, first=False):
    first = True
    full_paths = []
    missing_binaries = []
    included = []
    for prefix in prefixes:
        # look for bin files starting with current prefix:
        if suff is not None:
            results = [f for f in bins if op.basename(f) == prefix + suff]
        else:
            results = [f for f in bins if op.basename(f).startswith(prefix)]
        # make sure there is exactly one such bin file
        if not results:
            missing_binaries.append(prefix)
        elif len(results) > 1 and not first:
            eprint(f'Found multiple matches for prefix {prefix}:')
            eprint(results, sep='\n')
            raise IllegalArgumentError('Invalid binary files or group file')
        else:
            full_paths.append(results[0])
            included.append(op.basename(results[0]))

    # print information about missing binary files
    if missing_binaries:
        eprint(f'Error: {len(missing_binaries)} prefixes from groups file were not found in input bins:')
        for p in missing_binaries:
            eprint(p)
        raise IllegalArgumentError('groups file mismatch binary files')

    # check how many bins were not included in prefixes:
    excluded_bins = []
    for b in bins:
        if op.basename(b) not in included:
            excluded_bins.append(b)
    if excluded_bins:
        eprint(f'warning: ignoring {len(excluded_bins)} binary files not included in groups file')
        if len(excluded_bins) < 5:
            eprint(*excluded_bins, sep='\n')

    return full_paths
