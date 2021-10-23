import argparse
import subprocess
import warnings
from genomic_region import GenomicRegion
from utils_wgbs import add_GR_args, MAX_PAT_LEN
import numpy as np
from scipy import stats


# def region2sites(region):
#     path = Path(op.realpath(__file__))
#     refdir = op.join(op.join(path.parent.parent.parent, 'references'), 'default')
#     refdir = str(Path(refdir).resolve())
#     dictpath = op.join(refdir, 'CpG.bed.gz')
#     regions_from = regions.split(":")[-1].split("-")
#     # find CpG indexes in range of the region:
#     cmd = f'tabix {dictpath} {region}'
#     # cmd += 'awk \'(NR==1){first=$3} {lbp=$2} END{print first"-"$3+1}\''
#
#     # if bp_tuple[1] equals exactly a loci of a CpG site, this site is *not* included
#     # e.g., in hg19, chr6:71046415-71046562 is 9718430-9718435 in sites
#     cmd += f"awk -v b={self.bp_tuple[1]} "
#     cmd += '\'(NR==1){first=$3} END{if ($2<b) {r+=1}; print first"-"$3+r}\''
#     res = subprocess.check_output(cmd, shell=True).decode()
#     # eprint(cmd)
#
#     if len(set(res.strip().split('-'))) == 1:
#         res = '-1'
#
#     # throw error if there are no CpGs in range
#     if res.strip() == '-1':
#         raise IllegalArgumentError(f'Invalid genomic region: {self.region_str}. No CpGs in range')
#
#     s1, s2 = self._sites_str_to_tuple(res)
#     # s2 += 1     # non-inclusive
#     return s1, s2


def read_pat_vis(pat_text, start):
    first_ind = 0
    max_ind = 0

    reads_list = []
    ind_list = []
    counts_list = []
    is_first = True
    # with open(filename, 'r') if filename is not "-" else sys.stdin as f:
    for line in pat_text.split("\n"):
        # line = line.strip()
        if len(line) > 1:
            tokens = line.split("\t")
            cur_ind = int(tokens[1])
            cur_end = cur_ind + len(tokens[-2])
            if cur_end <= start:
                continue
            if is_first:
                first_ind = cur_ind
                is_first = False
            last_ind = cur_ind
            ind_list.append(last_ind)
            reads_list.append(tokens[-2])
            counts_list.append(int(tokens[-1]))
            if cur_end > max_ind:
                max_ind = cur_end
    num_cpgs = max_ind - first_ind
    pat_matrix = np.chararray((sum(counts_list), num_cpgs))
    pat_matrix[:] = b'.'
    read_index = 0
    for ((read, read_ind), count) in zip(zip(reads_list, ind_list), counts_list):
        col_ind = read_ind - first_ind
        for _ in range(count):
            pat_matrix[read_index, col_ind:(col_ind + len(read))] = list(read)
            read_index += 1
    return pat_matrix


def em_pat_matrix(pat_matrix):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c_char = b'C'
        t_char = b'T'

        c_per_col = 1e-3 + sum(pat_matrix == b'C')
        t_per_col = 1e-3 + sum(pat_matrix == b'T')
        n_per_col = c_per_col + t_per_col

        num_reads = pat_matrix.shape[0]
        num_cpgs = pat_matrix.shape[1]

        p_c = np.zeros([2, num_cpgs])
        p_c[0, :] = 0.9
        p_c[1, :] = 0.1
        p_t = 1 - p_c
        l_p_c = np.log2(p_c)
        l_p_t = np.log2(p_t)
        p_alleles = np.array([0.5, 0.5])
        l_p_alleles = np.log2(p_alleles)
        ll = np.infty * -1
        # jabba = np.zeros([2, num_cpgs])
        # zabba = np.zeros([2, num_cpgs])

        ll_alleles = np.zeros([2, num_reads])
        # ll_alleles_c = np.zeros([2, num_reads])
        ll_delta = 100
        while ll_delta > 0:
            # for row_ind in range(num_reads):
            #     #l_p_alleles[0]*np.ones(num_reads) + np.matmul(l_p_c[0,:], (pat_matrix == c_char).T) + np.matmul(l_p_t[0, :], (pat_matrix == t_char).T)
            #     ll_alleles_c[0, row_ind] = l_p_alleles[0] + sum(l_p_c[0, pat_matrix[row_ind, :] == c_char]) \
            #                              + sum(l_p_t[0, pat_matrix[row_ind, :] == t_char])
            #     ll_alleles_c[1, row_ind] = l_p_alleles[1] + sum(l_p_c[1, pat_matrix[row_ind, :] == c_char]) + \
            #                              sum(l_p_t[1, pat_matrix[row_ind, :] == t_char])
            ll_alleles[0, :] = l_p_alleles[0] * np.ones(num_reads) + \
                               np.matmul(l_p_c[0, :], (pat_matrix == c_char).T) + \
                               np.matmul(l_p_t[0, :], (pat_matrix == t_char).T)
            ll_alleles[1, :] = l_p_alleles[1] * np.ones(num_reads) + \
                               np.matmul(l_p_c[1, :], (pat_matrix == c_char).T) + \
                               np.matmul(l_p_t[1, :], (pat_matrix == t_char).T)
            # read_assignments_o = np.argmax(ll_alleles_c, axis=0)
            read_assignments = np.argmax(ll_alleles, axis=0)
            new_ll = sum(ll_alleles[0, read_assignments == 0]) + sum(ll_alleles[1, read_assignments == 1])
            ll_delta = new_ll - ll
            ll = new_ll

            p_alleles = np.array([sum(read_assignments == 0), sum(read_assignments == 1)]) / num_reads
            l_p_alleles = np.log2(p_alleles)
            # for j in range(2):
            #     for i in range(num_cpgs):
            #         p_c[j, i] = 1e-3 + sum(pat_matrix[read_assignments == j, i] == c_char)
            #         p_t[j, i] = 1e-3 + sum(pat_matrix[read_assignments == j, i] == t_char)
            for z in range(2):
                p_c[z, :] = 1e-3 + sum(pat_matrix[read_assignments == z, :] == c_char)
                p_t[z, :] = 1e-3 + sum(pat_matrix[read_assignments == z, :] == t_char)
            totals = p_c + p_t
            l_p_c = np.log2(p_c / totals)
            l_p_t = np.log2(p_t / totals)

        bpi = 2 ** (ll / sum(n_per_col))
        print(f"LL: {ll} | {num_reads} reads | {int(round(sum(n_per_col)))} observed | BPI: {bpi}")
        return ll


def calc_initial_liklihood(pat_matrix):
    c_per_col = 1e-3 + sum(pat_matrix == b'C')
    t_per_col = 1e-3 + sum(pat_matrix == b'T')
    n_per_col = c_per_col + t_per_col
    l_p_c = np.log2(c_per_col / n_per_col)
    l_p_t = np.log2(t_per_col / n_per_col)
    ll0 = sum(sum(pat_matrix == b'C') * l_p_c + sum(pat_matrix == b'T') * l_p_t)
    bpi = 2 ** (ll0 / sum(n_per_col))
    print(f"LL: {ll0} | {pat_matrix.shape[0]} reads | {int(round(sum(n_per_col)))} observed | BPI: {bpi}")
    return ll0


def pull_pat_file(region, file_name):
    cmd = f'tabix {file_name} {region}'
    return subprocess.check_output(cmd, shell=True).decode()


def add_args():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('pat', help="The input pat file")
    add_GR_args(parser)
    # add_multi_thread_args(parser)
    return parser


def parse_args(parser):
    args = parser.parse_args()
    return args


def main():
    """
    Test whether region is bimodal
    """
    parser = add_args()
    args = parse_args(parser)

    gr = GenomicRegion(region=args.region)
    s1, s2 = gr.sites
    chrom = args.region.split(":")[0]
    new_start = max(1, s1 - MAX_PAT_LEN)
    pat_region = f"{chrom}:{new_start}-{s2}"
    pat_text = pull_pat_file(pat_region, args.pat)
    pat_mat = read_pat_vis(pat_text, s1)
    ll0 = calc_initial_liklihood(pat_mat)
    post_em_ll = em_pat_matrix(pat_mat)
    test_stat = 2 * np.log(2) * (post_em_ll - ll0)
    pv = 1 - stats.chi2.cdf(test_stat, pat_mat.shape[1] - 1)
    pv = "{:,.3e}".format(pv)
    print(f"pvalue: {pv}")


if __name__ == '__main__':
    main()