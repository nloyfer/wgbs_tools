//
// Created by nloyfer on 11/2/18.
//
#include <memory>
#include <stdexcept>
#include <limits>
#include "segmentor.h"
/**
 * prints
 */

//void print_mem(double *mem, int nr_sites, int max_cpg) {
    //std::cerr << "\ncost" << std::endl;

    //for (int i = 0; i < nr_sites; i++) {
        //for (int j = 0; j < max_cpg; j++) {
            //std::cerr << mem[i * max_cpg + j] << '\t';
        //}
        //std::cerr << std::endl;
    //}
//}

//void print_MT(double *M, int *T, int nr_sites) {
    //std::cerr << "\nM" << '\t' << 'T' << std::endl;
    //for (int i = 0; i < nr_sites + 1; i++) {
        //std::cerr << M[i] << '\t' << T[i] << std::endl;
    //}
//}

void print_borders(std::vector<int> borders){
    for (auto i = borders.rbegin(); i != borders.rend(); ++i )
        printf("%d ", *i);
    printf("\n");
}

void segmentor::load_dists(uint32_t *dists) {

    if (params.max_bp == 0) { return; }

    int j = 0;
    for (std::string line; std::getline(std::cin, line); j++) {
        dists[j] = std::stoi(line);
    }
    if (j != nr_sites) {
        std::cerr << "Error: nr_sites != number of loci: " << nr_sites << " != " << j << ". Try different chunck size!"<< std::endl;
        throw 0;
    }
}

std::vector<int> segmentor::traceback(const int *T) {
    std::vector<int> borders;
    int i = nr_sites;
    borders.push_back(i);
    while (i > 0) {
        borders.push_back(i = std::max(0, T[i]));
    }
    return borders;
}

void segmentor::dp(std::vector<float*> &all_data){
    /**
     * Compute each row of the cost matrix on the fly into a small ring
     * buffer, then fold it directly into the DP transition.
     *
     * Previously a full (nr_sites x max_cpg) cost matrix was materialized
     * ahead of time in cost_memoization(). With defaults (chunk_size=60000,
     * max_cpg=1000) that was ~480 MB per worker.
     *
     * The DP at step i needs block costs for the last `max_cpg` rows only
     * (rows in [i+1-max_cpg, i]), so we keep exactly that many rows and
     * map row k to slot (k & ring_mask). Ring size is rounded up to the
     * next power of two so the per-access index is `& mask` rather than
     * `% ring_size` (the modulo was a measurable slowdown in the DP inner
     * loop). After: ring_size * max_cpg doubles, ~8 MB.
     *
     * NOTE: the intermediate `ll_k` must stay `float` (not `double`) to
     * match the pre-refactor numerics exactly — even a tiny precision
     * difference here can nudge marginal block-merge decisions at the DP
     * boundary. Verified byte-identical on whole-chromosome regressions.
     *
     * Cost formula:
     *   p_mle = (#meth + pseudo_count) / (#total + 2*pseudo_count)
     *   score = #meth * log2(p_mle) + #unmeth * log2(1-p_mle)
     * summed over datasets.
     */
    int nr_dsets = (int) all_data.size();
    std::unique_ptr<float[]>    nmeth(new float[nr_dsets]);
    std::unique_ptr<float[]>    ntotal(new float[nr_dsets]);
    std::unique_ptr<uint32_t[]> dists(new uint32_t[nr_sites]);
    load_dists(dists.get());

    int ring_size = 1;
    while (ring_size < max_cpg) ring_size <<= 1;
    int ring_mask = ring_size - 1;
    std::unique_ptr<double[]> mem_ring(new double[(size_t)ring_size * max_cpg]());

    std::unique_ptr<double[]> M(new double[nr_sites + 1]());
    std::unique_ptr<int[]>    T(new int[nr_sites + 1]());

    const double NEG_INF_D = -std::numeric_limits<double>::infinity();
    const double NEG_INF_F = -std::numeric_limits<float>::infinity();

    for (int i = 0; i < (int)nr_sites; i++) {
        // Compute row i of the cost matrix into its ring slot.
        double *row = &mem_ring[(size_t)(i & ring_mask) * max_cpg];
        std::fill(row, row + max_cpg, 0.0);

        memset(nmeth.get(),  0, (size_t)nr_dsets * sizeof(float));
        memset(ntotal.get(), 0, (size_t)nr_dsets * sizeof(float));

        int window = std::min((int)nr_sites - i, max_cpg);
        for (int j = 0; j < window; j++) {

            if ( (dists[i + j] - dists[i] > params.max_bp) || (dists[i + j] < dists[i]) ) {
                row[j] = NEG_INF_F;
                continue;
            }

            double ll_sum = 0;
            for (int k = 0; k < nr_dsets; k++) {
                nmeth[k]  += all_data[k][(i + j) * 2];
                ntotal[k] += all_data[k][((i + j) * 2) + 1];
                float ntotal_k = ntotal[k];
                float nmeth_k  = nmeth[k];
                if (!ntotal_k) { continue; }

                float p_mle_k = (nmeth_k + pseudo_count) / (ntotal_k + (2 * pseudo_count));
                float ll_k = 0;
                if (p_mle_k > 0.0) {
                    ll_k += (nmeth_k * log2(p_mle_k));
                }
                if (p_mle_k < 1.0) {
                    ll_k += (ntotal_k - nmeth_k) * log2(1.0 - p_mle_k);
                }
                ll_sum += ll_k;
            }
            if (ll_sum) row[j] = ll_sum;
        }

        // DP step i: M[i+1] = max over k in [max(0, i+1-max_cpg), i+1)
        //   of M[k] + mem[k][i - k].
        double best_score = NEG_INF_F;
        int best_ind = -1;
        int start_k = std::max(0, i + 1 - max_cpg);
        for (int k = start_k; k < i + 1; k++) {
            double *kth = &mem_ring[(size_t)(k & ring_mask) * max_cpg];
            double tmp = M[k] + kth[i - k];
            if (tmp > best_score) {
                best_score = tmp;
                best_ind = k;
            }
        }
        M[i+1] = best_score;
        T[i+1] = best_ind;
    }

    std::vector<int> borders = traceback(T.get());
    print_borders(borders);
}

/**
 *  Reading from files
 */
void segmentor::read_beta_file(const char *beta_path, float *data){
    // load beta files section to a temporary array cdata
    auto cdata = new char[nr_sites * 2];
    std::ifstream infile;
    infile.open(beta_path, std::ios::binary | std::ios::in);
    infile.seekg(start * 2);
    infile.read(cdata, nr_sites * 2);
    infile.close();

    // cast values to float(uint8) and set data variable
    for (int i = 0; i < nr_sites * 2; i++) {
        data[i] = (float) (uint8_t) cdata[i];
    }
    delete [] cdata;

    // validate data: in a beta file, right column must be greater
    // or equal the left column
    for (int i = 0; i < nr_sites; i++) {
        if (data[2 * i] > data[2 * i + 1]) {
            std::cerr << "invalid data, i = " << i << ". data: ";
            std::cerr << data[2 * i] << ", " << data[2 * i + 1] << std::endl;
            std::cerr << "beta path: " << beta_path << std::endl;
            throw 0;
        }
    }

}


void segmentor::dp_wrapper(){
    start = (uint) params.start;
    nr_sites = (uint) params.nr_sites;
    max_cpg = params.max_cpg;
    pseudo_count = params.pseudo_count;

    // Load beta files:
    // all_data is vector of floats pointers of size (nr_betas * nr_sites * 2)
    std::vector<float*> all_data;
    for (const auto &beta_path: beta_paths) {
        auto data = new float[nr_sites * 2];
        read_beta_file(beta_path.c_str(), data);
        all_data.emplace_back(data);
    }

    dp(all_data);

    // Delete data from memory
    for (auto dataset: all_data) {
        delete [] dataset;
    }
}

