//
// Created by nloyfer on 11/2/18.
//
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

void segmentor::cost_memoization(std::vector<float*> &all_data){
    /**
     * Create a cost matrix of size (nr_sites * nr_sites)
     * mem[i, j] = j >= i: cost of the block [i,...,j]
     *             i > j: cost(i,j) = -inf
     * The cost of a block is calculated as follows:
     * First, find p_mle = (#meth + pseudo_count) / (#total + 2*pseudo_count)
     * Then the score is the log likelihood:
     *              (#meth) * log(p_mle) + (#unmeth) * log(1-p_mle)
     * */
    auto nr_dsets = (int) all_data.size();
    auto *nmeth = new float[nr_dsets];
    auto *ntotal = new float[nr_dsets];
    auto *dists = new uint32_t[nr_sites];
    double ll_sum, z;
    float p_mle_k, ll_k, ntotal_k, nmeth_k;
    
    load_dists(dists);

    for (int i = 0; i < nr_sites; i++) {
        memset(nmeth, 0, nr_dsets * sizeof(*nmeth));
        memset(ntotal, 0, nr_dsets * sizeof(*ntotal));

        int window = std::min((int)nr_sites - i, max_cpg);
        for (int j = 0; j < window; j++) {


            if ( (dists[i + j] - dists[i] > params.max_bp) || (dists[i + j] < dists[i]) ) {
                mem[i * max_cpg + j] = -std::numeric_limits<float>::infinity();
                continue;
            }

            ll_sum = 0;

            for (int k = 0; k < nr_dsets; k++) {

                nmeth[k] += all_data[k][(i + j) * 2];
                ntotal[k] += all_data[k][((i + j) * 2) + 1];
                ntotal_k = ntotal[k];
                nmeth_k = nmeth[k];
                if (!ntotal_k) { continue; }

                //log likelihood = nmeth * log(p_mle) + (ntotal - nmeth) * log(1 - p_mle)
                p_mle_k = (nmeth_k + pseudo_count) / (ntotal_k + (2 * pseudo_count));
                ll_k = 0;
                if (p_mle_k > 0.0) {
                    ll_k += (nmeth_k * log2(p_mle_k));
                }
                if (p_mle_k < 1.0) {
                    ll_k += (ntotal_k - nmeth_k) * log2(1.0 - p_mle_k);
                }
                // weighted average of the log likelihoods
                ll_sum += ll_k;
            }

            if (ll_sum) {
                mem[i * max_cpg + j] = ll_sum;
            }
        }
    }
    delete[] nmeth;
    delete[] ntotal;
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

    // cost memoization:
    mem = new double[nr_sites * max_cpg]();
    cost_memoization(all_data);
    /**
     * init M, T, all to zeros.
     *
     *  M[i] - The maximal score for segmenting sites 1,...,i
     *
     *  T - Traceback table.
     *  T[i] is the index of the last border we opened in the current segmentation
     */
    auto M = new double[nr_sites + 1]();
    auto T = new int[nr_sites + 1]();

    for (int i = 0; i < nr_sites; i++) {
        double best_score = -std::numeric_limits<float>::infinity();
        int best_ind = -1;
        int start = std::max(0, i + 1 - max_cpg);
        for (int k = start; k < i + 1; k++) {
            double tmp = M[k] + mem[k * max_cpg + i - k];
            if (tmp > best_score) {
                best_score = tmp;
                best_ind = k;
            }
        }
        M[i+1] = best_score;
        T[i+1] = best_ind;
    }
    std::vector<int> borders = traceback(T);
   // print_mem(mem, nr_sites, max_cpg);
    //print_MT(M, T, nr_sites);
    print_borders(borders);

    delete [] mem; delete [] M; delete [] T;
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

