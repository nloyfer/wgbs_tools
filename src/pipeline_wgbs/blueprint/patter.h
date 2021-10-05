//
// Created by nloyfer on 9/14/18.
//

#ifndef FAST_PAT_PATTER_H
#define FAST_PAT_PATTER_H

#include <iostream>
#include <ctime>
//#include <string>
#include <vector>
#include <algorithm>    // std::sort
#include <fstream>
#include <regex>
#include <unordered_map>
#include <sstream>      // std::stringstream
#include <iomanip>      // std::setprecision
//#include <string_view>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <array>        // std::array

#define MAX_PAT_LEN 300
#define MAX_READ_LEN 2000

struct reads_stats {
    int nr_pairs = 0;
    int nr_empty = 0;
    int nr_short = 0;
    int nr_invalid = 0;
    int nr_bad_conv = 0;
};


struct mbias_ss {
    int meth[MAX_READ_LEN] = {0};
    int unmeth[MAX_READ_LEN] = {0};
};

class patter {
public:
    /** path to reference FASTA */
    std::string ref_path;
 
    // Current chromosome
    std::string chr;

    /** CpG-Index offset of the current chromosome.
     * i.e., How many CpGs there are before the first CpG in the current chromosome */
    int chrom_offset;

    std::string mbias_path;
    int min_cpg = 0;
    std::unordered_map<int, int> dict;
    std::string genome_ref;
    reads_stats readsStats;
    int line_i = 0;
    clock_t tick = clock();
    bool is_paired_end = false;
    bool blueprint = false;
    bool first_line(std::string &line);

    mbias_ss mbias[2];

    patter(std::string refpath, int coff, bool bp, std::string mb, int mc):
            ref_path(refpath), chrom_offset(coff), blueprint(bp), mbias_path(mb), min_cpg(mc) {}

    void load_genome_ref();
    int find_cpg_inds_offset();
    std::vector<long> fasta_index();

    int compareSeqToRef(std::string &seq, std::string &ref, bool reversed, std::string &meth_pattern);
    void print_stats_msg();
    void dump_mbias();
    void print_progress();
    int locus2CpGIndex(int locus);
    std::string clean_CIGAR(std::string seq, std::string CIGAR);
    std::vector<std::string> samLineToPatVec(std::vector<std::string> tokens);
    void proc2lines(std::vector<std::string> tokens1, std::vector<std::string> tokens2);
    void action();

    void initialize_patter(std::string &line_str);

};

std::vector<std::string> line2tokens(std::string &line);
void print_vec(std::vector<std::string> &vec);

#endif //FAST_PAT_PATTER_H
