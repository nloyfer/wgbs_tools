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
#include <fstream>
//#include <string_view>
#include <cstdio>
#include <memory>
#include <stdexcept>

#define MAX_PAT_LEN 300
#define MAX_READ_LEN 2000

struct reads_stats {
    int nr_pairs = 0;
    int nr_empty = 0;
    int nr_short = 0;
    int nr_invalid = 0;
    int nr_bad_conv = 0;
};

class snp_patter {
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
    long snp_pos;
    char snp_let1;
    char snp_let2;
    std::unordered_map<int, int> dict;
    std::string genome_ref;
    reads_stats readsStats;
    int line_i = 0;
    clock_t tick = clock();
    bool is_paired_end = false;
    bool blueprint = false;
    bool first_line(std::string &line);

    snp_patter(long sp, char sl1, char sl2):
            snp_pos(sp), snp_let1(sl1), snp_let2(sl2) {}

    char compareSeqToRef(std::string &seq, std::string &ref, bool reversed, std::string &meth_pattern, int start_pos);
    void print_stats_msg();
    void print_progress();
    std::string clean_CIGAR(std::string seq, std::string CIGAR);
    char samLineToPatVec(std::vector<std::string> tokens);
    char proc2lines(std::vector<std::string> tokens1, std::vector<std::string> tokens2);
    void action(std::string samFilepath);
    void proc_sam_in_stream(std::istream& in);

    void initialize_patter(std::string &line_str);

};

std::vector<std::string> line2tokens(std::string &line);
void print_vec(std::vector<std::string> &vec);

#endif //FAST_PAT_PATTER_H
