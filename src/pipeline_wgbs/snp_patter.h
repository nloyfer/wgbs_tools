//
// Created by nloyfer on 9/14/18.
//

#ifndef FAST_PAT_PATTER_H
#define FAST_PAT_PATTER_H

#include <iostream>
#include <ctime>
#include <vector>
#include <algorithm>    // std::sort
#include <fstream>
#include <regex>
#include <unordered_map>
#include <sstream>      // std::stringstream
#include <iomanip>      // std::setprecision
#include <fstream>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <array>

struct reads_stats {
    int nr_pairs = 0;
    int nr_empty = 0;
    int nr_short = 0;
    int nr_invalid = 0;
    int nr_bad_conv = 0;
};

class snp_patter {
public:
    std::string chr;
    long snp_pos;
    char snp_let1;
    char snp_let2;
    int qual_filter;
    reads_stats readsStats;
    int line_i = 0;
    clock_t tick = clock();
    bool is_paired_end = false;
    std::vector <std::string> dummy_tokens;
    bool first_line(std::string &line);

    snp_patter(long sp, char sl1, char sl2, int qual_filter):
            snp_pos(sp), snp_let1(sl1), snp_let2(sl2), qual_filter(qual_filter) {}

    char compareSeqToRef(std::string &seq, bool bottom, std::string &qual_str, int start_pos);
    void print_stats_msg();
    void print_progress();
    char samLineToSNP(std::vector<std::string> tokens);
    char proc2lines(std::vector<std::string> tokens1, std::vector<std::string> tokens2);
    char proc1line(std::vector<std::string> &tokens1);
    void action(std::string samFilepath);
    void proc_sam_in_stream(std::istream& in);

    void initialize_patter(std::string &line_str);
};

#endif //FAST_PAT_PATTER_H
