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
#include <array>        // std::array
#include <cstdio>
#include <memory>
#include <stdexcept>

#include "patter_utils.h"

struct reads_stats {
    int nr_pairs = 0;
    int nr_empty = 0;
    int nr_short = 0;
    int nr_invalid = 0;
    int nr_bad_conv = 0;
};

struct ReadOrient { // OT or OB
    char ref_chr;
    char unmeth_seq_chr;
    int shift;
    int mbias_ind;
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
    std::string region;

    bool* conv;

    ReadOrient OT{'C', 'T', 0, 0};
    ReadOrient OB{'G', 'A', 1, 1};

    std::string mbias_path;
    int min_cpg = 0;
    int clip_size = 0;
    int bsize = 0;  // position of the last CpG in the current chromosome
    std::unordered_map<int, int> dict;
    std::string genome_ref;
    reads_stats readsStats;
    int line_i = 0;
    clock_t tick = clock();
    bool is_long = false;
    bool is_paired_end = false;
    bool is_nanopore = false;
    bool np_dot = false; // Nanopore: Does MM field starts with "C+m." or "C+m?"?
    float np_thresh = 0.667;
    std::vector <std::string> dummy_tokens;
    void first_line(std::string &line);

    mbias_ss mbias_OT[2];
    mbias_ss mbias_OB[2];

    patter(std::string refpath, std::string rgn, std::string mb, int mc, int clip, bool is_np, float np_th, bool is_lng):
            ref_path(refpath), region(rgn), mbias_path(mb), min_cpg(mc), 
            clip_size(clip), is_nanopore(is_np), np_thresh(np_th), is_long(is_lng) {}
    ~patter() {delete[] conv;}
    void load_genome_ref();
    std::vector<long> fasta_index();

    int compareSeqToRef(std::string &seq, int start_locus, int samflag, std::string &meth_pattern);
    int np_call_meth(std::string &seq, int start_locus, int samflag, std::string &meth_pattern);
    int np_call_meth(std::string &seq, std::vector<std::string> &np_fields, 
                     int start_locus, int samflag, std::string &meth_pattern);
    void parse_np_fields(std::vector<std::string> &np_fields, 
                                 std::vector<int> &MM_vals, 
                                 std::vector<int> &ML_vals);
    std::string parse_ONT(std::vector <std::string> tokens);
    void print_stats_msg();
    void dump_mbias();
    void print_progress();
    int locus2CpGIndex(int locus);
    std::vector<std::string> samLineToPatVec(std::vector<std::string> tokens);
    std::vector<std::string> np_samLineToPatVec(std::vector<std::string> tokens);
    void proc2lines(std::vector<std::string> &tokens1, std::vector<std::string> &tokens2);
    void proc1line(std::vector <std::string> &tokens1);
    void parse_reads_from_stdin();

    void initialize_patter(std::string &line_str);

};

std::vector<std::string> get_np_fields(std::vector <std::string> &tokens);
std::vector <std::string> make_pat_vec(std::string chrom, int start_site, 
                                       std::string meth_pattern);
std::string parse_ONT(std::vector <std::string> tokens);
void parse_np_fields(std::vector<std::string> &np_fields, 
                    std::vector<int> &MM_vals, 
                    std::vector<int> &ML_vals);

#endif //FAST_PAT_PATTER_H
