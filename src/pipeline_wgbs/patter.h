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

#define MAX_PAT_LEN 300
#define MAX_READ_LEN 1000

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
    int offset = 0;
    //bool conv[1000000] = {false};
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
    bool is_paired_end = false;
    std::vector <std::string> dummy_tokens;
    bool first_line(std::string &line);

    mbias_ss mbias_OT[2];
    mbias_ss mbias_OB[2];

    patter(std::string refpath, std::string rgn, std::string mb, int mc, int clip):
            ref_path(refpath), region(rgn), mbias_path(mb), min_cpg(mc), clip_size(clip) {}
    ~patter() {delete[] conv;}
    void load_genome_ref();
    int find_cpg_inds_offset();
    std::vector<long> fasta_index();

    int compareSeqToRef(std::string &seq, int start_locus, int samflag, std::string &meth_pattern);
    void print_stats_msg();
    void dump_mbias();
    void print_progress();
    int locus2CpGIndex(int locus);
    std::string clean_CIGAR(std::string seq, std::string CIGAR);
    std::vector<std::string> samLineToPatVec(std::vector<std::string> tokens);
    void proc2lines(std::vector<std::string> &tokens1, std::vector<std::string> &tokens2);
    void proc1line(std::vector <std::string> &tokens1);
    void action();

    void initialize_patter(std::string &line_str);

};

std::vector<std::string> line2tokens(std::string &line);
void print_vec(std::vector<std::string> &vec);

#endif //FAST_PAT_PATTER_H
