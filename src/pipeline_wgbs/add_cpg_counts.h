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
#include <iomanip>      // std::setprecision
#include <unordered_map>
#include <sstream>      // std::stringstream
#include <array>        // std::array
#include "patter_utils.h"


#define MAX_PAT_LEN 300
#define MAX_READ_LEN 1000

struct reads_stats {
    int nr_pairs = 0;
    int nr_empty = 0;
    int nr_invalid = 0;
    int nr_short = 0;
};

struct ReadOrient { // OT or OB
    char ref_chr;
    char unmeth_seq_chr;
    int shift;
    int mbias_ind;
};

class patter {
public:
    std::string ref_path;
    std::string chr;
    int bsize = 0;  // position of the last CpG in the current chromosome
    bool* conv;
    clock_t tick = clock();
    std::string region;
    std::unordered_map<int, int> dict;
    std::string genome_ref;
    reads_stats readsStats;
    std::string TAGNAMETYPE = "YI:Z:";
    std::string PATTAGNAMETYPE = "XP:Z:";
    std::vector <std::string> dummy_tokens;
    int clip_size = 0;
    bool print_pat =false;
    std::string bed_file;
    std::string bed_ext_file;
    int min_cpg = 0;
    int line_i = 0;
    bool is_paired_end = false;

    ReadOrient OT{'C', 'T', 0, 0};
    ReadOrient OB{'G', 'A', 1, 1};

    struct MethylData {
        int countMethyl; int countUnmethyl;
        int originalIndex;
        std::string pattern;
    };


    patter(std::string refpath,  std::string rgn, int min_len, int clip, bool print_pat, std::string b_file, std::string b_ext_file): 
        ref_path(refpath), region(rgn), min_cpg(min_len), clip_size(clip), print_pat(print_pat), bed_file(b_file), bed_ext_file(b_ext_file) {}

    void load_genome_ref();
    std::vector<int> load_genome_helper(std::string r, std::string cmd);


    bool first_line(std::string &line);
    void print_stats_msg();
    void print_progress();
    int locus2CpGIndex(int locus);
    std::vector<std::string> samLineToPatVec(std::vector<std::string> tokens);
    void action_sam(std::string samFilePath);
    std::string samLineMethyldataMakeString(std::string originalLine, patter::MethylData md);
    void proc_sam_in_stream(std::istream& in);
    patter::MethylData merge_and_count_methyl_data(std::vector <std::string> l1, std::vector <std::string> l2);
    void proc1line(std::vector<std::string> &tokens1, std::string line1);
    void procPairAddMethylData(std::vector<std::string> tokens1, std::vector<std::string> tokens2,
                               std::string line1, std::string line2);
    int compareSeqToRef(std::string &seq, int start_locus, int samflag, std::string &meth_pattern);

    void initialize_patter(std::string &line_str);

    void proc_pair_sam_lines(std::string &line1,
                             std::string &line2);
        
    int update_conv(std::vector<int> loci, int start, int end, int counter);
};


#endif //FAST_PAT_PATTER_H

