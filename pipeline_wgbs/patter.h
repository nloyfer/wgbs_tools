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

#define MAX_PAT_LEN 300

struct reads_stats {
    int nr_pairs = 0;
    int nr_empty = 0;
    int nr_invalid = 0;
};

class patter {
public:
    std::string ref_path;
    std::string chrom_sz_path;
    std::string chr;
    std::unordered_map<int, int> dict;
    std::string genome_ref;
    reads_stats readsStats;
    int line_i = 0;
    bool paired_end = false;
    bool first_line(std::string &line);


    patter(std::string refpath, std::string cspath): ref_path(refpath), chrom_sz_path(cspath) {}

    void load_genome_ref();
    int find_cpg_inds_offset();
    std::vector<long> fasta_index();


    void print_stats_msg();
    int locus2CpGIndex(int locus);
    std::string clean_seq(std::string seq, std::string CIGAR);
    std::vector<std::string> samLineToPatVec(std::vector<std::string> tokens);
    void merge_and_print(std::vector<std::string> l1, std::vector<std::string> l2);
    void proc2lines(std::vector<std::string> tokens1, std::vector<std::string> tokens2);
    void action();
};

std::vector<std::string> line2tokens(std::string &line);
void print_vec(std::vector<std::string> &vec);

#endif //FAST_PAT_PATTER_H
