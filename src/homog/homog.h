//
// Created by nloyfer on 5/27/19.
//

#ifndef PATS_READS_LENS_PAT2RLEN_H
#define PATS_READS_LENS_PAT2RLEN_H

#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <stdexcept>    // std::invalid_argument
#include <algorithm>    // std::find
#include <array>        // std::array
#include <memory>       // std::unique_ptr
#include "../pipeline_wgbs/patter_utils.h"

#define SEP "\t"
#define DEFAULT_LEN "5"

class Homog {
    int32_t *counts;


    bool debug;
    bool inclusive;
    std::string blocks_path;
    std::string name;
    std::string sname;
    std::string chrom;
    std::vector<int> borders_starts;
    std::vector<int> borders_ends;
    std::vector<int> borders_counts;
    std::vector<std::string> coords;
    std::vector<float> range;
    int nr_blocks = 0;
    int cur_block_ind = 0;
    int min_cpgs = 0;

    int nr_bins;

    int read_blocks();

    void dump(int *data, int width);
    void update_m2(int block_ind, std::string pat, int count);

    int proc_line(std::vector <std::string> tokens);

    void update_block(int ind, std::string pat, std::string orig_pattern, int count);

    int32_t* init_array(int len);

    int blocks_helper(std::istream &instream);

public:
    Homog(std::string blocks_path, std::vector<float> range,
            int m_order, bool deb, std::string name,
            std::string chrom, bool inclusive);

    ~Homog();

    void parse();
};

#endif //PATS_READS_LENS_PAT2RLEN_H
