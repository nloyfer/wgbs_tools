//
// Created by nloyfer on 5/27/19.
//

#ifndef CVIEW_H
#define CVIEW_H

#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <array>
#include <algorithm> // find in main
#include "../pipeline_wgbs/patter_utils.h"

#define SEP "\t"

struct Block {
    int start;
    int end;
    int count;
};

class Cview {
    int32_t *counts;

    bool debug;
    bool verbose;
    std::vector<Block> borders;
    int nr_blocks = 0;
    int cur_block_ind = 0;
    int min_cpgs = 0;
    bool strip;
    bool strict;
    std::string blocks_path;
    std::string sites;

    int read_blocks();

    void dump(int *data, int width, std::string out_path);

    int proc_line(std::vector <std::string> tokens);

    int blocks_helper(std::istream &instream);
    void output_vec(std::vector <std::string> &tokens);
    void pass_read(std::vector <std::string> &tokens);
    void strip_read(std::vector <std::string> &tokens);

public:
    Cview(std::string bpath, std::string in_sites, bool in_strict, 
            bool in_strip, int in_min_cpgs, bool deb, bool ver):
        blocks_path(bpath), sites(in_sites), strict(in_strict), 
        strip(in_strip), min_cpgs(in_min_cpgs), debug(deb), verbose(ver) {}

    ~Cview() {}

    void parse();
};

#endif //CVIEW_H
