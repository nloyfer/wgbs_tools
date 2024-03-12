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
#include <string>
#include <stdlib.h>
#include <array>
#include "../pipeline_wgbs/patter_utils.h"

struct Block {
    int start;
    int end;
    int count;
};

int read_blocks(std::string block_data, std::vector<Block> &borders);

class Cview {
    std::vector<Block> borders;
    int nr_blocks = 0;
    int cur_block_ind = 0;
    int min_cpgs = 0;
    bool strip;
    bool strict;
    std::string blocks_path;
    std::string sites;
    bool verbose;

    int proc_line(std::vector <std::string> tokens);
    void pass_read(std::vector <std::string> &tokens);

public:
    Cview(std::string bpath, std::string in_sites, bool in_strict, 
            bool in_strip, int in_min_cpgs, bool ver):
        blocks_path(bpath), sites(in_sites), strict(in_strict), 
        strip(in_strip), min_cpgs(in_min_cpgs), verbose(ver) {}

    ~Cview() {}

    void parse();
};

#endif //CVIEW_H
