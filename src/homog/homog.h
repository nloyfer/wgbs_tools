//
// Created by nloyfer on 5/27/19.
//

#ifndef PATS_READS_LENS_PAT2RLEN_H
#define PATS_READS_LENS_PAT2RLEN_H

#include <vector>
#include <deque>
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
#include <climits>      // INT_MAX
#include "../pipeline_wgbs/patter_utils.h"

#define SEP "\t"
#define DEFAULT_LEN "5"

struct BlockData {
    int start;                    // startCpG
    int end;                      // endCpG
    std::string coords;           // "chr\tstart\tend"
    std::vector<int32_t> counts;  // bin counts, size = nr_bins
};

class Homog {

    bool debug;
    bool inclusive;
    std::string blocks_path;
    std::string name;
    std::string sname;
    std::string chrom;
    std::vector<float> range;
    int cur_block_ind = 0;
    int min_cpgs = 0;
    int nr_bins;

    // Sliding window of active blocks
    std::deque<BlockData> blocks;
    int flush_offset = 0;       // number of blocks already flushed
    FILE *blocks_pipe = nullptr;
    bool pipe_exhausted = false;
    int blocks_loaded = 0;      // total blocks loaded so far

    // Block streaming
    void open_blocks_pipe();
    void load_blocks_until(int max_start);

    // Flushing
    void flush_block(BlockData &b);
    void flush_through(int block_ind);
    void flush_remaining();

    // Deque accessors (absolute index → deque position)
    inline BlockData& block_at(int i) { return blocks[i - flush_offset]; }
    inline int block_start(int i) { return blocks[i - flush_offset].start; }
    inline int block_end(int i)   { return blocks[i - flush_offset].end; }
    inline int loaded_end()       { return flush_offset + (int)blocks.size(); }

    void update_m2(int block_ind, const std::string &pat, int count);
    int proc_line(const std::vector <std::string> &tokens);
    void update_block(int ind, const std::string &pat, const std::string &orig_pattern, int count);

public:
    Homog(std::string blocks_path, std::vector<float> range,
            int m_order, bool deb, std::string name,
            std::string chrom, bool inclusive);

    ~Homog();

    void parse();
};

#endif //PATS_READS_LENS_PAT2RLEN_H
