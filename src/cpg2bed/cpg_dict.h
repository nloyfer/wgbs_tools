
#ifndef CPG_DICT_H
#define CPG_DICT_H

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>      // std::stringstream
//#include <string>
//#include <iomanip>      // std::setprecision
//#include <string_view>
#include <array>        // std::array
//#include <cstdio>
//#include <memory>
//#include <stdexcept>
#include <algorithm>
#include <thread>

class cpg_dict {
public:
 
    // Current chromosome
    std::vector<std::string> chroms;
    std::vector<int> loci;
    std::vector<int> borders;

    cpg_dict(std::string dict_path, std::string chrom_size_path) {
        load_chroms(chrom_size_path); // load chrom sizes
        load_dict(dict_path); // load CpG->locus dict
    }
    std::string loc2chrom(int loc);
    void chrom_exists(std::string chrom);
    bool is_border(std::string chrom, int cpg_pos);
private:
    std::vector<int> sizes;
    int load_chroms(std::string chrom_size_path);
    void load_dict(std::string dict_path);
};


std::vector<std::string> line2tokens(std::string &line);
void print_vec(std::vector<std::string> &vec);
std::string exec(const char* cmd);


#endif //CPG_DICT_H
