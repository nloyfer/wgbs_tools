//
// Created by nloyfer on 9/14/18.
//

#ifndef PATTER_UTILS_H
#define PATTER_UTILS_H

#include <iostream>
//#include <string>
#include <vector>
#include <algorithm>    // std::sort
#include <fstream>
#include <regex>
#include <unordered_map>
#include <sstream>      // std::stringstream
#include <iomanip>      // std::setprecision
#include <array>        // std::array
#include <cstdio>
#include <memory>
#include <stdexcept>

#define MAX_PE_PAT_LEN 300
#define MAX_READ_LEN 1000


std::vector<std::string> line2tokens(const std::string &line);
void print_vec(std::vector<std::string> &vec);
void print_vec(std::vector<int> &vec);
std::string addCommas(int num);
std::vector<int> split_by_comma(std::string str_line);
std::vector<std::string> split_by_semicolon(std::string str_line);

std::string exec(const char* cmd);

bool are_paired(std::vector <std::string> tokens1,
                std::vector <std::string> tokens2);
std::string clean_CIGAR(std::string seq, std::string CIGAR);
int strip_pat(std::string &pat);
std::vector <std::string> pack_pat(std::string chrom, 
                                   int start_site, 
                                   std::string meth_pattern);
std::vector <std::string> merge_PE(std::vector<std::string> l1, 
                                std::vector<std::string> l2);
std::string reverse_comp(std::string seq);

const char UNKNOWN = '.';
#endif //PATTER_UTILS_H

