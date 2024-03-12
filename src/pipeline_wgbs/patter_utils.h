//
// Created by nloyfer on 9/14/18.
//

#ifndef PATTER_UTILS_H
#define PATTER_UTILS_H

#include <iostream>
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
void print_vec(const std::vector<std::string> &vec);
void print_vec(const std::vector<int> &vec);
void output_vec(const std::vector <std::string> &tokens);
std::string addCommas(const int num);


std::vector<int> split_by_comma(std::string str_line);
std::vector<float> split_float_by_comma(std::string str_line);
std::vector<std::string> split_by_semicolon(std::string str_line);

bool hasEnding(std::string const &fullString, std::string const &suffix);
bool is_number(const std::string& s);

std::string exec(const char* cmd);

bool are_paired(std::vector <std::string> tokens1,
                std::vector <std::string> tokens2);
std::string reverse_comp(std::string seq);
bool is_bottom(int samflag, bool is_paired_end);
std::string clean_CIGAR(std::string seq, std::string CIGAR);

void strip_read(std::vector <std::string> &tokens);
int strip_pat(std::string &pat);
std::vector <std::string> pack_pat(std::string chrom, 
                                   int start_site, 
                                   std::string meth_pattern);
std::vector <std::string> merge_PE(std::vector<std::string> l1, 
                                std::vector<std::string> l2);


/***************************************************
 *                                                 *
 *                 Useful symbols                  *
 *                                                 *
 **************************************************/

const char UNKNOWN = '.';
const char METH = 'C';
const char UNMETH = 'T';
const std::string TAB = "\t";

/***************************************************
 *                                                 *
 *        Input Arguments Parsering class          *
 *                                                 *
 **************************************************/
class InputParser{
public:
    InputParser (int &argc, char **argv){
        for (int i=1; i < argc; ++i)
            this->tokens.emplace_back(std::string(argv[i]));
    }
    // Get option 
    const std::string& getCmdOption(const std::string &option) const{
        std::vector<std::string>::const_iterator itr;
        itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
        if (itr != this->tokens.end() && ++itr != this->tokens.end()){
            return *itr;
        }
        static const std::string empty_string;
        return empty_string;
    }
    // Get boolean option 
    bool cmdOptionExists(const std::string &option) const{
        return std::find(this->tokens.begin(), this->tokens.end(), option)
               != this->tokens.end();
    }
    // Get option and set default if does not exist
    std::string getOptionWithDefault(const std::string &option, const std::string &defval) {
        std::string param_str = this->getCmdOption(option);
        if (param_str.empty()) { 
            return defval;
        }
        return param_str; 
    }
    // Get a required option
    std::string getRequiredOption(const std::string &option) {
        std::string param_str = this->getCmdOption(option);
        if (param_str.empty()) { 
            throw std::invalid_argument("Missing required argument: " + option);
        }
        return param_str;
    }
private:
    std::vector <std::string> tokens;
};


#endif //PATTER_UTILS_H

