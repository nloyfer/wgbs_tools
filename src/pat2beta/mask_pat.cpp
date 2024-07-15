#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <set>          // std::set
#include <string>
#include <iostream>
#include <stdexcept>      // std::invalid_argument
#include "../pipeline_wgbs/patter_utils.h"

// Load blocks
int read_bed(std::string bed_path, std::set<int> &sites_to_hide) {
    /**
     * Load bed (blocks) file  vector<int> borders_starts, borders_ends.
     */

    std::cerr << "[ mask_pat ] loading sites..." << std::endl;

    // choose the read command (cat/gunzip)
    std::string cmd = "cat " + bed_path;
    if ((bed_path.length() > 3) && (bed_path.substr(bed_path.length() - 3) == ".gz")) {
        cmd = "gunzip -c " + bed_path;
    }
    // load whole bed file to memory
    std::string bed_data = exec(cmd.c_str());
    if (bed_data.length() == 0) {
        throw std::invalid_argument("[ mask_pat ] Error: Unable to load bed:" + bed_path);
    }
    std::stringstream ss(bed_data);

    //Iterate lines
    std::vector <std::string> tokens;
    std::string line;
    int cur_start = 0, cur_end = 0;
    int bi = 0;
    while (std::getline(ss, line)) {

        // skip empty lines and comments
        if (line.empty() || (!(line.rfind("#", 0)))) { continue; }

        tokens = line2tokens(line);
        if (tokens.size() < 5) {
            std::cerr << "[ mask_pat ] Invalid blocks file format. ";
            std::cerr << "[ mask_pat ] Should be: chr start end startCpG endCpG\n";
            print_vec(tokens);
            throw std::invalid_argument("Invalid block format");
        }

        // skip header, if exists
        if ((bi == 0) && (tokens[0] == "chr")) { continue; }

        cur_start = std::stoi(tokens[3]);
        cur_end = std::stoi(tokens[4]);

        // If block is invalid, abort:
        if ((cur_end <= cur_start)) {
            std::cerr << "[ mask_pat ] Invalid block: " << cur_start << "\t" << cur_end << std::endl;
            throw std::invalid_argument("Invalid block: endCpG <= startCpG");
        } else if (cur_start < 1) {
            throw std::invalid_argument("Invalid block: startCpG < 1");
        }

        for (int i = cur_start; i < cur_end; i++) {
            sites_to_hide.insert(i);
        }
        bi++;
    }

    if (sites_to_hide.size() == 0) {
        std::cerr << "[ mask_pat ] Error while loading bed file. 0 sites found.\n";
        throw std::invalid_argument("");
     }
    std::cerr << "[ mask_pat ] loaded " << sites_to_hide.size() << " sites.\n";
    return 0;
}


int proc_line(std::vector<std::string> tokens, std::set<int> &sites_to_hide) {
    /** */
    if (tokens.size() < 4) {
        throw std::invalid_argument( "Invalid site in input file. too few columns" );
    }
    
    int site = std::stoi(tokens[1]);
    std::string pattern = tokens[2];
    int count = std::stoi(tokens[3]);
    auto read_len = (int) pattern.length();

    for(int i = 0; i < read_len; i++ ) {
        if (pattern[i] == UNKNOWN) { continue; }
        int cur_ind = site + i;
        if (sites_to_hide.find(site + i) != sites_to_hide.end()) {
            pattern[i] = UNKNOWN;
        }
    }
    tokens[2] = pattern;
    // strip pat and ignore this read if it's all dots
    int pos = strip_pat(tokens[2]);
    if (pos < 0) { return 0; }

    tokens[1] = std::to_string(pos + std::stoi(tokens[1]));
    output_vec(tokens);
    return 0;
}

void parse(std::string bed_path){

    try {
        std::set<int> sites_to_hide;
        read_bed(bed_path, sites_to_hide);

        int line_i = 0;
        for (std::string line_str; std::getline(std::cin, line_str);) {
            if (line_str.empty()) { continue; } // skip empty lines

            std::vector<std::string> tokens = line2tokens(line_str);
            if (tokens.size() < 4) {
                throw std::invalid_argument("too few columns in file, line " + std::to_string(line_i));
            } else if (!(tokens.empty())) {
                proc_line(tokens, sites_to_hide);
            } else {
                std::cerr << "something went wrong... tokens is empty" << std::endl;
            }
            line_i++;
        }
    }

    catch(std::exception &e) {
        std::cerr << "failed in mask_pat" << std::endl;
        std::cerr << e.what() << std::endl;
        return;
    }
}


/** main - mask out sites from pat reads */
int main( int argc, char *argv[])
{
    if (argc != 2){
        std::cerr << "Usage: mask_pat bed_path" << std::endl;
        return -1;
    }
    try {
        parse(std::string(argv[1]));
    }
    catch(std::exception &e) {
        std::cerr << e.what() << std::endl;
        return -1;
    }
}
