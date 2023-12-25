#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <iostream>
#include <stdexcept>      // std::invalid_argument


std::vector<std::string> line2tokens(std::string &line) {
    /** * Break string line to tokens, return it as a vector of strings */
    std::vector<std::string> result;
    std::string cell;
    std::stringstream lineStream(line);
    while(getline(lineStream, cell, '\t'))
        result.push_back(cell);
    return result;
}


class Pat2Beta {
    int *meth, *cover;
    int nr_sites = 0;
    int start = 1;
    int end = 1;
    void dumpbin();
    int proc_line(std::vector<std::string> tokens);

public:
    Pat2Beta(int in_start, int in_end);
    ~Pat2Beta();
    void parse();
};


Pat2Beta::~Pat2Beta() {
    delete[] meth;
    delete[] cover;
}

Pat2Beta::Pat2Beta(int in_start, int in_end) {

    start = in_start;
    end = in_end;
    nr_sites = end - start;

    meth = new int[nr_sites];
    cover = new int[nr_sites];
}

void Pat2Beta::dumpbin() {
    for ( int i = 0; i < nr_sites; i++ ) {
        std::cout << meth[i] << " " << cover[i] << " "; 
    }
    std::cout << std::endl;
}

int Pat2Beta::proc_line(std::vector<std::string> tokens) {
    /**
     * Given one line of the form "chr1 1245345 CCT 3", update the meth and coverage arrays -
     * Increase them in the corresponding sites.
     */
    if (tokens.size() < 4) {
        throw std::invalid_argument( "Invalid site in input file. too few columns" );
    }
    
    int site = std::stoi(tokens[1]);
    std::string pattern = tokens[2];
    int count = std::stoi(tokens[3]);
    auto read_len = (int) pattern.length();

    if (( site + read_len - 1 < start ) || (site >= end)){
        // skip this read - ends too soon or starts too late
        return 0;
    }

    for(int i = 0; i < read_len; i++ ) {
        int cur_char = pattern[i];
        int cur_ind = site - start  + i;
        if ((cur_ind >= nr_sites) || (cur_ind < 0)) {
            continue;
        }
        if (! ( (cur_char == 'T') || (cur_char == 'C') ) ) {
            continue;
        }
        cover[cur_ind] += count;
        if (cur_char == 'C') {
            meth[cur_ind] += count; 
        }
    }
    return 0;
}

void Pat2Beta::parse(){

    try {

        int line_i = 0;
        for (std::string line_str; std::getline(std::cin, line_str);) {
            if (line_str.empty()) { continue; } // skip empty lines

            std::vector<std::string> tokens = line2tokens(line_str);
            if (tokens.size() < 4) {
                throw std::invalid_argument("too few columns in file, line " + std::to_string(line_i));
            }
            else if (!(tokens.empty())) {
                proc_line(tokens);
            }
            else {
                std::cerr << "something went wrong... tokens is empty" << std::endl;
            }
            line_i++;
        }

        dumpbin();
    }
    catch(std::exception &e) {
        std::cerr << "failed calculating beta" << std::endl;
        std::cerr << e.what() << std::endl;
        return;
    }
}


/** main - generate and dump a beta file from stdin pat */
int main( int argc, char *argv[])
{
    if (argc != 3){
        std::cerr << "Usage: stdin2beta startCpG endCpG" << std::endl;
        return -1;
    }

    try {
        int start = std::stoi(argv[1]);
        int end = std::stoi(argv[2]);

        Pat2Beta *p = new Pat2Beta(start, end);
        p->parse();
        delete p;
    }
    catch(std::exception &e) {
        std::cerr << e.what() << std::endl;
        return -1;
    }
}
