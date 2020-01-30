#include <boost/algorithm/string/predicate.hpp>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <iostream>

#define NR_SITES 28217448


//  g++ -std=c++11 stdin2beta.cpp -o stdin2beta


std::vector<std::string> line2tokens(std::string &line) {
    /**
     * Break string line to tokens, return it as a vector of strings
     */
    std::vector<std::string> result;
    std::string cell;
    std::stringstream lineStream(line);
    while(getline(lineStream, cell, '\t'))
        result.push_back(cell);
    return result;
}


class Pat2Beta {
    int *meth, *cover;
    int nr_sites = NR_SITES;
    std::string beta_path;
    std::string pat_path;
    void dumpbin();
    int proc_line(std::vector<std::string> tokens);

public:
    Pat2Beta(std::string pat_path, int in_nr_sites);
    ~Pat2Beta();
    void parse();
};


Pat2Beta::~Pat2Beta() {
    delete[] meth;
    delete[] cover;
}

Pat2Beta::Pat2Beta(std::string beta_output_path, int in_nr_sites) {

    if ((!boost::algorithm::ends_with(beta_output_path, ".beta"))) {
        throw std::invalid_argument( "Ouput path must end with .beta\n" + beta_output_path );
    }
    beta_path = beta_output_path;
    nr_sites = in_nr_sites;

    meth = new int[nr_sites];
    cover = new int[nr_sites];
}

void Pat2Beta::dumpbin() {
    /**
     * Assuming meth and cover arrays are already filled,
     * dump them in the binary beta file.
     * avoid overflow.
     */
    std::ofstream bofs;
    bofs.open(beta_path, std::ios::binary);
    if (!(bofs.is_open())) {
        std::cerr << "could not open binary output file " << beta_path << std::endl;
        return;
    }

    if (!bofs.is_open()) {
        return;
    }
    std::cerr << "dumping to binary file " << std::endl;
    int j = 0;
    for ( int i = 0; i < nr_sites; i++ ) {
        int meth_nr = meth[i];
        int cov_nr = cover[i];

        // squeeze results to range [0,255]:
        if (cov_nr > 255) {
            meth_nr = (int) round(((float) meth_nr / (float) cov_nr) * 255.0);
            cov_nr = 255;
        }

        char bin_line[2];
        bin_line[0] = (char) meth_nr;
        bin_line[1] = (char) cov_nr;
        bofs.write((const char *) bin_line, 2);
        j++;
    }
    bofs.close();
}

int Pat2Beta::proc_line(std::vector<std::string> tokens) {
    /**
     * Given one line of the form "chr1 1245345 CCT 3", update the meth and coverage arrays -
     * Increase them in the corresponding sites.
     */
    if (tokens.size() < 4) {
        throw std::invalid_argument( "Invalid site in input file. too few columns" );
    }
    else {
        int site = std::stoi(tokens[1]);
        std::string pattern = tokens[2];
        int count = std::stoi(tokens[3]);
        auto read_len = (int) pattern.length();

        if ( site + read_len - 1 > nr_sites ){
            std::cerr << "Invalid site: " << site << " + " << read_len << std::endl;
            throw std::invalid_argument( "Invalid site in input file" );
        }
        if ( site < 1 ){
            std::cerr << "Invalid site: " << site << std::endl;
            throw std::invalid_argument( "Invalid site in input file" );
        }

        for(int i = 0; i < read_len; i++ ) {
            if ( (pattern[i] == 'T') | (pattern[i] == 'C') ) {
                cover[site - 1 + i] += count;
                meth[site - 1 + i] += count * (( pattern[i] == 'C' ) ? 1 : 0);
            }
        }

        return 0;
    }
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

        // dump_bed
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
    int nr_sites = NR_SITES;
    if (argc < 2){
        std::cerr << "Usage: stdin2beta OUTPUT_BETA_PATH [NR_SITES]" << std::endl;
        return -1;
    }
    if (argc == 3) {
        nr_sites = std::stoi(argv[2]);
    }

    try {
        Pat2Beta *p = new Pat2Beta(std::string(argv[1]), nr_sites);
        p->parse();
        delete p;
    }
    catch(std::exception &e) {
        std::cerr << e.what() << std::endl;
        return -1;
    }
}
