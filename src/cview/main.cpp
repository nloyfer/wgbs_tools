//
// Created by nloyfer on 5/27/19.
//

#include "cview.h"


/**
 *  Main
 *
 */

/*
 * Input Arguments Parsering class
 */
class InputParser {
public:
    InputParser(int &argc, char **argv) {
        for (int i = 1; i < argc; ++i)
            this->tokens.emplace_back(std::string(argv[i]));
    }

    const std::string &getCmdOption(const std::string &option) const {
        std::vector<std::string>::const_iterator itr;
        itr = std::find(this->tokens.begin(), this->tokens.end(), option);
        if (itr != this->tokens.end() && ++itr != this->tokens.end()) {
            return *itr;
        }
        static const std::string empty_string;
        return empty_string;
    }

    bool cmdOptionExists(const std::string &option) const {
        return std::find(this->tokens.begin(), this->tokens.end(), option)
               != this->tokens.end();
    }

private:
    std::vector <std::string> tokens;
};


std::string get_param_str(InputParser &input, std::string name, std::string defval) {
    std::string param_str = input.getCmdOption(name);
    if (!param_str.empty()) {
        return param_str;
    }
    return defval;
}


//void print_help() {
//    std::cerr << "" << std::endl;
//}

int main(int argc, char *argv[]) {

    //if ((argc < 3)) {
        //std::cerr << "Usage: cview [(-s SITES) | () | ()] [-d] [-b BLOCKS_PATH] [-l MIN_LEN] [--gzip]\n"
                  //<< "-l       Minimal sites in read to consider. Default is l=5.\n"
                  //<< "-r       Ranges of methylation average, in [0,1]. For example: 0,0.2001,0.8,1.\n"
                  //<< "-b       Blocks file. No header. Columns are:\n"
                  //<< "         <chr, start, end, startCpG, endCpG>. Default is s150.\n"
                  //<< "         File may be gzipped. Note sometimes bgzip causes problems.\n"
                  //<< "-d       Debug mode. Only use the first 2,500 blocks in the blocks file."
                  //<< std::endl;
        //return -1;
    //}

    InputParser input(argc, argv);

//    bool help = input.cmdOptionExists("-h");
//    if (help) {
//        print_help();
//        return 1;
//    }
    //std::string pat_path = std::string(argv[1]);
    //std::string output_path = std::string(argv[2]);
    //std::string blocks_path = get_param_str(input, "-b", DEFAULT_BLOCKS);
    //std::string range_str = get_param_str(input, "-r", "");
    //if (!range_str.size()) {
        //std::cout << "You must provide a range" << std::endl;
        //return 1;
    //}
    int min_cpgs = std::stoi(get_param_str(input, "--min_cpgs", "1"));

    bool debug = input.cmdOptionExists("-d");
    bool verbose = input.cmdOptionExists("-v");
    bool strict = input.cmdOptionExists("--strict");
    bool strip = input.cmdOptionExists("--strip");

    try {
        Cview(strict, strip, min_cpgs, debug, verbose).parse();
    }
    catch (std::exception &e) {
        std::cerr << e.what() << std::endl;
        return -1;
    }
    return 0;
}

