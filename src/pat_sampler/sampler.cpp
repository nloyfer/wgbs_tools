#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <chrono>
#include <random>
#include <stdexcept>      // std::invalid_argument


//  g++ -std=c++11 sampler.cpp -o pat_sampler


std::vector<std::string> line2tokens(std::string &line) {
    /** Break string line to tokens, return it as a vector of strings */
    std::vector<std::string> result;
    std::string cell;
    std::stringstream lineStream(line);
    while(getline(lineStream, cell, '\t'))
        result.push_back(cell);

    // validate line
    if (result.size() < 4)
        throw std::invalid_argument("too few columns in file, line: " + line);
    return result;
}

void output_read(std::vector <std::string> &tokens) {
    unsigned int i = 0;
    for (; i < tokens.size() - 1; i++) {
        std::cout << tokens[i] << "\t";
    }
    std::cout << tokens[i] << std::endl;
}


void consider_line(std::vector<std::string> &tokens, double rate, int reps) {
    int read_count = std::stoi(tokens[3]) * reps;

    // set seed
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);

    std::binomial_distribution<int> distribution (read_count, rate);

    int newcount = distribution(generator);
    if (newcount) {
        tokens[3] = std::to_string(newcount);
        output_read(tokens);
    }
}


/** main - sample reads from stdin (pat format). Output them to stdout */
int main( int argc, char *argv[])
{
    if (argc < 2){
        std::cerr << "Usage: pat_sampler RATE [REPS]" << std::endl;
        return -1;
    }

    try {
        double rate = std::stof(argv[1]);
        int reps = 1;
        if (argc > 2) {
            reps = std::stoi(argv[2]);
        }

        // Sample from stdin
        for (std::string line_str; std::getline(std::cin, line_str);) {
            if (line_str.empty()) { continue; } // skip empty lines
            std::vector<std::string> tokens = line2tokens(line_str);
            consider_line(tokens, rate, reps);
        }
    }
    catch(std::exception &e) {
        std::cerr << "failed sampling from pat" << std::endl;
        std::cerr << e.what() << std::endl;
        return -1;
    }
}
