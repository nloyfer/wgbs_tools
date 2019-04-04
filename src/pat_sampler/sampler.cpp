#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <chrono>
#include <random>

//  g++ -std=c++11 sampler.cpp -o pat_sampler


std::vector<std::string> line2tokens(std::string &line) {
    /** Break string line to tokens, return it as a vector of strings */
    std::vector<std::string> result;
    std::string cell;
    std::stringstream lineStream(line);
    while(getline(lineStream, cell, '\t'))
        result.push_back(cell);

    // validate line
    if (result.size() != 4)
        throw std::invalid_argument("too many/few columns in file, line: " + line);
    return result;
}


void consider_line(std::vector<std::string> &tokens, double rate) {
    int count = std::stoi(tokens[3]);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);

    std::binomial_distribution<int> distribution (count, rate);
    int newcount = distribution(generator);
    if (newcount) {
        std::cout << tokens[0] << "\t" << tokens[1] << "\t" << tokens[2] << "\t" << newcount << std::endl;
    }
}


/** main - sample reads from stdin (pat format). Output them to stdout */
int main( int argc, char *argv[])
{
    if (argc != 2){
        std::cerr << "Usage: pat_sampler RATE" << std::endl;
        return -1;
    }

    try {
        double rate = std::stof(argv[1]);

        // Sample from stdin
        for (std::string line_str; std::getline(std::cin, line_str);) {
            if (line_str.empty()) { continue; } // skip empty lines
            std::vector<std::string> tokens = line2tokens(line_str);
            consider_line(tokens, rate);
        }
    }
    catch(std::exception &e) {
        std::cerr << "failed sampling from pat" << std::endl;
        std::cerr << e.what() << std::endl;
        return -1;
    }
}
