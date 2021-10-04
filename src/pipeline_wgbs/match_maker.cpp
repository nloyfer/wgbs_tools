#include <iostream>
#include <ctime>
#include <vector>
#include <algorithm>    // std::sort
#include <fstream>
#include <regex>
#include <sstream>


std::vector<std::string> line2tokens(const std::string &line) {
    /** Break string line to words (a std::vector of string tokens) */
    std::vector<std::string> result;
    std::string cell;
    std::stringstream lineStream(line);
    while(getline(lineStream, cell, '\t'))
        result.push_back(cell);
    if (result.empty()) { throw std::runtime_error("line2tokens: tokens shouldn't be empty!"); }
    return result;
}

void print_vec(std::vector<std::string> &vec){
    /** print a vector to stderr, tab separated */
    for (auto &j: vec)
        std::cerr << j << std::endl;
}

std::string addCommas(int num) {
    /** add commas to an integer, and return it as a string */
    std::string s = std::to_string(num);
    int n = s.length() - 3;
    while (n > 0) {
        s.insert(n, ",");
        n -= 3;
    }
    return s;
}

int read_pos(std::string &read) { return std::stoi(line2tokens(read)[3]); }
int read_mate_pos(std::string &read) { return std::stoi(line2tokens(read)[7]); }

struct PairedEnd
{
    /**  This struct holds paired reads (read1, read2) 
     *   with all their fields from the BAM file.
     *   and a key, which is the start position of the read that starts first.
     */
    int key;    // starting position of the read that starts first.
    std::string read1, read2; // paired reads.

    PairedEnd(const std::string& r1, const std::string& r2): read1(r1), read2(r2) {
        std::vector<std::string> tokens = line2tokens(read1);
        key = std::stoi(tokens[3]);
        if (read2.empty()) { return; }
        
        // swap reads order if needed
        int k2 = std::stoi(tokens[7]);
        if ((k2 < key)) {
            key = k2;
            read1 = r2;
            read2 = r1;
        }
    }
};

void resortByPos(std::vector<PairedEnd> &data, std::ostream &myfile) {
    std::sort(data.begin(), data.end(),
              [](const PairedEnd& lhs, const PairedEnd& rhs) { return lhs.key < rhs.key; });
    for (auto &j: data) {
        myfile << j.read1 << std::endl;
        if (!(j.read2.empty())) myfile << j.read2 << std::endl;
    }
}

std::vector<std::string> flush_data(std::vector<std::string> &data,
                                    std::ostream &myfile, 
                                    bool last_chunk,
                                    bool output_singles) {
    /* sort lines in data by read name. The pairs should then be in adjacent lines.
     * go over the lines one by one, and flush / stream forward adjacent pairs.
     * return the single lines - the ones with no adjacent match. They'll probably find
     * their match on the next data vector */

    if (data.empty()) return data;

    std::sort (data.begin(), data.end());
    std::string dummy = "";
    std::vector<bool> flushed(data.size()); // mark which lines are flushed. The others will be returned as singles.
    std::vector<PairedEnd> pairs_vec;
    int last_line_pos = read_pos(data.at(data.size() - 1));
    for (unsigned int i = 0; i < data.size() - 1; i++) {
        std::string r1 = data.at(i);
        std::string r2 = data.at(i + 1); 

        // if r1 and r2 have the same name - they are a pair
        if (r1.substr(0, r1.find('\t')) == r2.substr(0, r2.find('\t'))) {
            pairs_vec.push_back(PairedEnd(r1, r2));
            flushed[i] = true;
            flushed[i + 1] = true;
            i++;
        }
        else {
            if (output_singles) {
                if ((read_mate_pos(r1) < last_line_pos) ||
                        (last_chunk)) {    // read i has no hope finding a mate
                    pairs_vec.push_back(PairedEnd(r1, dummy));
                    flushed[i] = true;
                } else {  // else, r1 is tagged as single
                    flushed[i] = false;
                }
            } else {
                std::cerr << "skipping " << line2tokens(r1)[0] << std::endl;
                flushed[i] = false;
            }
        }
    } 

    std::vector<std::string> optimistics;
    for (unsigned int i = 0; i < data.size(); i++) {
        if (!(flushed[i])) {
            if (output_singles) {
                pairs_vec.push_back(PairedEnd(data.at(i), dummy));
            } else {
                optimistics.push_back(data.at(i));
            }
        }
    }
    resortByPos(pairs_vec, myfile);

    data.clear();
    return optimistics;
}

void filter_singles(std::vector<std::string> &data, 
                    std::vector<std::string> &singles,
                    std::string last_line) {
    /* move from data to singles all lines that definitely don't have a mate in the input bam
     * let r be a read. assuming the bam is sorted by position, 
     * if the last seen line starts in a position larger than the position of r's expected mate,
     * then r's mate is undoubtedly missing. r is moved the the singles vector */
    if (data.empty()) return;

    std::vector<std::string> optimistics;
    int cur_pos = read_pos(last_line);
    for (unsigned int i = 0; i < data.size(); i++) {
        if (read_mate_pos(data[i]) < cur_pos)    // read i has no hope finding a mate
            singles.push_back(data[i]);
        else
            optimistics.push_back(data[i]);
    }
    data = optimistics;
}


void action(bool output_singles) {
    std::vector<std::string> data;      // a chunck of 2,000 input lines

    std::ostream &outfile(std::cout);

    int line_i = 0;
    //clock_t begin = clock();
    std::string log_pref = "[match maker] ";
//    std::ios_base::sync_with_stdio(false);  // improves reading speed by x70
    for (std::string line; std::getline(std::cin, line) && (line_i > -1); line_i++) {
        
        // add chromosome to log prefix
        if (line_i == 0) { log_pref += "[ " + line2tokens(line)[2] + " ] "; }


        data.push_back(line);
        //if (line_i && (line_i % 10000 == 0)) {
        if (line_i && (line_i % 2000 == 0)) {

            // flush data only if positions of first and last reads are far apart
            // This is rarely an issue. But sometimes in extremely high mappability 
            // regions (e.g. mm9 chr14:100734330), the data buffer could grow to 4M reads
            if (read_pos(line) - read_pos(data.at(0)) > 160) {
                data = flush_data(data, outfile, false, output_singles);
            }
        }
    }
    data = flush_data(data, outfile, true, output_singles);

    if (line_i == 0) {
        std::cerr << log_pref << "no reads found " << std::endl;
        return;
    }
    std::cerr << log_pref << "finished " << addCommas(line_i) << " lines." << std::endl;
    if (!(data.empty())) {
        std::cerr << log_pref << "Filtered " << addCommas(data.size()) << " unpaired reads" << std::endl;
        //std::cerr << data[0];
    }
}


int main(int argc, char **argv) {
    bool output_singles = (argc > 1) && (std::string(argv[1]) == "-s");
    // if output_singles is not set, output only pairs of reads, 
    // and ignore reads whose mate is absent.
    action(output_singles);
    return 0;
}
