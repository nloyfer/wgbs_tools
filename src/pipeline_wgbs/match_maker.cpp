#include <iostream>
#include <ctime>
#include <vector>
#include <algorithm>    // std::sort
#include <fstream>
#include <regex>
#include <sstream>
#include "patter_utils.h"


int read_pos(std::string &read) { return std::stoi(line2tokens(read)[3]); }
int read_mate_pos(std::string &read) { return std::stoi(line2tokens(read)[7]); }
std::string read_chrom(std::string &read) { return line2tokens(read)[2]; }

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

    std::string last_chrom = read_chrom(data.at(data.size() - 1)); // get last chrom before sort
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
                flushed[i] = false;
            }
        }
    } 

    // Go over the non flushed reads.
    // If it's last chunk and we output singles, flush them
    // otherwise, push them to optimistics for the next round (or to oblivion)
    std::vector<std::string> optimistics;
    for (unsigned int i = 0; i < data.size(); i++) {
        if (!(flushed[i])) {
            if (last_chunk && output_singles) {
                pairs_vec.push_back(PairedEnd(data.at(i), dummy));
            } else {
                if (read_chrom(data.at(i)) == last_chrom){
                    optimistics.push_back(data.at(i));
                }
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

    long long int line_i = 0;
    std::string log_pref = "[match maker] ";
    bool is_chrom_set = false;
    for (std::string line; std::getline(std::cin, line) && (line_i > -1); line_i++) {

        // skip header lines
        if ((!is_chrom_set) && (!line.empty()) && (line[0] == '@')){
            std::cout << line << std::endl;
            continue;
        }
        
        // add chromosome to log prefix
        if (!is_chrom_set) { 
            log_pref += "[ " + line2tokens(line)[2] + " ] "; 
            is_chrom_set = true;
        }


        data.push_back(line);
        // If we moved to the next chromosome, flush data
        std::string chrom = read_chrom(line);
        if (line_i && (chrom != read_chrom(data.at(0)))){
            data = flush_data(data, outfile, true, output_singles);
        }
        if (line_i && (line_i % 50000 == 0)) {
        //if (line_i && (line_i % 2000 == 0)) {

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
    }
}


int main(int argc, char **argv) {
    bool output_singles = (argc < 2) || (std::string(argv[1]) != "--drop_singles");
    // if output_singles is not set, output only pairs of reads, 
    // and ignore reads whose mate is absent.
    action(output_singles);
    return 0;
}
