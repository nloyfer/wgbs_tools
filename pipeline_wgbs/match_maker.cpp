#include <iostream>
#include <ctime>
#include <vector>
#include <algorithm>    // std::sort
#include <fstream>
#include <regex>
#include <sstream>


/*
 * Interesting flags:
 * read1, read2 = (99, 147)
 * read1, read2 = (83, 163)
 */

// awk '{if ($2==99 || $2==147 || $2==163 || $2==83) print}'

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

struct PairedEnd
{
    int key;
    std::string read1, read2;

    PairedEnd(const std::string& r1, const std::string& r2) {
        std::vector<std::string> tokens = line2tokens(r1);
        int k1 = stoi(tokens[3]);
        int k2 = stoi(tokens[7]);
        if (k2 > k1) {
            read1 = r1;
            read2 = r2;
            key = k1;
        } else {
            read1 = r2;
            read2 = r1;
            key = k2;
        }
    }
};

void resortByPos(std::vector<PairedEnd> &data, std::ostream &myfile) {
    std::sort(data.begin(), data.end(),
              [](const PairedEnd& lhs, const PairedEnd& rhs) { return lhs.key < rhs.key; });
//    for (int i = 0; i < data.size(); i++) {
    for (auto &j: data) {
        myfile << j.read1 << std::endl;
        myfile << j.read2 << std::endl;
    }
}


std::vector<std::string> flush_data(std::vector<std::string> &data, std::ostream &myfile) {
    /* sort lines in data by read name. The pairs should then be in adjacent lines.
     * go over the lines one by one, and flush / stream forward adjacent pairs.
     * return the single lines - the ones with no adjacent match. They'll probably find
     * their match on the next data vector */

    if (data.empty())
        return data;

    std::sort (data.begin(), data.end());
    std::vector<bool> flushed(data.size()); // mark which lines are flushed. The others will be returned as singles.
    std::vector<PairedEnd> pairs_vec;
    for (unsigned int i = 0; i < data.size() - 1; i++) {
        std::string r1 = data.at(i);
        std::string r2 = data.at(i + 1);

        // if r1 and r2 have the same name - they are a pair
        if (r1.substr(0, r1.find('\t')) == r2.substr(0, r2.find('\t'))) {
//            myfile << r1 << std::endl;
//            myfile << r2 << std::endl;
            pairs_vec.push_back(PairedEnd(r1, r2));
            flushed[i] = true;
            flushed[i + 1] = true;
            i++;
        }
        else {  // else, r1 is tagged as single
            flushed[i] = false;
        }
    }

    std::vector<std::string> singles;
    for (unsigned int i = 0; i < data.size(); i++) {
        if (!(flushed[i])) {
            singles.push_back(data.at(i));
        }
    }

    resortByPos(pairs_vec, myfile);

    data.clear();
    return singles;
}


void filter_singles(std::vector<std::string> &data, std::vector<std::string> &singles, std::string last_line) {
    /* move from data to singles all lines that definitely don't have a mate in the input bam
     * let r be a read. assuming the bam is sorted by position, if the last seen line starts in a position larger than
     * the position of r's expected mate, than r's mate is undoubtedly missing. r is moved the the singles vector */
    if (data.empty())
        return;

    std::vector<std::string> optimistics;
    int cur_pos = stoi(line2tokens(last_line)[3]);
    for (unsigned int i = 0; i < data.size(); i++) {
        if (stoi(line2tokens(data[i])[7]) < cur_pos)    // read i has no hope finding a mate
            singles.push_back(data[i]);
        else
            optimistics.push_back(data[i]);
    }
    data = optimistics;
}

std::string addCommas(int num) {
    auto s = std::to_string(num);
    int n = s.length() - 3;
    while (n > 0) {
        s.insert(n, ",");
        n -= 3;
    }
    return s;
}

void action() {
    std::vector<std::string> data;      // a chunck of 2,000 input lines
    std::vector<std::string> singles;   // lines with no current match

    std::ostream &outfile(std::cout);

    int line_i = 0;
    std::string log_pref = "[match maker] ";
//    std::ios_base::sync_with_stdio(false);  // improves reading speed by x70
    for (std::string line; std::getline(std::cin, line) && (line_i > -1); line_i++) {
        if (line_i == 0) {
             log_pref += "[ " + line2tokens(line)[2] + " ] ";
        }
        data.push_back(line);
        if (line_i && (line_i % 10000 == 0)) {

            // flush data only if positions of first and last reads are far apart
            // This is rarely an issue. But sometimes in extremely high mappability 
            // regions (e.g. mm9 chr14:100734330), the data buffer could grow to 4M reads
            if (stoi(line2tokens(line)[3]) - stoi(line2tokens(data[0])[3]) > 160) {
                data = flush_data(data, outfile);
                filter_singles(data, singles, line);
            }

            //if (line_i % 1000000 == 0) {
                //std::cerr  << log_pref << " line " << addCommas(line_i) << std::endl;
            //}
        }
    }
    data = flush_data(data, outfile);

    singles.insert( singles.end(), data.begin(), data.end() );  // concatenate singles and data
//    std::cerr << "singles:" << std::endl;
//    print_vec(singles);
    std::cerr << log_pref << "Number of unpaired reads: " << addCommas(singles.size()) << std::endl;
    if (line_i > 0) {
        std::cerr << log_pref << "finished " << addCommas(line_i) << " lines." << std::endl;
    } else {
        std::cerr << log_pref << "no pairs found " << std::endl;
    }
//    outfile.close();
}


int main() {
    clock_t begin = clock();

    action();

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //std::cerr << elapsed_secs << std::endl;

    return 0;
}
