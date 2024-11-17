//
// Created by nloyfer on 9/14/18.
//

#include "snp_patter.h"
#include "patter_utils.h"
#include <set>


/***************************************************************
 *                                                             *
 *             Process single / pair of reads                  *
 *                                                             *
 ***************************************************************/

char snp_patter::compareSeqToRef(std::string &seq,
                                bool bottom,
                                std::string &qual_str,
                                int start_pos) {

    char snp_letter = 'Z';
    int snp_index = snp_pos - start_pos;
    if (snp_index >= seq.length()) {
        return snp_letter;
    }
    int qual = int(qual_str[snp_index]) - 33;
    if (qual < qual_filter){
        return 'Z';
    }
    char snp_val = seq[snp_index];
    if (((snp_let1 == 'C' && snp_let2 == 'T') || (snp_let2 == 'C' && snp_let1 == 'T')) && !bottom){
        return snp_letter;
    }
    if (((snp_let1 == 'G' && snp_let2 == 'A') || (snp_let2 == 'G' && snp_let1 == 'A')) && bottom){
        return snp_letter;
    }
    std::set<char> allowed_letters1;
    if (snp_let1 == 'C' && snp_let2 != 'T' && !bottom) {
        allowed_letters1 = {'C', 'T'};
    } else if (snp_let1 == 'G' && snp_let2 != 'A' && bottom) {
        allowed_letters1 = {'G', 'A'};
    } else {
        allowed_letters1 = {snp_let1};
    }
    std::set<char> allowed_letters2;
    if (snp_let2 == 'C' && snp_let1 != 'T' && !bottom) {
        allowed_letters2 = {'C', 'T'};
    } else if (snp_let2 == 'G' && snp_let1 != 'A' && bottom) {
        allowed_letters2 = {'G', 'A'};
    } else {
        allowed_letters2 = {snp_let2};
    }
    if (allowed_letters1.find(snp_val) != allowed_letters1.end()){
        snp_letter = snp_let1;
    } else if (allowed_letters2.find(snp_val) != allowed_letters1.end()) {
        snp_letter = snp_let2;
    }
    return snp_letter;
}

char snp_patter::samLineToSNP(std::vector <std::string> tokens) {
    /** Given tokens of a sam line:
     *       QNAME, FLAG, RNAME (chrom), POS, MAPQ, CIGAR, RNEXT,
     *       PNEXT, TLEN, SEQ, QUAL, and possibly more.
     *  return a new vector with the following fields:
     *  [chr, first_CpG_ind, meth_pattern, start_loc, length]
     *
     *  In case line is empty or invalid, return an empty vector,
     *  and update the corresponding counter
     *  */

    if (tokens.empty()) { return 'Z'; }
    try {

        if (tokens.size() < 11) {
            throw std::invalid_argument("too few arguments in line");
        }
        unsigned long start_locus = stoul(tokens[3]);
        int samflag = stoi(tokens[1]);
        std::string seq = tokens[9];
        std::string qual_str = tokens[10];
        std::string CIGAR = tokens[5];

        // skip duplicated reads
        if ((samflag & 0x400) == 1024){ return 'Z'; }

        seq = clean_CIGAR(seq, CIGAR);
        qual_str = clean_CIGAR(qual_str, CIGAR);
        bool bottom = is_bottom(samflag, is_paired_end);

        return compareSeqToRef(seq, bottom, qual_str, start_locus);
    }
    catch (std::exception &e) {
        std::string msg = "[ " + chr + " ] " + "Exception while processing line "
                          + std::to_string(line_i) + ". Line content: \n";
        std::cerr << "[ patter ] " << msg;
        print_vec(tokens);
        std::cerr << "[ patter ] " << e.what() << std::endl;
        readsStats.nr_invalid++;
    }
    return 'Z'; // return empty vector
}


char snp_patter::proc2lines(std::vector <std::string> tokens1,
                            std::vector <std::string> tokens2) {

    /** print result to stdout */
    try {
        // sanity check: two lines must have the same QNAME
        if ((!(tokens2.empty())) && (!(tokens1.empty()))
            && (tokens1[0] != tokens2[0])) {
            readsStats.nr_invalid += 2;
            throw std::invalid_argument("lines are not complements!");
        }

        // Merge 2 complementary lines to a single output.
        char snp_read1 = samLineToSNP(tokens1);
        char snp_read2 = samLineToSNP(tokens2);
        if (snp_read1 == 'Z') {
            readsStats.nr_empty++;
            return snp_read2;
        }
        if (snp_read2 == 'Z') {
            if (!(tokens2.empty())){
                readsStats.nr_empty++;
            }
            return snp_read1;
        }
        if (snp_read1 != snp_read2) {
            readsStats.nr_empty++;
            return 'Z';
        }
        return snp_read1;
    }
    catch (std::exception &e) {
        return 'Z';
    }
}


/***************************************************************
 *                                                             *
 *                     print stats                             *
 *                                                             *
 ***************************************************************/

void snp_patter::print_stats_msg() {
    /** print informative summary message */

    int sucess = line_i ? int((1.0 - ((double) readsStats.nr_invalid / line_i)) * 100.0) : 0;

    std::string msg = "[ " + chr + " ] ";
    msg += "finished " + addCommas(line_i) + " lines. ";
    if (is_paired_end) {
        msg += "(" + addCommas(readsStats.nr_pairs) + " pairs). ";
    }
    msg += addCommas(line_i - readsStats.nr_empty - readsStats.nr_invalid) + " good, ";
    msg += addCommas(readsStats.nr_empty) + " empty, ";
    msg += addCommas(readsStats.nr_invalid) + " invalid. ";
    msg += "(success " + std::to_string(sucess) + "%)\n";
    std::cerr << "[ patter ] " << msg;
}

/***************************************************************
 *                                                             *
 *                Init                                         *
 *                                                             *
 ***************************************************************/

bool snp_patter::first_line(std::string &line) {
    // find out if the input is paired- or single-end
    try {
        // find out if single- or paired-end
        std::vector <std::string> tokens = line2tokens(line);
        auto flag = (uint16_t) stoi(tokens.at(1));
        chr = tokens.at(2);
        return (flag & 1);     // file is paired end iff flag 0x1 is on
    }
    catch (std::exception &e) {
        std::cerr << "[ patter ] " << "[ " + chr + " ]" << "Invalid first line: \n" << line;
        std::cerr << "\nexception: " << e.what() << std::endl;
        throw e;
    }
}

void snp_patter::initialize_patter(std::string &line_str) {
    // first line based initializations
    if (chr == "") {
        is_paired_end = first_line(line_str);
    }
}

/***************************************************************
 *                                                             *
 *                     Parse bam                               *
 *                                                             *
 ***************************************************************/

void snp_patter::print_progress(){
    if (line_i && !(line_i % 5000000)){
        clock_t tock = clock();
        double elapsed_secs = double(tock - tick) / (CLOCKS_PER_SEC * 60);
        tick = tock;
        std::cerr << "[ patter ] [ " + chr + " ]" << " line " << addCommas(line_i)
                  << " in " << std::setprecision(2) << elapsed_secs << " minutes." << std::endl;
    }
}

char snp_patter::proc1line(std::vector <std::string> &tokens1) {
    return proc2lines(tokens1, dummy_tokens);
}

void snp_patter::proc_sam_in_stream(std::istream& in){
    /** parse stdin for sam format lines (single- or pired-end).
     * Translate them to pat format, and output to stdout */

    bool first_in_pair = true;
    std::vector <std::string> tokens1, tokens2;
    std::string read1, read2;
    for (std::string line_str; std::getline(in, line_str); line_i++){

        print_progress();

        // skip empty lines
        if (line_str.empty()) { continue; }

        initialize_patter(line_str);

        if (tokens1.empty()) {
            tokens1 = line2tokens(line_str);
            read1 = line_str;
            if (!is_paired_end) {
                char snp_let = proc1line(tokens1);

                tokens1.clear();
                if (snp_let == snp_let1){
                    std::cout << read1 << "\n";
                }
            }
            continue;
        }
        tokens2 = line2tokens(line_str);
        read2 = line_str;
        if (are_paired(tokens1, tokens2)) {
            char snp_let = proc2lines(tokens1, tokens2);
            if (snp_let == snp_let1){
                std::cout << read1 << "\n" << read2 << "\n";
                readsStats.nr_pairs++;
            }
            tokens1.clear();
        } else {
            char snp_let = proc1line(tokens1);
            if (snp_let == snp_let1){
                std::cout << read1 << "\n";
            }
            tokens1 = tokens2;
            read1 = read2;
        }
    }
    if (! tokens1.empty()) {
        char snp_let = proc1line(tokens1);
        if (snp_let == snp_let1){
            std::cout << read1 << "\n";
        }
        tokens1 = tokens2;
        read1 = read2;
    }

    print_stats_msg();
}

void snp_patter::action(std::string samFilePath) {
    std::ifstream samFile(samFilePath, std::ios::in);

    if (!(samFile)){
        proc_sam_in_stream(std::cin);
    } else if (samFile.is_open()) {
            proc_sam_in_stream(samFile);
            samFile.close();
    }
}


/***************************************************************
 *                                                             *
 *                     Main                                    *
 *                                                             *
 ***************************************************************/

int main(int argc, char **argv) {
    try {
        InputParser input(argc, argv);
        long snp_pos = -1;
        char snp_let1 = 'Z';
        char snp_let2 = 'Z';
        int quality_filter = 0;
        if (input.cmdOptionExists("--snp_pos")){
            std::string snp_pos_string = input.getCmdOption("--snp_pos");

            if ( !is_number(snp_pos_string) )
                throw std::invalid_argument("invalid snp_pos argument. snp_pos should be a non-negative integer.");
            snp_pos = std::stol(snp_pos_string);
        } else {
            throw std::invalid_argument("must provide snp position.");
        }
        if (input.cmdOptionExists("--snp_let1")){
            std::string min_cpg_string = input.getCmdOption("--snp_let1");

            snp_let1 = min_cpg_string[0];
        } else {
            throw std::invalid_argument("must provide snp letter one.");
        }
        if (input.cmdOptionExists("--snp_let2")){
            std::string min_cpg_string = input.getCmdOption("--snp_let2");

            snp_let2 = min_cpg_string[0];
        } else {
            throw std::invalid_argument("must provide snp letter two.");
        }
        if (input.cmdOptionExists("--qual_filter")){
            std::string qual_str = input.getCmdOption("--qual_filter");

            if ( !is_number(qual_str) )
                throw std::invalid_argument("invalid qual_filter argument. qual_filter should be a non-negative integer.");
            quality_filter = std::stoi(qual_str);
        }
        if (argc < 3) {
            throw std::invalid_argument("Usage: patter GENOME_PATH CPG_CHROM_SIZE_PATH [--bam]");
        }
        snp_patter p(snp_pos, snp_let1, snp_let2, quality_filter);
        p.action("");

    }
    catch (std::exception &e) {
        std::cerr << "[ patter ] Failed! exception:" << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    }
    return 0;
}
