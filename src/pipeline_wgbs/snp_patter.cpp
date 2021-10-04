//
// Created by nloyfer on 9/14/18.
//

#include "snp_patter.h"
#include <set>

char METH = 'C';
char UNMETH = 'T';
char UNKNOWN = '.';
std::string TAB = "\t";

/***************************************************************
 *                                                             *
 *               Print methods                                 *
 *                                                             *
 ***************************************************************/

std::vector <std::string> line2tokens(std::string &line) {
    /** Break string line to words (a vector of string tokens) */
    std::vector <std::string> result;
    std::string cell;
    std::stringstream lineStream(line);
    while (getline(lineStream, cell, '\t')) {
        result.push_back(cell);
    }
    return result;
}

void print_vec(std::vector <std::string> &vec) {
    /** print a vector to stderr, tab separated */
    std::string sep = "";
    for (auto &j: vec) {
        std::cerr << sep << j;
        sep = TAB;
    }
    std::cerr << std::endl;
}


std::string addCommas(int num) {
    /** convert integer to string with commas */
    auto s = std::to_string(num);
    int n = s.length() - 3;
    while (n > 0) {
        s.insert(n, ",");
        n -= 3;
    }
    return s;
}

/***************************************************************
 *                                                             *
 *               Load FASTA                                    *
 *                                                             *
 ***************************************************************/

std::string exec(const char* cmd) {
    /** Execute a command and load output to string */
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("[ patter ] popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

/***************************************************************
 *                                                             *
 *               Niche                                         *
 *                                                             *
 ***************************************************************/

bool validate_seq_bp(std::string &seq, std::string &ref, bool reversed, int margin) {
    /** if blueprint filter is set, check for bisulfit conversionrate.
     * Return False if this rate is too low, or if there are not enough CH's
     */
    int nr_conv = 0;
    int nr_non_conv = 0;
    // case 1 - read 1 (not reversed)
    if (! reversed) {
        for (unsigned long j = 0; j < ref.length() - 1; j++) {
            // skip margins
            if ((j < margin) || (j >= seq.size() - margin)) { continue; }
            if ((ref[j] == 'C') && (ref[j + 1] != 'G')) {
                if (seq[j] == 'C') {
                    nr_non_conv++;
                } else if (seq[j] == 'T') {
                    nr_conv++;
                }
            }
        }
    } else { // case 2 - read 2 (reversed)
        for (unsigned long j = 1; j < ref.length(); j++) {
            // skip margins
            if ((j < margin) || (j >= seq.size() - margin)) { continue; }
            if ((ref[j] == 'G') && (ref[j + 1] != 'C')) {
                if (seq[j] == 'G') {
                    nr_non_conv++;
                } else if (seq[j] == 'A') {
                    nr_conv++;
                }
            }
        }
    }
    int nr_ch = nr_conv + nr_non_conv;
    if (nr_ch < 3) {
        return false;
    }
    return (((float) nr_conv / (float) nr_ch) >= 0.9) ? true : false;
}


/***************************************************************
 *                                                             *
 *             Process single / pair of reads                  *
 *                                                             *
 ***************************************************************/

std::string snp_patter::clean_CIGAR(std::string seq, std::string CIGAR) {

    /** use CIGAR string to adjust 'seq' so it will be comparable to the reference.
     * e.g, remove false letters ('I'), insert fictive letters ('D') etc. */

    // parse CIGAR and convert it to a couple of vectors: chars, nums.
    // e.g, '2S9M' will become ['S', 'M'] and [2, 9]include <set>
    std::vector<char> chars;
    std::vector<unsigned long> nums;
    std::string cur_num;
    for (auto c: CIGAR) {
        if (isdigit(c)) {
            cur_num.push_back(c);
        } else {
            nums.push_back(stoul(cur_num));
            cur_num.clear();
            chars.push_back(c);
        }
    }

    // build the adjusted seq, using original seq and the vectors:
    std::string adjusted_seq;
    for (int i = 0; i < (int) chars.size(); i++) {
        if (chars[i] == 'M') {
            adjusted_seq += seq.substr(0, nums[i]);
            seq = seq.substr(nums[i], seq.length() - nums[i]);
        } else if (chars[i] == 'D') {
            for (unsigned long j = 0; j < nums[i]; j++)
                adjusted_seq += 'N';
        } else if ((chars[i] == 'I') || (chars[i] == 'S')) {
            seq = seq.substr(nums[i], seq.length() - nums[i]);
        } else {
            throw std::invalid_argument("[ patter ] Unknown CIGAR character: " +
                                        std::string(1, chars[i]));
        }
    }

    return adjusted_seq;
}


int strip_pat(std::string &pat) {
    // remove dots from the tail (e.g. CCT.C.... -> CCT.C)
    pat = pat.substr(0, pat.find_last_not_of(UNKNOWN) + 1);
    if (pat == "") { return -1; }
    // remove dots from the head (..CCT -> CCT)
    int pos = pat.find_first_not_of(UNKNOWN);
    if (pos > 0) {
        pat = pat.substr(pos, pat.length() - pos);
    }
    return pos;
}

char snp_patter::compareSeqToRef(std::string &seq,
                                std::string &ref,
                                bool reversed,
                                std::string &meth_pattern,
                                int start_pos) {
    /** compare seq string to ref string. generate the methylation pattern, and return
     * the CpG index of the first CpG site in the seq (or -1 if there is none) */

    char snp_letter = 'Z';
    if (start_pos == 57418067){
        int x = 0;
    }
    int snp_index = snp_pos - start_pos;
    if (snp_index >= seq.length()) {
        return snp_letter;
    }
    char snp_val = seq[snp_index];
    if (((snp_let1 == 'C' && snp_let2 == 'T') || (snp_let2 == 'C' && snp_let1 == 'T')) && !reversed){
        return snp_letter;
    }
    if (((snp_let1 == 'G' && snp_let2 == 'A') || (snp_let2 == 'G' && snp_let1 == 'A')) && reversed){
        return snp_letter;
    }
    std::set<char> allowed_letters1;
    if (snp_let1 == 'C' && snp_let2 != 'T') {
        allowed_letters1 = {'C', 'T'};
    } else if (snp_let1 == 'G' && snp_let2 != 'A') {
        allowed_letters1 = {'G', 'A'};
    } else {
        allowed_letters1 = {snp_let1};
    }
    std::set<char> allowed_letters2;
    if (snp_let2 == 'C' && snp_let1 != 'T') {
        allowed_letters2 = {'C', 'T'};
    } else if (snp_let2 == 'G' && snp_let1 != 'A') {
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



std::vector <std::string> merge(std::vector <std::string> l1, std::vector <std::string> l2) {
    /** Merge 2 complementary lines to a single output.
     * each line has the following fields: [chr, startCpG, snp_let, pat]
     * One or more of the lines may be empty */

    // if one of the lines is empty - return the other
    if (l1.empty()) { return l2; }
    if (l2.empty()) { return l1; }

    // Swap lines s.t l1 starts before l2
    if (stoi(l1[1]) > stoi(l2[1])) {
        std::vector <std::string> tmp = l1;
        l1 = l2;
        l2 = tmp;
    }

    int start1 = stoi(l1[1]), start2 = stoi(l2[1]);
    std::string pat1 = l1[3], pat2 = l2[3];

    std::string merged_pat;  // output pattern
    int last_site = std::max(start1 + pat1.length(), start2 + pat2.length()); // location of last CpG from both reads

    if (last_site - start1 > MAX_PAT_LEN) // sanity check: make sure the two reads are not too far apart
        throw std::invalid_argument("invalid pairing. merged read is too long");

    // init merged_pat with missing values
    for (int i = start1; i < last_site; i++)
        merged_pat += ".";

    // set merged_pat head with pat1
    for (unsigned long i = 0; i < pat1.length(); i++)
        merged_pat[i] = pat1[i];

    // set pat2 in the adjusted position
    for (unsigned long i = 0; i < pat2.length(); i++) {
        int adj_i = i + start2 - start1;
        if (merged_pat[adj_i] == UNKNOWN) {   // this site was missing from read1
            merged_pat[adj_i] = pat2[i];
        } else if ((pat2[i] != UNKNOWN) && (merged_pat[adj_i] != pat2[i])) {
            // read1 and read2 disagree, and none of them is missing ('.').
            // treat this case as a missing value for now
            // future work: consider only the read with the higher quality.
            merged_pat[adj_i] = UNKNOWN;
        }
    }
    // strip merged pat:
    int pos = strip_pat(merged_pat);
    if (pos < 0 ) { return {}; }
    l1[1] = std::to_string(start1 + pos);
//    l1[]
    l1[3] = merged_pat;
    return l1;
}


char snp_patter::samLineToPatVec(std::vector <std::string> tokens) {
    /** Given tokens of a sam line:
     *       QNAME, FLAG, RNAME (chrom), POS, MAPQ, CIGAR, RNEXT,
     *       PNEXT, TLEN, SEQ, QUAL, and possibly more.
     *  return a new vector with the following fields:
     *  [chr, first_CpG_ind, meth_pattern, start_loc, length]
     *
     *  In case line is empty or invalid, return an empty vector,
     *  and update the corresponding counter
     *  */

    std::vector <std::string> res;
    if (tokens.empty()) {
        return 'Z';
    }
    try {

        if (tokens.size() < 11) {
            throw std::invalid_argument("too few arguments in line");
        }
        unsigned long start_locus = stoul(tokens[3]);   // Fourth field from bam file
        int samflag = stoi(tokens[1]);
        std::string seq = tokens[9];
        //std::string bp_qual = tokens[10];
        std::string CIGAR = tokens[5];

        seq = clean_CIGAR(seq, CIGAR);

//        unsigned long seq_len = seq.length();   // We may need the original length later.
        std::string ref = "";//genome_ref.substr(start_locus - 1, seq_len);

        // build methylation pattern:
        std::string meth_pattern;
        bool reversed;
        if (is_paired_end) {
            reversed = ( ((samflag & 0x53) == 83) || ((samflag & 0xA3) == 163) );
        } else {
            reversed = ((samflag & 0x0010) == 16);
        }
        if ((samflag & 0x400) == 1024){
            return 'Z';
        }

        return compareSeqToRef(seq, ref, reversed, meth_pattern, start_locus);
    }
    catch (std::exception &e) {
        std::string msg = "[ " + chr + " ] " + "Exception while processing line "
                          + std::to_string(line_i) + ". Line content: \n";
        std::cerr << "[ patter ] " << msg;
        print_vec(tokens);
        std::cerr << "[ patter ] " << e.what() << std::endl;
        readsStats.nr_invalid++;
    }
//    res.clear();
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
        char snp_read1 = samLineToPatVec(tokens1);
        char snp_read2 = samLineToPatVec(tokens2);
        if (snp_read1 == 'Z') {
            readsStats.nr_empty++;
            if (snp_read2 == 'Z'){
                readsStats.nr_empty++;
            }
            return snp_read2;
        }
        if (snp_read2 == 'Z') {
            readsStats.nr_empty++;
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
    if (min_cpg > 1) {
        msg += addCommas(readsStats.nr_short) + " with too few CpGs. ";
    }
    msg += addCommas(readsStats.nr_invalid) + " invalid. ";
    if (blueprint) {
        msg += addCommas(readsStats.nr_bad_conv) + " bad conversion. ";
    }
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
    if (genome_ref.empty()) {
        is_paired_end = first_line(line_str);
//        load_genome_ref();
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

void snp_patter::proc_sam_in_stream(std::istream& in){
    /** parse stdin for sam format lines (single- or pired-end).
     * Translate them to pat format, and output to stdout */

    bool first_in_pair = true;
//    std::ios_base::sync_with_stdio(false);
    std::vector <std::string> tokens1, tokens2;
    std::string read1, read2;
    for (std::string line_str; std::getline(in, line_str); line_i++){

        print_progress();

        // skip empty lines
        if (line_str.empty())
            continue;

        initialize_patter(line_str);

        // paired-end file, and current row is first out of a couple of rows
        if (first_in_pair && is_paired_end) {
            read1 = line_str;
            tokens1 = line2tokens(line_str);
            first_in_pair = false;
            readsStats.nr_pairs++;
            continue;
        }// otherwise (second row in couple, or not-paired-end), line is processed:

        tokens2 = line2tokens(line_str);

        if ((!(tokens2.empty())) && (!(tokens1.empty()))
            && (tokens1[0] != tokens2[0])) {
            // process couple of lines. write to stdout
            // in case of single-end input file, tokens1 will be empty
            // in case of non-paired read, use emtpy tokens for second param
            read2 = line_str;
            std::vector <std::string> empty_tokens;
            char snp_let = proc2lines(tokens1, empty_tokens);
            if (snp_let == snp_let1){
                std::cout << read1 << "\n";
            }
            tokens1 = tokens2;
            read1 = line_str;
        } else {
            first_in_pair = true;   // next line will be first of couple.
            read2 = line_str;
            char snp_let = proc2lines(tokens1, tokens2);
            if (snp_let == snp_let1){
                std::cout << read1 << "\n" << read2 << "\n";
            }
        }



    }
    print_stats_msg();
}

void snp_patter::action(std::string samFilePath) {
    std::ifstream samFile(samFilePath, std::ios::in);

    if (!(samFile)){
        proc_sam_in_stream(std::cin);
    } else {
        if (samFile.is_open()) {
            proc_sam_in_stream(samFile);
            samFile.close();
        }
    }
}


/***************************************************************
 *                                                             *
 *                     Main                                    *
 *                                                             *
 ***************************************************************/

/*
 * Input Arguments Parsering class
 */
class InputParser{
public:
    InputParser (int &argc, char **argv){
        for (int i=1; i < argc; ++i)
            this->tokens.emplace_back(std::string(argv[i]));
    }
    const std::string& getCmdOption(const std::string &option) const{
        std::vector<std::string>::const_iterator itr;
        itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
        if (itr != this->tokens.end() && ++itr != this->tokens.end()){
            return *itr;
        }
        static const std::string empty_string;
        return empty_string;
    }
    bool cmdOptionExists(const std::string &option) const{
        return std::find(this->tokens.begin(), this->tokens.end(), option)
               != this->tokens.end();
    }
private:
    std::vector <std::string> tokens;
};

bool is_number(const std::string& s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}


int main(int argc, char **argv) {
    clock_t begin = clock();
    try {
        InputParser input(argc, argv);
        long snp_pos = -1;
        char snp_let1 = 'Z';
        char snp_let2 = 'Z';
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
        if (argc < 3) {
            throw std::invalid_argument("Usage: patter GENOME_PATH CPG_CHROM_SIZE_PATH [--bam] [--mbias MBIAS_PATH]");
        }
        snp_patter p(snp_pos, snp_let1, snp_let2);
        p.action("");

    }
    catch (std::exception &e) {
        std::cerr << "[ patter ] Failed! exception:" << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    }
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//    std::cout << elapsed_secs << std::endl;
//    std::cerr << elapsed_secs << std::endl;
    return 0;
}
