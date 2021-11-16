//
// Created by jrosensk on 14/11/2021.
//

#include "add_cpg_counts.h"
//#include <string_view>

char METH = 'C';
char UNMETH = 'T';
char UNKNOWN = '.';
std::string TAB = "\t";


//bool DEBUG = true;
bool DEBUG = false;

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

/*
 * Interesting flags (paird-end):
 * read1, read2 = (99, 147)
 * read1, read2 = (83, 163)
 *
 * for single end:
 * 0 or 16
 * 2 or 18
 */

std::vector <std::string> line2tokens(std::string &line) {
    /** Break string line to words (a vector of string tokens) */
    std::vector <std::string> result;
    std::string cell;
    std::stringstream lineStream(line);
    while (getline(lineStream, cell, '\t')) {
        result.push_back(cell);
    }
//    if (result.empty()) { throw std::runtime_error("line2tokens: tokens shouldn't be empty!"); }
    return result;
}

void print_vec(std::vector <std::string> &vec) {
    /** print a vector to stderr, tab separated */
    for (auto &j: vec)
        std::cerr << j << TAB;
    std::cerr << std::endl;
}

void vec2string(std::vector <std::string> &vec) {
    /** print a 5 items vector to stdout, tab separated */
    if (vec.size() == 5)  // vec length must be either 5 or 0.
        std::cout << vec[0] + TAB + vec[1] + TAB + vec[2] + TAB + vec[3] + TAB + vec[4] << std::endl;
}


int patter::find_cpg_inds_offset() {
    /**
     * Find the CpG-Index offset of the current chromosome.
     * i.e., How many CpGs there are before the first CpG in the current chromosome
     */

    // Open file CpG.chrome.size
    std::ifstream index_file(chrom_sz_path, std::ios::in);
    if (!(index_file)) {
        throw std::invalid_argument(" Error: Unable to read chrome_sizes path: " + chrom_sz_path);
    }

    // parse chrome size file, accumulate offset
    int offset = 0;
    for (std::string line_str; std::getline(index_file, line_str);) {
        std::vector <std::string> tokens = line2tokens(line_str);
        if (tokens[0] == chr) {
            return offset;
        } else {
            offset += std::stoi(tokens[1]);
        }
    }
    throw std::invalid_argument("Invalid chromosome " + chr);

}

std::vector<long> patter::fasta_index() {
    /**
     * Read the genome.fa.fai file
     * It has 5 columns: NAME, LENGTH, OFFSET, LINEBASES, LINEWIDTH, QUALOFFSET
     * Find the line corresponding to 'chr', and return LENGTH and OFFSET.
     */
    std::ifstream index_file(ref_path + ".fai", std::ios::in);
    if (!(index_file)) {
        throw std::invalid_argument(" Error: Unable to read index (fai) path: " + ref_path + ".fai");
    }
    std::vector<long> fai_numbers;
    // parse index file, load the right chromosome
    for (std::string line_str; std::getline(index_file, line_str);) {
        std::vector <std::string> tokens = line2tokens(line_str);
        if (tokens[0] == chr) {
            fai_numbers.push_back(std::stol(tokens[1]));
            fai_numbers.push_back(std::stol(tokens[2]));
            return fai_numbers;
        }
    }
    throw std::invalid_argument(" Error: chromosome not found in fai file: " + chr);
}

void patter::load_genome_ref() {
    /** Load the genome reference, corresponding to chromosome */

    // Open genome.fa file:
    std::ifstream ref_file(ref_path, std::ios::in);
    if (!(ref_file)) {
        throw std::invalid_argument(" Error: Unable to read reference path: " + ref_path);
    }

    // Read fai:
    std::vector<long> fai_numbers = fasta_index();
    // Jump to the beginning of the chromosome:
    ref_file.seekg(fai_numbers[1]);

    // parse *.fa file, load the lines (of 'chr') to one long string
    for (std::string line_str; std::getline(ref_file, line_str);) {
        if (line_str[0] == '>') {
            break;
        }
        for (auto &c: line_str) c = (char) toupper(c);  // change all characters to Uppercase
        genome_ref += line_str;
    }
    // Was the reference sequence loaded? (not empty)
    if ((long) genome_ref.length() != fai_numbers[0]) {
        std::cerr << chr << "'s length is " << genome_ref.length() << " != " << fai_numbers[0] << std::endl;
        throw std::invalid_argument(" Error: Reference genome's length is wrong. chromosome " + chr);
    }

    // parse genome to generate a CpG indexes map: <locus, CpG-Index>
    int cpg_ind_offset = find_cpg_inds_offset();
    for (int loc = 1, IlmnID = 1 + cpg_ind_offset; loc < (int) genome_ref.length(); loc++) {
        if ((genome_ref[loc] == 'G') && (genome_ref[loc - 1] == 'C'))
            dict.insert(std::make_pair(loc, IlmnID++));
    }
}


std::string patter::clean_seq(std::string seq, std::string CIGAR) {

    /** use CIGAR string to adjust 'seq' so it will be comparable to the reference.
     * e.g, remove false letters ('I'), insert fictive letters ('D') etc. */

    // parse CIGAR and convert it to a couple of vectors: chars, nums.
    // e.g, '2S9M' will become ['S', 'M'] and [2, 9]
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
            throw std::invalid_argument("Unknown CIGAR character: " +
                                        std::string(1, chars[i]));
        }
    }

    return adjusted_seq;
}


int patter::locus2CpGIndex(int locus) {
    /** translate genomic locus to CpG index (in range 1,...,28M~) */
    int start_site = 0;
    auto search = dict.find(locus);
    if (search != dict.end()) {
        start_site = search->second;
    } else {
        // This is an internal error - should never happen.
        throw std::logic_error("Internal Error. Unknown CpG locus: " + std::to_string(locus));
    }
    return start_site;
}



patter::MethylData compareSeqToRef2(std::string &seq,
                                    std::string &ref,
                                    bool reversed) {
    /** compare seq string to ref string. generate the methylation pattern, and return
     * the CpG index of the first CpG site in the seq (or -1 if there is none) */

    // ignore first/last 'margin' characters, since they are often biased
    size_t margin = 3;

    int countMethyl = 0, countUnmethyl = 0;
    patter::MethylData curRes;

    // find CpG indexes on reference sequence
    std::vector<int> cpg_inds;
    for (unsigned long j = 0; j < ref.length() - 1; j++) {
        if ((ref[j] == 'C') && (ref[j + 1] == 'G')) {
            cpg_inds.push_back(j);
        }
    }
    if (cpg_inds.empty()) {
        curRes.countUnmethyl = 0;
        curRes.countMethyl = 0;
        return curRes;
    }

    // generate the methylation pattern (e.g 'CC.TC'),
    // by comparing the given sequence to reference at the CpG indexes
    char REF_CHAR = reversed ? 'G' : 'C';
    char UNMETH_SEQ_CHAR = reversed ? 'A' : 'T';
    int shift = reversed ? 1 : 0;

    for (unsigned long j: cpg_inds) {
        j += shift;
        char s = seq[j];
        if ((j >= margin) && (j < seq.size() - margin)) {
            if (s == UNMETH_SEQ_CHAR) {
                countUnmethyl++;
            } else if (s == REF_CHAR) {
                countMethyl++;
            }
        }
    }
    curRes.countUnmethyl = countUnmethyl;
    curRes.countMethyl = countMethyl;

    return curRes;
}


int compareSeqToRef(std::string &seq,
                    std::string &ref,
                    bool reversed,
                    std::string &meth_pattern) {
    /** compare seq string to ref string. generate the methylation pattern, and return
     * the CpG index of the first CpG site in the seq (or -1 if there is none) */

    // ignore first/last 'margin' characters, since they are often biased
    size_t margin = 3;

    // find CpG indexes on reference sequence
    std::vector<int> cpg_inds;
    for (unsigned long j = 0; j < ref.length() - 1; j++) {
        if ((ref[j] == 'C') && (ref[j + 1] == 'G')) {
            cpg_inds.push_back(j);
        }
    }
    if (cpg_inds.empty()) {
        return -1;
    }

    // generate the methylation pattern (e.g 'CC.TC'),
    // by comparing the given sequence to reference at the CpG indexes
    char REF_CHAR = reversed ? 'G' : 'C';
    char UNMETH_SEQ_CHAR = reversed ? 'A' : 'T';
    int shift = reversed ? 1 : 0;

    char cur_status;
    for (unsigned long j: cpg_inds) {
        j += shift;
        char s = seq[j];
        cur_status = UNKNOWN;
        if ((j >= margin) && (j < seq.size() - margin)) {
            if (s == UNMETH_SEQ_CHAR)
                cur_status = UNMETH;
            else if (s == REF_CHAR)
                cur_status = METH;
        }
        meth_pattern.push_back(cur_status);
    }

    // an all-dots pattern is equivalent to empty pattern.
    if (meth_pattern.find_first_not_of('.') == std::string::npos)
        meth_pattern = "";


//    std::cerr << seq << std::endl;
//    std::cerr << ref << std::endl;
//    for (auto c: cpg_inds){
//        std::cerr << c << ", ";
//    }
//    std::cerr << std::endl;

    return cpg_inds[0];
}

std::vector <std::string> merge(std::vector <std::string> l1, std::vector <std::string> l2) {
    /** Merge 2 complementary lines to a single output.
     * each line has the following fields: [chr, startCpG, pat, start_bp, read_len_bp]
     * One or more of the lines may be empty */

    // Swap lines s.t l1 starts before l2
    if (stoi(l1[1]) > stoi(l2[1])) {
        std::vector <std::string> tmp = l1;
        l1 = l2;
        l2 = tmp;
    }

    int start1 = stoi(l1[1]), start2 = stoi(l2[1]);
    std::string pat1 = l1[2], pat2 = l2[2];

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
        if (merged_pat[adj_i] == '.') {   // this site was missing from read1
            merged_pat[adj_i] = pat2[i];
        } else if ((pat2[i] != '.') && (merged_pat[adj_i] != pat2[i])) {
            // read1 and read2 disagree. treat this case as a missing value for now
            // future work: consider only the read with the higher quality.
            merged_pat[adj_i] = '.';
        }
    }
    l1[1] = std::to_string(start1);
    l1[2] = merged_pat;
    return l1;
}

void patter::merge_and_print(std::vector <std::string> l1, std::vector <std::string> l2) {
    /** Merge 2 complementary lines to a single output.
     * each line has the following fields: [chr, startCpG, pat, start_bp, read_len_bp]
     * One or more of the lines may be empty */

    //  First line is empty
    if (l1.empty()) {
        return vec2string(l2);
    }
    // Second line is empty
    if (l2.empty()) {
        return vec2string(l1);
    }
    l1 = merge(l1, l2);
    vec2string(l1);
}

patter::MethylData meth_pattern_count(std::string meth_pattern) {
    patter::MethylData res;
    int countMethyl = 0;
    int countUnmethyl = 0;
    for (int i = 0; i < meth_pattern.length(); i++){
        char cur_char = meth_pattern[i];
        if(cur_char == METH){
            countMethyl++;
        } else if (cur_char == UNMETH){
            countUnmethyl++;
        }
    }
    res.countMethyl = countMethyl;
    res.countUnmethyl = countUnmethyl;
    return res;
}

patter::MethylData patter::merge_and_count_methyl_data(std::vector <std::string> l1, std::vector <std::string> l2) {
    /** Merge 2 complementary lines to a single output.
     * each line has the following fields: [chr, startCpG, pat, start_bp, read_len_bp]
     * One or more of the lines may be empty */


    if (l1.empty() && l2.empty()){
        patter::MethylData res;
        res.countMethyl = 0;
        res.countUnmethyl = 0;
        return res;
    }
    if (l1.empty()) {
        return meth_pattern_count(l2[2]);
    }
    if (l2.empty()) {
        return meth_pattern_count(l1[2]);
    }

    l1 = merge(l1, l2);
    patter::MethylData res = meth_pattern_count(l1[2]);
    //2 is the index of the meth_pattern
    return res;
}

std::string patter::samLineMethyldataMakeString(std::string originalLine, patter::MethylData md) {
    return '\t' + TAGNAMETYPE + std::to_string(md.countMethyl) + ',' + std::to_string(md.countUnmethyl) + '\n';
}

std::string patter::samLineToSamLineWithMethCounts(std::vector <std::string> tokens, std::string originalLine) {
    try {
        patter::MethylData curRow = samLineToMethCounts(tokens, originalLine);
        return samLineMethyldataMakeString(originalLine, curRow);
    }
    catch (std::exception &e) {
        std::string msg = "[ " + chr + " ] " + "Exception while processing line "
                          + std::to_string(line_i) + ". Line content: \n";
        std::cerr << msg;
        print_vec(tokens);
        std::cerr << e.what() << std::endl;
        readsStats.nr_invalid++;
    }
    return originalLine; //TODO flag for throw error or print original line
}

patter::MethylData patter::samLineToMethCounts(std::vector <std::string> tokens, std::string originalLine) {
    patter::MethylData errorres;
    if (tokens.empty()) {
        return errorres; //todo throw error
    }
    try {

        if (tokens.size() < 11) {
            throw std::invalid_argument("too few arguments in line");
        }
        unsigned long start_locus = stoul(tokens[3]);   // Fourth field from bam file
        int samflag = stoi(tokens[1]);
        std::string seq = tokens[9];
        std::string bp_qual = tokens[10];
        std::string CIGAR = tokens[5];

        seq = clean_seq(seq, CIGAR);

        unsigned long seq_len = seq.length();   // We may need the original length later.
        std::string ref = genome_ref.substr(start_locus - 1, seq_len);

        bool reversed;
        if (is_paired_end) {
            reversed = ((samflag == 83) || (samflag == 163));
        } else {
            reversed = ((samflag & 0x0010) == 16);
        }
        patter::MethylData curRow = compareSeqToRef2(seq, ref, reversed);

        return curRow;
    }
    catch (std::exception &e) {
        std::string msg = "[ " + chr + " ] " + "Exception while processing line "
                          + std::to_string(line_i) + ". Line content: \n";
        std::cerr << msg;
        print_vec(tokens);
        std::cerr << e.what() << std::endl;
        readsStats.nr_invalid++;
    }
    return errorres; //TODO throw error
}


std::vector <std::string> patter::samLineToPatVec(std::vector <std::string> tokens) {
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
        return res;
    }
    try {

        if (tokens.size() < 11) {
            throw std::invalid_argument("too few arguments in line");
        }
        unsigned long start_locus = stoul(tokens[3]);   // Fourth field from bam file
        int samflag = stoi(tokens[1]);
        std::string seq = tokens[9];
        std::string bp_qual = tokens[10];
        std::string CIGAR = tokens[5];

        seq = clean_seq(seq, CIGAR);

        unsigned long seq_len = seq.length();   // We may need the original length later.
        std::string ref = genome_ref.substr(start_locus - 1, seq_len);
//        std::cerr << "s " << seq << std::endl;
//        std::cerr << "r " << ref << std::endl;
//        print_vec(tokens);

        // build methylation pattern:
        std::string meth_pattern;
        bool reversed;
        if (is_paired_end) {
            reversed = ((samflag == 83) || (samflag == 163));
        } else {
            reversed = ((samflag & 0x0010) == 16);
        }
        int first_ind = compareSeqToRef(seq, ref, reversed, meth_pattern);

        // in case read contains no CpG sites:
        if (meth_pattern.empty()) {
            readsStats.nr_empty++;
            return res;     // return empty vector
        }

        // translate first CpG locus to CpG index
        int start_site = locus2CpGIndex((int) start_locus + first_ind);

        // Push results into res vector and return:
        res.push_back(tokens[2]);                                  // chr
        res.push_back(std::to_string(start_site));                 // start site
        res.push_back(meth_pattern);                               // meth pattern
        res.push_back(std::to_string(stoi(tokens[3])));            // start locus

        // Read length:
        int read_len = std::abs(stoi(tokens[8]));
        if (!read_len) {        // Happens in single-end data
            read_len = int(seq_len);
        }
        res.push_back(std::to_string(read_len));                    // read length

        return res;
    }
    catch (std::exception &e) {
        std::string msg = "[ " + chr + " ] " + "Exception while processing line "
                          + std::to_string(line_i) + ". Line content: \n";
        std::cerr << msg;
        print_vec(tokens);
        std::cerr << e.what() << std::endl;
        readsStats.nr_invalid++;
    }
    res.clear();
    return res; // return empty vector
}


void patter::proc2lines(std::vector <std::string> tokens1,
                        std::vector <std::string> tokens2) {

    std::vector <std::string> l1, l2;

    /** print result to stdout */
    try {
        // sanity check: two lines must have the same QNAME
        if ((!(tokens2.empty())) && (!(tokens1.empty()))
            && (tokens1[0] != tokens2[0])) {
            readsStats.nr_invalid += 2;
            throw std::invalid_argument("lines are not complements!");
        }
        l1 = samLineToPatVec(tokens1);
        l2 = samLineToPatVec(tokens2);

        merge_and_print(l1, l2);
    }
    catch (std::exception &e) {
        std::string msg = "[ " + chr + " ] Exception while merging. lines ";
        msg += std::to_string(line_i) + ". Line content: \n";
        std::cerr << msg;
        print_vec(tokens1);
        print_vec(tokens2);
        std::cerr << e.what() << std::endl;
        return;
    }
}

void patter::procPairAddMethylData(std::vector <std::string> tokens1,
                                   std::vector <std::string> tokens2, std::string line1, std::string line2) {

    std::vector <std::string> l1, l2;

    /** print result to stdout */
    try {
        // sanity check: two lines must have the same QNAME
        if ((!(tokens2.empty())) && (!(tokens1.empty()))
            && (tokens1[0] != tokens2[0])) {
            readsStats.nr_invalid += 2;
            throw std::invalid_argument("lines are not complements!");
        }
        l1 = samLineToPatVec(tokens1);
        l2 = samLineToPatVec(tokens2);

        patter::MethylData res = merge_and_count_methyl_data(l1, l2);
        if (res.countMethyl + res.countUnmethyl < min_cpg) {
            readsStats.nr_short++;
            return;
        }
        if (!tokens1.empty()) {
            std::string toPrint1 = samLineMethyldataMakeString(line1, res);
            std::cout << line1 + toPrint1;
        }
        std::string toPrint2 = samLineMethyldataMakeString(line2, res);
        std::cout << line2 + toPrint2;
    }
    catch (std::exception &e) {
        std::string msg = "[ " + chr + " ] Exception while merging. lines ";
        msg += std::to_string(line_i) + ". Line content: \n";
        std::cerr << msg;
        print_vec(tokens1);
        print_vec(tokens2);
        std::cerr << e.what() << std::endl;
        return;
    }
}


bool patter::first_line(std::string &line) {
    // find out if the input is paired- or single-end
    try {
        // find out if single- or paired-end
        std::vector <std::string> tokens = line2tokens(line);
        auto flag = (uint16_t) stoi(tokens.at(1));
        chr = tokens.at(2);
        return (flag & 1);     // file is paired end iff flag 0x1 is on
    }
    catch (std::exception &e) {
        std::cerr << "[ " + chr + " ]" << "Invalid first line: \n" << line;
        std::cerr << "\nexception: " << e.what() << std::endl;
        throw e;
    }
}


void patter::print_stats_msg() {
    /** print informative summary message */

    int sucess = line_i ? int((1.0 - ((double) readsStats.nr_invalid / line_i)) * 100.0) : 0;

    std::string msg = "[ " + chr + " ] ";
    msg += "finished " + std::to_string(line_i) + " lines. ";
    if (is_paired_end) {
        msg += "(" + std::to_string(readsStats.nr_pairs) + " pairs). ";
    }
    msg += std::to_string(line_i - readsStats.nr_empty - readsStats.nr_short - readsStats.nr_invalid) + " good, ";
    msg += std::to_string(readsStats.nr_empty) + " empty, ";
    msg += std::to_string(readsStats.nr_short) + " short, ";
    msg += std::to_string(readsStats.nr_invalid) + " invalid. ";
    msg += "(sucess " + std::to_string(sucess) + "%)\n";
    std::cerr << msg;
}

void patter::handlyMethylCountSamLine(std::string line) {
    std::vector <std::string> tokens;
    if(line.at(0) == '@'){
        std::cout << line + "\n";
    } else {
        if (genome_ref.empty()) {
            is_paired_end = first_line(line);
            load_genome_ref();
        }
        tokens = line2tokens(line);
        std::string resString = samLineToSamLineWithMethCounts(tokens, line);
        std::cout << line + resString;
    }
}

void patter::print_progress(){
    if (line_i && !(line_i % 5000000))
        std::cerr << "[ " + chr + " ]" << "patter, line " << line_i << std::endl;
}

std::string get_first_non_empty_line(std::istream& in){
    std::string line;
    while (line.empty()){
        std::getline(in, line);
    }
    return line;
}

void patter::proc_sam_in_stream(std::istream& in) {
    std::vector <std::string> tokens1, tokens2;
    std::string line1, line2;
    bool first_in_pair = true;
    for (std::string line_str; std::getline(in, line_str); line_i++) {
        if(line_str.at(0) == '@'){
            std::cout << line_str + "\n";
        } else {
//            // skip empty lines
//            if (line_str.empty())
//                continue;
            initialize_patter(line_str);
            if (first_in_pair && is_paired_end) {
                line1 = line_str;
                first_in_pair = false;
                readsStats.nr_pairs++;
                continue;
            }
            print_progress();
            line2 = line_str;
            first_in_pair = true;
            line_i++;

            proc_pair_sam_lines(line1, line2);
        }

    }
}

void patter::action_sam(std::string samFilePath) {
    /** parse stdin for sam format lines (single- or pired-end).
     * Translate them to pat format, and output to stdout */

    std::ifstream samFile(samFilePath, std::ios::in);

    if (!(samFile)){
        proc_sam_in_stream(std::cin);
    } else {
        if (samFile.is_open()) {
            proc_sam_in_stream(samFile);
            samFile.close();
        }
    }
    print_stats_msg();
}

void patter::proc_pair_sam_lines(std::string &line1, std::string &line2) {
    std::vector<std::string> tokens1 = line2tokens(line1);
    std::vector<std::string> tokens2 = line2tokens(line2);
    print_progress();
    readsStats.nr_pairs++;
    procPairAddMethylData(tokens1, tokens2, line1, line2);
}

void patter::initialize_patter(std::string &line_str) {
    // first line based initializations
    if (genome_ref.empty()) {
        is_paired_end = first_line(line_str);
        load_genome_ref();
    }
}

void patter::action() {
    /** parse stdin for sam format lines (single- or pired-end).
     * Translate them to pat format, and output to stdout */

    if (DEBUG)
        std::cerr << "DEBUG mode ON" << std::endl;

    bool first_in_pair = true;
//    std::ios_base::sync_with_stdio(false);
    std::vector <std::string> tokens1, tokens2;
    for (std::string line_str; std::getline(std::cin, line_str); line_i++) {

        print_progress();

        // DEBUG - break after a few lines
        if (DEBUG && (line_i > 10000))
            break;

        // skip empty lines
        if (line_str.empty())
            continue;

        initialize_patter(line_str);

        // paired-end file, and current row is first out of a couple of rows
        if (first_in_pair && is_paired_end) {
            tokens1 = line2tokens(line_str);
            first_in_pair = false;
            readsStats.nr_pairs++;
            continue;
        }

        // otherwise (second row in couple, or not-paired-end), line is processed:
        tokens2 = line2tokens(line_str);
        first_in_pair = true;   // next line will be first of couple.

        // process couple of lines. write to stdout
        // in case of single-end input file, tokens1 will be empty
        proc2lines(tokens1, tokens2);

    }
    print_stats_msg();
}

bool is_number(const std::string& s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}


int main(int argc, char **argv) {
    clock_t begin = clock();
    try {
//        std::string genome_name = "/cs/cbio/netanel/tools/wgbs_tools/references/hg19/genome.fa";
//        std::string chrom_size_path = "/cs/cbio/netanel/tools/wgbs_tools/references/hg19/CpG.chrome.size";
//        std::string bam_path = "/cs/cbio/jon/projects/PyCharmProjects/wgbs_tools/pipeline_wgbs/check_check.sam";
//        patter p(genome_name, chrom_size_path);
//        p.action_sam(bam_path);
        InputParser input(argc, argv);
        int min_cpg = 1;
        if (input.cmdOptionExists("--min_cpg")){
            std::string min_cpg_string = input.getCmdOption("--min_cpg");
            if ( !is_number(min_cpg_string) )
                throw std::invalid_argument("Invalid min_cpg argument. Should be a non-negative integer.");
            min_cpg = std::stoi(min_cpg_string);
        }
        if (argc == 5) {
            patter p(argv[1], argv[2], min_cpg);
            p.action_sam("");
        } else
            throw std::invalid_argument("Usage: add_cpg_counts GENOME_PATH CPG_CHROM_SIZE_PATH bam --bam --min_cpg");
    }//                std::cerr  << "match maker, line_i " << line_i << std::endl;
    catch (std::exception &e) {
        std::cerr << "Failed! exception:" << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    }
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cerr << elapsed_secs << std::endl;
    return 0;
}

