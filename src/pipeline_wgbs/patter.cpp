//
// Created by nloyfer on 9/14/18.
//

#include "patter.h"

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

void patter::load_genome_ref() {
    /** Load the genome reference, corresponding to chromosome */

    // load the current chromosome sequence from the FASTA file 
    std::string cmd = "tabix " + ref_path + " " + region + " | cut -f2-3";
    std::string cur_chr = exec(cmd.c_str());
    if (cur_chr.length() == 0) {
        // If the region was empty due to lack of CpGs in range, bam2pat.py would have catched that earlier.
        throw std::invalid_argument("[ patter ] Error: Unable to read reference path: " + ref_path);
    }
    std::stringstream ss(cur_chr);

    std::vector<std::string> tokens;
    std::vector<int> loci;
    for (std::string line_str; std::getline(ss, line_str, '\n');) {
        tokens = line2tokens(line_str);
        int locus = stoi(tokens[0]);
        int cpg_ind = stoi(tokens[1]);
        dict.insert(std::make_pair(locus, cpg_ind));
        //std::cerr << locus << "  " << cpg_ind << std::endl;
        loci.push_back(locus);
    }
    offset = loci.at(0);
    //std::cerr << "offset: " << offset << std::endl;
    int bsize = loci.at(loci.size() - 1);

    conv = new bool[bsize]();
    for (int locus: loci) {
        conv[locus] = true;
    }
}

/***************************************************************
 *                                                             *
 *               Niche                                         *
 *                                                             *
 ***************************************************************/

void dump_mbias_helper(mbias_ss *mb, std::string outpath);
void patter::dump_mbias() {
    /** dump mbias info. Not finished or tested. */
    if (mbias_path == "") { return; }
    dump_mbias_helper(mbias_OT, mbias_path + ".OT.txt");
    dump_mbias_helper(mbias_OB, mbias_path + ".OB.txt");

}

void dump_mbias_helper(mbias_ss *mb, std::string outpath) {
    std::ofstream mbias_stream(outpath);
    std::string sep = "";
    mbias_stream << "r1m1\tr1u1\tr2m2\tr2u2\n";
    for (int pos = 0; pos < MAX_READ_LEN; pos++){
        for (int i = 0; i < 2; i++ ) {
            mbias_stream << sep << mb[i].meth[pos] ;
            sep = "\t";
            mbias_stream << sep << mb[i].unmeth[pos];
        }
        sep = "";
        mbias_stream << "\n";
    }
    mbias_stream.close();
}


/***************************************************************
 *                                                             *
 *             Process single / pair of reads                  *
 *                                                             *
 ***************************************************************/

std::string patter::clean_CIGAR(std::string seq, std::string CIGAR) {

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
            throw std::invalid_argument("[ patter ] Unknown CIGAR character: " +
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
        // Should never happen. Probably means a reference mismatch
        throw std::logic_error("[ patter ] Reference Error. Unknown CpG locus: " +
                                std::to_string(locus));
    }
    return start_site;
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

int patter::compareSeqToRef(std::string &seq,
                            int start_locus,
                            int samflag,
                            std::string &meth_pattern) {
    /** compare seq string to ref string. generate the methylation pattern, and return
     * the CpG index of the first CpG site in the seq (or -1 if there is none) */

    // get orientation 
    bool bottom;
    if (is_paired_end) {
        bottom = ( ((samflag & 0x53) == 83) || ((samflag & 0xA3) == 163) );
    } else {
        bottom = ((samflag & 0x10) == 16);
    }
    ReadOrient ro = bottom ? OB : OT;

    // get flag for mbias
    int mbias_ind;
    mbias_ss *mb;
    bool skip_mbias = false;
    if ((samflag & 0x53) == 0x53) {mb = mbias_OB; mbias_ind = 0;} // 83
    else if ((samflag & 0xA3) == 0xA3) {mb = mbias_OB; mbias_ind = 1;} // 163
    else if ((samflag & 0x63) == 0x63) {mb = mbias_OT; mbias_ind = 0;} // 99
    else if ((samflag & 0x93) == 0x93) {mb = mbias_OT; mbias_ind = 1;} // 147
    else { skip_mbias = true; }
    
    // generate the methylation pattern (e.g 'CC.TC'),
    // by comparing the given sequence to reference at the CpG indexes
    char cur_status;
    int j;
    int nr_cpgs = 0;
    int first_ind = -1;
    for (unsigned long i = 0; i < seq.length(); i++) {
        if (i >= MAX_READ_LEN) {skip_mbias = false;}
        if (conv[start_locus + i]) {
            j = i + ro.shift;
            char s = seq[j];
            cur_status = UNKNOWN;
            if (s == ro.unmeth_seq_chr) {
                cur_status = UNMETH;
                if (!skip_mbias) mb[mbias_ind].unmeth[j]++;
            } else if (s == ro.ref_chr) {
                cur_status = METH;
                if (!skip_mbias) mb[mbias_ind].meth[j]++;
            }
            // ignore first/last 'clip_size' characters, since they are often biased
            if (!((j >= clip_size) && (j < seq.size() - clip_size))) {
                cur_status = UNKNOWN;
            }
            if ((first_ind < 0) && (cur_status != UNKNOWN)) {
                first_ind = locus2CpGIndex(start_locus + i); 
            }
            if (first_ind > 0) {
                meth_pattern.push_back(cur_status);
            }
        }
    }
    if (strip_pat(meth_pattern)) {return -1;}
    return first_ind;
}

std::vector <std::string> merge(std::vector<std::string> l1, 
                                std::vector<std::string> l2) {
    /** Merge 2 complementary lines to a single output.
     * each line has the following fields: [chr, startCpG, pat]
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
        if (merged_pat[adj_i] == UNKNOWN) {   // this site was missing from read1
            merged_pat[adj_i] = pat2[i];
        } else if ((pat2[i] != UNKNOWN) && (merged_pat[adj_i] != pat2[i])) {
            // read1 and read2 disagree, and none of them is missing ('.').
            // treat this case as a missing value for now
            // future work: consider only the read with the higher quality.
            merged_pat[adj_i] = UNKNOWN;
        }
    }
    // strip merged pat (remove trailing dots):
    int pos = strip_pat(merged_pat);
    if (pos < 0 ) { return {}; }
    l1[1] = std::to_string(start1 + pos);
    l1[2] = merged_pat;
    return l1;
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
        //std::string bp_qual = tokens[10];
        std::string CIGAR = tokens[5];

        seq = clean_CIGAR(seq, CIGAR);

        // build methylation pattern:
        std::string meth_pattern;

        int start_site = compareSeqToRef(seq, start_locus, samflag, meth_pattern);

        // in case read contains no CpG sites:
        if (start_site < 1) {
            readsStats.nr_empty++;
            return res;     // return empty vector
        }

        // Push results into res vector and return:
        res.push_back(tokens[2]);                                  // chr
        res.push_back(std::to_string(start_site));                 // start site
        res.push_back(meth_pattern);                               // meth pattern

        return res;
    }
    catch (std::exception &e) {
        std::string msg = "[ " + chr + " ] " + "Exception while processing line "
                          + std::to_string(line_i) + ". Line content: \n";
        std::cerr << "[ patter ] " << msg;
        print_vec(tokens);
        std::cerr << "[ patter ] " << e.what() << std::endl;
        readsStats.nr_invalid++;
    }
    res.clear();
    return res; // return empty vector
}

bool are_paired(std::vector <std::string> tokens1,
                std::vector <std::string> tokens2) {
    // return true iff the reads are non empty and paired
    return ((!(tokens2.empty())) && 
            (!(tokens1.empty())) &&
            (tokens1[0] == tokens2[0]));
}

void patter::proc2lines(std::vector <std::string> &tokens1,
                        std::vector <std::string> &tokens2) {

    try {
        // sanity check: two lines must have the same QNAME
        if ((!(tokens2.empty())) && 
            (!(tokens1.empty())) &&
            (tokens1[0] != tokens2[0])){
            readsStats.nr_invalid += 2;
            throw std::invalid_argument("lines are not complements!");
        }

        // Merge 2 complementary lines to a single output. 
        std::vector <std::string> res = merge(samLineToPatVec(tokens1), samLineToPatVec(tokens2));
        if (res.empty()) { return; } 
        if ((res[2]).length() < min_cpg) {
            readsStats.nr_short++;
            return;
        }
        // print to stdout
        std::cout << res[0] + TAB + res[1] + TAB + res[2] + "\n";
    }
    catch (std::exception &e) {
        std::string msg = "[ " + chr + " ] Exception while merging. lines ";
        msg += std::to_string(line_i) + ". Line content: \n";
        std::cerr << "[ patter ] " << msg;
        print_vec(tokens1);
        print_vec(tokens2);
        std::cerr <<  "[ patter ] " << e.what() << std::endl;
        return;
    }
}

/***************************************************************
 *                                                             *
 *                     print stats                             *
 *                                                             *
 ***************************************************************/

void patter::print_stats_msg() {
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
    msg += "(success " + std::to_string(sucess) + "%)\n";
    std::cerr << "[ patter ] " << msg;
}

/***************************************************************
 *                                                             *
 *                Init                                         *
 *                                                             *
 ***************************************************************/

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
        std::cerr << "[ patter ] " << "[ " + chr + " ]" << "Invalid first line: \n" << line;
        std::cerr << "\nexception: " << e.what() << std::endl;
        throw e;
    }
}

void patter::initialize_patter(std::string &line_str) {
    // first line based initializations
    if (!(dict.empty())) { return; }
    is_paired_end = first_line(line_str);
    load_genome_ref();
}

/***************************************************************
 *                                                             *
 *                     Parse bam                               *
 *                                                             *
 ***************************************************************/

void patter::print_progress(){
    if (line_i && !(line_i % 5000000)){
        dump_mbias();
        clock_t tock = clock();
        double elapsed_secs = double(tock - tick) / (CLOCKS_PER_SEC * 60);
        tick = tock;
        std::cerr << "[ patter ] [ " + chr + " ]" << " line " << addCommas(line_i) 
            << " in " << std::setprecision(2) << elapsed_secs << " minutes." << std::endl;
    }
}

void patter::proc1line(std::vector <std::string> &tokens1) {
    proc2lines(tokens1, dummy_tokens);
    tokens1.clear();
}

void patter::action() {
    /** parse stdin for sam format lines (single- or pired-end).
     * Translate them to pat format, and output to stdout */

//    std::ios_base::sync_with_stdio(false);
    std::vector <std::string> tokens1, tokens2;
    for (std::string line_str; std::getline(std::cin, line_str); line_i++) {

        print_progress();

        // skip empty lines
        if (line_str.empty()) { continue; }

        initialize_patter(line_str);

        if (tokens1.empty()) {
            tokens1 = line2tokens(line_str);
            if (!is_paired_end) { 
                proc1line(tokens1);
                tokens1.clear();
            }
            continue;
        }
        tokens2 = line2tokens(line_str);
        if (are_paired(tokens1, tokens2)) {
            proc2lines(tokens1, tokens2);
            readsStats.nr_pairs++;
            tokens1.clear();
        } else {
            proc1line(tokens1);
            tokens1 = tokens2;
        }
    }
    if (! tokens1.empty()) { proc1line(tokens1); }
    print_stats_msg();
    dump_mbias();
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
        int min_cpg = 1;
        if (input.cmdOptionExists("--min_cpg")){
            std::string min_cpg_string = input.getCmdOption("--min_cpg");
            if ( !is_number(min_cpg_string) )
                throw std::invalid_argument("Invalid min_cpg argument. Should be a non-negative integer.");
            min_cpg = std::stoi(min_cpg_string);
        }
        int clip = 0;
        if (input.cmdOptionExists("--clip")){
            std::string clip_str = input.getCmdOption("--clip");
            if ( !is_number(clip_str) )
                throw std::invalid_argument("Invalid clip argument. Should be a non-negative integer.");
            clip = std::stoi(clip_str);
        }
        if (argc < 3) {
            throw std::invalid_argument("Usage: patter CPG_DICT REGION [--mbias MBIAS_PATH] [--clip CLIP]");
        }
        std::string mbias_path = input.getCmdOption("--mbias");
        patter p(argv[1], argv[2], mbias_path, min_cpg, clip);
        p.action();

    }
    catch (std::exception &e) {
        std::cerr << "[ patter ] Failed! exception:" << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    }
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //std::cerr << elapsed_secs << std::endl;
    return 0;
}
