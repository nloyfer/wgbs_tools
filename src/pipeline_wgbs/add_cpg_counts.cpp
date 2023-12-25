//
// Created by jrosensk on 14/11/2021.
//

#include "add_cpg_counts.h"


/***************************************************
 *                                                 *
 *         Region string parsing                   *
 *                                                 *
 ***************************************************/

std::vector<int> get_loci(std::string r){
    /** Take region string and output start, end vector
     *  Example: Input: chr1:1200-4500
     *           Output: [1200, 4500]
     */
    std::vector<int> v;
    std::string loci_str = r.substr(r.find(":") + 1);
    int start = stoi(loci_str.substr(0, loci_str.find('-')));
    v.push_back(start);
    int end = stoi(loci_str.substr(loci_str.find("-") + 1));
    v.push_back(end);
    return v;
}

std::string extend_region(std::string region){
    /** Take region string and extend it by 1,000bp for each side
     *  Example: Input:  chr1:1200-4500
     *           Output: chr1:200-5500
     */
    std::string chrom = region.substr(0, region.find(":"));
    std::vector<int> start_end_loci = get_loci(region);
    int start = std::max(1, start_end_loci.at(0) - 1000);
    int end = start_end_loci.at(1) + 1000;
    return chrom + ":" + std::to_string(start) + "-" + std::to_string(end);
}

/***************************************************
 *                                                 *
 *         Loading reference                       *
 *                                                 *
 ***************************************************/

std::vector<int> patter::load_genome_helper(std::string region, std::string cmd){
    // load the current chromosome sequence from the FASTA file
    std::string cur_chr = exec(cmd.c_str());
    if (cur_chr.length() == 0) {
        // If the region was empty due to lack of CpGs in range, bam2pat.py would have catched that earlier.
        std::cerr << "[add_cpg_counts] Failed command: " << cmd << std::endl;
        throw std::invalid_argument("Error: Unable to read reference path: " + ref_path + " at " + region + ".");
    }
    std::stringstream ss(cur_chr);

    std::vector<std::string> tokens;
    std::vector<int> loci;
    for (std::string line_str; std::getline(ss, line_str, '\n');) {
        tokens = line2tokens(line_str);
        int locus = stoi(tokens[0]);
        int cpg_ind = stoi(tokens[1]);
        dict.insert(std::make_pair(locus, cpg_ind));
        loci.push_back(locus);
    }
    bsize = loci.at(loci.size() - 1);
    conv = new bool[bsize]();
    return loci;
}

int patter::update_conv(std::vector<int> loci, int start, int end, int counter) {
    for (int i = counter; i < loci.size() ;i++) {
        counter = i; // update counter to be same as i for case when we break the loop and want to start from last index
        int locus = loci.at(i);
        if (locus < start) { continue; } // ignore all locus that no in region [start, end]
        if (locus >= start && locus <= end){
            conv[locus] = true;
        }
        if (locus > end){ break; } // break loop in order to get update [start, end]
    }
    return counter;
}

void patter::load_genome_ref() {
    /** Load the genome reference, corresponding to chromosome */

    if (bed_file.empty()) { // CASE region
        if (region.find(':') != std::string::npos) { // CASE for region=chr:start-end

            std::string extend_r = extend_region(region);
            std::string cmd = "tabix " + ref_path + " " + extend_r + " | cut -f2-3";
            std::vector <int> loci = load_genome_helper(extend_r, cmd); // init dict base on extended region (for reads with more cpgs than the region)

            // init conv base on original region
            std::vector<int> start_end_loci = get_loci(region);
            int start = start_end_loci.at(0);
            int end = start_end_loci.at(1);
            update_conv(loci, start, end, 0);
        } else { // CASE for region=chr
            std::string cmd = "tabix " + ref_path + " " + region + " | cut -f2-3";
            std::vector<int> loci = load_genome_helper(region, cmd);
            for (int locus: loci) {
                conv[locus] = true;
            }
        }
    } else { // CASE for bed file
        std::string chrom = region;  // in bed file case: region=chr
        std::string cmd = "tabix " + ref_path + " -R " + bed_ext_file + " | grep -w " + chrom  + " | cut -f2-3";
        // init dict base on extended bed file (for reads with more cpgs than the original region)
        std::vector<int> loci = load_genome_helper(chrom, cmd);

        // init conv base on original bed file
        cmd = "cat " + bed_file + " | grep -w " + chrom;
        std::string bed_chrs = exec(cmd.c_str());
        if (bed_chrs.length() == 0) {
            std::cerr << "no markers for chrom " << chrom << std::endl;
        }
        std::stringstream bed_chrs_ss(bed_chrs);
        int counter = 0;
        int prev_end = 0;
        int prev_start = 0;
        std::vector<std::string> tokens;
        for (std::string line_str; std::getline(bed_chrs_ss, line_str, '\n');) {
            tokens = line2tokens(line_str);
            int start = stoi(tokens[1]);
            int end = stoi(tokens[2]);
            // validate bed regions
            if (start >= end){
                throw std::invalid_argument("[add_cpg_count] [" + chrom + "] Error: bed file is not valid, start larger than end: " + line_str);
            }
            if (start < prev_start){
                throw std::invalid_argument("[add_bed_count] [" + chrom + "] Error: bed file is not sorted");
            }
            counter = update_conv(loci, start, end, counter);
            prev_end = end;
            prev_start = start;
        }
    }

}

/***************************************************
 *                                                 *
 *         patter stuff                            *
 *                                                 *
 ***************************************************/


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

int patter::compareSeqToRef(std::string &seq,
                            int start_locus,
                            int samflag,
                            std::string &meth_pattern) {
    /** compare seq string to ref string. generate the methylation pattern, and return
     * the CpG index of the first CpG site in the seq (or -1 if there is none) */

    // get orientation
    bool bottom = is_bottom(samflag, is_paired_end);
    ReadOrient ro = bottom ? OB : OT;

    // generate the methylation pattern (e.g 'CC.TC'),
    // by comparing the given sequence to reference at the CpG indexes
    char cur_status;
    int j;
    int nr_cpgs = 0;
    int first_ind = -1;
    for (unsigned long i = 0; i < seq.length(); i++) {

        // this deals with the case where a read exceeds
        // the last CpG of the chromosome. Ignore the rest of the read.
        if (start_locus + 1 > (bsize - 1)) {
            continue;
        }
        if (conv[start_locus + i]) {
            j = i + ro.shift;
            char s = seq[j];
            cur_status = UNKNOWN;
            if (s == ro.unmeth_seq_chr) {
                cur_status = UNMETH;
            } else if (s == ro.ref_chr) {
                cur_status = METH;
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
        if (merged_pat[adj_i] == UNKNOWN) {   // this site was missing from read1
            merged_pat[adj_i] = pat2[i];
        } else if ((pat2[i] != UNKNOWN) && (merged_pat[adj_i] != pat2[i])) {
            // read1 and read2 disagree. treat this case as a missing value for now
            // future work: consider only the read with the higher quality.
            merged_pat[adj_i] = UNKNOWN;
        }
    }
    int pos = strip_pat(merged_pat);
    if (pos < 0 ) { return {}; }
    l1[1] = std::to_string(start1 + pos);
    l1[2] = merged_pat;
    return l1;
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
    res.pattern = meth_pattern;
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
        res.pattern = "";
        return res;
    }
    if (l1.empty()) {
        return meth_pattern_count(l2[2]);
    }
    if (l2.empty()) {
        return meth_pattern_count(l1[2]);
    }
    std::vector <std::string> l1a = merge(l1, l2);
    if (l1a.empty()){
        patter::MethylData res;
        res.countMethyl = 0;
        res.countUnmethyl = 0;
        res.pattern = "";
        return res;
    }
    patter::MethylData res = meth_pattern_count(l1a[2]);
    //2 is the index of the meth_pattern
    return res;
}

std::string patter::samLineMethyldataMakeString(std::string originalLine, patter::MethylData md) {
    std::string intermediate = '\t' + TAGNAMETYPE + std::to_string(md.countMethyl) + ',' + std::to_string(md.countUnmethyl);
    if (print_pat){
        return  intermediate + '\t' + PATTAGNAMETYPE + md.pattern + '\n';
    } else {
        return intermediate + '\n';
    }
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
        res.push_back(std::to_string(stoi(tokens[3])));            // start locus

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

void patter::proc1line(std::vector <std::string> &tokens1, std::string line1) {
    procPairAddMethylData(tokens1, dummy_tokens, line1, "");
    tokens1.clear();
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
        std::string toPrint1 = samLineMethyldataMakeString(line1, res);
        std::cout << line1 + toPrint1;
        if (!tokens2.empty()) {
            std::string toPrint2 = samLineMethyldataMakeString(line2, res);
            std::cout << line2 + toPrint2;
        }
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

    if (chr == "") {return;}
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

void patter::print_progress(){
    if (line_i && !(line_i % 5000000)){
        clock_t tock = clock();
        double elapsed_secs = double(tock - tick) / (CLOCKS_PER_SEC * 60);
        tick = tock;
        std::cerr << "[ " + chr + " ]" << " line " << addCommas(line_i)
                  << " in " << std::setprecision(2) << elapsed_secs << " minutes." << std::endl;
    }
}

void patter::proc_sam_in_stream(std::istream& in) {
    std::vector <std::string> tokens1, tokens2;
    std::string line1, line2;
    for (std::string line_str; std::getline(in, line_str); line_i++) {
        if(line_str.at(0) == '@'){
            std::cout << line_str + "\n";
        } else {
            print_progress();
            // skip empty lines
            if (line_str.empty()) { continue; }
            initialize_patter(line_str);
            if (tokens1.empty()) {
                tokens1 = line2tokens(line_str);
                line1 = line_str;
                if (!is_paired_end) {
                    proc1line(tokens1, line1);
                    tokens1.clear();
                    line1 = "";
                }
                continue;
            }
            tokens2 = line2tokens(line_str);
            line2 = line_str;
            if (are_paired(tokens1, tokens2)) {
                procPairAddMethylData(tokens1, tokens2, line1, line2);
                readsStats.nr_pairs++;
                tokens1.clear();
                line1 = "";
            } else {
                proc1line(tokens1, line1);
                tokens1 = tokens2;
                line1 = line2;
            }
        }

    }
    if (! tokens1.empty()) { proc1line(tokens1, line1); }
}

void patter::action_sam(std::string samFilePath) {
    /** parse stdin for sam format lines (single- or pired-end).
     * Translate them to pat format, and output to stdout */

    std::ifstream samFile(samFilePath, std::ios::in);

    if (!(samFile)){
        proc_sam_in_stream(std::cin);
    } else if (samFile.is_open()) {
        proc_sam_in_stream(samFile);
        samFile.close();
    }
    print_stats_msg();
}

void patter::proc_pair_sam_lines(std::string &line1, std::string &line2) {
    std::vector<std::string> tokens1 = line2tokens(line1);
    std::vector<std::string> tokens2 = line2tokens(line2);
    readsStats.nr_pairs++;
    procPairAddMethylData(tokens1, tokens2, line1, line2);
}

void patter::initialize_patter(std::string &line_str) {
    // first line based initializations
    if (!(dict.empty())) { return; }
    is_paired_end = first_line(line_str);
    load_genome_ref();
}

/***************************************************
 *                                                 *
 *                  Main                           *
 *                                                 *
 ***************************************************/

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
        bool print_pat = false;
        if (input.cmdOptionExists("--pat")){
            print_pat = true;
        }
        std::string bed_file = input.getCmdOption("--bed_file");
        std::string bed_ext_file = input.getCmdOption("--bed_ext_file");
        if (argc >= 4) {
            patter p(argv[1], argv[2], min_cpg, clip, print_pat, bed_file, bed_ext_file);
            p.action_sam("");
        } else{

            throw std::invalid_argument("Usage: add_cpg counts CPG_CHROM_SIZE_PATH bam_input --bed_file BED_FILE --bed_ext_file BED_EXT_FILE [--min_cpg MIN_CPG_PER_READ] [--clip CLIP] ");
        }
    }
    catch (std::exception &e) {
        std::cerr << "Failed! exception:" << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    }
    return 0;
}

