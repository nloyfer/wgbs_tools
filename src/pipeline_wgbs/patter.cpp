//
// Created by nloyfer on 9/14/18.
//

#include "patter.h"

char METH = 'C';
char UNMETH = 'T';
//char UNKNOWN = '.';
std::string TAB = "\t";

// TODO: smarter threshold. 
// See https://github.com/nanoporetech/modkit/blob/master/filtering.md

/***************************************************************
 *                                                             *
 *               Load FASTA                                    *
 *                                                             *
 ***************************************************************/

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
        loci.push_back(locus);
    }

    bsize = loci.at(loci.size() - 1);

    conv = new bool[bsize]();
    for (int locus: loci) {
        conv[locus] = true;
    }
}

/***************************************************************
 *                                                             *
 *               M-bias                                        *
 *                                                             *
 ***************************************************************/

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

void patter::dump_mbias() {
    /** dump mbias info. Not finished or tested. */
    if (mbias_path == "") { return; }
    dump_mbias_helper(mbias_OT, mbias_path + ".OT.txt");
    dump_mbias_helper(mbias_OB, mbias_path + ".OB.txt");

}



/***************************************************************
 *                                                             *
 *             Process single / pair of reads                  *
 *                                                             *
 ***************************************************************/

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
    if (is_paired_end) {
        if ((samflag & 0x53) == 0x53) {mb = mbias_OB; mbias_ind = 0;} // 83
        else if ((samflag & 0xA3) == 0xA3) {mb = mbias_OB; mbias_ind = 1;} // 163
        else if ((samflag & 0x63) == 0x63) {mb = mbias_OT; mbias_ind = 0;} // 99
        else if ((samflag & 0x93) == 0x93) {mb = mbias_OT; mbias_ind = 1;} // 147
        else { skip_mbias = true; }
    } else {
        mbias_ind = 0;
        mb = bottom ? mbias_OB : mbias_OT;
    }
    
    // generate the methylation pattern (e.g 'CC.TC'),
    // by comparing the given sequence to reference at the CpG indexes
    char cur_status;
    int j;
    int mj;
    int first_ind = -1;
    for (unsigned long i = 0; i < seq.length(); i++) {
        
        // this deals with the case where a read exceeds
        // the last CpG of the chromosome. Ignore the rest of the read.
        if ((start_locus + i) > (bsize - 1)) { continue;  }
        mj = bottom ? (seq.length() - i - 1) : i;
        if (mj >= MAX_READ_LEN) {skip_mbias = true;} // read is too long. skip it to avoid seg. fault
        if (conv[start_locus + i]) {
            j = i + ro.shift;  // add a 1-pos shift for bottom strands
            char s = seq[j];
            cur_status = UNKNOWN;
            if (s == ro.unmeth_seq_chr) {
                cur_status = UNMETH;
                if (!skip_mbias) {
                    mb[mbias_ind].unmeth[mj]++;
                }
            } 
            else if (s == ro.ref_chr) {
                cur_status = METH;
                if (!skip_mbias) {
                    mb[mbias_ind].meth[mj]++;
                }
            }

            // ignore first/last 'clip_size' characters, since they are often biased
            // TODO: allow 4 different clip sizes (OT,OB,CTOT,CTOB)
            if (!((j >= clip_size) && (j < seq.size() - clip_size))) {
            //if (!(j < seq.size() - clip_size)) {
                cur_status = UNKNOWN;
            }
            // find first CpG index
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

bool is_this_CG(std::string &seq, int pos, bool bottom) {
    // make sure the current letter is a "C", and the next one is a "G"
    if (seq[pos] != 'C') {
        //std::cerr << "first is not C:" << seq[pos] << ". " << std::endl;
        return false;
    }
    if (pos >= seq.length() - 1) {
        //std::cerr << "pos is to big" << pos << ". " << std::endl;
        return false;
    }
    pos++;
    for (; pos < seq.length(); pos++) {
        //std::cerr << "pos: " << pos << ". seq[pos]: " << seq[pos] << std::endl;
        if (seq[pos] == 'G') {
            return true;
        } else if (seq[pos] == 'N') {
            //if (bottom) {
                //return false;
            //}
            continue;
        } else {
            return false;
        }

    }
    return false;
}

std::vector <std::string> patter::np_samLineToPatVec(std::vector <std::string> tokens) {
    /** same as samLineToPatVec, but for nanopore format
     *  */
    std::vector<std::string> empty_vec;

    // get MM and ML fields
    std::vector<std::string> np_fields = get_np_fields(tokens);
    std::string valMM = np_fields[0];
    std::string valML = np_fields[1];

    if ((valMM == "") || (valML == "") || (tokens[9] == "*")) {
        readsStats.nr_empty++;
        return empty_vec;
    }
    //std::cerr << tokens[9] << std::endl;

    unsigned long start_locus = stoul(tokens[3]);   // Fourth field from bam file
    int samflag = stoi(tokens[1]);
    std::string seq = tokens[9];
    std::string CIGAR = tokens[5];
    bool bottom = ((samflag & 0x10) == 16);

    seq = clean_CIGAR(seq, CIGAR);

    std::string np_mask = parse_ONT(tokens);
    np_mask = clean_CIGAR(np_mask, CIGAR);

    std::vector<int> MM_vals;
    std::vector<int> ML_vals;
    parse_np_fields(np_fields, MM_vals, ML_vals);

    if (MM_vals.empty()) { 
        readsStats.nr_empty++;
        return empty_vec; 
    }
    if (bottom) {
        std::reverse(ML_vals.begin(),ML_vals.end());
    }
    //std::cerr << "cseq: " << seq << std::endl;
    // build methylation pattern:
    std::string meth_pattern;
    char cur_status;
    int start_site = -1;
    int MM_vals_ind = 0;
    for (int i = 0; i < seq.length(); i++ ) {
        // this deals with the case where a read exceeds
        // the last CpG of the chromosome. Ignore the rest of the read.
        if ((start_locus + i) > (bsize - 1)) { continue;  }

        // Only consider cytosines in a CpG context
        if (!conv[start_locus + i]) { 
            //std::cerr << start_locus + i << " (" << i << " ) not CpG locus" << std::endl;
            continue;
        }

        // The current C is deleted according to the CIGAR
        if (np_mask[i] == 'N') { 
            cur_status = UNKNOWN;
        }
        // The current C is not mentioned in MM as modified
        else if (np_mask[i] == 'E') { 
            if (np_dot) { cur_status = UNMETH; } 
            else { cur_status = UNKNOWN; }
        } 
        // The current C is mentioned in MM as modified
        //else if (np_mask[i] == '1') { 
        else { 
            cur_status = UNKNOWN;

            if (np_mask[i] == 'M') {
                cur_status = METH;
            } else if (np_mask[i] == 'U') {   
                cur_status = UNMETH;
            }
            MM_vals_ind++;
        }

        // make sure the current letter is a "C", and the next one is a "G"
        if (!(is_this_CG(seq, i, bottom))) {
            cur_status = UNKNOWN;
        }
        //if (!(((seq[i] == 'C') && (i < (seq.length() - 1)) && (seq[i + 1] == 'G')) ||
                    //((i < (seq.length() - 2)) && (seq.substr(i, 3) == "CNG")))) { 
            //if (seq[i + 1] != 'G') {
                //std::cerr << "CNG: " << seq.substr(i, 3) << std::endl;
            //}
            //cur_status = UNKNOWN;
        //}

        if ((i>1000) && (i<1100)){
            //std::cerr << "i= " << i << ". loc: " << start_locus + i << ". mask: " << np_mask[i];
            //std::cerr << ". seq: " << seq[i] << seq[i + 1];
            //std::cerr << " status: " << cur_status << "\n";
        }

        // ignore first/last 'clip_size' characters, since they are often biased
        if (!((i >= clip_size) && (i < seq.size() - clip_size))) {
            cur_status = UNKNOWN;
        }
        // find first cpg index
        if ((start_site < 0) && (cur_status != UNKNOWN)) {
            start_site = locus2CpGIndex(start_locus + i); 
        }
        if (start_site > 0) {
            meth_pattern.push_back(cur_status);
        }
    }
    if (strip_pat(meth_pattern)) {
        readsStats.nr_empty++;
        return empty_vec;
    }

    // in case read contains no CpG sites:
    if (start_site < 1) {
        readsStats.nr_empty++;
        return empty_vec;     // return empty vector
    }
    return pack_pat(tokens[2], start_site, meth_pattern);
}

void patter::parse_np_fields(std::vector<std::string> &np_fields, 
                    std::vector<int> &MM_vals, 
                    std::vector<int> &ML_vals) {
    // validate fields
    std::string MM_str = np_fields[0];
    std::string ML_str = np_fields[1];

    if ((MM_str == "") || (ML_str == "")) {
        throw std::invalid_argument("Missing MM or ML fields");
    }

    // validate MM and trim head
    if (!((MM_str.substr(0, 5) == "C+m?,") ||
          (MM_str.substr(0, 5) == "C+m.,") ||
          (MM_str.substr(0, 4) == "C+m,"))) {
        throw std::invalid_argument("Unsupported MM field:" + MM_str);
    }
    np_dot = (MM_str[3] != '?');
    //np_dot = false; // TODO change
    // differentiate "C+m." and "C+m" from "C+m?"
    MM_str = MM_str.substr(MM_str.find_first_of(',') + 1, MM_str.length());

    // validate ML and trim leading comma
    if (!(ML_str.substr(0, 1) == ",")) {
        throw std::invalid_argument("Unsupported ML field:" + ML_str);
    }
    ML_str = ML_str.substr(1, ML_str.length());

    // split by commas
    std::vector<int> orig_MM_vals = split_by_comma(MM_str);
    ML_vals = split_by_comma(ML_str);

    if (orig_MM_vals.size() != ML_vals.size()) {
        throw std::invalid_argument("Error parsing MM and ML fields.\n" + MM_str + "\n" + ML_str);
    }
    // aggregate MM vals
    int pos = 0;
    for (int i = 0; i < orig_MM_vals.size(); i++) {
        pos += orig_MM_vals[i];
        MM_vals.push_back(pos++);
    }
}

std::string patter::parse_ONT(std::vector <std::string> tokens) {

    // get MM and ML fields
    std::vector<std::string> np_fields = get_np_fields(tokens);
    std::vector<int> MM_vals;
    std::vector<int> ML_vals;
    parse_np_fields(np_fields, MM_vals, ML_vals);

    int samflag = stoi(tokens[1]);
    bool bottom = ((samflag & 0x10) == 16);
    std::string work_seq = bottom ? reverse_comp(tokens[9]) : tokens[9];
    std::string np_mask;
    int C_counter = 0;
    int MM_vals_ind = 0;
    for (int i = 0; i < work_seq.length(); i++ ) {

        // not a C, not interesting
        if (work_seq[i] != 'C') { 
            np_mask += "E";
            continue;
        }

        if (C_counter == MM_vals[MM_vals_ind]) {
            //np_mask += "1";
            if (ML_vals[MM_vals_ind] > (255 * 0.66)) {   // TODO: hard-coded
                np_mask += "M";
            } else if (ML_vals[MM_vals_ind] < (255 * (1 - 0.66))) {
                np_mask  += "U";
            } else {
                np_mask  += "N";
            }
            //np_mask += std::to_string(int(float(ML_vals[MM_vals_ind]) / 25.5));
            MM_vals_ind++;
        } else {
            np_mask += "E";
        }
        C_counter++;
    }
    if (bottom) {
        std::reverse(np_mask.begin(), np_mask.end()); 
        // shift left by 1 pos
        np_mask = np_mask.substr(1, np_mask.length()) + "E";
    }
    return np_mask;
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
        if (is_nanopore) {
            return np_samLineToPatVec(tokens);
        }
        unsigned long start_locus = stoul(tokens[3]);   // Fourth field from bam file
        int samflag = stoi(tokens[1]);
        std::string seq = tokens[9];
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
        return pack_pat(tokens[2], start_site, meth_pattern);
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
        std::vector <std::string> res = merge_PE(samLineToPatVec(tokens1), 
                                                 samLineToPatVec(tokens2));
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

std::vector<std::string> get_np_fields(std::vector <std::string> &tokens) {
    std::string valMM = "";
    std::string valML = "";
    for (auto &j: tokens) {  // TODO: start looking from the 12'th field
        if (("MM:Z:" == j.substr(0, 5)) || ("Mm:Z:" == j.substr(0, 5))) {
            valMM = j.substr(5, j.length());
        }
        else if (("ML:B:C" == j.substr(0, 6)) || ("Ml:B:C" == j.substr(0, 6))) {
            valML = j.substr(6, j.length());
        }
    }
    // return results
    std::vector<std::string> res;
    res.push_back(valMM);
    res.push_back(valML);
    return res;
}

void patter::first_line(std::string &line) {
    // find out if the input is paired- or single-end
    // And if it's a nanopore bam file
    try {
        // find out if single- or paired-end
        std::vector <std::string> tokens = line2tokens(line);
        auto flag = (uint16_t) stoi(tokens.at(1));
        chr = tokens.at(2);
        is_paired_end = (flag & 1);     // file is paired end iff flag 0x1 is on

        // find out if this is a nanopore format (iff MM and ML fields are present)
        std::vector<std::string> npvec = get_np_fields(tokens);
        is_nanopore = (is_nanopore || ((npvec[0] != "") && (npvec[1] != "")));
        //std::cerr << "is nanopore: " << is_nanopore << std::endl;

        if (is_paired_end && is_nanopore) {
            throw std::invalid_argument("Unrecognized bam format: paired end and nanopore");
        }
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
    first_line(line_str);
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

void patter::parse_reads_from_stdin() {
    /** parse stdin for sam format lines (single- or pired-end).
     * Translate them to pat format, and output to stdout */

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

