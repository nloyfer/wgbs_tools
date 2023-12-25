
#include "patter.h"


bool is_this_CG(std::string &seq, int pos) {
    /* check if current letter is a "C", and the next one is a "G"
       But allow one or more "N"s to be inserted inbetween */
    // is the first letter a 'C'?
    if (seq[pos] != 'C') { return false; }
    // does the sequence end with that C?
    if (pos >= seq.length() - 1) { return false; }
    // check if the next letter is 'G'. Allow one or more N's inbetween
    // i.e., CG, CNG, CNNNNNG are all valid CpG sites
    pos++;
    for (; pos < seq.length(); pos++) {
        if      (seq[pos] == 'G') { return true; }
        else if (seq[pos] == 'N') { continue; }
        else                      { return false; }
    }
    return false;
}

void find_Cm_section(std::string &MM_str, std::string &ML_str) {
    std::vector<std::string> sMM = split_by_semicolon(MM_str);
    int i = 0;
    bool found = false;
    for (i = 0; i < sMM.size(); i++) {
        std::string sv_field = sMM[i];
        if (!((sv_field.substr(0, 5) == "C+m?,") ||
              (sv_field.substr(0, 5) == "C+m.,") ||
              (sv_field.substr(0, 4) == "C+m,"))) {
            continue;
        }
        else {
            found = true;
            break;
        }
    }
    if (!found) {
        throw std::invalid_argument("Unsupported MM field:" + MM_str);
    }
    MM_str = sMM.at(i);
    // trim beginning of MM_str
    std::string MM_str_s = MM_str.substr(MM_str.find_first_of(',') + 1, MM_str.length());
    std::vector<int> MM_vals = split_by_comma(MM_str_s);
    // trim beginning of ML_str
    ML_str = ML_str.substr(ML_str.find_first_of(',') + 1, ML_str.length());
    std::vector<int> ML_vals = split_by_comma(ML_str);
    int nr_MM_vals = MM_vals.size();
    int nr_ML_vals = ML_vals.size();
    if (!(nr_ML_vals % nr_MM_vals == 0)) {
        std::cerr << "Error in find_Cm: not modulu 0" << std::endl;;
        throw std::invalid_argument("Unsupported MM field:" + MM_str);
    }
    //std::cerr << "i:\n" << i << std::endl;;
    std::vector<int> sliced_ML_vals(ML_vals.begin() + (i * nr_MM_vals), ML_vals.begin() + ((i + 1) * nr_MM_vals));

    // vector of ints to string:
    std::ostringstream oss;
    // Iterate through the vector
    for (size_t j = 0; j < sliced_ML_vals.size(); ++j) {
        // Convert each integer to string and append to the stream
        oss << sliced_ML_vals[j];

        // Add a comma if it's not the last element
        if (j < sliced_ML_vals.size() - 1) {
            oss << ",";
        }
    }
    // Convert the stringstream to a string
    ML_str = "," + oss.str();
}



std::vector <std::string> patter::np_samLineToPatVec(std::vector <std::string> tokens) {
    /** same as samLineToPatVec, but for nanopore format  */
    std::vector<std::string> empty_vec;

    // get MM and ML fields
    std::vector<std::string> np_fields = get_np_fields(tokens);
    std::string valMM = np_fields[0];
    std::string valML = np_fields[1];

    if ((valMM == "") || (valML == "") || (tokens[9] == "*")) {
        readsStats.nr_empty++;
        return empty_vec;
    }

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
        else { 
            cur_status = UNKNOWN;

            if (np_mask[i] == 'M') {
                cur_status = METH;
            } else if (np_mask[i] == 'U') {   
                cur_status = UNMETH;
            }
            MM_vals_ind++;
        }

        // make sure the current letter is a "C", and the next one is a "G" (allow N's inbetween)
        if (!(is_this_CG(seq, i))) {
            cur_status = UNKNOWN;
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


std::string patter::parse_ONT(std::vector <std::string> tokens) {

    // get MM and ML fields
    std::vector<std::string> np_fields = get_np_fields(tokens);
    std::vector<int> MM_vals;
    std::vector<int> ML_vals;
    parse_np_fields(np_fields, MM_vals, ML_vals);

    int samflag = stoi(tokens[1]);
    bool bottom = ((samflag & 0x10) == 16);
    std::string work_seq = bottom ? reverse_comp(tokens[9]) : tokens[9];
    int C_counter = 0;
    int MM_vals_ind = 0;
    std::string np_mask(work_seq.length(), 'E'); // np_mask == "EEE..."
    for (int i = 0; i < work_seq.length(); i++ ) {

        // not a C, not interesting
        if (work_seq[i] == 'C') {
            if ( MM_vals_ind >= MM_vals.size()) { break;}
            char cur_status = 'N';
            if (C_counter == MM_vals.at(MM_vals_ind)) {
                if (ML_vals[MM_vals_ind] > (255 * np_thresh)) { cur_status = 'M'; }
                else if (ML_vals[MM_vals_ind] < (255 * (1 - np_thresh))) { cur_status = 'U'; }
                np_mask[i] = cur_status;
                MM_vals_ind++;
            }
            C_counter++;
            if (C_counter > MM_vals.at(MM_vals.size() - 1)) { break; }
        }
    }
    if (bottom) {
        std::reverse(np_mask.begin(), np_mask.end()); 
        // shift left by 1 pos
        np_mask = np_mask.substr(1, np_mask.length()) + "E";
    }
    return np_mask;
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
    find_Cm_section(MM_str, ML_str);
    np_dot = (MM_str[3] != '?'); // differentiate "C+m." and "C+m" from "C+m?"
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
