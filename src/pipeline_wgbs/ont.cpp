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

std::string patter::make_meth_mask(std::string &work_seq) {
    /** make a mask of the methylation status of the C's in the read
     * E - not a C
     * M - methylated        (p > 67%)
     * H - hydroxymethylated (p > 67%)
     * U - unmodified        (p < 33%)
     * N - ambiguous         (33% < p < 67%)
     */

    int C_counter = 0;
    int MM_vals_ind = 0;
    int MM_vals_h_ind = 0;
    std::string np_mask(work_seq.length(), 'E'); // np_mask == "EEE..."
    for (int i = 0; i < work_seq.length(); i++ ) {

        // not a C, not interesting
        if (work_seq[i] != 'C') { continue; }

        char cur_status = 'N';
        if (combine_mods) {
            // Combine 5hmC + 5mC: sum their probabilities and threshold as a single modification
            int h_prob = 0, m_prob = 0;
            bool has_h = (MM_vals_h_ind < MM_vals_h.size()) && (C_counter == MM_vals_h.at(MM_vals_h_ind));
            bool has_m = (MM_vals_ind < MM_vals.size()) && (C_counter == MM_vals.at(MM_vals_ind));
            if (has_h) { h_prob = ML_vals_h[MM_vals_h_ind]; MM_vals_h_ind++; }
            if (has_m) { m_prob = ML_vals[MM_vals_ind]; MM_vals_ind++; }
            if (has_h || has_m) {
                int combined = std::min(h_prob + m_prob, 255);
                if (combined > (255 * np_thresh)) {
                    cur_status = 'M';
                } else if (combined < (255 * (1 - np_thresh))) {
                    cur_status = 'U';
                }
                np_mask[i] = cur_status;
            }
        } else {
            // update 5hmc
            if ((MM_vals_h_ind < MM_vals_h.size()) && (C_counter == MM_vals_h.at(MM_vals_h_ind))) {
                if (ML_vals_h[MM_vals_h_ind] > (255 * np_thresh)) {
                    cur_status = 'H';
                } else if (ML_vals_h[MM_vals_h_ind] < (255 * (1 - np_thresh))) {
                    cur_status = 'U';
                }
                np_mask[i] = cur_status;
                MM_vals_h_ind++;
            }
            // update 5mc (overwrite 5hmc if both are present with high confidence)
            if ((MM_vals_ind < MM_vals.size()) && (C_counter == MM_vals.at(MM_vals_ind))) {
                if (ML_vals[MM_vals_ind] > (255 * np_thresh)) {
                    cur_status = 'M';
                } else if (ML_vals[MM_vals_ind] < (255 * (1 - np_thresh))) {
                    // make sure this is not a 5hmc
                    if (cur_status != 'H') {
                        cur_status = 'U';
                    }
                } else if (cur_status != 'H') {
                    cur_status = 'N';
                }
                np_mask[i] = cur_status;
                MM_vals_ind++;
            }
        }
        C_counter++;
    }
    return np_mask;
}


std::vector <std::string> patter::np_samLineToPatVec(std::vector <std::string> tokens) {
    /** same as samLineToPatVec, but for nanopore format  */
    std::vector<std::string> empty_vec;

    // get MM and ML fields
    parse_np_fields(tokens);

    if (((MM_vals.empty()) && (MM_vals_h.empty() && (!(np_dot)))) || (tokens[9] == "*")) {
        readsStats.nr_empty++;
        return empty_vec;
    }

    unsigned long start_locus = stoul(tokens[3]);   // Fourth field from bam file
    int samflag = stoi(tokens[1]);
    std::string seq = tokens[9];
    std::string CIGAR = tokens[5];
    bool bottom = ((samflag & 0x10) == 16);

    std::string orig_seq = std::string(tokens[9]);
    seq = clean_CIGAR(seq, CIGAR);

    if (bottom) {
        orig_seq = reverse_comp(orig_seq);
    }

    // Create mask (M/U/H/N)
    std::string np_mask = make_meth_mask(orig_seq);
    // debug print
    //std::cout << "[np_samLineToPatVec] Read is bottom strand: " << (bottom ? "true" : "false") << std::endl;

    // If bottom, orig_seq was RC, so np_mask is RC.
    // Flip it to Forward (Ref) orientation before cleaning CIGAR.
    if (bottom) { flip_string(np_mask); }

    np_mask = clean_CIGAR(np_mask, CIGAR);
    //if (bottom) { flip_string(np_mask); }
    // debug print
    //std::cout << "[np_samLineToPatVec] Cleaned np_mask: " << np_mask << std::endl;
    //std::cout << "mask: " << np_mask << std::endl;
    //std::cout << "seq:  " << seq << std::endl;

    
    // build methylation pattern:
    std::string meth_pattern;
    char cur_status;
    int start_site = -1;
    int MM_vals_ind = 0;
    int di = 0;
    // For bottom-strand reads, start at i=-1 to handle the edge case where
    // the read begins at the G position of a CpG (POS = CpG_C + 1). In that
    // case conv[start_locus - 1] == true (the C), but i=0 maps to the G.
    // di = i+1 = 0 is valid; conv[start_locus + (-1)] = conv[CpG_C] = true.
    int loop_start = bottom ? -1 : 0;
    for (int i = loop_start; i < (int)seq.length(); i++ ) {
        di = i;
        if (bottom) { di = i + 1; }
        // this deals with the case where a read exceeds
        // the last CpG of the chromosome. Ignore the rest of the read.
        if ((start_locus + i) > (bsize - 1)) { continue;  }
        // bounds check for negative i (bottom-strand edge case)
        if ((long long)start_locus + i < 0) { continue; }

        // Only consider cytosines in a CpG context
        if (!conv[start_locus + i]) {
            continue;
        }

        // bounds check: di = i+1 for bottom strand can exceed np_mask length
        if (di >= (int)np_mask.size()) { continue; }

        // The current C is deleted according to the CIGAR
        if (np_mask[di] == 'N') {
            cur_status = UNKNOWN;
        }
        // No explicit modification call at this position in MM.
        // With dot convention, unlisted C bases are implicitly unmethylated.
        // But if the read has N or a mismatch (not C/G) at the CpG, there is
        // no base to call — treat as unknown rather than unmethylated.
        else if (np_mask[di] == 'E') {
            bool has_base = (di >= 0 && di < (int)seq.size()) &&
                            (bottom ? seq[di] == 'G' : seq[di] == 'C');
            if (np_dot && has_base) { cur_status = UNMETH; }
            else { cur_status = UNKNOWN; }
        }
        // The current C is mentioned in MM as modified
        else {
            cur_status = UNKNOWN;

            if (np_mask[di] == 'M') {
                cur_status = METH;
            } else if (np_mask[di] == 'U') {
                cur_status = UNMETH;
            } else if (np_mask[di] == 'H') {
                cur_status = HYDROXY;
            }
            MM_vals_ind++;
        }

        // make sure the current letter is a "C", and the next one is a "G" (allow N's inbetween)
        // Skip in nanopore mode: reads are not bisulfite-converted, so conv[] is authoritative
        // Also skip when i<0 (bottom-strand edge case: can't index seq at -1)
        if (i >= 0 && !is_nanopore && !(is_this_CG(seq, i))) {
            cur_status = UNKNOWN;
        }

        // ignore first/last 'clip_size' characters, since they are often biased
        // Use di for bottom strand (di is the actual read position used)
        int clip_pos = bottom ? di : i;
        if (!((clip_pos >= clip_size) && (clip_pos < (int)seq.size() - clip_size))) {
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
    // debug print
    return pack_pat(tokens[2], start_site, meth_pattern);
}

void patter::parse_np_fields(std::vector<std::string> &tokens) {
    /* parse the MM and ML fields for 5hmc&5mc modifications in current read.
     * The MM and ML values for 5hmc are stored in MM_vals_h and ML_vals_h,
     * while the MM and ML values for 5mc are stored in MM_vals and ML_vals.
     * Also parses C+C? section (Biomodal-specific) and merges those positions
     * into MM_vals so they are treated as methylated.
     */
    np_dot = false;  // reset per-read; set only if MM field with '.' convention is found
    parse_np_fields_by_mod(tokens, "h");
    MM_vals_h = MM_vals;
    ML_vals_h = ML_vals;
    parse_np_fields_by_mod(tokens, "m");

    // Save C+m results (np_dot must reflect C+m convention, not C+C)
    bool np_dot_m = np_dot;
    std::vector<int> MM_vals_m = MM_vals;
    std::vector<int> ML_vals_m = ML_vals;

    // Parse C+C? section (Biomodal uses this for additional methylation calls)
    // cpc_call controls how C+C? positions are encoded:
    //   'C' (default): treat as 5mC → merge into MM_vals_m
    //   'H'          : treat as 5hmC → merge into MM_vals_h
    //   '.'          : treat as unknown → skip entirely
    parse_np_fields_by_mod(tokens, "C");
    if (!MM_vals.empty() && cpc_call != '.') {
        // Target vector: MM_vals_m (for 'C') or MM_vals_h (for 'H')
        std::vector<int> &target_vals = (cpc_call == 'H') ? MM_vals_h : MM_vals_m;
        std::vector<int> &target_mls  = (cpc_call == 'H') ? ML_vals_h : ML_vals_m;
        std::set<int> existing(target_vals.begin(), target_vals.end());
        for (int p : MM_vals) {
            if (existing.find(p) == existing.end()) {
                auto it = std::lower_bound(target_vals.begin(), target_vals.end(), p);
                int idx = it - target_vals.begin();
                target_vals.insert(it, p);
                target_mls.insert(target_mls.begin() + idx, 255);
            }
        }
    }

    // Restore C+m results (with C+C? merged in if cpc_call == 'C')
    np_dot = np_dot_m;
    MM_vals = MM_vals_m;
    ML_vals = ML_vals_m;
    return;
}

void patter::parse_np_fields_by_mod(std::vector<std::string> &tokens, 
                                    std::string mod_char) {
    /* parse the MM and ML fields for a specific modification 
     * ("C+m" or "C+h") and update the MM_vals and ML_vals vectors accordingly. 
     * This function is called twice in parse_np_fields: 
     * first for "h" (5hmc) and then for "m" (5mc). 
     * The results for 5hmc are stored in MM_vals_h and ML_vals_h, 
     * while the results for 5mc are stored in MM_vals and ML_vals.
     */

    MM_vals.clear();
    ML_vals.clear();

    // get MM and ML fields from the tokens
    std::string MM_str = "";
    std::string ML_str = "";
    if (!get_np_tags(tokens, MM_str, ML_str)) { return; }

    // subset MM and ML to the requested subsection (C+m or C+h)
    subset_to_Cm_section(MM_str, ML_str, np_dot, mod_char);

    // split MM and ML by commas
    std::vector<int> orig_MM_vals = split_by_comma(trim_from_first_comma(MM_str));
    if (ML_str == "") {
        // make a dummy ML vector of 255's (BioModal logic):
        ML_vals = std::vector<int>(orig_MM_vals.size(), 255);
    } else {
        ML_vals = split_by_comma(trim_from_first_comma(ML_str));
        if (orig_MM_vals.size() != ML_vals.size()) {
            throw std::invalid_argument("Error parsing MM and ML fields.\n" + MM_str + "\n" + ML_str);
        }
    }

    // aggregate MM vals
    int pos = 0;
    for (int i = 0; i < orig_MM_vals.size(); i++) {
        pos += orig_MM_vals[i];
        MM_vals.push_back(pos++);
    }
}

std::string find_Cm_substring(std::string &MM_str, int &i,
                              std::string mod_char) {
    /* Given the MM string (e.g. "C+m,0,1,2,3;C+h,0,2"), 
     * find the substring for the specific mod (e.g. "C+m") 
     * and return it (e.g. "C+m,0,1,2,3").
     * If the specific mod is not found, return empty string.
     */
    std::vector<std::string> sMM = split_by_semicolon(MM_str);
    bool found = false;
    for (i = 0; i < sMM.size(); i++) {
        std::string sv_field = sMM[i];
        if (!(sv_field.substr(0, 3) == "C+" + mod_char)) {
            continue;
        }
        else {
            found = true;
            break;
        }
    }
    if (!found) {
        return "";
    }
    return sMM[i];
}

std::string trim_from_first_comma(std::string &str) {
    /* Given a string like "C+m,0,1,2,3", return the substring after the first comma, i.e. "0,1,2,3"
     * If there is no comma, return empty string. If the input string is empty, return empty string.
     */
    if (str == "") { return str; }
    size_t comma_pos = str.find_first_of(',');
    if (comma_pos == std::string::npos) return "";
    str = str.substr(comma_pos + 1, str.length());
    return str;
}

std::string int_vec_to_str(std::vector<int> &vec) {
    // vector of ints to string:
    std::ostringstream oss;
    for (size_t j = 0; j < vec.size(); ++j) {
        // Convert each integer to string and append to the stream
        oss << vec[j];
        // Add a comma if it's not the last element
        if (j < vec.size() - 1) {
            oss << ",";
        }
    }
    // Convert the stringstream to a string
    return oss.str();
}

void subset_to_Cm_section(std::string &MM_str, std::string &ML_str, 
                          bool &np_dot, std::string mod_char) {
    /* subset the MM and ML strings to the section corresponding to the
     * specific mod ("C+m" or "C+h") and update the np_dot.
     *
     * For example, if mod_char == "m",
     *    MM_str: "C+m,0,1,2,3;C+h,0,2" -> "0,1,2,3"
     *    ML_str: "ML:B:C,255,255,255,255;C+h,255,255" -> "255,255,255,255"
     *    np_dot: false
     *
     *  and if mod_char == "h",
     *   MM_str: "C+m,0,1,2,3;C+h,0,2" -> "0,2"
     *   ML_str: "ML:B:C,255,255,255,255;C+h,255,255" -> "255,255"
     *   np_dot: false
     */

    // find MM_str for C+m
    int MM_pos = 0;
    std::string sub_MM = find_Cm_substring(MM_str, MM_pos, mod_char);
    if (sub_MM == "") {
        MM_str = "";
        ML_str = "";
        return;
    }
    MM_str = sub_MM;
    np_dot = (!((MM_str.size() > 3) && (MM_str[3] == '?'))); // differentiate "C+m." and "C+m" from "C+m?"
    
    // Use local copies so MM_str/ML_str are not modified here;
    // parse_np_fields_by_mod will call trim_from_first_comma on them afterward.
    std::string mm_copy = MM_str;
    std::vector<int> MM_vec = split_by_comma(trim_from_first_comma(mm_copy));
    std::string ml_copy = ML_str;
    std::vector<int> ML_vec = split_by_comma(trim_from_first_comma(ml_copy));

    if (ML_str == "") { return; } // No ML field (e.g., BioModal)
    
    int nr_MM_vals = MM_vec.size();
    if (nr_MM_vals == 0) {
        ML_str = "";
        return;
    }
    // If ML exists, check modulus
    if (!(ML_vec.size() % nr_MM_vals == 0) && ML_vec.size() > 0) {
        // Only throw if structure is fundamentally broken
        std::cerr << "Error in find_Cm: not modulu 0" << std::endl;;
        throw std::invalid_argument("Unsupported MM field:" + MM_str);
    }
    
    // slice the ML vector
    // Note: This assumes ML contains concatenated blocks for all mods in MM
    if (ML_vec.size() >= (MM_pos + 1) * nr_MM_vals) {
        std::vector<int> sliced_ML_vals(ML_vec.begin() + (MM_pos * nr_MM_vals), 
                                        ML_vec.begin() + ((MM_pos + 1) * nr_MM_vals));
        ML_str = "," + int_vec_to_str(sliced_ML_vals);
    }
}

bool get_np_tags(std::vector <std::string> &tokens,
                 std::string &MM_str,
                 std::string &ML_str) {
    /* Find the MM and ML fields and return their values as strings.
     * e.g., "MM:Z:C+m,0,1,2,3;C+h,0,2" -> "C+m,0,1,2,3;C+h,0,2"
     *       "ML:B:C,255,255,255,255;C+h,255,255" -> "C+m,255,255,255,255;C+h,255,255"
     */

    int fieldCount = 0;
    for (auto &j: tokens) { 
        if (fieldCount++ < 11) continue; // Skip the first 11 fields
        if (("MM:Z:" == j.substr(0, 5)) || ("Mm:Z:" == j.substr(0, 5))) {
            MM_str = j.substr(5, j.length());
        }
        else if (("ML:B:C" == j.substr(0, 6)) || ("Ml:B:C" == j.substr(0, 6))) {
            ML_str = j.substr(6, j.length());
        }
    }
    if (MM_str == "") { return false; }
    return true;
}
