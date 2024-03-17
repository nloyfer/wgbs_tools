
#include "patter_utils.h"

/***************************************************
 *                                                 *
 *                 Print methods                   *
 *                                                 *
 **************************************************/

std::vector <std::string> line2tokens(const std::string &line) {
    /** Break string line to words (a vector of string tokens) */
    std::vector <std::string> result;
    std::string cell;
    std::stringstream lineStream(line);
    while (getline(lineStream, cell, '\t')) {
        result.push_back(cell);
    }
    return result;
}

void print_vec(const std::vector <std::string> &vec) {
    /** print a vector to stderr, tab separated */
    const auto separator = "\t";
    auto sep = "";
    for (const auto &j: vec) {
        std::cerr << sep << j;
        sep = separator;
    }
    std::cerr << std::endl;
}

void print_vec(const std::vector <int> &vec) {
    /** print a vector to stderr, tab separated */
    const auto separator = "\t";
    auto sep = "";
    for (const auto &j: vec) {
        std::cerr << sep << j;
        sep = separator;
    }
    std::cerr << std::endl;
}

void output_vec(const std::vector <std::string> &vec) {
    /** print a vector to stdout, tab separated */
    const auto separator = "\t";
    auto sep = "";
    for (const auto &j: vec) {
        std::cout << sep << j;
        sep = separator;
    }
    std::cout << std::endl;
}

std::string addCommas(const int num) {
    /** convert integer to string with commas */
    auto s = std::to_string(num);
    int n = s.length() - 3;
    while (n > 0) {
        s.insert(n, ",");
        n -= 3;
    }
    return s;
}

/***************************************************
 *                                                 *
 *                 Split strings                   *
 *                                                 *
 **************************************************/

std::vector<float> split_float_by_comma(std::string str_line) {
    /** split a comma separated list of floats,
     * and output it as vector of floats */
    std::vector<float> prob_vec;
    std::stringstream ss(str_line);

    for (float i; ss >> i;) {
        prob_vec.push_back(i);
        if (ss.peek() == ',') { ss.ignore(); }
    }
    return prob_vec;
}

std::vector<int> split_by_comma(std::string str_line) {
    /** split a comma separated list of ints,
     * and output it as vector of ints */
    std::vector<int> prob_vec;
    std::stringstream ss(str_line);

    for (int i; ss >> i;) {
        prob_vec.push_back(i);
        if (ss.peek() == ',') { ss.ignore(); }
    }
    return prob_vec;
}

std::vector<std::string> split_by_semicolon(std::string str_line) {
    /** split a semicolon separated list of strings,
     * and output it as vector of strings */
    std::vector<std::string> str_vec;
    std::istringstream ss(str_line);
    std::string token;

    while (std::getline(ss, token, ';')) {
        str_vec.push_back(token);
    }
    return str_vec;
}

/***************************************************
 *                                                 *
 *                 General utils                   *
 *                                                 *
 **************************************************/

bool hasEnding(std::string const &fullString, std::string const &suffix) {
    /** return true iff "fullString" ends with "suffix" */
    if (fullString.length() >= suffix.length()) {
        return (0 == fullString.compare(fullString.length() - suffix.length(), suffix.length(), suffix));
    } else {
        return false;
    }
}

bool is_number(const std::string& s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

/***************************************************
 *                                                 *
 *                 Load from exec                  *
 *                                                 *
 **************************************************/

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


/***************************************************
 *                                                 *
 *                sam read methods                 *
 *                                                 *
 **************************************************/

bool is_bottom(int samflag, bool is_paired_end) {
    if (is_paired_end) { 
        return (((samflag & 0x53) == 83) || ((samflag & 0xA3) == 163));
    };
    return ((samflag & 0x10) == 16);
}

bool are_paired(std::vector <std::string> tokens1,
                std::vector <std::string> tokens2) {
    // return true iff the reads are non empty and paired
    return ((!(tokens2.empty())) && 
            (!(tokens1.empty())) &&
            (tokens1[0] == tokens2[0]));
}


std::string reverse_comp(std::string seq) {
    std::string revcomp = "";
    char oc, nc;
    for (int i = seq.length() - 1; i > -1; i--) {
        oc = seq[i];
        if (oc == 'A') {
            nc = 'T';
        } else if (oc == 'C') {
            nc = 'G';
        } else if (oc == 'G') {
            nc = 'C';
        } else if (oc == 'T') {
            nc = 'A';
        } else {
            throw std::runtime_error("[ patter ] Unsupported base");
        }
        revcomp += nc;
    }
    return revcomp;
}

/***************************************************
 *                                                 *
 *                     CIGAR                       *
 *                                                 *
 **************************************************/

std::string clean_CIGAR(std::string seq, std::string CIGAR) {

    /** use CIGAR string to adjust 'seq' so it will be comparable to the reference.
     * e.g, remove false letters ('I'), insert fictive letters ('D') etc. */

    // parse CIGAR and convert it to a couple of vectors: chars, nums.
    // e.g, '2S9M' will become ['S', 'M'] and [2, 9]
    std::vector<char> chars;
    std::vector<int> nums;
    std::string cur_num;
    for (auto c: CIGAR) {
        if (isdigit(c)) {
            cur_num.push_back(c);
        } else {
            nums.push_back(stoi(cur_num));
            cur_num.clear();
            chars.push_back(c);
        }
    }

    // build the adjusted seq
    std::string adjusted_seq;
    for (int i = 0; i < (int) chars.size(); i++) {
        char ch = chars[i]; // CIGAR character
        int num = nums[i];  // corresponding integer
        if ((ch == 'M') || (ch == '=') || (ch == 'X')) {
            adjusted_seq += seq.substr(0, num);
            seq = seq.substr(num, seq.length() - num);
        } else if ((ch == 'D') || (ch == 'N')) {
            for (unsigned long j = 0; j < num; j++)
                adjusted_seq += 'N';
        } else if ((ch == 'I') || (ch == 'S')) {
            seq = seq.substr(num, seq.length() - num);
        } else if (ch == 'H') {
            continue;
        } else {
            throw std::invalid_argument("[ patter ] Unknown CIGAR character: " +
                                        std::string(1, ch));
        }
    }

    return adjusted_seq;
}


/***************************************************
 *                                                 *
 *                pat read methods                 *
 *                                                 *
 **************************************************/

void strip_read(std::vector <std::string> &tokens) {
    // strip pat (remove trailing dots), and update the start position
    int pos = strip_pat(tokens[2]);
    if (pos < 0 ) { 
        tokens.clear();
        return;
    }
    tokens[1] = std::to_string(pos + std::stoi(tokens[1]));
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

std::vector <std::string> pack_pat(std::string chrom, int start_site, 
                                       std::string meth_pattern) {
    // Push results into res vector and return
    std::vector <std::string> res;
    res.push_back(chrom);                       // chr
    res.push_back(std::to_string(start_site));  // start site
    res.push_back(meth_pattern);                // meth pattern
    return res;
}

std::vector <std::string> merge_PE(std::vector<std::string> l1, 
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

    if (last_site - start1 > MAX_PE_PAT_LEN) // sanity check: make sure the two reads are not too far apart
        throw std::invalid_argument("invalid pairing. merged read is too long");

    // init merged_pat with missing values
    for (int i = start1; i < last_site; i++)
        merged_pat += UNKNOWN;

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
    l1[2] = merged_pat;
    strip_read(l1);
    return l1;
}

