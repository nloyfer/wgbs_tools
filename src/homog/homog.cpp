
#include "homog.h"


std::vector <std::string> line2tokens(std::string &line) {
    /** Break string line to tokens, return it as a vector of strings */
    std::vector <std::string> result;
    std::string cell;
    std::stringstream lineStream(line);
    while (getline(lineStream, cell, '\t'))
        result.push_back(cell);
    return result;
}

void print_vec(std::vector <std::string> &vec) {
    /** print a vector to stderr, tab separated */
    for (auto &j: vec)
        std::cerr << j << "\t";
    std::cerr << std::endl;
}

void print_vec(std::vector <int> &vec) {
    /** print a vector to stderr, tab separated */
    for (auto &j: vec)
        std::cerr << j << "\t";
    std::cerr << std::endl;
}


int32_t *Homog::init_array(int len) {
    int *arr = new int32_t[nr_blocks * len];
    std::fill_n(arr, nr_blocks * len, 0);
    return arr;
}

std::string exec(const char* cmd) {
    /** Execute a command and load output to string */
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("[ homog ] popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

// Constructor
Homog::Homog(std::string in_blocks_path, std::vector<float> in_range,
             int in_min_cpgs, bool deb, std::string in_name, std::string in_chrom) {
    min_cpgs = in_min_cpgs;
    blocks_path = in_blocks_path;
    range = in_range;
    chrom = in_chrom;
    debug = deb;
    name = in_name;
    sname = "[ homog " + name + " ] ";

    nr_bins = range.size() - 1;

    // load blocks file
    int r = read_blocks();
    nr_blocks = borders_starts.size();

    // Init arrays to zeros
    counts = init_array(nr_bins);
}

Homog::~Homog() {
    //delete[] counts;
}


int Homog::blocks_helper(std::istream &instream) {
    //Iterate lines
    std::vector <std::string> tokens;
    std::string line;
    int cur_start = 0, cur_end = 0;
    int bi = 0;
    int abc = 0;
    while (std::getline(instream, line)) {

        // skip empty lines and comments
        if (line.empty() || (!(line.rfind("#", 0)))) { continue; }

        tokens = line2tokens(line);
        if (tokens.size() < 5) {
            std::cerr << sname + "Invalid blocks file format. ";
            std::cerr << sname + "Should be: chr start end startCpG endCpG\n";
            print_vec(tokens);
            throw std::invalid_argument("Invalid block format");
        }

        // skip header, if exists
        if ((bi == 0) && (tokens[0] == "chr")) { continue; }

        // skip other chromosoms, if chrom!=""
        if ((chrom != "") && (tokens[0] != chrom)) { continue; }

        cur_start = std::stoi(tokens[3]);
        cur_end = std::stoi(tokens[4]);

        // If block is invalid, abort:
        if ((cur_end <= cur_start)) {
            std::cerr << sname + "Invalid block: " << cur_start << "\t" << cur_end << std::endl;
            throw std::invalid_argument("Invalid block: endCpG <= startCpG");
        } else if (cur_start < 1) {
            throw std::invalid_argument("Invalid block: startCpG < 1");
        }

        // skip borders with <min_cpgs CpGs
        //if ((cur_end - cur_start < min_cpgs)) {
            //continue;
        //}

        // if block is duplicated, continue
        if ((!borders_starts.empty()) &&    // Can't be dup if it's first 
                (borders_starts.back() == cur_start) &&
                (borders_ends.back() == cur_end)) {
            // only update the count and continue:
            borders_counts[bi - 1]++;
            continue;
        }
        
        // block isn't dup:
        borders_counts.push_back(1);
        borders_starts.push_back(cur_start);
        borders_ends.push_back(cur_end);
        coords.push_back(tokens[0] + "\t" + tokens[1] + "\t" + tokens[2]);
        

        // make sure there's no overlap or non-monotonic blocks
        if ( (bi > 0) && (borders_starts[bi] < borders_ends[bi - 1]) ) {
            std::cerr << sname + "Invalid block start: " << cur_start << std::endl;
            std::cerr << sname + "Make sure blocks are sorted by CpG-Index and monotonic (blockEnd > blockStart).\n";
            throw std::invalid_argument("Invalid blocks");
        }

        if (debug && (bi >= 2500)) {
            break;
        }
        bi++;
    }
    return 0;
}

int Homog::read_blocks() {
    /**
     * Load blocks gzipped file into vector<int> borders_starts, borders_ends.
     */


    std::cerr << sname + "loading blocks..." << std::endl;

    if (hasEnding(blocks_path, ".gz")) {
        std::string cmd = "gunzip -c " + blocks_path;
    else {
        std::string cmd = "cat " + blocks_path;
    }
    std::string blocks_data = exec(cmd.c_str());
    if (blocks_data.length() == 0) {
        throw std::invalid_argument("[ homog ] Error: Unable to blocks path:" + blocks_path);
    }
    std::stringstream ss(blocks_data);

    blocks_helper(ss);

    if (borders_starts.size() == 0) {
        std::cerr << sname + "Error while loading blocks. 0 blocks found.\n";
        //std::cerr << "Error while loading blocks. 0 blocks found.\nMaybe all blocks are too short.\n";
        throw std::invalid_argument("");
     }
    std::cerr << sname + "loaded " << borders_starts.size() << " blocks. last block: " ;
    std::cerr << borders_starts[borders_starts.size() - 1] << "-" << borders_ends[borders_ends.size() - 1] << std::endl;
    return 0;
}


void Homog::dump(int32_t *data, int width) {
    /**
     */
    for (int i = 0; i < nr_blocks; i++) {
        for (int d = 0; d < borders_counts[i]; d++) {
            for (int j = 0; j < width; j++) {
                std::cout << data[i * width + j];
                if (j < width - 1)
                    std::cout << SEP;
            }
            std::cout << std::endl;
        }
    }
}


void Homog::update_m2(int block_ind, std::string pat, int count) {

    int nrC = 0;
    int nrT = 0;
    for (int i = 0; i < pat.length(); i++) {
        if (pat[i] == 'C')
            nrC++;
        else if (pat[i] == 'T')
            nrT++;
    }

    if (nrC + nrT < min_cpgs) {
        return;
    }

    float meth = (float) nrC / (float) (nrC + nrT);
    if (meth < range[0])
        return;

    int bin_ind = 0;
    for (bin_ind = 0; bin_ind < range.size() - 1; bin_ind++) {
        if ((meth >= range[bin_ind]) && (meth < range[bin_ind + 1])) {
            break;
        }
    }
    if (bin_ind == range.size() - 1) {
        bin_ind--;
    }
//    std::cout << "block:" << block_ind << ", pat: " <<  pat << ", meth: " << meth << ", bin_ind:" <<  bin_ind << std::endl;
    counts[block_ind * nr_bins + bin_ind] += count;
}

void Homog::update_block(int block_ind, std::string pat, int32_t count) {

//    std::cerr << block_ind << ") updating: " << pat << std::endl;

    int len = pat.length();
    // skip reads with less then min_cpgs sites:
    if (len < min_cpgs) {
        return;
    }

    update_m2(block_ind, pat, count);
}

int Homog::proc_line(std::vector <std::string> tokens) {
    /**
     * Given one line of the form "chr1 1245345 CCT 3", update counts array -
     * Increase in the corresponding cells.
     */
    if (tokens.size() < 4) {
        throw std::invalid_argument("Invalid site in input file. too few columns");
    }

    int read_start = std::stoi(tokens[1]);
    std::string pattern = tokens[2];
    int count = std::stoi(tokens[3]);
    int read_end = read_start + (int) pattern.length() - 1;

//    print_vec(tokens);


    // read starts later than the last border - finish
    //                                          CTCTCT
    // |---|   |-----| |--|  ...  |-----| |-|
    if (read_start >= borders_ends[borders_ends.size() - 1]) {
        std::cerr << sname + "Reads reached the end of blocks: " << read_start << " > "
                  << borders_ends[borders_ends.size() - 1] << ". Breaking.\n";
        return 1;
    }
    // read starts after current block ends - move on to next relevant block
    //                   CTCTCT
    //         |-----|
    while ((read_start >= borders_ends[cur_block_ind]) && (cur_block_ind < nr_blocks)) {
        cur_block_ind++;
    }
    // todo: dump block line on the fly
    // read ends before current block starts: skip read.
    //  CTCTCT
    //         |-----|
    if (read_end < borders_starts[cur_block_ind])
        return 0;
    // read starts before current block, but continues to the block (and possibly beyond) - clip beginning of read
    //      CTCTCT
    //         |-----|
    if (read_start < borders_starts[cur_block_ind]) {
        pattern = pattern.substr(borders_starts[cur_block_ind] - read_start);
        read_start = borders_starts[cur_block_ind];
    }
    // If we reached here, current reads starts in cur_block_ind
    //        CTCTCT..CCTTTCTCT
    //     |-----|   |--|    |---|
    // tmp_block_ind index is running for the current read, from cur_block_ind until the read ends.
    int tmp_block_ind = cur_block_ind;
    while ((!pattern.empty()) && (tmp_block_ind < nr_blocks)) {
        if ((read_start >= borders_starts[tmp_block_ind]) && (read_start < borders_ends[tmp_block_ind])) {
            int head_size = std::min(borders_ends[tmp_block_ind] - read_start, (int) pattern.length());
            update_block(tmp_block_ind, pattern.substr(0, head_size), (int32_t) count);
            pattern = pattern.substr(head_size);
            read_start += head_size;
            tmp_block_ind++;
        } else if (read_end < borders_starts[tmp_block_ind]) {
            break;
        } else {
            pattern = pattern.substr(borders_starts[tmp_block_ind] - read_start);
            read_start = borders_starts[tmp_block_ind];
        }
    }
    return 0;
}

void Homog::debug_print(int ind, std::vector <std::string> tokens) {
    std::cerr << ind << ") " << borders_starts[ind] << "-" << borders_ends[ind] << std::endl;
    print_vec(tokens);
}


bool hasEnding(std::string const &fullString, std::string const &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

void Homog::parse() {

    try {
        int line_i = 0;
        for (std::string line_str; std::getline(std::cin, line_str);) {
            if (line_str.empty()) { continue; } // skip empty lines

            std::vector <std::string> tokens = line2tokens(line_str);
//            std::cerr << line_i << std::endl;
//	    print_vec(tokens); // TODO: del
            if (tokens.size() < 4) {
                throw std::invalid_argument("too few columns in file, line " + std::to_string(line_i));
            } else if (!(tokens.empty())) {
                int r = proc_line(tokens);
                if (r == 1) {
                    std::cerr << sname + "breaking" << std::endl;
                    break;
                }
            } else {
                std::cerr << "something went wrong... tokens is empty" << std::endl;
            }
            line_i++;
//            std::cerr << line_i << std::endl;
            if (line_i % 10000000 == 0) {
                std::cerr << sname << line_i / 1000000 << "M" << std::endl;
            }
        }

        // dump reads lengths file:
        dump(counts, nr_bins);
        std::cerr << sname + "finished " << line_i << " reads." << std::endl;

    }
    catch (std::exception &e) {
        std::cerr << sname + "failed calculating homog" << std::endl;
        std::cerr << sname << e.what() << std::endl;
        return;
    }
}


