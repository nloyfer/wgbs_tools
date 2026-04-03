#include "cview.h"

bool does_read_contain_dots(std::vector <std::string> &tokens) {
    // return true iff the read contains dots (unknown methylation).
    return (tokens[2].find('.') != std::string::npos);
}

void Cview::pass_read(std::vector <std::string> &tokens) {
    /** print a vector to stdout, tab separated */
    if (strip) {
        strip_read(tokens);
        if (tokens.empty()) { return; }
    }
    if (min_cpgs > tokens[2].length()) { return; }
    if (no_gaps && does_read_contain_dots(tokens)) { return; }
    output_vec(tokens);
}

std::string block_path_to_string_data(const std::string &blocks_path) {
    std::string cmd = hasEnding(blocks_path, ".gz") ? "gunzip -c " : "cat ";
    cmd += blocks_path + " | cut -f4-5 | sort -k1,1n";
    FILE *pipe = popen(cmd.c_str(), "r");
    if (pipe == nullptr) {
        throw std::invalid_argument("[ cview ] Error: Unable to open blocks file: " + blocks_path);
    }
    std::string block_data;
    char buf[4096];
    while (fgets(buf, sizeof(buf), pipe) != nullptr) {
        block_data += buf;
    }
    int status = pclose(pipe);
    if (block_data.empty()) {
        throw std::invalid_argument("[ cview ] Error: Unable to read blocks file: " + blocks_path);
    }
    if (status != 0) {
        std::cerr << "[ cview ] WARNING: blocks pipe exited with status " << status << std::endl;
    }
    return block_data;
}

int read_blocks(const std::string &block_data, std::vector<Block> &borders) {
    // Load blocks to string
    std::stringstream ss(block_data);

    //Iterate lines
    std::vector <std::string> tokens;
    int cur_start = 0, cur_end = 0, bi = 0;
    for (std::string line_str; std::getline(ss, line_str);) {

        // skip empty lines
        if (line_str.empty() || (!(line_str.rfind("#", 0)))) { continue; }
        tokens = line2tokens(line_str);

        if (tokens.size() != 2) {
            std::cerr << "Invalid blocks format. Should be: startCpG\tendCpG\n";
            std::cerr << tokens.size() << std::endl;
            print_vec(tokens);
            throw std::invalid_argument("Invalid block file");
        }
        // skip header line if exists
        if ((bi == 0) && (!is_number(tokens[0]))) { continue; }

        // skip blocks not covering any CpGs
        if ((tokens[0] == "NA") || (tokens[0] == "NaN")) { continue; }
        cur_start = std::stoi(tokens[0]);
        cur_end = std::stoi(tokens[1]);

        // If block is invalid, abort:
        if ((cur_end <= cur_start)) {
            std::cerr << "Invalid block: " << cur_start << "\t" << cur_end << std::endl;
            throw std::invalid_argument("Invalid block: endCpG <= startCpG");
        } else if (cur_start < 1) {
            throw std::invalid_argument("Invalid block: startCpG < 1");
        }

        Block b = {cur_start, cur_end};
        borders.push_back(b);

        bi++;
    }
    if (borders.empty()) {
        throw std::invalid_argument("Error while loading blocks. 0 blocks found.\n");
     }
    return borders.size();
}

int Cview::proc_line(std::vector <std::string> tokens) {
    /**
     * Given one line of the form "chr1 1245345 CCT 3", 
     * print it iff it overlaps the input blocks.
     * Apply the filters --strict, --strip, --no_gaps, and --min_cpg if requested
     */
    if (tokens.size() < 4) {
        throw std::invalid_argument("Invalid site in input file. too few columns");
    }

    int read_start = std::stoi(tokens[1]);
    std::string pattern = tokens[2];
    int read_end = read_start + (int) pattern.length() - 1;

    /** read starts later than the last border - finish
                                                CTCTCT
       |---|   |-----| |--|  ...  |-----| |-|           */
    if (read_start >= borders.at(borders.size() - 1).end) {
        if (verbose) {
            std::cerr << "Reads reached the end of blocks: " << read_start << " > "
                      << borders.at(borders.size() - 1).end << ". Breaking.\n";
        }
        return 1;
    }
    /** read starts after current block ends - move on to next relevant block
                       CTCTCT
             |-----|                                                         */
    while ((cur_block_ind < nr_blocks) && (read_start >= borders.at(cur_block_ind).end)) {
        cur_block_ind++;
    }
    if (cur_block_ind >= nr_blocks)
        return 1;

    /** read ends before current block starts: skip read.
        CTCTCT
               |-----|                                                      */
    if (read_end < borders.at(cur_block_ind).start)
        return 0;

    if (!strict) {
        pass_read(tokens);
        return 0;
    }
    /** read starts before current block, 
        but continues to the block (and possibly beyond) - clip beginning of read
            CTCTCT
               |-----|                                          */
    if (read_start < borders.at(cur_block_ind).start) {
        pattern = pattern.substr(borders.at(cur_block_ind).start - read_start);
        read_start = borders.at(cur_block_ind).start;
    }
    /** If we reached here, then strict==true and current read starts within cur_block_ind
              CTCTCT..CCTTTCTCT
           |-----|   |--|    |---|
              ^
       tmp_bi index is running for the current read, from cur_block_ind until the read ends. */
    int tmp_bi = cur_block_ind;
    while ((!pattern.empty()) && (tmp_bi < nr_blocks)) {
        /* case current tmp_bi has overlap with the read
                     .CCTTTCTCT
           |-----|   |--|    |---|
                      ^
        */
        if ((read_start >= borders.at(tmp_bi).start) && (read_start < borders.at(tmp_bi).end)) {
            int head_size = std::min(borders.at(tmp_bi).end - read_start, (int) pattern.length());
            std::vector<std::string> partial_read(tokens);
            partial_read[1] = std::to_string(read_start);
            partial_read[2] = pattern.substr(0, head_size);
            pass_read(partial_read);
            pattern = pattern.substr(head_size);
            read_start += head_size;
            tmp_bi++;
        } else if (read_end < borders.at(tmp_bi).start) {
            break;
        } else {
            pattern = pattern.substr(borders.at(tmp_bi).start - read_start);
            read_start = borders.at(tmp_bi).start;
        }
    }
    return 0;
}

void Cview::parse() {

    try {
        // load blocks:
        std::string block_data = (blocks_path == "") ? sites : block_path_to_string_data(blocks_path);
        nr_blocks = read_blocks(block_data, borders);

        // parse reads:
        int line_i = 0;
        for (std::string line; std::getline(std::cin, line);) {
            if (line.empty()) { continue; } // skip empty lines
            if (proc_line(line2tokens(line))) { break; }
            line_i++;
        }
        if (verbose) fprintf(stderr, "finished %s reads\n", addCommas(line_i).c_str());
    }
    catch (std::exception &e) {
        std::cerr << "failed in cview" << std::endl;
        std::cerr << e.what() << std::endl;
        return;
    }
}

int main(int argc, char *argv[]) {

    // input arguments
    InputParser input(argc, argv);
    int min_cpgs = std::stoi(input.getOptionWithDefault("--min_cpgs", "1"));
    std::string blocks_path = input.getCmdOption("--blocks_path");
    std::string sites = input.getCmdOption("--sites");
    bool verbose = input.cmdOptionExists("-v");
    bool strict = input.cmdOptionExists("--strict");
    bool strip = input.cmdOptionExists("--strip");
    bool no_gaps = input.cmdOptionExists("--no_gaps");

    try {
        Cview(blocks_path, sites, strict, strip, no_gaps, min_cpgs, verbose).parse();
    }
    catch (std::exception &e) {
        std::cerr << e.what() << std::endl;
        return -1;
    }
    return 0;
}
