

#include "cview.h"


void Cview::strip_read(std::vector <std::string> &tokens) {
    // strip pat (remove trailing dots):
    int pos = strip_pat(tokens[2]);
    if (pos < 0 ) { 
        tokens.clear();
        return;
    }
    tokens[1] = std::to_string(pos + std::stoi(tokens[1]));
}

void Cview::pass_read(std::vector <std::string> &tokens) {
    /** print a vector to stdout, tab separated */
    if (strip) {
        strip_read(tokens);
        if (tokens.empty()) { return; }
    }
    if (min_cpgs > tokens[2].length()) { return; }
    output_vec(tokens);
}

void Cview::output_vec(std::vector <std::string> &tokens) {
    const auto separator = "\t";
    auto sep = "";
    for (const auto &j: tokens) {
        std::cout << sep << j;
        sep = separator;
    }
    std::cout << std::endl;
}

std::string load_blocks(std::string blocks_path, std::string sites) {
    // Load the CpG blocks

    std::string block_data = "";
    if (blocks_path == "") {
        block_data = sites;
    }
    else {
        std::string cmd = "cat ";
        if ((blocks_path.length() > 3) && (blocks_path.substr(blocks_path.length() - 3) == ".gz")) {
            cmd = "gunzip -c ";
        }
        cmd += blocks_path + " " + " | grep -v ^# | cut -f4-5 | sort -k1,1n";
        block_data = exec(cmd.c_str());
        if (block_data.length() == 0) {
            throw std::invalid_argument("[ cview ] Error: Unable to read blocks file: " + blocks_path);
        }
    }
    return block_data;
}

int Cview::read_blocks() {
    // Load blocks to string
    std::string block_data = load_blocks(blocks_path, sites);
    std::stringstream ss(block_data);
    //
    //Iterate lines
    std::vector <std::string> tokens;
    int cur_start = 0, cur_end = 0, bi = 0;
    for (std::string line_str; std::getline(ss, line_str);) {

        // skip empty lines
        if (line_str.empty() || (!(line_str.rfind("#", 0)))) { continue; }
        tokens = line2tokens(line_str);

        if (tokens.size() != 2) {
            std::cerr << "Invalid blocks format. Should be: startCpG\tendCpG\n";
            print_vec(tokens);
            abort();
        }

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

        // if block is duplicated, continue
        if ((bi > 0) && (borders[bi - 1].start == cur_start) && (borders[bi - 1].end == cur_end)) {
            borders[bi - 1].count++;
            continue;
        }
        
        // block isn't dup:
        Block b = {cur_start, cur_end, 1};
        borders.push_back(b);

        // make sure there's no overlap or non-monotonic blocks
        if ( (bi > 0) && (borders[bi].start < borders[bi - 1].end) ) {
            std::cerr << "Invalid block start: " << cur_start << std::endl;
            std::cerr << "Make sure blocks are sorted by CpG-Index and monotonic (blockEnd > blockStart).\n";
            throw std::invalid_argument("Invalid blocks");
        }

        if (debug && (bi >= 2500)) {
            break;
        }
        bi++;
    }
    if (borders.size() == 0) {
        throw std::invalid_argument("Error while loading blocks. 0 blocks found.\n");
     }
    nr_blocks = borders.size();
    if (verbose)
        std::cerr << nr_blocks << " blocks are done" << std::endl;
    return 0;
}

int Cview::proc_line(std::vector <std::string> tokens) {
    /**
     * Given one line of the form "chr1 1245345 CCT 3", 
     * print it iff it overlaps the input blocks.
     * Apply the filters --strict, --strip and --min_cpg if requested
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
    if (read_start >= borders[borders.size() - 1].end) {
        if (verbose) {
            std::cerr << "Reads reached the end of blocks: " << read_start << " > "
                      << borders[borders.size() - 1].end << ". Breaking.\n";
        }
        return 1;
    }
    /** read starts after current block ends - move on to next relevant block
                       CTCTCT
             |-----|                                                         */
    while ((read_start >= borders[cur_block_ind].end) && (cur_block_ind < nr_blocks)) {
        cur_block_ind++;
    }
    // todo: dump block line on the fly
    
    /** read ends before current block starts: skip read.
        CTCTCT
               |-----|                                                      */
    if (read_end < borders[cur_block_ind].start)
        return 0;

    if (!strict) {
        pass_read(tokens);
        return 0;
    }
    /** read starts before current block, 
        but continues to the block (and possibly beyond) - clip beginning of read
            CTCTCT
               |-----|                                          */
    if  (read_start < borders[cur_block_ind].start) {
        pattern = pattern.substr(borders[cur_block_ind].start - read_start);
        read_start = borders[cur_block_ind].start;
    }
    /** If we reached here, current read starts in cur_block_ind (if strict==true)
        Or it starts before cur_block_ind but ends after cur_block start (if strict==false)
              CTCTCT..CCTTTCTCT
           |-----|   |--|    |---|
       tmp_bi index is running for the current read, from cur_block_ind until the read ends. */
    int tmp_bi = cur_block_ind;
    while ((!pattern.empty()) && (tmp_bi < nr_blocks)) {
        if ((read_start >= borders[tmp_bi].start) && (read_start < borders[tmp_bi].end)) {
            int head_size = std::min(borders[tmp_bi].end - read_start, (int) pattern.length());
            std::vector<std::string> partial_read(tokens);
            partial_read[1] = std::to_string(read_start);
            partial_read[2] = pattern.substr(0, head_size);
            pass_read(partial_read);
            pattern = pattern.substr(head_size);
            read_start += head_size;
            tmp_bi++;
        } else if (read_end < borders[tmp_bi].start) {
            break;
        } else {
            pattern = pattern.substr(borders[tmp_bi].start - read_start);
            read_start = borders[tmp_bi].start;
        }
    }
    return 0;
}


void Cview::parse() {

    try {
        // load blocks:
        read_blocks();

        // parse reads:
        int line_i = 0;
        for (std::string line; std::getline(std::cin, line);) {
            if (line.empty()) { continue; } // skip empty lines
            if (proc_line(line2tokens(line))) { break; }
            line_i++;
        }
        if (verbose) fprintf(stderr, "finished %d reads\n", line_i);
    }
    catch (std::exception &e) {
        std::cerr << "failed in cview" << std::endl;
        std::cerr << e.what() << std::endl;
        return;
    }
}



int main(int argc, char *argv[]) {

    InputParser input(argc, argv);

    int min_cpgs = std::stoi(input.getOptionWithDefault("--min_cpgs", "1"));
    std::string blocks_path = input.getCmdOption("--blocks_path");
    std::string sites = input.getCmdOption("--sites");

    bool debug = input.cmdOptionExists("-d");
    bool verbose = input.cmdOptionExists("-v");
    bool strict = input.cmdOptionExists("--strict");
    bool strip = input.cmdOptionExists("--strip");

    try {
        Cview(blocks_path, sites, strict, strip, min_cpgs, debug, verbose).parse();
    }
    catch (std::exception &e) {
        std::cerr << e.what() << std::endl;
        return -1;
    }
    return 0;
}
