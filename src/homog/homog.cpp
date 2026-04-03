
#include "homog.h"

Homog::Homog(std::string in_blocks_path, std::vector<float> in_range,
             int in_min_cpgs, bool deb, std::string in_name,
             std::string in_chrom, bool incl, bool in_sort_blocks) {
    min_cpgs = in_min_cpgs;
    blocks_path = in_blocks_path;
    range = in_range;
    chrom = in_chrom;
    debug = deb;
    inclusive = incl;
    sort_blocks = in_sort_blocks;
    name = in_name;
    sname = "[ homog " + name + " ] ";
    nr_bins = range.size() - 1;
}

Homog::~Homog() {
    close_blocks_pipe();
}


/***************************************************
 *                                                 *
 *          Block streaming & flushing             *
 *                                                 *
 **************************************************/

void Homog::open_blocks_pipe() {
    std::string cmd = "cat ";
    if (hasEnding(blocks_path, ".gz")) {
        cmd = "gunzip -c ";
    }
    cmd += blocks_path;
    if (sort_blocks) {
        cmd += " | sort -k4,4n -k5,5n";
        std::cerr << sname + "WARNING: blocks file is not sorted. "
                  << "Sorting on the fly (requires buffering the entire file in memory)." << std::endl;
    }
    blocks_pipe = popen(cmd.c_str(), "r");
    if (blocks_pipe == nullptr) {
        throw std::invalid_argument(sname + "Error: Unable to open blocks path: " + blocks_path);
    }
}

void Homog::close_blocks_pipe() {
    if (blocks_pipe != nullptr) {
        int status = pclose(blocks_pipe);
        blocks_pipe = nullptr;
        if (status != 0) {
            std::cerr << sname + "WARNING: blocks pipe exited with status "
                      << status << " (possible truncated or corrupt file)" << std::endl;
        }
    }
}

void Homog::load_blocks_until(int max_start) {
    /** Load blocks from pipe until next block's startCpG > max_start or EOF.
     *  This ensures all blocks that could overlap a read ending at max_start are loaded. */
    if (pipe_exhausted) return;

    char buf[4096];
    while (fgets(buf, sizeof(buf), blocks_pipe) != nullptr) {
        std::string line(buf);
        // strip trailing newline
        if (!line.empty() && line.back() == '\n') line.pop_back();

        // skip empty lines and comments
        if (line.empty() || (!(line.rfind("#", 0)))) { continue; }

        std::vector<std::string> tokens = line2tokens(line);
        if (tokens.size() < 5) {
            std::cerr << sname + "Invalid blocks file format. ";
            std::cerr << sname + "Should be: chr start end startCpG endCpG\n";
            print_vec(tokens);
            throw std::invalid_argument("Invalid block format");
        }

        // skip header, if exists
        if ((blocks_loaded == 0) && (tokens[0] == "chr")) { continue; }

        // skip other chromosomes, if chrom!=""
        if ((chrom != "") && (tokens[0] != chrom)) { continue; }

        int cur_start = std::stoi(tokens[3]);
        int cur_end = std::stoi(tokens[4]);

        // validate block
        if (cur_end <= cur_start) {
            std::cerr << sname + "Invalid block: " << cur_start << "\t" << cur_end << std::endl;
            throw std::invalid_argument("Invalid block: endCpG <= startCpG");
        } else if (cur_start < 1) {
            throw std::invalid_argument("Invalid block: startCpG < 1");
        }

        // push new block with zeroed counts
        BlockData b;
        b.start = cur_start;
        b.end = cur_end;
        b.counts.assign(nr_bins, 0);
        blocks.push_back(std::move(b));
        blocks_loaded++;

        if (debug && (blocks_loaded >= 2500)) {
            pipe_exhausted = true;
            return;
        }

        // stop loading if this block starts beyond what we need
        if (cur_start > max_start) {
            return;
        }
    }
    // reached EOF
    pipe_exhausted = true;
}

void Homog::flush_block(BlockData &b) {
    /** Output one block's counts as a tab-separated line to stdout */
    for (int j = 0; j < nr_bins; j++) {
        std::cout << b.counts[j];
        if (j < nr_bins - 1)
            std::cout << SEP;
    }
    std::cout << std::endl;
}

void Homog::flush_through(int block_ind) {
    /** Flush all blocks in [flush_offset, block_ind) — they are complete. */
    while (flush_offset < block_ind && !blocks.empty()) {
        flush_block(blocks.front());
        blocks.pop_front();
        flush_offset++;
    }
}

void Homog::flush_remaining() {
    /** Flush everything still in the deque (blocks after the last read). */
    while (!blocks.empty()) {
        flush_block(blocks.front());
        blocks.pop_front();
        flush_offset++;
    }
}


/***************************************************
 *                                                 *
 *             Read processing                     *
 *                                                 *
 **************************************************/

void Homog::update_m2(int block_ind, const std::string &pat, int offset, int length, int count) {

    int nrC = 0;
    int nrT = 0;
    for (int i = offset; i < offset + length; i++) {
        if ((pat[i] == METH) || (pat[i] == HYDROXY))
            nrC++;
        else if (pat[i] == UNMETH)
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
    block_at(block_ind).counts[bin_ind] += count;
}

void Homog::update_block(int block_ind, const std::string &pat, int offset, int length,
                         const std::string &orig_pat, int32_t count) {

    if (inclusive) {
        // inclusive mode: use the full original pattern
        if ((int)orig_pat.length() < min_cpgs) { return; }
        update_m2(block_ind, orig_pat, 0, (int)orig_pat.length(), count);
    } else {
        if (length < min_cpgs) { return; }
        update_m2(block_ind, pat, offset, length, count);
    }
}

int Homog::proc_line(const std::vector <std::string> &tokens) {
    /**
     * Given one line of the form "chr1 1245345 CCT 3", update counts array -
     * Increase in the corresponding cells.
     * Each overlapping block independently clips the original pattern.
     */
    if (tokens.size() < 4) {
        throw std::invalid_argument("Invalid site in input file. too few columns");
    }

    int read_start = std::stoi(tokens[1]);
    const std::string &pattern = tokens[2];
    int count = std::stoi(tokens[3]);
    int read_end = read_start + (int) pattern.length() - 1;
    int pat_len = (int) pattern.length();

    /** read starts later than the last loaded border — finish
        (if pipe is exhausted, no more blocks to load)
                                                    CTCTCT
       |---|   |-----| |--|  ...  |-----| |-|               */
    if (!blocks.empty() && pipe_exhausted &&
            read_start >= blocks.back().end) {
        std::cerr << sname + "Reads reached the end of blocks: " << read_start << " > "
                  << blocks.back().end << ". Breaking.\n";
        return 1;
    }
    /** advance past blocks whose end <= read_start
        (reads are sorted, so no future read can overlap these blocks)
                       CTCTCT
             |-----|                                         */
    while ((cur_block_ind < loaded_end()) && (read_start >= block_end(cur_block_ind))) {
        cur_block_ind++;
    }
    if (cur_block_ind >= loaded_end() && pipe_exhausted)
        return 1;
    /** read ends before current block starts: skip read.
        CTCTCT
               |-----|                                       */
    if (cur_block_ind < loaded_end() && read_end < block_start(cur_block_ind))
        return 0;

    /** For each block that could overlap this read, independently clip.
        Blocks may overlap, so each block gets its own clip of the original pattern.
              CTCTCTCTCTCTCTCT
           |------|                          block A: clip to A
               |------|                      block B: clip to B (independent)
                          |------|           block C: clip to C
                                   |---|     block D: no overlap, stop          */
    for (int bi = cur_block_ind; bi < loaded_end(); bi++) {
        if (block_start(bi) > read_end)
            break;  // no more blocks can overlap

        int overlap_start = std::max(read_start, block_start(bi));
        int overlap_end = std::min(read_start + pat_len, block_end(bi));
        if (overlap_start >= overlap_end)
            continue;

        int offset = overlap_start - read_start;
        int length = overlap_end - overlap_start;
        update_block(bi, pattern, offset, length, pattern, (int32_t) count);
    }
    return 0;
}

void Homog::parse() {

    try {
        open_blocks_pipe();

        int line_i = 0;
        for (std::string line_str; std::getline(std::cin, line_str);) {
            if (line_str.empty()) { continue; } // skip empty lines

            std::vector <std::string> tokens = line2tokens(line_str);
            if (tokens.size() < 4) {
                throw std::invalid_argument("too few columns in file, line " + std::to_string(line_i));
            } else if (!(tokens.empty())) {
                int read_end = std::stoi(tokens[1]) + (int) tokens[2].length() - 1;

                // ensure all blocks that could overlap this read are loaded
                load_blocks_until(read_end);

                // flush blocks that are now complete
                flush_through(cur_block_ind);

                int r = proc_line(tokens);
                if (r == 1) {
                    std::cerr << sname + "breaking" << std::endl;
                    break;
                }
            } else {
                std::cerr << "something went wrong... tokens is empty" << std::endl;
            }
            line_i++;
            if (line_i % 10000000 == 0) {
                std::cerr << sname << line_i / 1000000 << "M" << std::endl;
            }
        }

        // load any remaining blocks from the pipe
        load_blocks_until(INT_MAX);
        // flush everything still in the deque
        flush_remaining();
        // check pipe exit status
        close_blocks_pipe();

        std::cerr << sname + "finished " << line_i << " reads, "
                  << blocks_loaded << " blocks." << std::endl;

    }
    catch (std::exception &e) {
        std::cerr << sname + "failed calculating homog" << std::endl;
        std::cerr << sname << e.what() << std::endl;
        return;
    }
}


/***************************************************
 *                                                 *
 *                      Main                       *
 *                                                 *
 **************************************************/

std::vector<float> parse_range(std::string &range_str) {
    std::vector<float> vect = split_float_by_comma(range_str);

    float tmp = -1;
    for (std::size_t i = 0; i < vect.size(); i++) {
        if (vect[i] <= tmp) {
            std::cerr << "Invalid range - non monotonic: " << range_str << std::endl;
            throw std::invalid_argument("Invalid range");
        }
        if ((vect[i] < 0) || (vect[i] > 1)) {
            std::cerr << "Invalid range - must be in [0,1]: " << vect[i] << std::endl;
            throw std::invalid_argument("Invalid range");
        }
        tmp = vect[i];
    }
    if ((vect[0] > 0) || (vect[vect.size() - 1] < 1)) {
        std::cerr << "Invalid range - must start with 0 and end with 1."  << std::endl;
        throw std::invalid_argument("Invalid range");
    }
    return vect;
}

int main(int argc, char *argv[]) {

    if ((argc < 3)) {
        std::cerr << "Usage: EXEC -r RANGE -b BLOCKS_PATH [-l MIN_LEN] [-n NAME] [-d] [--sort_blocks]\n"
                  << "-l            Minimal sites in read to consider. Default is l=5.\n"
                  << "-r            Ranges of methylation average, in [0,1]. For example: 0,0.2001,0.8,1.\n"
                  << "-b            Blocks file. No header. First 5 columns are:\n"
                  << "              <chr, start, end, startCpG, endCpG>.\n"
                  << "              File may be gzipped. Note sometimes bgzip causes problems.\n"
                  << "-d            Debug mode. Only use the first 2,500 blocks in the blocks file.\n"
                  << "--sort_blocks Sort blocks by CpG index before processing.\n"
                  << std::endl;
        return -1;
    }

    InputParser input(argc, argv);
    std::string name = input.getCmdOption("-n");
    std::string blocks_path = input.getRequiredOption("-b");
    std::string range_str = input.getRequiredOption("-r");
    int min_cpgs = std::stoi(input.getOptionWithDefault("-l", DEFAULT_LEN));
    std::string chrom = input.getCmdOption("--chrom");
    bool debug = input.cmdOptionExists("-d");
    bool inclusive = input.cmdOptionExists("--inclusive");
    bool sort_blocks = input.cmdOptionExists("--sort_blocks");

    try {
        std::vector<float> range = parse_range(range_str);
        Homog(blocks_path, range, min_cpgs, debug, name, chrom, inclusive, sort_blocks).parse();
    }
    catch (std::exception &e) {
        std::cerr << e.what() << std::endl;
        return -1;
    }
    return 0;
}
