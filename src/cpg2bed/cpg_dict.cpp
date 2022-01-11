#include "cpg_dict.h"

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

std::string exec(const char* cmd) {
    /** Execute a command and load output to string */
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("[ add_loci ] popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

/***************************************************************
 *                                                             *
 *               Class methods                                 *
 *                                                             *
 ***************************************************************/

void load_dict_chr(std::string dict_path, std::vector<int> &locs, std::string chrom) {
    //[>* Load the CpG->Loci dict <]
    
    std::string cmd = "tabix " + dict_path + " " + chrom + " | cut -f2";
    std::string ref_data = exec(cmd.c_str());
    if (ref_data.length() == 0) {
        throw std::invalid_argument("[ cpg_dict ] Error: Unable to read reference path: " + dict_path);
    }
    std::stringstream ss(ref_data);

    int i = 0;
    for (std::string line; std::getline(ss, line, '\n');i++) {
        locs[i] = std::stoi(line);
    }
    //std::cerr << chrom << " ) dict is loaded. " << std::to_string(i) << std::endl;
}


void cpg_dict::load_dict(std::string dict_path) {
    
    // init vector of vectors by chromosome
    std::vector <std::vector<int>> v(chroms.size());
    for (int i = 0; i < chroms.size(); i++) {
        v[i].resize(sizes[i], 0);
    }

    // send threads to load chrom sections from CpG.bed.gz:
    std::vector<std::thread> myThreads;
    for (int i = 0; i < chroms.size(); i++) {
        myThreads.push_back(std::thread(load_dict_chr, dict_path, std::ref(v[i]), chroms[i]));
    }
    // join all threads
    for (int i = 0; i < chroms.size(); i++) { myThreads[i].join(); }
    // concatenate vectors from all chromosomes
    for (int i = 0; i < chroms.size(); i++) { loci.insert(loci.end(), v[i].begin(), v[i].end()); }
}


int cpg_dict::load_chroms(std::string chrom_size_path) {
    std::ifstream infile(chrom_size_path);
    if (!(infile.is_open())) { 
        throw std::runtime_error("[ cpg_dict ] Failed to open CpG chrom dict! " + chrom_size_path);
    }

    int t = 0;
    std::vector <std::string> tokens;
    for (std::string line_str; std::getline(infile, line_str);) {
        tokens = line2tokens(line_str);
        chroms.push_back(tokens[0]);
        int csize = std::stoi(tokens[1]);
        t += csize;
        borders.push_back(t);
        sizes.push_back(csize);
    }
    infile.close();
    return t;
}

bool cpg_dict::is_border(std::string chrom, int cpg_pos){
    // return true iif cpg_pos equals to the last CpG in chrom + 1
    // This is useful for the case where we have a region 
    // that ends exactly at the end of the chromosome (convert --sites_file)
    chrom_exists(chrom);
    for (int i = 0; i < chroms.size(); i++) {
        if (chrom == chroms[i]) {
            return (cpg_pos == borders[i]);
        }
    }
    throw std::runtime_error("[ cpg_dict ] Invalid chromosome: " + chrom + "\n");
}

void cpg_dict::chrom_exists(std::string chrom) {
    if (std::find(chroms.begin(), chroms.end(), chrom) != chroms.end()) {
        return;
    }
    throw std::runtime_error("[ cpg_dict ] Invalid chromosome: " + chrom + "\n");
}

std::string cpg_dict::loc2chrom(int loc) {
    int i = 0;
    for (i = 0; i < chroms.size(); i++) {
        if (loc <= borders[i]) {
            return chroms[i];
        }
    }
    // If loc == NR_CPGS + 1, we allow it 
    // (for endCpG, which is non-inclusive)
    if (loc == borders[i - 1] + 1) {
        return chroms[i - 1];
    }
    throw std::runtime_error("[ cpg_dict ] Could not find chromosome for site: " + std::to_string(loc));
}
