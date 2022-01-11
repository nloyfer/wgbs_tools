//
// Created by nloyfer on 22/12/21
//

#include "add_loci.h"



void cexception(std::string msg, int line_i) {
    std::string e = "[wt add_loci] line " + std::to_string(line_i) + ": " + msg;
    throw std::runtime_error(e);
}


/***************************************************************
 *                                                             *
 *               Parse input file                              *
 *                                                             *
 ***************************************************************/


void add_loci_main(std::string dpath, std::string cpath) {

    cpg_dict cd(dpath, cpath);

    // parse input site file
    std::string chrom1 = "";
    std::string chrom2 = "";
    int start = 0, end = 0, startCpG = 0, endCpG = 0, line_i = 0;
    std::vector <std::string> tokens;
    for (std::string line_str; std::getline(std::cin, line_str); line_i++) {

        tokens = line2tokens(line_str);
        int startCpG = std::stoi(tokens[0]);
        int endCpG = (tokens.size() == 1) ? startCpG + 1 : std::stoi(tokens[1]);

        // validate input
        if (endCpG < startCpG) { cexception("endCpG < startCpG", line_i); }
        if (startCpG < 1) { cexception("startCpG < 1", line_i); }
        if (endCpG < 1) { cexception("endCpG < 1", line_i); }
        
        chrom1 = cd.loc2chrom(startCpG);
        // validate chromosome
        chrom2 = cd.loc2chrom(endCpG);
        if (chrom1 != chrom2) { 
            if (!(cd.is_border(chrom1, endCpG - 1))) {
                cexception("Cross chromosomes", line_i); 
            }
        }
        
        start = cd.loci[startCpG - 1];
        end = (endCpG == startCpG) ? start + 2 : cd.loci[endCpG - 2] + 1;
        std::string TAB = "\t";
        std::cout << chrom1 << TAB << start << TAB << end << TAB << startCpG << TAB << endCpG << std::endl;

    }
}


/***************************************************************
 *                                                             *
 *                     Main                                    *
 *                                                             *
 ***************************************************************/


int main(int argc, char **argv) {
    try {
        if (argc != 3) {
            throw std::invalid_argument("Usage: add_loci CPG_DICT CHR_DICT");
        }
        add_loci_main(argv[1], argv[2]);

    }
    catch (std::exception &e) {
        std::cerr << "[ add_loci ] Failed! exception:" << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    }
    return 0;
}
