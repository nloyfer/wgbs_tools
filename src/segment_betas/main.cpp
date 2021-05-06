
#include "segmentor.h"


// g++ -std=c++11 *.cpp -o segment

/**
 *  Main
 *
 */

/*
 * Input Arguments Parsering class
 */
class InputParser{
public:
    InputParser (int &argc, char **argv){
        for (int i=1; i < argc; ++i)
            this->tokens.emplace_back(std::string(argv[i]));
    }
    const std::string& getCmdOption(const std::string &option) const{
        std::vector<std::string>::const_iterator itr;
        itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
        if (itr != this->tokens.end() && ++itr != this->tokens.end()){
            return *itr;
        }
        static const std::string empty_string;
        return empty_string;
    }
    bool cmdOptionExists(const std::string &option) const{
        return std::find(this->tokens.begin(), this->tokens.end(), option)
               != this->tokens.end();
    }
private:
    std::vector <std::string> tokens;
};

Params parse_params(InputParser &input){
    /*
     * fill the params struct
     * start and nr_sites are mandatory arguments
     */
    Params pams;
    pams.start = 0;
    pams.nr_sites = 0;
    pams.max_cpg = 1000;
    pams.max_bp = 0;
    pams.pseudo_count = 1;

    // parse start, nr_sites, max_cpg and pseudo_count

    // start site
    std::string start_str = input.getCmdOption("-s");
    if (!start_str.empty()) {
        pams.start = std::stoul(start_str);
    } else {
        throw "start sites (-s) must be provided\n";
    }

    // number of sites
    std::string nr_sites_str = input.getCmdOption("-n");
    if (!nr_sites_str.empty()) {
        pams.nr_sites = std::stoul(nr_sites_str);
    } else {
        throw "number of sites (-n) must be provided\n";
    }

    // max_cpg
    std::string max_cpg_str = input.getCmdOption("-max_cpg");
    if (!max_cpg_str.empty()) {
        pams.max_cpg = std::stoul(max_cpg_str);
    }

    // max_bp
    std::string max_bp_str = input.getCmdOption("-max_bp");
    if (!max_bp_str.empty()) {
        pams.max_bp = std::stoul(max_bp_str);
    }

    // pseudo count
    std::string psud_str = input.getCmdOption("-ps");
    if (!psud_str.empty()) {
        pams.pseudo_count = std::stof(psud_str);
    }

    return pams;
}

int main(int argc, char *argv[]) {

    InputParser input(argc, argv);
    if (argc < 6){
        std::cerr << "Usage: segment BETA_PATH [BETA_PATH...] -s START -n NR_SITES ";
        std::cerr << " [-m max_cpg] [-ps PSEUDO_COUNT]" << std::endl;
        return -1;
    }

    struct Params pams = parse_params(input);

    // parse beta files paths
    std::vector<std::string> beta_paths;
    for (int i = 1; i < argc; i++) {
        std::string s(argv[i]);
        if ((s.length() >= 6) && !(s.compare(s.length() - 5, 5, ".beta"))) {
            beta_paths.emplace_back(s);
        }
    }

    segmentor segmentor1(pams, beta_paths);
    segmentor1.dp_wrapper();
    return 0;
}

