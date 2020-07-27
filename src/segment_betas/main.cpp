
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

void parse_params(InputParser &input, std::vector<float> &params){
    /*
     * fill the params vector with:
     * [start, nr_sites, max_size, pseudo_count]
     * start and nr_sites are mandatory arguments
     * max_size anad pseudo_count have default values.
     */
    float max_size = 1000;
    float pseudo_count = 1;
    float start = 0, nr_sites = 0;

    // parse start, nr_sites, max_size and pseudo_count

    // start site
    std::string start_str = input.getCmdOption("-s");
    if (!start_str.empty()) {
        start = std::stoi(start_str);
    } else {
        throw "start sites (-s) must be provided\n";
    }

    // number of sites
    std::string nr_sites_str = input.getCmdOption("-n");
    if (!nr_sites_str.empty()) {
        nr_sites = std::stoi(nr_sites_str);
    } else {
        throw "number of sites (-n) must be provided\n";
    }

    // max_size
    std::string max_sz_str = input.getCmdOption("-m");
    if (!max_sz_str.empty()) {
        max_size = std::stoi(max_sz_str);
    }

    // pseudo count
    std::string psud_str = input.getCmdOption("-ps");
    if (!psud_str.empty()) {
        pseudo_count = std::stoi(psud_str);
    }

    params.emplace_back(start);
    params.emplace_back(nr_sites);
    params.emplace_back(max_size);
    params.emplace_back(pseudo_count);
}

int main(int argc, char *argv[]) {


    std::vector<float> params;

    InputParser input(argc, argv);
    if (argc < 6){
        std::cerr << "Usage: segment BETA_PATH [BETA_PATH...] -s START -n NR_SITES ";
        std::cerr << " [-m MAX_SIZE] [-ps PSEUDO_COUNT]" << std::endl;
        return -1;
    }

    parse_params(input, params);

    // parse beta files paths
    std::vector<std::string> beta_paths;
    for (int i = 1; i < argc; i++) {
        std::string s(argv[i]);
        if ((s.length() >= 6) && !(s.compare(s.length() - 5, 5, ".beta"))) {
            beta_paths.emplace_back(s);
        }
    }
    // print beta files paths
//    std::cerr << "beta files:\n";
//    for (const auto &b: beta_paths) {
//        std::cerr << b << "\n";
//    }

    segmentor segmentor1(params, beta_paths);
    segmentor1.dp_wrapper();
    return 0;
}