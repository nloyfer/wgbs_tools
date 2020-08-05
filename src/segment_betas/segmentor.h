//
// Created by nloyfer on 11/2/18.
//

#ifndef CSEGMENT_SEGMENTOR_H
#define CSEGMENT_SEGMENTOR_H

#include <iostream>

#include <iostream>
#include <ctime>
//#include <string>
#include <vector>
#include <algorithm>    // std::sort
#include <fstream>
#include <regex>
#include <unordered_map>
#include <cmath>
#include <math.h>


struct Params
{
    unsigned int start;
    unsigned int nr_sites;
    float pseudo_count;
    unsigned int max_size;
};


class segmentor {
    std::vector<std::string> beta_paths;
    Params params;
    uint start = 0;
    uint nr_sites = 0;
    float pseudo_count = 0;
    int max_size = 0;
    double *mem = nullptr;


    void read_beta_file(const char *beta_path, float *data);
    void cost_memoization(std::vector<float*> &all_data);
    std::vector<int> traceback(const int *T);
public:
    segmentor(Params &iparams, std::vector<std::string> &input_beta_paths):
            beta_paths(input_beta_paths), params(iparams) {}
    void dp_wrapper();
    void dp(std::vector<float*> &all_data);

};


#endif //CSEGMENT_SEGMENTOR_H
