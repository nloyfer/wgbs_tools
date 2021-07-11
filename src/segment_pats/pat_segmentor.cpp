//
// Created by jrosensk on 10/08/2020.
//

#include "pat_segmentor.h"
#include <vector>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <iostream>
#include <deque>

#define NR_SITES 28217448


//  g++ -std=c++11 stdin2beta.cpp -o stdin2beta
struct Params
{
    unsigned int start;
    unsigned int nr_sites;
    unsigned int max_cpg;
    unsigned int min_cpg_sites_per_read;
    unsigned int min_reads_per_window;
    float read_homog_cutoff;
    unsigned int min_bimodal_window_len;
    float bimodal_win_threshold;
};

struct ScoreObject {
    int count;
    int site;
    std::string pattern;
    int end_site;
};

struct Block {
    int start = -1;
    int end = -1;
    bool is_empty = true;
};

struct UXM {
    int num_u = 0;
    int num_m = 0;
    int num_x = 0;
};


std::vector<std::string> line2tokens(std::string &line) {
    /**
     * Break string line to tokens, return it as a vector of strings
     */
    std::vector<std::string> result;
    std::string cell;
    std::stringstream lineStream(line);
    while(getline(lineStream, cell, '\t'))
        result.push_back(cell);
    return result;
}


class PatSegmentor {
    int nr_sites = NR_SITES;
    int start_index = 1;
    int first_site = -1;
    int end_index = 20000;
    int max_window_size = 501;
    int min_cpg_sites_per_read = 5;
    int min_reads_per_window = 10;
    float read_homog_cutoff = 0.8;
    int min_bimodal_window_len = 80;
    float bimodal_win_threshold = 0.25;
    double** mem = nullptr;
    Params params;

public:
    PatSegmentor(Params params);
    ~PatSegmentor();
    void memoize_costs(std::vector<std::string> data_paths);
    void process_cur_window(std::deque<ScoreObject> prev_lines);
    void process_one_step(std::deque<ScoreObject> prev_lines);
    void scan_for_bimodal(std::string data_path);
    void scan_for_bimodal_2(std::vector<std::string> data_paths);
    void dp();
    std::vector<int> traceback(const int *T);
};


PatSegmentor::~PatSegmentor() {
}

PatSegmentor::PatSegmentor(Params in_params) {

    start_index = in_params.start;
    end_index = in_params.start + in_params.nr_sites - 1;
    nr_sites = in_params.nr_sites;
    max_window_size = in_params.max_cpg;
    first_site = -1;
    params = in_params;
    min_cpg_sites_per_read = in_params.min_cpg_sites_per_read;
    read_homog_cutoff = in_params.read_homog_cutoff;
    min_bimodal_window_len = in_params.min_bimodal_window_len;
    bimodal_win_threshold = in_params.bimodal_win_threshold;
    min_reads_per_window = in_params.min_reads_per_window;
}

struct ProportionCount {
    int min_size;
    int total;
    double min_proportion;
};

double process_score_object(int lower_bound, int upper_bound, ScoreObject scoreObject){
    int cur_site = std::max(lower_bound, scoreObject.site);
    int num_c = 0;
    int num_t = 0;
    while (cur_site < upper_bound and (cur_site - scoreObject.site) < scoreObject.pattern.size()){
        char cur_letter = scoreObject.pattern[cur_site - scoreObject.site];
        if (cur_letter == 'C'){
            num_c += 1;
        } else if (cur_letter == 'T'){
            num_t += 1;
        }
        cur_site++;
    }
    int total = num_c + num_t;
    int min_state = std::min(num_c, num_t);
    double min_proportion = total == 0 ? 0 : min_state / (double) total;
    return min_proportion;
}


void PatSegmentor::process_cur_window(std::deque<ScoreObject> prev_lines) {
    int cur_line_index = 0;
    int start_site = prev_lines[cur_line_index].site;
    int cur_end_index = prev_lines.size() - 1;
    int available_size = prev_lines[cur_end_index].site  - prev_lines[cur_line_index].end_site;

    while(available_size >= max_window_size){
        int i = start_site - first_site;
        for (int j = 1; j <= max_window_size; j++){
            double min_proportion_sum = 0;
            double total = 0;
            int cpg_j = start_site + j;
            for (ScoreObject scoreObject : prev_lines){
                if (scoreObject.end_site < start_site){
                    continue;
                }
                if (scoreObject.site < cpg_j){
                    double min_proportion = process_score_object(start_site, cpg_j, scoreObject);
                    min_proportion_sum += min_proportion * scoreObject.count;
                    total += 1 * scoreObject.count;
                } else {
                    break;
                }
            }
            double to_add = total > 0 ? min_proportion_sum / total : 0;
            mem[i][j] = mem[i][j] == -1 ? to_add : mem[i][j] + to_add;
        }
        int new_start_site = start_site;
        while (new_start_site == start_site and cur_line_index < prev_lines.size()){
            cur_line_index++;
            new_start_site = prev_lines[cur_line_index].site;
        }
        if (cur_line_index >= prev_lines.size()){
            break;
        }
        start_site = new_start_site;
        available_size = prev_lines[cur_end_index].site  - prev_lines[cur_line_index].end_site;
    }
}

void PatSegmentor::process_one_step(std::deque<ScoreObject> prev_lines) {
    int cur_line_index = 0;
    int start_site = prev_lines[cur_line_index].site;
    int end_index = prev_lines.size() - 1;
    int available_size = prev_lines[end_index].site  - prev_lines[cur_line_index].end_site;;
    int cur_size = prev_lines[end_index].site  - prev_lines[cur_line_index].end_site;;

    while(cur_size >= available_size){
        int i = start_site - first_site;
        if (i >= nr_sites){
            break;
        }
        for (int j = 1; j <= max_window_size; j++){
            double min_proportion_sum = 0;
            int total_Count = 0;
            int cpg_j = start_site + j;
            for (ScoreObject scoreObject : prev_lines){
                if (scoreObject.end_site < start_site){
                    continue;
                }
                if (scoreObject.site < cpg_j){
                    double min_proportion = process_score_object(start_site, cpg_j, scoreObject);
                    min_proportion_sum += min_proportion * scoreObject.count;
                    int x = scoreObject.count;
                    total_Count = total_Count + x;
                } else {
                    break;
                }
            }
            double to_add = total_Count > 0 ? (float)(min_proportion_sum) / (float)(total_Count) : 0;//total_Count

            if (mem[i][j] == -1){
                mem[i][j] = to_add;
            } else {
                mem[i][j] = mem[i][j] + to_add;
            }
            mem[i][j] = mem[i][j] == -1 ? to_add : mem[i][j] + to_add;
        }
        int new_start_site = start_site;
        while (new_start_site == start_site and cur_line_index < prev_lines.size()){
            cur_line_index++;
            new_start_site = prev_lines[cur_line_index].site;
        }
        if (cur_line_index >= prev_lines.size()){
            break;
        }
        start_site = new_start_site;
        cur_size = prev_lines[end_index].site  - prev_lines[cur_line_index].end_site;
    }
}

ScoreObject get_score_object_from_tokens(std::vector<std::string> tokens){
    if (tokens.size() < 4) {
        throw std::invalid_argument( "Invalid site in input file. too few columns" );
    }

    int site = std::stoi(tokens[1]);
    std::string pattern = tokens[2];
    int count = std::stoi(tokens[3]);
    auto read_len = (int) pattern.length();

    ScoreObject scoreObject;
    scoreObject.count = count;
    scoreObject.pattern = pattern;
    scoreObject.site = site;
    scoreObject.end_site = site + read_len - 1;
    return scoreObject;
}

void calc_score(std::vector<std::string> data_paths, int block_start, int block_end){
    try {

        int line_i = 0;
        std::deque<ScoreObject> lines_in_window;

        int num_files = data_paths.size();
        double total_score = 0;
        double total_count = 0;
        for (std::string data_path : data_paths){
            //"/cs/cbio/jon/projects/PyCharmProjects/wgbs_tools/src/segment_betas/trial_3.pat"
            std::ifstream input( data_path );

            for (std::string line_str; std::getline(input, line_str);) { //std::cin
                if (line_str.empty()) { continue; } // skip empty lines

                std::vector<std::string> tokens = line2tokens(line_str);

                if (tokens.size() < 4) {
                    throw std::invalid_argument("too few columns in file, line " + std::to_string(line_i));
                }
                else if (!(tokens.empty())) {

                    ScoreObject score_object = get_score_object_from_tokens(tokens);
                    if(score_object.site > block_end){
                        break;
                    }
                    if (score_object.end_site < block_start){
                        continue;
                    }
                    double min_proportion = (process_score_object(block_start, block_end, score_object)) * score_object.count;
                    total_count += score_object.count;
                    total_score += min_proportion;
                }
                else {
                    std::cerr << "something went wrong... tokens is empty" << std::endl;
                }
                line_i++;
            }
        }
        double avg_score = total_score / total_count;
        printf("%f ", avg_score);
        printf("\n");
    }
    catch(std::exception &e) {
        std::cerr << "failed segmenting pat" << std::endl;
        std::cerr << e.what() << std::endl;
        return;
    }
}

void PatSegmentor::scan_for_bimodal_2(std::vector<std::string> data_paths){

    try {
        UXM* uxMCounts = new UXM[nr_sites];
        for(int i = 0; i < nr_sites; ++i){
            UXM aUxm;
            aUxm.num_m = 0;
            aUxm.num_u = 0;
            aUxm.num_x = 0;
            uxMCounts[i] = aUxm;
        }
        int line_i = 0;
        std::vector<Block> blocks_list;
        for (std::string data_path : data_paths) {
            std::ifstream input(data_path);

            for (std::string line_str; std::getline(input, line_str);) { //std::cin
                if (line_str.empty()) { continue; } // skip empty lines

                std::vector<std::string> tokens = line2tokens(line_str);

                if (tokens.size() < 4) {
                    throw std::invalid_argument("too few columns in file, line " + std::to_string(line_i));
                } else if (!(tokens.empty())) {

                    ScoreObject score_object = get_score_object_from_tokens(tokens);
                    if (score_object.site > end_index) {
                        break;
                    }
                    if (score_object.end_site < start_index) {
                        continue;
                    }
                    if (first_site == -1) {
                        first_site = score_object.site;
                    }
                    int cur_site = score_object.site;
                    int cur_ind = cur_site - first_site;
                    if (score_object.pattern.length() > min_cpg_sites_per_read) {
                        int num_u = 0;
                        int num_m = 0;
                        for (int i = 0; i < score_object.pattern.length(); i++) {
                            char cur_letter = score_object.pattern[i];
                            if (cur_letter == 'C') {
                                num_m += 1;
                            } else if (cur_letter == 'T') {
                                num_u += 1;
                            }
                        }
                        if (num_u + num_m >= min_cpg_sites_per_read) {
                            float prop_u = (float) num_u / (float) (num_u + num_m);
                            float prop_m = 1 - prop_u;
                            if (prop_u >= read_homog_cutoff) { //0.8
                                for (int i = 0; i < score_object.pattern.length(); i++) {
                                    if (cur_ind + i >= 0){
                                        uxMCounts[cur_ind + i].num_u += (1 * score_object.count);
                                    }
                                }
                            } else if (prop_m >= read_homog_cutoff) {
                                for (int i = 0; i < score_object.pattern.length(); i++) {
                                    if (cur_ind + i >= 0){
                                        uxMCounts[cur_ind + i].num_m += (1 * score_object.count);
                                    }

                                }
                            } else {
                                for (int i = 0; i < score_object.pattern.length(); i++) {
                                    if (cur_ind + i >= 0){
                                        uxMCounts[cur_ind + i].num_x += (1 * score_object.count);
                                    }
                                }
                            }
                        }
                    }

                } else {
                    std::cerr << "something went wrong... tokens is empty" << std::endl;
                }
                line_i++;
            }
        }
        int end_site = 0;
        int start_site = 0;
        while (start_site < nr_sites - min_bimodal_window_len){ // 80
            end_site = start_site;
            int num_u_reads = 0;
            int num_x_reads = 0;
            int num_m_reads = 0;

            num_u_reads = uxMCounts[end_site].num_u;
            num_x_reads = uxMCounts[end_site].num_x;
            num_m_reads = uxMCounts[end_site].num_m;

            float total_reads = num_u_reads + num_x_reads + num_m_reads;
            float prop_u = total_reads >= min_reads_per_window ? (float) num_u_reads / total_reads : 0;
            float prop_m = total_reads >= min_reads_per_window ? (float) num_m_reads / total_reads : 0;
            Block cur_block;
            int num_fails = 1000;
            bool is_first = true;
            int length = 0;
            int continue_count = 0;
            while (((prop_m >= bimodal_win_threshold and prop_u >= bimodal_win_threshold) or num_fails < 3) and end_site < nr_sites - 1){
                if (is_first){
                    num_fails = 0;
                    is_first = false;
                }
                if (length >= 4){
                    cur_block.is_empty = false;
                    bimodal_win_threshold = 0.15;
                }
                cur_block.start = start_site;
                cur_block.end = end_site - num_fails - continue_count;
                end_site += 1;
                num_u_reads = uxMCounts[end_site].num_u;
                num_x_reads = uxMCounts[end_site].num_x;
                num_m_reads = uxMCounts[end_site].num_m;
                total_reads = num_u_reads + num_x_reads + num_m_reads;
                if (total_reads < min_reads_per_window){
                    continue_count += 1;
                    continue;
                }
                prop_u = (float) num_u_reads / total_reads;
                prop_m = (float) num_m_reads / total_reads;
                if (not (prop_m >= bimodal_win_threshold and prop_u >= bimodal_win_threshold)){
                    num_fails += 1;
                } else {
                    continue_count = 0;
                    length += 1;
                    num_fails = 0;
                }
                if (continue_count > 30){
                    break;
                }
            }
            if (!cur_block.is_empty){
                blocks_list.emplace_back(cur_block);
                start_site = end_site;
            } else {
                start_site += 1;
            }
        }
        delete uxMCounts;
        for (Block block : blocks_list){
            printf("%d ", block.start + first_site);
            printf("%d ", block.end + first_site - 1);
            printf("\n");
        }

    }
    catch(std::exception &e) {
        std::cerr << "failed segmenting pat" << std::endl;
        std::cerr << e.what() << std::endl;
        return;
    }
}


void PatSegmentor::scan_for_bimodal(std::string data_path){

    try {
        UXM* uxMCounts = new UXM[nr_sites];
        for(int i = 0; i < nr_sites; ++i){
            UXM aUxm;
            aUxm.num_m = 0;
            aUxm.num_u = 0;
            aUxm.num_x = 0;
            uxMCounts[i] = aUxm;
        }
        int line_i = 0;
        std::vector<Block> blocks_list;

        std::ifstream input( data_path );

        for (std::string line_str; std::getline(input, line_str);) { //std::cin
            if (line_str.empty()) { continue; } // skip empty lines

            std::vector<std::string> tokens = line2tokens(line_str);

            if (tokens.size() < 4) {
                throw std::invalid_argument("too few columns in file, line " + std::to_string(line_i));
            } else if (!(tokens.empty())) {

                ScoreObject score_object = get_score_object_from_tokens(tokens);
                if (score_object.site > end_index) {
                    break;
                }
                if (score_object.end_site < start_index) {
                    continue;
                }
                if (first_site == -1) {
                    first_site = score_object.site;
                }
                int cur_site = score_object.site;
                int cur_ind = cur_site - first_site;
                if (score_object.pattern.length() > min_cpg_sites_per_read) {
                    int num_u = 0;
                    int num_m = 0;
                    for (int i = 0; i < score_object.pattern.length(); i++) {
                        char cur_letter = score_object.pattern[i];
                        if (cur_letter == 'C') {
                            num_m += 1;
                        } else if (cur_letter == 'T') {
                            num_u += 1;
                        }
                    }
                    if (num_u + num_m > min_cpg_sites_per_read) {
                        float prop_u = (float) num_u / (float) (num_u + num_m);
                        float prop_m = 1 - prop_u;
                        if (prop_u >= read_homog_cutoff) { //0.8
                            uxMCounts[cur_ind].num_u += (1 * score_object.count);
                        } else if (prop_m >= read_homog_cutoff) {
                            uxMCounts[cur_ind].num_m += (1 * score_object.count);
                        } else {
                            uxMCounts[cur_ind].num_x += (1 * score_object.count);
                        }
                    }
                }

            } else {
                std::cerr << "something went wrong... tokens is empty" << std::endl;
            }
            line_i++;
        }
        int end_site = 0;
        int start_site = 0;
        while (start_site < nr_sites - min_bimodal_window_len){ // 80
            end_site = start_site;
            int num_u_reads = 0;
            int num_x_reads = 0;
            int num_m_reads = 0;
//            if (end_site + first_site >= 17135846){
//                int x = 0;
//            }
//            if (start_site + first_site >= 17135846){
//                int x = 0;
//            }
            while (end_site - start_site < min_bimodal_window_len){
                num_u_reads += uxMCounts[end_site].num_u;
                num_x_reads += uxMCounts[end_site].num_x;
                num_m_reads += uxMCounts[end_site].num_m;
                end_site += 1;
            }
            float total_reads = num_u_reads + num_x_reads + num_m_reads;
            float prop_u = total_reads > min_reads_per_window ? (float) num_u_reads / total_reads : 0;
            float prop_m = total_reads > min_reads_per_window ? (float) num_m_reads / total_reads : 0;
//            int num_fails = 100;
            Block cur_block;
            int cur_ind = start_site;
            while (prop_m >= bimodal_win_threshold and prop_u >= bimodal_win_threshold){ // 0.25
//                if (num_fails > 50){
//                    num_fails = 0;
//                }
                cur_block.is_empty = false;
                cur_block.start = start_site;
                cur_block.end = end_site;
                end_site += 1;
                num_u_reads += uxMCounts[end_site].num_u;
                num_x_reads += uxMCounts[end_site].num_x;
                num_m_reads += uxMCounts[end_site].num_m;
                num_u_reads -= uxMCounts[cur_ind].num_u;
                num_x_reads -= uxMCounts[cur_ind].num_x;
                num_m_reads -= uxMCounts[cur_ind].num_m;
                cur_ind += 1;
                total_reads = num_u_reads + num_x_reads + num_m_reads;
                prop_u = (float) num_u_reads / total_reads;
                prop_m = (float) num_m_reads / total_reads;
//                if (!(prop_m >= 0.25 and prop_u >= 0.25)){
//                    num_fails += 1;
//                } else {
//                    num_fails = 0;
//                }
            }
//            if (end_site + first_site >= 17135846){
//                int x = 0;
//            }
            if (!cur_block.is_empty){
                blocks_list.emplace_back(cur_block);
                start_site = end_site;
            } else {
                start_site += 1;
            }

        }
        delete uxMCounts;
        for (Block block : blocks_list){
//            if (block.end - block.start > 250){
//                printf("imp *** ");
//            }
            printf("%d ", block.start + first_site);
            printf("%d ", block.end + first_site);
            printf("\n");
        }

    }
    catch(std::exception &e) {
        std::cerr << "failed segmenting pat" << std::endl;
        std::cerr << e.what() << std::endl;
        return;
    }
}

void PatSegmentor::memoize_costs(std::vector<std::string> data_paths){

    try {

        int line_i = 0;
        std::deque<ScoreObject> lines_in_window;
        mem = new double*[nr_sites];
        for(int i = 0; i < nr_sites; ++i){
            mem[i] = new double[max_window_size];
            mem[i][0] = 1;
            for (int j = 1; j < max_window_size; j++){
                mem[i][j] = -1;
            }
        }
        int num_files = data_paths.size();
        for (std::string data_path : data_paths){
            //"/cs/cbio/jon/projects/PyCharmProjects/wgbs_tools/src/segment_betas/trial_3.pat"
            std::ifstream input( data_path );

            for (std::string line_str; std::getline(input, line_str);) { //std::cin
                if (line_str.empty()) { continue; } // skip empty lines

                std::vector<std::string> tokens = line2tokens(line_str);

                if (tokens.size() < 4) {
                    throw std::invalid_argument("too few columns in file, line " + std::to_string(line_i));
                }
                else if (!(tokens.empty())) {

                    ScoreObject score_object = get_score_object_from_tokens(tokens);
                    if(score_object.site > end_index){
                        break;
                    }
                    if (score_object.end_site < start_index){
                        continue;
                    }
                    if (first_site == -1){
                        first_site = score_object.site;
                    }
                    lines_in_window.push_back(score_object);
                    int available_window_size = score_object.site - lines_in_window[0].end_site;
                    if (available_window_size >= max_window_size){
                        process_cur_window(lines_in_window);
                    }
                    int cur_size = score_object.site - lines_in_window[0].end_site;
                    while (cur_size >= max_window_size){
                        lines_in_window.pop_front();
                        cur_size = score_object.site - lines_in_window[0].end_site;
                    }
                }
                else {
                    std::cerr << "something went wrong... tokens is empty" << std::endl;
                }
                line_i++;
            }
        }
        while (lines_in_window.size() > 0){
            process_one_step(lines_in_window);
            int cur_size = lines_in_window[lines_in_window.size() - 1].site - lines_in_window[0].end_site;
            int prev_size = cur_size;
            while (cur_size == prev_size){
                lines_in_window.pop_front();
                cur_size = lines_in_window[lines_in_window.size() - 1].site - lines_in_window[0].end_site;
            }
        }

        for(int i = 0; i < nr_sites; ++i){
            for (int j = 1; j < max_window_size; j++){
                mem[i][j] = mem[i][j] == -1 ? num_files : mem[i][j];
            }
        }
    }
    catch(std::exception &e) {
        std::cerr << "failed segmenting pat" << std::endl;
        std::cerr << e.what() << std::endl;
        return;
    }
}

std::vector<int> PatSegmentor::traceback(const int *T) {
    std::vector<int> borders;
    int i = nr_sites;
    borders.push_back(i + first_site);
    int border = std::max(0, T[i]);
    while (i > 0 and border == 0){
        i--;
        border = std::max(0, T[i]);
    }
    borders.push_back(border);
    i = std::max(0, T[i] - first_site);
    while (i > 0) {
        border = std::max(0, T[i]);
        i = std::max(0, T[i] - first_site);
        borders.push_back(border);
    }
    return borders;
}

void print_borders(std::vector<int> borders){
    for (auto i = borders.rbegin(); i != borders.rend(); ++i )
        printf("%d ", *i);
    printf("\n");
}

void PatSegmentor::dp(){
    /**
     * init M, T, all to zeros.
     *
     *  M[i] - The maximal score for segmenting sites 1,...,i
     *
     *  T - Traceback table.
     *  T[i] is the index of the last border we opened in the current segmentation
     */
    if (first_site > 0) {
        auto M = new double[nr_sites + 1]();
        auto T = new int[nr_sites + 1]();

        for (int j = 1; j < nr_sites; j++) {
            double best_score = std::numeric_limits<double>::infinity();
            int best_ind = -1;
            int start = std::max(0, j + 1 - max_window_size);
            for (int i = start; i < j; i++) {
                double tmp = M[i] + mem[i][j - i] + 0.0065;
                if (tmp < best_score) {
                    best_score = tmp;
                    best_ind = i + first_site;
                }
            }
            M[j] = best_score;
            T[j] = best_ind;
        }
        std::vector<int> borders = traceback(T);
        print_borders(borders);
        delete [] M; delete [] T;
    } else {
        for (int i = start_index; i < start_index + nr_sites; i += max_window_size - 1 )
            printf("%d ", i);
        printf("\n");
    }

    delete [] mem;
}

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
    pams.min_cpg_sites_per_read = 3;
    pams.read_homog_cutoff = 0.75;
    pams.min_bimodal_window_len = 50;
    pams.bimodal_win_threshold = 0.20;
    pams.min_reads_per_window = 5;

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

//    // max_bp
//    std::string max_bp_str = input.getCmdOption("-max_bp");
//    if (!max_bp_str.empty()) {
//        pams.max_bp = std::stoul(max_bp_str);
//    }

//    // rev_dict
//    std::string revd_str = input.getCmdOption("-rd");
//    if (!revd_str.empty()) {
//        pams.revdict = revd_str;
//    } else if (!max_bp_str.empty()) {
//        throw "Please provide path to revdict\n";
//    }

//    // pseudo count
//    std::string psud_str = input.getCmdOption("-ps");
//    if (!psud_str.empty()) {
//        pams.pseudo_count = std::stof(psud_str);
//    }

    std::string min_cpg_pre_read = input.getCmdOption("-cpgs_per_read");
    if (!min_cpg_pre_read.empty()) {
        pams.min_cpg_sites_per_read = std::stof(min_cpg_pre_read);
    }

    std::string read_homog_cutoff = input.getCmdOption("-read_homog_cutoff");
    if (!read_homog_cutoff.empty()) {
        pams.read_homog_cutoff = std::stof(read_homog_cutoff);
    }

    std::string bimodal_win_threshold = input.getCmdOption("-min_u_m_proportion");
    if (!bimodal_win_threshold.empty()) {
        pams.bimodal_win_threshold = std::stof(bimodal_win_threshold);
    }

    std::string min_bimodal_window_len = input.getCmdOption("-window_size");
    if (!min_bimodal_window_len.empty()) {
        pams.min_bimodal_window_len = std::stof(min_bimodal_window_len);
    }
    std::string min_reads_per_window = input.getCmdOption("-reads_per_window");
    if (!min_reads_per_window.empty()) {
        pams.min_reads_per_window = std::stof(min_reads_per_window);
    }

    return pams;
}



int main( int argc, char *argv[])
{
    InputParser input(argc, argv);
    if (argc < 6){
        std::cerr << "Usage: segment PAT_PATH [PAT_PATH...] -s START -n NR_SITES ";
        std::cerr << " [-m max_cpg] [-ps PSEUDO_COUNT]" << std::endl;
        return -1;
    }

    struct Params pams = parse_params(input);
//    struct Params pams = Params();

    // parse beta files paths
    std::vector<std::string> beta_paths;
    for (int i = 1; i < argc; i++) {
        std::string s(argv[i]);

        if ((s.length() >= 6) && !(s.compare(s.length() - 4, 4, ".pat"))) {
            beta_paths.emplace_back(s);
        }
    }
//    beta_paths.emplace_back("/cs/cbio/jon/check_prostate.tsv");

    try {
//        int start = 6072452;//std::stoi(argv[1]);
//        int end = 6231195;//std::stoi(argv[2]);

//        pams.start = 11160512; //6192452;
//        pams.nr_sites = 60000000;
//        pams.bimodal_win_threshold = 0.20;
//        pams.min_bimodal_window_len = 50;
//        pams.min_cpg_sites_per_read = 3;
//        pams.read_homog_cutoff = 0.75;
//        pams.min_reads_per_window = 5;

        PatSegmentor *p = new PatSegmentor(pams);
//        p->scan_for_bimodal(beta_paths.front());
        p->scan_for_bimodal_2(beta_paths);
//        p->memoize_costs(beta_paths);
//        p->dp();
        delete p;
    }
    catch(std::exception &e) {
        std::cerr << e.what() << std::endl;
        return -1;
    }
}
