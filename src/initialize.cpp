/* Data initializtion routines */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>
#include "global.h"
#include "rand.h"
#include <misc.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <fmt/core.h>
#include <nlohmann/json.hpp>
#include <map>
#include <filesystem>
#include <random>
#include <algorithm>

using json = nlohmann::json;
namespace fs = std::filesystem;

std::string prefix_init_file;
int sols_init_file;
int n_injected_points;
std::string injected_points_filename;


std::vector<std::string> my_s_h_u_vec;
struct var_t {
    std::string name;
    std::vector<int> location;
    double amount;
    int bmp;
};

void to_json(json &j, const var_t &p) {
    j = json{{"name",  p.name},
             {"location", p.location},
             {"amount",   p.amount},
             {"bmp", p.bmp}};
}

void from_json(const json &j, var_t &p) {
    j.at("name").get_to(p.name);
    j.at("location").get_to(p.location);
    j.at("amount").get_to(p.amount);
    j.at("bmp").get_to(p.bmp);
}
void read_json_lc_file(std::vector<std::vector<double> > &init_pop, std::string json_filename) {
    std::ifstream ifs_lc(json_filename);
    json jf_lc = json::parse(ifs_lc);
    std::unordered_map<std::string, double> loaded_lc_map;

    for (json j_lc: jf_lc) {
        for (json bmp_lc: j_lc) {
            var_t my_bmp = bmp_lc;
            auto s_lc = my_bmp.location[0];
            auto h_lc = my_bmp.location[1];
            auto u_lc = my_bmp.location[2];
            auto t_bmp_lc = my_bmp.location[3];
            std::string curr_parcel_lc = fmt::format("{}_{}_{}", s_lc, h_lc, u_lc, t_bmp_lc);
            loaded_lc_map[curr_parcel_lc] = my_bmp.amount;
        }
    }
}

void read_json_file(std::vector<std::vector<double> > &init_pop, std::string json_filename) {
    if(!fs::exists(json_filename)) return;
    fmt::print("jsonfilename: {}\n ", json_filename);

    std::ifstream ifs(json_filename);
    json jf = json::parse(ifs);
    int counter = 0;
    init_pop.resize(jf.size());

    for (json j: jf){
        std::vector<double> xreal(nreal, 0.0);
        for(json bmp: j) {
            var_t my_bmp = bmp;
            auto s = my_bmp.location[0];
            auto h = my_bmp.location[1];
            auto u = my_bmp.location[2];
            std::string curr_parcel = fmt::format("{}_{}_{}", s,h,u);

            if (std::find(my_s_h_u_vec.begin(), my_s_h_u_vec.end(), curr_parcel) == my_s_h_u_vec.end()) {
                my_s_h_u_vec.push_back(curr_parcel);
            }


            std::string s_tmp = fmt::format("{}_{}_{}_{}", s, h, u, my_bmp.bmp);
            fmt::print("stmp: {}", s_tmp);
            int idx = get_s_h_u_b(s_tmp);
            if (idx < 0){
                std::cout<<" READ_json Error, no s h u b in the file\n";
                exit(0);
            }
            else if(idx >nreal) {
                std::cout<<"Error idx out of bound";
                exit(0);
            }

            xreal[idx] = my_bmp.amount;

        }


        if (is_lc_enabled == true) {

            int lc_idx = lc_begin_;
            for (const auto &key: lc_keys_) {
                auto bmp_group = land_conversion_from_bmp_to[key];
                xreal[lc_idx] = 1.0;
                ++lc_idx;

                bool flag = false;
                if (std::find(my_s_h_u_vec.begin(), my_s_h_u_vec.end(), key) != my_s_h_u_vec.end()) {
                    //flag = true;
                }
                for (const auto &bmp: bmp_group) {
                    if (flag) {
                        xreal[lc_idx] = 0;
                    } else {
                        xreal[lc_idx] = rndreal(0.0, 0.05);
                    }
                    ++lc_idx;
                }
                /*
                for (const auto &bmp: bmp_group) {
                    xreal[lc_idx] = rndreal(0.0, 0.05);
                    ++lc_idx;
                }*/
            }
        }

        if (is_animal_enabled == true) {
            int animal_idx = animal_begin_;
            for (const auto &key: animal_keys_) {
                auto bmp_group = animal_complete_[key];
                xreal[animal_idx] = 1.0;
                ++animal_idx;
                for (const auto &bmp: bmp_group) {
                    xreal[animal_idx] = rndreal(0.0, 0.05);
                    ++animal_idx;
                }
            }
        }
        init_pop[counter] = xreal;
        counter++;
    }

    std::cout<<"Read solutions: "<<counter<<std::endl;

}
int get_random_number(int n) {
    std::random_device rd;                          // Obtain a random seed from the hardware
    std::mt19937 gen(rd());                         // Seed the random number generator
    std::uniform_int_distribution<int> dist(0, n);   // Define the range of the distribution

    return dist(gen);                               // Generate a random number within the defined range
}

void initialize_pop(population *pop)
{
    int i;
    std::vector<std::vector<double> > init_pop;
        read_json_file(init_pop, injected_points_filename);

        n_injected_points = (int) init_pop.size();
        for (i = 0; i < n_injected_points; i++) {
            if (i >= popsize)
                break;

            //std::copy(init_pop[i].begin(), init_pop[i].end(), pop->ind[i].xreal);
            for (int j = 0; j < nreal; j++) {
                pop->ind[i].xreal[j] = init_pop[i][j];
            }
        }
    for (i = n_injected_points; i < popsize; i++)
    {
        if (n_injected_points != 0) {
            int injected = get_random_number(n_injected_points-1);
            for (int j = 0; j < nreal; j++) {
                pop->ind[i].xreal[j] = init_pop[injected][j];
            }
        }
        initialize_ind(&(pop->ind[i]), n_injected_points);
    }
    return;
}

/* Function to initialize an individual randomly */
void initialize_ind(individual *ind, int n_injected_points)
{
    int j, k;
    if (nreal != 0)
    {
        for (j = 0; j < ef_end_; j++)
        {
            if(n_injected_points == 0) {
                ind->xreal[j] = rndreal(min_realvar[j], max_realvar[j]);
            }
        }
        if (is_lc_enabled == true) {
            int lc_idx = lc_begin_;
            for (const auto &key: lc_keys_) {
                auto bmp_group = land_conversion_from_bmp_to[key];
                ind->xreal[lc_idx] = 1.0;
                ++lc_idx;
                bool flag = false;
                if (std::find(my_s_h_u_vec.begin(), my_s_h_u_vec.end(), key) != my_s_h_u_vec.end()) {
                    flag = true;
                }
                for (const auto &bmp: bmp_group) {
                    if (flag) {
                        ind->xreal[lc_idx] = 0;
                    } else {
                        ind->xreal[lc_idx] = rndreal(0.0, 0.05);
                    }
                    ++lc_idx;
                }
            }
        }


        if (is_animal_enabled == true) {
            int animal_idx = animal_begin_;
            for (const auto &key: animal_keys_) {
                auto bmp_group = animal_complete_[key];
                ind->xreal[animal_idx] = 1.0;
                ++animal_idx;
                for (const auto &bmp: bmp_group) {
                    ind->xreal[animal_idx] = rndreal(0.0, 0.05);
                    ++animal_idx;
                }
            }
        }

    }

    if (nbin != 0)
    {
        for (j = 0; j < nbin; j++)
        {
            for (k = 0; k < nbits[j]; k++)
            {
                if (randomperc() <= 0.5)
                {
                    ind->gene[j][k] = 0;
                }
                else
                {
                    ind->gene[j][k] = 1;
                }
            }
        }
    }
    return;
}

