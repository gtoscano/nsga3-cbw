/* NSGA-III routine (implementation of the 'main' function) */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "global.h"
#include "rand.h"
#include <iostream>
#include <string>
#include <misc.hpp>
#include <vector>
#include <fmt/core.h>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
/*
 * For: Test Evaluation
 * Program Arguments: 0.7946199 be2e77b4-da21-4d4d-ba3f-bf88bae73618 corecast 0 ./ 20 1 3000 /tmp/be2e77b4-da21-4d4d-ba3f-bf88bae73618.csv 0 0
 * Redirect: /opt/opt4cast/alg/nsga34cast/input_data/be2e77b4-da21-4d4d-ba3f-bf88bae73618_cbw.in
 */

/*
 * For: Test NSGA-III (no injection)
 * Program Arguments: 0.7946199 be2e77b4-da21-4d4d-ba3f-bf88bae73618 corecast+smartinit 10 /opt/opt4cast/output/nsga3/be2e77b4-da21-4d4d-ba3f-bf88bae73618/config/ipopt.json 20 10 0 ./ 0.3 10
 * Redirect: /opt/opt4cast/alg/nsga34cast/input_data/be2e77b4-da21-4d4d-ba3f-bf88bae73618_cbw.in
 */
//run_exec 3511ee35-f9cf-4540-8356-f044b3e6f9d9 acbae501-814c-4133-9608-7f67c56ad970
// empty_38_6611_256_6_8_59_1_6608_36_2_31_8_406
// 20 2 2 0.2 20 0 0 8
//
// /home/gtoscano/django/api4opt4/optimization/nsga3/build/nsga3 0.98108506 3511ee35-f9cf-4540-8356-f044b3e6f9d9 corecast+smartinit 20 ./ 20 2 0 /tmp/3511ee35-f9cf-4540-8356-f044b3e6f9d9.csv 0.2 20 2 8 < /opt/opt4cast/output/nsga3/3511ee35-f9cf-4540-8356-f044b3e6f9d9/config/3511ee35-f9cf-4540-8356-f044b3e6f9d9_cbw.in

namespace fs = std::filesystem;

using json = nlohmann::json;
void split_str(std::string const &str, const char delim,
               std::vector<std::string> &out);

extern std::vector<double> alpha_stored_vec;
int nreal;
int nbin;
int nobj;
int ncon;
int popsize;
double pcross_real;
double pcross_bin;
double pmut_real;
double pmut_bin;
double eta_c;
double eta_m;
int ngen;
int nbinmut;
int nrealmut;
int nbincross;
int nrealcross;
int *nbits;
double *min_realvar;
double *max_realvar;


double *min_binvar;
double *max_binvar;
int bitlength;
int choice;
int obj1;
int obj2;
int obj3;
int angle1;
int angle2;

/* Extra Vaiables for MNSGA-II*/

double *ref_pt[NMAX];
double *ideal_point;
double **plane_point;
double scaling;

int steps;
int nref;
int curr_gen;
int onref;
int created_around[NMAX];
int adaptive_increment;
int active_ref_pts[10];
int start_incr;

/* Extra varaibles for CBW */

/* end CBW */
namespace my_execute{

    struct CommandResult {
        std::string output;
        int exitstatus;
        friend std::ostream &operator<<(std::ostream &os, const CommandResult &result) {
            os << "command exitstatus: " << result.exitstatus << " output: " << result.output;
            return os;
        }
        bool operator==(const CommandResult &rhs) const {
            return output == rhs.output &&
                   exitstatus == rhs.exitstatus;
        }
        bool operator!=(const CommandResult &rhs) const {
            return !(rhs == *this);
        }
    };

    class Command {
    public:
        /**
             * Execute system command and get STDOUT result.
             * Regular system() only gives back exit status, this gives back output as well.
             * @param command system command to execute
             * @return commandResult containing STDOUT (not stderr) output & exitstatus
             * of command. Empty if command failed (or has no output). If you want stderr,
             * use shell redirection (2&>1).
             */
        static CommandResult exec(const std::string &command) {
            int exitcode = 0;
            std::array<char, 1048576> buffer {};
            std::string result;
            FILE *pipe = popen(command.c_str(), "r");
            if (pipe == nullptr) {
                throw std::runtime_error("popen() failed!");
            }
            try {
                std::size_t bytesread;
                while ((bytesread = std::fread(buffer.data(), sizeof(buffer.at(0)), sizeof(buffer), pipe)) != 0) {
                    result += std::string(buffer.data(), bytesread);
                }
            } catch (...) {
                pclose(pipe);
                throw;
            }
            exitcode = WEXITSTATUS(pclose(pipe));
            return CommandResult{result, exitcode};
        }
    };
}

void scenario__compute_lc_keys() {

    for (const auto& pair : land_conversion_from_bmp_to) {
        lc_keys_.push_back(pair.first);
    }

    // Sort the vector of keys
    std::sort(lc_keys_.begin(), lc_keys_.end());
}

void scenario__compute_animal_keys() {

    for (const auto& pair : animal_complete_) {
        animal_keys_.push_back(pair.first);
    }

    // Sort the vector of keys
    std::sort(animal_keys_.begin(), animal_keys_.end());
}

void scenario__load(const std::string& filename) {

    // Open the JSON file
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        exit(-1);
        return;
    }

    // Parse the JSON file directly into a nlohmann::json object
    json json_obj = json::parse(file);

    // Access the JSON data
    amount_ = json_obj["amount"].get<std::unordered_map<std::string, double>>();
    land_conversion_from_bmp_to = json_obj["land_conversion_to"].get<std::unordered_map<std::string, std::vector<std::string>>>();
    bmp_cost_ = json_obj["bmp_cost"].get<std::unordered_map<std::string, double>>();
    animal_complete_ = json_obj["animal_complete"].get<std::unordered_map<std::string, std::vector<int>>>();
    animal_ = json_obj["animal_unit"].get<std::unordered_map<std::string, double>>();
    scenario__compute_lc_keys();
    scenario__compute_animal_keys();
}


int scenario__compute_lc_size() {
    int counter = 0;
    for (const auto &[key, bmp_group]: land_conversion_from_bmp_to) {
        ++counter;
        for (const auto &bmp: bmp_group) {
            ++counter;
        }
    }
    return counter;
}

int scenario__compute_animal_size() {
    int counter = 0;
    for (const auto &[key, bmp_group]: animal_complete_) {
        ++counter;
        for (const auto &bmp: bmp_group) {
            ++counter;
        }
    }
    return counter;
}


double scenario__alpha_minus2(const std::string& key, double original_amount) {
    double minus = 0.0;
    double original = original_amount;
    /*
    if (amount_.contains(key)) {
        original = amount_.at(key);
        std::cout<<"Original "<<original<<std::endl;
    }
     */

    if (amount_minus_.contains(key)) {
        minus = amount_minus_.at(key);
        //std::cout<<"My Minus "<<minus<<std::endl;
    }

    return original - minus;
}

double scenario__alpha_plus_minus2(const std::string& key, double original_amount) {

    double plus = 0.0;

    if (amount_plus_.contains(key)) {
        plus = amount_plus_.at(key);
    }

    return plus + scenario__alpha_minus2(key, original_amount);
}

void scenario__compute_lc2() {
    for (const auto& entry : lc_x_) {
        auto [s, h, u, bmp, amount] = entry;
        std::cout<<fmt::format("S: {}, h: {}, u: {}, bmp: {}, amount: {}\n", s, h, u, bmp, amount);
    }
}

void scenario__sum_alpha2(std::unordered_map<std::string, double>& am,  const std::string& key, double acreage) {
    double tmp = 0.0;
    if (am.contains(key)) {
        tmp = am.at(key);
    }
    tmp += acreage;
    am[key] = tmp;
}

void scenario__load_alpha2(std::vector<double>& my_alpha) {
    size_t parcel_idx = 0;
    my_alpha.clear();
    for (auto &&bmp_per_parcel: bmp_grp_src_links_vec) {
        int s = s_h_u_vec[parcel_idx][0];
        int h = s_h_u_vec[parcel_idx][1];
        int u = s_h_u_vec[parcel_idx][2];
        auto from_key = fmt::format("{}_{}_{}", s, h, u);
        //auto new_alpha = scenario__alpha_minus(from_key, alpha_stored_vec[parcel_idx]);
        double new_alpha = scenario__alpha_plus_minus2(from_key, alpha_stored_vec[parcel_idx]);
        if (alpha_vec[parcel_idx] != new_alpha) {
            //std::cout << "New Alpha: " << new_alpha << " Previos Alpha: " << alpha_vec[parcel_idx] << " My other alpha " << amount_[from_key] << std::endl;
        }

        my_alpha.push_back(new_alpha);
        ++parcel_idx;
    }
}

double scenario__normalize_animal2(const std::vector<double>& x) {
    int counter = animal_begin_;
    double total_cost = 0.0;
    animal_x_.clear();
    //std::vector<std::string> lc_altered_keys;
    for (const std::string& key : animal_keys_) {
        std::vector<int> bmp_group =  animal_complete_[key];

        std::vector<std::pair<double, int>> grp_tmp;
        std::vector <std::string> key_split;
        split_str(key, '_', key_split);
        std::string base_condition = key_split[0];
        std::string county = key_split[1];
        std::string load_source = key_split[2];
        std::string animal_id = key_split[3];

        double sum = x[counter];
        ++counter;

        for (int bmp : bmp_group) {
            grp_tmp.push_back({x[counter], bmp});
            sum += x[counter];
            ++counter;
        }

        double pct_accum = 0.0;
        fmt::print("grp_tmp: {}\n", grp_tmp.size());
        for (auto [pct, bmp]: grp_tmp) {
            double norm_pct =  (double) pct / sum;
            if (norm_pct * animal_[key] > 1.0) {
                double amount = (norm_pct * animal_[key]);
                std::string key_bmp_cost = fmt::format("{}_{}", 8, bmp);
                double cost = amount * bmp_cost_[key_bmp_cost];
                total_cost += cost;
                animal_x_.push_back({std::stoi(base_condition), std::stoi(county), std::stoi(load_source), std::stoi(animal_id), bmp, amount});
            }
        }
    }

    return total_cost;
}

double scenario__normalize_land_conversion2(const std::vector<double>& x) {
    int counter = lc_begin_;
    amount_minus_.clear();
    amount_plus_.clear();
    lc_x_.clear();
    double total_cost = 0.0;
    //std::vector<std::string> lc_altered_keys;
    for (const auto& key : lc_keys_) {
        std::vector<std::string> bmp_group =  land_conversion_from_bmp_to[key];
        std::vector<std::pair<double, std::string>> grp_tmp;
        std::vector <std::string> key_split;
        double sum = x[counter];
        split_str(key, '_', key_split);
        ++counter;

        for (std::string bmp : bmp_group) {
            grp_tmp.push_back({x[counter], bmp});
            sum += x[counter];
            ++counter;
        }

        double pct_accum = 0.0;
        std::vector<std::pair<double, std::string>> grp_pct_tmp;
        for (auto line : grp_tmp) {
            double pct = line.first;
            std::string to = line.second;
            double norm_pct =  pct / sum;
            std::vector <std::string> out_to;
            split_str(to, '_', out_to);
            auto bmp = std::stoi(out_to[0]);
            auto key_to = fmt::format("{}_{}_{}", key_split[0], key_split[1], out_to[1]);
            scenario__sum_alpha2(amount_plus_, key_to, norm_pct * amount_[key]);
            pct_accum += norm_pct;
            if (norm_pct * amount_[key] > 1.0) {
                double amount = (norm_pct * amount_[key]);
                auto key_bmp_cost = fmt::format("{}_{}", 8, bmp);
                double cost = amount * bmp_cost_[key_bmp_cost];
                total_cost += cost;
                lc_x_.push_back({std::stoi(key_split[0]), std::stoi(key_split[1]), std::stoi(key_split[2]), bmp, norm_pct * amount_[key]});
            }
        }
        scenario__sum_alpha2(amount_minus_, key, pct_accum * amount_[key]);
    }
    return total_cost;
}


int main(int argc, char **argv)
{
    int i, j;
    FILE *fpt1;
    FILE *fpt2;
    FILE *fpt3;
    FILE *fpt4;
    FILE *fpt5;
    FILE *fpt6;
    FILE *gp=NULL;
    population *parent_pop;
    population *child_pop;
    population *mixed_pop;
    std::string algorithm = "nsga2";
    n_injected_points = 0;
    if (argc < 2)
    {
      std::clog<<"Usage ./nsga3r random_seed UUID  [corecast | mathmodel | corecast+smartinit | mathmodel+smartinit] #injected_points injected_points_filename\n";
      std::clog<<"Usage Example: ./nsga3r 0.1 00000000-0000-0000-0000-000000000075 mathmodel 11 injected_points/tucker_vw_nds.in < ../input_data/tucker_vw.in\n";
        exit(1);
    }
    seed = (double) std::stof(argv[1]);
    emo_uuid = argv[2];
    std::string eval_options[] = {"corecast", "mathmodel", "corecast+smartinit", "mathmodel+smartinit"};
    opt4cast_evaluation = argv[3];
    n_injected_points = std::stoi(argv[4]);
    injected_points_filename = argv[5];
    popsize = std::stoi(argv[6]);//4
    //popsize=100;
    ngen = std::stoi(argv[7]); //1
    //ngen=20;
    mode = std::stoi(argv[8]);//0 is for NSGA-II >= 1 is for evaluate
    std::string scenario_filename = argv[9]; // /tmp/...csv
    double ipopt_reduction = std::stof(argv[10]); //0
    double ipopt_popsize = std::stoi(argv[11]); //0
    int corecast_gen = std::stoi(argv[12]); //2
    int cost_profile_idx = std::stoi(argv[13]); //2

    int pollutant_idx = 0;
    double limit_alpha = 1.0;

    scenario__load("/opt/opt4cast/csvs/prueba.json");
    std::string msu_cbpo_path = getEnvVar("MSU_CBPO_PATH", "/opt/opt4cast");
    if (find(cbegin(eval_options), cend(eval_options), opt4cast_evaluation) == cend(eval_options)) {
       std::cerr << "Invalid option: " << opt4cast_evaluation << ", please use [send | retrieve]\n";
       return -1;
    }
    if (opt4cast_evaluation =="corecast+smartinit" || opt4cast_evaluation == "mathmodel+smartinit") {
        if (injected_points_filename == "./") {
            injected_points_filename = fmt::format("{}/output/nsga3/{}/config/ipopt.json", msu_cbpo_path, emo_uuid);
        }
        if (!fs::exists(injected_points_filename)) {

            //be2e77b4-da21-4d4d-ba3f-bf88bae73618 0.3 0 1 1.0 4 10
            // 1        2                                 3                  4           5          6        7    8     9         10         11
            // RND      EMO-UUID opt4cast_evaluation n_injected  n_file    popsize   ngen  mode scen_file                                                                   reduction  cc_gen
            //0.7946199 be2e77b4 corecast+smartinit 10          ipopt.json 20        1     0    /config/  0.5   1

            int ipopt = 1;/*
            std::string exec_string = fmt::format("/home/gtoscano/CLionProjects/msu_cbpo_opt/epa {} {} {} {} {} {} {}", emo_uuid, ipopt_reduction, pollutant_idx, ipopt, limit_alpha, cost_idx, ipopt_popsize);
            */

            std::string env_var = "OPT4CAST_EPS_CNSTR_PATH";
            std::string EPS_CNSTR_PATH = getEnvVar(env_var);
            std::string exec_string = fmt::format("{} {} {} {} {} {} {} {}",
                                                  EPS_CNSTR_PATH, emo_uuid, ipopt_reduction, pollutant_idx, ipopt, limit_alpha, cost_profile_idx, ipopt_popsize);
            using namespace my_execute;
            CommandResult nullbyteCommand = Command::exec(exec_string); // NOLINT(bugprone-string-literal-with-embedded-nul)
            std::ofstream ofile("/tmp/filename2.txt");
            ofile<<exec_string<<std::endl;
            ofile << "Output using fread: " << nullbyteCommand << std::endl;
            ofile.close();


        }

        std::string dir_path = fmt::format("{}/output/nsga3/{}/front", msu_cbpo_path, emo_uuid);

        if(!fs::exists(dir_path)) //Does not exist
            fs::create_directories(dir_path);
        std::string filename_src = fmt::format("{}/output/nsga3/{}/config/ipopt.json", msu_cbpo_path, emo_uuid);
        std::string filename_dst = fmt::format("{}/output/nsga3/{}/front/epsilon.json", msu_cbpo_path, emo_uuid);

        if(fs::exists(filename_src)) //Does it exist
            fs::copy(filename_src, filename_dst, fs::copy_options::update_existing);
    }
    else {
        n_injected_points = 0;
        injected_points_filename = "";
    }

    int my_new_variables = 0;
    int animal_variables = 0;
    if(is_lc_enabled) {
        my_new_variables = scenario__compute_lc_size();
    }
    if(is_animal_enabled) {
        animal_variables = scenario__compute_animal_size();
    }

    std::cout<<"My new variables: "<<my_new_variables<<std::endl;
    if (seed <= 0.0 || seed >= 1.0)
    {
        printf("\n Entered seed value is wrong, seed value must be in (0,1) \n");
        exit(1);
    }
  try {
    std::string dir_path = fmt::format("{}/output/nsga3/{}", msu_cbpo_path, emo_uuid);
    std::filesystem::create_directories(dir_path);
  }
  catch(const std::exception & error) {
      std::cout<<"error\n"<<error.what()<<"\n";
      std::clog << fmt::format("when saving files last generation")<<error.what();
  }
  try {
    std::string dir_path = fmt::format("{}/output/nsga3/{}/front", msu_cbpo_path, emo_uuid);
    std::filesystem::create_directories(dir_path);
  }
  catch(const std::exception & error) {
      std::cout<<"error\n"<<error.what()<<"\n";
      std::clog << fmt::format("when saving files last generation")<<error.what();
  }
    std::string stored_path = fmt::format("{}/output/nsga3/{}/front/", msu_cbpo_path, emo_uuid);
    std::string specific_path = stored_path + "initial_pop.out";
    fpt1 = fopen(specific_path.c_str(), "w");
    specific_path = stored_path + "final_pop.out";
    fpt2 = fopen(specific_path.c_str(), "w");
    specific_path = stored_path + "best_pop.out";
    fpt3 = fopen(specific_path.c_str(), "w");
    specific_path = stored_path + "all_pop.out";
    fpt4 = fopen(specific_path.c_str(), "w");
    specific_path = stored_path + "params.out";
    fpt5 = fopen(specific_path.c_str(), "w");
    specific_path = stored_path + "ref.out";
    fpt6 = fopen(specific_path.c_str(), "w");
    fprintf(fpt1, "# This file contains the data of initial population\n");
    fprintf(fpt2, "# This file contains the data of final population\n");
    fprintf(fpt3, "# This file contains the data of final feasible population (if found)\n");
    fprintf(fpt4, "# This file contains the data of all generations\n");
    fprintf(fpt5, "# This file contains information about inputs as read by the program\n");

    //printf("\n Enter the number of objectives : ");
    scanf("%d", &nobj);
    printf("nobj: %d\n", nobj);

    if (nobj < 1)
    {
        printf("\n number of objectives entered is : %d", nobj);
        printf("\n Wrong number of objectives entered, hence exiting \n");
        exit(1);
    }
    //printf("\n Enter the number of steps in structured reference points (0 if preferential reference points are supplied): ");
    scanf("%d", &steps);
    printf("steps: %d\n", steps);
    steps = (int) popsize * 0.90;

    if (steps < 0)
    {
        printf("\n number of steps entered is : %d", steps);
        printf("\n Wrong number of steps entered, hence exiting \n");
        exit(1);
    }
    scaling = (double)steps;
    if (steps == 0)
    {
        printf("\n Enter the number of preferential reference points : ");
        scanf("%d", &nref);
        onref = nref;
        scaling = 20.0;
    }

    //printf("\n Do you want reference points to be adaptive (0 for NO) (1 for yes) : ");
    scanf("%d", &adaptive_increment);
    printf("adaptive_increment: %d\n", adaptive_increment);
    adaptive_increment = 1;

    //printf("\n Enter the population size (a multiple of 4) : ");
    int my_popsize;
    scanf("%d", &my_popsize);
    printf("my_popsize: %d\n", my_popsize);
    if (popsize < 4 || (popsize % 4) != 0)
    {
      std::clog<<fmt::format("\n population size read is : {}\n", popsize);
      std::clog<<"\n Wrong population size entered, hence exiting \n";
        exit(1);
    }
//printf("\n Enter the number of generations : ");
    int my_ngen;
    scanf("%d", &my_ngen);
    printf("my_ngen: %d\n", my_ngen);

    //printf("\n Enter the number of constraints : ");
    scanf("%d", &ncon);
    printf("ncon: %d\n", ncon);
    ncon++;
    if (ncon < 0)
    {
      std::clog<<fmt::format("\n number of constraints entered is : {}\n", ncon);
      std::clog<<"Wrong number of constraints enetered, hence exiting \n";
        exit(1);
    }
    //printf("\n Enter the number of real variables : ");
    scanf("%d", &nreal);
    printf("nreal: %d\n", nreal);

        ef_size_ = nreal;
        ef_begin_ = 0;
        ef_end_ = ef_begin_ + ef_size_ - 1;

    if (is_lc_enabled) {
        lc_size_ = my_new_variables;
        lc_begin_ = nreal;
        lc_end_ = lc_begin_ + lc_size_ - 1;
    }

    if(is_animal_enabled) {
        animal_size_ = animal_variables;
        animal_begin_ = ef_size_ + lc_size_;
        animal_end_ = animal_begin_ + animal_size_ - 1;
    }

    nreal += (my_new_variables + animal_variables);

    if (nreal < 0)
    {
      std::clog<<fmt::format("\n number of real variables entered is : {}\n", nreal);
      std::clog<<"Wrong number of variables entered, hence exiting \n";
        exit(1);
    }

    if (nreal != 0)
    {
        min_realvar = (double *)malloc(nreal * sizeof(double));
        max_realvar = (double *)malloc(nreal * sizeof(double));
        for (i = 0; i < nreal; i++)
        {
            //printf("\n Enter the lower limit of real variable %d : ", i + 1);
            //scanf("%lf", &min_realvar[i]);
            min_realvar[i] = 0.0;
            //printf("\n Enter the upper limit of real variable %d : ", i + 1);
            //scanf("%lf", &max_realvar[i]);
            //max_realvar[i]);
            max_realvar[i] = 1.0;
            if (max_realvar[i] <= min_realvar[i])
            {
              std::clog<<"\n Wrong limits entered for the min and max bounds of real variable, hence exiting \n";
                exit(1);
            }
        }
        //printf("\n Enter the probability of crossover of real variable (0.6-1.0) : ");
        scanf("%lf", &pcross_real);
        printf("pcross_real: %lf\n", pcross_real);
        if (pcross_real < 0.0 || pcross_real > 1.0)
        {
          std::clog<<fmt::format("\n Probability of crossover entered is : {}\n", pcross_real);
          std::clog<<"Entered value of probability of crossover of real variables is out of bounds, hence exiting \n";
            exit(1);
        }
        //printf("\n Enter the probablity of mutation of real variables (1/nreal) : ");
        /* scanf ("%lf",&pmut_real);*/
        pmut_real = 1.0 / (double)nreal;
        if (pmut_real < 0.0 || pmut_real > 1.0)
        {
          std::clog<<fmt::format("\n Probability of mutation entered is : {}\n", pmut_real);
          std::clog<<"Entered value of probability of mutation of real variables is out of bounds, hence exiting \n";
            exit(1);
        }
        //printf("\n Enter the value of distribution index for crossover (5-20): ");
        scanf("%lf", &eta_c);
        printf("eta_c: %lf\n", eta_c);
        if (eta_c <= 0)
        {
          std::clog<<fmt::format("\n The value entered is : {}\n", eta_c);
          std::clog<<"Wrong value of distribution index for crossover entered, hence exiting \n";
            exit(1);
        }
        //printf("\n Enter the value of distribution index for mutation (5-50): ");
        scanf("%lf", &eta_m);
        printf("eta_m: %lf\n", eta_m);
        if (eta_m <= 0)
        {
          std::clog<<fmt::format("\n The value entered is : %e\n", eta_m);
          std::clog<<"Wrong value of distribution index for mutation entered, hence exiting \n";
            exit(1);
        }
    }
    //printf("\n Enter the number of binary variables : ");
    scanf("%d", &nbin);
    printf("nbin: %d\n", nbin);
    if (nbin != 0)
    {
      std::clog<<fmt::format("Number of binary variables entered is : {}\n", nbin);
      std::clog<<"Wrong number of binary variables entered, hence exiting \n";
        exit(1);
    }
    if (nbin != 0)
    {
        nbits = (int *)malloc(nbin * sizeof(int));
        min_binvar = (double *)malloc(nbin * sizeof(double));
        max_binvar = (double *)malloc(nbin * sizeof(double));
        for (i = 0; i < nbin; i++)
        {
            //printf("\n Enter the number of bits for binary variable %d : ", i + 1);
            scanf("%d", &nbits[i]);
            if (nbits[i] < 1)
            {
              std::clog<<"\n Wrong number of bits for binary variable entered, hence exiting\n";
                exit(1);
            }
            //printf("\n Enter the lower limit of binary variable %d : ", i + 1);
            scanf("%lf", &min_binvar[i]);
            //printf("\n Enter the upper limit of binary variable %d : ", i + 1);
            scanf("%lf", &max_binvar[i]);
            if (max_binvar[i] <= min_binvar[i])
            {
              std::clog<<"Wrong limits entered for the min and max bounds of binary variable entered, hence exiting\n";
                exit(1);
            }
        }
        //printf("\n Enter the probability of crossover of binary variable (0.6-1.0): ");
        scanf("%lf", &pcross_bin);
        if (pcross_bin < 0.0 || pcross_bin > 1.0)
        {
          std::clog<<fmt::format("Probability of crossover entered is : {}\n", pcross_bin);
            std::clog<<"Entered value of probability of crossover of binary variables is out of bounds, hence exiting \n";
            exit(1);
        }
        //printf("\n Enter the probability of mutation of binary variables (1/nbits): ");
        scanf("%lf", &pmut_bin);
        if (pmut_bin < 0.0 || pmut_bin > 1.0)
        {
          std::clog<<fmt::format("Probability of mutation entered is : %e\n", pmut_bin);
          std::clog<<"Entered value of probability  of mutation of binary variables is out of bounds, hence exiting \n";
            exit(1);
        }
    }
    if (nreal == 0 && nbin == 0)
    {
        std::clog<<"\n Number of real as well as binary variables, both are zero, hence exiting\n";
        exit(1);
    }
    choice = 0;
    //std::clog<<"\n Do you want to use gnuplot to display the results realtime (0 for NO) (1 for yes) :\n";
    scanf("%d", &choice);
    if (choice != 0 && choice != 1)
    {
      std::clog<<fmt::format("Entered the wrong choice, hence exiting, choice entered was {}\n", choice);
        exit(1);
    }
    if (choice == 1)
    {
        gp = popen(GNUPLOT_COMMAND, "w");
        if (gp == NULL)
        {
          std::clog<<"Could not open a pipe to gnuplot, check the definition of GNUPLOT_COMMAND in file global.h\n";
          std::clog<<"Edit the string to suit your system configuration and rerun the program\n";
            exit(1);
        }
        if (nobj == 2)
        {
            //printf("\n Enter the objective for X axis display : ");
            scanf("%d", &obj1);
            if (obj1 < 1 || obj1 > nobj)
            {
              std::clog<<fmt::format("Wrong value of X objective entered, value entered was {}\n", obj1);
                exit(1);
            }
            //printf("\n Enter the objective for Y axis display : ");
            scanf("%d", &obj2);
            if (obj2 < 1 || obj2 > nobj)
            {
              std::clog<<fmt::format("\n Wrong value of Y objective entered, value entered was {}\n", obj2);
                exit(1);
            }
            obj3 = -1;
        }
        else
        {
            //printf("\n #obj > 2, 2D display or a 3D display ?, enter 2 for 2D and 3 for 3D :");
            scanf("%d", &choice);
            if (choice != 2 && choice != 3)
            {
              std::clog<<fmt::format("Entered the wrong choice, hence exiting, choice entered was {}\n", choice);
                exit(1);
            }
            if (choice == 2)
            {
                //printf("\n Enter the objective for X axis display : ");
                scanf("%d", &obj1);
                if (obj1 < 1 || obj1 > nobj)
                {
                    std::clog<<fmt::format("\n Wrong value of X objective entered, value entered was {}\n", obj1);
                    exit(1);
                }
                //printf("\n Enter the objective for Y axis display : ");
                scanf("%d", &obj2);
                if (obj2 < 1 || obj2 > nobj)
                {
                    std::clog<<fmt::format("\n Wrong value of Y objective entered, value entered was {}\n", obj2);
                    exit(1);
                }
                obj3 = -1;
            }
            else
            {
                //printf("\n Enter the objective for X axis display : ");
                scanf("%d", &obj1);
                if (obj1 < 1 || obj1 > nobj)
                {
                    std::clog<<fmt::format("\n Wrong value of X objective entered, value entered was {}\n", obj1);
                    exit(1);
                }
                //printf("\n Enter the objective for Y axis display : ");
                scanf("%d", &obj2);
                if (obj2 < 1 || obj2 > nobj)
                {
                    std::clog<<fmt::format("Wrong value of Y objective entered, value entered was {}\n", obj2);
                    exit(1);
                }
                //printf("\n Enter the objective for Z axis display : ");
                scanf("%d", &obj3);
                if (obj3 < 1 || obj3 > nobj)
                {
                    std::clog<<fmt::format("Wrong value of Z objective entered, value entered was {}\n", obj3);
                    exit(1);
                }
                //printf("\n You have chosen 3D display, hence location of eye required \n");
                //printf("\n Enter the first angle (an integer in the range 0-180) (if not known, enter 60) :");
                scanf("%d", &angle1);
                if (angle1 < 0 || angle1 > 180)
                {
                    std::clog<<"Wrong value for first angle entered, hence exiting\n";
                    exit(1);
                }
                //printf("\n Enter the second angle (an integer in the range 0-360) (if not known, enter 30) :");
                scanf("%d", &angle2);
                if (angle2 < 0 || angle2 > 360)
                {
                    std::clog<<"Wrong value for second angle entered, hence exiting\n";
                    exit(1);
                }
            }
        }
    }

    if (mode > 0) {
      popsize = 1;
      ngen = 1;
    }
    //printf("\n Input data successfully entered, now performing initialization \n");
    fprintf(fpt5, "\n Population size = %d", popsize);
    fprintf(fpt5, "\n Number of generations = %d", ngen);
    fprintf(fpt5, "\n Number of objective functions = %d", nobj);
    fprintf(fpt5, "\n Number of constraints = %d", ncon);
    fprintf(fpt5, "\n Number of real variables = %d", nreal);
    if (nreal != 0)
    {
        for (i = 0; i < nreal; i++)
        {
            fprintf(fpt5, "\n Lower limit of real variable %d = %e", i + 1, min_realvar[i]);
            fprintf(fpt5, "\n Upper limit of real variable %d = %e", i + 1, max_realvar[i]);
        }
        fprintf(fpt5, "\n Probability of crossover of real variable = %e", pcross_real);
        fprintf(fpt5, "\n Probability of mutation of real variable = %e", pmut_real);
        fprintf(fpt5, "\n Distribution index for crossover = %e", eta_c);
        fprintf(fpt5, "\n Distribution index for mutation = %e", eta_m);
    }
    fprintf(fpt5, "\n Number of binary variables = %d", nbin);
    if (nbin != 0)
    {
        for (i = 0; i < nbin; i++)
        {
            fprintf(fpt5, "\n Number of bits for binary variable %d = %d", i + 1, nbits[i]);
            fprintf(fpt5, "\n Lower limit of binary variable %d = %e", i + 1, min_binvar[i]);
            fprintf(fpt5, "\n Upper limit of binary variable %d = %e", i + 1, max_binvar[i]);
        }
        fprintf(fpt5, "\n Probability of crossover of binary variable = %e", pcross_bin);
        fprintf(fpt5, "\n Probability of mutation of binary variable = %e", pmut_bin);
    }
    fprintf(fpt5, "\n Seed for random number generator = %e", seed);
    fprintf(fpt5, "\n Approach %d used for reference points", adaptive_increment);
    bitlength = 0;
    if (nbin != 0)
    {
        for (i = 0; i < nbin; i++)
        {
            bitlength += nbits[i];
        }
    }
    fprintf(fpt1, "# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n", nobj, ncon, nreal, bitlength);
    fprintf(fpt2, "# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n", nobj, ncon, nreal, bitlength);
    fprintf(fpt3, "# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n", nobj, ncon, nreal, bitlength);
    fprintf(fpt4, "# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n", nobj, ncon, nreal, bitlength);
    nbinmut = 0;
    nrealmut = 0;
    nbincross = 0;
    nrealcross = 0;

    start_incr = 0;
    for (i = 0; i < 10; i++)
        active_ref_pts[i] = -1000000;

    parent_pop = (population *)malloc(sizeof(population));
    child_pop = (population *)malloc(sizeof(population));
    mixed_pop = (population *)malloc(sizeof(population));
    allocate_memory_pop(parent_pop, popsize);
    allocate_memory_pop(child_pop, popsize);
    allocate_memory_pop(mixed_pop, 2 * popsize);

    ideal_point = (double *)malloc(nobj * sizeof(double));
    plane_point = (double **)malloc(nobj * sizeof(double *));

    for (i = 0; i < nobj; i++)
        plane_point[i] = (double *)malloc(nobj * sizeof(double));
    for (i = 0; i < NMAX; i++)
        ref_pt[i] = (double *)malloc(nobj * sizeof(double));
    create_ref_points(steps);

    randomize();
    /* Aditional Code CBW */

    //prefix_init_file = std::to_string(prefix);
    //double max_constr = 0.90;
    int pollutant_id = 0;
    //int limit_alpha = 1;
    //int constr_eval = 0;
    //std::string fileName("/home/gtoscano/output/");
    //fileName += prefix_init_file;
    //fileName += "_output.txt";

    read_global_files(emo_uuid.c_str(), pollutant_id, cost_profile_idx);

    int nvars, nconstraints;
    get_info(nvars, nconstraints, fmt::format("{}/output/nsga3/config/", msu_cbpo_path));
    ++nconstraints;
    

    std::cout<<fmt::format("{}/output/nsga3/config/", msu_cbpo_path)<<emo_uuid<<"\n";
    //load_info(nvars, nconstraints, fmt::format("{}/output/nsga3/{}/config/", msu_cbpo_path, emo_uuid), emo_uuid);

    /* End CBW */

  /* gtp: loki
    std::string routing_name = "opt4cast_log_emo";
    std::string msg = fmt::format("{}_{}_{}", "init", emo_uuid, ngen);
    send_message(routing_name, msg);


    routing_name = "opt4cast_log_gen";
    msg = fmt::format("{}_{}_{}", "init", emo_uuid, popsize);
    send_message(routing_name, msg);
    */

    initialize_pop(parent_pop);

    if (mode >0) {
      for (int _i(0); _i < nvars; ++_i) {
        parent_pop->ind[0].xreal[_i] = 0.0;

      }

      CSVReader my_reader(scenario_filename);
      std::vector<std::vector<std::string>> dataList = my_reader.getData();
      int _s, _h, _u, _b;
      double _a;
      for (int _i(0); _i < (int) dataList.size(); ++_i) {
        _s = std::stoi(dataList[_i][0]);
        _h = std::stoi(dataList[_i][1]);
        _u = std::stoi(dataList[_i][2]);
        _b = std::stoi(dataList[_i][3]);
        _a = std::stod(dataList[_i][4]);
        std::string s_tmp = fmt::format("{}_{}_{}_{}", _s, _h, _u, _b);
        std::cout<<"s_tmp is "<<s_tmp<<"\n";
        int idx = get_s_h_u_b(s_tmp);

        if (idx < 0){
          std::cout<<"GTP Error, no s h u b in the file\n";
        }
        parent_pop->ind[0].xreal[idx] = _a;

      }
    }

    std::clog<<"Initialization done, now performing first generation\n";
    curr_gen = 1;

    decode_pop(parent_pop);
    evaluate_pop(parent_pop, curr_gen, ngen, corecast_gen, true, curr_gen);

  /* gtp: loki
    routing_name = "opt4cast_log_gen";
    msg = fmt::format("{}_{}_{}", "terminate", emo_uuid, popsize);
    send_message(routing_name, msg);
    */
      assign_rank(parent_pop);
      find_ideal_point(parent_pop);
      elitist_sorting(parent_pop, parent_pop);

    report_pop(parent_pop, fpt1);
    fprintf(fpt4, "# gen = 1\n");
    report_pop(parent_pop, fpt4);
    std::cout<<"gen = 1\n";
    fflush(stdout);
    /*if (choice!=0)    onthefly_display (parent_pop,gp,1);*/
    fflush(fpt1);
    fflush(fpt2);
    fflush(fpt3);
    fflush(fpt4);
    fflush(fpt5);
    sleep(1);
    bool reevaluate_flag = true;

    for (i = 2; i <= ngen; i++)
    {
        curr_gen = i;

        if (curr_gen >= corecast_gen && reevaluate_flag == true) {
          reevaluate_flag = false;
          evaluate_pop(parent_pop, curr_gen, ngen, corecast_gen, true, curr_gen-1);
        }
        selection(parent_pop, child_pop);
        mutation_pop(child_pop);
        decode_pop(child_pop);

      /* gtp: loki
            routing_name = "opt4cast_log_gen";
            msg = fmt::format("{}_{}_{}", "init", emo_uuid, popsize);
            send_message(routing_name, msg);
            */

        evaluate_pop(child_pop, curr_gen, ngen, corecast_gen, false, curr_gen);

      /* gtp: loki
            routing_name = "opt4cast_log_gen";
            msg = fmt::format("{}_{}_{}", "terminate", emo_uuid, popsize);
            send_message(routing_name, msg);
            */

        merge(parent_pop, child_pop, mixed_pop);
        assign_rank_mixedpop(mixed_pop);
        update_ideal_point(child_pop);
        if (algorithm == "nsga3") {
          elitist_sorting(mixed_pop, parent_pop);
        }
        else {
          fill_nondominated_sort(mixed_pop, parent_pop);
        }

        /* Comment following four lines if information for all
        generations is not desired, it will speed up the execution */
        fprintf(fpt4, "# gen = %d\n", i);
        /* report_pop(parent_pop,fpt4);
         fflush(fpt4);*/
        if (choice != 0 && i == ngen)
            onthefly_display(parent_pop, gp, i);
        onthefly_display2(parent_pop, i);
        std::cout<<fmt::format("gen = {}\n", i);
    }

    if (mode > 0) { // >0 evaluation
        onthefly_display(parent_pop, gp, 1);
        std::string outcome = fmt::format("{}/output/nsga3/{}/front/outcome.txt", msu_cbpo_path, emo_uuid);
        std::ofstream outcome_f(outcome);
        outcome_f.precision(10);
        outcome_f << fmt::format("{:.2f},{:.2f}\n", parent_pop->ind[0].obj[0], parent_pop->ind[0].obj[1]);
        outcome_f.close();
    }

  /* gtp: loki
    routing_name = "opt4cast_log_emo";
    msg = fmt::format("{}_{}_{}", "terminate", emo_uuid, ngen);
    send_message(routing_name, msg);
    */
    //std::cout<<fmt::format("Original count of ref points was : {}\t while new count is :{} and scaling is {}\n", onref, nref, scaling);
    std::cout<<"Generations finished, now reporting solutions\n";
    for (i = 0; i < nref; i++)
    {
        for (j = 0; j < nobj; j++)
            fprintf(fpt6, "%f\t", ref_pt[i][j]);
        fprintf(fpt6, "\n");
    }


    report_pop(parent_pop, fpt2);
    report_feasible(parent_pop, fpt3);
    if (nreal != 0)
    {
        fprintf(fpt5, "\n Number of crossover of real variable = %d", nrealcross);
        fprintf(fpt5, "\n Number of mutation of real variable = %d", nrealmut);
    }
    if (nbin != 0)
    {
        fprintf(fpt5, "\n Number of crossover of binary variable = %d", nbincross);
        fprintf(fpt5, "\n Number of mutation of binary variable = %d", nbinmut);
    }
    fprintf(fpt5, "\n Number of reference points initially = %d", onref);
    fprintf(fpt5, "\n Number of reference points finally = %d", nref);

    fflush(stdout);
    fflush(fpt1);
    fflush(fpt2);
    fflush(fpt3);
    fflush(fpt4);
    fflush(fpt5);
    fflush(fpt6);
    fclose(fpt1);
    fclose(fpt2);
    fclose(fpt3);
    fclose(fpt4);
    fclose(fpt5);
    fclose(fpt6);
    if (choice != 0)
    {
        pclose(gp);
    }
    if (nreal != 0)
    {
        free(min_realvar);
        free(max_realvar);
    }
    if (nbin != 0)
    {
        free(min_binvar);
        free(max_binvar);
        free(nbits);
    }
    for (i = 0; i < nobj; i++)
    {
        free(plane_point[i]);
    }
    for (i = 0; i < NMAX; i++)
    {
        free(ref_pt[i]);
    }
    free(plane_point);
    free(ideal_point);

    deallocate_memory_pop(parent_pop, popsize);
    deallocate_memory_pop(child_pop, popsize);
    deallocate_memory_pop(mixed_pop, 2 * popsize);
    free(parent_pop);
    free(child_pop);
    free(mixed_pop);


    std::cout<<"\n Routine successfully exited \n";
    return (0);
}
