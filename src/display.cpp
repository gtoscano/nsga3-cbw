/* Routines to display the population information using gnuplot */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <iostream>
#include <string>
#include <fmt/core.h>
#include <filesystem>


#include "global.h"
#include "rand.h"
#include <fstream>

#include "misc.hpp"

extern int mode;
namespace fs = std::filesystem;

std::vector<double> true_pf_lb;
std::vector<double> true_pf_ub;
double hv_raw_pf;
double pct_hv_to_exit;
/* Function to display the current population for the subsequent generation */
void onthefly_display(population *pop, FILE *gp, int ii)
{
    int i, j;
    int flag;
    FILE *fpt, *fp;

    std::string msu_cbpo_path = getEnvVar("MSU_CBPO_PATH", "/opt/opt4cast");
    std::string stored_path = fmt::format("{}/output/nsga3/{}/front/",
                                          msu_cbpo_path, emo_uuid);
    std::string specific_path = stored_path + "plot.out";
    fpt = fopen(specific_path.c_str(), "w");
    specific_path = stored_path + "pareto_front.out";
    fp = fopen(specific_path.c_str(), "a");
    specific_path = stored_path + "pareto_front_uuid.out";
    std::ofstream ofile(specific_path);
//    specific_path = stored_path + "pareto_front_indexes_fx.out";
 //   std::ofstream pf_file(specific_path);

    specific_path = fmt::format("{}{}", stored_path,"pareto_front.out");
    //fp_gen = fopen(specific_path.c_str(), "w");
    /*fprintf(fp,"#\n");*/
    flag = 0;
    std::string dir_path = fmt::format("{}/output/nsga3/{}/front",
                                       msu_cbpo_path, emo_uuid);
    if(!fs::exists(dir_path)) //Does not exist
        fs::create_directories(dir_path);

    specific_path = stored_path + "exec_uuid_registry.out";
    std::ofstream ofile_registry(specific_path);
    for (auto const&[key, value] : uuid_registry)
        ofile_registry << key << " " << value << std::endl;
    ofile_registry.close();

    specific_path = stored_path + "uuid_surviving_pop_registry.out";
    std::ofstream ofile_surviving_pop_registry(specific_path);
    for (auto const&[key, value] : uuid_surviving_pop_registry)
        ofile_surviving_pop_registry << key << " " << value << std::endl;
    ofile_surviving_pop_registry.close();

    for (i = 0; i < popsize; i++) {
        if (pop->ind[i].constr_violation >= 0.0 && pop->ind[i].rank == 1 )
        //if (pop->ind[i].constr_violation == 0)
        {
            auto uuid_registry_idx = fmt::format("{}_{}", ii, i);
            std::string tmp_obj = fmt::format("{}", pop->ind[i].obj[0]);
            for (int j(1); j<nobj; ++j){
                tmp_obj = fmt::format("{},{}",tmp_obj, pop->ind[i].obj[j]);
            }
            std::string exec_uuid;

            if (ii == 1)
                exec_uuid = uuid_registry[uuid_registry_idx];
            else
                exec_uuid = uuid_surviving_pop_registry[uuid_registry_idx];


            ofile<<fmt::format("{},{},{}\n", i, exec_uuid, tmp_obj);
                std::string filename_src = fmt::format("{}/output/nsga3/{}/{}_output_t.csv",
                                                       msu_cbpo_path, emo_uuid, exec_uuid);
                std::string filename_dst = fmt::format("{}/output/nsga3/{}/front/{}_output_t.csv",
                                                       msu_cbpo_path, emo_uuid, i);
                if(fs::exists(filename_src)) //Does it exist
                    fs::copy(filename_src, filename_dst, fs::copy_options::update_existing);
                filename_src = fmt::format("{}/output/nsga3/{}/{}_output_t.txt",
                                           msu_cbpo_path, emo_uuid, exec_uuid);
                filename_dst = fmt::format("{}/output/nsga3/{}/front/{}_output_t.txt",
                                           msu_cbpo_path, emo_uuid, i);
                if(fs::exists(filename_src)) //Does it exist
                    fs::copy(filename_src, filename_dst, fs::copy_options::update_existing);
                filename_src = fmt::format("{}/output/nsga3/{}/{}_reportloads.csv",
                                           msu_cbpo_path, emo_uuid, exec_uuid);
                filename_dst = fmt::format("{}/output/nsga3/{}/front/{}_reportloads.csv",
                                           msu_cbpo_path, emo_uuid, i);
                if(fs::exists(filename_src)) //Does it exist
                    fs::copy(filename_src, filename_dst, fs::copy_options::update_existing);
                filename_src = fmt::format("{}/output/nsga3/{}/{}_reportloads.parquet",
                                           msu_cbpo_path, emo_uuid, exec_uuid);
                filename_dst = fmt::format("{}/output/nsga3/{}/front/{}_reportloads.parquet", msu_cbpo_path,
                                           emo_uuid, i);
                if(fs::exists(filename_src)) //Does it exist
                    fs::copy(filename_src, filename_dst, fs::copy_options::update_existing);
                filename_src = fmt::format("{}/output/nsga3/{}/{}_impbmpsubmittedland.parquet",
                                           msu_cbpo_path, emo_uuid, exec_uuid);
                filename_dst = fmt::format("{}/output/nsga3/{}/front/{}_impbmpsubmittedland.parquet",
                                           msu_cbpo_path, emo_uuid, i);
                if(fs::exists(filename_src)) //Does it exist
                    fs::copy(filename_src, filename_dst, fs::copy_options::update_existing);
            filename_src = fmt::format("{}/output/nsga3/{}/{}_impbmpsubmittedanimal.parquet",
                                       msu_cbpo_path, emo_uuid, exec_uuid);
            filename_dst = fmt::format("{}/output/nsga3/{}/front/{}_impbmpsubmittedanimal.parquet",
                                       msu_cbpo_path, emo_uuid, i);
            if(fs::exists(filename_src)) //Does it exist
                fs::copy(filename_src, filename_dst, fs::copy_options::update_existing);
                filename_src = fmt::format("{}/output/nsga3/{}/{}_bmp_summary_t.txt",
                                           msu_cbpo_path, emo_uuid, exec_uuid);
                filename_dst = fmt::format("{}/output/nsga3/{}/front/{}_bmp_summary_t.txt",
                                           msu_cbpo_path, emo_uuid, i);
                if(fs::exists(filename_src)) //Does it exist
                    fs::copy(filename_src, filename_dst, fs::copy_options::update_existing);

                //get_outcome_by_lrs(emo_uuid, exec_uuid);
                filename_src = fmt::format("{}/output/nsga3/{}/{}_outcome.csv", msu_cbpo_path, emo_uuid, exec_uuid);
                filename_dst = fmt::format("{}/output/nsga3/{}/front/outcome_lrs.csv", msu_cbpo_path, emo_uuid);
                if(fs::exists(filename_src)) //Does it exist
                    fs::copy(filename_src, filename_dst, fs::copy_options::update_existing);


            if (nobj < 4)
            {
                if (choice != 3)
                {
                    if (pop->ind[i].constr[0] >= 0) {
                        fprintf(fpt, "%e\t%e\n", pop->ind[i].obj[obj1 - 1], pop->ind[i].obj[obj2 - 1]);
                        fprintf(fp, "%e\t%e\n", pop->ind[i].obj[obj1 - 1], pop->ind[i].obj[obj2 - 1]);
                        //fprintf(fp_gen, "%e\t%e\n", pop->ind[i].obj[obj1 - 1], pop->ind[i].obj[obj2 - 1]);
                    }
                    else{
                     std::cout<<pop->ind[i].constr[0]<<"\n";

                    }
                }
                else
                {
                    fprintf(fpt, "%e\t%e\t%e\n", pop->ind[i].obj[obj1 - 1], pop->ind[i].obj[obj2 - 1], pop->ind[i].obj[obj3 - 1]);
                    fprintf(fp, "%e\t%e\t%e\n", pop->ind[i].obj[obj1 - 1], pop->ind[i].obj[obj2 - 1], pop->ind[i].obj[obj3 - 1]);
                    //fprintf(fp_gen, "%e\t%e\t%e\n", pop->ind[i].obj[obj1 - 1], pop->ind[i].obj[obj2 - 1], pop->ind[i].obj[obj3 - 1]);
                }
                fflush(fpt);
                fflush(fp);
                //fflush(fp_gen);
                flag = 1;
            }
            else
            {
                for (j = 0; j < nobj; j++)
                {
                    fprintf(fpt, "%e\t", pop->ind[i].obj[j]);
                    fprintf(fp, "%e\t", pop->ind[i].obj[j]);
                    //fprintf(fp_gen, "%e\t", pop->ind[i].obj[j]);
                }
                fprintf(fpt, "\n");
                fprintf(fp, "\n");
                //fprintf(fp_gen, "\n");
            }
        }
    }
    fprintf(fp, "#\n");
    if (flag == 0 || nobj > 3)
    {
        printf("\n No feasible soln in this pop, hence no display");
    }
    else
    {
      /*
        if (choice != 3)
            fprintf(gp, "set title 'Generation #%d'\n unset key\n plot 'plot.out' w points pointtype 6 pointsize 1\n", ii);
        else
            fprintf(gp, "set title 'Generation #%d'\n set view %d,%d\n unset key\n splot 'plot.out' w points pointtype 6 pointsize 1\n set ticslevel 0\n", ii, angle1, angle2);
            */
        fflush(gp);
    }
    fclose(fpt);
    fclose(fp);
    ofile.close();
    //fclose(fp_gen);
    return;
}
void onthefly_display2(population *pop, int ii)
{
    int i, j;
    int flag;
    FILE *fp;

    std::string msu_cbpo_path = getEnvVar("MSU_CBPO_PATH", "/opt/opt4cast");
    std::string stored_path = fmt::format("{}/output/nsga3/{}/front/gen_{}_sol.out", msu_cbpo_path, emo_uuid, ii);
    fp = fopen(stored_path.c_str(), "a");

    flag = 0;
    for (i = 0; i < popsize; i++)
    {
      if(true)
        //if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank == 1)
        //if (pop->ind[i].constr_violation == 0)
        {
            if (nobj < 4)
            {
                if (choice != 3)
                {
                    if (pop->ind[i].constr[0] >= 0) {
                        fprintf(fp, "%e\t%e\n", pop->ind[i].obj[obj1 - 1], pop->ind[i].obj[obj2 - 1]);
                    }
                }
                else
                {
                    fprintf(fp, "%e\t%e\t%e\n", pop->ind[i].obj[obj1 - 1], pop->ind[i].obj[obj2 - 1], pop->ind[i].obj[obj3 - 1]);
                }
                fflush(fp);
                flag = 1;
            }
            else
            {
                for (j = 0; j < nobj; j++)
                {
                    fprintf(fp, "%e\t", pop->ind[i].obj[j]);
                }
                fprintf(fp, "\n");
            }
        }
    }
    if (flag == 0 || nobj > 3)
    {
        printf("\n No feasible soln in this pop, hence no display");
    }
    else
    {
      /*
        if (choice != 3)
            fprintf(gp, "set title 'Generation #%d'\n unset key\n plot 'plot.out' w points pointtype 6 pointsize 1\n", ii);
        else
            fprintf(gp, "set title 'Generation #%d'\n set view %d,%d\n unset key\n splot 'plot.out' w points pointtype 6 pointsize 1\n set ticslevel 0\n", ii, angle1, angle2);
            */
        //fflush(gp);
    }
    fclose(fp);
    return;
}

/* Function to display the current population for the subsequent generation */
int onthefly_display2(population *pop)
{
    int i, j;
    FILE *fpt; 
    fpt = fopen("current_pareto_front.out", "w");
    /*fprintf(fp,"#\n");*/
    fprintf(fpt,"#\n");
    int counter= 0;
    for (i = 0; i < popsize; i++)
    {
      
        if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank == 1)
        {
            if (nobj < 4)
            {
                if (choice != 3)
                {

                    double norm_obj_1 = (pop->ind[i].obj[obj1 - 1] - true_pf_lb[obj1-1]) / (true_pf_ub[obj1-1]-true_pf_lb[obj1-1]);
                    double norm_obj_2 = (pop->ind[i].obj[obj2 - 1] - true_pf_lb[obj2-1]) / (true_pf_ub[obj2-1]-true_pf_lb[obj2-1]);
                    
                    if (pop->ind[i].constr[0] >= 0) {
                        if ( pop->ind[i].obj[obj1-1] <= true_pf_ub[obj1-1] && pop->ind[i].obj[obj2 - 1] <= true_pf_ub[obj2-1]){
                            fprintf(fpt, "%e %e\n", norm_obj_1, norm_obj_2);
                            ++counter;
                        }
                    }
                }
                else
                {
                    double norm_obj_1 = (pop->ind[i].obj[obj1 - 1] - true_pf_lb[obj1-1]) / (true_pf_ub[obj1-1]-true_pf_lb[obj1-1]);
                    double norm_obj_2 = (pop->ind[i].obj[obj2 - 1] - true_pf_lb[obj2-1]) / (true_pf_ub[obj2-1]-true_pf_lb[obj2-1]);
                    double norm_obj_3 = (pop->ind[i].obj[obj3 - 1] - true_pf_lb[obj3-1]) / (true_pf_ub[obj3-1]-true_pf_lb[obj3-1]);
                    if (pop->ind[i].obj[obj1 - 1] <= true_pf_ub[obj1 - 1] && pop->ind[i].obj[obj2 - 1] <= true_pf_ub[obj2 - 1] && pop->ind[i].obj[obj3 - 1] <= true_pf_ub[obj3-1]){
                        fprintf(fpt, "%e %e %e\n", norm_obj_1, norm_obj_2, norm_obj_3);
                        ++counter;
                    }
                }
                fflush(fpt);
            }
            else
            {
                bool flag = false;
                for (j = 0; j < nobj; j++)
                {
                    if (pop->ind[i].obj[j] > true_pf_ub[j]) {
                        flag = true;
                        break;
                    }
                }
                if (flag == false){
                for (j = 0; j < nobj; j++)
                {
                    double norm_obj = (pop->ind[i].obj[j] - true_pf_lb[j]) / (true_pf_ub[j]-true_pf_lb[j]);
                    fprintf(fpt, "%e ", norm_obj);
                    ++counter;
                }
                fprintf(fpt, "\n");
            }
            }
        }
    }
    fprintf(fpt,"#\n");
    fclose(fpt);
    return counter;
}
