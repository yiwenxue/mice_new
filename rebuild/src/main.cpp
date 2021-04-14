#include <cstdio>
#include <cstdlib>
#include <gsl/gsl_statistics_double.h>
#include <iostream>
#include <ostream>
#include <random>
#include <string>
#include "utility.h"

#define MIN(x,y) (((x)>(y))?(y):(x))

using namespace std;

int loadmice(std::string infile){
    FILE *config = fopen(infile.c_str(), "r");

    if (config == NULL){
        printf("Error!");
        exit(-1);
    }

    char _name[255];
    char _path[255];
    char _mutant[10];
    int batch;
    int n,m;
    int start, length;

    std::vector<std::string> script;

    TMicemem mice;
    while(fscanf(config,"%s %s %s %d %d %d %d %d\n",_path,_name,_mutant,&batch,&start,&length,&n,&m) == 8){
        if (_path[0] == '#')
            continue;

        struct period_dfa newperiod;

        std::string path(_path);
        std::string mutant(_mutant);
        std::string name(_name);

        mice.name = name;
        mice.path = path;
        mice.batch = batch;
        mice.ifmutant = (mutant == "mutant"? true:false);

        TMice Mice(mice);
        Mice.nullize();
        
        /* std::cout << "Next loop" << std::endl; */
        Mice.SetRange(start*24,length*24); 
        /* std::cout << "Range seted" << std::endl; */
        Mice.ReadMice(mice);
        /* std::cout << "Readed" << std::endl; */
        Mice.SetwinSize(90);

        Mice.Average();

        /* Draw the overview of data | Only for debug. */
        /* TCanvas c1(""); */
        /* c1.cd(); */
        /* TGraph gr_data_ave_temp(length * 24 * (360/Mice.GetwinSize()), Mice.tempAve.time, Mice.tempAve.data); */
        /* gr_data_ave_temp.GetXaxis()->SetTitle("time(hour)"); */
        /* gr_data_ave_temp.GetXaxis()->CenterTitle(); */
        /* gr_data_ave_temp.GetYaxis()->SetTitle("Temperature(C)"); */
        /* gr_data_ave_temp.GetYaxis()->CenterTitle(); */
        /* gr_data_ave_temp.Draw(); */

        /* c1.Update(); */
        /* c1.Print(("data_ave_temp" + mice.name+ ".pdf").c_str()); */
        /* c1.Clear(); */
        /* TGraph gr_data_ave_act(length * 24 * (360/Mice.GetwinSize()), Mice.actAve.time, Mice.actAve.data); */
        /* gr_data_ave_act.GetXaxis()->SetTitle("time(hour)"); */
        /* gr_data_ave_act.GetXaxis()->CenterTitle(); */
        /* gr_data_ave_act.GetYaxis()->SetTitle("Temperature(C)"); */
        /* gr_data_ave_act.GetYaxis()->CenterTitle(); */
        /* gr_data_ave_act.Draw(); */
        /* c1.Update(); */
        /* c1.Print(("data_ave_act" + mice.name + ".pdf").c_str()); */
        /* c1.Clear(); */

        /* int lenth = 10; */
        /* data_overview(mice.name + "MA.pdf", Mice.tempAve, Mice.actAve, 0, 1000); */

        /* Sync mice_sync(Mice.actAve.data + start*24*4, Mice.tempAve.data+start*24*4, length*24*4, n, m); */

        /* mice_sync.PSI_decay(mice.name, 4 * 24 * 15, 1); */

        /* mice_sync.plot((mice.name + "sync_dist_" + std::to_string(90)), 90); */
        /* mice_sync.PSI_rho(0); */
        /* /1* mice_sync.PSI_lambda(0); *1/ */

        /* /1* if (batch == 3) *1/ */
        /*     mice_sync.phasePlot(mice.name); */
        /*     /1* mice_sync.phasePlot(mice.name, 1, 1); *1/ */
        /* /1* else if (batch == 5) *1/ */
        /*     /1* mice_sync.phasePlot(mice.name, 2, 3); *1/ */
        /* /1* else if (batch == 6) *1/ */
        /*     /1* mice_sync.phasePlot(mice.name, 3, 2); *1/ */

        /* FILE *phase_out = fopen((mice.name + "phase.txt").c_str(), "w"); */
        /* for (int i = 0; i < mice_sync.size; ++i) { */
        /*     double x = 0; */
        /*     double y = 0; */
        /*     int integer = 0; */
        /*     fprintf(phase_out, "%8.3lf\t%8.4lf\n", mice_sync.phasex[i], mice_sync.phasey[i]); */
        /* } */
        /* fclose(phase_out); */

        /* /1* script.push_back(mice.name + "\t" + std::to_string(mice_sync.rho_index) + "\t" + std::to_string(mice_sync.lambda_index) + "\t" + std::to_string(mice_sync.decay_radio)); *1/ */
        /* script.push_back(mice.name + "\t" + std::to_string(mice_sync.rho_index) + "\t" + std::to_string(mice_sync.decay_radio)); */

        /* char filename[255]; */
        /* sprintf(filename, "../causality/batch%d/%s.Activity.txt", batch, _name); */
        /* printf("filename1 : %s\n", filename); */
        /* FILE *file1 = fopen(filename, "w"); */
        /* sprintf(filename, "../causality/batch%d/%s.Temperature.txt", batch, _name); */
        /* printf("filename2 : %s\n", filename); */
        /* FILE *file2 = fopen(filename, "w"); */

        /* int tstart = start * 24 * 360; */
        /* int tend = tstart + length * 24 * 360; */

        /* if(file1 != NULL && file2 != NULL){ */
        /*     for(int i=tstart; i< tend; i ++) */
        /*     { */
        /*         fprintf(file1, "%.7lf\t%.7lf\n", Mice.act.time[i], Mice.act.data[i]); */
        /*         fprintf(file2, "%.7lf\t%.7lf\n", Mice.temp.time[i], Mice.temp.data[i]); */
        /*     } */
        /* } else { */
        /*     fprintf(stdout, "../../causality/batch%d/%s.Activity.txt", batch, _name); */
        /*     fprintf(stdout, "../../causality/batch%d/%s.Activity.txt", batch, _name); */
        /* } */

        /* fclose(file1); */
        /* fclose(file2); */

        /* std::cout << "Averaged" << std::endl; */
        double Parmax[4] = {5, 0.32, 24, 50};
        double Parmin[4] = {0, 0.20, -24, 30};
        Mice.SetfitSize(48);
        Mice.SetfitStride(4);
        Mice.SetParRange(Parmin,Parmax);
        Mice.PrintDetails();
        /* std::cout << "Details printed" << Mice.tempAve.size<< std::endl; */

        Mice.GetPeriod();
        /* /1* std::cout << "Period got" << std::endl; *1/ */
        /* /1* Mice.DrawOverview(); *1/ */
        /* Mice.DrawPeriodDist(); */
        /* script.push_back(mice.name +"\tCircadian Rhythm\t" + std::to_string(Mice.periodl) + "\t" + std::to_string(Mice.periodr) + "\t" + std::to_string(Mice.periodall)); */
        /* Mice.PeriodDFA(1); */
        /* script.push_back(mice.name +"\tDFA1\t"+ std::to_string(Mice.period_dfa.lbp) + "\t" + std::to_string(Mice.period_dfa.rbp) + "\t" + std::to_string(Mice.period_dfa.index_l) + "\t" + std::to_string(Mice.period_dfa.index_r)); */
        Mice.PeriodDFA(2);
        /* /1* script.push_back(mice.name +"\tDFA2\t"+ std::to_string(Mice.period_dfa.lbp) + "\t" + std::to_string(Mice.period_dfa.rbp) + "\t" + std::to_string(Mice.period_dfa.index_l) + "\t" + std::to_string(Mice.period_dfa.index_r)); *1/ */
        Mice.PeriodDFA(3);
        /* script.push_back(mice.name +"\tDFA3\t"+ std::to_string(Mice.period_dfa.lbp) + "\t" + std::to_string(Mice.period_dfa.rbp) + "\t" + std::to_string(Mice.period_dfa.index_l) + "\t" + std::to_string(Mice.period_dfa.index_r)); */
        /* Mice.PeriodDFA(4); */
        /* /1* script.push_back(mice.name +"\tDFA4\t"+ std::to_string(Mice.period_dfa.lbp) + "\t" + std::to_string(Mice.period_dfa.rbp) + "\t" + std::to_string(Mice.period_dfa.index_l) + "\t" + std::to_string(Mice.period_dfa.index_r)); *1/ */
        Mice.RhythmRemap();

        /* Mice.calculate_residual(); */

        /* std::cout << "The evil part" << std::endl; */

        for (int i = 1; i < 5; ++i) {
            if (i != 2 && i != 3)
                continue;
            // DFA with freq filter
            DFA_filter_plot(Mice.actAve.data, length * 24 * (360/Mice.GetwinSize()), i, mice.name + " activity", &newperiod, Mice.actRem);
            DFA_filter_plot(Mice.tempAve.data, length * 24 * (360/Mice.GetwinSize()), i, mice.name + " temperature", &newperiod, Mice.tempRem);

            /* DFA_plot(Mice.residuals_act, length * 24 * (360/Mice.GetwinSize()), i, mice.name + " act_detrend", &newperiod,0); */

            DFA_plot(Mice.actAve.data, length * 24 * (360/Mice.GetwinSize()), i, mice.name + " activity_1day", &newperiod, 1);
            /* DFA_dual_plot(Mice.actAve.data, length * 24 * (360/Mice.GetwinSize()), i, mice.name + " activity", &newperiod); */

            /* DFA_plot(Mice.residuals_temp, length * 24 * (360/Mice.GetwinSize()), i, mice.name + " temp_detrend", &newperiod,0); */
            /* DFA_plot(Mice.residuals_act, length * 24 * (360/Mice.GetwinSize()), i, mice.name + " act_detrend", &newperiod,0); */
            script.push_back(mice.name +"\tDFA" + std::to_string(i)+ "\tactivity 1day\t"+ std::to_string(newperiod.lbp) + "\t" + std::to_string(newperiod.rbp) + "\t" + std::to_string(newperiod.index_l) + "\t" + std::to_string(newperiod.index_r));

            DFA_plot(Mice.tempAve.data, length * 24 * (360/Mice.GetwinSize()), i, mice.name + " temperature_1day", &newperiod, 1);
            /* DFA_dual_plot(Mice.tempAve.data, length * 24 * (360/Mice.GetwinSize()), i, mice.name + " temperature", &newperiod); */
            script.push_back(mice.name +"\tDFA" + std::to_string(i)+ "\ttemperature 1day\t"+ std::to_string(newperiod.lbp) + "\t" + std::to_string(newperiod.rbp) + "\t" + std::to_string(newperiod.index_l) + "\t" + std::to_string(newperiod.index_r));
        }

        /* Mice.residuals_cleanup(); */

        /* std::cout << "Remaped" << std::endl; */
        /* Mice.DrawHeatmap(); */
        /* std::cout << "Drawn" << std::endl; */
        /* Mice.GetErrorPercentage(0); */
    }

    for (std::string a: script){
        std::cout << a << std::endl;
    }
    fclose(config);
    return 0;
}

int main(int argc, char *argv[])
{
    loadmice("./mice");
    return 0;
}
