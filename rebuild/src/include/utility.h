#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <ostream>
#include <vector>

#ifndef ROOT_TGraph
#include <TGraph.h>
#endif

#ifndef ROOT_TAxis
#include <TAxis.h>
#endif

#ifndef ROOT_TLegend
#include <TLegend.h>
#endif

#ifndef ROOT_TF1
#include <TF1.h>
#endif

#ifndef ROOT_TCanvas 
#include <TCanvas.h>
#endif

#ifndef ROOT_TH2F
#include <TH2F.h>
#endif

#include <TStyle.h>
#include <TGraphErrors.h>

#ifdef __GSL_STATISTICS_H__
#include <gsl/gsl_statistics_double.h>
#endif

#ifndef MATHEMATICS
#include <mathematics.h>
#endif

/// @details A structure to mane the objects in the class more arrangable.	
/// @author Yiwen Xue
typedef struct TimeSeq{
    /// The size of data.
    int size;
    /// Store the data value.
    double *data;
    /// Store the time. The unit is hour.
    double *time;
} TimeSeq;

typedef struct FitCosinor{
    double *p;
    double chisq;
    double ndf;
} FitCosinor;

int Micedb(std::string infile);

typedef struct TMicemem{
    std::string name;
    std::string path;
    int  batch;
    bool ifmutant;
} TMicemem;

typedef enum RangeIndex{
    rSTART,
    rEND,
    rLENGTH,
} RangeIndex;

int fileLines(std::string path);

class TMice{
    public:
        TMice(TMicemem _mice, int _winSize=90);
        ~TMice();

        TMicemem mice;

        /// This array is used to store the effective period, it is used for the dfa analysis.
        TimeSeq period;
        /// @details This function returns a TimeSeq object which contains the efficient period 
        TimeSeq *GetPeriodEfficient(void ){
            return &this->period;
        }
        /// @details This function is used to draw a figure about how does the period distrubuted.
        int DrawPeriodDist(void );
        int GetPeriod(void );

        int nullize();
        int PeriodDFA(int dfa_order);

        /// The structures to store the data read from file.
        TimeSeq actOri, act, actStd, actAve, actDfa;
        TimeSeq tempOri, temp, tempAve, tempDfa;
        TimeSeq tempRem, actRem, dfaRem;

        int SetRange(int _start, int _length); // public
        /// @details return the range, rSTART means start, rEND means end, rLENGTH means the length.
        int GetDataRange(RangeIndex index); 
        
        /// @details Function to set the fit parameter range, which is super important for nonlinear fit.
        int SetParRange(double *_parMin, double *_parMax);
        /// @details return the Paramenter range, i=0,1 means min or max, j=0,1,2,3 represent p0,p1,p2,p3
        double GetParRange(int i, int j);
        int PrintParRange();

        /// @details Function to Readin the data from file. It will update the mice information, the data, and do a base calculation(Data remap included).
        bool ReadMice(TMicemem _mice);

        int SetwinSize(int _size);
        int GetwinSize();
        int SetwinNum(int _size);
        int GetwinNum();

        /// Remember that in this program the unit of TimeSeq is hour, dont use second any more.
        int SetfitSize(int _size);
        int GetfitSize();
        int SetfitStride(int _step);
        int GetfitStride();

        double GetfitChilimit();
        int SetfitChilimit(double _chi);

        /// @details 
        bool RhythmRemap(); 

        /// @details 
        bool Average();

        /// @details Print the details of this mice.
        int PrintDetails();

        /// @details Draw the heat map. 
        int DrawHeatmap();

        /// draw an over view to select data.
        int DrawOverview();

        /// @details Count the number of errors and voids in original data and calculate the ratio errors divided by total numbers
        double GetErrorPercentage(int );
        int errors_temp;
        int errors_act;

        /// structure to store the period dfa information.
        struct {
            double lbp;
            double index_l;
            double index_r;
            double rbp;
            double midbp;
        } period_dfa;
    private:
        TCanvas *canvas;
        TGraph *grmice;
        double fitChilimit;

        /// A flag to determine if the range was seted.
        int ifRangeSeted; // private
        /// @brief Define the range of data to use. And the unite is hours.
        int start, length; // private 
        int start_index;
        
        double parMax[4]; // private
        double parMin[4]; // private

        /// the number of all windows.
        int winNum;  
        /// size of single window segment (How many points in every window).
        int winSize;

        /// Parameter to determine the fit length. (unit hour)
        int fitSize;
        /// Parameter to determine the Stride of every fit(The shift between every two fit).
        int fitStride; 

    public:
        /* /// used to calculate range of the dfa exponent. */
        /* double dfa_max; */
        /* /// used to calculate range of the dfa exponent. */
        /* double dfa_min; */

        /* double ave_rsq; */
        /* double std_rsq; */
        /* double dfa_rsq; */
        /* void print_data(); */
        /* void print_details(); */
        /* void plot_all(); */
        /* void plot_overview(); */
        /* void plot_ave_rhy(); */
        /* void plot_std_rhy(); */
        /* void plot_dfa_rhy(); */
        /* void plot_heatmap(); */
        /* void plot_powerspec(); */
        /* void mice_info();  // extract informations from the file name */ 
        /* void file_lines(void);  // find the totle lines in the data file */ 
        /* void read_mice_data(void); // read in the data with optimize (Nan, and lack of some points) */
        /* void solve_windows(); // help to solve the windows numbers and windows size */
        /* void solve_init(); // calculate the average value for each window. */
};

bool printMice(TMicemem mice);

int mktable(TMice *mice1, TMice *mice2, double sync_index, std::string batch, bool ifmutant);

int DFA_plot(double *data, int data_size, int dfa_order, std::string name,struct period_dfa *input);

struct period_dfa{
    double lbp;
    double index_l;
    double index_r;
    double rbp;
    double midbp;
};

int loadmice(std::string infile){
    FILE *config = fopen(infile.c_str(), "r");

    char _name[255];
    char _path[255];
    char _mutant[10];
    int batch;
    int start, length;

    std::vector<std::string> script;


    TMicemem mice;
    while(fscanf(config,"%s %s %s %d %d %d\n",_path,_name,_mutant,&batch,&start,&length) == 6){

        struct period_dfa newperiod;

        /* if(batch != 3 && batch != 5) */
        /* if(batch != 3) */
        /* if (name != "12Otx2") */
            /* continue; */
        TMice Mice(mice);
        Mice.nullize();

        std::string path(_path);
        std::string mutant(_mutant);
        std::string name(_name);

        mice.name = name;
        mice.path = path;
        mice.batch = batch;
        mice.ifmutant = (mutant == "mutant"? true:false);

        std::cout << "Next loop" << std::endl;
        Mice.SetRange(0,length*24); 
        std::cout << "Range seted" << std::endl;
        Mice.ReadMice(mice);
        std::cout << "Readed" << std::endl;
        Mice.SetwinSize(90);


        Mice.Average();

        char filename[255];
        sprintf(filename, "../causality/batch%d/%s.Activity.txt", batch, _name);
        printf("filename1 : %s\n", filename);
        FILE *file1 = fopen(filename, "w");
        sprintf(filename, "../causality/batch%d/%s.Temperature.txt", batch, _name);
        printf("filename2 : %s\n", filename);
        FILE *file2 = fopen(filename, "w");

        if(file1 != NULL && file2 != NULL){
            for(int i=0; i< Mice.act.size; i ++)
            {
                fprintf(file1, "%.7lf\t%.7lf\n", Mice.act.time[i], Mice.act.data[i]);
                fprintf(file2, "%.7lf\t%.7lf\n", Mice.temp.time[i], Mice.temp.data[i]);
            }
        } else {
            fprintf(stdout, "../../causality/batch%d/%s.Activity.txt", batch, _name);
            fprintf(stdout, "../../causality/batch%d/%s.Activity.txt", batch, _name);
        }

        fclose(file1);
        fclose(file2);

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
        /* Mice.PeriodDFA(1); */
        /* script.push_back(mice.name +"\tDFA1\t"+ std::to_string(Mice.period_dfa.lbp) + "\t" + std::to_string(Mice.period_dfa.rbp) + "\t" + std::to_string(Mice.period_dfa.index_l) + "\t" + std::to_string(Mice.period_dfa.index_r)); */
        /* Mice.PeriodDFA(2); */
        /* /1* script.push_back(mice.name +"\tDFA2\t"+ std::to_string(Mice.period_dfa.lbp) + "\t" + std::to_string(Mice.period_dfa.rbp) + "\t" + std::to_string(Mice.period_dfa.index_l) + "\t" + std::to_string(Mice.period_dfa.index_r)); *1/ */
        /* Mice.PeriodDFA(3); */
        /* script.push_back(mice.name +"\tDFA3\t"+ std::to_string(Mice.period_dfa.lbp) + "\t" + std::to_string(Mice.period_dfa.rbp) + "\t" + std::to_string(Mice.period_dfa.index_l) + "\t" + std::to_string(Mice.period_dfa.index_r)); */
        /* Mice.PeriodDFA(4); */
        /* /1* script.push_back(mice.name +"\tDFA4\t"+ std::to_string(Mice.period_dfa.lbp) + "\t" + std::to_string(Mice.period_dfa.rbp) + "\t" + std::to_string(Mice.period_dfa.index_l) + "\t" + std::to_string(Mice.period_dfa.index_r)); *1/ */
        /* Mice.RhythmRemap(); */


        for (int i = 1; i < 6; ++i) {
            if (i != 2)
                continue;
            /* DFA_plot(Mice.actAve.data, Mice.actAve.size, i, mice.name + "_activity_" + std::to_string(i), &newperiod); */
            /* script.push_back(mice.name +"\tDFA" + std::to_string(i)+ "\tactivity\t"+ std::to_string(newperiod.lbp) + "\t" + std::to_string(newperiod.rbp) + "\t" + std::to_string(newperiod.index_l) + "\t" + std::to_string(newperiod.index_r)); */
            DFA_plot(Mice.tempAve.data, Mice.tempAve.size, i, mice.name + "_temperature_" + std::to_string(i), &newperiod);
            script.push_back(mice.name +"\tDFA" + std::to_string(i)+ "\ttemperature\t"+ std::to_string(newperiod.lbp) + "\t" + std::to_string(newperiod.rbp) + "\t" + std::to_string(newperiod.index_l) + "\t" + std::to_string(newperiod.index_r));
        }

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


double fit_crossover(double *x, double *par);
double draw_crossover(double *x, double *par);
