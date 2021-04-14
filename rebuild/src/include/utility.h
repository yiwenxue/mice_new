#include <Rtypes.h>
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

        /// Draw the heatmap of a time series.

        /// @details Count the number of errors and voids in original data and calculate the ratio errors divided by total numbers
        double GetErrorPercentage(int );
        int errors_temp;
        int errors_act;

        double periodl;
        double periodr;
        double periodall;

        /// structure to store the period dfa information.
        struct {
            double lbp;
            double index_l;
            double index_r;
            double rbp;
            double midbp;
        } period_dfa;
        double *residuals_temp;
        double *residuals_act;
        void calculate_residual();
        void residuals_cleanup();
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

int DFA_dual_plot(double *data, int data_size, int dfa_order, std::string name,struct period_dfa *input);

int DFA_plot(double *data, int data_size, int dfa_order, std::string name,struct period_dfa *input, double crossover = 0);

int DFA_filter_plot(double *data, int data_size, int dfa_order, std::string name,struct period_dfa *input, TimeSeq dataRem);

struct period_dfa{
    double lbp;
    double index_l;
    double index_r;
    double rbp;
    double midbp;
};

void data_overview(std::string name, TimeSeq seriesx, TimeSeq seriesy, int rangel = 0, int ranger = 0){
    auto cv = new TCanvas("super cv");
    cv->SetCanvasSize(480, 270);
    auto pad = new TPad("pad", "", 0, 0, 1 ,1);

    TGraph *gr1, *gr2;

    if (ranger == 0){
        gr1 = new TGraph(seriesx.size, seriesx.time, seriesx.data);
        gr2 = new TGraph(seriesy.size, seriesy.time, seriesy.data);
    }else if (ranger > seriesx.size || ranger > seriesy.size || rangel < 0){
        printf("[!!] On plotting overview: out of range\n");
        return;
    } else {
        gr1 = new TGraph(ranger-rangel, seriesx.time+rangel, seriesx.data+rangel);
        gr2 = new TGraph(ranger-rangel, seriesy.time+rangel, seriesy.data+rangel);
    }
    gr1->SetTitle("Act");
    gr2->SetTitle("Temp");

    double *x;
    int N;
    x = gr1->GetX();
    N = gr1->GetN();
    for(int i=0;i<N;i++){
        x[i] /= 24;
    }
    /* gr1->GetXaxis()->SetRangeUser(x[0], x[N-1]); */
    gr1->GetXaxis()->SetRange(0, N-1);
    gr1->GetYaxis()->SetTickLength(0.01);
    x = gr2->GetX();
    N = gr2->GetN();
    for(int i=0;i<N;i++){
        x[i] /= 24;
    }
    /* gr2->GetXaxis()->SetRangeUser(x[0], x[N-1]); */
    gr2->GetYaxis()->SetTickLength(0.01);
    gr2->GetXaxis()->SetRange(0, N-1);

    gr1->GetXaxis()->SetTitle("Time (day)");
    gr1->GetXaxis()->CenterTitle(1);
    gr2->GetXaxis()->SetTitle("Time (day)");
    gr2->GetXaxis()->CenterTitle(1);

    cv->SetFillStyle(4000);
    pad->Divide(1,2, 0,0);
    pad->SetFillStyle(4000);
    pad->Draw();
    pad->cd(1);
    gr1->Draw();
    pad->cd(2);
    gr2->Draw();

    cv->Update();
    cv->Print(name.c_str(), "Title:title");

    delete gr1;
    delete gr2;
    delete pad;
    delete cv;
}

void residuals_test(){
    TMicemem micemem;
    micemem.name = "12Otx2";
    micemem.batch = 3;
    micemem.ifmutant = true;
    micemem.path = "../../data/";

    TMice mice(micemem);
    mice.nullize();
    mice.SetRange(0, 36*24);
    mice.ReadMice(micemem);
    mice.SetwinSize(90);
    mice.Average();

    auto cv = new TCanvas("cv");
    cv->cd();

    auto grmice = new TGraph(mice.tempAve.size, mice.tempAve.time, mice.tempAve.data);
    grmice->GetXaxis()->SetRangeUser(300, 348);
    grmice->Draw();

    auto fit1 = new TF1("fit1", "[0]*sin([1]*x+[2])+[3]", 300, 348);
    fit1->SetParLimits(0, 0, 10);
    fit1->SetParLimits(1, 0.20, 0.32);
    fit1->SetParLimits(2, -24, 24);
    /* fit1->SetParLimits(3, 0, 20); */
    fit1->SetParLimits(3, 30, 50);
    grmice->Fit(fit1, "R");

    fit1->Draw("SAME");

    int start_idx, end_idx=0; 
    for (start_idx = end_idx; start_idx < mice.tempAve.size; ++start_idx) {
       if (mice.tempAve.time[start_idx] >= 300)
           break;
    }
    for (end_idx = start_idx; end_idx < mice.tempAve.size; ++end_idx) {
       if (mice.tempAve.time[end_idx+1] >= 348)
           break;
    }

    printf("time[%d,%d] = [%f,%f]\n", start_idx, end_idx, mice.tempAve.time[start_idx],
            mice.tempAve.time[end_idx]);

    int resize = end_idx - start_idx + 1;
    double *residuals = new double[resize];
    for (int i = 0; i < resize; ++i) {
        residuals[i] = mice.tempAve.data[i + start_idx] - fit1->Eval(mice.tempAve.time[i + start_idx]) + 35;
    }

    auto gresidual = new TGraph(resize, mice.tempAve.time+start_idx, residuals);
    gresidual->SetLineColor(kGreen);
    gresidual->Draw("SAME");
}



double fit_crossover(double *x, double *par);
double fit_crossover_dual(double *x, double *par);
double dfa_functions(double *x, double *par);
double draw_crossover(double *x, double *par);

struct dfa_data {
    double scal;
    double fluc;
};

enum mice_type {
    mice_type_control = 0x01,
    mice_type_mutant = 0x02,
    mice_type_mutant_hyper = 0x04,
};

void normalize_double(double *data, int size);

void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny,
                     Float_t lMargin, Float_t rMargin,
                     Float_t bMargin, Float_t tMargin);

void get_waveform(double *data, int size, double *waveform, int len);