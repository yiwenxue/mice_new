// C++ headers
#include <random>
#include <iostream>
#include <ctime>
#include <fftw3.h>
#include <cmath>
#include <string>


// ROOT headers
#include <RtypesCore.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h> 
#include <TGaxis.h> 
#include <RConfigure.h>
#include <TColor.h>
#include <TROOT.h>
#include <TMultiGraph.h>
#include <TLegend.h>

#ifndef __GSL_HISTOGRAM_H_
#include <gsl/gsl_histogram.h>
#endif

#ifndef __GSL_STATISTICS_H__
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_statistics_double.h>
#endif 

#ifndef __GSL_MULTIFIT_H__
#include <gsl/gsl_multifit.h>
#endif

#ifndef  __GSL_SF_H__ 
#include <gsl/gsl_sf.h>
#endif 

#ifndef GNUPLOT 
#include <plot.h>
#endif

#ifndef MATHEMATICS 
#define MATHEMATICS 
#endif 


#define NMICE 1000
#define RE 0
#define IM 1

/// Mathematical function for least squre polynominal fit 
/// @details This class is will perform a polynominal fit (for any order)
class Polyfit{
    public:
        Polyfit(double *datax, double *datay, int _size, int _order);
        virtual ~Polyfit();
        double fitfunc(double input);
        inline int get_degree(){
            return degree;
        }
        inline double get_rsq(){
            return coeff;
        }
        void resize(int new_size);
        void init(double *datax, double *datay);
        void recal(double *datax, double *datay); // size of the new seq should have the same size as the setting one;
    public:
        std::string func; 
        double *store;
    private:
        gsl_multifit_linear_workspace *ws;
        gsl_matrix *cov, *X;
        gsl_vector *y, *c;
        int size;
        int degree;
        double coeff;
        double meanx;
        double tss;
        int i, j;
    private:
        bool polyfit();
};

class Cosinor{
    public:
        Cosinor(double *x, double *y, int size, int degree);
        virtual ~Cosinor();
        double fitfunc(double input);
        inline int get_degree()
        {
            return degree;
        }
        inline double get_rsq()
        {
            return coeff;
        }
        void init(double *datax, double *datay);
        void recal(double *datax, double *datay);
    public:
        std::string func; 
    private:
        double omega;
        int size;
        int degree;
        double *store;
        double coeff;
        double tss;

        gsl_multifit_linear_workspace *ws;
        gsl_matrix *X, *cov;
        gsl_vector *Y, *c;
        int i,j;
    private:
        bool cosinor();
};


class DFA{
    public:
        DFA(double *data, int data_size, int dfa_order);
        virtual ~DFA();

        inline double get_rsq()
        {
            return coeff;
        }
        inline int get_order()
        {
            return order;
        }
        inline int get_size_o()
        {
            return size_o;
        }
        void init(double *data);        // init function;
        void recal(double *data); // reuse object without delete;
        double fitfunc(double input); // the final function; 
        void plot_over(int point); // plot out 
    public:
        double index; // hurst index
        double coeff; // quality of the final fit 
        std::string func; // the final function for gnuplot 
        double store[2]; // to store the final func;
        double *diff; // 
        double *dfax; // 
        double *dfay; // 
        double *data_x; // 
        double *data_y; // 

    private:
        int i, j, k;
        int size_o; // size of fn ~ n 
        int size; // size of data 
        int order; // dfa order
        double rito;
        int minbox;
        int maxbox;
    private:
        int *seg_size;
        void mkbox();
        void solver();
};
class Powerspec{
    public:
        Powerspec(double *data, int size);
        ~Powerspec();
    public: 
        double maximum;
        int maxindex;
        double maxfreq;
        double maxperi;

        double rito;

        double *result;
        int size;
        int halfsize;
    private:
        void solver(double *data);
};

class Sync{
    public:
        Sync(double *in_datax, double *in_datay, int in_size);
        ~Sync();
    public:
        double *datax;
        double *dataxh;
        double *phasex;
        double *ampx;

        double *datay;
        double *datayh;
        double *phasey;
        double *ampy;

        int size;
        double sync_index;
    public:
        void plot(std::string name, int seg_size);

    private:
        int i;
        void solver();
};

void hilbert_trans(double *in,double *output,int num);

int fft(fftw_complex *in, fftw_complex *out, int num);

int ifft(fftw_complex *in, fftw_complex *out, int num);

class DFAI{
    public:
        DFAI(double *data, int data_size, int dfa_order);
        virtual ~DFAI();

        inline double get_rsq()
        {
            return coeff;
        }
        inline int get_order()
        {
            return order;
        }
        inline int get_size_o()
        {
            return size_o;
        }
        void init(double *data);        // init function;
        void recal(double *data); // reuse object without delete;
        double fitfunc(double input); // the final function; 
        void plot_over(int point); // plot out 
    public:
        double index; // hurst index
        double coeff; // quality of the final fit 
        std::string func; // the final function for gnuplot 
        double store[2]; // to store the final func;
        double *diff; // 
        double *dfax; // 
        double *dfay; // 
        double *data_x; // 
        double *data_xi;
        double *data_y; // 
        double *data_yi; 

    private:
        int i, j, k;
        int size_o; // size of fn ~ n 
        int size; // size of data 
        int order; // dfa order
        double rito;
        int minbox;
        int maxbox;
    private:
        int *seg_size;
        void mkbox();
        void solver();
};
