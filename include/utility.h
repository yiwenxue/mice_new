#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <ostream>
#include <gsl/gsl_statistics_double.h>

#include <mathematics.h>

#define BATCH3 

class mice{
    public:
        std::string mice_name;
        std::string filename;
        std::string fullpath;
        std::string path;
        std::string type;
    public:
        double *data_original; // all data read from the original file;
        double *data_original_t; // all data read from the original file;
        double *data;     // Only the usable data will be stored in this array;
        double *data_t;     // Only the usable data will be stored in this array;
        double *average;  // the average data 
        double *average_days;
        double *dfa;      // the dfa array containsthe index for every windows;
        double *dfa_days;
        double *dfa_days_err;
        double *std;      // the std data;
        double *std_days;

        int day_num;
        int ave_days_num;
        int seg_size; // size of single window segment
        int seg_num;  // the number of all windows
        int start;
        int length;

        double dfa_max;
        double dfa_min;

        double ave_rsq;
        double std_rsq;
        double dfa_rsq;
    public:
        mice(std::string filename, int windows_size,int start, int length);
        ~mice();
        void print_data();
        void print_details();
        inline int get_size(){
            return this->size_data;
        };
        inline int get_size_origin(){
            return this->size_original;
        };
        void plot_all();
        void plot_overview();
        void plot_ave_rhy();
        void plot_std_rhy();
        void plot_dfa_rhy();
        void plot_heatmap();
        void plot_powerspec();
        void detail_vision(int start, int length);
    protected:
        int d_seg_size;
        int d_seg_num;
    private:
        int size_original;    // original number of data, errors included 
        int size_data;   // the size of the data 
        std::ifstream infile;
    private:
        void mice_path_name();  // extract informations from the file name 
        void file_lines(void);  // find the totle lines in the data file 
        void read_mice_data(void); // read in the data with optimize (Nan, and lack of some points)
        void solve_windows(); // help to solve the windows numbers and windows size
        void solve_init(); // calculate the average value for each window.
};

int mktable(mice *mice1, mice *mice2, double sync_index, std::string batch, bool ifmutant);
