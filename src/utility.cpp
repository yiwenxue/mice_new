#include <cmath>
#include <cstdlib>
#include <gsl/gsl_statistics_double.h>
#include <math.h>
#include <string>
#include <utility.h>

mice::mice(std::string filename, int windows_size, int start_in, int length_in): fullpath(filename), seg_size(windows_size), start(start_in), length(length_in){
    mice_path_name();
    read_mice_data();
    solve_windows();
    solve_init();
}

mice::~mice(){
    /* delete data_ori; delete data; delete average; delete dfa; delete std; */
    delete[] data; delete[] data_t;
    delete[] data_original_t; delete[] data_original;
    delete[] average; delete[] std;
    delete[] average_days; delete[] std_days;
    delete[] dfa; delete[] dfa_days;
    delete[] dfa_days_err;
    std::cout << "Mice " + mice_name + "_" + type + " deleted" << std::endl;
}

// tempt function for debug
void mice::print_data(){
    for (int i = 0; i < size_original; ++i) {
        std::cout << "Line: " << i << "\tval:" << data[i] << std::endl;
    }
}

// minimum function to get filename, micename, file path, from full path.
void mice::mice_path_name(){
#ifdef __WIN32
    char sep = '\\';
#else
    char sep = '/';
#endif
    auto mark = fullpath.rfind(sep);
    path = fullpath.substr(0,mark+1);
    filename = fullpath.substr(mark+1,fullpath.size()-mark-1); // minus 1 because I don't want the '.' char in my name(../data/filename.activity.txt).
    sep = '.';
    mark = filename.find(sep);
    mice_name = filename.substr(0,mark);
    type = filename.substr(mark+1,filename.rfind(sep)-mark-1);
}

// easy function to output the general details of the mice.
void mice::print_details(){
    std::cout << "file:          " << filename << std::endl;
    std::cout << "path:          " << path << std::endl;
    std::cout << "name:          " << mice_name << std::endl;
    std::cout << "type:          " << type << std::endl;
    std::cout << "totle lines:   " << size_original << std::endl;
    std::cout << "reduced lines: " << size_data << std::endl;
    std::cout << "seg_num:       " << seg_num << std::endl;
    std::cout << "seg_size:      " << seg_size << std::endl;
}

void mice::file_lines(void){
    FILE *input = fopen(fullpath.c_str(),"r");
    if(input == NULL){
        fprintf(stderr,"(File Lines)Failed to read %s\n",filename.c_str());
        return;
    }

    //Count the totla lines;
    char chr;
    size_original = 0;
    chr = fgetc(input);
    while(chr != EOF){
        if(chr == '\n')
            size_original ++;
        chr = fgetc(input);
    }
    fclose(input);
}

void mice::read_mice_data(void){
    file_lines();
    FILE* temp = fopen(fullpath.c_str(),"r");
    //Check if one can read this file;
    if(temp == NULL){
        fprintf(stderr,"(Read ALL)Failed to read %s\n",filename.c_str());
        return;
    }

    data_original = new double[size_original];
    data_original_t = new double[size_original];
    for (int i = 0; i < size_original; ++i) {
        data_original[i] = 0.0;
        data_original_t[i] = 0.0;
    }

    double tempx = 0, tempy = 0;
    /* // read in all data into data_original */
    for (int i = 1; i < size_original; ++i) {
        if(fscanf (temp,"%lf\t%lf\n",&tempx, &tempy) == 2){
            /* // remove all nan numbers */ 
            if(!std::isnan(tempy))
                data_original[i] = tempy;
            else 
                data_original[i] = data_original[i-1];
            data_original_t[i] = tempx;
        }
    }
    fclose(temp);

    size_data = data_original_t[size_original-1]/10;
    double *temp_data = new double[size_data];
    double *temp_data_t = new double[size_data];
    int datamrk = 0, diff_t = 0;
    for (int i=1;i<size_original;i++){
        diff_t = (data_original_t[i] - data_original_t[i-1])/10;
        if( diff_t == 1 ){
            temp_data_t[datamrk] = data_original_t[i];
            temp_data[datamrk] = data_original[i];
            datamrk ++;
        }else if (diff_t < 100){
            for (int j = 1; j <= diff_t; ++j) {
                temp_data[datamrk] = data_original[i-1] + 
                    (data_original[i]-data_original[i-1])/diff_t * j;
                temp_data_t[datamrk] = data_original_t[i-1] + 
                    (data_original_t[i]-data_original_t[i-1])/diff_t * j;
                datamrk ++;
            }
        }
        if(diff_t > 100.0){
#ifdef BATCH3
            {datamrk = 0; std::cout << "for Batch3 or 5" << std::endl;}
#else 
            {break; std::cout << "for Batch 6" << std::endl; } 
#endif
        }
    }
    size_data = (datamrk/(48*360)) * 48*360; // datamrk is the reduced size, and then the size_data is integer multiple 2 days.
    int tmp =( start+length ) * 360 * 24;
    int ttmp = 0;
    if(tmp > size_data){
        std::cerr << "Error, data not enough" << std::endl;
        data = new double[size_data];
        data_t = new double[size_data];
        for (int i = 0; i < size_data; ++i) {
            data[i] = temp_data[i];
            data_t[i] = temp_data_t[i];
        }
    }else {size_data = tmp;
        ttmp = start * 360 * 24;
        data = new double[size_data];
        data_t = new double[size_data];
        for (int i = ttmp; i < tmp; ++i) {
            data[i] = temp_data[i];
            data_t[i] = temp_data_t[i];
        }
    }
    delete[] temp_data_t; delete[] temp_data;
}

void mice::solve_windows(){
    seg_num = size_data / seg_size;
    d_seg_size = 1080;
    d_seg_num = size_data / d_seg_size;
}

void mice::solve_init(){
    int count = 48 * 360 / seg_size; 
    ave_days_num = count;

    average = new double[seg_num];
    std = new double[seg_num];
    average_days = new double[count];
    std_days = new double[count];

    for (int i = 0; i < seg_num; ++i) {
        average[i] = gsl_stats_mean(data+(i*seg_size), 1, seg_size);
        std[i] = gsl_stats_sd(data+(i*seg_size), 1, seg_size);
    }
    for (int i = 0; i < count; ++i) {
        average_days[i] = gsl_stats_mean(average+i, count, seg_num/count);
        std_days[i] = gsl_stats_mean(std+i, count, seg_num/count);
    }

    dfa = new double[d_seg_num];
    dfa_days = new double[48 * 360 / d_seg_size];
    dfa_days_err = new double[48 * 360 / d_seg_size];
}

void mice::plot_all(){
    std::cout << "Drawing over view..." << std::flush;
    this->plot_overview();
    std::cout << " Finished!" << std::endl;
    std::cout << "Drawing heatmap..." << std::flush;
    this->plot_heatmap();
    std::cout << " Finished!" << std::endl;
    std::cout << "Drawing std rhythm..." << std::flush;
    this->plot_std_rhy();
    std::cout << " Finished!" << std::endl;
    std::cout << "Drawing average rhythm..." << std::flush;
    this->plot_ave_rhy();
    std::cout << " Finished!" << std::endl;
    std::cout << "Drawing dfa rhythm..." << std::flush;
    this->plot_dfa_rhy();
    std::cout << " Finished!" << std::endl;
    std::cout << "Drawing powerspectrum..." << std::flush;
    this->plot_powerspec();
    std::cout << " Finished!" << std::endl;
}

void mice::plot_overview(){
    GnuplotPipe gp;
    gp.sendLine("set terminal pdf size 5,4");
    gp.sendLine("set title \'mice:" + mice_name + " " + type + " overview(" + std::to_string(seg_size) + ")\'" );
    gp.sendLine("set xlabel \'days\'");
    gp.sendLine("set xtics 2");
    gp.sendLine("set grid");
    gp.sendLine("set ylabel '"+ type + "'");
    gp.sendLine("set output '" + mice_name + "_" + type + "_overview.pdf'");
    gp.sendLine("plot '-' w lines lw 1 title 'Average:" + 
            std::to_string(gsl_stats_mean(data, 1, size_data)).substr(0,4) + "'");
    for (int i = 0; i < seg_num; ++i) {
        gp.sendLine( std::to_string(i*(seg_size/8640.0)) + "\t" + std::to_string(average[i]));
    }
    gp.sendLine("e");
    gp.sendLine("set out");
};

void mice::plot_ave_rhy(){
    GnuplotPipe gp;
    gp.sendLine("set terminal pdf size 5,4");
    gp.sendLine("set xrange [0:48]");
    gp.sendLine("set grid");
    // send data to gnuplot
    int count = 48*360/seg_size; 
    {   
        gp.sendLine("set ylabel '"+ type + "'");
        gp.sendLine("set title 'mice:" + mice_name + " " + type + " rhythm(" + std::to_string(seg_size) + ")'" );
        gp.sendLine("set xlabel 'hours'");
        gp.sendLine("set output '" + mice_name + "_" + type + "_rhythm.pdf'");
        double *data_x = new double[count];
        for (int i = 0; i < count; ++i) {
            data_x[i] = i * seg_size;
        }
        Cosinor new_cos(data_x, average_days, count, 5);
        ave_rsq = new_cos.get_rsq();
        gp.sendLine("plot '-' w points ps 1 notitle ," 
                + new_cos.func + "title 'R^2=" + std::to_string(ave_rsq) + "'" );
        /* std::cout << new_cos.func << std::endl; */
        for (int i = 0; i < count; ++i) {
            gp.sendLine(std::to_string((i)*seg_size/360.0) + "\t" + std::to_string(average_days[i]));
        }
        gp.sendLine("e");
        gp.sendLine("set out");
        delete[] data_x;
    }
}

void mice::plot_std_rhy(){
    if(type == "Activity"){
        GnuplotPipe gp;
        gp.sendLine("set terminal pdf size 5,4");
        gp.sendLine("set xrange [0:48]");
        gp.sendLine("set ylabel '"+ type + "'");
        gp.sendLine("set xlabel 'hours'");
        gp.sendLine("set grid");
        // send data to gnuplot
        int count = 48*360/seg_size; 

        {   double *data_x = new double[count];
            for (int i = 0; i < count; ++i) {
                data_x[i] = i * seg_size;
            }
            Cosinor new_cos(data_x, std_days, count, 5);
            std_rsq = new_cos.get_rsq();
            gp.sendLine("set title 'mice:" + mice_name + " " + type + " Std rhythm(" + std::to_string(seg_size) + ")'" );
            gp.sendLine("set output '" + mice_name + "_" + type + "_std_rhythm.pdf'");
            gp.sendLine("plot '-' w points ps 1 notitl , "
                    + new_cos.func + "title 'R^2=" + std::to_string(std_rsq) + "'" );
            /* std::cout << new_cos.func << std::endl; */

            for (int i = 0; i < count; ++i) {
                gp.sendLine(std::to_string((i)*seg_size/360.0) + "\t" + std::to_string(std_days[i]));
            }
            gp.sendLine("e");
            gp.sendLine("set out");
            delete[] data_x;
        }
    }
}

void mice::plot_dfa_rhy(){
    int count = 48 * 360 / d_seg_size; 
    for (int i = 0; i < d_seg_num; ++i) {      
        DFA new_dfa(data+(i*d_seg_size), d_seg_size, 1);
        dfa[i] = new_dfa.index;
        if(isnan(dfa[i]) || std::isinf(dfa[i])){
            new_dfa.plot_over(i);
            std::cout << "DFA[" << i << "] = " << dfa[i] << " dfax[last]: " << new_dfa.dfax[new_dfa.get_size_o()-1] << std::endl;
        }
    }
    for (int i = 0; i < count; ++i) {
        dfa_days[i] = gsl_stats_mean(dfa+i, count, d_seg_num/count);
    }

    GnuplotPipe gp;
    gp.sendLine("set terminal pdf size 5,4");
    gp.sendLine("set xrange [0:48]");
    gp.sendLine("set grid");
    // send data to gnuplot

    {   double *data_x = new double[count];
        for (int i = 0; i < count; ++i) {
            data_x[i] = i * d_seg_size;
        }
        Cosinor new_cos(data_x, dfa_days, count, 5);
        dfa_rsq = new_cos.get_rsq();
        dfa_max = gsl_stats_max(dfa_days, 1, count);
        dfa_min = gsl_stats_min(dfa_days, 1, count);
        gp.sendLine("set title \'mice:" + mice_name + " " + type + " dfa rhythm(" + std::to_string(d_seg_size) + ")\'" );
        gp.sendLine("set output '" + mice_name + "_" + type + "_dfa_rhythm.pdf'");
        gp.sendLine("set xlabel 'hours'");
        gp.sendLine("set ylabel 'α'");
        gp.sendLine("set grid");
        gp.sendLine("plot '-' w points ps 1 title 'Range:" 
                + std::to_string(dfa_min).substr(0,4) + "-" 
                + std::to_string(dfa_max).substr(0,4) + "',"
                + new_cos.func + "title 'R^2=" + std::to_string(dfa_rsq) + "'" );
        for (int i = 0; i < count; ++i) {
            gp.sendLine(std::to_string((i)*d_seg_size/360.0) + "\t" + std::to_string(dfa_days[i]));
        }
        gp.sendLine("e");
        gp.sendLine("set out");
        delete[] data_x;
    }
}

void mice::plot_powerspec(){
    Powerspec new_power(data,size_data);
    GnuplotPipe gp;
    gp.sendLine("set terminal pdf size 5,4");
    gp.sendLine("set xrange [" + std::to_string(-1/1000.0/20.0) + ":" + std::to_string(1/1000.0) + "]");
    gp.sendLine("set ylabel '"+ type + "'");
    gp.sendLine("set xlabel 'frequency'");
    gp.sendLine("set grid");
    // send data to gnuplot
    {
        gp.sendLine("set title 'mice:" + mice_name + " " + type + " powerspectrum'" );
        gp.sendLine("set output '" + mice_name + "_" + type + "_powerspec.pdf'");
        gp.sendLine("plot '-' w lines lw 1 title 'Maxindex:" + std::to_string(new_power.maxperi).substr(0,4) + "(h)'  ");
        for (int i = 0; i < new_power.halfsize ; ++i) {
            gp.sendLine(std::to_string(i*new_power.rito) + "\t" + std::to_string(new_power.result[i]));
        }
        gp.sendLine("e");
        gp.sendLine("set out");
    }
}

void mice::plot_heatmap(){
    GnuplotPipe gp;
    gp.sendLine("set terminal pdf size 5,4");
    gp.sendLine("set title \'mice:" + mice_name + " " + type + " heatmap(" + std::to_string(seg_size) + ")\'" );
    gp.sendLine("set output '" + mice_name + "_" + type + "_heatmap.pdf'");
    gp.sendLine("set tics nomirror out scale 0.75");
    gp.sendLine("set xlabel 'hours'");
    gp.sendLine("set ylabel 'days'");

    gp.sendLine("$map3 << EOF");
    // send data to gnuplot
    int raws = size_data / 360 / 48; 
    day_num = size_data / 360 / 24;
    int count = 48 * 360 / seg_size;
    std::string temp = "" ;
    int xtics = count / 12;
    for (int j = 0; j < count; ++j) {
        if((j+xtics/2)%xtics == 0)
            temp += ", " + std::to_string(j*48.0/count).substr(0,4);
        else 
            temp += ", ";
    }
    gp.sendLine(temp);
    for (int i = 0; i < raws; i++) {
        temp = std::to_string(i*2+2) ;
        for (int j = 0; j < count; ++j) {
               temp += ", " + std::to_string(average[i*count+j]);
        }
        gp.sendLine(temp);
    }
    gp.sendLine("EOF");
    gp.sendLine("set datafile separator comma");
    gp.sendLine("plot '$map3' matrix rowheaders columnheaders using 1:2:3 with image notitle ");
    gp.sendLine("set datafile separator");
    gp.sendLine("set out");
}

void mice::detail_vision(int start, int length){
    GnuplotPipe gp;
    gp.sendLine("set terminal pdf size 5,4");
    gp.sendLine("set title \'mice:" + mice_name + " " + type + "detail");
    gp.sendLine("set xlabel 'sec'");
    gp.sendLine("set grid");
    gp.sendLine("set output '" + mice_name + "_" + type + "_det.pdf'");
    gp.sendLine("plot '-' w lines lw 1 title 'Average:" + 
            std::to_string(gsl_stats_mean(data, 1, size_data)).substr(0,4) + "'");
    length = start + length;
    for (int i = start; i < length; ++i) {
        gp.sendLine( std::to_string(data_t[i]) + "\t" + std::to_string(data[i]));
    }
    gp.sendLine("e");
    gp.sendLine("set out");
};

int mktable(mice *mice1,mice *mice2, double sync_index, std::string batch, bool ifmutant){
    std::string output;
    if(mice1->type == "Activity"){
        output = 
            "\"" + batch + "\"," + // mice batch
            "\"" + mice1->mice_name + "\"," + // mice name
            (ifmutant?("\"mutant\""):("\"WT\"")) + "," + // mutant?
            std::to_string(mice1->day_num) + "," + // Num od days

            std::to_string(gsl_stats_mean(mice1->average, 1, mice1->seg_num) ) + "," + // Activity average rsq
            std::to_string(mice1->ave_rsq) + "," + // Activity average rsq
            std::to_string(mice1->std_rsq) + "," + // Activity Std rsq
            "\"" + std::to_string(mice1->dfa_min).substr(0,4) + "~" + std::to_string(mice1->dfa_max).substr(0,4) + "\"," + // DFA range
            std::to_string(mice1->dfa_rsq) + "," + // DFA rhythm rsq 
            // Then Temperature 
            std::to_string(gsl_stats_mean(mice2->average, 1, mice2->seg_num) ) + "," + // Activity average rsq
            std::to_string(mice2->ave_rsq) + "," + // Activity average rsq
            "\"" + std::to_string(mice2->dfa_min).substr(0,4) + "~" + std::to_string(mice2->dfa_max).substr(0,4) + "\"," + // DFA range
            std::to_string(mice2->dfa_rsq) + "," + // DFA rhythm rsq 
            std::to_string(sync_index);
    } else if(mice1->type == "Temperature"){
        output = 
            "\"" + batch + "\"," + // mice batch
            "\"" + mice1->mice_name + "\"," + // mice name
            (ifmutant?("\"mutant\""):("\"WT\"")) + "," + // mutant?
            std::to_string(mice2->day_num) + "," + // Num od days

            std::to_string(gsl_stats_mean(mice2->average, 1, mice2->seg_num) ) + "," + // Activity average rsq
            std::to_string(mice2->ave_rsq) + "," + // Activity average rsq
            std::to_string(mice2->std_rsq) + "," + // Activity Std rsq
            "\"" + std::to_string(mice2->dfa_min).substr(0,4) + "~" + std::to_string(mice2->dfa_max).substr(0,4) + "\"," + // DFA range
            std::to_string(mice2->dfa_rsq) + "," + // DFA rhythm rsq 
            // Then Temperature 
            std::to_string(gsl_stats_mean(mice1->average, 1, mice1->seg_num) ) + "," + // Activity average rsq
            std::to_string(mice1->ave_rsq) + "," + // Activity average rsq
            "\"" + std::to_string(mice1->dfa_min).substr(0,4) + "~" + std::to_string(mice1->dfa_max).substr(0,4) + "\"," + // DFA range
            std::to_string(mice1->dfa_rsq) + "," + // DFA rhythm rsq 
            std::to_string(sync_index);
    }

    std::cout << output << std::endl;
    system(("echo -e "+output+" >> table.txt").c_str());

    return 0;
}
