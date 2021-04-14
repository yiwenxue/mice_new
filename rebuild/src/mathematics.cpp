#include <Rtypes.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <gsl/gsl_statistics_double.h>
#include <iterator>
#include <math.h>
#include <string>
#ifndef MATHEMATICS
#include <mathematics.h>
#endif

DFA_MAG::DFA_MAG(double *data, int data_size, int dfa_order) : size(data_size), order(dfa_order)
{
    rito = pow(2.0, 1.0 / 5);
    minbox = 2 * (order + 1);
    maxbox = size / 4;
    size_o = log10(maxbox / (double)minbox) / log10(rito) + 1.5;

    // allocated memory
    data_x = new double[size];
    data_y = new double[size];
    diff = new double[size];
    dfax = new double[size_o];
    dfay = new double[size_o];
    seg_size = new int[size_o];

    mkbox();
    init(data);
    solver();

    //  For Debug
    /* TCanvas c1(""); */
    /* c1.cd(); */
    /* TGraph gr_original(size, data_x, data); */
    /* gr_original.Draw(); */
    /* c1.Update(); */
    /* c1.Print("dfa_test_originaldata.pdf"); */
    /* c1.Clear(); */

    /* TGraph gr_profile(size, data_x, data_y); */
    /* gr_profile.Draw(); */
    /* c1.Update(); */
    /* c1.Print("dfa_test_profile.pdf"); */
    /* c1.Clear(); */
};

DFA_MAG::~DFA_MAG()
{
    delete[] dfax;
    delete[] dfay;
    delete[] data_x;
    delete[] data_y;
    delete[] diff;
    delete[] seg_size;
};

void DFA_MAG::init(double *data)
{
    // First step: To get the profile
    // remove and integral
    double rm = gsl_stats_mean(data, 1, size);
    diff[0] = 0.0;
    data_y[0] = data[0] - rm;
    data_x[0] = 0;
    for (i = 1; i < size; ++i)
    {
        diff[i] = 0.0;
        data_y[i] = fabs(data[i] - rm) + data_y[i - 1];
        data_x[i] = i;
    }
    for (i = 0; i < size_o; ++i)
    {
        dfax[i] = 0.0;
        dfay[i] = 0.0;
    }
}

void DFA_MAG::mkbox()
{
    for (i = 1, j = 1, seg_size[0] = minbox; j < size_o && seg_size[j - 1] < maxbox; ++i)
        if ((seg_size[j] = minbox * pow(rito, i) + 0.5) > seg_size[j - 1])
            j++;
    size_o = j;
}

void DFA_MAG::recal(double *data)
{
    init(data);
    solver();
}

void DFA_MAG::solver()
{
    // init the window info
    int end = 0, tail = 0;
    Polyfit new_fit(data_x, data_y, 10, order + 1);

    for (i = 0; i < size_o; i++)
    {
        end = size - seg_size[i];
        new_fit.resize(seg_size[i]);
        for (j = 0; j < end; j += seg_size[i])
        {
            new_fit.recal(data_x + j, data_y + j);
            tail = j + seg_size[i];
            for (k = j; k < tail; k++)
            {
                diff[k] = new_fit.fitfunc(data_x[k]) - data_y[k];
                diff[k] = diff[k] * diff[k];
            }
        }
        dfay[i] = gsl_stats_mean(diff, 1, end);
        dfay[i] = sqrt(dfay[i]);
        dfax[i] = seg_size[i];
        // log lize
        dfay[i] = log10(dfay[i]);
        dfax[i] = log10(dfax[i]);
    }

    {
        // fit loged function with least linear regression
        Polyfit new_fit(dfax, dfay, size_o, 2);
        coeff = new_fit.get_rsq();
        store[0] = new_fit.store[0];
        store[1] = new_fit.store[1];
        index = store[1];
    }

    for (int i = 0; i < size_o; ++i)
    {
        dfay[i] = pow(10.0, dfay[i]);
        dfax[i] = pow(10.0, dfax[i]);
    }

    /* func = std::to_string(store[0]) + "+x*" + std::to_string(store[1]); */
    func = "pow(10.0," + std::to_string(store[0]) + ") *pow(x," + std::to_string(store[1]) + ")";
}

void DFA_MAG::plot_over(int point)
{
    GnuplotPipe gp;
    gp.sendLine("set terminal pdf size 5,4");
    gp.sendLine("set output 'DFA_MAG_" + std::to_string(point) + ".pdf'");
    gp.sendLine("set title 'DFA_MAG'");
    gp.sendLine("set logscale x 10");
    gp.sendLine("plot '-' w points ps 1.5 notitle");
    for (i = 0; i < size_o; ++i)
    {
        gp.sendLine(std::to_string(dfax[i]) + "\t" +
                    std::to_string(dfay[i]));
    }
    gp.sendLine("e");
    gp.sendLine("set output");
}

Polyfit::Polyfit(double *datax, double *datay, int _size, int _degree)
    : size(_size), degree(_degree)
{
    store = new double[degree];
    X = gsl_matrix_alloc(size, degree);
    y = gsl_vector_alloc(size);
    c = gsl_vector_alloc(degree);
    cov = gsl_matrix_alloc(degree, degree);
    ws = gsl_multifit_linear_alloc(size, degree);

    init(datax, datay);
    polyfit();
};

Polyfit::~Polyfit()
{
    delete[] store;
    gsl_multifit_linear_free(ws);
    gsl_matrix_free(X);
    gsl_matrix_free(cov);
    gsl_vector_free(y);
    gsl_vector_free(c);
}

void Polyfit::resize(int new_size)
{
    gsl_multifit_linear_free(ws);
    gsl_matrix_free(X);
    gsl_vector_free(y);

    size = new_size;
    X = gsl_matrix_alloc(size, degree);
    y = gsl_vector_alloc(size);
    ws = gsl_multifit_linear_alloc(size, degree);
}

double Polyfit::fitfunc(double input)
{
    double temp = store[0];
    for (i = 1; i < degree; ++i) {
        temp += store[i] * pow(input,i);
    }
    return temp;
};

void Polyfit::init(double *datax,double *datay)
{
    meanx = gsl_stats_mean(datax,1,size); 

    for(i=0; i < size; i++) {
        gsl_matrix_set(X, i, 0, 1.0);
        for(j=1; j < degree; j++) {
            gsl_matrix_set(X, i, j, pow(datax[i]-meanx, j));
        }
        gsl_vector_set(y, i, datay[i]);
    }
    tss = gsl_stats_tss(datay,1,size);
}

void Polyfit::recal(double *datax, double *datay)
{
    init(datax, datay);
    polyfit();
}

bool Polyfit::polyfit()
{
    // calculate
    gsl_multifit_linear(X, y, c, cov, &coeff, ws);
    coeff = 1.0 - coeff/tss;

    for(i=0; i < degree; i++){
        store[i] = gsl_vector_get(c, i);
    }
    for(i=0;i<degree;i++){
        for(int j=i+1;j<degree;j++){
            int order = gsl_sf_fact(j)/gsl_sf_fact(i)/gsl_sf_fact(j-i);
            store[i]+= store[j]*order*pow(-1*meanx,j-i);
        }
    }
    // this will give you a string to describe the regression function
    func = std::to_string(store[0]);
    for (i = 1; i < degree; ++i) {
        func += " + " + std::to_string(store[i]) + "*x**" + std::to_string(i);
    }
    return true;
};

Cosinor::Cosinor(double *data_x, double *data_y, int _size, int _degree): size(_size), degree(_degree)
{
    omega = M_PI/360.0/12.0;
    store = new double[degree];
    X = gsl_matrix_alloc (size,degree);
    Y = gsl_vector_alloc (size);
    c = gsl_vector_alloc (degree);
    cov = gsl_matrix_alloc (degree,degree);
    ws = gsl_multifit_linear_alloc (size,degree);

    init(data_x, data_y);
    cosinor();
};

Cosinor::~Cosinor()
{
    delete[] store;
    gsl_multifit_linear_free (ws);
    gsl_matrix_free (X);
    gsl_vector_free (Y);
    gsl_vector_free (c);
    gsl_matrix_free (cov);
}

double Cosinor::fitfunc(double input)
{
    double temp = store[0];
    for(j=1;j<degree;j+=2){
        temp += store[j] * cos( omega * input) +
            store[j+1] * sin( omega * input); 
    }
    return temp;
};

void Cosinor::init(double *datax, double *datay)
{
    for(i=0;i<size;i++){
        gsl_matrix_set (X, i, 0, 1.0);
        gsl_vector_set (Y, i, datay[i]);
    }
    for(j=1;j<degree;j+=2){
        for(i=0;i<size;i++){
            gsl_matrix_set (X, i, j, cos(omega * (j+1)/2 * datax[i]));
            gsl_matrix_set (X, i, j+1, sin(omega * (j+1)/2 * datax[i]));
        }
    }
    tss = gsl_stats_tss(datay,1,size);
}

void Cosinor::recal(double *datax, double *datay)
{
    init(datax, datay);
    cosinor();
}

bool Cosinor::cosinor()
{
    gsl_multifit_linear (X, Y, c, cov, &coeff, ws);
    coeff = 1.0 - coeff/tss;
    for(i=0;i<degree;i++){
        store[i] = gsl_vector_get (c,i);
    }
    // this will give you a string to describe the regression function
    func = std::to_string(store[0]);
    for (i = 1; i < degree; i+=2 ) {
        func += " + " + std::to_string(store[i]) + "*cos(" + std::to_string(omega *(i+1)/2 *360) + "*x)" 
            + " + " + std::to_string(store[i+1]) + "*sin(" + std::to_string(omega *(i+1)/2 *360) + "*x)";
    }
    return true;
};

DFA::DFA(double *data, int data_size, int dfa_order): size(data_size), order(dfa_order)
{
    rito = pow(2.0, 1.0 / 5);
    minbox = 2*(order+1);
    maxbox = size/4;
    size_o = log10(maxbox / (double)minbox)/log10(rito) + 1.5;

    // allocated memory
    data_x = new double[size];
    data_y = new double[size];
    diff = new double[size];
    dfax = new double[size_o];
    dfay = new double[size_o];
    seg_size = new int[size_o];

    mkbox();
    init(data);
    solver();

    //  For Debug
    /* TCanvas c1(""); */
    /* c1.cd(); */
    /* TGraph gr_original(size, data_x, data); */
    /* gr_original.Draw(); */
    /* c1.Update(); */
    /* c1.Print("dfa_test_originaldata.pdf"); */
    /* c1.Clear(); */

    /* TGraph gr_profile(size, data_x, data_y); */
    /* gr_profile.Draw(); */
    /* c1.Update(); */
    /* c1.Print("dfa_test_profile.pdf"); */
    /* c1.Clear(); */
};

DFA::~DFA()
{
    delete[] dfax;
    delete[] dfay;
    delete[] data_x;
    delete[] data_y;
    delete[] diff;
    delete[] seg_size;
};

void DFA::init(double *data)
{
    // First step: To get the profile
    // remove and integral 
    double rm = gsl_stats_mean(data, 1, size);
    diff[0] = 0.0;
    data_y[0] = data[0] - rm;
    data_x[0] = 0;
    for (i = 1; i < size; ++i) {
        diff[i] = 0.0;
        data_y[i] = data[i] - rm + data_y[i-1];
        data_x[i] = i;
    }
    for (i = 0; i < size_o; ++i) {
        dfax[i] = 0.0;
        dfay[i] = 0.0;
    }
}

void DFA::mkbox()
{
    for(i = 1, j = 1, seg_size[0] = minbox; j < size_o && seg_size[j-1] < maxbox; ++i)
        if((seg_size[j] = minbox * pow(rito,i) + 0.5) > seg_size[j-1])
            j++;
    size_o = j;
}

void DFA::recal(double *data)
{
    init(data);
    solver();
}

void DFA::solver()
{
    // init the window info
    int end = 0, tail = 0; 
    Polyfit new_fit( data_x, data_y, 10, order+1);

    for (i = 0; i < size_o; i++) { 
        end = size - seg_size[i];
        new_fit.resize(seg_size[i]);
        for (j=0; j < end; j+= seg_size[i]){
            new_fit.recal(data_x+j, data_y+j);
            tail = j+seg_size[i];
            for(k=j; k< tail; k++){
                diff[k] = new_fit.fitfunc(data_x[k]) - data_y[k];
                diff[k] = diff[k]*diff[k];
            }
        }
        dfay[i] = gsl_stats_mean(diff, 1, end);
        dfay[i] = sqrt(dfay[i]);
        dfax[i] = seg_size[i];
        // log lize
        dfay[i] = log10(dfay[i]);
        dfax[i] = log10(dfax[i]);
    }

    {
        // fit loged function with least linear regression
        Polyfit new_fit(dfax,dfay,size_o,2);
        coeff = new_fit.get_rsq();
        store[0] = new_fit.store[0]; 
        store[1] = new_fit.store[1];
        index = store[1];
    }

    for (int i = 0; i < size_o; ++i) {
        dfay[i] = pow(10.0,dfay[i]);
        dfax[i] = pow(10.0,dfax[i]);
    }

    /* func = std::to_string(store[0]) + "+x*" + std::to_string(store[1]); */
    func = "pow(10.0,"+std::to_string(store[0]) + ") *pow(x,"+std::to_string(store[1]) +")"  ;
}

void DFA::plot_over(int point)
{
    GnuplotPipe gp; 
    gp.sendLine("set terminal pdf size 5,4");
    gp.sendLine("set output 'DFA_" + std::to_string(point) + ".pdf'");
    gp.sendLine("set title 'DFA'");
    gp.sendLine("set logscale x 10");
    gp.sendLine("plot '-' w points ps 1.5 notitle" );
    for (i = 0; i < size_o; ++i) {
        gp.sendLine(std::to_string(dfax[i]) + "\t" + 
                std::to_string(dfay[i]) );
    }
    gp.sendLine("e");
    gp.sendLine("set output");

}

double DFA::fitfunc(double input) 
{
    return store[0] + store[1]*input;
}

Powerspec::Powerspec(double *data, int new_size)
{
    size = new_size;
    halfsize = size / 2;
    result = new double[halfsize];
    
    solver(data); 
}

Powerspec::~Powerspec()
{
    delete[] result;
}

void Powerspec::solver(double *data)
{
    double mean = gsl_stats_mean(data, 1, size);
    // process fft 
    fftw_complex * input, * output;
    input = fftw_alloc_complex(size);
    output = fftw_alloc_complex(size);

    for (int i = 0; i < size; ++i) {
       input[i][RE] = data[i] - mean; input[i][IM] = 0.0; 
       output[i][RE] = 0.0; output[i][IM] = 0.0; 
    }
    fftw_plan p;
    p=fftw_plan_dft_1d(size,input,output,-1,FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);


    // set the information 
    rito = 1 / 10.0 / size; // the interval between ~; 

    // save the amplitude info into amp array;
    double total = 0; 
    for (int i = 0; i < halfsize; ++i) {
        result[i] = output[i][RE]*output[i][RE] + output[i][IM]*output[i][IM];
        total += result[i];
    }
    double max = gsl_stats_max(result, 1, halfsize);
    for (int i = 0; i < halfsize; ++i) {
        result[i] /= total;
    }

    // find the largest peek in the array;
    max = gsl_stats_max(result, 1, halfsize);
    maxindex = gsl_stats_max_index(result, 1, halfsize);
    maxfreq = maxindex*rito; 
    maxperi = (1.0/maxfreq)/3600.0;
    
    fftw_free(input);
    fftw_free(output);
}

Sync::Sync(double *in_datax, double *in_datay, int in_size)
{
    size = in_size;
    datax = new double[size];
    dataxh = new double[size];
    ampx = new double[size];
    phasex = new double[size];
    datay = new double[size];
    datayh = new double[size];
    ampy = new double[size];
    phasey = new double[size];

    phase_diffx = new double[size - 1];
    phase_diffy = new double[size - 1];
    pd = new double[size];
    pdcirx = new double[size];
    pdciry = new double[size];

    double meanx = gsl_stats_mean(in_datax, 1, size);
    double meany = gsl_stats_mean(in_datay, 1, size);
    for(i = 0;i < size; ++i){
        datax[i] = in_datax[i] - meanx;
        datay[i] = in_datay[i] - meany;

        dataxh[i] = 0.0;
        ampx[i] = 0.0;
        phasex[i] = 0.0;

        datayh[i] = 0.0;
        ampy[i] = 0.0;
        phasey[i] = 0.0;
    }
    SetIndex(1, 1);
    SetRange(0, size);
    solver();
}
Sync::Sync(double *in_datax, double *in_datay, int in_size, int _start, int _end, int type)
{
    size = in_size;
    datax = new double[size];
    dataxh = new double[size];
    ampx = new double[size];
    phasex = new double[size];
    datay = new double[size];
    datayh = new double[size];
    ampy = new double[size];
    phasey = new double[size];

    phase_diffx = new double[size - 1];
    phase_diffy = new double[size - 1];
    pd = new double[size];
    pdcirx = new double[size];
    pdciry = new double[size];

    double meanx = gsl_stats_mean(in_datax, 1, size);
    double meany = gsl_stats_mean(in_datay, 1, size);
    for (i = 0; i < size; ++i)
    {
        datax[i] = in_datax[i] - meanx;
        datay[i] = in_datay[i] - meany;

        dataxh[i] = 0.0;
        ampx[i] = 0.0;
        phasex[i] = 0.0;

        datayh[i] = 0.0;
        ampy[i] = 0.0;
        phasey[i] = 0.0;
    }
    SetIndex(1, 1);
    SetRange(_start, _end);
    solver();
}
Sync::Sync(double *in_datax, double *in_datay, int in_size, int n, int m)
{
    size = in_size;
    datax = new double[size];
    dataxh = new double[size];
    ampx = new double[size];
    phasex = new double[size];
    datay = new double[size];
    datayh = new double[size];
    ampy = new double[size];
    phasey = new double[size];

    phase_diffx = new double[size - 1];
    phase_diffy = new double[size - 1];
    pd = new double[size];
    pdcirx = new double[size];
    pdciry = new double[size];

    double meanx = gsl_stats_mean(in_datax, 1, size);
    double meany = gsl_stats_mean(in_datay, 1, size);
    for (i = 0; i < size; ++i)
    {
        datax[i] = in_datax[i] - meanx;
        datay[i] = in_datay[i] - meany;

        dataxh[i] = 0.0;
        ampx[i] = 0.0;
        phasex[i] = 0.0;

        datayh[i] = 0.0;
        ampy[i] = 0.0;
        phasey[i] = 0.0;
    }
    SetIndex(n, m);
    SetRange(0, size);
    solver();
}

Sync::~Sync()
{
    delete[] datax;
    delete[] dataxh;
    delete[] ampx;
    delete[] phasex;
    delete[] datay;
    delete[] datayh;
    delete[] ampy;
    delete[] phasey;

    delete[] phase_diffx;
    delete[] phase_diffy;
    delete[] pd;
    delete[] pdcirx;
    delete[] pdciry;
}

void Sync::GetPhase()
{
    hilbert_trans(datax, dataxh, size);
    hilbert_trans(datay, datayh, size);

    size = end;

    /* double minimumx = gsl_stats_mean(datax, 1, size); */
    /* double minimumy = gsl_stats_mean(datay, 1, size); */
    /* double meanxh =  gsl_stats_mean(dataxh, 1, size); */
    /* double meanyh =  gsl_stats_mean(datayh, 1, size); */

    /* for (i = 0; i < size; ++i) { */
    /*     /1* datax[i] += -3*(minimumx); *1/ */
    /*     /1* datay[i] += -2*(minimumy); *1/ */
    /*     datax[i] += 0.1 * meanxh; */
    /*     datay[i] += 0.1 * meanyh; */
    /*     dataxh[i] -= 1.1*meanxh; */
    /*     datayh[i] -= 1.1*meanyh; */
    /* } */

    // Windows wise
    /* for (i = 0; i < size; ++i) { */
    /*     ampx[i] = sqrt( datax[i] * datax[i] + dataxh[i] * dataxh[i]); */
    /*     phasex[i] = atan2(dataxh[i],datax[i]); */
    /*     ampy[i] = sqrt( datay[i] * datay[i] + datayh[i] * datayh[i]); */
    /*     phasey[i] = atan2(datayh[i],datay[i]); */
    /* } */
    /* printf("Start to get phase!!!\n"); */
    int winSize = 36 * 24;
    int steps = size / winSize;
    for (int i = 0; i < steps; ++i)
    {
        int offset = i * winSize;

        double meanxh = gsl_stats_max(dataxh + offset, 1, winSize) + gsl_stats_min(dataxh + offset, 1, winSize);
        double meanyh = gsl_stats_max(datayh + offset, 1, winSize) + gsl_stats_min(datayh + offset, 1, winSize);
        double meanx = gsl_stats_max(datax + offset, 1, winSize) + gsl_stats_min(datax + offset, 1, winSize);
        double meany = gsl_stats_max(datay + offset, 1, winSize) + gsl_stats_min(datay + offset, 1, winSize);

        meanyh /= 2.0;
        meanxh /= 2.0;
        meany /= 2.0;
        meanx /= 2.0;

        /* printf("mid: %lf\t median: %lf\n", meanx, gsl_stats_median(datax+offset, 1, winSize)); */

        /* printf("mid: %lf\t median: %lf\n", meany, gsl_stats_median_from_sorted_data(datay+offset, 1, winSize)); */
        double rate = 0.4;
        meanxh = meanxh * (1 + rate) - rate * gsl_stats_median_from_sorted_data(dataxh + offset, 1, winSize);
        meanyh = meanyh * (1 + rate) - rate * gsl_stats_median_from_sorted_data(datayh + offset, 1, winSize);
        meanx = meanx * (1 + rate) - rate * gsl_stats_median_from_sorted_data(datax + offset, 1, winSize);
        meany = meany * (1 + rate) - rate * gsl_stats_median_from_sorted_data(datay + offset, 1, winSize);

        for (int j = 0; j < winSize; ++j)
        {
            offset = i * winSize + j;
            /* printf("meanx: %lf\n", meanx); */
            datax[offset] -= 1 * meanx;
            datay[offset] -= 1 * meany;
            dataxh[offset] -= 1 * meanxh;
            datayh[offset] -= 1 * meanyh;
            ampx[offset] = sqrt(datax[offset] * datax[offset] + dataxh[offset] * dataxh[offset]);
            phasex[offset] = atan2(dataxh[offset], datax[offset]);
            ampy[offset] = sqrt(datay[offset] * datay[offset] + datayh[offset] * datayh[offset]);
            phasey[offset] = atan2(datayh[offset], datay[offset]);
            /* offset ++; */
        }
        /* for (int j = offset+1; j < size; ++j) { */
        /*     phasex[j] += */
        /* } */
    }

    for (int i = 0; i < size - 1; ++i)
    {
        if (phasey[i] - phasey[i + 1] > 1.5 * M_PI)
        {
            for (int j = i + 1; j < size; ++j)
            {
                phasey[j] += 2 * M_PI;
            }
        }
        if (phasex[i] - phasex[i + 1] > 1.5 * M_PI)
        {
            for (int j = i + 1; j < size; ++j)
            {
                phasex[j] += 2 * M_PI;
            }
        }
    }

    for (int i = 0; i < size; ++i)
    {
        phasex[i] -= phasex[0];
        phasey[i] -= phasey[0];
    }
    for (int i = 1; i < size; ++i)
    {
        if (phasey[i] - phasey[i - 1] > 1.5 * M_PI)
        {
            for (int j = i; j < size; ++j)
            {
                phasey[j] -= 2 * M_PI;
            }
        }
        if (phasex[i] - phasex[i - 1] > 1.5 * M_PI)
        {
            for (int j = i; j < size; ++j)
            {
                phasex[j] -= 2 * M_PI;
            }
        }
    }

    /* for (int i = 0; i < size; ++i) { */
    /*     if (phasex[i] > M_PI ){ */
    /*         for (int j = i; j < size; ++j) { */
    /*             phasex[j] -= 2*M_PI; */
    /*         } */
    /*     } */
    /*     if (phasex[i] < -1*M_PI ){ */
    /*         for (int j = i; j < size; ++j) { */
    /*             phasex[j] += 2*M_PI; */
    /*         } */
    /*     } */
    /*     if (phasey[i] > M_PI ){ */
    /*         for (int j = i; j < size; ++j) { */
    /*             phasey[j] -= 2*M_PI; */
    /*         } */
    /*     } */
    /*     if (phasey[i] < -1*M_PI ){ */
    /*         for (int j = i; j < size; ++j) { */
    /*             phasey[j] += 2*M_PI; */
    /*         } */
    /*     } */
    /* } */
    /* for (int i = 0; i < size; ++i) { */
    /*     double txx = phasex[i]; */
    /*     txx = fmod(txx, 2* M_PI); */
    /*     if (txx > M_PI) */
    /*         txx -= M_PI; */
    /*     else if (txx < -1 *M_PI) */
    /*         txx += M_PI; */
    /*     phasex[i] = txx; */
    /*     txx = phasey[i]; */
    /*     txx = fmod(txx, 2* M_PI); */
    /*     if (txx > M_PI) */
    /*         txx -= M_PI; */
    /*     else if (txx < -1 *M_PI) */
    /*         txx += M_PI; */
    /*     phasey[i] = txx; */
    /* } */

    double maxpx = gsl_stats_max(phasex, 1, size);
    double minpx = gsl_stats_min(phasex, 1, size);
    double maxpy = gsl_stats_max(phasey, 1, size);
    double minpy = gsl_stats_min(phasey, 1, size);
    double mid, rate;

    mid = (maxpx + minpx) / 2.0;
    rate = 2 * M_PI / (maxpx - minpx);
    /* std::cout << "Rate: " << rate << std::endl; */
    for (int i = 0; i < size; ++i)
    {
        /* phasex[i] -= mid; */
        /* phasex[i] *= rate; */
    }

    mid = (maxpy + minpy) / 2.0;
    rate = 2 * M_PI / (maxpy - minpy);
    /* std::cout << "Rate: " << rate << std::endl; */
    for (int i = 0; i < size; ++i)
    {
        /* phasey[i] -= mid; */
        /* phasey[i] *= rate; */
    }

    for (int i = 0; i < size; ++i)
    {
        pd[i] = phasex[i] * n - phasey[i] * m;
    }
    for (int i = 0; i < size - 2; i++)
    {
        phase_diffx[i] = phasex[i + 1] - phasex[i];
        phase_diffy[i] = phasey[i + 1] - phasey[i];
    }
    phase_diffx[size - 2] = phasex[size - 1] - phasex[size - 2];
    phase_diffy[size - 2] = phasey[size - 1] - phasey[size - 2];
}

void Sync::SetIndex(int n, int m)
{
    this->n = n;
    this->m = m;
}

void Sync::phasePlot(std::string name)
{

    double scale = 0.7;

    /* double nsize = 4 * 24 * 6; */
    double nsize = size;
    double *tx, *ty;
    double *x = new double[size];

    for (int i = 0; i < size; i++)
    {
        x[i] = i;
    }

    auto cv = new TCanvas("c0");
    cv->SetTitle(name.c_str());
    cv->SetCanvasSize(900 * scale, 600 * scale);
    /* cv->SetFillColor(0); */
    cv->SetFillStyle(0);
    cv->SetFrameFillStyle(0);

    auto gr1 = new TGraph(size, datax, dataxh);
    gr1->SetTitle("Activity");
    gr1->GetXaxis()->SetTitle("S(t)");
    gr1->GetXaxis()->CenterTitle(1);
    gr1->GetYaxis()->SetTitle("S_{H}(t)");
    gr1->GetYaxis()->CenterTitle(1);

    auto gr2 = new TGraph(size, datay, datayh);
    gr2->SetTitle("Temperature");
    gr2->GetXaxis()->SetTitle("S(t)");
    gr2->GetXaxis()->CenterTitle(1);
    gr2->GetYaxis()->SetTitle("S_{H}(t)");
    gr2->GetYaxis()->CenterTitle(1);

    /* double minx = gsl_stats_min(phasex, 1, size); */
    /* double maxx = gsl_stats_max(phasex, 1, size); */
    /* double miny = gsl_stats_min(phasey, 1, size); */
    /* double maxy = gsl_stats_max(phasey, 1, size); */
    /* maxx = fmax(maxx, maxy); */
    /* minx = fmin(minx, miny); */

    auto gr3 = new TGraph(size, phasex, phasey);
    gr3->SetTitle("Phase");
    tx = gr3->GetX();
    ty = gr3->GetY();

    for (int i = 0; i < size; ++i)
    {
        tx[i] /= M_PI;
        ty[i] /= M_PI;
    }

    gr3->GetXaxis()->SetTitle("#phi_{1} (#pi)");
    gr3->GetXaxis()->CenterTitle(1);
    gr3->GetYaxis()->SetTitle("#phi_{2} (#pi)");
    gr3->GetYaxis()->CenterTitle(1);

    TH1F hist("hist1", "Phase diff dist", 90, -1, 1);
    double txx;
    for (int i = 0; i < size; ++i)
    {
        txx = fmod(pd[i], 2 * M_PI);
        if (txx > M_PI)
            txx -= 2 * M_PI;
        else if (txx < -1 * M_PI)
            txx += 2 * M_PI;
        pd[i] = txx;

        pdcirx[i] = cos(txx);
        pdciry[i] = sin(txx);

        hist.Fill(txx / M_PI, 1);
    }

    hist.Scale(1 / hist.GetEntries());
    hist.GetXaxis()->SetTitle(("#Phi_{" + std::to_string(n) + "," + std::to_string(m) + "} (#pi)").c_str());
    hist.GetXaxis()->CenterTitle(1);

    auto grcir = new TGraph(size, pdcirx, pdciry);
    grcir->SetMarkerStyle(4);
    grcir->SetMarkerSize(1);
    grcir->SetMarkerColorAlpha(kRed, 0.02);

    TH1F shift("", "", 1, -1.2, 1.2);
    shift.SetTitle(("#Phi_{" + std::to_string(n) + "," + std::to_string(m) + "}").c_str());
    shift.SetStats(0);
    shift.GetYaxis()->SetRangeUser(-1.2, 1.2);
    shift.GetXaxis()->SetTitle("x");
    shift.GetXaxis()->CenterTitle(1);
    shift.GetYaxis()->SetTitle("y");
    shift.GetYaxis()->CenterTitle(1);

    auto f1 = new TF1("f1", (std::to_string(((double)n) / m) + "*x").c_str(), -100, 1000);

    auto pad = new TPad("pad", "", .0, .0, 1, 1);
    /* pad->SetFillStyle(4000); */

    int orix = 0, oriy = 0;
    TGraph pointgr(1, &orix, &oriy);
    pad->Draw();

    pad->SetGrid();
    pad->SetFillStyle(4000);

    pad->Divide(3, 2, 0.01, 0.01);

    pad->cd(1)->SetGrid();
    pad->cd(1)->SetLeftMargin(1.1e-1);
    gr1->Draw("AL");
    pointgr.SetMarkerStyle(8);
    pointgr.SetMarkerSize(1);
    pointgr.SetMarkerColor(kRed);
    pointgr.Draw("SAME P");
    pad->cd(2)->SetGrid();
    gr2->Draw("AL");
    pointgr.Draw("SAME P");
    pad->cd(3)->SetGrid();
    pad->cd(3)->SetLeftMargin(1.1e-1);
    gr3->Draw("AP");
    f1->Draw("LSAME");
    pad->cd(4)->SetGrid();
    shift.Draw("p");
    grcir->Draw("P SAME");
    pointgr.Draw("SAME P");
    pad->Modified();
    pad->cd(5);
    hist.Draw("HIST");
    pad->Update();

    cv->Update();
    cv->Print((name + "_hilbert.pdf").c_str(), "Title:Hilbert");
    /* cv->Clear(); */

    delete gr1;
    delete gr2;
    delete gr3;
    delete f1;
    delete pad;
    delete cv;

    cv = new TCanvas("c0");
    cv->SetTitle(name.c_str());
    cv->SetCanvasSize(700 * scale, 1000 * scale);
    pad = new TPad("pad", "", 0, 0, 1, 1);

    int tsize;
    int unit = 36;
    auto gr4 = new TGraph(nsize, x, phasex);
    tsize = gr4->GetN();
    ty = gr4->GetY();
    tx = gr4->GetX();
    for (int i = 0; i < tsize; ++i)
    {
        ty[i] /= M_PI;
        tx[i] /= unit;
    }

    gr4->SetTitle("Activity Phase");
    gr4->GetXaxis()->SetTitle("Time (hour)");
    gr4->GetXaxis()->CenterTitle(1);
    /* gr4->GetXaxis()->SetRangeUser(0, nsize-1); */
    gr4->GetYaxis()->SetTitle("#phi_{1} (#pi)");
    gr4->GetYaxis()->SetTickLength(0.01);
    gr4->GetYaxis()->CenterTitle(1);

    auto gr5 = new TGraph(nsize, x, phasey);
    tsize = gr5->GetN();
    ty = gr5->GetY();
    tx = gr5->GetX();
    for (int i = 0; i < tsize; ++i)
    {
        ty[i] /= M_PI;
        tx[i] /= unit;
    }
    gr5->SetTitle("Temperature Phase");
    gr5->GetXaxis()->SetTitle("Time (hour)");
    gr5->GetXaxis()->CenterTitle(1);
    /* gr5->GetXaxis()->SetRangeUser(0, nsize-1); */
    gr5->GetYaxis()->SetTitle("#phi_{1} (#pi)");
    gr5->GetYaxis()->SetTickLength(0.01);
    gr5->GetYaxis()->CenterTitle(1);

    auto gr6 = new TGraph(nsize - 1, x, phase_diffx);
    ty = gr6->GetY();
    tx = gr6->GetX();
    tsize = gr6->GetN();
    for (int i = 0; i < tsize; ++i)
    {
        ty[i] /= M_PI;
        tx[i] /= unit;
    }
    gr6->SetTitle("Activity Phase incerments");
    gr6->GetXaxis()->SetTitle("Time (hour)");
    gr6->GetXaxis()->CenterTitle(1);
    /* gr6->GetXaxis()->SetRangeUser(-20, nsize-1); */
    gr6->GetYaxis()->SetTitle("#phi_{1} (#pi)");
    gr6->GetYaxis()->SetTickLength(0.01);
    gr6->GetYaxis()->CenterTitle(1);

    auto gr7 = new TGraph(nsize - 1, x, phase_diffy);
    ty = gr7->GetY();
    tx = gr7->GetX();
    tsize = gr7->GetN();
    for (int i = 0; i < tsize; ++i)
    {
        ty[i] /= M_PI;
        tx[i] /= unit;
    }
    gr7->SetTitle("Temperature Phase incerments");
    gr7->GetXaxis()->SetTitle("Time (hour)");
    gr7->GetXaxis()->CenterTitle(1);
    /* gr7->GetXaxis()->SetRangeUser(0, nsize-1); */
    gr7->GetYaxis()->SetTitle("#phi_{1} (#pi)");
    gr7->GetYaxis()->SetTickLength(0.01);
    gr7->GetYaxis()->CenterTitle(1);

    auto gr8 = new TGraph(nsize - 1, x, pd);
    ty = gr8->GetY();
    tx = gr8->GetX();
    tsize = gr8->GetN();
    for (int i = 0; i < tsize; ++i)
    {
        ty[i] /= M_PI;
        tx[i] /= unit;
    }
    gr8->SetTitle("Phase Difference");
    gr8->GetXaxis()->SetTitle("Time (hour)");
    gr8->GetXaxis()->CenterTitle(1);
    /* gr8->GetXaxis()->SetRangeUser(0, nsize-1); */
    gr8->GetYaxis()->SetTitle(("#phi_{" + std::to_string(n) + "," + std::to_string(m) + "} (#pi)").c_str());
    gr8->GetYaxis()->SetTickLength(0.01);
    gr8->GetYaxis()->CenterTitle(1);

    pad->Divide(1, 5, 0.01, 0.01);
    pad->Draw();
    pad->cd(1)->SetGrid(); //->SetTopMargin(0.1);
    gr4->Draw("AL");
    pad->cd(2)->SetGrid(); //->SetBottomMargin(0.15);
    gr5->Draw("AL");
    pad->cd(3)->SetGrid(); //->SetTopMargin(0.1);
    gr6->Draw("AL");
    pad->cd(4)->SetGrid(); //->SetTopMargin(1.1);
    gr7->Draw("AL");
    pad->cd(5)->SetGrid(); //->SetTopMargin(1.1);
    gr8->Draw("AL");

    cv->Update();
    cv->Print((name + "_phase.pdf").c_str(), "Title:Hilbert");
    /* cv->Clear(); */

    delete[] x;

    delete gr4;
    delete gr5;
    delete gr6;
    delete gr7;
    delete gr8;
    delete grcir;
    delete pad;
    delete cv;
}

void Sync::PSI_lambda(int delay)
// PSI ased on MSC:
{
    if (phasex == NULL || phasey == NULL || ampx == NULL || ampy == NULL)
    {
        GetPhase();
    }
    // divide windows
    int phaseWinNum = 90;
    double *lambdar = new double[phaseWinNum];
    double *lambdai = new double[phaseWinNum];
    int *lambdanum = new int[phaseWinNum];

    for (int i = 0; i < phaseWinNum; ++i)
    {
        lambdai[i] = 0;
        lambdar[i] = 0;
        lambdanum[i] = 0;
        assert(lambdai[i] == 0 && "It should be 0");
        assert(lambdar[i] == 0 && "It should be 0");
        assert(lambdanum[i] == 0 && "It should be 0");
    }

    int tm = this->m,
        tn = this->n;
    double gapx = 2 * M_PI * tn;
    double gapy = 2 * M_PI * tm;
    double sgapx = gapx / phaseWinNum;

    double *px = NULL, *py = NULL;
    int tsize = 0;
    if (delay > 0)
    {
        px = phasex + delay;
        py = phasey;
        tsize = size - delay;
    }
    else
    {
        px = phasex;
        py = phasey - delay;
        tsize = size + delay;
    }

    for (int i = 0; i < tsize; ++i)
    {
        double tmpx = px[i];
        double tmpy = py[i];

        // calculate the extra phase to 2pi*m.
        tmpx = fmod(tmpx, gapx);
        if (tmpx < 0)
            tmpx += gapx;
        tmpy = fmod(tmpy, gapy);
        if (tmpy < 0)
            tmpy += gapy;
        assert(tmpx >= 0 && tmpy >= 0 && "Phase should be larger than 0");
        assert(tmpx <= gapx && tmpy <= gapy);

        // calculate in which bin.
        int tmp = (int)(tmpx / sgapx);
        assert(tmp < phaseWinNum && "index should be less than the array size!");

        lambdar[tmp] += cos(tmpy / tm);
        lambdai[tmp] += sin(tmpy / tm);
        lambdanum[tmp] += 1;
    }

    int tmpnum = 0;
    for (int i = 0; i < phaseWinNum; ++i)
    {
        tmpnum += lambdanum[i];
    }

    /* for (int i = 0; i < phaseWinNum; ++i) { */
    /*     printf("%5.2f~%5.2f: [", i*sgap/M_PI, atan2(lambdai[i],lambdar[i])/M_PI); */
    /*     int len = (int )(4000.0 * lambdanum[i] / tmpnum); */
    /*     for (int j = 0; j < len; ++j) { */
    /*         printf("#"); */
    /*     } */
    /*     printf("\n"); */
    /* } */
    if (tmpnum != tsize)
        printf("something disappeared!: tmpnum: %d\t total: %d\n", tmpnum, tsize);
    assert(tmpnum == tsize && "the number in histgram not match the sum");

    double lambda = 0;

    for (int i = 0; i < phaseWinNum; i++)
    {
        // before divide, remember to check if 0
        if (lambdanum[i] == 0)
            continue;
        double ttmp = 1 / (double)lambdanum[i];
        //calculate the average, absolute should be less than 1
        ttmp *= sqrt(lambdar[i] * lambdar[i] + lambdai[i] * lambdai[i]);

        assert(ttmp <= 1.0 && "abs should be less than 1");
        lambda += ttmp;
    }
    lambda /= phaseWinNum;
    assert(lambda <= 1 && "This value thould in [0,1]");
    lambda_index = lambda;
    delete[] lambdai;
    delete[] lambdar;
    delete[] lambdanum;
}

void Sync::PSI_gamma(int delay)
{
    if (phasex == NULL || phasey == NULL || ampx == NULL || ampy == NULL)
    {
        GetPhase();
    }
    gamma_index = 0;

    double *px = NULL, *py = NULL;
    int tsize = 0;
    if (delay > 0)
    {
        px = phasex + delay;
        py = phasey;
        tsize = size - delay;
    }
    else
    {
        px = phasex;
        py = phasey - delay;
        tsize = size + delay;
    }

    double tcos = 0, tsin = 0;
    double diff = 0;
    for (int i = 0; i < tsize; ++i)
    {
        diff = px[i] * n - py[i] * m;
        diff = fmod(diff, 2 * M_PI);
        if (diff < -1 * M_PI)
            diff += 2 * M_PI;
        else if (diff > M_PI)
            diff -= 2 * M_PI;
        tcos += cos(diff);
        tsin += sin(diff);
    }
    tsin /= tsize;
    tcos /= tsize;
    gamma_index = sqrt(tsin * tsin + tcos * tcos);
}

void Sync::PSI_rho(int delay)
// PSI based on Shannon Entropy
{
    if (phasex == NULL || phasey == NULL || ampx == NULL || ampy == NULL)
    {
        GetPhase();
    }
    double *shannon = new double[50];
    gsl_histogram *h = gsl_histogram_alloc(50);
    gsl_histogram_set_ranges_uniform(h, -M_PI, M_PI);
    double diff;
    double *px = NULL, *py = NULL;

    int tsize = 0;
    if (delay > 0)
    {
        px = phasex + delay;
        py = phasey;
        tsize = size - delay;
    }
    else
    {
        px = phasex;
        py = phasey - delay;
        tsize = size + delay;
    }
    for (int i = 0; i < tsize; ++i)
    {
        diff = px[i] * n - py[i] * m;
        diff = fmod(diff, 2 * M_PI);
        if(diff < -1*M_PI)
            diff += 2 * M_PI;
        else if(diff > M_PI)
            diff -= 2 * M_PI;
        gsl_histogram_increment(h, diff);
    }
    for (i = 0; i < 50; ++i) {
        if (delay > 0)
            shannon[i] = gsl_histogram_get(h, i) / (double)(size - delay);
        else
            shannon[i] = gsl_histogram_get(h, i) / (double)(size + delay);
    }
    rho_index = 0;
    for (i = 0; i < 50; ++i) {
        if(shannon[i] == 0)
            continue;
        else
            rho_index -= shannon[i] * log(shannon[i]);
    }
    rho_index = (log(50) - rho_index) / log(50);
    delete[] shannon;
}

double shannonET(double *data, int size)
{
    double SET = 0;
    assert(data != NULL && "On calculate shannon Entropy. Input data points to NULL");
    for (int i = 0; i < size; ++i)
    {
        assert(data[i] > 0 && "Log of negative!");
        if (data[i] == 0)
            continue;
        SET += data[i] * log(data[i]);
    }
    return -1 * SET;
}

void Sync::PSI_decay(std::string name, int t_max, int t_step, std::string method)
{
    int loops = t_max * 2 / t_step;
    double *time = (double *)malloc(loops * sizeof(double));
    double *psidx = (double *)malloc(loops * sizeof(double));
    int timei = -loops / 2.0 * t_step;

    if (method == "rho")
    {
        for (int i = 0; i < loops; i += 1)
        {
            time[i] = timei;
            PSI_rho(timei);
            timei += t_step;
            psidx[i] = rho_index;
        }
    }
    else if (method == "gamma")
    {
        for (int i = 0; i < loops; i += 1)
        {
            time[i] = timei;
            PSI_gamma(timei);
            timei += t_step;
            psidx[i] = gamma_index;
        }
    }
    else if (method == "lambda")
    {
        for (int i = 0; i < loops; i += 1)
        {
            time[i] = timei;
            PSI_lambda(timei);
            timei += t_step;
            psidx[i] = lambda_index;
        }
    }
    else
    {
        printf("PSI_decay(): Error: Wrong method \"%s\", select one from [\"rho\", \"gamma\", \"lambda\"]", method.c_str());
    }

    double unit = 36.0 * 24;
    for (int i = 0; i < loops; i += 1)
    {
        time[i] /= unit;
    }

    TCanvas cv("c0");
    cv.SetCanvasSize(1000, 500);
    cv.SetGrid(20, 10);

    TGraph gr1(loops, time, psidx);
    gr1.SetName("Sync decay");
    gr1.SetTitle(("PSI_{" + std::to_string(n) + "," + std::to_string(m) + "} decay").c_str());
    gr1.GetXaxis()->SetTitle("Time lag (day)");
    gr1.GetXaxis()->CenterTitle();
    gr1.GetYaxis()->SetTitle("PSI (rho)");
    gr1.GetYaxis()->CenterTitle();
    /* gr1.GetXaxis()->SetRangeUser(-15.0, 15.0); */
    gr1.SetLineColor(kRed);
    gr1.Draw("LA");

    // Func to find peaks
    /* int dsize = t_max/2/24+1; */
    /* double *peak = new double[dsize]; */
    /* double *ptime = new double[dsize]; */

    /* int delaay= 24; */

    /* peak[0] = gsl_stats_max((psidx), 1, delaay); */
    /* int idx = gsl_stats_max_index((psidx), 1, delaay); */
    /* ptime[0] = time[idx]; */
    /* for (int i = 0; i < dsize-2; ++i) { */
    /*     peak[i+1] = gsl_stats_max((psidx + i*24*4 + delaay), 1, 24*4); */
    /*     int idx = gsl_stats_max_index((psidx + i*24*4 + delaay), 1, 24*4) + i*24*4 + delaay; */
    /*     ptime[i+1] = time[idx]; */
    /* } */
    /* peak[dsize-1] = gsl_stats_max((psidx + (dsize-2)*24*4 + delaay), 1, 24*4-delaay); */
    /* idx = gsl_stats_max_index((psidx + (dsize-2)*24*4 + delaay), 1, 24*4-delaay) + (dsize-2)*24*4 + delaay; */
    /* ptime[dsize-1] = time[idx]; */

    double mean = 0, max = 0, std = 0;
    max = gsl_stats_max(psidx, 1, loops);
    mean = gsl_stats_mean(psidx, 1, loops);
    std = gsl_stats_sd(psidx, 1, loops);

    decay_radio = max - mean;
    decay_radio /= std;

    decay_max_time = time[(int)gsl_stats_max_index(psidx, 1, loops)];

    /* TGraph gr2(31, ptime, peak); */
    /* gr2.SetName("local max"); */
    /* gr2.SetMarkerStyle(8); */
    /* gr2.SetMarkerSize(1); */
    /* gr2.Draw("PSAME"); */

    auto legend = new TLegend(0.75, 0.8, 0.9, 0.9);
    legend->AddEntry("Sync decay", ("psi_{max}=" + std::to_string(gsl_stats_max(psidx, 1, loops))).c_str(), "l");
    legend->AddEntry("", ("  t_{lag}=" + std::to_string(decay_max_time)).c_str(), "");
    legend->AddEntry("", ("W=" + std::to_string(decay_radio)).c_str(), "");
    /* legend->AddEntry("local max", "maximum", "p"); */
    legend->Draw();

    cv.Update();
    cv.Print((name + "_psi_" + method + "_decay.pdf").c_str(), "Title:PSI decay");
    cv.Clear();

    free(time);
    free(psidx);
}

void Sync::solver()
{
    GetPhase();
}

void Sync::plot(std::string name, int seg_size)
{
    double diff;
    GnuplotPipe gp; 
    gp.sendLine("set terminal pdf size 5,4");
    gp.sendLine("set grid");

    gp.sendLine("set output '"+name+"_phase_diff.pdf'");
    /* gp.sendLine("set xrange [0:"+ std::to_string(size) +"]"); */
    gp.sendLine("plot '-' w points ps 0.5 title 'x'"); //" , '-' w lines lw 0.5 title 'y'");
    /* for(i =0;i<size;i++){ */
    /*     gp.sendLine(std::to_string(i) + "\t" + std::to_string(phasex[i])); */
    /* } */
    /* gp.sendLine("e"); */
    for(i =0;i<size;i++){
        diff = phasex[i]-phasey[i];
        if(diff < -1*M_PI)
            diff += M_PI;
        else if(diff > M_PI)
            diff -= M_PI;
        gp.sendLine(std::to_string(phasex[i]) + "\t" + std::to_string(phasey[i]));
    }
    gp.sendLine("e");
    gp.sendLine("set output");

    gp.sendLine("set output '"+name+"_Sync.pdf'");
    gp.sendLine("set autoscale");
    gp.sendLine("bin(x,s) = s*int(x/s)");
    gp.sendLine("set boxwidth 0.1");
    gp.sendLine("set xrange [-pi:pi]");
    gp.sendLine("set style fill solid 0.5");
    gp.sendLine("set xlabel 'phase difference'");
    gp.sendLine("plot '-'  u (bin($1-$2,0.1)) smooth frequency w boxes title 'index = " + std::to_string(rho_index).substr(0, 5) + "'");
    for(i = 0;i < size ; i++){
        diff = phasex[i]-phasey[i];
        if(diff < -1*M_PI)
            diff += M_PI;
        else if(diff > M_PI)
            diff -= M_PI;
        gp.sendLine(std::to_string(diff) + "\t0");
    }
    gp.sendLine("e");
    gp.sendLine("set output");
}

void hilbert_trans(double *in, double *output,int num)
{
    int i=0;
    fftw_complex *out;
    out = fftw_alloc_complex(num);
    for(i=0;i<num;i++){
        out[i][RE] = in[i];
        out[i][IM] = 0.0;
    }
    fftw_plan plan;
    plan = fftw_plan_dft_1d(num,out,out,FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    int hN = num>> 1;
    int numRem = hN; 
    for(i=1;i<hN;i++){
        out[i][RE] *= 2;
        out[i][IM] *= 2;
    }

    if(NMICE%2 == 0)
        numRem --;
    else if(num>1){
        out[hN][RE] *= 2;
        out[hN][IM] *= 2;
    }
    for(i=hN+1;i<num;i++){
        out[i][RE] = 0;
        out[i][IM] = 0;
    }

    plan = fftw_plan_dft_1d(num,out,out,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    for(i=0;i<num;i++){
        output[i] = out[i][IM]/num;
    }
    fftw_free(out);
}

int fft(fftw_complex *in, fftw_complex *out, int num)
{
    fftw_plan p; 
    p = fftw_plan_dft_1d(num,in,out,-1,FFTW_ESTIMATE); 
    fftw_execute(p);
    fftw_destroy_plan(p);
    return 1;
}

int ifft(fftw_complex *in, fftw_complex *out, int num)
{
    fftw_plan p; 
    p = fftw_plan_dft_1d(num,in,out,+1,FFTW_ESTIMATE); 
    fftw_execute(p);
    fftw_destroy_plan(p);
    for(int i=0;i<num;i++){
        out[i][RE] /= num;
        out[i][IM] /= num;
    }
    return 1;
}

DFAI::DFAI(double *data, int data_size, int dfa_order): size(data_size), order(dfa_order)
{
    rito = pow(2.0, 1.0 / 5);
    minbox = 2*(order+1);
    maxbox = size / 2;
    size_o = log10(maxbox / (double)minbox)/log10(rito) + 1.5;

    // allocated memory
    data_xi = new double[size];
    data_yi = new double[size];
    data_x = new double[size];
    data_y = new double[size];
    diff = new double[size];
    dfax = new double[size_o];
    dfay = new double[size_o];
    seg_size = new int[size_o];

    mkbox();
    init(data);
    solver();
};

DFAI::~DFAI()
{
    delete[] dfax;
    delete[] dfay;
    delete[] data_xi;
    delete[] data_yi;
    delete[] data_x;
    delete[] data_y;
    delete[] diff;
    delete[] seg_size;
};

void DFAI::init(double *data)
{
    // First step: To get the profile
    // remove and integral 
    double rm = gsl_stats_mean(data, 1, size);
    diff[0] = 0.0;
    data_y[0] = data[0] - rm;
    data_x[0] = 0;
    for (i = 1; i < size; ++i) {
        diff[i] = 0.0;
        data_y[i] = data[i] - rm + data_y[i-1];
        data_x[i] = i;
    }
    for (int i = 0; i < size; ++i) {
        data_yi[i] = data_y[size-i-1];
        data_xi[i] = data_x[size-i-1];
    }
    for (i = 0; i < size_o; ++i) {
        dfax[i] = 0.0;
        dfay[i] = 0.0;
    }
}

void DFAI::mkbox()
{
    for(i = 1, j = 1, seg_size[0] = minbox; j < size_o && seg_size[j-1] < maxbox; ++i)
        if((seg_size[j] = minbox * pow(rito,i) + 0.5) > seg_size[j-1])
            j++;
    size_o = j;
}

void DFAI::recal(double *data)
{
    init(data);
    solver();
}

void DFAI::solver()
{
    // init the window info
    int end = 0, tail = 0; 
    Polyfit new_fit( data_x, data_y, 10, order+1);

    double diff_temp = 0;

    for (i = 0; i < size_o; i++) { 
        end = size - seg_size[i];
        new_fit.resize(seg_size[i]);
        for (j=0; j < end; j+= seg_size[i]){
            new_fit.recal(data_x+j, data_y+j);
            tail = j+seg_size[i];
            for(k=j; k< tail; k++){
                diff_temp = new_fit.fitfunc(data_x[k]) - data_y[k];
                /* diff[k] = new_fit.fitfunc(data_x[k]) - data_y[k]; */
                diff[k] = diff_temp*diff_temp;
            }
        }
        dfay[i] = gsl_stats_mean(diff, 1, end);
        for (j=0; j < end; j+= seg_size[i]){
            new_fit.recal(data_xi+j, data_yi+j);
            tail = j+seg_size[i];
            for(k=j; k< tail; k++){
                diff_temp = new_fit.fitfunc(data_xi[k]) - data_yi[k];
                diff[k] = diff_temp*diff_temp;
                /* diff[k] /= 2.0; */
                /* diff[k] = new_fit.fitfunc(data_xi[k]) - data_yi[k]; */
                /* diff[k] = diff[k]*diff[k]; */
            }
        }
        dfay[i] += gsl_stats_mean(diff, 1, end);
        dfay[i] /= 2;
        dfay[i] = sqrt(dfay[i]);
        dfax[i] = seg_size[i];
        // log lize
        dfay[i] = log10(dfay[i]);
        dfax[i] = log10(dfax[i]);
    }

    {
        // fit loged function with least linear regression
        Polyfit new_fit(dfax,dfay,size_o,2);
        coeff = new_fit.get_rsq();
        store[0] = new_fit.store[0]; 
        store[1] = new_fit.store[1];
        index = store[1];
    }

    for (int i = 0; i < size_o; ++i) {
        dfay[i] = pow(10.0,dfay[i]);
        dfax[i] = pow(10.0,dfax[i]);
    }

    /* func = std::to_string(store[0]) + "+x*" + std::to_string(store[1]); */
    func = "pow(10.0,"+std::to_string(store[0]) + ") *pow(x,"+std::to_string(store[1]) +")"  ;
}

void DFAI::plot_over(int point)
{
    GnuplotPipe gp; 
    gp.sendLine("set terminal pdf size 5,4");
    gp.sendLine("set output 'DFA_" + std::to_string(point) + ".pdf'");
    gp.sendLine("set title 'DFA'");
    gp.sendLine("set logscale x 10");
    gp.sendLine("plot '-' w points ps 1.5 notitle" );
    for (i = 0; i < size_o; ++i) {
        gp.sendLine(std::to_string(dfax[i]) + "\t" + 
                std::to_string(dfay[i]) );
    }
    gp.sendLine("e");
    gp.sendLine("set output");
}

double DFAI::fitfunc(double input) 
{
    return store[0] + store[1]*input;
}
