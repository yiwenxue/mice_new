#include <mathematics.h>

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
    rito = pow(2.0,1.0/8);
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
    func = std::to_string(store[0]) + " *x** " + std::to_string(store[1]) ;
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
}

void Sync::solver()
{
    hilbert_trans(datax, dataxh, size);
    hilbert_trans(datay, datayh, size);


    double *shannon = new double[50];
    gsl_histogram *h = gsl_histogram_alloc(50);
    gsl_histogram_set_ranges_uniform(h, -M_PI, M_PI);
    double diff;
    for (i = 0; i < size; ++i) {
        ampx[i] = sqrt( datax[i] * datax[i] + dataxh[i] * dataxh[i]);
        phasex[i] = atan2(dataxh[i],datax[i]);
        ampy[i] = sqrt( datay[i] * datay[i] + datayh[i] * datayh[i]);
        phasey[i] = atan2(datayh[i],datay[i]);
        diff = phasex[i] - phasey[i];
        /* std::cout << ampx[i] <<"\t"<< ampy[i] <<"\t"<< diff << std::endl; */
        if(diff < -1*M_PI)
            diff += M_PI;
        else if(diff > M_PI)
            diff -= M_PI;
        gsl_histogram_increment(h, diff);
    }

    for (i = 0; i < 50; ++i) {
        shannon[i] = gsl_histogram_get(h,i)/(double)size;
    }

    sync_index = 0;
    for (i = 0; i < 50; ++i) {
        if(shannon[i] == 0)
            continue;
        else 
            sync_index -= shannon[i] * log(shannon[i]);
    }
    sync_index = (log(50) - sync_index ) / log(50);

    delete[] shannon;
}

void Sync::plot(std::string name, int seg_size)
{
    double diff;
    GnuplotPipe gp; 
    gp.sendLine("set terminal pdf size 5,4");
    gp.sendLine("set grid");

    gp.sendLine("set output '"+name+"_phase_diff.pdf'");
    gp.sendLine("set xrange [0:"+ std::to_string(size) +"]");
    gp.sendLine("plot '-' w lines lw 0.5 title 'x'");//" , '-' w lines lw 0.5 title 'y'");
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
        gp.sendLine(std::to_string(i) + "\t" + std::to_string(diff));
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
    gp.sendLine("plot '-'  u (bin($1-$2,0.1)) smooth frequency w boxes title 'index = "+ std::to_string(sync_index).substr(0,5) +"'");
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

