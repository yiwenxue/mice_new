#include <TCanvas.h>
#include <TGraph.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <utility.h>

/// @details A function to get number of a file.
int fileLines(std::string path)
{
    FILE *input = fopen(path.c_str(),"r");
    if(input == NULL)
    {
        fprintf(stderr,"(File Lines)Failed to read %s\n",path.c_str());
        return -1;
    }

    //Count the totla lines;
    char chr;
    int lines = 0;
    chr = fgetc(input);
    while(chr != EOF)
    {
        if(chr == '\n')
            lines ++;
        chr = fgetc(input);
    }
    // Don't forget the last line.
    lines ++; 
    fclose(input);
    return lines;
}

/// @details 
int TMice::GetDataRange(RangeIndex _index)
{
    if( _index == rSTART )
        return this->start;
    else if( _index == rEND )
        return (this->start + this->length );
    else if ( _index == rLENGTH ) 
        return this->length;
    else return -1;
}

int  TMice::SetParRange(double *_parMin, double *_parMax)
{
    for(int i=0;i<4;i++)
    {
        this->parMax[i] = _parMax[i];
        this->parMin[i] = _parMin[i];
    }
    return 0;
}

double TMice::GetParRange(int i, int j)
{
    double tmp;
    if(i==0)
        tmp = parMin[j];
    else if(i==1)
        tmp = parMax[j];
    else 
        tmp = NAN;
    return tmp;
}

int TMice::PrintParRange()
{
    for (int i = 0; i < 4; ++i) {
        std::cout << "p" << i << "Min = " << this->parMin[i] 
            << "\tp" << i << "Max = " << this->parMax[i] 
            << std::endl;
    }
    return 0;
}

bool TMice::Average()
{
    if(actAve.data == nullptr)
    {
        delete[]  actAve.data;
        delete[]  actAve.time;
    }
    if(actStd.data == nullptr)
    {
        delete[]  actStd.data;
        delete[]  actStd.time;
    }
    if(tempAve.data == nullptr)
    {
        delete[] tempAve.data;
        delete[] tempAve.time;
    }
    actAve.size = winNum;
    actAve.data = new double[winNum];
    actAve.time = new double[winNum];
    actStd.size = winNum;
    actStd.data = new double[winNum];
    actStd.time = new double[winNum];
    tempAve.size = winNum;
    tempAve.data = new double[winNum];
    tempAve.time = new double[winNum];

    for (int i = 0; i < winNum; ++i) {
        actAve.data[i] = gsl_stats_mean(act.data+i*winSize, 1, winSize);
        actStd.data[i] = gsl_stats_variance(act.data+i*winSize, 1, winSize);
        tempAve.data[i] = gsl_stats_mean(temp.data+i*winSize, 1, winSize);
        actAve.time[i] = act.time[i*winSize];
        actStd.time[i] = actAve.time[i];
        tempAve.time[i] = temp.time[winSize*i];
    }
    return true;
}

int TMice::PrintDetails()
{
    std::cout << "=====Mice details======" << std::endl;
    std::cout << "Name:\t" << mice.name << std::endl;
    std::cout << "Path:\t" << mice.path + "/batch" << mice.batch << "/"+mice.name << ".Temperature.txt" << std::endl;
    std::cout << "Path:\t" << mice.path + "/batch" << mice.batch << "/"+mice.name << ".Activity.txt" << std::endl;
    std::cout << "Batch:\tBatch" << mice.batch << std::endl;
    std::cout << "Mutant:\t" << (mice.ifmutant?"True":"False") << std::endl;
    std::cout << "=====Window details====" << std::endl;
    std::cout << "Window Size:\t"   << this->winSize << std::endl;
    std::cout << "Window Num:\t"      << this->winNum << std::endl;
    std::cout << "temp size:\t" << temp.size << std::endl;
    std::cout << "tempAve size:\t" << tempAve.size << std::endl;
    std::cout << "act size:\t" << act.size << std::endl;
    std::cout << "actAve size:\t" << actAve.size << std::endl;
    std::cout << "=====Fit Details======="     << std::endl;
    std::cout << "Fit Window Size:\t" << this->fitSize << std::endl;
    std::cout << "Fit Window Stride:\t"<<this->fitStride  << std::endl;
    PrintParRange();
    std::cout << "=====Data errors=======" << std::endl;
    std::cout << "Act  Errors           : " << this->errors_act  << std::endl;
    std::cout << "Act  Total Number     : " << this->actOri.size << std::endl;
    std::cout << "Act  Errors Percentage: " << GetErrorPercentage(0) << std::endl;
    std::cout << "Temp Errors           : " << this->errors_temp  << std::endl;
    std::cout << "Temp Total Number     : " << this->tempOri.size << std::endl;
    std::cout << "Temp Errors Percentage: " << GetErrorPercentage(1) << std::endl;
    std::cout << "======================\n" << std::endl;

    return 0;
}

/// @brief A function to read in the mice tag db(A index file).
/// @details The db stores the name, batch, path, ifmutant, start, length, and winSize;
int Micedb(std::string infile)
{
    int num =0;
    return num;
}

int TMice::SetfitSize(int _size)
{
    this->fitSize = _size;
    return 0;
}
int TMice::GetfitSize()
{
    return this->fitSize;
}

int TMice::SetfitStride(int _step)
{
    this->fitStride = _step;
    return 0;
}
int TMice::GetfitStride()
{
    return this->fitStride;
}

double TMice::GetfitChilimit()
{
    return this->fitChilimit;
}
int TMice::SetfitChilimit(double _chi)
{
    fitChilimit = _chi;
    return 0;
}

int TMice::GetPeriod()
{
    TF1 *fit1 = nullptr;
    int itmax = (tempAve.time[tempAve.size-1] - fitSize - tempAve.time[0])/fitStride; // how many times will the fit takeplace;
    itmax ++; // Must plus 1 because from 0 to itmax is itmax +1 

    FitCosinor *result = new FitCosinor[itmax];
    grmice = new TGraph(tempAve.size, tempAve.time, tempAve.data); 
    canvas->cd();

    /// @details Perform the nonlinear fit at all regions
    for(int i=0;i<itmax;i++)
    {
        fit1 = new TF1("fit1", "[0]*sin([1]*x+[2])+[3]",tempAve.time[0] + i*fitStride 
                , tempAve.time[0] + fitSize + i*fitStride);
        for(int j=0;j<4;j++)
        {
            fit1->SetParLimits(j, parMin[j],parMax[j]);
            fit1->SetParameter(j,(parMin[j]+parMax[j])/2.0);
        }
        for (int j = 0; j < 10; ++j) 
        {
            /// the fit details are disabled, if you want to see the details, please just remove "Q" 
            grmice->Fit(fit1,"RQ");
            result[i].chisq = fit1->GetChisquare();
            if(result[i].chisq < fitChilimit)
                break;
        }
        result[i].p = new double[4];
        for (int j = 0; j < 4; ++j) 
        {
            result[i].p[j] = fit1->GetParameter(j);
        }
        result[i].ndf = fit1->GetNDF();
        if(result[i].chisq > fitChilimit)
        {
            grmice->GetXaxis()->SetRangeUser(tempAve.time[0] + i*fitStride,tempAve.time[0] + fitSize + i*fitStride);
            grmice->Draw("");
            canvas->Update();
            canvas->Print((mice.name+"_badpoint_batch"+std::to_string(mice.batch)+"["+ std::to_string(tempAve.time[0] + i*fitStride)
                        +","+ std::to_string(tempAve.time[0] + fitSize + i*fitStride) +"].pdf").c_str(),"Title:errors");
            canvas->Clear();
            std::cout << "Attention! Bad fit appear at ["
                << tempAve.time[0] + i*fitStride << "," << tempAve.time[0] + fitSize + i*fitStride
                << "].\n###################################################" << std::endl;
            /* return -1; */
        }
        delete fit1;
    }
    /* std::cout << "Fit finished" << std::endl; */

    /// Arrays to store the period information temporarily
    if(period.data != nullptr){
        delete [] period.data;
        delete [] period.time;
    }
    period.data = new double[itmax];
    period.time = new double[itmax];
    period.size = itmax;

    for (int i = 0; i < itmax; ++i) {
        period.data[i] = 2*M_PI/result[i].p[1];
        period.time[i] = tempAve.time[0] + fitSize/2 + i*fitStride;
        /* std::cout << "index: " << i << "\ttime: " << period.time[i] << "\tperiod:" << period.data[i] << std::endl; */
    }

    /// @details This loop will get a average period. 
    /* int head, it, top; */
    /* for (int i = 0; i < tempAve.size ; i+= fitStride*360/winSize) // improvable, i is the index of */ 
    /* { */
    /*     head = int(tempAve.time[i] - tempAve.time[0])/fitStride; */
    /*     // Use the head minus head or tail minus tail. you must use the same part. */
    /*     top = head - int(tempAve.time[tempAve.size-1] - tempAve.time[0] - fitSize + 0.5)/fitStride; */

    /*     // Top is used to detect the end part. */ 
    /*     // it is used to solve the start condation, and can also record the number of periods to average */ 
    /*     if(top < 0) */
    /*         top = 0; */ 
    /*     period.data[i] = 0; */
    /*     it=top; */
    /*     while(it < (fitSize/fitStride) && head-it >= 0) */ 
    /*     { */
    /*         period.data[i] += 2*M_PI/result[i/(fitStride*360/winSize) - it].p[1]; */ 
    /*         it++; */
    /*     } */
    /*     period.data[i] /= (it -top); */ 
    /*     for (int j = 0; j < fitStride*360/winSize ; ++j) */ 
    /*     { */
    /*         period.data[i + j] = period.data[i]; */
    /*         /1* std::cout << "index: " << i+j << "\t period: " << period.data[i] << std::endl; *1/ */
    /*     } */
    /* } */
    /* std::cout << "Average the period" << std::endl; */

    for (int i = 0; i < itmax; ++i) 
    {
        delete[] result[i].p; 
    }
    delete[] result;
    delete grmice;
    return 0;
}

bool TMice::RhythmRemap()
{

    /* // This will turn the time into phase. */
    /* tempAve.time[0] = 0;// this->winSize/360.0/period[0]*2*M_PI; */
    /* actAve.time[0] = 0;// this->winSize/360.0/period[0]*2*M_PI; */

    std::cout << "=====Remap the data to phase space====="  << std::endl;
    int remapdays, phaseDsize;
    double total_phase = 0;
    for (int i = 0; i < period.size; ++i) {
        total_phase += fitStride * 2 *M_PI / period.data[i];
    }
    remapdays = int(total_phase / 2 / M_PI);
    phaseDsize = remapdays * (24 * 360 / winSize);
    std::cout << "Total phase is: " << total_phase << std::endl;
    std::cout << "Total days is: " << remapdays << std::endl;
    std::cout << "Total phase numbers: " << phaseDsize << std::endl;

    /* for (int i = 1; i < tempAve.size; ++i) */ 
    /* { */
    /*     tempAve.time[i] = this->winSize/360.0/period.data[i]*2*M_PI + this->tempAve.time[i-1]; */
    /*     actAve.time[i] = this->winSize/360.0/period.data[i]*2*M_PI + this->actAve.time[i-1]; */
    /*     period.time[i] = tempAve.time[i] / 2.0 /M_PI; */
    /* } */
    /* std::cout << "Get the phase" << std::endl; */

    /// Because the time in tempave is not eaqual distributed, I need a resample.
    if (tempRem.data != nullptr ) 
    {
        delete [] tempRem.data;
        delete [] tempRem.time;
    }
    if (actRem.data != nullptr ) 
    {
        delete [] actRem.data;
        delete [] actRem.time;
    }

    tempRem.size=phaseDsize;
    tempRem.time = new double[tempRem.size];
    tempRem.data = new double[tempRem.size];
    actRem.size=phaseDsize;
    actRem.time = new double[actRem.size];
    actRem.data = new double[actRem.size];

    /* for (int i = 1, j = 1; i < actAve.size;) */
    /* { */
    /*     while( j *2*M_PI/(360/this->winSize)/24 < actAve.time[i]) */
    /*     { */
    /*         actRem.time[j] = j; */ 
    /*         actRem.data[j] = actAve.data[i-1]+(actAve.data[i]-actAve.data[i-1])/(actAve.time[i]-actAve.time[i-1])*(j*2*M_PI/(360/this->winSize)/24 - actAve.time[i-1]); */  
    /*         j++; */ 
    /*     } */
    /*     i++; */
    /* } */
    /* int lasti =0; */ 
    /* for (int i = 1, j = 1; i < tempAve.size;) */ 
    /* { */
    /*     while( j *2*M_PI/(360/this->winSize)/24 < tempAve.time[i]) */
    /*     { */
    /*         tempRem.time[j] = j; */ 
    /*         tempRem.data[j] = tempAve.data[i-1]+(tempAve.data[i]-tempAve.data[i-1])/(tempAve.time[i]-tempAve.time[i-1])*(j*2*M_PI/(360/this->winSize)/24 - tempAve.time[i-1]); */  
    /*         /1* std::cout << tempRem.size <<"\t"<< i << "\t" << j <<"\t"<< tempAve.time[i-1] <<"\t"<< j*2*M_PI/(360/this->winSize)/24 <<"\t"<< tempAve.time[i] << std::endl; *1/ */
    /*         j++; */ 
    /*         lasti = i; */
    /*     } */
    /*     i++; */
    /* } */

    std::cout << std::endl;
    std::cout << "======================================="  << std::endl;
    return true;
}

//
int TMice::DrawHeatmap()
{
    auto *hist = new TH2F("hits",("Heatmap of Temperature("+mice.name+")").c_str()
            ,360/this->winSize*48, 0, 360/this->winSize*48, this->length/48-1,0,this->length/24-2);
    for (int i = 0; i < this->tempRem.size; ++i)
    {
        hist->Fill(int(tempRem.time[i])%(360/this->winSize*48), i/(360/this->winSize*48)*2, tempRem.data[i]);
    }
    canvas->Clear();
    hist->SetStats(kFALSE);
    hist->Draw("COLZ");
    canvas->Print((mice.name+"_Temperature.pdf").c_str(),"Title:Heatmap of Temperature");
    delete hist;
    canvas->Clear();
    hist = new TH2F("hits",("Heatmap of Activity("+mice.name+")").c_str()
            ,360/this->winSize*24, 0, 360/this->winSize*48, this->length/48,0,this->length/24);
    for (int i = 0; i < this->actRem.size; ++i)
    {
        hist->Fill(int(actRem.time[i])%(360/this->winSize*48), i/(360/this->winSize*48)*2, actRem.data[i]);
    }
    canvas->Clear();
    hist->SetStats(kFALSE);
    hist->Draw("COLZ");
    canvas->Print((mice.name+"_Activity.pdf").c_str(),"Title:Heatmap of Activity");
    delete hist;
    canvas->Clear();
    return 0;
}


int TMice::DrawPeriodDist(void )
{

    int temp_size;
    temp_size = period.size; // the length for all data
    std::cout << "period length: " << temp_size ;
    temp_size = length / fitStride - 1; // the length for selected
    std::cout << "Selected length: " << temp_size << std::endl;
    

    double *avep = new double[temp_size];
    double *avet = new double[temp_size];
    for (int i = 0; i < temp_size; ++i) {
        avep[i] = 0;
        avet[i] = 0;
    }

    int j=0;
    for (int i = 0; i < temp_size ; ++i) {
        if( i+1 < fitSize/fitStride )
            j = i+1;
        else if( i+fitSize/fitStride > temp_size )
            j = temp_size - i;
        else j=fitSize/fitStride;
        avep[i] = gsl_stats_mean(period.data+i-j/2, 1, j);
        avet[i] = gsl_stats_mean(period.time+i-j/2, 1, j);
        avet[i] /= 24;
    }

    TF1 *f1 = new TF1("f1","24", 0 , period.time[temp_size-1]);

    /* grmice = new TGraph(temp_size, period.time, period.data); */
    grmice = new TGraph(temp_size, avet, avep);
    grmice->SetTitle(("Period distribution(" +mice.name+ ")").c_str());
    grmice->GetYaxis()->CenterTitle();
    grmice->GetYaxis()->SetTitle("Period(hours)");
    grmice->GetXaxis()->CenterTitle();
    grmice->GetXaxis()->SetTitle("Time(days)");
    /* canvas->cd(); */
    grmice->Draw("ALP");
    f1->Draw("SAME");
    /* std::cout << "Plot diagram to canvas" << std::endl; */
    /* f1->Draw("SAME"); */
    canvas->Update();
    canvas->Print((mice.name+"_period.pdf").c_str(),"Title:Period");
    canvas->Clear();
    /* std::cout << "Period diagram saved" << std::endl; */

    delete f1;
    delete grmice;
    delete [] avep;
    delete [] avet;

    return 0;
}

int TMice::DrawOverview()
{
    grmice = new TGraph(tempAve.size, tempAve.time, tempAve.data);
    grmice->SetTitle(("Temperature overvire(" +mice.name+ ")" ).c_str());
    grmice->GetYaxis()->CenterTitle();
    grmice->GetYaxis()->SetTitle("Temperatue(C)");
    grmice->GetXaxis()->CenterTitle();
    grmice->GetXaxis()->SetTitle("Time(hours)");

    grmice->Draw();
    canvas->Update();
    canvas->Print((mice.name+"_temperature_Overview.pdf").c_str(),"Title:temperature");
    canvas->Clear();
    delete grmice;

    grmice = new TGraph(actAve.size, actAve.time, actAve.data);
    grmice->SetTitle(("Activity overvire(" +mice.name+ ")").c_str());
    grmice->GetYaxis()->CenterTitle();
    grmice->GetYaxis()->SetTitle("Activity");
    grmice->GetXaxis()->CenterTitle();
    grmice->GetXaxis()->SetTitle("Time(hours)");

    grmice->Draw();
    canvas->Update();
    canvas->Print((mice.name+"_activity_Overview.pdf").c_str(),"Title:activity");
    canvas->Clear();
    delete grmice;
    return 0;
}


int DFA_plot(double *data,
        int data_size,
        int dfa_order,
        std::string name,
        struct period_dfa *input)
{
    auto newdfa = new DFA(data, data_size, dfa_order);
    int num = newdfa->get_size_o();
    double *xp = new double[num];
    double *yp = new double[num];
    double pars[5];

    for (int i = 0; i < num; ++i) {
        xp[i] = newdfa->dfax[i];
        yp[i] = newdfa->dfay[i];
    }

    yp[0] = log10(newdfa->dfay[0]);
    xp[0] = log10(newdfa->dfax[0]);
    for (int i = 1; i < num; ++i) {
        yp[i] = log10(newdfa->dfay[i]);
        xp[i] = log10(newdfa->dfax[i]);
        std::cout << xp[i] - xp[i-1] << std::endl;
    }
    auto *canvas = new TCanvas("c0");
    auto f1 = new TF1("f1",fit_crossover,xp[0],xp[num-1],5);
    auto *plot = new TGraph(num, xp, yp);

    f1->SetParLimits(0, 0, 2);
    f1->SetParameter(0, 1.0);
    f1->SetParLimits(2, 0, 2);
    f1->SetParameter(2, 0.5);
    f1->SetParLimits(4, xp[4], xp[num-4]);
    f1->SetParameter(4, (xp[(int)(num * 0.6)]) );
    plot->Fit(f1,"QR");
    delete plot;

    std::cout << "==DFA"<< dfa_order <<" F~s with piecewise fit" << std::endl;
    for (int i = 0; i < 5; ++i) {
        pars[i] = f1->GetParameter(i);
        std::cout << "Pars"<< i<< ": " << pars[i] << std::endl;
    }

    auto f2 = new TF1("f2",draw_crossover,newdfa->dfax[0], pow(10,pars[4]),5);
    auto f3 = new TF1("f3",draw_crossover,pow(10,pars[4]),newdfa->dfax[num-1],5);
    for (int i = 0; i < 4; ++i) {
        f2->SetParameter(i, pars[i]);
        f3->SetParameter(i, pars[i]);
    }
    f2->SetParameter(4, pow(10, pars[4]));
    f3->SetParameter(4, pow(10, pars[4]));

    double rvalue, lvalue;
    lvalue = pow(10, pars[4]);
    for (int i = 0; i < num; ++i) {
        if ( lvalue > newdfa->dfax[i] )
            continue;
        else
        {
            lvalue = newdfa->dfax[i-1];
            rvalue = newdfa->dfax[i];
            break;
        }
    }
    
    char crosslabel[255]; 
    sprintf(crosslabel, "Crossover:%.2fd-%.2fd-%.2fd",lvalue * 0.010416667,pow(10.0,pars[4]) * 0.010416667,rvalue * 0.010416667);
    input->lbp = lvalue*0.010416667;
    input->rbp = rvalue*0.010416667;
    input->midbp = pow(10.0,pars[4])*0.010416667;
    input->index_l = pars[0];
    input->index_r = pars[2];

    plot = new TGraph(num, newdfa->dfax, newdfa->dfay);
    plot->SetName("period_dfa");
    plot->SetMarkerStyle(8);
    plot->SetMarkerSize(0.3);
    plot->SetTitle(("DFA("+std::to_string(dfa_order)+") of " + name).c_str());
    plot->GetYaxis()->CenterTitle();
    plot->GetYaxis()->SetTitle("Fn");
    plot->GetXaxis()->CenterTitle();
    plot->GetXaxis()->SetTitle("N");
    canvas->SetLogx(1);
    canvas->SetLogy(1);

    f2->SetLineColor(kRedBlue);
    f3->SetLineColor(kPink);

    plot->Draw("AP");
    f2->Draw("SAME");
    f3->Draw("SAME");

    auto legend = new TLegend(0.1,0.8,0.4,0.9);
    legend->AddEntry("period_dfa","DFA data","l");
    legend->AddEntry("f2",("Exp: " + std::to_string(pars[0])).c_str(),"l");
    legend->AddEntry("f3",("Exp: " + std::to_string(pars[2])).c_str(),"l");
    legend->AddEntry("",crosslabel,"l");
    legend->Draw();

    canvas->Update();
    canvas->Print(("dfa_sketch_" + name + ".pdf").c_str(),"Title:sketch");
    canvas->Clear();
    canvas->SetLogx(0);
    canvas->SetLogy(0);

    delete[] yp;
    delete[] xp;
    delete f1;
    delete f2;
    delete f3;

    delete plot;
    delete canvas;

    return 0;
}

int TMice::PeriodDFA(int dfa_order)
{

    std::cout << "=====DFA"<< dfa_order <<" for periods=======" << std::endl;
    start_index = start;
    auto newdfa = new DFAI(period.data+start_index,(length-fitStride)/fitStride,dfa_order);

    int num = newdfa->get_size_o();
    double *yp = new double[num];
    double *xp = new double[num];
    double pars[5];

    for (int i = 0; i < num; ++i) {
        yp[i] = log10(newdfa->dfay[i]);
        xp[i] = log10(newdfa->dfax[i]);
    }
    auto f1 = new TF1("f1",fit_crossover,xp[0],xp[num-1],5);
    grmice = new TGraph(num, xp, yp);

    f1->SetParLimits(4, xp[4], xp[num-4]);
    f1->SetParameter(4, (xp[num-1]) );
    grmice->Fit(f1,"QR");
    delete grmice;

    std::cout << "==DFA"<< dfa_order <<" F~s with piecewise fit" << std::endl;
    for (int i = 0; i < 5; ++i) {
        pars[i] = f1->GetParameter(i);
        std::cout << "Pars"<< i<< ": " << pars[i] << std::endl;
    }

    auto f2 = new TF1("f2",draw_crossover,newdfa->dfax[0], pow(10,pars[4]),5);
    auto f3 = new TF1("f3",draw_crossover,pow(10,pars[4]),newdfa->dfax[num-1],5);
    for (int i = 0; i < 4; ++i) {
        f2->SetParameter(i, pars[i]);
        f3->SetParameter(i, pars[i]);
    }
    f2->SetParameter(4, pow(10, pars[4]));
    f3->SetParameter(4, pow(10, pars[4]));

    double rvalue, lvalue;
    lvalue = pow(10, pars[4]);
    for (int i = 0; i < num; ++i) {
        if ( lvalue > newdfa->dfax[i] )
            continue;
        else
        {
            lvalue = newdfa->dfax[i-1];
            rvalue = newdfa->dfax[i];
            break;
        }
    }
    /* std::cout << "Crossover-> l: " << lvalue << " c: " << pow(10, pars[4]) << " r: " << rvalue << std::endl; */
    /* std::cout << "Crossover-> l: " <<lvalue*fitStride/24<< "d c: " << pow(10, pars[4])*fitStride/24 << "d r: " << rvalue*fitStride/24 << "d"<< std::endl; */
    char crosslabel[255]; 
    period_dfa.lbp = lvalue*fitStride/24;
    period_dfa.rbp = rvalue*fitStride/24;
    period_dfa.midbp = pow(10.0,pars[4])*fitStride/24;
    period_dfa.index_l = pars[0];
    period_dfa.index_r = pars[2];
    sprintf(crosslabel, "Crossover:%.2fd-%.2fd-%.2fd",lvalue*fitStride/24,pow(10.0,pars[4])*fitStride/24,rvalue*fitStride/24);

    grmice = new TGraph(num, newdfa->dfax, newdfa->dfay);
    grmice->SetName("period_dfa");
    grmice->SetMarkerStyle(8);
    grmice->SetMarkerSize(1);
    grmice->SetTitle(("DFA("+std::to_string(dfa_order)+") of period distribute(" +mice.name+ ")").c_str());
    grmice->GetYaxis()->CenterTitle();
    grmice->GetYaxis()->SetTitle("Fn");
    grmice->GetXaxis()->CenterTitle();
    grmice->GetXaxis()->SetTitle("N");
    canvas->SetLogx(1);
    canvas->SetLogy(1);

    f2->SetLineColor(kRedBlue);
    f3->SetLineColor(kPink);

    grmice->Draw("AP");
    f2->Draw("SAME");
    f3->Draw("SAME");
    /* f1->Draw("SAME"); */


    auto legend = new TLegend(0.1,0.8,0.4,0.9);
    legend->AddEntry("period_dfa","DFA data","l");
    legend->AddEntry("f2",("Exp: " + std::to_string(pars[0])).c_str(),"l");
    legend->AddEntry("f3",("Exp: " + std::to_string(pars[2])).c_str(),"l");
    legend->AddEntry("",crosslabel,"l");
    legend->Draw();

    canvas->Update();
    canvas->Print((mice.name + "_DFA" + std::to_string(dfa_order) + ".pdf").c_str(),"Title:dfa");
    canvas->Clear();
    canvas->SetLogx(0);
    canvas->SetLogy(0);

    delete[] yp;
    delete[] xp;
    delete f1;
    delete f2;
    delete f3;
    delete grmice;
    std::cout << "============================" << std::endl;

    return 0;
}

double fit_crossover(double *x, double *par)
{
    if( x[0] < par[4])
    {
        return par[0] * x[0] + par[1];
    }
    else 
    {
        return par[2] * x[0] + par[3];
    }
}

double draw_crossover(double *x, double *par)
{
    if( x[0] < par[4])
    {
        return pow(10, par[1]) * pow(x[0], par[0]);
    }
    else
    {
        return pow(10, par[3]) * pow(x[0], par[2]);
    }
}

/// @brief A class to solve do the rhythm analysis for our mice.
/// @details A class to solve do the rhythm analysis for our mice.
TMice::TMice(TMicemem _mice, int _winSize)
{   
    mice.name = _mice.name;
    mice.path = _mice.path;
    mice.batch = _mice.batch; 
    mice.ifmutant = _mice.ifmutant;
    winSize = _winSize;
    fitChilimit = 200;
    canvas = new TCanvas("canvas");
}


bool printMice(TMicemem mice)
{
    std::cout << "\nMice details" << std::endl;
    std::cout << "Name:\t" << mice.name << std::endl;
    std::cout << "Path:\t" << mice.path + "/batch" << mice.batch << "/"+mice.name << ".Temperature.txt" << std::endl;
    std::cout << "Path:\t" << mice.path + "/batch" << mice.batch << "/"+mice.name << ".Activity.txt" << std::endl;
    std::cout << "Batch:\tBatch" << mice.batch << std::endl;
    std::cout << "Mutant:\t" << mice.ifmutant << std::endl;
    std::cout << std::endl;
    return 1;
}

int TMice::nullize()
{
    tempOri.data = NULL;
    tempOri.time = NULL;
    temp.data = NULL;
    temp.time = NULL;
    actOri.data = NULL;
    actOri.time = NULL;
    act.data = NULL;
    act.time = NULL;
    dfaRem.data = NULL;
    dfaRem.time = NULL;
    tempRem.data = NULL;
    tempRem.time = NULL;
    actRem.data = NULL;
    actRem.time = NULL;
    period.data = NULL;
    period.time = NULL;

    return 0;
}

TMice::~TMice()
{
    if(tempOri.data != NULL)
    {
        delete[] tempOri.data;
        delete[] tempOri.time;
    }
    if(temp.data != NULL)
    {
        delete[] temp.data;
        delete[] temp.time;
    }
    if(actOri.data != NULL)
    {
        delete[] actOri.data;
        delete[] actOri.time;
    }
    if(act.data != NULL)
    {
        delete[] act.data;
        delete[] act.time;
    }
    if(dfaRem.data != nullptr){
        delete [] dfaRem.data;
        delete [] dfaRem.time;
    }
    if(tempRem.data != nullptr){
        delete [] tempRem.data;
        delete [] tempRem.time;
    }
    if(actRem.data != nullptr){
        delete [] actRem.data;
        delete [] actRem.time;
    }
    if(period.data != nullptr){
        delete [] period.data;
        delete [] period.time;
    }
    if(this->canvas == nullptr)
        delete canvas;
}

/// This function set the window size and number
int TMice::SetwinSize(int _size)
{
    this->winSize = _size;
    this->winNum = act.size / this->winSize;
    return 0;
}

/// This function set the window number and size
int TMice::SetwinNum(int _num)
{
    this->winNum = _num;
    this->winSize = act.size / this->winNum;
    return 0;
}

int TMice::GetwinNum()
{
    return this->winNum;
}

int TMice::GetwinSize()
{
    return this->winSize;
}

double TMice::GetErrorPercentage(int i)
{
    static double percentage;
    if( i == 0 )
        percentage = (double)errors_act/actOri.size;
    else if( i == 1 )
        percentage = (double)errors_temp/tempOri.size;
    return percentage;
}

/// The function to select which part you want to use.
/// @details The input is start and length, in the unit of hours.
int TMice::SetRange(int _start, int _length)
{
    start = _start;
    length = _length;
    ifRangeSeted = 1; 
    return 0;
}

/// @brief The function to read in data according to the TMicemem mice. 
/// @details First, get the number of lines of these data. and assign memory according to this value. Then read in all data into origin array, and from the last value of origin array, we can get the biggest array size we need(Though we not need this big).  Then the real data array.  Also the different batches are operated with different method here.
bool TMice::ReadMice(TMicemem _mice)
{
    mice.name = _mice.name;
    mice.path = _mice.path;
    mice.batch = _mice.batch; 
    mice.ifmutant = _mice.ifmutant;

    errors_act = 0;
    errors_temp = 0;
    
    std::string filepath = mice.path + "/batch" + std::to_string(mice.batch) +"/"+ mice.name + ".Activity.txt";
    int lines = fileLines(filepath.c_str());
    int datamrk = 0, diff_t = 0;
    double tempx = 0, tempy = 0;

    double *data_original = new double[lines];
    double *data_original_t = new double[lines];
    double *temp_data, *temp_data_t;
    int size_data = 0;

    actOri.size = lines;
    if(actOri.data != NULL)
    {
        delete[] actOri.data;
        delete[] actOri.time;
    }
    actOri.data = new double[lines];
    actOri.time = new double[lines];

    FILE* infile = fopen(filepath.c_str(),"r");
    //Check if one can read this file;
    if(infile == NULL)
    {
        fprintf(stderr,"(Read ALL)Failed to read %s\n",filepath.c_str());
        return -1;
    }

    for (int i = 0; i < lines; ++i)
    {
        data_original[i] = 0.0;
        data_original_t[i] = 0.0;
    }

    // read in all Act data into data_original
    for (int i = 1; i < lines; ++i) {
        if(fscanf (infile,"%lf\t%lf\n",&tempx, &tempy) == 2)
        {
            // remove all nan numbers
            if(!std::isnan(tempy))
                data_original[i] = tempy;
            else{
                data_original[i] = data_original[i-1];
                this->errors_act ++;
            }
            // unit of time should be second 
            data_original_t[i] = tempx;
        }
    }
    fclose(infile);
    for (int i = 0; i < lines; ++i) {
        actOri.data[i] = data_original[i];
        actOri.time[i] = data_original_t[i];
    }

    size_data = data_original_t[lines-1]/10;
    temp_data = new double[size_data];
    temp_data_t = new double[size_data];

    for (int i=1;i<lines;i++)
    {
        diff_t = (data_original_t[i] - data_original_t[i-1])/10;
        if( diff_t == 1 )
        {
            temp_data_t[datamrk] = data_original_t[i];
            temp_data[datamrk] = data_original[i];
            datamrk ++;
        }
        else if (diff_t < 100)
        {
            for (int j = 1; j <= diff_t; ++j)
            {
                temp_data[datamrk] = data_original[i-1] + 
                    (data_original[i]-data_original[i-1])/diff_t * j;
                temp_data_t[datamrk] = data_original_t[i-1] + 
                    (data_original_t[i]-data_original_t[i-1])/diff_t * j;
                datamrk ++;
            }
        }
        if(diff_t > 100.0)
        {    
            if(mice.batch == 3 || mice.batch == 5)
            {
                datamrk = 0; 
                /* std::cout << "for Batch3 or 5" << std::endl; */
            }
            if(mice.batch == 6)
            {
                datamrk--; 
                /* std::cout << "for Batch 6" << std::endl; */
                break;
            }
        }
    }
    // datamrk is the reduced size, and then the size_data is integer multiple 2 days.
    size_data = (datamrk/(48*360)) * 48*360;

    if(act.data != nullptr)
    {
        delete[] act.data;
        delete[] act.time;
    }
    act.size = size_data;
    act.data = new double[size_data];
    act.time = new double[size_data];

    for (int i = 0; i < size_data; ++i) {
        act.data[i] = temp_data[i];
        act.time[i] = temp_data_t[i]/3600.0;
    }

    delete[] data_original; delete[] data_original_t;
    delete[] temp_data_t; delete[] temp_data;

// Then start to read in the temperature data;
    filepath = mice.path + "/batch" + std::to_string(mice.batch) +"/"+ mice.name + ".Temperature.txt";
    lines = fileLines(filepath);
    datamrk = 0; diff_t = 0;
    tempx = 0; tempy = 0;

    data_original = new double[lines];
    data_original_t = new double[lines];
    size_data = 0;

    tempOri.size = lines;
    if(tempOri.data != NULL)
    {
        delete[] tempOri.data;
        delete[] tempOri.time;
    }
    tempOri.data = new double[lines];
    tempOri.time = new double[lines];

    infile = fopen(filepath.c_str(),"r");
    //Check if one can read this file;
    if(infile == NULL)
    {
        fprintf(stderr,"(Read ALL)Failed to read %s\n",filepath.c_str());
        return -1;
    }

    for (int i = 0; i < lines; ++i)
    {
        data_original[i] = 0.0;
        data_original_t[i] = 0.0;
    }

    // read in all Temp data into data_original
    for (int i = 1; i < lines; ++i) {
        if(fscanf (infile,"%lf\t%lf\n",&tempx, &tempy) == 2)
        {
            // remove all nan numbers 
            if(!std::isnan(tempy))
                data_original[i] = tempy;
            else {
                data_original[i] = data_original[i-1];
                this->errors_temp ++;
            }
            // unit of time should be second
            data_original_t[i] = tempx;
        }
    }
    fclose(infile);
    for (int i = 0; i < lines; ++i) {
        tempOri.data[i] = data_original[i];
        tempOri.time[i] = data_original_t[i];
    }

    size_data = data_original_t[lines-1]/10;
    temp_data = new double[size_data];
    temp_data_t = new double[size_data];

    for (int i=1;i<lines;i++)
    {
        diff_t = (data_original_t[i] - data_original_t[i-1])/10;
        if( diff_t == 1 )
        {
            temp_data_t[datamrk] = data_original_t[i];
            temp_data[datamrk] = data_original[i];
            datamrk ++;
        }
        else if (diff_t < 100)
        {
            for (int j = 1; j <= diff_t; ++j)
            {
                temp_data[datamrk] = data_original[i-1] + 
                    (data_original[i]-data_original[i-1])/diff_t * j;
                temp_data_t[datamrk] = data_original_t[i-1] + 
                    (data_original_t[i]-data_original_t[i-1])/diff_t * j;
                datamrk ++;
            }
        }
        if(diff_t > 100.0)
        {    
            if(mice.batch == 3 || mice.batch == 5)
            {
                datamrk = 0; 
                /* std::cout << "for Batch3 or 5" << std::endl; */
            }
            if(mice.batch == 6)
            {
                datamrk--; 
                /* std::cout << "for Batch 6" << std::endl; */
                break;
            }
        }
    }
    // datamrk is the reduced size, and then the size_data is integer multiple 2 days.
    size_data = (datamrk/(48*360)) * 48*360;

    if(temp.data != nullptr)
    {
        delete[] temp.data;
        delete[] temp.time;
    }
    temp.size = size_data;
    temp.data = new double[size_data];
    temp.time = new double[size_data];
    for (int i = 0; i < size_data; ++i)
    {
        temp.data[i] = temp_data[i];
        // Turn unit to hour
        temp.time[i] = temp_data_t[i]/3600.0;
    }

    delete[] data_original; delete[] data_original_t;
    delete[] temp_data_t; delete[] temp_data;

    return true;
}

