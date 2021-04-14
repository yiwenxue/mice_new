#include <Rtypes.h>
#include <bits/types/FILE.h>
#include <fstream>
#include <map>
#include <unordered_map>
#include <TCanvas.h>
#include <TGraph.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fftw3.h>
#include <gsl/gsl_statistics_double.h>
#include <iostream>
#include <math.h>
#include <string>
#include <utility.h>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>

int fuckmrzhang()
{
    printf("Fuck Mr.Zhang\n");
    printf("Version: 0.2\n");

    return 0;
}

int root_load()
{
    printf("-----------------------------\n");
    printf("Project: Mice data Analysis\n");
    printf("-----------------------------\n\n");
    printf("ReLoading Include Path...");
    gROOT->ProcessLine(".I src/include/");
    printf("Succeed!\nReLoading Sources...");
    gROOT->ProcessLine(".L src/utility.cpp");
    gROOT->ProcessLine(".L src/mathematics.cpp");
    printf("Succeed!\nLoading Library...");
    gROOT->ProcessLine(".L /usr/lib/libfftw3.so");
    gROOT->ProcessLine(".L /usr/lib/libgslcblas.so");
    gROOT->ProcessLine(".L /usr/lib/libgsl.so");
    gROOT->ProcessLine("usleep(100000)");
    printf("Succeed!\nReInitializing mice...");
    printf("Succees!\n\n");
    return 0;
}

/// @details A function to get number of a file.
int fileLines(std::string path)
{
    FILE *input = fopen(path.c_str(), "r");
    if (input == NULL)
    {
        fprintf(stderr, "(File Lines)Failed to read %s\n", path.c_str());
        return -1;
    }

    //Count the totla lines;
    char chr;
    int lines = 0;
    chr = fgetc(input);
    while (chr != EOF)
    {
        if (chr == '\n')
            lines++;
        chr = fgetc(input);
    }
    // Don't forget the last line.
    lines++;
    fclose(input);
    return lines;
}

/// @details
int TMice::GetDataRange(RangeIndex _index)
{
    if (_index == rSTART)
        return this->start;
    else if (_index == rEND)
        return (this->start + this->length);
    else if (_index == rLENGTH)
        return this->length;
    else
        return -1;
}

int TMice::SetParRange(double *_parMin, double *_parMax)
{
    for (int i = 0; i < 4; i++)
    {
        this->parMax[i] = _parMax[i];
        this->parMin[i] = _parMin[i];
    }
    return 0;
}

double TMice::GetParRange(int i, int j)
{
    double tmp;
    if (i == 0)
        tmp = parMin[j];
    else if (i == 1)
        tmp = parMax[j];
    else
        tmp = NAN;
    return tmp;
}

int TMice::PrintParRange()
{
    for (int i = 0; i < 4; ++i)
    {
        std::cout << "p" << i << "Min = " << this->parMin[i]
                  << "\tp" << i << "Max = " << this->parMax[i]
                  << std::endl;
    }
    return 0;
}

bool TMice::Average()
{
    if (actAve.data == nullptr)
    {
        delete[] actAve.data;
        delete[] actAve.time;
    }
    if (actStd.data == nullptr)
    {
        delete[] actStd.data;
        delete[] actStd.time;
    }
    if (tempAve.data == nullptr)
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

    for (int i = 0; i < winNum; ++i)
    {
        actAve.data[i] = gsl_stats_mean(act.data + i * winSize, 1, winSize);
        actStd.data[i] = gsl_stats_variance(act.data + i * winSize, 1, winSize);
        tempAve.data[i] = gsl_stats_mean(temp.data + i * winSize, 1, winSize);
        actAve.time[i] = act.time[i * winSize];
        actStd.time[i] = actAve.time[i];
        tempAve.time[i] = temp.time[winSize * i];
    }
    return true;
}

int TMice::PrintDetails()
{
    std::cout << "=====Mice details======" << std::endl;
    std::cout << "Name:\t" << mice.name << std::endl;
    std::cout << "Path:\t" << mice.path + "/batch" << mice.batch << "/" + mice.name << ".Temperature.txt" << std::endl;
    std::cout << "Path:\t" << mice.path + "/batch" << mice.batch << "/" + mice.name << ".Activity.txt" << std::endl;
    std::cout << "Batch:\tBatch" << mice.batch << std::endl;
    std::cout << "Mutant:\t" << (mice.ifmutant ? "True" : "False") << std::endl;
    std::cout << "Data Len:\t" << length << "hours|" << length / 24 << "days" << std::endl;
    std::cout << "=====Window details====" << std::endl;
    std::cout << "Window Size:\t" << this->winSize << std::endl;
    std::cout << "Window Num:\t" << this->winNum << std::endl;
    std::cout << "temp size:\t" << temp.size << std::endl;
    std::cout << "tempAve size:\t" << tempAve.size << std::endl;
    std::cout << "act size:\t" << act.size << std::endl;
    std::cout << "actAve size:\t" << actAve.size << std::endl;
    std::cout << "=====Fit Details=======" << std::endl;
    std::cout << "Fit Window Size:\t" << this->fitSize << std::endl;
    std::cout << "Fit Window Stride:\t" << this->fitStride << std::endl;
    PrintParRange();
    std::cout << "=====Data errors=======" << std::endl;
    std::cout << "Act  Errors           : " << this->errors_act << std::endl;
    std::cout << "Act  Total Number     : " << this->actOri.size << std::endl;
    std::cout << "Act  Errors Percentage: " << GetErrorPercentage(0) << std::endl;
    std::cout << "Temp Errors           : " << this->errors_temp << std::endl;
    std::cout << "Temp Total Number     : " << this->tempOri.size << std::endl;
    std::cout << "Temp Errors Percentage: " << GetErrorPercentage(1) << std::endl;
    std::cout << "======================\n"
              << std::endl;

    return 0;
}

/// @brief A function to read in the mice tag db(A index file).
/// @details The db stores the name, batch, path, ifmutant, start, length, and winSize;
int Micedb(std::string infile)
{
    int num = 0;
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
    int itmax = (tempAve.time[tempAve.size - 1] - fitSize - tempAve.time[0]) / fitStride; // how many times will the fit takeplace;
    itmax++;                                                                              // Must plus 1 because from 0 to itmax is itmax +1, if you don't understand, draw it and then count.
    residuals_temp = new double[tempAve.size];

    FitCosinor *result = new FitCosinor[itmax];
    grmice = new TGraph(tempAve.size, tempAve.time, tempAve.data);
    canvas->cd();

    int start_idx, end_idx = 0;

    /// @details Perform the nonlinear fit at all regions
    for (int i = 0; i < itmax; i++)
    {
        fit1 = new TF1("fit1", "[0]*sin([1]*x+[2])+[3]", tempAve.time[0] + i * fitStride, tempAve.time[0] + fitSize + i * fitStride);
        for (int j = 0; j < 4; j++)
        {
            fit1->SetParLimits(j, parMin[j], parMax[j]);
            fit1->SetParameter(j, (parMin[j] + parMax[j]) / 2.0);
        }
        for (int j = 0; j < 10; ++j)
        {
            /// the fit details are disabled, if you want to see the details, please just remove "Q"
            grmice->Fit(fit1, "RQ");
            result[i].chisq = fit1->GetChisquare();
            if (result[i].chisq < fitChilimit)
                break;
        }

        /* for (start_idx = end_idx; start_idx < tempAve.size; ++start_idx) { */
        /*     if (tempAve.time[start_idx] >= tempAve.time[0] + i*fitStride) */
        /*         break; */
        /* } */
        /* for (end_idx = start_idx; end_idx < tempAve.size; ++end_idx) { */
        /*     if (tempAve.time[end_idx+1] >= tempAve.time[0] + fitSize + i*fitStride) */
        /*         break; */
        /* } */
        /* std::cout << "Done !" << std::endl; */
        /* printf("time[%d,%d] = [%f,%f]\n", start_idx, end_idx, tempAve.time[start_idx], */
        /* tempAve.time[end_idx]); */
        /* for (int i = start_idx; i <= end_idx; ++i) { */
        /*     residuals_temp[i] = tempAve.data[i] - fit1->Eval(tempAve.time[i]); */
        /* } */

        result[i].p = new double[4];
        for (int j = 0; j < 4; ++j)
        {
            result[i].p[j] = fit1->GetParameter(j);
        }
        result[i].ndf = fit1->GetNDF();

        if (result[i].chisq > fitChilimit)
        {
            grmice->GetXaxis()->SetRangeUser(tempAve.time[0] + i * fitStride, tempAve.time[0] + fitSize + i * fitStride);
            grmice->Draw("");
            canvas->Update();
            canvas->Print((mice.name + "_badpoint_batch" + std::to_string(mice.batch) + "[" + std::to_string(tempAve.time[0] + i * fitStride) + "," + std::to_string(tempAve.time[0] + fitSize + i * fitStride) + "].pdf").c_str(), "Title:errors");
            canvas->Clear();
            std::cout << "Attention! Bad fit appear at ["
                      << tempAve.time[0] + i * fitStride << "," << tempAve.time[0] + fitSize + i * fitStride
                      << "].\n###################################################" << std::endl;
            /* return -1; */
        }
        delete fit1;
    }
    /* std::cout << "Fit finished" << std::endl; */

    /// Arrays to store the period information temporarily
    if (period.data != nullptr)
    {
        delete[] period.data;
        delete[] period.time;
    }
    period.data = new double[itmax];
    period.time = new double[itmax];
    period.size = itmax;

    for (int i = 0; i < itmax; ++i)
    {
        period.data[i] = 2 * M_PI / result[i].p[1];
        period.time[i] = tempAve.time[0] + fitSize / 2 + i * fitStride;
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

    FILE *period_output = fopen(("period_output/" + mice.name + "_" + std::to_string(mice.batch) + ".txt").c_str(), "w");

    for (int i = 0; i < itmax; i++)
    {
        fprintf(period_output, "%lf\t%lf\n", period.time[i], period.data[i]);
    }

    fclose(period_output);

    for (int i = 0; i < itmax; ++i)
    {
        delete[] result[i].p;
    }
    delete[] result;
    delete grmice;
    return 0;
}

void TMice::calculate_residual()
{
    TF1 *fit1 = nullptr;
    int itmax = (tempAve.time[tempAve.size - 1] - fitSize - tempAve.time[0]) / fitStride; // how many times will the fit takeplace;
    itmax++;                                                                              // Must plus 1 because from 0 to itmax is itmax +1, if you don't understand, draw it and then count.
    residuals_temp = new double[tempAve.size];

    FitCosinor *result = new FitCosinor[itmax];
    grmice = new TGraph(tempAve.size, tempAve.time, tempAve.data);
    canvas->cd();

    int start_idx, end_idx = 0;

    /// @details Perform the nonlinear fit at all regions
    for (int i = 0; i < itmax; i++)
    {
        fit1 = new TF1("fit1", "[0]*sin([1]*x+[2])+[3]", tempAve.time[0] + i * fitStride, tempAve.time[0] + fitSize + i * fitStride);
        for (int j = 0; j < 4; j++)
        {
            fit1->SetParLimits(j, parMin[j], parMax[j]);
            fit1->SetParameter(j, (parMin[j] + parMax[j]) / 2.0);
        }
        for (int j = 0; j < 10; ++j)
        {
            /// the fit details are disabled, if you want to see the details, please just remove "Q"
            grmice->Fit(fit1, "RQ");
            result[i].chisq = fit1->GetChisquare();
            if (result[i].chisq < fitChilimit)
                break;
        }

        for (start_idx = end_idx; start_idx < tempAve.size; ++start_idx)
        {
            if (tempAve.time[start_idx] >= tempAve.time[0] + i * fitStride)
                break;
        }
        for (end_idx = start_idx; end_idx < tempAve.size; ++end_idx)
        {
            if (tempAve.time[end_idx + 1] >= tempAve.time[0] + fitSize + i * fitStride)
                break;
        }
        /* std::cout << "Done !" << std::endl; */
        /* printf("time[%d,%d] = [%f,%f]\n", start_idx, end_idx, tempAve.time[start_idx], */
        /* tempAve.time[end_idx]); */
        for (int i = start_idx; i <= end_idx; ++i)
        {
            residuals_temp[i] = tempAve.data[i] - fit1->Eval(tempAve.time[i]);
        }

        result[i].p = new double[4];
        for (int j = 0; j < 4; ++j)
        {
            result[i].p[j] = fit1->GetParameter(j);
        }
        result[i].ndf = fit1->GetNDF();

        if (result[i].chisq > fitChilimit)
        {
            grmice->GetXaxis()->SetRangeUser(tempAve.time[0] + i * fitStride, tempAve.time[0] + fitSize + i * fitStride);
            grmice->Draw("");
            canvas->Update();
            canvas->Print((mice.name + "_badpoint_batch" + std::to_string(mice.batch) + "[" + std::to_string(tempAve.time[0] + i * fitStride) + "," + std::to_string(tempAve.time[0] + fitSize + i * fitStride) + "].pdf").c_str(), "Title:errors");
            canvas->Clear();
            std::cout << "Attention! Bad fit appear at ["
                      << tempAve.time[0] + i * fitStride << "," << tempAve.time[0] + fitSize + i * fitStride
                      << "].\n###################################################" << std::endl;
            /* return -1; */
        }
        delete fit1;
    }

    delete[] result;
    delete grmice;

    fit1 = nullptr;
    itmax = (actAve.time[actAve.size - 1] - fitSize - actAve.time[0]) / fitStride; // how many times will the fit takeplace;
    itmax++;                                                                       // Must plus 1 because from 0 to itmax is itmax +1, if you don't understand, draw it and then count.
    residuals_act = new double[actAve.size];

    result = new FitCosinor[itmax];
    grmice = new TGraph(actAve.size, actAve.time, actAve.data);
    canvas->cd();

    start_idx = end_idx = 0;
    /// @details Perform the nonlinear fit at all regions
    for (int i = 0; i < itmax; i++)
    {
        fit1 = new TF1("fit1", "[0]*sin([1]*x+[2])+[3]", actAve.time[0] + i * fitStride, actAve.time[0] + fitSize + i * fitStride);
        for (start_idx = end_idx; start_idx < actAve.size; ++start_idx)
        {
            if (actAve.time[start_idx] >= actAve.time[0] + i * fitStride)
                break;
        }
        for (end_idx = start_idx; end_idx < actAve.size; ++end_idx)
        {
            if (actAve.time[end_idx + 1] >= actAve.time[0] + fitSize + i * fitStride)
                break;
        }
        fit1->SetParLimits(0, 0, 10);
        fit1->SetParLimits(1, 0.20, 0.32);
        fit1->SetParLimits(2, -24, 24);
        fit1->SetParLimits(3, gsl_stats_min(actAve.data + start_idx, 1, end_idx - start_idx),
                           gsl_stats_min(actAve.data + start_idx, 1, end_idx - start_idx));
        for (int j = 0; j < 10; ++j)
        {
            /// the fit details are disabled, if you want to see the details, please just remove "Q"
            grmice->Fit(fit1, "RQ");
            result[i].chisq = fit1->GetChisquare();
            if (result[i].chisq < fitChilimit)
                break;
        }
        /* std::cout << "Done !" << std::endl; */
        /* printf("time[%d,%d] = [%f,%f]\n", start_idx, end_idx, actAve.time[start_idx], */
        /* actAve.time[end_idx]); */
        for (int i = start_idx; i <= end_idx; ++i)
        {
            residuals_act[i] = actAve.data[i] - fit1->Eval(actAve.time[i]);
        }

        result[i].p = new double[4];
        for (int j = 0; j < 4; ++j)
        {
            result[i].p[j] = fit1->GetParameter(j);
        }
        result[i].ndf = fit1->GetNDF();

        if (result[i].chisq > fitChilimit)
        {
            grmice->GetXaxis()->SetRangeUser(actAve.time[0] + i * fitStride, actAve.time[0] + fitSize + i * fitStride);
            grmice->Draw("");
            canvas->Update();
            canvas->Print((mice.name + "_badpoint_batch" + std::to_string(mice.batch) + "[" + std::to_string(actAve.time[0] + i * fitStride) + "," + std::to_string(actAve.time[0] + fitSize + i * fitStride) + "].pdf").c_str(), "Title:errors");
            canvas->Clear();
            std::cout << "Attention! Bad fit appear at ["
                      << actAve.time[0] + i * fitStride << "," << actAve.time[0] + fitSize + i * fitStride
                      << "].\n###################################################" << std::endl;
            /* return -1; */
        }
        delete fit1;
    }

    delete[] result;
    delete grmice;
}

void TMice::residuals_cleanup()
{
    delete[] residuals_act;
    delete[] residuals_temp;
}

bool TMice::RhythmRemap()
{

    /* // This will turn the time into phase. */
    /* tempAve.time[0] = 0;// this->winSize/360.0/period[0]*2*M_PI; */
    /* actAve.time[0] = 0;// this->winSize/360.0/period[0]*2*M_PI; */

    std::cout << "=====Remap the data to phase space=====" << std::endl;
    int remapdays, phaseDsize;
    double total_phase = 0;
    for (int i = 0; i < period.size; ++i)
    {
        total_phase += fitStride * 2 * M_PI / period.data[i];
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

    double *phase = new double[tempAve.size];
    phase[0] = 0;
    for (int i = 0; i < period.size; ++i)
    {
        int gap = fitStride * 360 / winSize;
        double dphase = winSize * 2 * M_PI / 360 / period.data[i];
        for (int j = 0; j < gap; ++j)
        {
            phase[i * gap + j] = dphase;
        }
    }
    int itsize = period.size * fitStride * 360 / winSize;
    phase[0] = 0;
    for (int i = 1; i < itsize; ++i)
    {
        phase[i] = phase[i - 1] + phase[i];
    }

    std::cout << "itsize: " << itsize << "\tphase size: " << tempAve.size << std::endl;
    std::cout << "Total phase is: " << phase[itsize - 1] << std::endl;

    /// Because the time in tempave is not eaqual distributed, I need a resample.
    if (tempRem.data != nullptr)
    {
        delete[] tempRem.data;
        delete[] tempRem.time;
    }
    if (actRem.data != nullptr)
    {
        delete[] actRem.data;
        delete[] actRem.time;
    }

    tempRem.size = phaseDsize;
    tempRem.time = new double[tempRem.size];
    tempRem.data = new double[tempRem.size];

    actRem.size = phaseDsize;
    actRem.time = new double[actRem.size];
    actRem.data = new double[actRem.size];

    tempRem.data[0] = tempAve.data[0];
    tempRem.time[0] = 0;
    actRem.data[0] = actAve.data[0];
    actRem.time[0] = 0;

    int actend = 0, tempend = 0;

    for (int i = 1, j = 1; i < actRem.size;)
    {
        while (j * 2 * M_PI / (360.0 / this->winSize) / 24 < phase[i] && j < actRem.size)
        {
            actRem.time[j] = j * 2 * M_PI / (360.0 / this->winSize) / 24;
            actRem.data[j] = actAve.data[i - 1] + (actAve.data[i] - actAve.data[i - 1]) / (phase[i] - phase[i - 1]) * (j * 2 * M_PI / (360.0 / this->winSize) / 24 - phase[i - 1]);
            j++;
        }
        actend = std::min(i, j);
        i++;
    }

    for (int i = 1, j = 1; i < tempRem.size;)
    {
        while (j * 2 * M_PI / (360.0 / this->winSize) / 24 < phase[i] && j < tempRem.size)
        {
            tempRem.time[j] = j * 2 * M_PI / (360.0 / this->winSize) / 24;
            tempRem.data[j] = tempAve.data[i - 1] + (tempAve.data[i] - tempAve.data[i - 1]) / (phase[i] - phase[i - 1]) * (j * 2 * M_PI / (360.0 / this->winSize) / 24 - phase[i - 1]);
            j++;
        }
        tempend = std::min(i, j);
        i++;
    }
    actRem.size = actend;
    tempRem.size = tempend;

    std::cout << "Ends: " << actRem.time[actRem.size - 1] << "\t" << tempRem.time[tempRem.size - 1] << std::endl;

    remapdays = (int)(tempRem.time[tempRem.size - 1] / 2 / M_PI);
    double stopphase = remapdays * 2 * M_PI;
    int stopidx = 0;
    for (int i = 0; i < tempRem.size; ++i)
    {
        if (tempRem.time[i] > stopphase)
        {
            stopidx = i;
            break;
        }
    }
    std::cout << "Stop at: " << stopidx << std::endl;

    std::cout << "Temp[0]: " << tempAve.data[0] << std::endl;

    int output_days = 0, data_len = 0;
    if (remapdays > 30)
        output_days = 30;
    data_len = stopidx / remapdays * output_days;
    FILE *file_output = fopen(("/home/yiwen/mice_physiology/micenew/rebuild/group_average/" + mice.name + "_temperature.txt").c_str(), "w");
    double normalize_mean = gsl_stats_mean(tempRem.data, 1, data_len);
    double normalize_stdev = gsl_stats_sd(tempRem.data, 1, data_len);
    for (int i = 0; i < data_len; ++i)
    {
        fprintf(file_output, "%lf\n", (tempRem.data[i]));
    }
    fclose(file_output);

    file_output = fopen(("/home/yiwen/mice_physiology/micenew/rebuild/group_average/" + mice.name + "_activiry.txt").c_str(), "w");
    normalize_mean = gsl_stats_mean(actRem.data, 1, data_len);
    normalize_stdev = gsl_stats_sd(actRem.data, 1, data_len);
    for (int i = 0; i < data_len; ++i)
    {
        fprintf(file_output, "%lf\n", (actRem.data[i]));
    }
    fclose(file_output);

    int daylen = stopidx / remapdays;
    double *waveform = new double[daylen];
    double *waveformx = new double[daylen];
    std::cout << "Daylen: " << daylen << std::endl;
    for (int i = 0; i < daylen; ++i)
    {
        waveform[i] = 0;
    }
    for (int i = 0; i < remapdays; ++i)
    {
        for (int j = 0; j < daylen; ++j)
        {
            waveform[j] += tempRem.data[i * daylen + j];
        }
    }
    for (int j = 0; j < daylen; ++j)
    {
        waveformx[j] = tempRem.time[j] / 2 / M_PI * 24;
        waveform[j] /= remapdays;
        if (std::isnan(waveform[j]))
        {
            std::cout << "Shitttttttttt" << std::endl;
        }
    }
    for (int i = 0; i < remapdays; ++i)
    {
        for (int j = 0; j < daylen; ++j)
        {
            tempRem.data[i * daylen + j] -= waveform[j];
        }
    }
    tempRem.size = daylen * remapdays;
    std::cout << "Daylen: " << daylen << std::endl;

    file_output = fopen(("/home/yiwen/mice_physiology/micenew/rebuild/group_average/" + mice.name + "_temperature_waveform.txt").c_str(), "w");
    normalize_mean = gsl_stats_mean(waveform, 1, daylen);
    normalize_stdev = gsl_stats_sd(waveform, 1, daylen);
    for (int i = 0; i < daylen; ++i)
    {
        fprintf(file_output, "%lf\n", (waveform[i] - normalize_mean) / normalize_stdev);
    }
    fclose(file_output);

    auto pad = new TPad("pad", "", .0, .0, 1, 1);
    pad->Draw();
    pad->SetGrid();
    pad->SetFillStyle(4000);
    pad->Divide(1, 2, 0.01, 0.01);

    auto tempWaveform = new TGraph(daylen, waveformx, waveform);
    pad->cd(1);
    tempWaveform->SetTitle("Temperature Waveform");
    tempWaveform->GetXaxis()->SetTitle("Time(h)");
    tempWaveform->GetXaxis()->CenterTitle();
    tempWaveform->Draw("AL");
    delete[] waveform;
    delete[] waveformx;

    remapdays = (int)(actRem.time[actRem.size - 1] / 2 / M_PI);
    stopphase = remapdays * 2 * M_PI;
    stopidx = 0;
    for (int i = 0; i < actRem.size; ++i)
    {
        if (actRem.time[i] > stopphase)
        {
            stopidx = i;
            break;
        }
    }
    std::cout << "Stop at: " << stopidx << std::endl;

    daylen = stopidx / remapdays;
    waveform = new double[daylen];
    waveformx = new double[daylen];
    for (int i = 0; i < daylen; ++i)
    {
        waveform[i] = 0;
    }
    for (int i = 0; i < remapdays; ++i)
    {
        for (int j = 0; j < daylen; ++j)
        {
            waveform[j] += actRem.data[i * daylen + j];
        }
    }
    for (int j = 0; j < daylen; ++j)
    {
        waveformx[j] = actRem.time[j] / 2 / M_PI * 24;
        waveform[j] /= remapdays;
        if (std::isnan(waveform[j]))
        {
            std::cout << "Shitttttttttt" << std::endl;
        }
    }
    for (int i = 0; i < remapdays; ++i)
    {
        for (int j = 0; j < daylen; ++j)
        {
            actRem.data[i * daylen + j] -= waveform[j];
        }
    }
    actRem.size = daylen * remapdays;

    file_output = fopen(("/home/yiwen/mice_physiology/micenew/rebuild/group_average/" + mice.name + "_activity_waveform.txt").c_str(), "w");
    normalize_mean = gsl_stats_mean(waveform, 1, daylen);
    normalize_stdev = gsl_stats_sd(waveform, 1, daylen);
    for (int i = 0; i < daylen; ++i)
    {
        fprintf(file_output, "%lf\n", (waveform[i] - normalize_mean) / normalize_stdev);
    }
    fclose(file_output);

    auto actWaveform = new TGraph(daylen, waveformx, waveform);
    pad->cd(2);
    actWaveform->SetTitle("Activity Waveform");
    actWaveform->GetXaxis()->SetTitle("Time(h)");
    actWaveform->GetXaxis()->CenterTitle();
    actWaveform->Draw("AL");
    delete[] waveform;
    delete[] waveformx;

    canvas->Print((mice.name + "_waveform.pdf").c_str());
    delete[] phase;
    delete actWaveform;
    delete tempWaveform;

    std::cout << std::endl;
    std::cout << "=======================================" << std::endl;
    return true;
}

//
int TMice::DrawHeatmap()
{
    auto *hist = new TH2F("hits", ("Heatmap of Temperature(" + mice.name + ")").c_str(), 360 / this->winSize * 48, 0, 360 / this->winSize * 48, this->length / 48 - 1, 0, this->length / 24 - 2);
    for (int i = 0; i < this->tempAve.size; ++i)
    {
        hist->Fill(int(tempAve.time[i]) % (360 / this->winSize * 48), i / (360 / this->winSize * 48) * 2, tempAve.data[i]);
    }
    canvas->Clear();
    hist->SetStats(kFALSE);
    hist->Draw("COLZ");
    canvas->Print((mice.name + "_Temperature.pdf").c_str(), "Title:Heatmap of Temperature");
    delete hist;
    canvas->Clear();
    hist = new TH2F("hits", ("Heatmap of Activity(" + mice.name + ")").c_str(), 360 / this->winSize * 24, 0, 360 / this->winSize * 48, this->length / 48, 0, this->length / 24);
    for (int i = 0; i < this->actAve.size; ++i)
    {
        hist->Fill(int(actAve.time[i]) % (360 / this->winSize * 48), i / (360 / this->winSize * 48) * 2, actAve.data[i]);
    }
    canvas->Clear();
    hist->SetStats(kFALSE);
    hist->Draw("COLZ");
    canvas->Print((mice.name + "_Activity.pdf").c_str(), "Title:Heatmap of Activity");
    delete hist;
    canvas->Clear();
    return 0;
}

int TMice::DrawPeriodDist(void)
{
    int temp_size;
    temp_size = period.size; // the length for all data
    std::cout << "period length: " << temp_size;
    temp_size = length / fitStride - 1; // the length for selected
    std::cout << "Selected length: " << temp_size << std::endl;

    double *avep = new double[temp_size];
    double *avet = new double[temp_size];
    for (int i = 0; i < temp_size; ++i)
    {
        avep[i] = 0;
        avet[i] = 0;
    }

    int j = 0;
    for (int i = 0; i < temp_size; ++i)
    {
        if (i + 1 < fitSize / fitStride)
            j = i + 1;
        else if (i + fitSize / fitStride > temp_size)
            j = temp_size - i;
        else
            j = fitSize / fitStride;
        avep[i] = gsl_stats_mean(period.data + i - j / 2, 1, j);
        avet[i] = gsl_stats_mean(period.time + i - j / 2, 1, j);
        avet[i] /= 24;
    }
    double ma = avet[0];
    for (int i = 0; i < temp_size; ++i)
    {
        avet[i] -= ma;
    }

    TF1 *f1 = new TF1("f1", "24", 0, period.time[temp_size - 1]);

    /* grmice = new TGraph(temp_size, period.time, period.data); */
    grmice = new TGraph(temp_size, avet, avep);
    double ave1 = gsl_stats_mean(avep, 1, temp_size / 2);
    double ave2 = gsl_stats_mean(avep + temp_size / 2, 1, temp_size / 2);
    periodl = ave1;
    periodr = ave2;
    periodall = gsl_stats_mean(avep, 1, temp_size);
    std::string label1("average before half: " + std::to_string(ave1));
    std::string label2("average after half: " + std::to_string(ave2));
    grmice->SetTitle(("Circadian Rhythm(" + mice.name + ")").c_str());
    grmice->GetYaxis()->CenterTitle();
    grmice->GetYaxis()->SetTitle("Rhythm (hours)");
    grmice->GetXaxis()->CenterTitle();
    grmice->GetXaxis()->SetTitle("Time (days)");
    /* canvas->cd(); */
    grmice->Draw("ALP");
    f1->Draw("SAME");

    auto legend = new TLegend(0.1, 0.8, 0.4, 0.9);
    legend->AddEntry("", label1.c_str(), "l");
    legend->AddEntry("", label2.c_str(), "l");
    legend->Draw();

    /* std::cout << "Plot diagram to canvas" << std::endl; */
    /* f1->Draw("SAME"); */
    canvas->Update();
    canvas->Print((mice.name + "_period.pdf").c_str(), "Title:Period");
    canvas->Clear();
    /* std::cout << "Period diagram saved" << std::endl; */

    delete f1;
    delete grmice;
    delete legend;
    delete[] avep;
    delete[] avet;

    return 0;
}

int TMice::DrawOverview()
{
    grmice = new TGraph(tempAve.size, tempAve.time, tempAve.data);
    grmice->SetTitle(("Temperature overvire(" + mice.name + ")").c_str());
    grmice->GetYaxis()->CenterTitle();
    grmice->GetYaxis()->SetTitle("Temperatue(C)");
    grmice->GetXaxis()->CenterTitle();
    grmice->GetXaxis()->SetTitle("Time(hours)");

    grmice->Draw();
    canvas->Update();
    canvas->Print((mice.name + "_temperature_Overview.pdf").c_str(), "Title:temperature");
    canvas->Clear();
    delete grmice;

    grmice = new TGraph(actAve.size, actAve.time, actAve.data);
    grmice->SetTitle(("Activity overvire(" + mice.name + ")").c_str());
    grmice->GetYaxis()->CenterTitle();
    grmice->GetYaxis()->SetTitle("Activity");
    grmice->GetXaxis()->CenterTitle();
    grmice->GetXaxis()->SetTitle("Time(hours)");

    grmice->Draw();
    canvas->Update();
    canvas->Print((mice.name + "_activity_Overview.pdf").c_str(), "Title:activity");
    canvas->Clear();
    delete grmice;
    return 0;
}
//
int DFA_dual_plot(double *data,
                  int data_size,
                  int dfa_order,
                  std::string name,
                  struct period_dfa *input)
{
    std::cout << "==DFA" << dfa_order << " " << name << " F~s with piecewise fit" << std::endl;
    auto newdfa = new DFA(data, data_size, dfa_order);
    int num = newdfa->get_size_o();
    double *xp = new double[num];
    double *yp = new double[num];
    double pars[8];

    for (int i = 0; i < num; ++i)
    {
        xp[i] = newdfa->dfax[i];
        yp[i] = newdfa->dfay[i];
    }

    yp[0] = log10(newdfa->dfay[0]);
    xp[0] = log10(newdfa->dfax[0]);
    for (int i = 1; i < num; ++i)
    {
        yp[i] = log10(newdfa->dfay[i]);
        xp[i] = log10(newdfa->dfax[i]);
        /* std::cout << xp[i] - xp[i-1] << std::endl; */
    }
    /* The preparation of data ends here */

    /* Prepare to fit the dfa points with a special function */
    auto *canvas = new TCanvas("c0");
    auto f1 = new TF1("f1", fit_crossover_dual, xp[0], xp[num - 1], 8);
    auto *plot = new TGraph(num, xp, yp);

    f1->SetParLimits(0, 0, 2);
    f1->SetParameter(0, 1.0);
    f1->SetParLimits(2, 0, 2);
    f1->SetParameter(2, 1);
    f1->SetParLimits(4, 0, 1);
    f1->SetParameter(4, 0.5);

    f1->SetParLimits(6, xp[6], xp[0] * 0.6 + xp[num - 1] * 0.4);
    f1->SetParameter(6, (xp[0] * 0.6 + xp[num - 1] * 0.4));
    f1->SetParLimits(7, xp[4] - xp[0], xp[6] - xp[0]);
    f1->SetParameter(7, xp[6] - xp[0]);

    double chisq = 0;
    for (int i = 0; i < 30; ++i)
    {
        plot->Fit(f1, "R");
        chisq = f1->GetChisquare();
        if (chisq < 0.001)
            break;
    }

    /* std::cout << "Fit done" << std::endl; */

    delete plot;

    for (int i = 0; i < 8; ++i)
    {
        pars[i] = f1->GetParameter(i);
        std::cout << "Pars" << i << ": " << pars[i] << std::endl;
    }

    auto f2 = new TF1("f2", dfa_functions, newdfa->dfax[0], pow(10, pars[6]), 2);
    auto f3 = new TF1("f3", dfa_functions, pow(10, pars[6]), pow(10, pars[6] + pars[7]), 2);
    auto f4 = new TF1("f4", dfa_functions, pow(10, pars[6] + pars[7]), newdfa->dfax[num - 1], 2);

    f2->SetParameter(0, pars[0]);
    f2->SetParameter(1, pars[1]);
    f3->SetParameter(0, pars[2]);
    f3->SetParameter(1, pars[3]);
    f4->SetParameter(0, pars[4]);
    f4->SetParameter(1, pars[5]);

    /* double rvalue, lvalue; */
    /* lvalue = pow(10, pars[4]); */
    /* for (int i = 0; i < num; ++i) { */
    /*     if ( lvalue > newdfa->dfax[i] ) */
    /*         continue; */
    /*     else */
    /*     { */
    /*         lvalue = newdfa->dfax[i-1]; */
    /*         rvalue = newdfa->dfax[i]; */
    /*         break; */
    /*     } */
    /* } */

    /* char crosslabel[255]; */
    /* sprintf(crosslabel, "Crossover:%.2fd-%.2fd-%.2fd",lvalue * 0.010416667,pow(10.0,pars[4]) * 0.010416667,rvalue * 0.010416667); */
    /* input->lbp = lvalue*0.010416667; */
    /* input->rbp = rvalue*0.010416667; */
    /* input->midbp = pow(10.0,pars[4])*0.010416667; */
    /* input->index_l = pars[0]; */
    /* input->index_r = pars[2]; */

    plot = new TGraph(num, newdfa->dfax, newdfa->dfay);
    plot->SetName("period_dfa");
    plot->SetMarkerStyle(8);
    plot->SetMarkerSize(0.6);
    plot->SetTitle(("DFA(" + std::to_string(dfa_order) + ") of " + name).c_str());
    plot->GetYaxis()->CenterTitle();
    plot->GetYaxis()->SetTitle("Fn");
    plot->GetXaxis()->CenterTitle();
    plot->GetXaxis()->SetTitle("N");
    canvas->SetLogx(1);
    canvas->SetLogy(1);

    f2->SetLineColor(kRedBlue);
    f3->SetLineColor(kPink);
    f3->SetLineColor(kGreen);

    plot->Draw("AP");
    f2->Draw("SAME");
    f3->Draw("SAME");
    f4->Draw("SAME");

    /* auto legend = new TLegend(0.1,0.8,0.4,0.9); */
    /* legend->AddEntry("period_dfa","DFA data","l"); */
    /* legend->AddEntry("f2",("Exp: " + std::to_string(pars[0])).c_str(),"l"); */
    /* legend->AddEntry("f3",("Exp: " + std::to_string(pars[2])).c_str(),"l"); */
    /* legend->AddEntry("",crosslabel,"l"); */
    /* legend->Draw(); */

    canvas->Update();
    name[name.find(" ")] = '\t';

    std::cout << name << "\tDFA" << dfa_order << " Crossovers\t" << pow(10.0, pars[6]) / 4.0 << "\t" << pow(10, pars[6] + pars[7]) / 4.0 << "\thour" << std::endl;
    std::cout << name << "\tDFA" << dfa_order << " Alpha\t" << pars[0] << "\t" << pars[2] << "\t" << pars[4] << std::endl;
    name[name.find("\t")] = '_';

    canvas->Print((name + "_DFA" + std::to_string(dfa_order) + ".pdf").c_str(), "Title:sketch");
    canvas->Clear();
    canvas->SetLogx(0);
    canvas->SetLogy(0);

    delete[] yp;
    delete[] xp;
    delete f1;
    delete f2;
    delete f3;
    delete f4;

    delete newdfa;
    delete plot;
    delete canvas;
    std::cout << "===========================================" << std::endl;

    return 0;
}

//
int DFA_plot(double *data,
             int data_size,
             int dfa_order,
             std::string name,
             struct period_dfa *input,
             double crossover)
{
    std::cout << "==DFA" << dfa_order << " " << name << " F~s with piecewise fit" << std::endl;
    auto newdfa = new DFA(data, data_size, dfa_order);
    int num = newdfa->get_size_o();
    double *xp = new double[num];
    double *yp = new double[num];
    double pars[5];

    for (int i = 0; i < num; ++i)
    {
        xp[i] = newdfa->dfax[i];
        yp[i] = newdfa->dfay[i];
    }

    yp[0] = log10(newdfa->dfay[0]);
    xp[0] = log10(newdfa->dfax[0]);
    for (int i = 1; i < num; ++i)
    {
        yp[i] = log10(newdfa->dfay[i]);
        xp[i] = log10(newdfa->dfax[i]);
        /* std::cout << xp[i] - xp[i-1] << std::endl; */
    }

    auto *canvas = new TCanvas("c0");
    auto f1 = new TF1("f1", fit_crossover, xp[0], xp[num - 1], 5);
    auto *plot = new TGraph(num, xp, yp);

    f1->SetParLimits(0, 0, 2);
    f1->SetParameter(0, 1.0);
    f1->SetParLimits(2, 0, 1);
    f1->SetParameter(2, 0.5);
    if (crossover == 0)
    {
        f1->SetParLimits(4, log10(0.5 / 0.010416667), xp[num - 1]);
        f1->SetParameter(4, (xp[0] + xp[num - 1]) / 2.0);
    }
    else
    {
        double day1 = log10(crossover / 0.010416667);
        f1->SetParLimits(4, day1, day1);
        f1->SetParameter(4, day1);
    }
    double chisq = 0;
    for (int i = 0; i < 20; ++i)
    {
        plot->Fit(f1, "Q");
        chisq = f1->GetChisquare();
        if (chisq < 0.01)
            break;
    }
    delete plot;

    for (int i = 0; i < 5; ++i)
    {
        pars[i] = f1->GetParameter(i);
        std::cout << "Pars" << i << ": " << pars[i] << std::endl;
    }

    auto f2 = new TF1("f2", draw_crossover, newdfa->dfax[0], pow(10, pars[4]), 5);
    auto f3 = new TF1("f3", draw_crossover, pow(10, pars[4]), newdfa->dfax[num - 1], 5);
    for (int i = 0; i < 4; ++i)
    {
        f2->SetParameter(i, pars[i]);
        f3->SetParameter(i, pars[i]);
    }
    f2->SetParameter(4, pow(10, pars[4]));
    f3->SetParameter(4, pow(10, pars[4]));

    double rvalue, lvalue;
    lvalue = pow(10, pars[4]);
    for (int i = 0; i < num; ++i)
    {
        if (lvalue > newdfa->dfax[i])
            continue;
        else
        {
            lvalue = newdfa->dfax[i - 1];
            rvalue = newdfa->dfax[i];
            break;
        }
    }

    char crosslabel[255];
    sprintf(crosslabel, "Crossover:%.2fd-%.2fd-%.2fd", lvalue * 0.010416667, pow(10.0, pars[4]) * 0.010416667, rvalue * 0.010416667);
    input->lbp = lvalue * 0.010416667;
    input->rbp = rvalue * 0.010416667;
    input->midbp = pow(10.0, pars[4]) * 0.010416667;
    input->index_l = pars[0];
    input->index_r = pars[2];

    plot = new TGraph(num, newdfa->dfax, newdfa->dfay);
    plot->SetName("period_dfa");
    plot->SetMarkerStyle(8);
    plot->SetMarkerSize(0.6);
    plot->SetTitle(("DFA(" + std::to_string(dfa_order) + ") of " + name).c_str());
    plot->GetYaxis()->CenterTitle();
    plot->GetYaxis()->SetTitle("Fn");
    plot->GetXaxis()->CenterTitle();
    plot->GetXaxis()->SetTitle("N");
    /* double *xx = plot->GetX(); */
    /* double *yy = plot->GetY(); */
    /* for (int i = 0; i < num; ++i) { */
    /*     yy[i] /= xx[i]; */
    /* } */
    canvas->SetLogx(1);
    canvas->SetLogy(1);

    f2->SetLineColor(kRedBlue);
    f3->SetLineColor(kPink);

    plot->Draw("AP");
    /* f2->Draw("SAME"); */
    /* f3->Draw("SAME"); */

    auto legend = new TLegend(0.1, 0.8, 0.4, 0.9);
    legend->AddEntry("period_dfa", "DFA data", "l");
    legend->AddEntry("f2", ("Exp: " + std::to_string(pars[0])).c_str(), "l");
    legend->AddEntry("f3", ("Exp: " + std::to_string(pars[2])).c_str(), "l");
    legend->AddEntry("", crosslabel, "l");
    /* legend->Draw(); */

    canvas->Update();
    name[name.find(" ")] = '_';
    canvas->Print((name + "_DFA" + std::to_string(dfa_order) + ".pdf").c_str(), "Title:sketch");
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
    std::cout << "===========================================" << std::endl;

    return 0;
}

int DFA_filter_plot(double *data,
                    int data_size,
                    int dfa_order,
                    std::string name,
                    struct period_dfa *input,
                    TimeSeq dataRem)
{

    FILE *dfa_storage = nullptr;
    std::string temp_name = name;
    temp_name[temp_name.find(" ")] = '_';

    auto c1 = new TCanvas("c0");
    double scale = 2;
    c1->SetTitle(name.c_str());
    c1->SetCanvasSize(1200 * scale, 700 * scale);
    c1->SetFillStyle(0);
    c1->SetFrameFillStyle(0);
    c1->cd();

    auto pad = new TPad("pad", "", .0, .0, 1, 1);
    /* gPad->RedrawAxis("g"); */
    pad->Draw();
    pad->SetGrid();
    pad->SetFillStyle(4000);
    pad->Divide(4, 3, 0.01, 0.01);

    TGraph *graph1 = nullptr;
    TGraph *graph2 = nullptr;
    TGraph *graph3 = nullptr;
    TGraph *graph4 = nullptr;
    TGraph *graph5 = nullptr;
    TGraph *graph6 = nullptr;

    int complex_size = data_size; //((int) (data_size *1.0 + 0.5))/2 +1;
    fftw_complex *middle = fftw_alloc_complex(complex_size);
    for (int i = 0; i < complex_size; ++i)
    {
        middle[i][0] = 0.0;
        middle[i][1] = 0.0;
    }

    double *x = new double[data_size];
    double mean = gsl_stats_mean(data, 1, data_size);
    for (int i = 0; i < data_size; ++i)
    {
        x[i] = i;
        data[i] -= mean;
    }

    fftw_plan plan = fftw_plan_dft_r2c_1d(complex_size, data, middle, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    double lim = 48;
    double lim2 = 96;
    double limit = log10(lim);
    // DFA plots
    std::cout << "==DFA" << dfa_order << " " << name << " F~s with piecewise fit" << std::endl;
    auto newdfa = new DFA(data, data_size, dfa_order);
    int num = newdfa->get_size_o();
    double *xp = new double[num];
    double *yp = new double[num];
    double pars[5];

    dfa_storage = fopen(("/home/yiwen/mice_physiology/micenew/rebuild/group_average/" + temp_name + "_dfa_data_ori.txt").c_str(), "w");
    for (int i = 0; i < num; ++i)
    {
        xp[i] = newdfa->dfax[i];
        yp[i] = newdfa->dfay[i];
        fprintf(dfa_storage, "%lf\t%lf\n", newdfa->dfax[i], newdfa->dfay[i]);
    }
    fclose(dfa_storage);

    yp[0] = log10(newdfa->dfay[0]);
    xp[0] = log10(newdfa->dfax[0]);
    for (int i = 1; i < num; ++i)
    {
        yp[i] = log10(newdfa->dfay[i]);
        xp[i] = log10(newdfa->dfax[i]);
        /* std::cout << xp[i] - xp[i-1] << std::endl; */
    }
    std::cout << "xp: " << xp[0] << "~" << xp[num - 1] << "\t" << limit << std::endl;
    auto f1 = new TF1("f1", fit_crossover, xp[0], limit, 5);
    auto *plot = new TGraph(num, xp, yp);
    f1->SetParLimits(0, 0, 2);
    f1->SetParameter(0, 1.0);
    /* f1->SetParLimits(2, 0, 1); */
    f1->SetParameter(2, 0.5);
    f1->SetParLimits(4, xp[2], limit);
    f1->SetParameter(4, (xp[2]));

    double chisq = 0;
    for (int i = 0; i < 20; ++i)
    {
        plot->Fit(f1, "R");
        chisq = f1->GetChisquare();
        if (chisq < 0.0001)
            break;
    }
    delete plot;

    for (int i = 0; i < 5; ++i)
    {
        pars[i] = f1->GetParameter(i);
        std::cout << "Pars" << i << ": " << pars[i] << std::endl;
    }

    auto f2 = new TF1("f2", draw_crossover, newdfa->dfax[0] / 4.0, pow(10, pars[4]) / 4.0, 5);
    auto f3 = new TF1("f3", draw_crossover, pow(10, pars[4]) / 4.0, lim / 4.0, 5);
    for (int i = 0; i < 4; ++i)
    {
        f2->SetParameter(i, pars[i]);
        f3->SetParameter(i, pars[i]);
    }
    f2->SetParameter(4, pow(10, pars[4]));
    f3->SetParameter(4, pow(10, pars[4]));
    double rvalue, lvalue;
    lvalue = pow(10, pars[4]);
    for (int i = 0; i < num; ++i)
    {
        if (lvalue > newdfa->dfax[i])
            continue;
        else
        {
            lvalue = newdfa->dfax[i - 1];
            rvalue = newdfa->dfax[i];
            break;
        }
    }

    char crosslabel[255];
    sprintf(crosslabel, "Crossover:%.2fh-%.2fh-%.2fh", lvalue * 0.010416667 * 24, pow(10.0, pars[4]) * 0.010416667 * 24, rvalue * 0.010416667 * 24);
    /* sprintf(crosslabel, "Crossover:%.2fd-%.2fd-%.2fd",lvalue * 0.010416667,pow(10.0,pars[4]) * 0.010416667,rvalue * 0.010416667); */

    std::cout << name << "\tDFA(ori) " << dfa_order << " Crossovers\t" << pow(10.0, pars[4]) * 0.010416667 * 24 << "\thour" << std::endl;
    std::cout << name << "\tDFA(ori) " << dfa_order << " Alpha\t" << pars[0] << "\t" << pars[2] << std::endl;

    auto legend = new TLegend(0.1, 0.8, 0.4, 0.9);
    legend->AddEntry("f2", ("Exp: " + std::to_string(pars[0])).c_str(), "l");
    legend->AddEntry("f3", ("Exp: " + std::to_string(pars[2])).c_str(), "l");
    legend->AddEntry("", crosslabel, "l");

    graph3 = new TGraph(num, newdfa->dfax, newdfa->dfay);
    graph3->SetName("period_dfa");
    graph3->SetMarkerStyle(8);
    graph3->SetMarkerSize(0.6);
    graph3->SetTitle(("DFA(" + std::to_string(dfa_order) + ") of " + name).c_str());
    graph3->GetYaxis()->CenterTitle();
    graph3->GetYaxis()->SetTitle("Fn");
    graph3->GetXaxis()->CenterTitle();
    graph3->GetXaxis()->SetTitle("N(hour)");
    double *dfahour = graph3->GetX();
    for (int i = 0; i < num; ++i)
    {
        dfahour[i] /= 4.0;
    }

    /* pad->cd(3)->SetGrid(); */
    pad->cd(3)->SetLeftMargin(1.1e-1);
    pad->cd(3)->SetLogx(1);
    pad->cd(3)->SetLogy(1);
    graph3->GetXaxis()->SetRangeUser(0, lim2 / 4.0);
    graph3->Draw("AP");
    f2->SetLineWidth(1);
    f3->SetLineWidth(1);

    f2->Draw("ALSAME");
    f3->Draw("ALSAME");
    legend->Draw("SAME");

    /* c1->Update(); */
    /* c1->Print("original_dfa.pdf"); */
    /* pad->cd(3)->SetLogx(0); */
    /* pad->cd(3)->SetLogy(0); */
    delete[] xp;
    delete[] yp;
    delete newdfa;

    // DFA plots
    std::cout << "==DFA_MAG" << dfa_order << " " << name << " F~s with piecewise fit" << std::endl;
    auto newdfamag = new DFA_MAG(data, data_size, dfa_order);
    num = newdfamag->get_size_o();
    xp = new double[num];
    yp = new double[num];

    dfa_storage = fopen(("/home/yiwen/mice_physiology/micenew/rebuild/group_average/" + temp_name + "_dfamag_data_ori.txt").c_str(), "w");

    for (int i = 0; i < num; ++i)
    {
        xp[i] = newdfamag->dfax[i];
        yp[i] = newdfamag->dfay[i];
        fprintf(dfa_storage, "%lf\t%lf\n", newdfamag->dfax[i], newdfamag->dfay[i]);
    }
    fclose(dfa_storage);
    yp[0] = log10(newdfamag->dfay[0]);
    xp[0] = log10(newdfamag->dfax[0]);
    for (int i = 1; i < num; ++i)
    {
        yp[i] = log10(newdfamag->dfay[i]);
        xp[i] = log10(newdfamag->dfax[i]);
        /* std::cout << xp[i] - xp[i-1] << std::endl; */
    }
    auto f10 = new TF1("f10", fit_crossover, xp[0], limit, 5);
    plot = new TGraph(num, xp, yp);
    f10->SetParLimits(0, 0, 2);
    f10->SetParameter(0, 1.0);
    /* f10->SetParLimits(2, 0, 1); */
    f10->SetParameter(2, 0.5);
    f10->SetParLimits(4, xp[2], limit);
    f10->SetParameter(4, (xp[2]));
    chisq = 0;
    for (int i = 0; i < 20; ++i)
    {
        plot->Fit(f10, "R");
        chisq = f10->GetChisquare();
        if (chisq < 0.0001)
            break;
    }
    delete plot;

    for (int i = 0; i < 5; ++i)
    {
        pars[i] = f10->GetParameter(i);
        std::cout << "Pars" << i << ": " << pars[i] << std::endl;
    }

    auto f11 = new TF1("f11", draw_crossover, newdfamag->dfax[0] / 4.0, pow(10, pars[4]) / 4.0, 5);
    auto f12 = new TF1("f12", draw_crossover, pow(10, pars[4]) / 4.0, lim / 4.0, 5);
    for (int i = 0; i < 4; ++i)
    {
        f11->SetParameter(i, pars[i]);
        f12->SetParameter(i, pars[i]);
    }
    f11->SetParameter(4, pow(10, pars[4]));
    f12->SetParameter(4, pow(10, pars[4]));
    lvalue = pow(10, pars[4]);
    for (int i = 0; i < num; ++i)
    {
        if (lvalue > newdfamag->dfax[i])
            continue;
        else
        {
            lvalue = newdfamag->dfax[i - 1];
            rvalue = newdfamag->dfax[i];
            break;
        }
    }

    sprintf(crosslabel, "Crossover:%.2fh-%.2fh-%.2fh", lvalue * 0.010416667 * 24, pow(10.0, pars[4]) * 0.010416667 * 24, rvalue * 0.010416667 * 24);
    /* sprintf(crosslabel, "Crossover:%.2fd-%.2fd-%.2fd",lvalue * 0.010416667,pow(10.0,pars[4]) * 0.010416667,rvalue * 0.010416667); */
    std::cout << name << "\tDFA(ori mag) " << dfa_order << " Crossovers\t" << pow(10.0, pars[4]) * 0.010416667 * 24 << "\thour" << std::endl;
    std::cout << name << "\tDFA(ori mag) " << dfa_order << " Alpha\t" << pars[0] << "\t" << pars[2] << std::endl;
    auto legend4 = new TLegend(0.1, 0.8, 0.4, 0.9);
    legend4->AddEntry("f11", ("Exp: " + std::to_string(pars[0])).c_str(), "l");
    legend4->AddEntry("f12", ("Exp: " + std::to_string(pars[2])).c_str(), "l");
    legend4->AddEntry("", crosslabel, "l");

    pad->cd(4);
    auto graph10 = new TGraph(num, newdfamag->dfax, newdfamag->dfay);
    graph10->SetName("period_dfa");
    graph10->SetMarkerStyle(8);
    graph10->SetMarkerSize(0.6);
    graph10->SetTitle(("DFA_{mag}(" + std::to_string(dfa_order) + ") of " + name).c_str());
    graph10->GetYaxis()->CenterTitle();
    graph10->GetYaxis()->SetTitle("Fn");
    graph10->GetXaxis()->SetTitle("N(hour)");
    graph10->GetXaxis()->CenterTitle();
    dfahour = graph10->GetX();
    for (int i = 0; i < num; ++i)
    {
        dfahour[i] /= 4.0;
    }

    pad->cd(4)->SetLeftMargin(1.1e-1);
    pad->cd(4)->SetLogx(1);
    pad->cd(4)->SetLogy(1);
    graph10->GetXaxis()->SetRangeUser(0, lim2 / 4.0);
    graph10->Draw("AP");
    f11->SetLineWidth(1);
    f12->SetLineWidth(1);

    f11->Draw("ALSAME");
    f12->Draw("ALSAME");
    legend4->Draw("SAME");
    delete newdfamag;

    // The original data;
    graph1 = new TGraph(data_size, x, data);
    pad->cd(1)->SetGrid();
    pad->cd(1)->SetLeftMargin(1.1e-1);
    graph1->SetTitle((name).c_str());
    dfahour = graph1->GetX();
    for (int i = 0; i < data_size; ++i)
    {
        dfahour[i] /= 96.0;
    }
    graph1->GetXaxis()->SetTitle("day");
    graph1->GetXaxis()->CenterTitle();
    graph1->SetLineWidth(1);
    graph1->Draw("AL");
    std::cout << "Size: " << complex_size << std::endl;

    double *powers = new double[complex_size];
    double *powerx = new double[complex_size];

    for (int i = 0; i < complex_size; ++i)
    {
        powerx[i] = complex_size / 4.0 / (i + 1.0);
        powers[i] = middle[i][0] * middle[i][0] + middle[i][1] * middle[i][1];
    }

    double xx[2], yy[2];
    yy[0] = gsl_stats_max(powers, 1, complex_size);
    yy[1] = gsl_stats_min(powers, 1, complex_size);
    xx[0] = gsl_stats_max(powerx, 1, complex_size);
    xx[1] = gsl_stats_min(powerx, 1, complex_size) * 2;
    yy[1] -= (yy[0] - yy[1]) * 0.05;
    TGraph setrange(2, xx, yy);
    setrange.SetTitle(("PowerSpec " + name).c_str());
    setrange.SetMarkerColorAlpha(kBlack, .0);
    setrange.GetXaxis()->CenterTitle();
    setrange.GetXaxis()->SetTitle("hour");

    // Powerspectrum
    graph2 = new TGraph(complex_size, powerx, powers);
    /* pad->cd(2)->SetGrid(); */
    pad->cd(2)->SetLeftMargin(1.1e-1);
    pad->cd(2)->SetLogx();
    graph2->SetLineWidth(1);
    setrange.Draw("AP");
    graph2->Draw("LSAME");

    // Filter stage
    for (int i = 0; i < 10; ++i)
    {
        int maxid = gsl_stats_max_index(powers, 1, complex_size);
        double oldp = powers[maxid];
        powers[maxid] = 0.5 * (powers[maxid - 1] + powers[maxid + 1]);
        double radio = powers[maxid] / oldp;
        radio = sqrt(radio);
        /* middle[maxid][0] *= radio; */
        /* middle[maxid][1] *= radio; */
        middle[maxid][0] *= 0;
        middle[maxid][1] *= 0;
    }

    yy[0] = gsl_stats_max(powers, 1, complex_size);
    yy[1] = gsl_stats_min(powers, 1, complex_size);
    xx[0] = gsl_stats_max(powerx, 1, complex_size);
    xx[1] = gsl_stats_min(powerx, 1, complex_size) * 2;
    yy[1] -= (yy[0] - yy[1]) * 0.05;
    TGraph setrange1(2, xx, yy);
    setrange1.SetTitle(("PowerSpec " + name + " filtered").c_str());
    setrange1.SetMarkerColorAlpha(kBlack, .0);
    setrange1.GetXaxis()->CenterTitle();
    setrange1.GetXaxis()->SetTitle("hour");

    graph5 = new TGraph(complex_size, powerx, powers);
    pad->cd(6)->SetLeftMargin(1.1e-1);
    pad->cd(6)->SetLogx();
    graph5->SetLineWidth(1);
    setrange1.Draw("AP");
    graph5->Draw("LSAME");

    double *result = new double[data_size];
    plan = fftw_plan_dft_c2r_1d(data_size, middle, result, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    double ratio = gsl_stats_max(data, 1, data_size);
    ratio /= gsl_stats_max(result, 1, data_size);
    for (int i = 0; i < data_size; ++i)
    {
        result[i] *= ratio;
    }

    graph4 = new TGraph(data_size, x, result);
    pad->cd(5)->SetGrid();
    pad->cd(5)->SetLeftMargin(1.1e-1);
    graph4->SetTitle((name + " filtered").c_str());
    dfahour = graph4->GetX();
    for (int i = 0; i < data_size; ++i)
    {
        dfahour[i] /= 96.0;
    }
    graph4->GetXaxis()->SetTitle("day");
    graph4->GetXaxis()->CenterTitle();
    graph4->SetLineWidth(1);
    graph4->Draw("AL");

    std::cout << "==DFA" << dfa_order << " " << name << " F~s with piecewise fit" << std::endl;
    newdfa = new DFA(result, data_size, dfa_order);
    num = newdfa->get_size_o();

    delete[] xp;
    delete[] yp;
    xp = new double[num];
    yp = new double[num];

    dfa_storage = fopen(("/home/yiwen/mice_physiology/micenew/rebuild/group_average/" + temp_name + "_dfa_data_detrend.txt").c_str(), "w");
    for (int i = 0; i < num; ++i)
    {
        xp[i] = newdfa->dfax[i];
        yp[i] = newdfa->dfay[i];
        fprintf(dfa_storage, "%lf\t%lf\n", newdfa->dfax[i], newdfa->dfay[i]);
    }
    fclose(dfa_storage);

    yp[0] = log10(newdfa->dfay[0]);
    xp[0] = log10(newdfa->dfax[0]);
    for (int i = 1; i < num; ++i)
    {
        yp[i] = log10(newdfa->dfay[i]);
        xp[i] = log10(newdfa->dfax[i]);
        /* std::cout << xp[i] - xp[i-1] << std::endl; */
    }
    auto f4 = new TF1("f1", fit_crossover, xp[0], limit, 5);
    plot = new TGraph(num, xp, yp);
    f4->SetParLimits(0, 0, 2);
    f4->SetParameter(0, 1.0);
    /* f4->SetParLimits(2, 0, 1); */
    f4->SetParameter(2, 0.5);
    f4->SetParLimits(4, xp[2], limit);
    f4->SetParameter(4, (xp[2]));
    chisq = 0;
    for (int i = 0; i < 20; ++i)
    {
        plot->Fit(f4, "R");
        chisq = f4->GetChisquare();
        if (chisq < 0.0001)
            break;
    }
    delete plot;

    for (int i = 0; i < 5; ++i)
    {
        pars[i] = f4->GetParameter(i);
        std::cout << "Pars" << i << ": " << pars[i] << std::endl;
    }

    auto f5 = new TF1("f2", draw_crossover, newdfa->dfax[0] / 4.0, pow(10, pars[4]) / 4.0, 5);
    auto f6 = new TF1("f3", draw_crossover, pow(10, pars[4]) / 4.0, lim / 4.0, 5);
    for (int i = 0; i < 4; ++i)
    {
        f5->SetParameter(i, pars[i]);
        f6->SetParameter(i, pars[i]);
    }
    f5->SetParameter(4, pow(10, pars[4]));
    f6->SetParameter(4, pow(10, pars[4]));
    rvalue = 0, lvalue = 0;
    lvalue = pow(10, pars[4]);
    for (int i = 0; i < num; ++i)
    {
        if (lvalue > newdfa->dfax[i])
            continue;
        else
        {
            lvalue = newdfa->dfax[i - 1];
            rvalue = newdfa->dfax[i];
            break;
        }
    }

    sprintf(crosslabel, "Crossover:%.2fh-%.2fh-%.2fh", lvalue * 0.010416667 * 24, pow(10.0, pars[4]) * 0.010416667 * 24, rvalue * 0.010416667 * 24);
    std::cout << name << "\tDFA(detrend) " << dfa_order << " Crossovers\t" << pow(10.0, pars[4]) * 0.010416667 * 24 << "\thour" << std::endl;
    std::cout << name << "\tDFA(detrend) " << dfa_order << " Alpha\t" << pars[0] << "\t" << pars[2] << std::endl;
    auto legend2 = new TLegend(0.1, 0.8, 0.4, 0.9);
    legend2->AddEntry("f2", ("Exp: " + std::to_string(pars[0])).c_str(), "l");
    legend2->AddEntry("f3", ("Exp: " + std::to_string(pars[2])).c_str(), "l");
    legend2->AddEntry("", crosslabel, "l");

    graph6 = new TGraph(num, newdfa->dfax, newdfa->dfay);
    graph6->SetName("period_dfa");
    graph6->SetMarkerStyle(8);
    graph6->SetMarkerSize(0.6);
    graph6->SetTitle(("DFA(" + std::to_string(dfa_order) + ") of " + name + " filtered").c_str());
    graph6->GetYaxis()->CenterTitle();
    graph6->GetYaxis()->SetTitle("Fn");
    graph6->GetXaxis()->CenterTitle();
    graph6->GetXaxis()->SetTitle("N(hour)");

    dfahour = graph6->GetX();
    for (int i = 0; i < num; ++i)
    {
        dfahour[i] /= 4.0;
    }
    /* pad->cd(6)->SetGrid(); */
    pad->cd(7)->SetLeftMargin(1.1e-1);
    pad->cd(7)->SetLogx(1);
    pad->cd(7)->SetLogy(1);
    graph6->GetXaxis()->SetRangeUser(0, lim2 / 4.0);
    graph6->Draw("AP");
    f5->SetLineWidth(1);
    f6->SetLineWidth(1);
    f5->Draw("ALSAME");
    f6->Draw("ALSAME");
    legend2->Draw("SAME");
    delete newdfa;

    // DFA plots
    std::cout << "==DFA_MAG" << dfa_order << " " << name << " F~s with piecewise fit" << std::endl;
    newdfamag = new DFA_MAG(result, data_size, dfa_order);
    num = newdfamag->get_size_o();
    xp = new double[num];
    yp = new double[num];

    dfa_storage = fopen(("/home/yiwen/mice_physiology/micenew/rebuild/group_average/" + temp_name + "_dfamag_data_detrend.txt").c_str(), "w");
    for (int i = 0; i < num; ++i)
    {
        xp[i] = newdfamag->dfax[i];
        yp[i] = newdfamag->dfay[i];
        fprintf(dfa_storage, "%lf\t%lf\n", newdfamag->dfax[i], newdfamag->dfay[i]);
    }
    fclose(dfa_storage);

    yp[0] = log10(newdfamag->dfay[0]);
    xp[0] = log10(newdfamag->dfax[0]);
    for (int i = 1; i < num; ++i)
    {
        yp[i] = log10(newdfamag->dfay[i]);
        xp[i] = log10(newdfamag->dfax[i]);
        /* std::cout << xp[i] - xp[i-1] << std::endl; */
    }
    auto f13 = new TF1("f12", fit_crossover, xp[0], limit, 5);
    plot = new TGraph(num, xp, yp);
    f13->SetParLimits(0, 0, 2);
    f13->SetParameter(0, 1.0);
    /* f13->SetParLimits(2, 0, 1); */
    f13->SetParameter(2, 0.5);
    f13->SetParLimits(4, xp[2], limit);
    f13->SetParameter(4, (xp[2]));
    chisq = 0;
    for (int i = 0; i < 20; ++i)
    {
        plot->Fit(f13, "R");
        chisq = f13->GetChisquare();
        if (chisq < 0.0001)
            break;
    }
    delete plot;

    for (int i = 0; i < 5; ++i)
    {
        pars[i] = f13->GetParameter(i);
        std::cout << "Pars" << i << ": " << pars[i] << std::endl;
    }

    auto f14 = new TF1("f14", draw_crossover, newdfamag->dfax[0] / 4.0, pow(10, pars[4]) / 4.0, 5);
    auto f15 = new TF1("f15", draw_crossover, pow(10, pars[4]) / 4.0, lim / 4.0, 5);
    for (int i = 0; i < 4; ++i)
    {
        f14->SetParameter(i, pars[i]);
        f15->SetParameter(i, pars[i]);
    }
    f14->SetParameter(4, pow(10, pars[4]));
    f15->SetParameter(4, pow(10, pars[4]));
    lvalue = pow(10, pars[4]);
    for (int i = 0; i < num; ++i)
    {
        if (lvalue > newdfamag->dfax[i])
            continue;
        else
        {
            lvalue = newdfamag->dfax[i - 1];
            rvalue = newdfamag->dfax[i];
            break;
        }
    }

    sprintf(crosslabel, "Crossover:%.2fh-%.2fh-%.2fh", lvalue * 0.010416667 * 24, pow(10.0, pars[4]) * 0.010416667 * 24, rvalue * 0.010416667 * 24);
    std::cout << name << "\tDFA(detrend mag) " << dfa_order << " Crossovers\t" << pow(10.0, pars[4]) * 0.010416667 * 24 << "\thour" << std::endl;
    std::cout << name << "\tDFA(detrend mag) " << dfa_order << " Alpha\t" << pars[0] << "\t" << pars[2] << std::endl;
    /* sprintf(crosslabel, "Crossover:%.2fd-%.2fd-%.2fd",lvalue * 0.010416667,pow(10.0,pars[4]) * 0.010416667,rvalue * 0.010416667); */
    auto legend5 = new TLegend(0.1, 0.8, 0.4, 0.9);
    legend5->AddEntry("f14", ("Exp: " + std::to_string(pars[0])).c_str(), "l");
    legend5->AddEntry("f15", ("Exp: " + std::to_string(pars[2])).c_str(), "l");
    legend5->AddEntry("", crosslabel, "l");

    pad->cd(8);
    auto graph11 = new TGraph(num, newdfamag->dfax, newdfamag->dfay);
    graph11->SetName("period_dfa");
    graph11->SetMarkerStyle(8);
    graph11->SetMarkerSize(0.6);
    graph11->SetTitle(("DFA_{mag}(" + std::to_string(dfa_order) + ") of " + name + " filtered").c_str());
    graph11->GetYaxis()->CenterTitle();
    graph11->GetYaxis()->SetTitle("Fn");
    graph11->GetXaxis()->CenterTitle();
    graph11->GetXaxis()->SetTitle("N(hour)");
    dfahour = graph11->GetX();
    for (int i = 0; i < num; ++i)
    {
        dfahour[i] /= 4.0;
    }

    pad->cd(8)->SetLeftMargin(1.1e-1);
    pad->cd(8)->SetLogx(1);
    pad->cd(8)->SetLogy(1);
    graph11->GetXaxis()->SetRangeUser(0, lim2 / 4.0);
    graph11->Draw("AP");
    f14->SetLineWidth(1);
    f15->SetLineWidth(1);

    f14->Draw("ALSAME");
    f15->Draw("ALSAME");
    legend5->Draw("SAME");
    delete newdfamag;

    auto graph7 = new TGraph(dataRem.size, dataRem.time, dataRem.data);
    /* for (int i = 0; i < dataRem.size; ++i) { */
    /*     std::cout << dataRem.time[i] << "\t" << dataRem.data[i] << std::endl; */
    /* } */
    pad->cd(9)->SetGrid();
    pad->cd(9)->SetLeftMargin(1.1e-1);
    graph7->SetTitle((name + " no waveform").c_str());
    dfahour = graph7->GetX();
    for (int i = 0; i < dataRem.size; ++i)
    {
        dfahour[i] /= 2 * M_PI;
    }
    graph7->GetXaxis()->SetTitle("day");
    graph7->GetXaxis()->CenterTitle();
    graph7->SetLineWidth(1);
    graph7->Draw("AL");

    complex_size = dataRem.size; //((int) (data_size *1.0 + 0.5))/2 +1;
    fftw_free(middle);
    middle = fftw_alloc_complex(complex_size);
    for (int i = 0; i < complex_size; ++i)
    {
        middle[i][0] = 0.0;
        middle[i][1] = 0.0;
    }

    delete[] x;
    x = new double[dataRem.size];
    mean = gsl_stats_mean(dataRem.data, 1, dataRem.size);

    plan = fftw_plan_dft_r2c_1d(complex_size, dataRem.data, middle, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    delete[] powerx;
    delete[] powers;
    powers = new double[complex_size];
    powerx = new double[complex_size];

    for (int i = 0; i < complex_size; ++i)
    {
        powerx[i] = complex_size / 4.0 / (i + 1.0);
        powers[i] = middle[i][0] * middle[i][0] + middle[i][1] * middle[i][1];
    }

    yy[0] = gsl_stats_max(powers, 1, complex_size);
    yy[1] = gsl_stats_min(powers, 1, complex_size);
    xx[0] = gsl_stats_max(powerx, 1, complex_size);
    xx[1] = gsl_stats_min(powerx, 1, complex_size) * 2;
    yy[1] -= (yy[0] - yy[1]) * 0.05;
    TGraph setrange2(2, xx, yy);
    setrange2.SetTitle(("PowerSpec " + name + " no waveform").c_str());
    setrange2.SetMarkerColorAlpha(kBlack, .0);
    setrange2.GetXaxis()->CenterTitle();
    setrange2.GetXaxis()->SetTitle("hour");

    auto graph8 = new TGraph(dataRem.size, powerx, powers);
    pad->cd(10)->SetLeftMargin(1.1e-1);
    pad->cd(10)->SetLogx();
    graph8->SetLineStyle(1);
    graph8->SetLineWidth(1);
    setrange2.Draw("AP");
    graph8->Draw("LSAME");

    std::cout << "==DFA" << dfa_order << " " << name << " F~s with piecewise fit" << std::endl;
    newdfa = new DFA(dataRem.data, dataRem.size, dfa_order);
    num = newdfa->get_size_o();

    std::cout << "waveform size: " << num << std::endl;
    delete[] xp;
    delete[] yp;
    xp = new double[num];
    yp = new double[num];

    dfa_storage = fopen(("/home/yiwen/mice_physiology/micenew/rebuild/group_average/" + temp_name + "_dfa_data_nowaveform.txt").c_str(), "w");
    for (int i = 0; i < num; ++i)
    {
        xp[i] = newdfa->dfax[i];
        yp[i] = newdfa->dfay[i];
        fprintf(dfa_storage, "%lf\t%lf\n", newdfa->dfax[i], newdfa->dfay[i]);
    }
    fclose(dfa_storage);

    yp[0] = log10(newdfa->dfay[0]);
    xp[0] = log10(newdfa->dfax[0]);
    for (int i = 1; i < num; ++i)
    {
        yp[i] = log10(newdfa->dfay[i]);
        xp[i] = log10(newdfa->dfax[i]);
        /* std::cout << xp[i] - xp[i-1] << std::endl; */
    }
    auto f7 = new TF1("f1", fit_crossover, xp[0], limit, 5);
    plot = new TGraph(num, xp, yp);
    f7->SetParLimits(0, 0, 2);
    f7->SetParameter(0, 1.0);
    /* f7->SetParLimits(2, 0, 1); */
    f7->SetParameter(2, 0.5);
    f7->SetParLimits(4, xp[2], limit);
    f7->SetParameter(4, (xp[2]));
    chisq = 0;
    for (int i = 0; i < 20; ++i)
    {
        plot->Fit(f7, "R");
        chisq = f7->GetChisquare();
        if (chisq < 0.0001)
            break;
    }
    delete plot;

    for (int i = 0; i < 5; ++i)
    {
        pars[i] = f7->GetParameter(i);
        std::cout << "Pars" << i << ": " << pars[i] << std::endl;
    }

    auto f8 = new TF1("f8", draw_crossover, newdfa->dfax[0] / 4.0, pow(10, pars[4]) / 4.0, 5);
    auto f9 = new TF1("f9", draw_crossover, pow(10, pars[4]) / 4.0, lim / 4.0, 5);
    for (int i = 0; i < 4; ++i)
    {
        f8->SetParameter(i, pars[i]);
        f9->SetParameter(i, pars[i]);
    }
    f8->SetParameter(4, pow(10, pars[4]));
    f9->SetParameter(4, pow(10, pars[4]));
    rvalue = 0, lvalue = 0;
    lvalue = pow(10, pars[4]);
    for (int i = 0; i < num; ++i)
    {
        if (lvalue > newdfa->dfax[i])
            continue;
        else
        {
            lvalue = newdfa->dfax[i - 1];
            rvalue = newdfa->dfax[i];
            break;
        }
    }

    sprintf(crosslabel, "Crossover:%.2fh-%.2fh-%.2fh", lvalue * 0.010416667 * 24, pow(10.0, pars[4]) * 0.010416667 * 24, rvalue * 0.010416667 * 24);
    std::cout << name << "\tDFA(no waveform) " << dfa_order << " Crossovers\t" << pow(10.0, pars[4]) * 0.010416667 * 24 << "\thour" << std::endl;
    std::cout << name << "\tDFA(no waveform) " << dfa_order << " Alpha\t" << pars[0] << "\t" << pars[2] << std::endl;
    auto legend3 = new TLegend(0.1, 0.8, 0.4, 0.9);
    legend3->AddEntry("f8", ("Exp: " + std::to_string(pars[0])).c_str(), "l");
    legend3->AddEntry("f9", ("Exp: " + std::to_string(pars[2])).c_str(), "l");
    legend3->AddEntry("", crosslabel, "l");

    pad->cd(11);
    auto graph9 = new TGraph(num, newdfa->dfax, newdfa->dfay);
    std::cout << "What the hell? num: " << num << std::endl;
    graph9->SetMarkerStyle(8);
    graph9->SetMarkerSize(0.6);
    graph9->SetTitle(("DFA(" + std::to_string(dfa_order) + ") of " + name + " no waveform").c_str());
    graph9->GetYaxis()->CenterTitle();
    graph9->GetYaxis()->SetTitle("Fn");
    graph9->GetXaxis()->CenterTitle();
    graph9->GetXaxis()->SetTitle("N");

    dfahour = graph9->GetX();
    for (int i = 0; i < num; ++i)
    {
        dfahour[i] /= 4.0;
        /* std::cout << "dfa hour: " << dfahour[i] << std::endl; */
    }
    pad->cd(11)->SetLeftMargin(1.1e-1);
    pad->cd(11)->SetLogx(1);
    pad->cd(11)->SetLogy(1);
    graph9->GetXaxis()->SetRangeUser(0, lim2 / 4.0);
    graph9->Draw("AP");
    f8->SetLineWidth(1);
    f9->SetLineWidth(1);
    f8->Draw("ALSAME");
    f9->Draw("ALSAME");
    legend3->Draw("SAME");
    delete newdfa;

    // DFA plots
    std::cout << "==DFA_MAG" << dfa_order << " " << name << " F~s with piecewise fit" << std::endl;
    newdfamag = new DFA_MAG(dataRem.data, dataRem.size, dfa_order);
    num = newdfamag->get_size_o();
    std::cout << "waveform MAG size: " << num << std::endl;
    xp = new double[num];
    yp = new double[num];

    dfa_storage = fopen(("/home/yiwen/mice_physiology/micenew/rebuild/group_average/" + temp_name + "_dfamag_data_nowaveform.txt").c_str(), "w");
    for (int i = 0; i < num; ++i)
    {
        xp[i] = newdfamag->dfax[i];
        yp[i] = newdfamag->dfay[i];
        fprintf(dfa_storage, "%lf\t%lf\n", newdfamag->dfax[i], newdfamag->dfay[i]);
    }
    fclose(dfa_storage);

    yp[0] = log10(newdfamag->dfay[0]);
    xp[0] = log10(newdfamag->dfax[0]);
    for (int i = 1; i < num; ++i)
    {
        yp[i] = log10(newdfamag->dfay[i]);
        xp[i] = log10(newdfamag->dfax[i]);
        /* std::cout << xp[i] - xp[i-1] << std::endl; */
    }
    auto f16 = new TF1("f16", fit_crossover, xp[0], limit, 5);
    plot = new TGraph(num, xp, yp);
    f16->SetParLimits(0, 0, 2);
    f16->SetParameter(0, 1.0);
    /* f16->SetParLimits(2, 0, 1); */
    f16->SetParameter(2, 0.5);
    f16->SetParLimits(4, xp[2], limit);
    f16->SetParameter(4, (xp[2]));
    chisq = 0;
    for (int i = 0; i < 20; ++i)
    {
        plot->Fit(f16, "R");
        chisq = f16->GetChisquare();
        if (chisq < 0.0001)
            break;
    }
    delete plot;

    for (int i = 0; i < 5; ++i)
    {
        pars[i] = f16->GetParameter(i);
        std::cout << "Pars" << i << ": " << pars[i] << std::endl;
    }

    auto f17 = new TF1("f14", draw_crossover, newdfamag->dfax[0] / 4.0, pow(10, pars[4]) / 4.0, 5);
    auto f18 = new TF1("f15", draw_crossover, pow(10, pars[4]) / 4.0, lim / 4.0, 5);
    for (int i = 0; i < 4; ++i)
    {
        f17->SetParameter(i, pars[i]);
        f18->SetParameter(i, pars[i]);
    }
    f17->SetParameter(4, pow(10, pars[4]));
    f18->SetParameter(4, pow(10, pars[4]));
    lvalue = pow(10, pars[4]);
    for (int i = 0; i < num; ++i)
    {
        if (lvalue > newdfamag->dfax[i])
            continue;
        else
        {
            lvalue = newdfamag->dfax[i - 1];
            rvalue = newdfamag->dfax[i];
            break;
        }
    }

    sprintf(crosslabel, "Crossover:%.2fh-%.2fh-%.2fh", lvalue * 0.010416667 * 24, pow(10.0, pars[4]) * 0.010416667 * 24, rvalue * 0.010416667 * 24);
    std::cout << name << "\tDFA(no waveform|mag) " << dfa_order << " Crossovers\t" << pow(10.0, pars[4]) * 0.010416667 * 24 << "\thour" << std::endl;
    std::cout << name << "\tDFA(no waveform|mag) " << dfa_order << " Alpha\t" << pars[0] << "\t" << pars[2] << std::endl;
    /* sprintf(crosslabel, "Crossover:%.2fd-%.2fd-%.2fd",lvalue * 0.010416667,pow(10.0,pars[4]) * 0.010416667,rvalue * 0.010416667); */
    auto legend6 = new TLegend(0.1, 0.8, 0.4, 0.9);
    legend6->AddEntry("f14", ("Exp: " + std::to_string(pars[0])).c_str(), "l");
    legend6->AddEntry("f15", ("Exp: " + std::to_string(pars[2])).c_str(), "l");
    legend6->AddEntry("", crosslabel, "l");

    pad->cd(12);
    auto graph12 = new TGraph(num, newdfamag->dfax, newdfamag->dfay);
    graph12->SetName("period_dfa");
    graph12->SetMarkerStyle(8);
    graph12->SetMarkerSize(0.6);
    graph12->SetTitle(("DFA_{mag}(" + std::to_string(dfa_order) + ") of " + name + " no waveform").c_str());
    graph12->GetYaxis()->CenterTitle();
    graph12->GetYaxis()->SetTitle("Fn");
    graph12->GetXaxis()->CenterTitle();
    graph12->GetXaxis()->SetTitle("N(hour)");
    dfahour = graph12->GetX();
    for (int i = 0; i < num; ++i)
    {
        dfahour[i] /= 4.0;
    }

    pad->cd(12)->SetLeftMargin(1.1e-1);
    pad->cd(12)->SetLogx(1);
    pad->cd(12)->SetLogy(1);
    graph12->GetXaxis()->SetRangeUser(0, lim2 / 4.0);
    graph12->Draw("AP");
    f17->SetLineWidth(1);
    f18->SetLineWidth(1);

    f17->Draw("ALSAME");
    f18->Draw("ALSAME");
    legend6->Draw("SAME");
    delete newdfamag;

    c1->Update();
    name[name.find(" ")] = '_';
    c1->Print((name + "_DFA" + std::to_string(dfa_order) + "_detrend").c_str(), "Title:detrend dfa");

    delete graph1;
    delete graph3;
    delete graph5;
    delete graph2;
    delete graph4;
    delete graph6;
    delete graph7;
    delete graph8;
    delete graph9;
    delete graph10;
    delete graph11;
    delete graph12;
    delete f1;
    delete f2;
    delete f3;
    delete f4;
    delete f5;
    delete f6;
    delete f7;
    delete f8;
    delete f9;
    delete f10;
    delete f11;
    delete f12;
    delete f13;
    delete f14;
    delete f16;
    delete f17;
    delete f18;
    delete legend;
    delete legend2;
    delete legend3;
    delete legend4;
    delete legend5;
    delete legend6;

    delete[] powers;
    delete[] powerx;
    delete[] result;
    delete[] x;

    fftw_free(middle);
    delete c1;

    return 0;
}

int TMice::PeriodDFA(int dfa_order)
{
    std::cout << "=====DFA" << dfa_order << " for periods=======" << std::endl;
    start_index = start;
    auto newdfa = new DFAI(period.data + start_index, (length - fitStride) / fitStride, dfa_order);

    int num = newdfa->get_size_o();
    double *yp = new double[num];
    double *xp = new double[num];
    double pars[5];

    for (int i = 0; i < num; ++i)
    {
        yp[i] = log10(newdfa->dfay[i]);
        xp[i] = log10(newdfa->dfax[i] * fitStride);
    }

    auto dfa_storage = fopen(("/home/yiwen/mice_physiology/micenew/rebuild/group_average/" + mice.name + "_dfa_circadian_data_ori.txt").c_str(), "w");
    for (int i = 0; i < num; ++i)
    {
        // xp[i] = newdfa->dfax[i];
        // yp[i] = newdfa->dfay[i];
        fprintf(dfa_storage, "%lf\t%lf\n", newdfa->dfax[i], newdfa->dfay[i]);
    }
    fclose(dfa_storage);

    std::cout << "Start fit" << std::endl;

    auto f1 = new TF1("f1", fit_crossover, xp[0], xp[num - 1], 5);
    grmice = new TGraph(num, xp, yp);

    f1->SetParLimits(0, 0, 2);
    f1->SetParameter(0, 1.0);
    /* f1->SetParLimits(2, 0, 2); */
    f1->SetParameter(2, 0.5);
    f1->SetParLimits(4, xp[3], xp[num - 2]);
    f1->SetParameter(4, (xp[0] + xp[num - 2]) / 2.0);

    double chisq = 0;
    for (int i = 0; i < 20; ++i)
    {
        grmice->Fit(f1, "R");
        chisq = f1->GetChisquare();
        if (chisq < 0.001)
            break;
    }
    delete grmice;

    std::cout << "==DFA" << dfa_order << " F~s with piecewise fit" << std::endl;
    for (int i = 0; i < 5; ++i)
    {
        pars[i] = f1->GetParameter(i);
        std::cout << "Pars" << i << ": " << pars[i] << std::endl;
    }

    auto f2 = new TF1("f2", draw_crossover, newdfa->dfax[0] * fitStride, pow(10, pars[4]), 5);
    auto f3 = new TF1("f3", draw_crossover, pow(10, pars[4]), newdfa->dfax[num - 1] * fitStride, 5);
    std::cout << xp[0] << "<<" << pars[4] << "<<" << xp[num - 1] << std::endl;
    for (int i = 0; i < 4; ++i)
    {
        f2->SetParameter(i, pars[i]);
        f3->SetParameter(i, pars[i]);
    }
    f2->SetParameter(4, pow(10, pars[4]));
    f3->SetParameter(4, pow(10, pars[4]));

    double rvalue, lvalue;
    lvalue = pow(10, pars[4]);
    for (int i = 0; i < num; ++i)
    {
        if (lvalue > newdfa->dfax[i])
            continue;
        else
        {
            lvalue = newdfa->dfax[i - 1];
            rvalue = newdfa->dfax[i];
            break;
        }
    }

    std::cout << "Fit Finished" << std::endl;

    std::cout << "Crossover-> l: " << lvalue << " c: " << pow(10, pars[4]) << " r: " << rvalue << std::endl;
    /* std::cout << "Crossover-> l: " <<lvalue*fitStride/24<< "d c: " << pow(10, pars[4])*fitStride/24 << "d r: " << rvalue*fitStride/24 << "d"<< std::endl; */
    char crosslabel[255];
    period_dfa.lbp = lvalue / 24;
    period_dfa.rbp = rvalue / 24;
    period_dfa.midbp = pow(10.0, pars[4]) * fitStride / 24;
    period_dfa.index_l = pars[0];
    period_dfa.index_r = pars[2];
    sprintf(crosslabel, "Crossover:%.2fd-%.2fd-%.2fd", lvalue / 24, pow(10.0, pars[4]) / 24, rvalue / 24);

    for (int i = 0; i < num; i++)
    {
        xp[i] = newdfa->dfax[i] * fitStride;
    }

    std::cout << "start to plot" << std::endl;

    grmice = new TGraph(num, xp, newdfa->dfay);
    grmice->SetName("cirdadian_rhythm_dfa");
    grmice->SetMarkerStyle(8);
    grmice->SetMarkerSize(1);
    grmice->SetTitle(("DFA(" + std::to_string(dfa_order) + ") of Circadian Rhythm(" + mice.name + ")").c_str());
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

    auto legend = new TLegend(0.1, 0.8, 0.4, 0.9);
    legend->AddEntry("period_dfa", "DFA data", "l");
    legend->AddEntry("f2", ("Exp: " + std::to_string(pars[0])).c_str(), "l");
    legend->AddEntry("f3", ("Exp: " + std::to_string(pars[2])).c_str(), "l");
    legend->AddEntry("", crosslabel, "l");
    legend->Draw();

    std::cout << "Plot done" << std::endl;
    canvas->Update();
    canvas->Print((mice.name + "_DFA" + std::to_string(dfa_order) + ".pdf").c_str(), "Title:dfa");
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
    if (x[0] < par[4])
    {
        return par[0] * x[0] + par[1];
    }
    else
    {
        return par[2] * x[0] + par[3];
    }
}

double fit_crossover_dual(double *x, double *par)
{
    if (x[0] < par[6])
    {
        return par[0] * x[0] + par[1];
    }
    else
    {
        if (x[0] < (par[6] + par[7]))
            return par[2] * x[0] + par[3];
        else
            return par[4] * x[0] + par[5];
    }
}

// this function need 2 parameter like y = a x + b
double dfa_functions(double *x, double *par)
{
    return pow(10, par[1]) * pow(x[0], par[0]);
}

double draw_crossover(double *x, double *par)
{
    if (x[0] * 4.0 < par[4])
    {
        return pow(10, par[1]) * pow(x[0] * 4.0, par[0]);
    }
    else
    {
        return pow(10, par[3]) * pow(x[0] * 4.0, par[2]);
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
    tempAve.data = NULL;
    tempAve.time = NULL;
    actAve.data = NULL;
    actAve.time = NULL;
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
    if (tempOri.data != NULL)
    {
        delete[] tempOri.data;
        delete[] tempOri.time;
    }
    if (temp.data != NULL)
    {
        delete[] temp.data;
        delete[] temp.time;
    }
    if (actOri.data != NULL)
    {
        delete[] actOri.data;
        delete[] actOri.time;
    }
    if (act.data != NULL)
    {
        delete[] act.data;
        delete[] act.time;
    }
    if (dfaRem.data != nullptr)
    {
        delete[] dfaRem.data;
        delete[] dfaRem.time;
    }
    if (tempRem.data != nullptr)
    {
        delete[] tempRem.data;
        delete[] tempRem.time;
    }
    if (actRem.data != nullptr)
    {
        delete[] actRem.data;
        delete[] actRem.time;
    }
    if (period.data != nullptr)
    {
        delete[] period.data;
        delete[] period.time;
    }
    if (this->canvas == nullptr)
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
    if (i == 0)
        percentage = (double)errors_act / actOri.size;
    else if (i == 1)
        percentage = (double)errors_temp / tempOri.size;
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

    std::string filepath = mice.path + "/batch" + std::to_string(mice.batch) + "/" + mice.name + ".Activity.txt";
    int lines = fileLines(filepath.c_str());
    int datamrk = 0, diff_t = 0;
    double tempx = 0, tempy = 0;

    double *data_original = new double[lines];
    double *data_original_t = new double[lines];
    double *temp_data, *temp_data_t;
    int size_data = 0;

    actOri.size = lines;
    if (actOri.data != NULL)
    {
        delete[] actOri.data;
        delete[] actOri.time;
    }
    actOri.data = new double[lines];
    actOri.time = new double[lines];

    FILE *infile = fopen(filepath.c_str(), "r");
    //Check if one can read this file;
    if (infile == NULL)
    {
        fprintf(stderr, "(Read ALL)Failed to read %s\n", filepath.c_str());
        return -1;
    }

    for (int i = 0; i < lines; ++i)
    {
        data_original[i] = 0.0;
        data_original_t[i] = 0.0;
    }

    // read in all Act data into data_original
    for (int i = 1; i < lines; ++i)
    {
        if (fscanf(infile, "%lf\t%lf\n", &tempx, &tempy) == 2)
        {
            // remove all nan numbers
            if (!std::isnan(tempy))
                data_original[i] = tempy;
            else
            {
                data_original[i] = data_original[i - 1];
                this->errors_act++;
            }
            // unit of time should be second
            data_original_t[i] = tempx;
        }
    }
    fclose(infile);
    for (int i = 0; i < lines; ++i)
    {
        actOri.data[i] = data_original[i];
        actOri.time[i] = data_original_t[i];
    }

    size_data = data_original_t[lines - 1] / 10;
    temp_data = new double[size_data];
    temp_data_t = new double[size_data];

    for (int i = 1; i < lines; i++)
    {
        diff_t = (data_original_t[i] - data_original_t[i - 1]) / 10;
        if (diff_t == 1)
        {
            temp_data_t[datamrk] = data_original_t[i];
            temp_data[datamrk] = data_original[i];
            datamrk++;
        }
        else if (diff_t < 100)
        {
            for (int j = 1; j <= diff_t; ++j)
            {
                temp_data[datamrk] = data_original[i - 1] +
                                     (data_original[i] - data_original[i - 1]) / diff_t * j;
                temp_data_t[datamrk] = data_original_t[i - 1] +
                                       (data_original_t[i] - data_original_t[i - 1]) / diff_t * j;
                datamrk++;
            }
        }
        if (diff_t > 100.0)
        {
            if (mice.batch == 3 || mice.batch == 5)
            {
                datamrk = 0;
                /* std::cout << "for Batch3 or 5" << std::endl; */
            }
            if (mice.batch == 6)
            {
                datamrk--;
                /* std::cout << "for Batch 6" << std::endl; */
                break;
            }
        }
    }
    // datamrk is the reduced size, and then the size_data is integer multiple 2 days.
    size_data = (datamrk / (48 * 360)) * 48 * 360;

    if (act.data != nullptr)
    {
        delete[] act.data;
        delete[] act.time;
    }
    act.size = size_data;
    act.data = new double[size_data];
    act.time = new double[size_data];

    for (int i = 0; i < size_data; ++i)
    {
        act.data[i] = temp_data[i];
        act.time[i] = temp_data_t[i] / 3600.0;
    }

    delete[] data_original;
    delete[] data_original_t;
    delete[] temp_data_t;
    delete[] temp_data;

    // Then start to read in the temperature data;
    filepath = mice.path + "/batch" + std::to_string(mice.batch) + "/" + mice.name + ".Temperature.txt";
    lines = fileLines(filepath);
    datamrk = 0;
    diff_t = 0;
    tempx = 0;
    tempy = 0;

    data_original = new double[lines];
    data_original_t = new double[lines];
    size_data = 0;

    tempOri.size = lines;
    if (tempOri.data != NULL)
    {
        delete[] tempOri.data;
        delete[] tempOri.time;
    }
    tempOri.data = new double[lines];
    tempOri.time = new double[lines];

    infile = fopen(filepath.c_str(), "r");
    //Check if one can read this file;
    if (infile == NULL)
    {
        fprintf(stderr, "(Read ALL)Failed to read %s\n", filepath.c_str());
        return -1;
    }

    for (int i = 0; i < lines; ++i)
    {
        data_original[i] = 0.0;
        data_original_t[i] = 0.0;
    }

    // read in all Temp data into data_original
    for (int i = 1; i < lines; ++i)
    {
        if (fscanf(infile, "%lf\t%lf\n", &tempx, &tempy) == 2)
        {
            // remove all nan numbers
            if (!std::isnan(tempy))
                data_original[i] = tempy;
            else
            {
                data_original[i] = data_original[i - 1];
                this->errors_temp++;
            }
            // unit of time should be second
            data_original_t[i] = tempx;
        }
    }
    fclose(infile);
    for (int i = 0; i < lines; ++i)
    {
        tempOri.data[i] = data_original[i];
        tempOri.time[i] = data_original_t[i];
    }

    size_data = data_original_t[lines - 1] / 10;
    temp_data = new double[size_data];
    temp_data_t = new double[size_data];

    for (int i = 1; i < lines; i++)
    {
        diff_t = (data_original_t[i] - data_original_t[i - 1]) / 10;
        if (diff_t == 1)
        {
            temp_data_t[datamrk] = data_original_t[i];
            temp_data[datamrk] = data_original[i];
            datamrk++;
        }
        else if (diff_t < 100)
        {
            for (int j = 1; j <= diff_t; ++j)
            {
                temp_data[datamrk] = data_original[i - 1] +
                                     (data_original[i] - data_original[i - 1]) / diff_t * j;
                temp_data_t[datamrk] = data_original_t[i - 1] +
                                       (data_original_t[i] - data_original_t[i - 1]) / diff_t * j;
                datamrk++;
            }
        }
        if (diff_t > 100.0)
        {
            if (mice.batch == 3 || mice.batch == 5)
            {
                datamrk = 0;
                /* std::cout << "for Batch3 or 5" << std::endl; */
            }
            if (mice.batch == 6)
            {
                datamrk--;
                /* std::cout << "for Batch 6" << std::endl; */
                break;
            }
        }
    }
    // datamrk is the reduced size, and then the size_data is integer multiple 2 days.
    size_data = (datamrk / (48 * 360)) * 48 * 360;

    if (temp.data != nullptr)
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
        temp.time[i] = temp_data_t[i] / 3600.0;
    }

    delete[] data_original;
    delete[] data_original_t;
    delete[] temp_data_t;
    delete[] temp_data;

    return true;
}

void group()
{
    int the_choosen_len = 2880;
    std::string path = "/home/yiwen/mice_physiology/micenew/rebuild/group_average/";
    std::unordered_map<std::string, int> data;
    std::ifstream file;
    file.open(path + "mice");
    if (!file.is_open())
        return;
    std::string name;
    std::string type;
    int i;
    while (!file.eof())
    {
        file >> name;
        file >> type;
        file >> i;
        data[name] = i;
    }
    int largest = 0;
    for (auto &a : data)
    {
        if (a.second > largest)
            largest = a.second;
    }
    std::cout << "largest: " << largest << std::endl;
    int reduced_len = the_choosen_len - largest;
    double *temp = new double[the_choosen_len];
    for (auto &a : data)
    {
        std::cout << "name: " << a.first << "\tdata:" << largest - a.second << "\tto\tend-" << a.second << std::endl;
        FILE *temp_data = fopen((path + a.first + "_activiry.txt").c_str(), "r");
        if (temp_data == NULL)
            continue;
        for (int i = 0; i < the_choosen_len; ++i)
        {
            fscanf(temp_data, "%lf\n", temp + i);
        }
        fclose(temp_data);
        temp_data = fopen((path + a.first + "_activity_reduce.txt").c_str(), "w");
        for (int j = 0; j < reduced_len; ++j)
        {
            fprintf(temp_data, "%lf\n", temp[j + largest - a.second]);
        }
        fclose(temp_data);
    }
    delete[] temp;
    file.close();
}

void get_activity_alignment(std::string name, std::string name2)
{
    int the_choosen_len = 2880;
    int draw_len = 4 * 24 * 4;
    /* FILE *config = fopen("/home/yiwen/mice_physiology/micenew/rebuild/group_average/mice", "r"); */
    std::string path = "/home/yiwen/mice_physiology/micenew/rebuild/group_average/";

    int len = fileLines(path + name + "_temperature.txt");
    len -= 1;
    FILE *data = fopen((path + name + "_temperature.txt").c_str(), "r");
    if (data == NULL)
        return;
    double *y = new double[len];
    for (int i = 0; i < len; ++i)
    {
        fscanf(data, "%lf\n", y + i);
    }
    if (len != the_choosen_len)
    {
        /* double *ny = new double[the_choosen_len]; */
        /* int idx = 0; */
        /* for (int j = 0; j < the_choosen_len; ++j) { */
        /*     if () */
        /* } */
        printf("Fuck u %d != %d\n", len, the_choosen_len);
    }
    double *x = new double[the_choosen_len];
    for (int i = 0; i < the_choosen_len; ++i)
    {
        x[i] = i;
    }
    fclose(data);

    data = fopen((path + name2 + "_temperature.txt").c_str(), "r");
    if (data == NULL)
    {
        delete[] x;
        delete[] y;
        return;
    }
    double *y2 = new double[the_choosen_len];
    for (int i = 0; i < the_choosen_len; ++i)
    {
        fscanf(data, "%lf\n", y2 + i);
    }
    fclose(data);

    auto cv = new TCanvas("cv");
    cv->cd();
    auto graph1 = new TGraph(draw_len, x, y);
    graph1->SetLineColor(kRedBlue);
    graph1->Draw("AL");

    auto graph2 = new TGraph(draw_len, x, y2);
    graph2->SetLineColor(kGreenRedViolet);
    graph2->Draw("SAMEL");
    cv->Modified();
    cv->Update();
    int op = 0;
    int offset = 0;
    while (1)
    {
        scanf("%d", &op);
        offset += op;
        if (op == 0)
            break;
        double *tmpx = graph2->GetX();
        for (int i = 0; i < draw_len; ++i)
        {
            tmpx[i] += op;
        }
        graph2->Draw("SAMEL");
        cv->Modified();
        cv->Update();
    }
    printf("Offset: %d\n", offset);

    cv->Print("Fck.pdf");

    delete graph1;
    delete graph2;

    delete cv;
    delete[] y;
    delete[] y2;
    delete[] x;
}

// calculate the group average of activity. done !
void group_average()
{
    int the_choosen_len = 2790;
    std::string path = "/home/yiwen/mice_physiology/micenew/rebuild/group_average/";
    std::unordered_map<std::string, mice_type> data;
    std::ifstream file;
    file.open(path + "mice");
    if (!file.is_open())
        return;
    std::string name;
    std::string type;
    int i;
    while (!file.eof())
    {
        file >> name;
        file >> type;
        file >> i;
        if (type == "mutant")
            data[name] = mice_type_mutant;
        else if (type == "hyper")
            data[name] = mice_type_mutant_hyper;
        else if (type == "control")
            data[name] = mice_type_control;
    }
    double *mutant_data = new double[the_choosen_len];
    double *control_data = new double[the_choosen_len];
    double *mutant_hyper_data = new double[the_choosen_len];
    // double *x = new double [the_choosen_len];

    for (int i = 0; i < the_choosen_len; ++i)
    {
        mutant_data[i] = 0;
        control_data[i] = 0;
        mutant_hyper_data[i] = 0;
        // x[i] = 0;
    }

    int mutant_count = 0, control_count = 0, hyper_count = 0;
    for (auto &a : data)
    {
        if (a.second & mice_type_mutant)
            mutant_count++;
        else if (a.second & mice_type_control)
            control_count++;
        else if (a.second & mice_type_mutant_hyper)
            hyper_count++;
        else
        {
            std::cout << "Error: wrong mice type" << std::endl;
            return;
        }

        std::cout << "name: " << a.first << "\ttype:" << (a.second == mice_type_mutant ? "mutant" : (a.second == mice_type_mutant_hyper ? "hyper" : "control")) << std::endl;

        double temp;
        FILE *temp_data = fopen((path + a.first + "_activity_reduce.txt").c_str(), "r");
        if (temp_data == NULL)
        {
            std::cerr << "File doesnot exist: " + (path + a.first + "_activity_reduce.txt") << std::endl;
        }

        if (a.second == mice_type_mutant)
            for (int i = 0; i < the_choosen_len; ++i)
            {
                fscanf(temp_data, "%lf\n", &temp);
                mutant_data[i] += temp;
                if (std::isnan(temp))
                    std::cout << a.first << " what the hell?" << std::endl;
            }
        else if (a.second == mice_type_control)
            for (int i = 0; i < the_choosen_len; ++i)
            {
                fscanf(temp_data, "%lf\n", &temp);
                control_data[i] += temp;
                if (std::isnan(temp))
                    std::cout << a.first << " what the hell?" << std::endl;
            }
        else if (a.second == mice_type_mutant_hyper)
            for (int i = 0; i < the_choosen_len; ++i)
            {
                fscanf(temp_data, "%lf\n", &temp);
                mutant_hyper_data[i] += temp;
                if (std::isnan(temp))
                    std::cout << a.first << " what the hell?" << std::endl;
            }
        fclose(temp_data);
    }
    FILE *mutant = fopen((path + "mutant_group_average_activiry.txt").c_str(), "w");
    FILE *control = fopen((path + "control_group_average_activiry.txt").c_str(), "w");
    FILE *mutant_hyper = fopen((path + "mutant_hyper_group_average_activiry.txt").c_str(), "w");
    for (int i = 0; i < the_choosen_len; ++i)
    {
        mutant_data[i] /= mutant_count;
        control_data[i] /= control_count;
        mutant_hyper_data[i] /= hyper_count;
        fprintf(mutant, "%lf\n", mutant_data[i]);
        fprintf(control, "%lf\n", control_data[i]);
        fprintf(mutant_hyper, "%lf\n", mutant_hyper_data[i]);
    }

    fclose(mutant);
    fclose(control);
    fclose(mutant_hyper);
    delete[] mutant_data;
    delete[] control_data;
    delete[] mutant_hyper_data;
}

// do some dfa for group average of activity
void group_average_dfa(std::string name)
{
    // calculate average then the dfa plot.
    std::string path = "/home/yiwen/mice_physiology/micenew/rebuild/group_average/";
    if (name != "mutant")
        name = "control";
    int lines = fileLines(path + name + "_group_average_activiry.txt");
    FILE *file = fopen((path + name + "_group_average_activiry.txt").c_str(), "r");
    double *data = new double[lines];
    for (int i = 0; i < lines; ++i)
    {
        fscanf(file, "%lf\n", data + i);
    }
    int days = lines / 96;
    int the_choosen_len = days * 96;

    TimeSeq timeseq;
    timeseq.size = the_choosen_len;
    timeseq.data = new double[the_choosen_len];
    timeseq.time = new double[the_choosen_len];

    for (int i = 0; i < the_choosen_len; ++i)
    {
        timeseq.data[i] = data[i];
        timeseq.time[i] = i * 2 * M_PI / 96;
    }

    double *waveform = new double[96];
    for (int i = 0; i < 96; ++i)
    {
        waveform[i] = 0;
    }
    for (int i = 0; i < the_choosen_len; ++i)
    {
        waveform[i % 96] += timeseq.data[i];
    }

    for (int i = 0; i < 96; ++i)
    {
        waveform[i] /= days;
    }
    for (int i = 0; i < the_choosen_len; ++i)
    {
        timeseq.data[i] -= waveform[i % 96];
    }

    struct period_dfa a;
    DFA_filter_plot(data, the_choosen_len, 2, name + " activity", &a, timeseq);

    delete[] data;
    delete[] timeseq.data;
    delete[] timeseq.time;
}

// calculate the group average of dfa data
void dfa_group_average(std::string midfix, std::string suffix)
{
    // do dfa and then average the dfa data.
    int the_choosen_len = 34;
    if (midfix == "dfa_circadian")
        the_choosen_len = 20;
    std::string path = "/home/yiwen/mice_physiology/micenew/rebuild/group_average/";
    std::unordered_map<std::string, mice_type> data;
    std::ifstream file;
    file.open(path + "mice");
    if (!file.is_open())
        return;
    std::string name;
    std::string type;
    int i;
    while (!file.eof())
    {
        file >> name;
        file >> type;
        file >> i;
        // std::cout << name << "\t" << type << "\t" << i << std::endl;
        if (type == "mutant")
            data[name] = mice_type_mutant;
        else if (type == "hyper")
            data[name] = mice_type_mutant_hyper;
        else if (type == "control")
            data[name] = mice_type_control;
    }
    double *mutant_data = new double[the_choosen_len];
    double *control_data = new double[the_choosen_len];
    double *mutant_hyper_data = new double[the_choosen_len];
    double *x = new double[the_choosen_len];

    for (int i = 0; i < the_choosen_len; ++i)
    {
        mutant_data[i] = 0;
        control_data[i] = 0;
        mutant_hyper_data[i] = 0;
        x[i] = 0;
    }

    std::string file_name;

    std::unordered_map<std::string, mice_type>::iterator it = data.begin();
    if (midfix != "dfa_circadian")
        file_name = path + it->first + "_activity_" + midfix + "_data_" + suffix + ".txt";
    else
        file_name = path + it->first + "_" + midfix + "_data_" + suffix + ".txt";

    FILE *temp_data = fopen(file_name.c_str(), "r");
    if (temp_data == NULL)
    {
        std::cerr << "File doesnot exist: " << file_name << std::endl;
        return;
    }

    double tempx, tempy;
    for (int i = 0; i < the_choosen_len; ++i)
    {
        fscanf(temp_data, "%lf\t%lf\n", &tempx, &tempy);
        x[i] = tempx;
    }

    int mutant_count = 0, control_count = 0, hyper_count = 0;
    for (auto &a : data)
    {
        if (a.second & mice_type_mutant)
            mutant_count++;
        else if (a.second & mice_type_control)
            control_count++;
        else if (a.second & mice_type_mutant_hyper)
            hyper_count++;
        else
        {
            std::cout << "Error: wrong mice type" << std::endl;
            return;
        }

        if (midfix != "dfa_circadian")
            file_name = path + a.first + "_activity_" + midfix + "_data_" + suffix + ".txt";
        else
            file_name = path + a.first + "_" + midfix + "_data_" + suffix + ".txt";

        printf("file name: %s\n", file_name.c_str());

        std::cout << "name: " << a.first << "\ttype:" << (a.second == mice_type_mutant ? "mutant" : (a.second == mice_type_mutant_hyper ? "hyper" : "control")) << std::endl;

        temp_data = fopen(file_name.c_str(), "r");
        if (temp_data == NULL)
        {
            std::cerr << "File doesnot exist: " + file_name << std::endl;
        }

        if (a.second == mice_type_mutant)
        {
            for (int i = 0; i < the_choosen_len; ++i)
            {
                fscanf(temp_data, "%lf\t%lf\n", &tempx, &tempy);
                mutant_data[i] += tempy;
                if (std::isnan(tempy))
                    std::cout << a.first << " what the hell?" << std::endl;
            }
        }
        else if (a.second == mice_type_control)
        {
            for (int i = 0; i < the_choosen_len; ++i)
            {
                fscanf(temp_data, "%lf\t%lf\n", &tempx, &tempy);
                control_data[i] += tempy;
                if (std::isnan(tempy))
                    std::cout << a.first << " what the hell?" << std::endl;
            }
        }
        else if (a.second == mice_type_mutant_hyper)
        {
            for (int i = 0; i < the_choosen_len; ++i)
            {
                fscanf(temp_data, "%lf\t%lf\n", &tempx, &tempy);
                mutant_hyper_data[i] += tempy;
                if (std::isnan(tempy))
                    std::cout << a.first << " what the hell?" << std::endl;
            }
        }
        else
        {
            std::cout << "Error: wrong mice type" << std::endl;
            return;
        }

        fclose(temp_data);
    }

    FILE *mutant = fopen((path + "mutant_group_average_" + midfix + "_data_" + suffix + ".txt").c_str(), "w");
    FILE *control = fopen((path + "control_group_average_" + midfix + "_data_" + suffix + ".txt").c_str(), "w");
    FILE *mutant_hyper = fopen((path + "mutant_hyper_group_average_" + midfix + "_data_" + suffix + ".txt").c_str(), "w");
    printf("Control\t: %d\n", control_count);
    printf("Mutant\t: %d\n", mutant_count);
    printf("Hyper\t: %d\n", hyper_count);

    for (int i = 0; i < the_choosen_len; ++i)
    {
        mutant_data[i] /= mutant_count;
        control_data[i] /= control_count;
        mutant_hyper_data[i] /= hyper_count;
        fprintf(mutant, "%lf\t%lf\n", x[i], mutant_data[i]);
        fprintf(control, "%lf\t%lf\n", x[i], control_data[i]);
        fprintf(mutant_hyper, "%lf\t%lf\n", x[i], mutant_hyper_data[i]);
        printf("%lf\t", x[i]);
        printf("%lf\t%lf\t%lf\n", control_data[i], mutant_data[i], mutant_hyper_data[i]);
    }

    fclose(mutant);
    fclose(control);
    fclose(mutant_hyper);
    delete[] mutant_data;
    delete[] control_data;
    delete[] mutant_hyper_data;
}

// Draw the group average of dfa
void dfa_group_average_plot(std::string midfix, std::string suffix = "ori", int dfa_order = 2, std::string plot_mark = "A")
{
    std::string type;

    // Draw the averaged dfa data into plots.
    std::string path = "/home/yiwen/mice_physiology/micenew/rebuild/group_average/";
    int the_choosen_len = 34;
    int reduced_len = 0;
    std::vector<double> dfax(the_choosen_len);
    std::vector<double> dfay(the_choosen_len);
    std::vector<double> xp(the_choosen_len);
    std::vector<double> yp(the_choosen_len);

    auto canvas = std::make_unique<TCanvas>("canvas");
    double scale = 0.7;
    canvas->SetTitle("Hello World");
    gStyle->SetOptStat(0);
    canvas->SetCanvasSize(450 * scale, 500 * scale);

    double lim = 48;
    double lim2 = 96;
    double limit = log10(lim);

    double *dfahour = nullptr;

    type = "control";
    FILE *temp_data = fopen((path + type + "_group_average_" + midfix + "_data_" + suffix + ".txt").c_str(), "r");
    if (temp_data == NULL)
    {
        std::cerr << "Cant open file: " << path + type + "_group_average_" + midfix + "_data_" + suffix + ".txt" << std::endl;
        return;
    }
    double tempx, tempy;
    for (int i = 0; i < the_choosen_len; ++i)
    {
        fscanf(temp_data, "%lf\t%lf\n", &tempx, &tempy);
        dfax[i] = tempx;
        dfay[i] = tempy;
        xp[i] = log10(dfax[i]);
        yp[i] = log10(dfay[i]);
    }
    for (; reduced_len < the_choosen_len; reduced_len++)
    {
        if (dfax[reduced_len] >= 4 * 36)
        {
            break;
        }
    }

    fclose(temp_data);
    auto graph_control = std::make_shared<TGraph>(reduced_len, dfax.data(), dfay.data());
    int ori_data_leng = 2790;
    std::vector<double> ori_data(ori_data_leng);
    temp_data = fopen((path + type + "_group_average_activiry.txt").c_str(), "r");
    if (temp_data == NULL)
    {
        std::cerr << "Cant open file: " << path + type + "_group_average_activiry.txt" << std::endl;
        return;
    }
    for (int i = 0; i < ori_data_leng; ++i)
    {
        fscanf(temp_data, "%lf\n", &tempy);
        ori_data[i] = tempy;
    }
    fclose(temp_data);
    unsigned num = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(ori_data.begin(), ori_data.end(), std::default_random_engine(num));

    auto newdfa = std::make_unique<DFA>(ori_data.data(), ori_data.size(), 2);
    auto graph_suffule = std::make_shared<TGraph>(newdfa->get_size_o(), newdfa->dfax, newdfa->dfay);

    type = "mutant";
    temp_data = fopen((path + type + "_group_average_" + midfix + "_data_" + suffix + ".txt").c_str(), "r");
    if (temp_data == NULL)
    {
        std::cerr << "Cant open file: " << path + type + "_group_average_" + midfix + "_data_" + suffix + ".txt" << std::endl;
        return;
    }
    for (int i = 0; i < the_choosen_len; ++i)
    {
        fscanf(temp_data, "%lf\t%lf\n", &tempx, &tempy);
        dfax[i] = tempx;
        dfay[i] = tempy;
        xp[i] = log10(dfax[i]);
        yp[i] = log10(dfay[i]);
    }
    fclose(temp_data);
    auto graph_mutant = std::make_shared<TGraph>(reduced_len, dfax.data(), dfay.data());

    type = "mutant_hyper";
    temp_data = fopen((path + type + "_group_average_" + midfix + "_data_" + suffix + ".txt").c_str(), "r");
    if (temp_data == NULL)
    {
        std::cerr << "Cant open file: " << path + type + "_group_average_" + midfix + "_data_" + suffix + ".txt" << std::endl;
        return;
    }
    for (int i = 0; i < the_choosen_len; ++i)
    {
        fscanf(temp_data, "%lf\t%lf\n", &tempx, &tempy);
        dfax[i] = tempx;
        assert(tempy > 0);
        dfay[i] = tempy;
        xp[i] = log10(dfax[i]);
        yp[i] = log10(dfay[i]);
    }
    fclose(temp_data);
    double pars[5];
    auto f1 = new TF1("f1", fit_crossover, xp[0], limit, 5);
    auto *plot = new TGraph(the_choosen_len, xp.data(), yp.data());
    f1->SetParLimits(0, 0, 2);
    f1->SetParameter(0, 1.0);
    f1->SetParLimits(2, 0, 1);
    f1->SetParameter(2, 0.5);
    f1->SetParLimits(4, xp[2], limit);
    f1->SetParameter(4, xp[8]);

    double chisq = 0;
    for (int i = 0; i < 20; ++i)
    {
        plot->Fit(f1, "R");
        chisq = f1->GetChisquare();
        if (chisq < 0.0001)
            break;
    }
    delete plot;

    for (int i = 0; i < 5; ++i)
    {
        pars[i] = f1->GetParameter(i);
        std::cout << "Pars" << i << ": " << pars[i] << std::endl;
    }

    auto f2 = new TF1("f2", draw_crossover, dfax[0] / 4.0, pow(10, pars[4]) / 4.0, 5);
    auto f3 = new TF1("f3", draw_crossover, pow(10, pars[4]) / 4.0, lim / 4.0, 5);
    for (int i = 0; i < 4; ++i)
    {
        f2->SetParameter(i, pars[i]);
        f3->SetParameter(i, pars[i]);
    }
    f2->SetParameter(4, pow(10, pars[4]));
    f3->SetParameter(4, pow(10, pars[4]));
    double rvalue, lvalue;
    lvalue = pow(10, pars[4]);
    for (int i = 0; i < the_choosen_len; ++i)
    {
        if (lvalue > dfax[i])
            continue;
        else
        {
            lvalue = dfax[i - 1];
            rvalue = dfax[i];
            break;
        }
    }
    auto graph_hyper = std::make_shared<TGraph>(reduced_len, dfax.data(), dfay.data());

    dfahour = graph_control->GetX();
    for (int i = 0; i < reduced_len; ++i)
    {
        dfahour[i] /= 4.0;
    }
    dfahour = graph_mutant->GetX();
    for (int i = 0; i < reduced_len; ++i)
    {
        dfahour[i] /= 4.0;
    }
    dfahour = graph_hyper->GetX();
    for (int i = 0; i < reduced_len; ++i)
    {
        dfahour[i] /= 4.0;
    }
    dfahour = graph_suffule->GetX();
    for (int i = 0; i < reduced_len; ++i)
    {
        dfahour[i] /= 4.0;
    }

    if (midfix == "dfa")
        graph_control->SetTitle(("Correlation analysis group average of " + (type == "mutant_hyper" ? "hyper" : type)).c_str());
    else if (midfix == "dfamag")
        graph_control->SetTitle(("Magnitude analysis group average of " + (type == "mutant_hyper" ? "hyper" : type)).c_str());
    graph_control->SetName("graph_detrend");
    graph_control->SetMarkerStyle(72);
    graph_control->SetMarkerColor(kGreen);
    graph_control->SetMarkerSize(1);
    graph_mutant->SetName("graph_detrend");
    graph_mutant->SetMarkerStyle(71);
    graph_mutant->SetMarkerColor(kBlue);
    graph_mutant->SetMarkerSize(1);
    graph_hyper->SetName("graph_detrend");
    graph_hyper->SetMarkerStyle(73);
    graph_hyper->SetMarkerColor(kRed);
    graph_hyper->SetMarkerSize(1);
    graph_suffule->SetMarkerStyle(74);
    graph_suffule->SetMarkerColor(kBlack);
    graph_suffule->SetMarkerSize(1);

    std::string name;
    if (type == "mutant")
    {
        name = "12Otx2";
    }
    else if (type == "mutant_hyper")
    {
        name = "51M";
    }
    else if (type == "control")
    {
        name = "53";
    }
    else
        exit(-1);

    name = "53";
    temp_data = fopen((path + name + "_activity_" + midfix + "_data_" + suffix + ".txt").c_str(), "r");
    for (int i = 0; i < the_choosen_len; ++i)
    {
        fscanf(temp_data, "%lf\t%lf\n", &tempx, &tempy);
        dfax[i] = tempx;
        dfay[i] = tempy;
        xp[i] = log10(dfax[i]);
        yp[i] = log10(dfay[i]);
    }
    fclose(temp_data);
    auto graph_inv_control = std::make_shared<TGraph>(reduced_len, dfax.data(), dfay.data());

    temp_data = fopen((path + name + "_activity_reduce.txt").c_str(), "r");
    if (temp_data == NULL)
    {
        std::cerr << "Cant open file: " << path + name + "_activity_reduce.txt" << std::endl;
        perror(("Cant open file: " + path + name + "_activity_reduce.txt").c_str());
        return;
    }
    for (int i = 0; i < ori_data_leng; ++i)
    {
        fscanf(temp_data, "%lf\n", &tempy);
        ori_data[i] = tempy;
    }
    fclose(temp_data);
    num = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(ori_data.begin(), ori_data.end(), std::default_random_engine(num));

    newdfa = std::make_unique<DFA>(ori_data.data(), ori_data.size(), 2);
    auto control_inv_suffule = std::make_shared<TGraph>(newdfa->get_size_o(), newdfa->dfax, newdfa->dfay);

    name = "12Otx2";
    temp_data = fopen((path + name + "_activity_" + midfix + "_data_" + suffix + ".txt").c_str(), "r");
    for (int i = 0; i < the_choosen_len; ++i)
    {
        fscanf(temp_data, "%lf\t%lf\n", &tempx, &tempy);
        dfax[i] = tempx;
        dfay[i] = tempy;
        xp[i] = log10(dfax[i]);
        yp[i] = log10(dfay[i]);
    }
    fclose(temp_data);
    auto graph_inv_mutant = std::make_shared<TGraph>(reduced_len, dfax.data(), dfay.data());

    name = "64M";
    temp_data = fopen((path + name + "_activity_" + midfix + "_data_" + suffix + ".txt").c_str(), "r");
    for (int i = 0; i < the_choosen_len; ++i)
    {
        fscanf(temp_data, "%lf\t%lf\n", &tempx, &tempy);
        dfax[i] = tempx;
        dfay[i] = tempy;
        xp[i] = log10(dfax[i]);
        yp[i] = log10(dfay[i]);
    }
    fclose(temp_data);
    auto graph_inv_hyper = std::make_shared<TGraph>(reduced_len, dfax.data(), dfay.data());

    dfahour = graph_inv_control->GetX();
    for (int i = 0; i < reduced_len; ++i)
    {
        dfahour[i] /= 4.0;
    }
    dfahour = graph_inv_mutant->GetX();
    for (int i = 0; i < reduced_len; ++i)
    {
        dfahour[i] /= 4.0;
    }
    dfahour = graph_inv_hyper->GetX();
    for (int i = 0; i < reduced_len; ++i)
    {
        dfahour[i] /= 4.0;
    }
    dfahour = control_inv_suffule->GetX();
    for (int i = 0; i < reduced_len; ++i)
    {
        dfahour[i] /= 4.0;
    }

    if (midfix == "dfa")
        graph_inv_control->SetTitle(("Correlation analysis individual of " + (type == "mutant_hyper" ? "hyper" : type)).c_str());
    else if (midfix == "dfamag")
        graph_inv_control->SetTitle(("Magnitude analysis individual of " + (type == "mutant_hyper" ? "hyper" : type)).c_str());

    graph_inv_control->SetName("graph_detrend");
    graph_inv_control->SetMarkerStyle(72);
    graph_inv_control->SetMarkerColor(kGreen);
    graph_inv_control->SetMarkerSize(1);
    graph_inv_mutant->SetName("graph_detrend");
    graph_inv_mutant->SetMarkerStyle(71);
    graph_inv_mutant->SetMarkerColor(kBlue);
    graph_inv_mutant->SetMarkerSize(1);
    graph_inv_hyper->SetName("graph_detrend");
    graph_inv_hyper->SetMarkerStyle(73);
    graph_inv_hyper->SetMarkerColor(kRed);
    graph_inv_hyper->SetMarkerSize(1);

    control_inv_suffule->SetMarkerStyle(74);
    control_inv_suffule->SetMarkerColor(kBlack);
    control_inv_suffule->SetMarkerSize(1);
    canvas->Divide(1, 2, 0, 0);

    canvas->cd(2)->SetLogx(1);
    canvas->cd(2)->SetLogy(1);
    std::cout << "Draw First" << std::endl;

    // graph_control->GetYaxis()->SetTitle("F(n)");
    // graph_control->GetYaxis()->CenterTitle();
    // graph_control->GetXaxis()->SetTitle("Time(h)");
    // graph_control->GetXaxis()->CenterTitle();

    auto gaussian_curve = std::make_shared<TF1>("gaussian", "0.2*x^0.5", 0, 1e10);
    auto pink_curve = std::make_shared<TF1>("pink", "4*x^1", 0, 1e10);

    gaussian_curve->SetLineColor(kBlack);
    gaussian_curve->SetLineStyle(kDashed);
    pink_curve->SetLineColor(kBlack);

    auto multigr2 = std::make_shared<TMultiGraph>();
    multigr2->Add(graph_control.get());
    multigr2->Add(graph_mutant.get());
    multigr2->Add(graph_hyper.get());
    multigr2->Add(graph_suffule.get());
    multigr2->Draw("AP");
    gaussian_curve->Draw("SAME");
    pink_curve->Draw("SAME");

    // multigr2->GetYaxis()->SetTitle("F(n)");
    multigr2->GetYaxis()->CenterTitle();
    // multigr2->GetXaxis()->SetTitle("Time(h)");
    multigr2->GetXaxis()->CenterTitle();
    multigr2->GetXaxis()->SetLimits(dfax[0] * 0.9 / 4, dfax[reduced_len - 1] / 4 * 1.1);
    auto legend2 = std::make_unique<TLegend>(0.1, 0.88, 0.34, 1);
    legend2->SetName("legend");

    legend2->AddEntry(pink_curve.get(), "alpha = 1.0", "l");
    legend2->AddEntry(gaussian_curve.get(), "alpha = 0.5", "l");
    // legend->SetBorderSize(0);
    legend2->Draw();
    // auto l = std::make_unique<TLine>(pow(10.0,pars[4]) * 0.010416667*24,-3,pow(10.0,pars[4]) * 0.010416667*24,100);
    // if (type != "mutant_hyper" && midfix != "dfamag")
    //     l->Draw("SAMEL");

    // Draw Second
    canvas->cd(1)->SetLogx(1);
    canvas->cd(1)->SetLogy(1);

    // graph_inv_control->GetYaxis()->SetTitle("F(n)");
    graph_inv_control->GetYaxis()->CenterTitle();
    // graph_inv_control->GetXaxis()->SetTitle("Time(h)");
    graph_inv_control->GetXaxis()->CenterTitle();

    auto multigr1 = std::make_shared<TMultiGraph>();
    multigr1->Add(graph_inv_control.get());
    multigr1->Add(graph_inv_mutant.get());
    multigr1->Add(graph_inv_hyper.get());
    multigr1->Add(control_inv_suffule.get());
    multigr1->Draw("AP");
    gaussian_curve->Draw("SAME");
    pink_curve->Draw("SAME");

    // multigr1->GetYaxis()->SetTitle("F(n)");
    multigr1->GetYaxis()->CenterTitle();
    multigr1->GetXaxis()->SetLimits(dfax[0] * 0.9 / 4, dfax[reduced_len - 1] / 4 * 1.1);

    auto legend = std::make_unique<TLegend>(0.1, 0.76, 0.47, 1);
    legend->SetName("legend");

    legend->AddEntry(graph_inv_control.get(), "control", "p");
    legend->AddEntry(graph_inv_mutant.get(), "mutant - normal activity", "p");
    legend->AddEntry(graph_inv_hyper.get(), "mutant - hyper activity", "p");
    legend->AddEntry(control_inv_suffule.get(), "shuffled", "p");
    // legend->SetBorderSize(0);
    legend->Draw();

    canvas->cd();

    // canvas->DrawFrame(0, 0, 1, 1);
    auto latex = std::make_shared<TLatex>();
    latex->SetTextFont(42);
    latex->SetTextAlign(13);
    latex->SetTextSize(0.045);
    latex->DrawLatex(0, 1, plot_mark.c_str());

    latex->SetTextAlign(11);
    latex->DrawLatex(0.44, 0, "Time(h)");

    latex->SetTextAngle(90);
    latex->SetTextAlign(13);
    latex->DrawLatex(0, 0.22, "Detrended Fluctuation Function F(n)");

    canvas->Print((midfix + "_group_average_" + suffix + "_data.pdf").c_str());

    dfax.clear();
    xp.clear();
    dfay.clear();
    yp.clear();
}

double *get_powerspec(double *data, int size)
{
    fftw_complex *middle = fftw_alloc_complex(size);
    for (size_t i = 0; i < size; i++)
    {
        middle[i][0] = 0.0;
        middle[i][1] = 0.0;
    }
    fftw_plan plan = fftw_plan_dft_r2c_1d(size, data, middle, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    auto powers = new double[size];

    for (size_t i = 0; i < size; i++)
    {
        powers[i] = middle[i][0] * middle[i][0] + middle[i][1] * middle[i][1];
    }

    fftw_free(middle);

    return powers;
}

// Draw the powerspec and waveform for grouopaverage and individuals
void powerspec_waveform_plot()
{
    std::string path = "/home/yiwen/mice_physiology/micenew/rebuild/group_average/";
    int the_choosen_len = 2790;
    // individual
    std::string control_mouse = "57";
    std::string mutant_mouse = "52M";
    std::string mutant_hyper_mouse = "51M";
    double *tempx = nullptr;

    // readin activity data
    std::vector<double> control(the_choosen_len);
    std::vector<double> mutant(the_choosen_len);
    std::vector<double> hyper(the_choosen_len);

    double alignment = 0;

    FILE *handle = fopen((path + control_mouse + "_activity_reduce.txt").c_str(), "r");
    if (handle == NULL)
    {
        std::cerr << "Cant open this file: " << path + control_mouse + "_activity_reduce.txt" << std::endl;
        return;
    }
    for (size_t i = 0; i < the_choosen_len; i++)
    {
        fscanf(handle, "%lf\n", &control.data()[i]);
    }
    fclose(handle);
    auto powers = get_powerspec(control.data(), the_choosen_len);
    normalize_double(powers, the_choosen_len);
    alignment = gsl_stats_max(powers, 1, the_choosen_len) * 1.1;
    // Set the proper frequency for the powerspec.
    auto x = new double[the_choosen_len];
    for (size_t i = 0; i < the_choosen_len; i++)
    {
        x[i] = the_choosen_len / (4.0) / (i + 1.0);
    }

    auto waveform = new double[96];
    auto waveformx = new double[96];
    for (size_t i = 0; i < 96; i++)
    {
        waveformx[i] = i;
    }

    auto control_plot = std::make_shared<TGraph>(the_choosen_len, x, powers);
    control_plot->SetName("control_plot");
    control_plot->SetLineColor(kGreen);
    get_waveform(control.data(), the_choosen_len, waveform, 96);
    auto control_waveform = std::make_unique<TGraph>(96, waveformx, waveform);
    control_waveform->SetLineColor(kGreen);
    control_waveform->SetTitle("");
    delete[] powers;

    handle = fopen((path + mutant_mouse + "_activity_reduce.txt").c_str(), "r");
    if (handle == NULL)
    {
        std::cerr << "Cant open this file: " << path + mutant_mouse + "_activity_reduce.txt" << std::endl;
        return;
    }
    for (size_t i = 0; i < the_choosen_len; i++)
    {
        fscanf(handle, "%lf\n", &mutant.data()[i]);
    }
    fclose(handle);
    powers = get_powerspec(mutant.data(), the_choosen_len);

    normalize_double(powers, the_choosen_len);
    for (size_t i = 0; i < the_choosen_len; i++)
    {
        powers[i] += alignment;
    }
    alignment += (gsl_stats_max(powers, 1, the_choosen_len) - alignment) * 1.1;

    auto mutant_plot = std::make_shared<TGraph>(the_choosen_len, x, powers);
    mutant_plot->SetName("mutant_plot");
    mutant_plot->SetLineColor(kBlue);
    get_waveform(mutant.data(), the_choosen_len, waveform, 96);
    auto mutant_waveform = std::make_unique<TGraph>(96, waveformx, waveform);
    mutant_waveform->SetLineColor(kBlue);
    mutant_waveform->SetTitle("");
    delete[] powers;

    handle = fopen((path + mutant_hyper_mouse + "_activity_reduce.txt").c_str(), "r");
    if (handle == NULL)
    {
        std::cerr << "Cant open this file: " << path + mutant_hyper_mouse + "_activity_reduce.txt" << std::endl;
        return;
    }
    for (size_t i = 0; i < the_choosen_len; i++)
    {
        fscanf(handle, "%lf\n", &hyper.data()[i]);
    }
    fclose(handle);
    powers = get_powerspec(hyper.data(), the_choosen_len);
    normalize_double(powers, the_choosen_len);
    for (size_t i = 0; i < the_choosen_len; i++)
    {
        powers[i] += alignment;
    }
    alignment = 0;

    auto hyper_plot = std::make_shared<TGraph>(the_choosen_len, x, powers);
    hyper_plot->SetName("hyper_plot");
    hyper_plot->SetLineColor(kRed);
    get_waveform(hyper.data(), the_choosen_len, waveform, 96);
    auto hyper_waveform = std::make_unique<TGraph>(96, waveformx, waveform);
    hyper_waveform->SetLineColor(kRed);
    hyper_waveform->SetTitle("");
    delete[] powers;

    tempx = mutant_waveform->GetX();
    for (size_t i = 0; i < 96; i++)
    {
        tempx[i] /= 4.0;
    }
    tempx = control_waveform->GetX();
    for (size_t i = 0; i < 96; i++)
    {
        tempx[i] /= 4.0;
    }
    tempx = hyper_waveform->GetX();
    for (size_t i = 0; i < 96; i++)
    {
        tempx[i] /= 4.0;
    }

    /**
     * Set the graph elements, and plot.
    */
    auto canvas = std::make_unique<TCanvas>("canvas", "canvas", 160, 90);
    canvas->SetFillStyle(0);
    canvas->SetFrameFillStyle(0);
    canvas->Divide(2, 1);
    canvas->cd(1)->SetPad("padl", "padl", 0, 0, 0.7, 1);
    canvas->cd(2)->SetPad("padr", "padr", 0.7, 0, 1, 1);

    auto padl = (TPad *)gROOT->FindObject("padl");
    auto padr = (TPad *)gROOT->FindObject("padr");
    padl->SetFillStyle(4000);
    padl->SetFrameFillStyle(4000);
    padr->SetFillStyle(4000);
    padr->SetFrameFillStyle(4000);

    auto legend = std::make_unique<TLegend>(0.15, 0.7, 0.35, 0.85);
    legend->SetName("legend");
    legend->AddEntry(control_plot.get(), "control", "l");
    legend->AddEntry(mutant_plot.get(), "mutant", "l");
    legend->AddEntry(hyper_plot.get(), "hyper", "l");
    legend->SetBorderSize(0);
    auto multigraph = std::make_unique<TMultiGraph>();

    multigraph->Add(control_plot.get());
    multigraph->Add(mutant_plot.get());
    multigraph->Add(hyper_plot.get());

    canvas->cd(1);
    padl->Draw();
    padl->cd()->SetMargin(0.1, 0.01, 0.1, 0.1);
    padl->cd()->SetLogx(1);
    multigraph->GetXaxis()->SetRangeUser(0.5, 36);
    multigraph->GetXaxis()->SetTitle("Time scale (h)");
    multigraph->GetXaxis()->CenterTitle();
    multigraph->GetYaxis()->SetTitle("Power spectrum");
    multigraph->GetYaxis()->CenterTitle();
    multigraph->Draw("AL");
    legend->Draw();

    canvas->cd(2);
    padr->Draw();
    // padr->SetMargin(0.1, 0.01, 0.1, 0.1);
    padr->Divide(1, 3, 0, 0);

    padr->cd(3)->SetFillStyle(4000);
    padr->cd(3)->SetFrameFillStyle(4000);
    padr->cd(3)->SetMargin(0.01, 0.01, 0.25, 0.0);
    control_waveform->GetXaxis()->SetLabelSize(0.06);

    control_waveform->Draw("AL");
    padr->cd(2)->SetFillStyle(4000);
    padr->cd(2)->SetFrameFillStyle(4000);
    padr->cd(2)->SetMargin(0.01, 0.01, 0.0, 0.0);
    mutant_waveform->GetXaxis()->SetLabelSize(0);

    mutant_waveform->Draw("AL");
    padr->cd(1)->SetFillStyle(4000);
    padr->cd(1)->SetFrameFillStyle(4000);
    padr->cd(1)->SetMargin(0.01, 0.01, 0.0, 0.2);
    hyper_waveform->GetXaxis()->SetLabelSize(0);

    hyper_waveform->Draw("AL");

    canvas->Update();
    canvas->Print("spectrum_individual.pdf");
    std::cout << "Individual Done!" << std::endl;

    // groupaverage plots
    handle = fopen((path + "control_group_average_activiry.txt").c_str(), "r");
    if (handle == NULL)
    {
        std::cerr << "Cant open this file: " << (path + "control_group_average_activiry.txt") << std::endl;
        return;
    }
    for (size_t i = 0; i < the_choosen_len; i++)
    {
        fscanf(handle, "%lf\n", &control.data()[i]);
    }
    fclose(handle);

    powers = get_powerspec(control.data(), the_choosen_len);
    normalize_double(powers, the_choosen_len);
    alignment = gsl_stats_max(powers, 1, the_choosen_len) * 1.1;
    // Set the proper frequency for the powerspec.
    x = new double[the_choosen_len];
    for (size_t i = 0; i < the_choosen_len; i++)
    {
        x[i] = the_choosen_len / (4.0) / (i + 1.0);
    }

    waveform = new double[96];
    waveformx = new double[96];
    for (size_t i = 0; i < 96; i++)
    {
        waveformx[i] = i;
    }
    control_plot.reset();
    control_plot = std::make_shared<TGraph>(the_choosen_len, x, powers);
    control_plot->SetName("control_plot");
    control_plot->SetLineColor(kGreen);
    get_waveform(control.data(), the_choosen_len, waveform, 96);
    control_waveform.reset();
    control_waveform = std::make_unique<TGraph>(96, waveformx, waveform);
    control_waveform->SetLineColor(kGreen);
    control_waveform->SetTitle("");
    delete[] powers;

    handle = fopen((path + "mutant_group_average_activiry.txt").c_str(), "r");
    if (handle == NULL)
    {
        std::cerr << "Cant open this file: " << (path + "mutant_group_average_activiry.txt") << std::endl;
        return;
    }
    for (size_t i = 0; i < the_choosen_len; i++)
    {
        fscanf(handle, "%lf\n", &mutant.data()[i]);
    }
    fclose(handle);

    powers = get_powerspec(mutant.data(), the_choosen_len);
    normalize_double(powers, the_choosen_len);
    for (size_t i = 0; i < the_choosen_len; i++)
    {
        powers[i] += alignment;
    }
    alignment += (gsl_stats_max(powers, 1, the_choosen_len) - alignment) * 1.1;
    mutant_plot.reset();
    mutant_plot = std::make_shared<TGraph>(the_choosen_len, x, powers);
    mutant_plot->SetName("mutant_plot");
    mutant_plot->SetLineColor(kBlue);
    get_waveform(mutant.data(), the_choosen_len, waveform, 96);
    mutant_waveform.reset();
    mutant_waveform = std::make_unique<TGraph>(96, waveformx, waveform);
    mutant_waveform->SetLineColor(kBlue);
    mutant_waveform->SetTitle("");
    delete[] powers;

    handle = fopen((path + "mutant_hyper_group_average_activiry.txt").c_str(), "r");
    if (handle == NULL)
    {
        std::cerr << "Cant open this file: " << (path + "mutant_hyper_group_average_activiry.txt") << std::endl;
        return;
    }
    for (size_t i = 0; i < the_choosen_len; i++)
    {
        fscanf(handle, "%lf\n", &hyper.data()[i]);
    }
    fclose(handle);
    powers = get_powerspec(hyper.data(), the_choosen_len);
    normalize_double(powers, the_choosen_len);
    for (size_t i = 0; i < the_choosen_len; i++)
    {
        powers[i] += alignment;
    }
    alignment = 0;

    hyper_plot.reset();
    hyper_plot = std::make_shared<TGraph>(the_choosen_len, x, powers);
    hyper_plot->SetName("hyper_plot");
    hyper_plot->SetLineColor(kRed);
    get_waveform(hyper.data(), the_choosen_len, waveform, 96);
    hyper_waveform.reset();
    hyper_waveform = std::make_unique<TGraph>(96, waveformx, waveform);
    hyper_waveform->SetLineColor(kRed);
    hyper_waveform->SetTitle("");
    delete[] powers;

    tempx = mutant_waveform->GetX();
    for (size_t i = 0; i < 96; i++)
    {
        tempx[i] /= 4.0;
    }
    tempx = control_waveform->GetX();
    for (size_t i = 0; i < 96; i++)
    {
        tempx[i] /= 4.0;
    }
    tempx = hyper_waveform->GetX();
    for (size_t i = 0; i < 96; i++)
    {
        tempx[i] /= 4.0;
    }

    std::cout << "Prepared" << std::endl;
    canvas.reset();
    canvas = std::make_unique<TCanvas>("canvas", "canvas", 160, 90);
    canvas->SetFillStyle(0);
    canvas->SetFrameFillStyle(0);
    canvas->Divide(2, 1);
    canvas->cd(1)->SetPad("padl", "padl", 0, 0, 0.7, 1);
    canvas->cd(2)->SetPad("padr", "padr", 0.7, 0, 1, 1);

    padl = (TPad *)gROOT->FindObject("padl");
    padr = (TPad *)gROOT->FindObject("padr");
    padl->SetFillStyle(4000);
    padl->SetFrameFillStyle(4000);
    padr->SetFillStyle(4000);
    padr->SetFrameFillStyle(4000);

    canvas->cd(1);
    padl->Clear();
    padl->cd()->SetMargin(0.1, 0.01, 0.1, 0.1);
    padl->cd()->SetLogx(1);

    multigraph.reset();
    multigraph = std::make_unique<TMultiGraph>();
    multigraph->Add(control_plot.get());
    multigraph->Add(mutant_plot.get());
    multigraph->Add(hyper_plot.get());

    multigraph->GetXaxis()->SetRangeUser(0.5, 36);
    multigraph->GetXaxis()->SetTitle("Time scale (h)");
    multigraph->GetXaxis()->CenterTitle();
    multigraph->GetYaxis()->SetTitle("Power spectrum");
    multigraph->GetYaxis()->CenterTitle();
    multigraph->Draw("AL");
    legend->Draw();

    canvas->cd(2);
    padr->Clear();
    // padr->SetMargin(0.1, 0.01, 0.1, 0.1);
    padr->Divide(1, 3, 0, 0);

    padr->cd(3)->SetFillStyle(4000);
    padr->cd(3)->SetFrameFillStyle(4000);
    padr->cd(3)->SetMargin(0.01, 0.01, 0.25, 0.0);
    control_waveform->GetXaxis()->SetLabelSize(0.06);

    // control_waveform->GetXaxis()->SetTitle("Time (h)");
    // control_waveform->GetXaxis()->CenterTitle();
    control_waveform->Draw("AL");
    padr->cd(2)->SetFillStyle(4000);
    padr->cd(2)->SetFrameFillStyle(4000);
    padr->cd(2)->SetMargin(0.01, 0.01, 0.0, 0.0);
    mutant_waveform->GetXaxis()->SetLabelSize(0);

    mutant_waveform->Draw("AL");
    padr->cd(1)->SetFillStyle(4000);
    padr->cd(1)->SetFrameFillStyle(4000);
    padr->cd(1)->SetMargin(0.01, 0.01, 0.0, 0.2);
    hyper_waveform->GetXaxis()->SetLabelSize(0);

    hyper_waveform->Draw("AL");
    std::cout << "Prepared" << std::endl;

    canvas->Update();
    canvas->Print("spectrum_group_average.pdf");
    std::cout << "Group Ave Done!" << std::endl;
}

void normalize_double(double *data, int size)
{
    double mean = gsl_stats_mean(data, 1, size);
    double sd = gsl_stats_sd(data, 1, size);
    for (size_t i = 0; i < size; i++)
    {
        data[i] -= mean;
        data[i] /= sd;
    }
}

void get_waveform(double *data, int size, double *waveform, int len)
{
    int lod = size / len;
    for (size_t i = 0; i < len; i++)
    {
        waveform[i] = 0;
    }
    for (size_t i = 0; i < lod; i++)
    {
        for (size_t j = 0; j < len; j++)
        {
            waveform[j] += data[i * len + j];
        }
    }
}

void draw_overview_reduced()
{
    std::string path = "/home/yiwen/mice_physiology/micenew/rebuild/group_average/";
    int the_choosen_len = 2790;
    // individual
    std::string control_mouse = "57";
    std::string mutant_mouse = "52M";
    std::string mutant_hyper_mouse = "51M";

    std::vector<double> control(the_choosen_len);
    std::vector<double> mutant(the_choosen_len);
    std::vector<double> hyper(the_choosen_len);

    FILE *handle = fopen((path + control_mouse + "_activity_reduce.txt").c_str(), "r");
    if (handle == NULL)
    {
        std::cerr << "Cant open this file: " << path + control_mouse + "_activity_reduce.txt" << std::endl;
        return;
    }
    for (size_t i = 0; i < the_choosen_len; i++)
    {
        fscanf(handle, "%lf\n", &control.data()[i]);
    }
    fclose(handle);

    handle = fopen((path + mutant_mouse + "_activity_reduce.txt").c_str(), "r");
    if (handle == NULL)
    {
        std::cerr << "Cant open this file: " << path + mutant_mouse + "_activity_reduce.txt" << std::endl;
        return;
    }
    for (size_t i = 0; i < the_choosen_len; i++)
    {
        fscanf(handle, "%lf\n", &mutant.data()[i]);
    }
    fclose(handle);

    handle = fopen((path + mutant_hyper_mouse + "_activity_reduce.txt").c_str(), "r");
    if (handle == NULL)
    {
        std::cerr << "Cant open this file: " << path + mutant_hyper_mouse + "_activity_reduce.txt" << std::endl;
        return;
    }
    for (size_t i = 0; i < the_choosen_len; i++)
    {
        fscanf(handle, "%lf\n", &hyper.data()[i]);
    }
    fclose(handle);

    auto canvas = std::make_unique<TCanvas>("canvas", "canvas", 140, 90);
    canvas->SetFillStyle(0);
    canvas->SetFrameFillStyle(0);
    canvas->Divide(3, 1);

    canvas->cd(1)->SetPad("padb", "padb", 0, 0, 1, 0.4);
    canvas->cd(2)->SetPad("padm", "padm", 0, 0.4, 1, 0.7);
    canvas->cd(3)->SetPad("padt", "padt", 0, 0.7, 1, 1);

    auto padt = (TPad *)gROOT->FindObject("padt");
    auto padm = (TPad *)gROOT->FindObject("padm");
    auto padb = (TPad *)gROOT->FindObject("padb");

    auto x = new double[the_choosen_len];
    for (size_t i = 0; i < the_choosen_len; i++)
    {
        x[i] = i / 4.0 / 24.0;
    }
    auto control_plot = std::make_shared<TGraph>(the_choosen_len, x, control.data());
    control_plot->SetName("control_plot");
    control_plot->SetTitle("");
    control_plot->SetLineColor(kGreen);
    auto mutant_plot = std::make_shared<TGraph>(the_choosen_len, x, mutant.data());
    mutant_plot->SetName("mutant_plot");
    mutant_plot->SetTitle("");
    mutant_plot->SetLineColor(kBlue);
    auto hyper_plot = std::make_shared<TGraph>(the_choosen_len, x, hyper.data());
    hyper_plot->SetName("hyper_plot");
    hyper_plot->SetTitle("");
    hyper_plot->SetLineColor(kRed);

    canvas->cd(3);
    padt->cd();
    padt->SetFillStyle(4000);
    padt->SetFrameFillStyle(4000);
    padt->SetMargin(0.1, 0.1, 0., 0.1);
    control_plot->GetXaxis()->SetRangeUser(0, 29);
    control_plot->GetYaxis()->SetTickLength(0.01);
    control_plot->GetXaxis()->SetTickLength(0.05);
    // control_plot->GetYaxis()->SetLabelSize(0);
    control_plot->GetXaxis()->SetLabelSize(0);
    control_plot->GetYaxis()->SetTitle("Activity (h)");
    control_plot->GetYaxis()->CenterTitle();
    control_plot->Draw("AL");

    canvas->cd(2);
    padm->cd();
    padm->SetFillStyle(4000);
    padm->SetFrameFillStyle(4000);
    padm->SetMargin(0.1, 0.1, 0., 0.);
    mutant_plot->GetXaxis()->SetRangeUser(0, 29);
    mutant_plot->GetYaxis()->SetTickLength(0.01);
    mutant_plot->GetXaxis()->SetTickLength(0.05);
    // mutant_plot->GetYaxis()->SetLabelSize(0);
    mutant_plot->GetXaxis()->SetLabelSize(0);

    mutant_plot->GetYaxis()->SetTitle("Activity (h)");
    mutant_plot->GetYaxis()->CenterTitle();
    mutant_plot->Draw("AL");

    canvas->cd(1);
    padb->cd();
    padb->SetFillStyle(4000);
    padb->SetFrameFillStyle(4000);
    padb->SetMargin(0.1, 0.1, 0.3, 0.);
    hyper_plot->GetXaxis()->SetRangeUser(0, 29);
    hyper_plot->GetYaxis()->SetTickLength(0.01);
    hyper_plot->GetXaxis()->SetTickLength(0.05);
    // hyper_plot->GetYaxis()->SetLabelSize(0);
    hyper_plot->GetXaxis()->SetTitle("Time (h)");
    hyper_plot->GetXaxis()->CenterTitle();
    hyper_plot->GetYaxis()->SetTitle("Activity (h)");
    hyper_plot->GetYaxis()->CenterTitle();
    hyper_plot->Draw("AL");

    auto legend = std::make_unique<TLegend>(0.15, 0.4, 0.45, 0.8);
    legend->SetName("legend");
    legend->AddEntry(control_plot.get(), "control", "l");
    legend->AddEntry(mutant_plot.get(), "mutant", "l");
    legend->AddEntry(hyper_plot.get(), "hyper", "l");
    legend->SetBorderSize(0);
    padt->cd();
    legend->Draw();

    canvas->Update();
    canvas->Print("individual_activity_overview.pdf");
}

void draw_overview_circadian()
{
    std::string path = "/home/yiwen/mice_physiology/micenew/rebuild/group_average/";
    int the_choosen_len = 216;
    // individual
    std::string control_mouse = "57";
    std::string mutant_mouse = "52M";
    std::string mutant_hyper_mouse = "51M";

    std::vector<double> control(the_choosen_len);
    std::vector<double> mutant(the_choosen_len);
    std::vector<double> hyper(the_choosen_len);

    FILE *handle = fopen((path + control_mouse + "_circadian.txt").c_str(), "r");
    if (handle == NULL)
    {
        std::cerr << "Cant open this file: " << path + control_mouse + "_circadian.txt" << std::endl;
        return;
    }
    for (size_t i = 0; i < the_choosen_len; i++)
    {
        fscanf(handle, "%*lf\t%lf\n", &control.data()[i]);
    }
    fclose(handle);

    handle = fopen((path + mutant_mouse + "_circadian.txt").c_str(), "r");
    if (handle == NULL)
    {
        std::cerr << "Cant open this file: " << path + mutant_mouse + "_circadian.txt" << std::endl;
        return;
    }
    for (size_t i = 0; i < the_choosen_len; i++)
    {
        fscanf(handle, "%*lf\t%lf\n", &mutant.data()[i]);
    }
    fclose(handle);

    handle = fopen((path + mutant_hyper_mouse + "_circadian.txt").c_str(), "r");
    if (handle == NULL)
    {
        std::cerr << "Cant open this file: " << path + mutant_hyper_mouse + "_circadian.txt" << std::endl;
        return;
    }
    for (size_t i = 0; i < the_choosen_len; i++)
    {
        fscanf(handle, "%*lf\t%lf\n", &hyper.data()[i]);
    }
    fclose(handle);

    float scale = 0.9;

    auto canvas = std::make_unique<TCanvas>("canvas", "canvas", 140 * scale, 90 * scale);
    canvas->SetFillStyle(0);
    canvas->SetFrameFillStyle(0);
    canvas->Divide(3, 1);

    canvas->cd(1)->SetPad("padb", "padb", 0, 0, 1, 0.4);
    canvas->cd(2)->SetPad("padm", "padm", 0, 0.4, 1, 0.7);
    canvas->cd(3)->SetPad("padt", "padt", 0, 0.7, 1, 1);

    auto padt = (TPad *)gROOT->FindObject("padt");
    auto padm = (TPad *)gROOT->FindObject("padm");
    auto padb = (TPad *)gROOT->FindObject("padb");

    auto x = new double[the_choosen_len];
    for (size_t i = 0; i < the_choosen_len; i++)
    {
        x[i] = i / 6.0;
    }
    auto control_plot = std::make_shared<TGraph>(the_choosen_len, x, control.data());
    control_plot->SetName("control_plot");
    control_plot->SetTitle("");
    control_plot->SetLineColor(kGreen);
    auto mutant_plot = std::make_shared<TGraph>(the_choosen_len, x, mutant.data());
    mutant_plot->SetName("mutant_plot");
    mutant_plot->SetTitle("");
    mutant_plot->SetLineColor(kBlue);
    auto hyper_plot = std::make_shared<TGraph>(the_choosen_len, x, hyper.data());
    hyper_plot->SetName("hyper_plot");
    hyper_plot->SetTitle("");
    hyper_plot->SetLineColor(kRed);

    canvas->cd(3);
    padt->cd();
    padt->SetFillStyle(4000);
    padt->SetFrameFillStyle(4000);
    padt->SetMargin(0.1, 0.1, 0., 0.1);
    control_plot->GetXaxis()->SetRangeUser(0, 29);
    control_plot->GetYaxis()->SetTickLength(0.01);
    control_plot->GetXaxis()->SetTickLength(0.05);
    // control_plot->GetYaxis()->SetLabelSize(0);
    control_plot->GetXaxis()->SetLabelSize(0);
    control_plot->GetXaxis()->SetTitle("Time (day)");
    control_plot->GetXaxis()->CenterTitle();
    control_plot->GetYaxis()->SetTitle("Circadian Rhythm (h)");
    control_plot->GetYaxis()->CenterTitle();
    control_plot->Draw("AL");

    canvas->cd(2);
    padm->cd();
    padm->SetFillStyle(4000);
    padm->SetFrameFillStyle(4000);
    padm->SetMargin(0.1, 0.1, 0., 0.);
    mutant_plot->GetXaxis()->SetRangeUser(0, 29);
    mutant_plot->GetYaxis()->SetTickLength(0.01);
    mutant_plot->GetXaxis()->SetTickLength(0.05);
    // mutant_plot->GetYaxis()->SetLabelSize(0);
    mutant_plot->GetXaxis()->SetLabelSize(0);
    mutant_plot->GetXaxis()->SetTitle("Time (h)");
    mutant_plot->GetXaxis()->CenterTitle();
    mutant_plot->GetYaxis()->SetTitle("Circadian Rhythm (h)");
    mutant_plot->GetYaxis()->CenterTitle();
    mutant_plot->Draw("AL");

    // gStyle->SetLabelSize(.1, "XY");
    // gStyle->SetTitleSize(.1, "XY");
    // gStyle->SetTitleFontSize(.1);

    canvas->cd(1);
    padb->cd();
    padb->SetFillStyle(4000);
    padb->SetFrameFillStyle(4000);
    padb->SetMargin(0.1, 0.1, 0.3, 0.);
    hyper_plot->GetXaxis()->SetRangeUser(0, 29);
    hyper_plot->GetYaxis()->SetTickLength(0.01);
    hyper_plot->GetXaxis()->SetTickLength(0.05);
    // hyper_plot->GetYaxis()->SetLabelSize(0);
    hyper_plot->GetXaxis()->SetTitle("Time (h)");
    hyper_plot->GetXaxis()->CenterTitle();
    hyper_plot->GetYaxis()->SetTitle("Circadian Rhythm (h)");
    hyper_plot->GetYaxis()->CenterTitle();
    // hyper_plot->GetXaxis()->SetLabelSize(0.1);
    // hyper_plot->GetXaxis()->SetTitleSize(0.1);
    // hyper_plot->GetYaxis()->SetLabelSize(0.07);
    // hyper_plot->GetYaxis()->SetTitleSize(0.07);
    hyper_plot->Draw("AL");

    auto legend = std::make_unique<TLegend>(0.67, 0.55, 0.92, 0.96);
    legend->SetName("legend");
    legend->AddEntry(control_plot.get(), "control", "l");
    legend->AddEntry(mutant_plot.get(), "mutant - normal activity", "l");
    legend->AddEntry(hyper_plot.get(), "hyper - hyper activity", "l");
    // legend->SetBorderSize(0);
    padt->cd();
    legend->Draw();

    canvas->Update();
    canvas->Print("individual_circadian_overview.pdf");
}

void draw_hist_bar()
{
    const Int_t nx = 3;

    const char *labels[] = {"", "control night",
                            "", "mutant night",
                            "", "hyper night"};
    auto canvas = std::make_unique<TCanvas>("canvas");
    // float scale = 0.4;
    // canvas->SetCanvasSize(400 * scale, 300 * scale);

    const double day_data[] = {1.03417599145397, 1.25890757696252, 1.64363785908487};
    const double day_error[] = {0.226751580570385 / 9.0, 0.563334663324478 / 6.0, 0.709046560898137 / 6.0};
    const double night_data[] = {4.18429687091677, 8.05761237283794, 25.2633722831102};
    const double night_error[] = {1.48140031926948 / 9.0, 1.14716264352397 / 6.0, 11.2908820923777 / 6.0};

    auto hist = std::make_unique<TH1F>("h", "", nx, 0, nx);

    hist->SetStats(0);
    hist->SetFillColor(4);
    hist->SetBarWidth(0.4);
    hist->SetBarOffset(0.1);
    for (size_t i = 0; i < nx; i++)
    {
        hist->SetBinContent(i + 1, day_data[i]);
        hist->SetBinError(i + 1, day_error[i]);
        hist->GetXaxis()->SetBinLabel(i + 1, labels[i * 2]);
    }

    auto hist2 = std::make_unique<TH1F>("h", "", nx, 0, nx);
    hist2->SetStats(0);
    hist2->SetFillColor(35);
    hist2->SetBarWidth(0.4);
    hist2->SetBarOffset(0.5);
    hist2->SetMinimum(0);
    for (size_t i = 0; i < nx; i++)
    {
        hist2->SetBinContent(i + 1, night_data[i]);
        hist2->SetBinError(i + 1, night_error[i]);
        hist2->GetXaxis()->SetBinLabel(i + 1, labels[i * 2]);
        // hist2->GetXaxis()->SetBinLabel(i*2+1, labels[i*2]);
    }

    auto legend = std::make_unique<TLegend>(0.2, 0.7, 0.45, 0.8);
    legend->SetBorderSize(0);
    // auto multi = std::make_unique<TMultiGraph>();

    canvas->cd();
    gStyle->SetErrorX(0);

    hist2->Draw("bar e1");
    hist->Draw("bar e1 same");

    legend->AddEntry(hist.get(), "day activity", "f");
    legend->AddEntry(hist2.get(), "night activity", "f");
    legend->Draw();

    auto latex = std::make_shared<TLatex>();
    latex->SetTextFont(42);
    latex->SetTextSize(0.045);
    latex->DrawLatex(0.4, -2, "control");
    latex->DrawLatex(1.4, -2, "mutant");
    latex->DrawLatex(2.4, -2, "hyper");

    canvas->Print("hist_bar_act.pdf");
}

void draw_hist_bar_temperature()
{
    const Int_t nx = 3;

    const char *labels[] = {"", "control night",
                            "", "mutant night",
                            "", "hyper night"};
    auto canvas = std::make_unique<TCanvas>("canvas");
    // float scale = 0.4;
    // canvas->SetCanvasSize(400 * scale, 300 * scale);

    const double day_data[] = {35.3427490522401, 34.9104196906129, 34.9626726710689};
    const double day_error[] = {0.226751580570385 / 9, 0.563334663324478 / 6, 0.709046560898137 / 6};
    const double night_data[] = {37.1248049807406, 36.6844434806201, 37.0944549264285};
    const double night_error[] = {0.665320842285663 / 9, 0.575850773424495 / 6, 0.460725058332665 / 6};

    auto hist = std::make_unique<TH1F>("h", "", nx, 0, nx);

    hist->SetStats(0);
    hist->SetFillColor(4);
    hist->SetBarWidth(0.4);
    hist->SetBarOffset(0.1);
    for (size_t i = 0; i < nx; i++)
    {
        hist->SetBinContent(i + 1, day_data[i]);
        hist->SetBinError(i + 1, day_error[i]);
        hist->GetXaxis()->SetBinLabel(i + 1, labels[i * 2]);
    }

    auto hist2 = std::make_unique<TH1F>("h", "", nx, 0, nx);
    hist2->SetStats(0);
    hist2->SetFillColor(35);
    hist2->SetBarWidth(0.4);
    hist2->SetBarOffset(0.5);
    hist2->SetMinimum(34);
    for (size_t i = 0; i < nx; i++)
    {
        hist2->SetBinContent(i + 1, night_data[i]);
        hist2->SetBinError(i + 1, night_error[i]);
        hist2->GetXaxis()->SetBinLabel(i + 1, labels[i * 2]);
        // hist2->GetXaxis()->SetBinLabel(i*2+1, labels[i*2]);
    }

    auto legend = std::make_unique<TLegend>(0.2, 0.75, 0.45, 0.85);
    legend->SetBorderSize(0);
    // auto multi = std::make_unique<TMultiGraph>();

    canvas->cd();
    gStyle->SetErrorX(0);

    hist2->Draw("bar e1");
    hist->Draw("bar e1 same");

    legend->AddEntry(hist.get(), "day temperature", "f");
    legend->AddEntry(hist2.get(), "night temperature", "f");
    legend->Draw();

    auto latex = std::make_shared<TLatex>();
    latex->SetTextFont(42);
    latex->SetTextSize(0.045);
    latex->DrawLatex(0.4, 33.75, "control");
    latex->DrawLatex(1.4, 33.75, "mutant");
    latex->DrawLatex(2.4, 33.75, "hyper");

    canvas->Print("hist_bar_temp.pdf");
}

void dfa_circadian_group_average_plot(int dfa_order)
{
    std::string suffix, midfix;
    std::string type;
    // Draw the averaged dfa data into plots.
    std::string path = "/home/yiwen/mice_physiology/micenew/rebuild/group_average/";
    int the_choosen_len = 20;
    std::vector<double> dfax(the_choosen_len);
    std::vector<double> dfay(the_choosen_len);
    std::vector<double> xp(the_choosen_len);
    std::vector<double> yp(the_choosen_len);

    auto canvas = std::make_unique<TCanvas>("canvas");
    double scale = 0.7;
    canvas->SetTitle("Hello World");
    gStyle->SetOptStat(0);
    canvas->SetCanvasSize(450 * scale, 500 * scale);

    double lim = 48;
    double lim2 = 96;
    double limit = log10(lim);

    double *dfahour = nullptr;

    suffix = "ori";
    midfix = "dfa_circadian";
    type = "control";
    FILE *temp_data = fopen((path + type + "_group_average_" + midfix + "_data_" + suffix + ".txt").c_str(), "r");
    if (temp_data == NULL)
    {
        std::cerr << "Cant open file: " << path + type + "_group_average_" + midfix + "_data_" + suffix + ".txt" << std::endl;
        return;
    }
    double tempx, tempy;
    for (int i = 0; i < the_choosen_len; ++i)
    {
        fscanf(temp_data, "%lf\t%lf\n", &tempx, &tempy);
        dfax[i] = tempx;
        dfay[i] = tempy;
        xp[i] = log10(dfax[i]);
        yp[i] = log10(dfay[i]);
    }
    fclose(temp_data);
    auto graph_control = std::make_shared<TGraph>(the_choosen_len, dfax.data(), dfay.data());

    type = "mutant";
    temp_data = fopen((path + type + "_group_average_" + midfix + "_data_" + suffix + ".txt").c_str(), "r");
    if (temp_data == NULL)
    {
        std::cerr << "Cant open file: " << path + type + "_group_average_" + midfix + "_data_" + suffix + ".txt" << std::endl;
        return;
    }
    for (int i = 0; i < the_choosen_len; ++i)
    {
        fscanf(temp_data, "%lf\t%lf\n", &tempx, &tempy);
        dfax[i] = tempx;
        dfay[i] = tempy;
        xp[i] = log10(dfax[i]);
        yp[i] = log10(dfay[i]);
    }
    fclose(temp_data);
    auto graph_mutant = std::make_shared<TGraph>(the_choosen_len, dfax.data(), dfay.data());

    type = "mutant_hyper";
    temp_data = fopen((path + type + "_group_average_" + midfix + "_data_" + suffix + ".txt").c_str(), "r");
    if (temp_data == NULL)
    {
        std::cerr << "Cant open file: " << path + type + "_group_average_" + midfix + "_data_" + suffix + ".txt" << std::endl;
        return;
    }
    for (int i = 0; i < the_choosen_len; ++i)
    {
        fscanf(temp_data, "%lf\t%lf\n", &tempx, &tempy);
        dfax[i] = tempx;
        assert(tempy > 0);
        dfay[i] = tempy;
        xp[i] = log10(dfax[i]);
        yp[i] = log10(dfay[i]);
    }
    fclose(temp_data);
    auto graph_mutant_hyper = std::make_shared<TGraph>(the_choosen_len, dfax.data(), dfay.data());

    dfahour = graph_control->GetX();
    for (int i = 0; i < the_choosen_len; ++i)
    {
        dfahour[i] *= 4.0;
    }
    dfahour = graph_mutant->GetX();
    for (int i = 0; i < the_choosen_len; ++i)
    {
        dfahour[i] *= 4.0;
    }
    dfahour = graph_mutant_hyper->GetX();
    for (int i = 0; i < the_choosen_len; ++i)
    {
        dfahour[i] *= 4.0;
    }

    // graph_mutant_hyper->GetXaxis()->SetRangeUser(0, lim2/4.0);

    std::string name;

    name = "57";
    temp_data = fopen((path + name + "_" + midfix + "_data_" + suffix + ".txt").c_str(), "r");
    for (int i = 0; i < the_choosen_len; ++i)
    {
        fscanf(temp_data, "%lf\t%lf\n", &tempx, &tempy);
        dfax[i] = tempx;
        dfay[i] = tempy;
        xp[i] = log10(dfax[i]);
        yp[i] = log10(dfay[i]);
    }
    fclose(temp_data);
    auto graph_inv_control = std::make_shared<TGraph>(the_choosen_len, dfax.data(), dfay.data());

    name = "12Otx2";
    temp_data = fopen((path + name + "_" + midfix + "_data_" + suffix + ".txt").c_str(), "r");
    for (int i = 0; i < the_choosen_len; ++i)
    {
        fscanf(temp_data, "%lf\t%lf\n", &tempx, &tempy);
        dfax[i] = tempx;
        dfay[i] = tempy;
        xp[i] = log10(dfax[i]);
        yp[i] = log10(dfay[i]);
    }
    fclose(temp_data);
    auto graph_inv_mutant = std::make_shared<TGraph>(the_choosen_len, dfax.data(), dfay.data());

    name = "51M";
    temp_data = fopen((path + name + "_" + midfix + "_data_" + suffix + ".txt").c_str(), "r");
    for (int i = 0; i < the_choosen_len; ++i)
    {
        fscanf(temp_data, "%lf\t%lf\n", &tempx, &tempy);
        dfax[i] = tempx;
        dfay[i] = tempy;
        xp[i] = log10(dfax[i]);
        yp[i] = log10(dfay[i]);
    }
    fclose(temp_data);
    // double pars[5];
    // auto f1 = new TF1("f1", fit_crossover, xp[0], limit, 5);
    // auto *plot = new TGraph(the_choosen_len, xp.data(), yp.data());
    // f1->SetParLimits(0, 0, 2);
    // f1->SetParameter(0, 1.0);
    // f1->SetParLimits(2, 0, 1);
    // f1->SetParameter(2, 0.5);
    // f1->SetParLimits(4, xp[2], limit);
    // f1->SetParameter(4, xp[8]);

    // double chisq = 0;
    // for (int i = 0; i < 20; ++i)
    // {
    //     plot->Fit(f1, "R");
    //     chisq = f1->GetChisquare();
    //     if (chisq < 0.0001)
    //         break;
    // }
    // delete plot;

    // for (int i = 0; i < 5; ++i)
    // {
    //     pars[i] = f1->GetParameter(i);
    //     std::cout << "Pars" << i << ": " << pars[i] << std::endl;
    // }

    // auto f2 = new TF1("f2", draw_crossover, dfax[0] / 4.0, pow(10, pars[4]) / 4.0, 5);
    // auto f3 = new TF1("f3", draw_crossover, pow(10, pars[4]) / 4.0, lim / 4.0, 5);
    // for (int i = 0; i < 4; ++i)
    // {
    //     f2->SetParameter(i, pars[i]);
    //     f3->SetParameter(i, pars[i]);
    // }
    // f2->SetParameter(4, pow(10, pars[4]));
    // f3->SetParameter(4, pow(10, pars[4]));
    // double rvalue, lvalue;
    // lvalue = pow(10, pars[4]);
    // for (int i = 0; i < the_choosen_len; ++i)
    // {
    //     if (lvalue > dfax[i])
    //         continue;
    //     else
    //     {
    //         lvalue = dfax[i - 1];
    //         rvalue = dfax[i];
    //         break;
    //     }
    // }
    // std::cout << name << "\tDFA(ori) " << dfa_order << " Crossovers\t" << pow(10.0, pars[4]) * 0.010416667 * 24 << "\thour" << std::endl;
    // std::cout << name << "\tDFA(ori) " << dfa_order << " Alpha\t" << pars[0] << "\t" << pars[2] << std::endl;
    // auto func = std::make_shared<TF1>("func", "x=5", 5, 5);

    auto graph_inv_mutant_hyper = std::make_shared<TGraph>(the_choosen_len, dfax.data(), dfay.data());

    dfahour = graph_inv_control->GetX();
    for (int i = 0; i < the_choosen_len; ++i)
    {
        dfahour[i] *= 4.0;
    }
    dfahour = graph_inv_mutant->GetX();
    for (int i = 0; i < the_choosen_len; ++i)
    {
        dfahour[i] *= 4.0;
    }
    dfahour = graph_inv_mutant_hyper->GetX();
    for (int i = 0; i < the_choosen_len; ++i)
    {
        dfahour[i] *= 4.0;
    }

    graph_control->SetTitle(("Correlation analysis group average of circadian rhythm"));
    graph_control->SetName("graph_detrend");
    graph_control->SetMarkerStyle(72);
    graph_control->SetMarkerColor(kGreen);
    graph_control->SetMarkerSize(1);
    graph_mutant->SetName("graph_detrend");
    graph_mutant->SetMarkerStyle(71);
    graph_mutant->SetMarkerColor(kBlue);
    graph_mutant->SetMarkerSize(1);
    graph_mutant_hyper->SetName("graph_detrend");
    graph_mutant_hyper->SetMarkerStyle(73);
    graph_mutant_hyper->SetMarkerColor(kRed);
    graph_mutant_hyper->SetMarkerSize(1);

    graph_inv_control->SetTitle(("Correlation analysis individial of circadian rhythm"));
    graph_inv_control->SetName("graph_detrend");
    graph_inv_control->SetMarkerStyle(72);
    graph_inv_control->SetMarkerColor(kGreen);
    graph_inv_control->SetMarkerSize(1);
    graph_inv_mutant->SetName("graph_detrend");
    graph_inv_mutant->SetMarkerStyle(71);
    graph_inv_mutant->SetMarkerColor(kBlue);
    graph_inv_mutant->SetMarkerSize(1);
    graph_inv_mutant_hyper->SetName("graph_detrend");
    graph_inv_mutant_hyper->SetMarkerStyle(73);
    graph_inv_mutant_hyper->SetMarkerColor(kRed);
    graph_inv_mutant_hyper->SetMarkerSize(1);

    canvas->Divide(1, 2, 0, 0);

    canvas->cd(1)->SetLogx(1);
    canvas->cd(1)->SetLogy(1);

    graph_inv_control->GetYaxis()->SetTitle("F(n)");
    graph_inv_control->GetYaxis()->CenterTitle();
    graph_inv_control->GetXaxis()->SetTitle("Time(h)");
    graph_inv_control->GetXaxis()->CenterTitle();
    // graph_inv_control->Draw("AP");
    // auto l = std::make_unique<TLine>(pow(10.0, pars[4]) * 0.010416667 * 24, -3, pow(10.0, pars[4]) * 0.010416667 * 24, 100);
    // graph_inv_mutant->Draw("SAMEP");
    // graph_inv_mutant_hyper->Draw("SAMEP");

    auto gaussian_curve = std::make_shared<TF1>("gaussian", "0.03*x^0.5", 0, 1e10);
    auto pink_curve = std::make_shared<TF1>("pink", "0.015*x^1", 0, 1e10);

    gaussian_curve->SetLineColor(kBlack);
    gaussian_curve->SetLineStyle(kDashed);
    pink_curve->SetLineColor(kBlack);

    auto legend = std::make_unique<TLegend>(0.15, 0.75, 0.55, 0.95);
    legend->SetName("legend");
    legend->AddEntry(graph_control.get(), "control", "p");
    legend->AddEntry(graph_mutant.get(), "mutant - normal activity", "p");
    legend->AddEntry(graph_mutant_hyper.get(), "mutant - hyper activity", "p");
    legend->SetBorderSize(0);

    auto multigr1 = std::make_shared<TMultiGraph>();
    multigr1->Add(graph_inv_control.get());
    multigr1->Add(graph_inv_mutant.get());
    multigr1->Add(graph_inv_mutant_hyper.get());

    multigr1->Draw("AP");
    gaussian_curve->Draw("SAME");
    pink_curve->Draw("SAME");
    legend->Draw();

    // multigr1->GetYaxis()->SetTitle("F(n)");
    multigr1->GetYaxis()->CenterTitle();
    // multigr1->GetXaxis()->SetTitle("Time(h)");
    multigr1->GetXaxis()->CenterTitle();
    multigr1->GetXaxis()->SetLimits(dfax[0] * 0.9 * 4, dfax[the_choosen_len - 1] * 4 * 1.1);

    /// Draw Second
    canvas->cd(2)->SetLogx(1);
    canvas->cd(2)->SetLogy(1);
    std::cout << "Draw Second" << std::endl;
    // graph_control->GetYaxis()->SetTitle("F(n)");
    graph_control->GetYaxis()->CenterTitle();
    // graph_control->GetXaxis()->SetTitle("Time(h)");
    graph_control->GetXaxis()->CenterTitle();

    auto multigr2 = std::make_shared<TMultiGraph>();
    multigr2->Add(graph_control.get());
    multigr2->Add(graph_mutant.get());
    multigr2->Add(graph_mutant_hyper.get());

    multigr2->Draw("AP");
    gaussian_curve->Draw("SAME");
    pink_curve->Draw("SAME");

    auto legend2 = std::make_unique<TLegend>(0.15, 0.83, 0.39, 0.95);
    legend2->SetBorderSize(0);
    legend2->SetName("legend");

    legend2->AddEntry(pink_curve.get(), "alpha = 1.0", "l");
    legend2->AddEntry(gaussian_curve.get(), "alpha = 0.5", "l");
    // legend->SetBorderSize(0);
    legend2->Draw();

    // multigr2->GetYaxis()->SetTitle("F(n)");
    multigr2->GetYaxis()->CenterTitle();
    // multigr2->GetXaxis()->SetTitle("Time(h)");
    multigr2->GetXaxis()->CenterTitle();
    multigr2->GetXaxis()->SetLimits(dfax[0] * 0.9 * 4, dfax[the_choosen_len - 1] * 4 * 1.1);

    canvas->cd();

    // canvas->DrawFrame(0, 0, 1, 1);
    auto latex = std::make_shared<TLatex>();
    latex->SetTextFont(42);
    latex->SetTextSize(0.045);
    latex->SetTextAlign(11);
    latex->DrawLatex(0.44, 0, "Time(h)");

    latex->SetTextAngle(90);
    latex->SetTextAlign(13);
    latex->DrawLatex(0, 0.22, "Detrended Fluctuation Function F(n)");

    canvas->Print(("group_average_" + midfix + "_data.pdf").c_str());

    dfax.clear();
    xp.clear();
    dfay.clear();
    yp.clear();
}

void run_all_post_plot()
{

    group();

    group_average();

    dfa_group_average("dfa", "ori");
    dfa_group_average("dfa", "nowaveform");
    dfa_group_average("dfa", "detrend");
    dfa_group_average("dfamag", "ori");
    dfa_group_average("dfamag", "detrend");
    dfa_group_average("dfamag", "nowaveform");

    dfa_group_average_plot("dfa", "ori", 2, "A");
    dfa_group_average_plot("dfamag", "ori", 2, "D");

    dfa_group_average_plot("dfa", "detrend", 2, "B");
    dfa_group_average_plot("dfamag", "detrend", 2, "E");

    dfa_group_average_plot("dfa", "nowaveform", 2, "C");
    dfa_group_average_plot("dfamag", "nowaveform", 2, "F");

    // dfa_group_average_plot("control", "dfa", 2);
    // dfa_group_average_plot("mutant", "dfa", 2);
    // dfa_group_average_plot("mutant_hyper", "dfa", 2);

    // dfa_group_average_plot("control", "dfamag", 2);
    // dfa_group_average_plot("mutant", "dfamag", 2);
    // dfa_group_average_plot("mutant_hyper", "dfamag", 2);

    dfa_group_average("dfa_circadian", "ori");

    dfa_circadian_group_average_plot(2);

    powerspec_waveform_plot();

    draw_overview_reduced();

    draw_overview_circadian();

    draw_hist_bar();

    draw_hist_bar_temperature();
}