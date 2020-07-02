/* #include <TPDF.h> */
#include <iostream>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1F.h>

int main()
{
    TCanvas* canvas = new TCanvas("canvas");
    TH1F* histo = new TH1F("histo","test 1",10,0.,10.);
    histo->SetFillColor(2);
    histo->Fill(2.);
    histo->Draw();
    canvas->Print("plots.pdf(","Title:One bin filled");
    histo->Fill(4.);
    histo->Draw();
    canvas->Print("plots.pdf","Title:Two bins filled");
    histo->Fill(6.);
    histo->Draw();
    canvas->Print("plots.pdf","Title:Three bins filled");
    histo->Fill(8.);
    histo->Draw();
    canvas->Print("plots.pdf","Title:Four bins filled");
    histo->Fill(8.);
    histo->Draw();
    canvas->Print("plots.pdf)","Title:The fourth bin content is 2");
    double x[10] ={1,2,3,4,5,6,7,8,9,10};
    double y[10] ={2,3,4,5,6,7,8,9,10,11};
    TGraph *gr = new TGraph(10,x,y);
    canvas->Clear();
    gr->Draw();
    sleep(1);
    delete gr;
    delete histo;
    delete canvas;
    return 0;
}
