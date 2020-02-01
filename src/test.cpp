/* #include <TPDF.h> */
#include <iostream>
#include <TF1.h>
#include <TCanvas.h>
#include <TH1F.h>

int main(int argc, char *argv[])
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
    return 0;
}
