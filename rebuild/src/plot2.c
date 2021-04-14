{
    float r = 0.3;
    float epsilon = 0.02;

    TH1F *h1 = new TH1F("h1", "test1", 100, -3, 3);
    h1->SetStats(0);
    h1->GetXaxis()->SetLabelSize(0.);
    h1->GetXaxis()->SetTitleSize(0.);
    h1->FillRandom("gaus", 200000);
    h1->SetMinimum(0.);

    TH1F *h2 = new TH1F("h2", "test2", 100 , -3, 3);
    h2->SetStats(0);
    h2->FillRandom("gaus", 100000);

    TCanvas *c1 = new TCanvas("c1", "example", 600, 700);
    TPad *pad1 = new TPad("pad1", "pad1", 0, r-epsilon, 1, 1);
    pad1->SetBottomMargin(epsilon);
    c1->cd();
    pad1->Draw();
    pad1->cd();
    h1->Draw();

    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, r*(1-epsilon));
    pad2->SetTopMargin(0);
    pad2->SetFillColor(0);
    pad2->SetFillStyle(0);
    c1->cd();
    pad2->Draw();
    pad2->cd();
    h2->Draw("ep");
}
