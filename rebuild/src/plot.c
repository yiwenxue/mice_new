{
    auto cv = new TCanvas("super_cv");
    cv->SetFillColor(0);
    cv->SetFrameFillStyle(0);

    auto f1 = new TF1("f1", "sin(x)", -50, 50);
    auto f3 = new TF1("f2", "cos(x)", -50, 50);
    auto f2 = new TF1("f3", "0.5*sin(x)+0.5", -50, 50);
    auto f4 = new TF1("f4", "0.5*cos(x)-0.5", -50, 50);

    f1->SetNpx(1000);
    f2->SetNpx(1000);
    f3->SetNpx(1000);
    f4->SetNpx(1000);

    auto pad = new TPad("pad", "", 0,0,1,1);
    pad->SetMargin(1e-3, 1e-3, 1e-3, 1e-3);

    pad->SetFillStyle(4000);
    pad->SetBorderMode (1);
    pad->SetBorderSize (1);
    pad->Divide(2,2, 0.01, 0.01);
    pad->Draw();

    pad->cd(1);
    f1->Draw();
    pad->cd(2);
    f2->Draw();
    pad->cd(3);
    f3->Draw();
    pad->cd(4);
    f4->Draw();


    cv->Print("Hello.pdf","Title:Hello");

    /* delete f1; */ 
    /* delete f2; */
    /* delete f3; */
    /* delete f4; */
    /* delete pad; */
    /* delete cv; */ 

}
