#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLatex.h"

#include "ActSRIM.h"

#include <iostream>


void stragglingCheck()
{
    // SRIM
    std::string particle {"11Li"};
    std::string path{"../SRIM files/"};
    std::string gas{"900mb_CF4_90-10"};
    std::string silicon{"silicon"};
    auto* srim {new ActPhysics::SRIM};
    srim->ReadTable("ParticleGas", path + particle + "_" + gas + ".txt");
    srim->ReadTable("ParticleInSil", path + particle + "_" + silicon + ".txt");

    auto hEloss1 {new TH1D("hEloss", "eLoss; Tini = 80 MeV", 100, 0, 20)};
    auto hEloss2 {new TH1D("hEloss2", "eLoss; Tini = 50 MeV", 100, 0, 20)};
    for (int i = 0; i < 10000; i++)
    {
        double Tini1 = 80; // MeV
        double Tini2 = 30; // MeV

        auto eAfter1 {srim->SlowWithStraggling("ParticleInSil", Tini1, 0.06, 0)};
        auto eAfter2 {srim->SlowWithStraggling("ParticleInSil", Tini2, 0.06, 0)};

        auto eLoss1 {Tini1 - eAfter1};
        auto eLoss2 {Tini2 - eAfter2};
        hEloss1->Fill(eLoss1);
        hEloss2->Fill(eLoss2);
    }

    // Fit the histograms to gausians
    auto f1 = new TF1("f1", "gaus", 0, 100);
    hEloss1->Fit(f1, "R");
    auto f2 = new TF1("f2", "gaus", 0, 100);
    hEloss2->Fit(f2, "R");
    f1->SetLineColor(kRed);
    f2->SetLineColor(kBlue);
    f1->SetTitle("Gaussian Fit 1");
    f2->SetTitle("Gaussian Fit 2");
    f1->SetLineWidth(2);    
    f2->SetLineWidth(2);

    // Crear canvas 
    auto c {new TCanvas("new", "Straggling Analysis", 1200, 800)};
    c->DivideSquare(2);
    // Plot the histograms, fit function and the values of fit parameters
    c->cd(1);
    hEloss1->DrawClone("colz");
    f1->Draw("same");
    auto mean1 = f1->GetParameter(1);
    auto sigma1 = f1->GetParameter(2);
    auto mean1Err = f1->GetParError(1);
    auto sigma1Err = f1->GetParError(2);
    auto mean1Text = new TLatex(10, 500, Form("Mean: %.2f #pm %.2f", mean1, mean1Err));
    mean1Text->SetTextColor(kRed);
    mean1Text->Draw();
    auto sigma1Text = new TLatex(10, 400, Form("Sigma: %.2f #pm %.2f", sigma1, sigma1Err));
    sigma1Text->SetTextColor(kRed);
    sigma1Text->Draw();
    c->cd(2);
    hEloss2->DrawClone("colz");
    f2->Draw("same");
    auto mean2 = f2->GetParameter(1);
    auto sigma2 = f2->GetParameter(2);
    auto mean2Err = f2->GetParError(1);
    auto sigma2Err = f2->GetParError(2);
    auto mean2Text = new TLatex(10, 500, Form("Mean: %.2f #pm %.2f", mean2, mean2Err));
    mean2Text->SetTextColor(kBlue);
    mean2Text->Draw();
    auto sigma2Text = new TLatex(10, 400, Form("Sigma: %.2f #pm %.2f", sigma2, sigma2Err));
    sigma2Text->SetTextColor(kBlue);
    sigma2Text->Draw();
    
}