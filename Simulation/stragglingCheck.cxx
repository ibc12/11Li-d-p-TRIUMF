#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLatex.h"
#include "TRandom3.h"

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

    auto hEloss1 {new TH1D("hEloss", "eLoss Tini = 80 MeV", 100, 0, 20)};
    auto hEloss2 {new TH1D("hEloss2", "eLoss Tini = 50 MeV", 100, 0, 20)};
    auto hEloss3 {new TH1D("hEloss3", "eLoss Tini = 30 MeV", 100, 0, 30)};
    for (int i = 0; i < 10000; i++)
    {
        double Tini1 = 80; // MeV
        double Tini2 = 50; // MeV
        double Tini3 = 30; // MeV

        auto eAfter1 {srim->SlowWithStraggling("ParticleInSil", Tini1, 0.1, 0)};
        auto eAfter2 {srim->SlowWithStraggling("ParticleInSil", Tini2, 0.1, 0)};
        auto eAfter3 {srim->SlowWithStraggling("ParticleInSil", Tini3, 0.1, 0)};

        auto eLoss1 {Tini1 - eAfter1};
        auto eLoss2 {Tini2 - eAfter2};
        auto eLoss3 {Tini3 - eAfter3};
        hEloss1->Fill(eLoss1);
        hEloss2->Fill(eLoss2);
        hEloss3->Fill(eLoss3);
    }

    // Fit the histograms to gausians
    auto f1 = new TF1("f1", "gaus", 0, 100);
    hEloss1->Fit(f1, "R");
    auto f2 = new TF1("f2", "gaus", 0, 100);
    hEloss2->Fit(f2, "R");
    auto f3 = new TF1("f3", "gaus", 0, 100);
    hEloss3->Fit(f3, "R");
    f1->SetLineColor(kRed);
    f2->SetLineColor(kBlue);
    f3->SetLineColor(kGreen);
    f1->SetTitle("Gaussian Fit 1");
    f2->SetTitle("Gaussian Fit 2");
    f3->SetTitle("Gaussian Fit 3");
    f1->SetLineWidth(2);    
    f2->SetLineWidth(2);
    f3->SetLineWidth(2);

    // Crear canvas 
    auto c {new TCanvas("new", "Straggling Analysis", 1200, 800)};
    c->DivideSquare(3);
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
    c->cd(3);
    hEloss3->DrawClone("colz");
    f3->Draw("same");
    auto mean3 = f3->GetParameter(1);
    auto sigma3 = f3->GetParameter(2);
    auto mean3Err = f3->GetParError(1);
    auto sigma3Err = f3->GetParError(2);    
    auto mean3Text = new TLatex(1, 500, Form("Mean: %.2f #pm %.2f", mean3, mean3Err));
    mean3Text->SetTextColor(kGreen);
    mean3Text->Draw();
    auto sigma3Text = new TLatex(1, 400, Form("Sigma: %.2f #pm %.2f", sigma3, sigma3Err));
    sigma3Text->SetTextColor(kGreen);
    sigma3Text->Draw();
    
}