
#include <string>
#include <vector>

#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TStyle.h"

#include "ActKinematics.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/HistoModels.hxx"

#include <iostream>

void plotExpectedKinematics()
{
    // Get the output files
    std::vector<std::string> outputFiles = {
        "../Outputs/7.5MeV/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_silspecs_spacer.root",
        "../Outputs/7.5MeV/2H_1H_TRIUMF_Eex_0.130_nPS_0_pPS_0_silspecs_spacer.root",
        "../Outputs/7.5MeV/2H_1H_TRIUMF_Eex_0.435_nPS_0_pPS_0_silspecs_spacer.root",
        "../Outputs/7.5MeV/2H_1H_TRIUMF_Eex_2.000_nPS_0_pPS_0_silspecs_spacer.root",
        "../Outputs/7.5MeV/2H_1H_TRIUMF_Eex_5.000_nPS_0_pPS_0_silspecs_spacer.root",
        "../Outputs/7.5MeV/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_silspecs_spacer.root",
        "../Outputs/7.5MeV/2H_2H_TRIUMF_Eex_0.818_nPS_0_pPS_0_silspecs_spacer.root",
        "../Outputs/7.5MeV/2H_3H_TRIUMF_Eex_0.000_nPS_0_pPS_0_silspecs_spacer.root",
    };

    std::vector<std::string> reactionChannel = {
        "11Li + d -> 12Li + p (g.s. s-wave)",
        "11Li + d -> 12Li + p (Ex = 0.130 MeV p-wave)",
        "11Li + d -> 12Li + p (Ex = 0.435 MeV p-wave)",
        "11Li + d -> 12Li + p (Ex = 2.000 MeV d-wave)",
        "11Li + d -> 12Li + p (Ex = 5.000 MeV d-wave)",
        "11Li + d -> 12Li + d (elastic scattering)",
        "11Li + d -> 12Li + d (Ex = 0.818 MeV 2 neutron phase space)",
        "11Li + d -> 12Li + t (g.s.)",
    };

    std::vector<ActPhysics::Kinematics> kinematics = {
        ActPhysics::Kinematics("11Li", "d", "p", 7.5 * 11, 0.0),
        ActPhysics::Kinematics("11Li", "d", "p", 7.5 * 11, 0.130),
        ActPhysics::Kinematics("11Li", "d", "p", 7.5 * 11, 0.435),
        ActPhysics::Kinematics("11Li", "d", "p", 7.5 * 11, 2.0),
        ActPhysics::Kinematics("11Li", "d", "p", 7.5 * 11, 5.0),
        ActPhysics::Kinematics("11Li", "d", "d", 7.5 * 11, 0.0),
        ActPhysics::Kinematics("11Li", "d", "d", 7.5 * 11, 0.818),
        ActPhysics::Kinematics("11Li", "d", "t", 7.5 * 11, 0.0),
    };

    gStyle->SetOptStat(0);

    std::vector<TH2D*> histosKin;
    std::vector<TH2D*> histosAngles;

    for (size_t i = 0; i < outputFiles.size(); ++i)
    {
        TFile *file = TFile::Open(outputFiles[i].c_str(), "READ");
        if (!file || file->IsZombie())
        {
            std::cerr << "Error opening file: " << outputFiles[i] << std::endl;
            continue;
        }

        TTree *tree = (TTree *)file->Get("SimulationTTree");
        TTree *treeHeavy = (TTree *)file->Get("SimulationTTreeHeavy");

        if (!tree || !treeHeavy)
        {
            std::cerr << "Missing TTree(s) in file: " << outputFiles[i] << std::endl;
            file->Close();
            continue;
        }

        double theta3Lab, EVertex, theta4Lab;
        tree->SetBranchAddress("theta3Lab", &theta3Lab);
        tree->SetBranchAddress("EVertex", &EVertex);
        treeHeavy->SetBranchAddress("theta4Lab", &theta4Lab);

        Long64_t nEntries = tree->GetEntries();
        if (nEntries != treeHeavy->GetEntries())
        {
            std::cerr << "Mismatch in entries between trees in file: " << outputFiles[i] << std::endl;
            file->Close();
            continue;
        }

        // Histogram models
        TH2D *hKin = new TH2D("hKin", "Kinematics;#theta_{light, Lab} [#circ];E_{light} [MeV]", 350, 0, 165, 350, 0, 60);
        TH2D *hTheta3Theta4 = new TH2D("hTheta3Theta4",
                                       "Theta3Lab vs Theta4Lab; #theta_{Light, Lab} [deg]; #theta_{Heavy, Lab} [deg]",
                                       180, 0, 180, 180, 0, 20);
        hKin->SetDirectory(0);

        for (Long64_t j = 0; j < nEntries; ++j)
        {
            tree->GetEntry(j);
            treeHeavy->GetEntry(j);

            hKin->Fill(theta3Lab, EVertex);
            hTheta3Theta4->Fill(theta3Lab, theta4Lab * TMath::RadToDeg());
        }

        // Store histograms in vectors
        histosKin.push_back(hKin);
        histosAngles.push_back(hTheta3Theta4);

        // Canvas and draw
        TCanvas *c = new TCanvas(Form("c%d", (int)i), reactionChannel[i].c_str(), 1200, 600);
        c->Divide(2, 1);

        c->cd(1);
        hKin->Draw("COLZ");
        auto g3{kinematics[i].GetKinematicLine3()};
        g3->Draw("same");
        TLatex latex;
        latex.SetNDC(true);      // NDC coordinates (normalized 0–1)
        latex.SetTextSize(0.04); // Adjust text size
        latex.DrawLatex(0.15, 0.85, reactionChannel[i].c_str());

        c->cd(2);
        auto* gthetas {kinematics[i].GetTheta3vs4Line()};
        gthetas->SetLineColor(kOrange);
        hTheta3Theta4->Draw("COLZ");
        gthetas->Draw("l");

        c->Update();
    }
    gROOT->SetSelectedPad(nullptr);

    auto cAll = new TCanvas("cAll", "All Kinematics", 1200, 600);
    cAll->Divide(2, 1);
    // Draw all in the same pad
    for (size_t i = 0; i < histosKin.size(); ++i)
    {
        std::cout<<i<<'\n';
        if (i == 0)
        {
            cAll->cd(1);
            histosKin[i]->DrawClone();
            cAll->cd(2);
            histosAngles[i]->DrawClone();
        }
        else
        {
            cAll->cd(1);
            histosKin[i]->DrawClone("SAME");
            cAll->cd(2);
            histosAngles[i]->DrawClone("SAME");
        }
    }
}