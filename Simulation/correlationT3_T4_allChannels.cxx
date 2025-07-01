#include "ActSilSpecs.h"

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TString.h"
#include "TMath.h"
#include "TPaveText.h"

#include "Math/Point3D.h"

#include <string>
#include <vector>
#include <cmath>
#include <iostream>

bool BothHitf0(double theta3Lab, double phi3, double theta4Lab, double phi4, ROOT::Math::XYZPoint vertex, 
               ActPhysics::SilSpecs *sils)
{
    ROOT::Math::XYZVector direction3{TMath::Cos(theta3Lab * TMath::DegToRad()), TMath::Sin(theta3Lab * TMath::DegToRad()) * TMath::Sin(phi3),
                            TMath::Sin(theta3Lab * TMath::DegToRad()) * TMath::Cos(phi3)};
    ROOT::Math::XYZVector direction4{TMath::Cos(theta4Lab), TMath::Sin(theta4Lab) * TMath::Sin(phi4),
                            TMath::Sin(theta4Lab) * TMath::Cos(phi4)};       


    int silIndex3{};
    ROOT::Math::XYZPoint silPoint3{};
    int silIndex4{};
    ROOT::Math::XYZPoint silPoint4{};                 
    // Check if both particles hit the f0 layer
    // std::tie(silIndex3, silPoint3) = sils->FindSPInLayer("f0", vertex, direction3);
    std::tie(silIndex4, silPoint4) = sils->FindSPInLayer("f0", vertex, direction4);
    return (silIndex3 != -1 && silIndex4 != -1);
}

void correlationT3_T4_allChannels()
{
    // Silicons
    auto *sils{new ActPhysics::SilSpecs};
    std::string silConfig("silicons_reverse");
    sils->ReadFile("../configs/" + silConfig + ".conf");
    for (auto &[name, layer] : sils->GetLayers())
    {
        if (name == "f0" || name == "f1")
            layer.MoveZTo(75, {3});
        if (name == "f2" || name == "f3")
            layer.MoveZTo(125, {0});
        if (name == "l0" || name == "r0")
            layer.MoveZTo(75, {3});
    }
    sils->DrawGeo();


    std::string path = "../Outputs/7.5MeV/";

    std::vector<std::pair<std::string, std::string>> files = {
        {"dp1NPS", "2H_1H_TRIUMF_Eex_0.000_nPS_1_pPS_0_silicons_reverse"},
        {"dp3NPS", "2H_1H_TRIUMF_Eex_0.000_nPS_3_pPS_0_silicons_reverse"},
        {"dp", "2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_silicons_reverse"},
        {"dt", "2H_3H_TRIUMF_Eex_0.000_nPS_0_pPS_0_silicons_reverse"},
        {"elastic", "2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_silicons_reverse"},
        {"inelastic", "2H_2H_TRIUMF_Eex_1.266_nPS_0_pPS_0_silicons_reverse"},
        {"inelastic2NPS", "2H_2H_TRIUMF_Eex_0.818_nPS_2_pPS_0_silicons_reverse"}
    };

    int n = files.size();
    int nCols = 3;
    int nRows = (n + nCols - 1) / nCols;

    TCanvas* c = new TCanvas("c", "All Histos", 1200, 400 * nRows);
    c->Divide(nCols, nRows);

    for (size_t i = 0; i < files.size(); ++i) {
        const auto& [name, filename] = files[i];
        std::string fullPath = path + filename + ".root";

        TFile* f = TFile::Open(fullPath.c_str(), "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "Error opening file: " << fullPath << std::endl;
            continue;
        }

        TTree* t3 = (TTree*)f->Get("SimulationTTreeNoCuts");
        TTree* t4 = (TTree*)f->Get("SimulationTTreeHeavy");
        if (!t3 || !t4) {
            std::cerr << "Missing trees in: " << fullPath << std::endl;
            f->Close();
            continue;
        }

        double theta3, theta4, phi3, phi4, T3, T4;
        ROOT::Math::XYZPoint vertex;
        t3->SetBranchAddress("T3Lab", &T3);
        t3->SetBranchAddress("theta3Lab", &theta3);
        t3->SetBranchAddress("phi3CM", &phi3);
        t3->SetBranchAddress("T4Lab", &T4);
        t3->SetBranchAddress("theta4Lab", &theta4);
        t3->SetBranchAddress("phi4CM", &phi4);
        t3->SetBranchAddress("RP", &vertex);

        int nEntries = t3->GetEntries();
        // int nEntries = std::min(t3->GetEntries(), t4->GetEntries());
        auto hist = new TH2D(("h_" + name).c_str(),
                             " ;T_{3}^{Lab} (MeV);T_{4}^{Lab} (MeV)",
                             100, 0, 80, 100, 0, 80);

        for (int j = 0; j < nEntries; ++j) 
        {
            t3->GetEntry(j);
            t4->GetEntry(j);
            if(BothHitf0(theta3, phi3, theta4, phi4, vertex, sils))
                hist->Fill(T3, T4);
        }

        hist->GetXaxis()->SetTitleSize(0.04);
        hist->GetYaxis()->SetTitleSize(0.04);
        hist->GetXaxis()->SetLabelSize(0.025);
        hist->GetYaxis()->SetLabelSize(0.025);


        c->cd(i + 1);
        hist->DrawClone("colz");
        // gPad->SetGrid();
        // gPad->SetMargin(0.12, 0.15, 0.15, 0.10);

        hist->SetTitle("");

        // Dibujar título manual grande arriba
        TPaveText* title = new TPaveText(0.1, 0.91, 0.9, 0.98, "NDC");
        title->AddText(name.c_str());
        title->SetFillColor(0);
        title->SetTextSize(0.08); // Aquí cambias el tamaño del título
        title->SetTextAlign(22);  // Centrado
        title->Draw();
        
    }

    // c->SaveAs("../Figures/AllChannels_T3_vs_T4.png");
}