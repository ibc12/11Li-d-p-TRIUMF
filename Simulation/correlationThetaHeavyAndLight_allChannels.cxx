#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TString.h"
#include "TMath.h"
#include "TPaveText.h"

#include <string>
#include <vector>
#include <iostream>

void correlationThetaHeavyAndLight_allChannels()
{
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

        TTree* t3 = (TTree*)f->Get("SimulationTTree");
        TTree* t4 = (TTree*)f->Get("SimulationTTreeHeavy");
        if (!t3 || !t4) {
            std::cerr << "Missing trees in: " << fullPath << std::endl;
            f->Close();
            continue;
        }

        double theta3, theta4;
        t3->SetBranchAddress("theta3Lab", &theta3);
        t4->SetBranchAddress("theta4Lab", &theta4);

        int nEntries = std::min(t3->GetEntries(), t4->GetEntries());
        auto hist = new TH2D(("h_" + name).c_str(),
                             " ;#theta_{3}^{Lab} (deg);#theta_{4}^{Lab} (deg)",
                             180, 0, 180, 40, 0, 20);

        for (int j = 0; j < nEntries; ++j) {
            t3->GetEntry(j);
            t4->GetEntry(j);
            hist->Fill(theta3, theta4 * TMath::RadToDeg());
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

    c->SaveAs("../Figures/AllChannels_theta3_vs_theta4.png");
}