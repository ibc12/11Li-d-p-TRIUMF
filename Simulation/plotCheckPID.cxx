#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"

#include <iostream>


void plotCheckPID()
{
    std::vector<std::string> particles{"1H", "2H", "3H"};
    std::vector<TFile*> inFiles;
    std::vector<TH2D*> hkin;
    std::vector<TH2D*> hPID;
    std::vector<TH2D*> hPIDLength;

    for (const auto& p : particles)
    {
        std::string fileName = "../DebugOutputs/checkPID_output" + p + ".root";
        TFile* file = TFile::Open(fileName.c_str(), "READ");
        if (!file || file->IsZombie())
        {
            std::cerr << "Error al abrir el archivo ROOT para " << p << std::endl;
            inFiles.push_back(nullptr);
            hkin.push_back(nullptr);
            hPID.push_back(nullptr);
            hPIDLength.push_back(nullptr);
            if (file) file->Close();
            continue;
        }

        auto hkinHist = (TH2D*)file->Get("hkinLi");
        auto hPIDHist = (TH2D*)file->Get("hPID");
        auto hPIDLengthHist = (TH2D*)file->Get("hPIDLength");

        if (!hkinHist || !hPIDHist || !hPIDLengthHist)
        {
            std::cerr << "Uno o más histogramas no se pudieron cargar para " << p << std::endl;
            file->Close();
            inFiles.push_back(nullptr);
            hkin.push_back(nullptr);
            hPID.push_back(nullptr);
            hPIDLength.push_back(nullptr);
            continue;
        }

        hPIDLengthHist->SetTitle(("PID Length for " + p).c_str());

        inFiles.push_back(file);
        hkin.push_back(hkinHist);
        hPID.push_back(hPIDHist);
        hPIDLength.push_back(hPIDLengthHist);
    }

    // Crear canvas para PID
    auto cPID = new TCanvas("cPID", "PID Analysis", 1200, 800);
    int nParticles = particles.size();
    cPID->DivideSquare(nParticles + 1);  // Uno por partícula + combinado final

    // Dibujar cada hPIDLength por separado
    for (size_t i = 0; i < particles.size(); ++i)
    {
        cPID->cd(i + 1);
        if (hPIDLength[i]) {
            hPIDLength[i]->DrawClone("colz");
            gPad->SetTitle(particles[i].c_str());
        }
    }

    // Plot combinado
    cPID->cd(nParticles + 1);
    bool firstDrawn = false;
    for (size_t i = 0; i < particles.size(); ++i)
    {
        if (!hPIDLength[i]) continue;
        auto clone = (TH2D*)hPIDLength[i]->Clone();
        clone->SetTitle("PID Length: All Particles");
        clone->SetLineColor(i + 1);  // Diferente color
        clone->SetMarkerColor(i + 1);
        clone->SetLineWidth(2);
        clone->SetStats(0);
        if (!firstDrawn)
        {
            clone->Draw("colz");
            firstDrawn = true;
        }
        else
        {
            clone->Draw("colz same");
        }
    }
}