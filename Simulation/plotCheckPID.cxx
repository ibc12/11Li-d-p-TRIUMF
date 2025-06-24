#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"

#include <iostream>
#include <vector>
#include <string>

#include "../Histos.h"

void plotCheckPID()
{
    bool isTelescope {true}; // Cambiar a true si es telescopio

    std::vector<std::string> particles{"7Li", "8Li", "9Li", "11Li"};
    // std::vector<std::string> particles{"1H", "2H", "3H"};
    std::vector<TFile*>       inFiles;
    std::vector<TH2D*>        hkin;
    std::vector<TH2D*>        hPIDfront;
    std::vector<TH2D*>        hPIDLengthfront;
    std::vector<TH2D*>        hPIDside;
    std::vector<TH2D*>        hPIDLengthside;

    // 1) Abrimos todos los archivos y extraemos los histogramas correspondientes
    for (const auto& p : particles)
    {
        std::string fileName;
        if (isTelescope)
        {
            fileName = "../DebugOutputs/checkPID_outputTelescope" + p + ".root";
        }
        else
        {
            fileName = "../DebugOutputs/checkPID_output" + p + ".root";
        }

        TFile* file = TFile::Open(fileName.c_str(), "READ");
        if (!file || file->IsZombie())
        {
            std::cout << "Error al abrir el archivo ROOT para " << p << std::endl;
            inFiles.push_back(nullptr);
            hkin.push_back(nullptr);
            hPIDfront.push_back(nullptr);
            hPIDLengthfront.push_back(nullptr);
            hPIDside.push_back(nullptr);
            hPIDLengthside.push_back(nullptr);
            if (file) file->Close();
            continue;
        }

        // Nombres exactos de los objetos dentro del .root
        auto hkinHist            = (TH2D*)file->Get("hkinLi");
        auto hPIDHistfront       = (TH2D*)file->Get("hPIDfront");
        auto hPIDLengthHistfront = (TH2D*)file->Get("hPIDLengthfront");
        auto hPIDHistside        = (TH2D*)file->Get("hPIDside");
        auto hPIDLengthHistside  = (TH2D*)file->Get("hPIDLengthside");

        if (!hkinHist || !hPIDHistfront || !hPIDLengthHistfront || !hPIDHistside || !hPIDLengthHistside)
        {
            std::cout << "Uno o más histogramas no se pudieron cargar para " << p << std::endl;
            file->Close();
            inFiles.push_back(nullptr);
            hkin.push_back(nullptr);
            hPIDfront.push_back(nullptr);
            hPIDLengthfront.push_back(nullptr);
            hPIDside.push_back(nullptr);
            hPIDLengthside.push_back(nullptr);
            continue;
        }

        // Ajustamos títulos de los histogramas “length”
        hPIDLengthHistfront->SetTitle(("PID Length for front " + p).c_str());
        hPIDLengthHistside->SetTitle(("PID Length for side  " + p).c_str());

        // Guardamos punteros en los vectores correspondientes
        inFiles.push_back(file);
        hkin.push_back(hkinHist);
        hPIDfront.push_back(hPIDHistfront);
        hPIDLengthfront.push_back(hPIDLengthHistfront);
        hPIDside.push_back(hPIDHistside);
        hPIDLengthside.push_back(hPIDLengthHistside);
    }

    int nParticles = particles.size();

    // 2) CANVAS para “front”
    auto cPIDfront = new TCanvas("cPID_front", "PID Analysis - Front Silicons", 1200, 800);
    cPIDfront->DivideSquare(nParticles + 1);  // Una columna por partícula + un plot combinado final

    if (isTelescope)
    {
        // Dibujar cada hPIDside (Telescopio) en su pad
        for (size_t i = 0; i < particles.size(); ++i)
        {
            cPIDfront->cd(i + 1);
            if (hPIDside[i])
            {
                hPIDfront[i]->DrawClone("colz");
                hPIDfront[i]->GetYaxis()->SetTitle("#DeltaE_{sil1}");
                hPIDfront[i]->GetXaxis()->SetTitle("#DeltaE_{sil2}");
                gPad->SetTitle(particles[i].c_str());
            }
        }

        // Plot combinado (todos los hPIDside juntos con suma de cuentas)
        cPIDfront->cd(nParticles + 1);
        auto hSum {Histos::PIDHeavy.GetHistogram()};

        for (size_t i = 0; i < particles.size(); ++i)
        {
            if (!hPIDfront[i]) continue;

            if (!hSum)
            {
                
            }
            else
            {
                hSum->Add(hPIDfront[i]);  // suma bin a bin
            }
        }
        hSum->GetYaxis()->SetTitle("#DeltaE_{sil1}");
        hSum->GetXaxis()->SetTitle("#DeltaE_{sil2}");
        hSum->SetTitle("PID for all Li isotopes");
        if (hSum)
            hSum->DrawClone("colz");
    }
    else
    {
        // Dibujar cada hPIDLengthfront (sin telescopio) en su pad
        for (size_t i = 0; i < particles.size(); ++i)
        {
            cPIDfront->cd(i + 1);
            if (hPIDLengthfront[i])
            {
                hPIDLengthfront[i]->DrawClone("colz");
                gPad->SetTitle(particles[i].c_str());
            }
        }

        // Plot combinado para “front length”
        cPIDfront->cd(nParticles + 1);
        bool firstDrawnFL = false;
        for (size_t i = 0; i < particles.size(); ++i)
        {
            if (!hPIDLengthfront[i]) continue;
            auto clone = (TH2D*)hPIDLengthfront[i]->Clone();
            clone->SetTitle("PID Length Front: All Particles");
            clone->SetLineColor(i + 1);
            clone->SetMarkerColor(i + 1);
            clone->SetLineWidth(2);
            clone->SetStats(0);

            if (!firstDrawnFL)
            {
                clone->DrawClone("colz");
                firstDrawnFL = true;
            }
            else
            {
                clone->DrawClone("colz same");
            }
        }
    }

    // 3) CANVAS para “side”
    auto cPIDside = new TCanvas("cPID_side", "PID Analysis - Side Silicons", 1200, 800);
    cPIDside->DivideSquare(nParticles + 1);

    if (isTelescope)
    {
        
        // Dibujar cada hPIDside (Telescopio) en su pad
        for (size_t i = 0; i < particles.size(); ++i)
        {
            cPIDside->cd(i + 1);
            if (hPIDside[i])
            {
                hPIDside[i]->DrawClone("colz");
                hPIDside[i]->GetYaxis()->SetTitle("#DeltaE_{sil1}");
                hPIDside[i]->GetXaxis()->SetTitle("#DeltaE_{sil2}");
                gPad->SetTitle(particles[i].c_str());
            }
        }

        // Plot combinado (todos los hPIDside juntos)
        cPIDside->cd(nParticles + 1);
        bool firstDrawnS = false;
        for (size_t i = 0; i < particles.size(); ++i)
        {
            if (!hPIDside[i]) continue;
            auto clone = (TH2D*)hPIDside[i]->Clone();
            clone->SetTitle("PID Side: All Particles");
            clone->GetYaxis()->SetTitle("#DeltaE_{sil1}");
            clone->GetXaxis()->SetTitle("#DeltaE_{sil2}");
            clone->SetLineColor(i + 1);
            clone->SetMarkerColor(i + 1);
            clone->SetLineWidth(2);
            clone->SetStats(0);

            if (!firstDrawnS)
            {
                clone->DrawClone("colz");
                firstDrawnS = true;
            }
            else
            {
                clone->DrawClone("colz same");
            }
        }
    }
    else
    {
        // Dibujar cada hPIDLengthside (sin telescopio) en su pad
        for (size_t i = 0; i < particles.size(); ++i)
        {
            cPIDside->cd(i + 1);
            if (hPIDLengthside[i])
            {
                hPIDLengthside[i]->DrawClone("colz");
                gPad->SetTitle(particles[i].c_str());
            }
        }

        // Plot combinado para “side length”
        cPIDside->cd(nParticles + 1);
        bool firstDrawnSL = false;
        for (size_t i = 0; i < particles.size(); ++i)
        {
            if (!hPIDLengthside[i]) continue;
            auto clone = (TH2D*)hPIDLengthside[i]->Clone();
            clone->SetTitle("PID Length Side: All Particles");
            clone->SetLineColor(i + 1);
            clone->SetMarkerColor(i + 1);
            clone->SetLineWidth(2);
            clone->SetStats(0);

            if (!firstDrawnSL)
            {
                clone->DrawClone("colz");
                firstDrawnSL = true;
            }
            else
            {
                clone->DrawClone("colz same");
            }
        }
    }
}
