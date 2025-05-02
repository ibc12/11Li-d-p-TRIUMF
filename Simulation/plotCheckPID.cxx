#include "TFile.h"
#include "TH2D.h"
#include "TCanvas.h"

#include <iostream>


void plotCheckPID()
{
    // Abrir el archivo ROOT para 9Li
    TFile* inFile9Li = TFile::Open("../DebugOutputs/checkPID_output9Li.root", "READ");
    if (!inFile9Li || inFile9Li->IsZombie())
    {
        std::cerr << "Error al abrir el archivo ROOT" << std::endl;
        return;
    }

    // Recuperar los histogramas
    auto hkinLi9Li = (TH2D*)inFile9Li->Get("hkinLi");
    auto hPID9Li = (TH2D*)inFile9Li->Get("hPID");
    auto hPIDLength9Li = (TH2D*)inFile9Li->Get("hPIDLength");
    hPIDLength9Li->SetTitle("PID Length 9Li");

    if (!hkinLi9Li || !hPID9Li || !hPIDLength9Li)
    {
        std::cerr << "Uno o más histogramas no se pudieron cargar." << std::endl;
        inFile9Li->Close();
        return;
    }

    // Abrir el archivo ROOT para 11Li
    TFile* inFile11Li = TFile::Open("../DebugOutputs/checkPID_output11Li.root", "READ");
    if (!inFile11Li || inFile11Li->IsZombie())
    {
        std::cerr << "Error al abrir el archivo ROOT" << std::endl;
        inFile9Li->Close();
        return;
    }
    // Recuperar los histogramas
    auto hPIDLength11Li = (TH2D*)inFile11Li->Get("hPIDLength");
    auto hPID11Li = (TH2D*)inFile11Li->Get("hPID");
    auto hkinLi11Li = (TH2D*)inFile11Li->Get("hkinLi");
    if (!hkinLi11Li || !hPID11Li || !hPIDLength11Li)
    {
        std::cerr << "Uno o más histogramas no se pudieron cargar." << std::endl;
        inFile11Li->Close();
        return;
    }

    // Abrir el archivo ROOT para 1H
    TFile* inFile1H = TFile::Open("../DebugOutputs/checkPID_output1H.root", "READ");
    if (!inFile1H || inFile1H->IsZombie())
    {
        std::cerr << "Error al abrir el archivo ROOT" << std::endl;
        inFile1H->Close();
        return;
    }
    // Recuperar los histogramas
    auto hPIDLength1H = (TH2D*)inFile1H->Get("hPIDLength");
    auto hPID1H = (TH2D*)inFile1H->Get("hPID");
    auto hkinLi1H = (TH2D*)inFile1H->Get("hkinLi");
    if (!hkinLi1H || !hPID1H || !hPIDLength1H)
    {
        std::cerr << "Uno o más histogramas no se pudieron cargar." << std::endl;
        inFile1H->Close();
        return;
    }

    // Canvas para PID
    auto cPID = new TCanvas("cPID", "PID Analysis", 800, 600);
    cPID->DivideSquare(4);
    cPID->cd(1);
    hPIDLength9Li->DrawClone("colz");
    cPID->cd(2);
    hPIDLength11Li->DrawClone("colz");
    cPID->cd(3);
    hPID1H->DrawClone("colz");
    cPID->cd(4);
    hPIDLength11Li->DrawClone("colz");
    hPIDLength9Li->DrawClone("same");
    hPIDLength1H->DrawClone("same");
    cPID->SetTitle("PID Length");

    inFile9Li->Close();
}