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
        inFile9Li->Close();
        inFile11Li->Close();
        return;
    }

    // Canvas para PID
    auto cPID = new TCanvas("cPID", "PID Analysis", 800, 600);
    cPID->DivideSquare(3);
    cPID->cd(1);
    hPIDLength9Li->DrawClone("colz");
    cPID->cd(2);
    hPIDLength11Li->DrawClone("colz");
    cPID->cd(3);
    hPIDLength11Li->DrawClone("colz");
    hPIDLength9Li->DrawClone("same");
    cPID->SetTitle("PID Length");

    inFile9Li->Close();
}