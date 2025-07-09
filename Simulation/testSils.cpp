#ifndef test_cxx
#define test_cxx
#include "ActKinematics.h"
#include "ActSRIM.h"
#include "ActSilSpecs.h"
#include "ActTPCParameters.h"
#include "ActCrossSection.h"
#include "ActColors.h"

#include "TCanvas.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TPolyLine3D.h"

#include <cmath>
#include <iostream>
#include <string>
#include <unordered_map>

#include "../Histos.h"

using XYZPoint = ROOT::Math::XYZPoint;
using XYZVector = ROOT::Math::XYZVector;

void testSils()
{
    ActRoot::TPCParameters tpc {"Actar"};
    std::cout << "TPC: " << tpc.X() << " " << tpc.Y() << " " << tpc.Z() << '\n';
    // Silicons
    auto* sils {new ActPhysics::SilSpecs};
    sils->ReadFile("../configs/silicons_reverse.conf");
    sils->Print();
    const double sigmaSil {0.060 / 2.355}; // Si resolution
    auto silRes = std::make_unique<TF1>(
        "silRes", [=](double* x, double* p) { return sigmaSil * TMath::Sqrt(x[0] / 5.5); }, 0.0, 100.0, 1);
    std::vector<std::string> silLayers {"f0", "f1", "f2", "f3", "l0", "r0"};
    // We have to centre the silicons with the beam input
    // In real life beam window is not at Z / 2
    for (auto &[name, layer] : sils->GetLayers())
    {
        if (name == "f0" || name == "f1")
            layer.MoveZTo(75, {3});
        if (name == "f2" || name == "f3")
            layer.MoveZTo(125, {0});
        if (name == "l0" || name == "r0")
            layer.MoveZTo(75, {3});
    }

    //sils->DrawGeo();

    double thetaLab = 90 * TMath::DegToRad();
    double phiLab1 = 90 * TMath::DegToRad();
    double phiLab2 = -90 * TMath::DegToRad(); 
    XYZVector dir1 = {TMath::Cos(thetaLab), TMath::Sin(thetaLab) * TMath::Sin(phiLab1),
                             TMath::Sin(thetaLab) * TMath::Cos(phiLab1)};
    XYZVector dir2 = {TMath::Cos(thetaLab), TMath::Sin(thetaLab) * TMath::Sin(phiLab2),
                             TMath::Sin(thetaLab) * TMath::Cos(phiLab2)};
    XYZPoint vertex {128, 128, 30}; // Mitad de la cámara



    int silIndex0 = -1;
    ROOT::Math::XYZPoint silPoint0;
    std::string layer0;
    for(auto layer : silLayers)
        {
            std::tie(silIndex0, silPoint0) = sils->FindSPInLayer(layer, vertex, dir1);
            if(silIndex0 != -1)
            {
                layer0 = layer;
                break;
            }                
        } 
    int silIndex1 = -1;
    ROOT::Math::XYZPoint silPoint1;
    std::string layer1;
    for(auto layer : silLayers)
        {
            std::tie(silIndex1, silPoint1) = sils->FindSPInLayer(layer, vertex, dir2);
            if(silIndex1 != -1)
            {
                layer1 = layer;
                break;
            }                
        } 

    std::cout<<"Distance 1: "<<(silPoint0 - vertex).R()<<std::endl;
    std::cout<<"Distance 2: "<<(silPoint1 - vertex).R()<<std::endl;
    std::cout<<"Layer 1: "<<layer0<<std::endl;  
    std::cout<<"Layer 2: "<<layer1<<std::endl;





    // Crear un TPolyLine3D para cada trayectoria
    TPolyLine3D* track1 = new TPolyLine3D(2);
    track1->SetPoint(0, vertex.X(), vertex.Y(), vertex.Z());
    track1->SetPoint(1, silPoint0.X(), silPoint0.Y(), silPoint0.Z());
    track1->SetLineColor(kRed);
    track1->SetLineWidth(2);

    TPolyLine3D* track2 = new TPolyLine3D(2);
    track2->SetPoint(0, vertex.X(), vertex.Y(), vertex.Z());
    track2->SetPoint(1, silPoint1.X(), silPoint1.Y(), silPoint1.Z());
    track2->SetLineColor(kBlue);
    track2->SetLineWidth(2);

    // Dibujar la geometría primero
    sils->DrawGeo();  // Asegúrate de que esta función dibuja en el mismo TCanvas

    // Dibujar las trazas después
    track1->Draw("same");
    track2->Draw("same");
}
#endif