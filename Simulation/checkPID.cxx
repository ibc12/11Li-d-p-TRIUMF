#ifndef checkPID_cxx
#define checkPID_cxx

#include "ActSRIM.h"
#include "ActSilSpecs.h"
#include "ActTPCParameters.h"
#include "ActParticle.h"
#include "ActLine.h"

#include "Math/Point3D.h"
#include "TMath.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TPolyLine3D.h"
#include "TFile.h"

#include "../Histos.h"

ROOT::Math::XYZPoint ComputeLimitPoint(ROOT::Math::XYZVector directionHeavy, ROOT::Math::XYZPoint RP)
{
    // Compute the distance from RP to end of pad plane 0<X<256 0<Y<256 0<Z<256

    // Define the direction given by theta and phi
    ROOT::Math::XYZPointF firstPoint {RP};
    ROOT::Math::XYZPointF secondPoint {RP + 100 * directionHeavy};

    // Let's create a line and then compute the point where x=256 (low angle in heavy, so x will be always the limit)
    auto* line {new ActRoot::Line(firstPoint, secondPoint)};

    auto pointLimit {line->MoveToX(256)};
    ROOT::Math::XYZPoint pointLimitD {pointLimit};

    double distance {(RP - pointLimit).R()}; // distance in mm
    delete line;
    return pointLimitD;
}

double ComputeDistancef0toPoint(ROOT::Math::XYZVector directionHeavy, ROOT::Math::XYZPoint SP)
{
    // Define the direction given by theta and phi
    ROOT::Math::XYZPointF firstPoint {SP};
    ROOT::Math::XYZPointF secondPoint {SP + 100 * directionHeavy};

    // Let's create a line and then compute the point where x=256 (low angle in heavy, so x will be always the limit)
    auto* line {new ActRoot::Line(firstPoint, secondPoint)};

    auto pointf0 {line->MoveToX(362)}; // position of f0, 181 in pad units

    double distancef0f2 {(pointf0 - SP).R()}; // distance in mm
    delete line;
    return distancef0f2;
}

void checkPID()
{
    // Simulate 11Li transport in chamber
    ActRoot::TPCParameters tpc {"Actar"};
    // Silicons
    auto* sils {new ActPhysics::SilSpecs};
    sils->ReadFile("../configs/silicons_reverse.conf");
    sils->Print();
    const double sigmaSil {0.060 / 2.355}; // Si resolution
    auto silRes = std::make_unique<TF1>(
        "silRes", [=](double* x, double* p) { return sigmaSil * TMath::Sqrt(x[0] / 5.5); }, 0.0, 100.0, 1);
    std::vector<std::string> silLayers {"f0", "f2", "l0", "r0"};
    // We have to centre the silicons with the beam input
    // In real life beam window is not at Z / 2
    for(auto& [name, layer] : sils->GetLayers())
    {   
        if(name == "f0" || name == "f1")
            layer.MoveZTo(75, {3});
        if(name == "f2")
            layer.MoveZTo(125, {0});
        if(name == "l0" || name == "r0")
            layer.MoveZTo(75, {3});
    }
    sils->DrawGeo();

    // SRIM
    auto* srim {new ActPhysics::SRIM};
    srim->ReadTable("11LiGas", "../SRIM files/11Li_900mb_CF4_90-10.txt");
    srim->ReadTable("11LiInSil", "../SRIM files/11Li_silicon.txt");

    ROOT::Math::XYZPoint vertex {128, 128, 128}; // center of TPC

    // Histos
    auto hSPf0 {Histos::SP.GetHistogram()};
    hSPf0->SetTitle("SP for f0");
    auto hSPf2 {Histos::SP.GetHistogram()};
    hSPf2->SetTitle("SP for f2");
    auto hTheta11LiOut {Histos::ThetaLabHeavy.GetHistogram()};
    hTheta11LiOut->SetTitle("Theta 11Li out");
    auto hkinLi {Histos::KinHeavy.GetHistogram()};
    // PIDs
    auto hPID {Histos::PID.GetHistogram()};
    hPID->SetTitle("PID for 11Li");
    auto hPIDLength {Histos::PIDlength.GetHistogram()};
    hPIDLength->SetTitle("PID for 11Li length");

    int counter {0};

    // 11Li
    for(int i = 0; i < 100000; i++)
    {
        double TLi11 {gRandom->Uniform(50, 80)}; // MeV
        
        vertex.SetX(gRandom->Uniform(0,256));
        
        double phiLi11 {gRandom->Uniform(0, 2 * TMath::Pi())};
        double thetaLi11 {gRandom->Uniform(3 * TMath::DegToRad(), 7 * TMath::DegToRad())};
        std::cout << "Theta: " << thetaLi11 * TMath::RadToDeg() << " Phi: " << phiLi11 * TMath::RadToDeg() << '\n';
        ROOT::Math::XYZVector directionLi11 {TMath::Cos(thetaLi11), TMath::Sin(thetaLi11) * TMath::Sin(phiLi11),
                              TMath::Sin(thetaLi11) * TMath::Cos(phiLi11)};


        // Check hit for the 11Li 
        int silIndex = -1;
        ROOT::Math::XYZPoint silPoint;
        std::string layerHit;
        for(auto layer : silLayers)
        {
            std::tie(silIndex, silPoint) = sils->FindSPInLayer(layer, vertex, directionLi11);
            if(silIndex != -1)
            {
                layerHit = layer;
                break;
            }                
        }
        if(silIndex == -1)
        {
            hTheta11LiOut->Fill(thetaLi11 * TMath::RadToDeg());
            if(counter < 50)
            {
                // Draw lines to check hits
                double x0 = vertex.X();
                double y0 = vertex.Y();
                double z0 = vertex.Z();
                // Definir la dirección escalada (para visualizar mejor)
                double scale = 200.0;  // Ajusta el tamaño de la línea
                double x1 = x0 + scale * directionLi11.X();
                double y1 = y0 + scale * directionLi11.Y();
                double z1 = z0 + scale * directionLi11.Z();
                // Crear la línea 3D
                double x[] = {x0, x1};
                double y[] = {y0, y1};
                double z[] = {z0, z1};
                TPolyLine3D* line = new TPolyLine3D(2, x, y, z);
                line->SetLineColor(kRed);
                line->SetLineWidth(2);
                line->DrawClone("same");
            }
            
            continue; // If a silicon is not reached, don't continue with punchthough calculation
        }

        // Calculation of DeltaE-E
        // First, slow inside detector volume
        auto limitPointGas {ComputeLimitPoint(directionLi11, vertex)};
        double distanceInside {(vertex - limitPointGas).R()};
        auto energyAfterInside {srim->SlowWithStraggling("11LiGas", TLi11, distanceInside)};
        // auto DeltaE {TLi11 - energyAfterInside};
        // Second, slow in gas before silicon
        auto distanceInterGas {(limitPointGas - silPoint).R()};
        auto energyAfterInterGas {srim->SlowWithStraggling("11LiGas", energyAfterInside, distanceInterGas)};
        auto DeltaE {TLi11 - energyAfterInside};
        // Finally, slow in silicon
        auto energyAfterSil {srim->SlowWithStraggling("11LiInSil", energyAfterInterGas, sils->GetLayer(layerHit).GetUnit().GetThickness(), thetaLi11)};

        double eLoss {energyAfterInterGas - energyAfterSil};
        //if(layerHit == "f0")
        //{
        //    eLoss = energyAfterInterGas - energyAfterSil;
        //}
        //if(layerHit == "f2")
        //{
        //    double distanceBorderTof0 {ComputeDistancef0toPoint(directionLi11, limitPointGas)};
        //    double distancef0f2 {ComputeDistancef0toPoint(directionLi11, silPoint)};
        //    
        //    double energyJustBeforef0 {srim->SlowWithStraggling("11LiGas", energyAfterInside, distanceBorderTof0)};
        //    double energyJustBeforef2 {srim->SlowWithStraggling("11LiGas", energyJustBeforef0, distancef0f2)};
        //    
        //    double eLossBetweenf0Andf2 = energyJustBeforef0 - energyJustBeforef2;
        //
        //    eLoss = energyAfterInterGas - energyAfterSil;
        //}
        


        // Fill Histos
        if(layerHit == "f0")
        {
            hSPf0->Fill(silPoint.Y(), silPoint.Z());
        }
        if(layerHit == "f2")
        {
            hSPf2->Fill(silPoint.Y(), silPoint.Z());
        }

        hkinLi->Fill(thetaLi11 * TMath::RadToDeg(), TLi11);

        hPID->Fill(eLoss, DeltaE);
        hPIDLength->Fill(eLoss, DeltaE / distanceInside);
    }
    
    auto c {new TCanvas("c", "c", 800, 600)};
    c->DivideSquare(2);
    c->cd(1);
    hSPf0->DrawClone("colz");
    sils->GetLayer("f0").GetSilMatrix()->Draw();
    c->cd(2);   
    hSPf2->DrawClone("colz");
    sils->GetLayer("f2").GetSilMatrix()->Draw();

    auto c1 {new TCanvas("c1", "c1", 800, 600)};
    c1->DivideSquare(2);
    c1->cd(1);
    hTheta11LiOut->DrawClone();
    c1->cd(2);
    hkinLi->DrawClone("colz");

    auto cPID {new TCanvas("cPID", "cPID", 800, 600)};
    cPID->DivideSquare(2);
    cPID->cd(1);
    hPID->DrawClone("colz");
    cPID->cd(2);
    hPIDLength->DrawClone("colz");

    
    // Save histos
    TFile* outFile = new TFile("../DebugOutputs/checkPID_output11Li.root", "RECREATE");
    hkinLi->Write("hkinLi");
    hPID->Write("hPID");
    hPIDLength->Write("hPIDLength");
    outFile->Close();
    delete outFile;


}
#endif