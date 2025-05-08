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
#include "TF1.h"
#include "TGraph.h"
#include "TSpline.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

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

void checkPIDtelescope()
{
    // Simulate 11Li transport in chamber
    ActRoot::TPCParameters tpc {"Actar"};
    // Silicons
    auto* sils {new ActPhysics::SilSpecs};
    sils->ReadFile("../configs/silicons_telescope.conf");
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
        if(name == "f2" || name == "f3")
            layer.MoveZTo(125, {0});
        if(name == "l0" || name == "r0")
            layer.MoveZTo(75, {3});
    }
    sils->DrawGeo();

    // SRIM
    std::string particle {"9Li"};
    std::string path{"../SRIM files/"};
    std::string gas{"900mb_CF4_90-10"};
    std::string silicon{"silicon"};
    auto* srim {new ActPhysics::SRIM};
    srim->ReadTable("ParticleGas", path + particle + "_" + gas + ".txt");
    srim->ReadTable("ParticleInSil", path + particle + "_" + silicon + ".txt");
    // LISE 
    std::string fileLISE {"../LISE files/11Li_silicon.txt"};
    // auto splineStragglingLISE {getSplineStragglingLISE(fileLISE)};
    srim->SetStragglingLISE("ParticleInSil", fileLISE);

    ROOT::Math::XYZPoint vertex {128, 128, 128}; // center of TPC

    // Histos
    auto hSPf0 {Histos::SP.GetHistogram()};
    hSPf0->SetTitle("SP for f0");
    auto hSPf2 {Histos::SP.GetHistogram()};
    hSPf2->SetTitle("SP for f2");
    auto hTheta11LiOut {Histos::ThetaLabHeavy.GetHistogram()};
    hTheta11LiOut->SetTitle("Theta 11Li out");
    auto hkinLi {Histos::KinHeavy.GetHistogram()};
    auto hEsilAftervsBefore {Histos::EsilAftervsBefore.GetHistogram()};
    hEsilAftervsBefore->SetTitle("EsilAfter vs Ebefore");
    // PIDs
    std::shared_ptr<TH2D> hPID = Histos::PIDLight.GetHistogram();
    std::shared_ptr<TH2D> hPIDLength = Histos::PIDLight.GetHistogram();
    // Straggling
    auto hStragglingBegining {Histos::Straggling.GetHistogram()};
    hStragglingBegining->SetTitle("Straggling begining");
    auto hStragglingEnd {Histos::Straggling.GetHistogram()};
    hStragglingEnd->SetTitle("Straggling end");
    auto hStragglingDistance {Histos::Straggling.GetHistogram()};
    hStragglingDistance->SetTitle("Straggling distance");
    auto hEnergyDiferenceWithoutStraggling {new TH1D("hEnergyDiferenceWithoutStraggling", "Energy difference without straggling", 100, -10, 10)};

    if (particle == "1H" || particle == "2H" || particle == "3H" || particle == "3He" || particle == "4He")
    {
        hPID = Histos::PIDLight.GetHistogram();
        hPID->SetTitle(("PID for " + particle).c_str());

        hPIDLength = Histos::PIDLightlength.GetHistogram();
        hPIDLength->SetTitle(("PID length for " + particle).c_str());
    }
    else
    {
        hPID = Histos::PIDHeavy.GetHistogram();
        hPID->SetTitle(("PID for " + particle).c_str());

        hPIDLength = Histos::PIDHeavylength.GetHistogram();
        hPIDLength->SetTitle(("PID length for " + particle).c_str());
    }

    int counter {0};
    bool stragglinSiliconsEnabled {true};


    // 11Li
    for(int i = 0; i < 100000; i++)
    {
        double udistRef {};
        double uRiniRef {}; 
        double uRafterRef {};

        double Tparticle {gRandom->Uniform(0,80)}; // MeV
        
        vertex.SetX(gRandom->Uniform(0,256));
        
        double phiParticle {gRandom->Uniform(0, 2 * TMath::Pi())};
        double thetaParticle {gRandom->Uniform(0, 8 * TMath::DegToRad())}; 
        // std::cout << "Theta: " << thetaLi11 * TMath::RadToDeg() << " Phi: " << phiLi11 * TMath::RadToDeg() << '\n';
        ROOT::Math::XYZVector directionParticle {TMath::Cos(thetaParticle), TMath::Sin(thetaParticle) * TMath::Sin(phiParticle),
                              TMath::Sin(thetaParticle) * TMath::Cos(phiParticle)};


        // Check hit for the 11Li 
        int silIndex = -1;
        ROOT::Math::XYZPoint silPoint;
        std::string layerHit;
        for(auto layer : silLayers)
        {
            std::tie(silIndex, silPoint) = sils->FindSPInLayer(layer, vertex, directionParticle);
            if(silIndex != -1)
            {
                layerHit = layer;
                break;
            }                
        }
        auto silNormal {sils->GetLayer(layerHit).GetNormal().Unit()};
        auto angleWithSil {TMath::ACos(directionParticle.Unit().Dot(silNormal))};
        if(silIndex == -1 || layerHit != "f2")
        {
            // hTheta11LiOut->Fill(thetaParticle * TMath::RadToDeg());
            // if(counter < 50)
            // {
            //     // Draw lines to check hits
            //     double x0 = vertex.X();
            //     double y0 = vertex.Y();
            //     double z0 = vertex.Z();
            //     // Definir la dirección escalada (para visualizar mejor)
            //     double scale = 200.0;  // Ajusta el tamaño de la línea
            //     double x1 = x0 + scale * directionParticle.X();
            //     double y1 = y0 + scale * directionParticle.Y();
            //     double z1 = z0 + scale * directionParticle.Z();
            //     // Crear la línea 3D
            //     double x[] = {x0, x1};
            //     double y[] = {y0, y1};
            //     double z[] = {z0, z1};
            //     TPolyLine3D* line = new TPolyLine3D(2, x, y, z);
            //     line->SetLineColor(kRed);
            //     line->SetLineWidth(2);
            //     line->DrawClone("same");
            // }
            continue; // If a silicon is not reached, don't continue with punchthough calculation
        }

        // Check if te particle hit also the second 0º silicon
        int silIndex1{};
        ROOT::Math::XYZPoint silPoint1{};
        double T3AfterSil1{-1};
        std::tie(silIndex1, silPoint1) = sils->FindSPInLayer("f3", vertex, directionParticle);
        if (silIndex1 == -1)
        {
            continue; // if don't hit second silicon continue, we would not have info for PID reconstruction
        }

        // Calculation of DeltaE-E
        // First, slow before silicon
        auto limitPointGas {ComputeLimitPoint(directionParticle, vertex)};
        double distanceInside {(vertex - limitPointGas).R()};
        double distanceGas {(vertex - silPoint).R()};
        auto energyBeforeSilicons {srim->SlowWithStraggling("ParticleGas", Tparticle, distanceGas)};
        // Second, slow in first silicon, DeltaE
        double energyAfterFirstSil {};
        double energyDiferenceWithoutStraggling {};
        if(stragglinSiliconsEnabled)
        {
            energyAfterFirstSil = srim->SlowWithStraggling("ParticleInSil", energyBeforeSilicons, sils->GetLayer(layerHit).GetUnit().GetThickness(), angleWithSil);
            if(energyAfterFirstSil == 0)
            {
                continue;
            }
            energyDiferenceWithoutStraggling = energyAfterFirstSil - srim->Slow("ParticleInSil", energyBeforeSilicons, sils->GetLayer(layerHit).GetUnit().GetThickness(), angleWithSil);
        }
        else
        {
            energyAfterFirstSil = srim->Slow("ParticleInSil", energyBeforeSilicons, sils->GetLayer(layerHit).GetUnit().GetThickness(), angleWithSil);
            if(energyAfterFirstSil == 0)
            {
                continue;
            }
        }
        auto DeltaE {energyBeforeSilicons - energyAfterFirstSil};
        if(!stragglinSiliconsEnabled)
        {
            DeltaE = gRandom->Gaus(DeltaE, silRes->Eval(DeltaE)); // after silicon resolution
        }
        // Third, slow in gas before second silicon
        auto distanceInterGas {(silPoint - silPoint1).R()};
        auto energyAfterInterGas {srim->SlowWithStraggling("ParticleGas", energyAfterFirstSil, distanceInterGas)};
        // Finally, slow in second silicon
        auto energyAfterSecondSil {srim->SlowWithStraggling("ParticleInSil", energyAfterInterGas, sils->GetLayer("f3").GetUnit().GetThickness(), angleWithSil)};

        double eLoss {energyAfterInterGas - energyAfterSecondSil};
        if(!stragglinSiliconsEnabled)
        {
            eLoss = gRandom->Gaus(eLoss, silRes->Eval(eLoss)); // after silicon resolution
        }
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

        hkinLi->Fill(thetaParticle * TMath::RadToDeg(), Tparticle);

        if(energyAfterSecondSil == 0 && eLoss >0)
        {
            hPID->Fill(eLoss, DeltaE);
            hPIDLength->Fill(eLoss, DeltaE / distanceInside);
            hEsilAftervsBefore->Fill(DeltaE, energyBeforeSilicons);
            // std::cout <<"Energy Before first silicon: " << energyBeforeSilicons << std::endl;
            // std::cout<< - energyBeforeSilicons + eLoss + DeltaE + (energyAfterFirstSil - energyAfterInterGas)  << std::endl;
        }
        hStragglingDistance->Fill(udistRef * 1000);
        hStragglingBegining->Fill(uRiniRef * 1000);
        hStragglingEnd->Fill(uRafterRef * 1000);
        hEnergyDiferenceWithoutStraggling->Fill(energyDiferenceWithoutStraggling);
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
    cPID->DivideSquare(3);
    cPID->cd(1);
    hPID->DrawClone("colz");
    cPID->cd(2);
    hPIDLength->DrawClone("colz");
    cPID->cd(3);
    hEsilAftervsBefore->DrawClone("colz");

    auto cStraggling {new TCanvas("cStraggling", "cStraggling", 800, 600)};
    cStraggling->DivideSquare(4);
    cStraggling->cd(1);
    hStragglingBegining->DrawClone("colz");
    cStraggling->cd(2);
    hStragglingEnd->DrawClone("colz");
    cStraggling->cd(3);
    hStragglingDistance->DrawClone("colz");
    cStraggling->cd(4);
    hEnergyDiferenceWithoutStraggling->DrawClone("colz");

    
    // Save histos
    // TFile* outFile = new TFile(("../DebugOutputs/checkPID_outputTelescope" + particle  + ".root").c_str(), "RECREATE");
    // hkinLi->Write("hkinLi");
    // hPID->Write("hPID");
    // hPIDLength->Write("hPIDLength");
    // outFile->Close();
    // delete outFile;
}
#endif