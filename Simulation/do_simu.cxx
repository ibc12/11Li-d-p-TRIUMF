#ifndef triumf_cxx
#define triumf_cxx
#include "ActKinematics.h"
#include "ActSRIM.h"
#include "ActSilSpecs.h"
#include "ActTPCParameters.h"
#include "ActCrossSection.h"
#include "ActColors.h"
#include "ActParticle.h"
#include "ActDecayGenerator.h"

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

#include <cmath>
#include <iostream>
#include <string>
#include <unordered_map>

#include "../Histos.h"

using XYZPoint = ROOT::Math::XYZPoint;
using XYZVector = ROOT::Math::XYZVector;

XYZPoint SampleVertex(ActRoot::TPCParameters* tpc)
{
    // Define sigmas along Y and Z
    double sigmaY {4};
    double sigmaZ {4};
    auto y {gRandom->Gaus(tpc->Y() / 2, sigmaY)};
    auto z {gRandom->Gaus(tpc->Z() / 2, sigmaZ)};
    auto x {gRandom->Uniform() * tpc->X()};
    return {x, y, z};
}

std::pair<double, double> SampleCM()
{
    auto theta {TMath::ACos(gRandom->Uniform(-1, 1))};
    auto phi {gRandom->Uniform(0, TMath::TwoPi())};
    return {theta, phi};
}

void ApplyNaN(double& e, double t = 0, const std::string& comment = "stopped")
{
    if(e <= t)
        e = std::nan(comment.c_str());
}

void ApplyThetaRes(double& theta)
{
    double sigma {0.95 / 2.355}; // FWHM to sigma
    theta = gRandom->Gaus(theta, sigma * TMath::DegToRad());
}

double RandomizeBeamEnergy(double Tini, double sigma)
{
    return gRandom->Gaus(Tini, sigma);
}



void do_simu(const std::string& beam, const std::string& target, const std::string& light, const std::string& heavy, int neutronPS, int protonPS, double Tbeam, double Ex, bool inspect)
{
    // Ex = 0.435;
    // Set number of iterations
    auto niter {static_cast<int>(1e5)};
    gRandom->SetSeed(0);
    // Initialize detectors
    // TPC
    ActRoot::TPCParameters tpc {"Actar"};
    std::cout << "TPC: " << tpc.X() << " " << tpc.Y() << " " << tpc.Z() << '\n';
    // Silicons
    auto* sils {new ActPhysics::SilSpecs};
    sils->ReadFile("../configs/silicons_reverse.conf");
    sils->Print();
    const double sigmaSil {0.060 / 2.355}; // Si resolution
    auto silRes = std::make_unique<TF1>(
        "silRes", [=](double* x, double* p) { return sigmaSil * TMath::Sqrt(x[0] / 5.5); }, 0.0, 100.0, 1);
    std::vector<std::string> silLayers {"f0", "l0", "r0"};
    // We have to centre the silicons with the beam input
    // In real life beam window is not at Z / 2
    for(auto& [name, layer] : sils->GetLayers())
    {   
        if(name == "f0" || name == "f1")
            layer.MoveZTo(75, {3});
        if(name == "l0" || name == "r0")
            layer.MoveZTo(75, {3});
    }


        
    std::cout << "Sils Z centred at : " << tpc.Z() / 2 << " mm" << '\n';
    sils->DrawGeo();
    // This means: make the Z of silicons {5,6,...} be that zOfBeam.
    // shift the others accordingly
    // sils->DrawGeo();

    // Sigmas
    const double sigmaPercentBeam {0.001122};
    // Flags for resolution
    bool exResolution {true};

    // SRIM
    auto* srim {new ActPhysics::SRIM};
    // Transfer d,p
    srim->ReadTable("light", "../SRIM files/1H_900mb_CF4_90-10.txt");
    srim->ReadTable("beam", "../SRIM files/11Li_900mb_CF4_90-10.txt");
    srim->ReadTable("heavy", "../SRIM files/12Li_900mb_CF4_90-10.txt");
    srim->ReadTable("lightInSil", "../SRIM files/1H_silicon.txt");
    srim->ReadTable("heavyInSil", "../SRIM files/12Li_silicon.txt");
    // Transfer d,t
    // srim->ReadTable("light", "../SRIM files/3H_900mb_CF4_90-10.txt");
    // srim->ReadTable("beam", "../SRIM files/11Li_900mb_CF4_90-10.txt");
    // srim->ReadTable("heavy", "../SRIM files/10Li_900mb_CF4_90-10.txt");
    // srim->ReadTable("lightInSil", "../SRIM files/3H_silicon.txt");
    // srim->ReadTable("heavyInSil", "../SRIM files/10Li_silicon.txt");

    // Kinematics
    auto* kin {new ActPhysics::Kinematics {beam, target, light, heavy, Tbeam, Ex}};

    // cross section sampler
    auto* xs {new ActSim::CrossSection()};
    if(Ex == 0.)
    {
        TString data_to_read {TString::Format("../Inputs/TheoXS/%.1fMeV/dp/angs12nospin.dat", Tbeam / 11)};
        xs->ReadFile(data_to_read.Data());
        std::cout<<xs->GetTotalXSmbarn()<<std::endl;
    }
    else if(Ex == 0.130)
    {
        TString data_to_read {TString::Format("../Inputs/TheoXS/%.1fMeV/dp/angp12nospin.dat", Tbeam / 11)};
        xs->ReadFile(data_to_read.Data());
        std::cout<<xs->GetTotalXSmbarn()<<std::endl;
    }
    else if(Ex == 0.435)
    {
        TString data_to_read {TString::Format("../Inputs/TheoXS/%.1fMeV/dp/angp32nospin.dat", Tbeam / 11)};
        xs->ReadFile(data_to_read.Data());
        std::cout<<xs->GetTotalXSmbarn()<<std::endl;
    }
    xs->DrawTheo();
    double x {xs->GetTotalXSmbarn()};
    std::cout << "Total XS = " << x << " mb" << '\n';

    // Declare histograms
    auto hKin {Histos::Kin.GetHistogram()};
    hKin->SetTitle("Sampled kinematics");
    auto hKinRec {Histos::Kin.GetHistogram()};
    hKinRec->SetTitle("Reconstructed kinematics");
    // Silicon hits
    auto hSPf0 {Histos::SP.GetHistogram()};
    hSPf0->SetTitle("SP for f0");
    auto hSPl0 {Histos::SP.GetHistogram()};
    hSPl0->SetTitle("SP for l0");
    auto hSPr0 {Histos::SP.GetHistogram()};
    hSPr0->SetTitle("SP for r0");
    // RP and angles
    auto hRP {Histos::RP.GetHistogram()};
    auto hRP_ZY {Histos::RP_ZY.GetHistogram()};
    auto hThetaCMAll {Histos::ThetaCM.GetHistogram()};
    auto hThetaCM {Histos::ThetaCM.GetHistogram()};
    auto hThetaLabAll {Histos::ThetaLab.GetHistogram()};
    auto hThetaLab {Histos::ThetaLab.GetHistogram()};
    auto hEx {Histos::Ex.GetHistogram()};
    auto hThetaCMThetaLab {Histos::ThetaCMThetaLab.GetHistogram()};
    auto hRPxEfirstSil {Histos::RP_E.GetHistogram()};
    hRPxEfirstSil->SetTitle("RPvsE if just 1 Sil");
    auto hRPxEbothSil {Histos::RP_E.GetHistogram()};
    hRPxEbothSil->SetTitle("RPvsE if 2 Sil");
    auto hRPxEStoppedGas {Histos::RP_E.GetHistogram()};
    auto hRPeffAll {Histos::RP_eff.GetHistogram()};
    auto hRPeffIn {Histos::RP_eff.GetHistogram()};
    auto hThetaVertexInGas {Histos::ThetaLabVertex.GetHistogram()};
    auto hRangeGas {Histos::RangeInGas.GetHistogram()};
    auto hPhiLab {Histos::PhiLab.GetHistogram()};
    auto hPhiCM {Histos::PhiLab.GetHistogram()};
    hPhiCM->SetTitle("PhiCM");
    // Partial Efficeincies
    auto hThetaLabNoReachSil {Histos::ThetaLab.GetHistogram()};
    auto hThetaCMNoReachSil {Histos::ThetaCM.GetHistogram()};
    auto hThetaLabDirectionSil {Histos::ThetaLab.GetHistogram()};
    auto hThetaCMDirectionSil {Histos::ThetaCM.GetHistogram()};
    auto hThetaLabMeassureSil {Histos::ThetaLab.GetHistogram()};
    auto hThetaCMMeassureSil {Histos::ThetaCM.GetHistogram()};
    // Heavy Particle
    auto hRangeHeavyInGas {Histos::RangeInGas.GetHistogram()};
    hRangeHeavyInGas->SetTitle("Range of heavy in gas");
    auto hRangeHeavyInSilAfterGas {Histos::RangeInSil.GetHistogram()};
    hRangeHeavyInSilAfterGas->SetAxisRange(0, 1);
    hRangeHeavyInSilAfterGas->SetTitle("Range of heavy in silicon after gas");
    auto hT4Lab {Histos::T4Lab.GetHistogram()};
    auto hT4LabAfterGas {Histos::T4Lab.GetHistogram()};
    hT4LabAfterGas->SetTitle("Kinetic Energy of heavy after gas");
    auto hTbeamAfterGas {Histos::T4Lab.GetHistogram()};
    hTbeamAfterGas->SetTitle("Kinetic Energy of beam after gas");
    auto hTheta4Lab {Histos::ThetaLabHeavy.GetHistogram()};
    // Lateral silicon plots, why less counting
    auto hSPl0noCut {Histos::SP.GetHistogram()};
    hSPl0noCut->SetTitle("SP for l0 no Cuts");
    auto hSPr0noCut {Histos::SP.GetHistogram()};
    hSPr0noCut->SetTitle("SP for r0 no Cuts");
    auto hTheta3LabReachLatSil {Histos::ThetaLab.GetHistogram()};
    hTheta3LabReachLatSil->SetTitle("Theta3Lab for particles that reach lateral silicons");
    auto hTheta3LabNoReachLatSil {Histos::ThetaLab.GetHistogram()};
    hTheta3LabNoReachLatSil->SetTitle("Theta3Lab for particles that don't reach lateral silicons");
    auto hT3LabNoReachLatSil {Histos::T3Lab.GetHistogram()};
    hT3LabNoReachLatSil->SetTitle("T3Lab for particles that reach lateral silicons");
    auto hT3LabReachLatSil {Histos::T3Lab.GetHistogram()};
    hT3LabReachLatSil->SetTitle("T3Lab for particles that reach lateral silicons");
    auto hTheta3LabLatSilPunch {Histos::ThetaLab.GetHistogram()};
    hTheta3LabLatSilPunch->SetTitle("Theta3Lab for particles that punchthrough lateral silicons");
    auto hRangeInGasSin {Histos::RangeInGasSin.GetHistogram()};
    auto hdistanceSPtol0 {Histos::DistanceSP.GetHistogram()};
    hdistanceSPtol0->SetTitle("Distance to SP for l0");
    auto hdistanceSPtor0 {Histos::DistanceSP.GetHistogram()};
    hdistanceSPtor0->SetTitle("Distance to SP for r0");
    auto hRange1 {Histos::RangeInGasSin.GetHistogram()};
    hRange1->SetTitle("Range");
    auto hT3test {Histos::T3Lab.GetHistogram()};
    hT3test->SetTitle("T3Lab");
    // Kin and eff for lat silicons
    auto hKinLatSil {Histos::Kin.GetHistogram()};
    hKinLatSil->SetTitle("Kinematics for particles that go to lateral silicons (impact + stop)");
    auto hKinLatSilReach {Histos::Kin.GetHistogram()};
    hKinLatSilReach->SetTitle("Kinematics for particles that reach lateral silicons");
    auto hKinLatStopInside {Histos::Kin.GetHistogram()};
    hKinLatStopInside->SetTitle("Kinematics for particles that stop inside chamber");
    auto hKinStopInside {Histos::Kin.GetHistogram()};
    hKinStopInside->SetTitle("Kinematics for particles that stop inside chamber");
    auto hKinNoReachSil {Histos::Kin.GetHistogram()};
    hKinNoReachSil->SetTitle("Kinematics for particles that don't have a SP");
    auto hKinStopBetweenSilAndGas {Histos::Kin.GetHistogram()};
    hKinStopBetweenSilAndGas->SetTitle("Kinematics for particles that stop between lateral silicons and gas");
    auto hThetaLabStopInside {Histos::ThetaLab.GetHistogram()};
    hThetaLabStopInside->SetTitle("ThetaLab for particles that stop inside chanber");
    auto hThetaCMStopInside {Histos::ThetaCM.GetHistogram()};
    hThetaCMStopInside->SetTitle("ThetaCM for particles that stop inside chanber");
    
    // End histograms declaration

    // Gas information and normalization factor
    double beamParticles {(2e3) * 6 * 24 * 3600}; // 3e5 pps 6 days
    double gasDensity {4.4807e-4}; // g/cm3
    double gasMolarDensity {0.9 * 4.0282 + 0.1 * 87.993612}; // g/mol 
    double gasParticles {(gasDensity/gasMolarDensity) * 6.022e23 * 25.6}; // particles/cm3 * ACTAR length
    double alpha {beamParticles * gasParticles * xs->GetTotalXScm2() / niter};
    std::cout<<"Norm value: "<<alpha<<std::endl;

    // File to save data
    TString fileName {
        TString::Format("../Outputs/%.1fMeV/transfer_TRIUMF_Eex_%.3f_nPS_%d_pPS_%d.root", Tbeam / 11, Ex, neutronPS, protonPS)};
    auto* outFile {new TFile(fileName, "recreate")};
    auto* outTree {new TTree("SimulationTTree", "A TTree containing only our Eex obtained by simulation")};
    double theta3CM_tree {};
    outTree->Branch("theta3CM", &theta3CM_tree);
    double Eex_tree {};
    outTree->Branch("Eex", &Eex_tree);
    double EVertex_tree {};
    outTree->Branch("EVertex", &EVertex_tree);
    double theta3Lab_tree {};
    outTree->Branch("theta3Lab", &theta3Lab_tree);
    double phi3CM_tree {};
    outTree->Branch("phi3CM", &phi3CM_tree);
    auto* outTreeHeavy {new TTree("SimulationTTreeHeavy", "A TTree containing heavy information")};
    double theta4Lab_tree {};
    outTreeHeavy->Branch("theta4Lab", &theta4Lab_tree);
    double phi4Lab_tree {};
    outTreeHeavy->Branch("phi4Lab", &phi4Lab_tree);
    double T4Lab_tree {};
    outTreeHeavy->Branch("T4Lab", &T4Lab_tree);
    ROOT::Math::XYZPoint RP_tree {};
    outTreeHeavy->Branch("RP", &RP_tree);

    // Set Random Ex if needed (no xs available, so will be uniform distributed)
    bool setRandomEx {false};
    if(setRandomEx)
    {
        Ex = 0.;
    }

    // RUN!
    // print fancy info (dont pay attention to this)
    std::cout << BOLDMAGENTA << "Running for Ex = " << Ex << " MeV" << RESET << '\n';
    std::cout << BOLDGREEN;
    const int percentPrint {5};
    int step {niter / (100 / percentPrint)};
    int nextPrint {step};
    int percent {};
    for(int it = 0; it < niter; it++)
    {
        // Print progress
        if(it >= nextPrint)
        {
            percent = 100 * it / niter;
            std::cout << "\r" << std::string(percent / percentPrint, '|') << percent << "%";
            std::cout.flush();
            nextPrint += step;
        }
        // Sample vertex position
        auto vertex {SampleVertex(&tpc)};
        hRPeffAll->Fill(vertex.X());
        // Randomize Ex in a BW distribution
        double randEx = Ex;
        if(exResolution)
        {
            if(Ex == 0)
            {
                randEx = gRandom->BreitWigner(Ex, 0.1);
            }
            else if(Ex == 0.130)
            {
                randEx = gRandom->BreitWigner(Ex, 0.015);
            }
            else if(Ex == 0.435)
            {
                randEx = gRandom->BreitWigner(Ex, 0.08);
            }
        }

        // Randomize beam energy, slow beam with straggling and check if reaction can happen
        auto TbeamRand = RandomizeBeamEnergy(Tbeam, sigmaPercentBeam * Tbeam);
        auto TbeamCorr {srim->SlowWithStraggling("beam", TbeamRand, vertex.X())};
        auto beamThreshold {ActPhysics::Kinematics(beam, target, light, heavy, -1, randEx).GetT1Thresh()};
        if(std::isnan(TbeamCorr) || TbeamCorr < beamThreshold){
            continue;
        }
        // Obtein angles
        double theta3CMBefore {-1};
        if(setRandomEx)
        {
            auto [sampledTheta3CM, phi3CM] {SampleCM()};
            theta3CMBefore = sampledTheta3CM * TMath::RadToDeg();
        }
        else
        {
            while(theta3CMBefore < 0){theta3CMBefore = xs->Sample(gRandom->Uniform());} // sample in deg
        }
        double phi3CM {gRandom -> Uniform(0, 2 * TMath::Pi())};
        kin->SetBeamEnergyAndEx(TbeamCorr, randEx);
        // Sample Ex/Ecm/thetaCM...
        // Generate lab kinematics for protons
        kin->ComputeRecoilKinematics(theta3CMBefore * TMath::DegToRad(), phi3CM);
        // Fill thetaCMall
        hThetaCMAll->Fill(theta3CMBefore);
        // Extract direction
        auto T3Lab {kin->GetT3Lab()};
        auto theta3Lab {kin->GetTheta3Lab()};
        // Save without resolution
        auto theta3LabSampled {theta3Lab};
        // Apply angle resolution
        ApplyThetaRes(theta3Lab);
        hThetaLabAll->Fill(theta3Lab * TMath::RadToDeg());
        double theta3CM {kin->ReconstructTheta3CMFromLab(Tbeam, theta3Lab)};
        auto phi3Lab {kin->GetPhi3Lab()};
        hPhiLab->Fill(phi3Lab * TMath::RadToDeg());
        hPhiCM->Fill(phi3CM * TMath::RadToDeg());
        XYZVector direction {TMath::Cos(theta3Lab), TMath::Sin(theta3Lab) * TMath::Sin(phi3Lab),
                             TMath::Sin(theta3Lab) * TMath::Cos(phi3Lab)};
        // Analysis of heavy particle range, 12Li decays to 11Li
        double theta4Lab {kin->GetTheta4Lab()};
        hTheta4Lab->Fill(theta4Lab * TMath::RadToDeg());
        auto phi4Lab {kin->GetPhi4Lab()};
        XYZVector directionHeavy {TMath::Cos(theta4Lab), TMath::Sin(theta4Lab) * TMath::Sin(phi4Lab),
                             TMath::Sin(theta4Lab) * TMath::Cos(phi4Lab)};
        double T4Lab {kin->GetT4Lab()};
        hT4Lab->Fill(T4Lab);
        double T4LabAfterGas {srim->Slow("heavy", T4Lab, (256 - vertex.X()) + 100)}; // Very low angle, so not took into acount
        hT4LabAfterGas->Fill(T4LabAfterGas);
        double rangeHeavyInGas {srim->EvalRange("heavy", T4Lab)};
        hRangeHeavyInGas->Fill(rangeHeavyInGas);
        double rangeHeavyInSilAfterGas {srim->EvalRange("heavyInSil", T4LabAfterGas)}; 
        hRangeHeavyInSilAfterGas->Fill(rangeHeavyInSilAfterGas);
        // Analysis of unreacted beam particle energy in silicons
        double TbeamAfterGas {srim->Slow("beam", TbeamCorr, 256 + 100)}; // ACTAR lenght + silicon separation
        hTbeamAfterGas->Fill(TbeamAfterGas);
        // Threshold L1, particles that stop in actar. Check before doing the continues
        double rangeInGas {srim->EvalRange("light", T3Lab)};
        ROOT::Math::XYZPoint finalPointGas {vertex + rangeInGas * direction.Unit()};
        if(0 <= finalPointGas.X() && finalPointGas.X() <= 256 && 0 <= finalPointGas.Y() && finalPointGas.Y() <= 256 && 0 <= finalPointGas.Z() && finalPointGas.Z() <= 256)
        {
            double distanceXY {TMath::Sqrt(pow(vertex.X() - finalPointGas.X(),2) + pow(vertex.Y() - finalPointGas.Y(),2))};
            if(distanceXY >= 0)
            {
                hThetaVertexInGas->Fill(theta3Lab * TMath::RadToDeg(), vertex.X());
                hRangeGas->Fill(rangeInGas);
                hKinStopInside->Fill(theta3Lab * TMath::RadToDeg(), T3Lab);

                hRangeInGasSin->Fill(theta3Lab * TMath::RadToDeg(), rangeInGas * TMath::Sin(theta3Lab));
                hRange1->Fill(theta3Lab * TMath::RadToDeg(), rangeInGas);
                hT3test->Fill(T3Lab);

                // Fill also eff hists
                hThetaCM->Fill(theta3CMBefore);
                hThetaLab->Fill(theta3Lab * TMath::RadToDeg());
                hThetaLabStopInside->Fill(theta3Lab * TMath::RadToDeg());
                hThetaCMStopInside->Fill(theta3CMBefore);
            }
        }
        // How to check whether tracks would read the silicons with new class:
        int silIndex0 = -1;
        ROOT::Math::XYZPoint silPoint0;
        std::string layer0;
        for(auto layer : silLayers)
        {
            std::tie(silIndex0, silPoint0) = sils->FindSPInLayer(layer, vertex, direction);
            if(silIndex0 != -1)
            {
                layer0 = layer;
                break;
            }                
        }            
        // "f0": key name of layer to check for SP
        // silIndex == -1 if NO SP
        // else, returns the silicon index
        // silPoint0: SP in millimetres. in standard ACTAR frame (same as analysis)
        if(silIndex0 == -1)
        {
            hKinNoReachSil->Fill(theta3Lab * TMath::RadToDeg(), T3Lab);
            hThetaLabNoReachSil->Fill(theta3Lab * TMath::RadToDeg());
            hThetaCMNoReachSil->Fill(theta3CMBefore);
            continue;
        }
        hThetaLabDirectionSil->Fill(theta3Lab * TMath::RadToDeg());
        hThetaCMDirectionSil->Fill(theta3CMBefore);
            
        // Fill hist to get counts in laterals with no cuts
        if(layer0 == "l0")
        {
            hSPl0noCut->Fill(silPoint0.X(), silPoint0.Z());
            hKinLatSil->Fill(theta3Lab * TMath::RadToDeg(), T3Lab);
            hdistanceSPtol0->Fill((silPoint0 - vertex).R());
        }
        if(layer0 == "r0")
        {
            hSPr0noCut->Fill(silPoint0.X(), silPoint0.Z());
            hKinLatSil->Fill(theta3Lab * TMath::RadToDeg(), T3Lab);
            hdistanceSPtor0->Fill((silPoint0 - vertex).R());
        }
        
        // Slow down light in gas
        auto T3AtSil {srim->SlowWithStraggling("light", T3Lab, (silPoint0 - vertex).R())};
        // Check if stopped
        ApplyNaN(T3AtSil);
        if(std::isnan(T3AtSil))
        {
            hThetaLabNoReachSil->Fill(theta3Lab * TMath::RadToDeg());
            hThetaCMNoReachSil->Fill(theta3CMBefore);
            hRPxEStoppedGas->Fill(vertex.X(), T3Lab);
            if(layer0 == "l0" || layer0 == "r0")
            {
                hTheta3LabNoReachLatSil->Fill(theta3Lab * TMath::RadToDeg());
                hT3LabNoReachLatSil->Fill(T3Lab);
                if(0 <= finalPointGas.X() && finalPointGas.X() <= 256 && 0 <= finalPointGas.Y() && finalPointGas.Y() <= 256 && 0 <= finalPointGas.Z() && finalPointGas.Z() <= 256)
                {
                    hKinLatStopInside->Fill(theta3Lab * TMath::RadToDeg(), T3Lab);
                }
                else
                {
                    hKinStopBetweenSilAndGas->Fill(theta3Lab * TMath::RadToDeg(), T3Lab);
                }
            }   
            continue;
        }
        
        // Slow down in silicon
        auto normal {sils->GetLayer(layer0).GetNormal()};
        auto angleWithNormal {TMath::ACos(direction.Unit().Dot(normal.Unit()))};
        auto T3AfterSil0 {srim->SlowWithStraggling("lightInSil", T3AtSil, sils->GetLayer(layer0).GetUnit().GetThickness(),
                                                   angleWithNormal)};
        auto eLoss0preSilRes {T3AtSil - T3AfterSil0};
        auto eLoss0 {gRandom->Gaus(eLoss0preSilRes, silRes->Eval(eLoss0preSilRes))}; // after silicon resolution
        ApplyNaN(eLoss0, sils->GetLayer(layer0).GetThresholds().at(silIndex0));
        int count = 0;
        if(std::isnan(eLoss0))
        {
            hThetaLabNoReachSil->Fill(theta3Lab * TMath::RadToDeg());
            hThetaCMNoReachSil->Fill(theta3CMBefore);

            hTheta3LabNoReachLatSil->Fill(theta3Lab * TMath::RadToDeg());
            hT3LabNoReachLatSil->Fill(T3Lab);
            
            if(layer0 == "l0" || layer0 == "r0")
                hKinStopBetweenSilAndGas->Fill(theta3Lab * TMath::RadToDeg(), T3Lab); // No having more energy than threshold = stop between gas and silicon
            continue;
        }
            

        // Analysis of punchthrough in lateral silicons
        if(T3AfterSil0 > 0 && (layer0 == "l0" || layer0 == "r0"))
        {
            hTheta3LabLatSilPunch->Fill(theta3Lab * TMath::RadToDeg());
        }


        // Apply 2nd layer of silicons
        double T3AfterInterGas {};
        int silIndex1 {};
        ROOT::Math::XYZPoint silPoint1 {};
        double eLoss1 {};
        double T3AfterSil1 {-1};
        if(T3AfterSil0 > 0. && layer0 == "f0") 
        {
            std::tie(silIndex1, silPoint1) = sils->FindSPInLayer("f1", vertex, direction);
            if(silIndex1 == -1)
            {} // If a silicon is not reached, don't continue with punchthough calculation
            else
            {
                T3AfterInterGas = {srim->SlowWithStraggling("light", T3AfterSil0, (silPoint0 - silPoint1).R())};
                if(T3AfterInterGas == 0)
                {} // If slow in gas don't continue with calculation
                else
                {
                    T3AfterSil1 =srim->SlowWithStraggling("lightInSil", T3AfterInterGas, sils->GetLayer("f1").GetUnit().GetThickness(),
                                                   angleWithNormal);
                    auto eLoss1preSilRes {T3AfterInterGas - T3AfterSil1};
                    eLoss1 = gRandom->Gaus(eLoss1preSilRes, silRes->Eval(eLoss1preSilRes)); // after silicon resolution
                    ApplyNaN(eLoss1, sils->GetLayer("f1").GetThresholds().at(silIndex1));
                    if(std::isnan(eLoss1))
                        eLoss1 = 0;
                }
            }
        }

        // Reconstruct!
        bool isOk {T3AfterSil0 == 0 || T3AfterSil1 == 0}; // no punchthrouhg
        if(isOk)
        {
            // Assuming no punchthrough!
            double T3Rec {};
            if(eLoss1 == 0)
            {
                T3Rec = srim->EvalInitialEnergy("light", eLoss0, (silPoint0 - vertex).R());
                hRPxEfirstSil->Fill(vertex.X(), T3Rec);
            }
            else
            {
                auto T3Rec0 {srim->EvalInitialEnergy("light", eLoss0, (silPoint0 - vertex).R())};
                hRPxEfirstSil->Fill(vertex.X(), T3Rec0); // This is to show what would be obtained if there was just one sil

                // Reconstruction ob T3 with 2 silicon layers
                auto T3Rec1 {srim->EvalInitialEnergy("light", eLoss1, (silPoint1 - silPoint0).R())};
                T3Rec = srim->EvalInitialEnergy("light", eLoss0+T3Rec1, (silPoint0 - vertex).R());
            }

            auto ExRec {kin->ReconstructExcitationEnergy(T3Rec, theta3Lab)};

            // Fill
            hKin->Fill(theta3LabSampled * TMath::RadToDeg(), T3Lab); // before uncertainties implemented
            hKinRec->Fill(theta3Lab * TMath::RadToDeg(), T3Rec);     // after reconstruction
            hEx->Fill(ExRec);
            hRP->Fill(vertex.X(), vertex.Y());
            hRP_ZY->Fill(vertex.Y(), vertex.Z());
            if(layer0 == "f0")
            {
                hSPf0->Fill(silPoint0.Y(), silPoint0.Z());
            }
            if(layer0 == "l0")
            {
                hSPl0->Fill(silPoint0.X(), silPoint0.Z());
                hTheta3LabReachLatSil->Fill(theta3Lab * TMath::RadToDeg());
                hKinLatSilReach->Fill(theta3Lab * TMath::RadToDeg(), T3Lab);
            }
            if(layer0 == "r0")
            {
                hSPr0->Fill(silPoint0.X(), silPoint0.Z());
                hTheta3LabReachLatSil->Fill(theta3Lab * TMath::RadToDeg());
                hKinLatSilReach->Fill(theta3Lab * TMath::RadToDeg(), T3Lab);
            }
            hThetaCM->Fill(theta3CMBefore); // only thetaCm that enter our cuts
            hThetaLab->Fill(theta3Lab * TMath::RadToDeg());
            hThetaCMThetaLab->Fill(theta3CMBefore, theta3LabSampled * TMath::RadToDeg());
            hRPxEbothSil->Fill(vertex.X(), T3Rec);
            hRPeffIn->Fill(vertex.X());
            
            hThetaLabMeassureSil->Fill(theta3Lab * TMath::RadToDeg());
            hThetaCMMeassureSil->Fill(theta3CMBefore);

            // write to TTree
            Eex_tree = ExRec;
            theta3CM_tree = theta3CM * TMath::RadToDeg();
            EVertex_tree = T3Rec;
            theta3Lab_tree = theta3Lab * TMath::RadToDeg();
            phi3CM_tree = phi3CM;
            outTree->Fill();
            theta4Lab_tree = theta4Lab;
            phi4Lab_tree = phi4Lab;
            T4Lab_tree = T4Lab;
            RP_tree = vertex;
            outTreeHeavy->Fill();
        }
    }

    // outFile->Write();
    outFile->Close();

    // Compute efficiency
    auto* effCM {new TEfficiency {*hThetaCM, *hThetaCMAll}};
    effCM->SetNameTitle("effCM", " #epsilon_{TOT} (#theta_{CM});#epsilon;#theta_{CM} [#circ]");
    auto* effLab {new TEfficiency {*hThetaLab, *hThetaLabAll}};
    effLab->SetNameTitle("effLab", "#epsilon_{TOT} (#theta_{Lab});#epsilon;#theta_{Lab} [#circ]");
    auto* effLatSil {new TEfficiency {*hTheta3LabReachLatSil, *hThetaLabAll}};
    effLatSil->SetNameTitle("effLatSil", "#epsilon_{TOT} (#theta_{Lab}) Lateral Silicon;#epsilon;#theta_{Lab} [#circ]");
    auto* effStopInside {new TEfficiency {*hThetaLabStopInside, *hThetaLabAll}};
    effStopInside->SetNameTitle("effStopInside", "#epsilon_{TOT} (#theta_{Lab}) events stopped inside;#epsilon;#theta_{Lab} [#circ]");
    auto effDetectionSilLab {new TEfficiency {*hThetaLabMeassureSil, *hThetaLabDirectionSil}};
    effDetectionSilLab->SetNameTitle("effDetectionSil", "#epsilon_{Sil} (#theta_{Lab}) events detected in silicons that go to them;#epsilon;#theta_{Lab} [#circ]");
    auto effDetectionSilCM {new TEfficiency {*hThetaCMMeassureSil, *hThetaCMDirectionSil}};
    effDetectionSilCM->SetNameTitle("effDetectionSilCM", "#epsilon_{Sil} (#theta_{CM}) events detected in silicons that go to them;#epsilon;#theta_{CM} [#circ]");
    auto effDetectionL1Lab {new TEfficiency {*hThetaLabStopInside, *hThetaLabNoReachSil}};
    effDetectionL1Lab->SetNameTitle("effDetectionL1Lab", "#epsilon_{L1} (#theta_{Lab}) events detected L1 (no Sil signal);#epsilon;#theta_{Lab} [#circ]");
    auto effDetectionL1CM {new TEfficiency {*hThetaCMStopInside, *hThetaCMNoReachSil}};
    effDetectionL1CM->SetNameTitle("effDetectionL1CM", "#epsilon_{L1} (#theta_{CM}) events detected L1 (no Sil signal);#epsilon;#theta_{CM} [#circ]");
    // Now efficiency in intervals
    auto* effIntervals {new TEfficiency {*hRPeffIn, *hRPeffAll}};
    effIntervals->SetNameTitle("effIntervals", "#RPx efficiency;#epsilon;#RPx [#mm]");


    // Draw if not running for multiple Exs
    if(inspect)
    {
        auto* cSP {new TCanvas {"cSP", "Sil Points"}};
        cSP->DivideSquare(3);
        cSP->cd(1);
        hSPf0->DrawClone("colz");
        sils->GetLayer("f0").GetSilMatrix()->Draw();
        cSP->cd(2);
        hSPl0->DrawClone("colz");
        sils->GetLayer("l0").GetSilMatrix()->Draw();
        cSP->cd(3);
        hSPr0->DrawClone("colz");
        sils->GetLayer("r0").GetSilMatrix()->Draw();

        auto* c0 {new TCanvas {"c0", "Sim inspect 0"}};
        c0->DivideSquare(6);
        c0->cd(1);
        hKin->DrawClone("colz");
        // Draw theo kin
        kin->SetBeamEnergyAndEx(Tbeam, Ex);
        auto* gtheo {kin->GetKinematicLine3()};
        gtheo->Draw("l");
        c0->cd(2);
        hKinRec->DrawClone("colz");
        gtheo->Draw("l");
        c0->cd(3);
        c0->cd(4);
        hRP->DrawClone("colz");
        c0->cd(5);
        hThetaCMAll->SetTitle("All #theta_{CM}");
        hThetaCMAll->DrawClone();
        c0->cd(6);
        hThetaCM->SetTitle("#theta_{CM} in cuts");
        hThetaCM->DrawClone();

        auto* c1 {new TCanvas {"c1", "Sim inspect 1"}};
        c1->DivideSquare(6);
        c1->cd(1);
        hPhiLab->DrawClone();
        c1->cd(2);
        hEx->DrawClone();
        c1->cd(3);
        hThetaCMThetaLab->DrawClone("colz");
        auto* gCMLab {kin->GetThetaLabvsThetaCMLine()};
        gCMLab->Draw("l");
        c1->cd(4);
        hRPxEfirstSil->DrawClone();
        c1->cd(5);
        hRPxEbothSil->DrawClone();
        c1->cd(6);
        hRPxEStoppedGas->DrawClone();

        auto* cEff {new TCanvas {"cEff", "Eff plots"}};
        cEff->DivideSquare(7);
        cEff->cd(1);
        effIntervals->Draw("apl");
        cEff->cd(2);
        effCM->Draw("apl");
        cEff->cd(3);
        effLab->Draw("apl");
        cEff->cd(4);
        hThetaCMAll->DrawClone();
        cEff->cd(5);
        hThetaCM->DrawClone();
        cEff->cd(6);
        hThetaLabAll->DrawClone();
        cEff->cd(7);
        hThetaLab->DrawClone();

        auto* cL1 {new TCanvas {"cL1", "Plots for L1 events"}};
        cL1->DivideSquare(6);
        cL1->cd(1);
        hThetaLabStopInside->DrawClone();
        cL1->cd(2);
        hThetaVertexInGas->DrawClone("colz");
        cL1->cd(3);
        hRangeGas->DrawClone();

        auto* cHeavy {new TCanvas {"cHeavy", "Heavy Particle Study"}};
        cHeavy->DivideSquare(6);
        cHeavy->cd(1);
        hT4Lab->DrawClone();
        cHeavy->cd(2);
        hRangeHeavyInGas->DrawClone();
        cHeavy->cd(3);
        hRangeHeavyInSilAfterGas->DrawClone();
        cHeavy->cd(4);
        hTheta4Lab->DrawClone();
        cHeavy->cd(5);
        hT4LabAfterGas->DrawClone();
        cHeavy->cd(6);
        hTbeamAfterGas->DrawClone();

        auto* cLat {new TCanvas {"cLatnoCut", "lateral Silicon Study"}};
        cLat->DivideSquare(7);
        cLat->cd(1);
        hSPl0noCut->DrawClone("colz");
        sils->GetLayer("l0").GetSilMatrix()->DrawClone();
        cLat->cd(2);
        hSPr0noCut->DrawClone("colz");
        sils->GetLayer("r0").GetSilMatrix()->DrawClone();
        cLat->cd(3);
        hTheta3LabReachLatSil->DrawClone();
        cLat->cd(4);
        hTheta3LabNoReachLatSil->DrawClone();
        cLat->cd(5);
        hT3LabNoReachLatSil->DrawClone();
        cLat->cd(6);
        hTheta3LabLatSilPunch->DrawClone();
        cLat->cd(7);

        auto cLatKin {new TCanvas {"cLatKin", "Kinematics for lateral silicons"}};
        cLatKin->DivideSquare(10);
        cLatKin->cd(1);
        hKinLatSil->DrawClone();
        cLatKin->cd(2);
        hKinLatSilReach->DrawClone();
        cLatKin->cd(3);
        hKinStopBetweenSilAndGas->DrawClone();
        cLatKin->cd(4);
        hKinLatStopInside->DrawClone();
        cLatKin->cd(5);
        hKinStopInside->DrawClone();
        cLatKin->cd(6);
        hKinNoReachSil->DrawClone();
        cLatKin->cd(7);
        hTheta3LabReachLatSil->DrawClone();
        cLatKin->cd(8);
        hThetaLabStopInside->DrawClone();
        cLatKin->cd(9);
        effLatSil->Draw("apl");
        cLatKin->cd(10);
        effStopInside->Draw("apl");

        auto cProtonStopped {new TCanvas {"cProtonStopped", "Protons stopped in gas"}};
        cProtonStopped->DivideSquare(4);
        cProtonStopped->cd(1);
        hRangeInGasSin->DrawClone();
        //hRange1->DrawClone("same");
        cProtonStopped->cd(2);
        hKinStopInside->DrawClone();
        gtheo->Draw("l");
        cProtonStopped->cd(3);
        hT3test->DrawClone();
        

        auto cCheckAsymmetry {new TCanvas {"cCheckAsymmetry", "Check Assymetry"}};
        cCheckAsymmetry->DivideSquare(6);
        cCheckAsymmetry->cd(1);
        hPhiCM->DrawClone();
        cCheckAsymmetry->cd(2);
        hPhiLab->DrawClone();
        cCheckAsymmetry->cd(3);
        hRP->DrawClone("colz");
        cCheckAsymmetry->cd(4);
        hRP_ZY->DrawClone("colz");
        cCheckAsymmetry->cd(5);
        hdistanceSPtol0->DrawClone();
        cCheckAsymmetry->cd(6);
        hdistanceSPtor0->DrawClone();

        auto cEfficiency {new TCanvas {"cEfficiency", "Efficiency"}};
        cEfficiency->DivideSquare(6);
        cEfficiency->cd(1);
        effLab->Draw("apl");
        cEfficiency->cd(2);
        effDetectionSilLab->Draw("apl");
        cEfficiency->cd(3);
        effDetectionL1Lab->Draw("apl");
        cEfficiency->cd(4);
        effCM->Draw("apl");
        cEfficiency->cd(5);
        effDetectionSilCM->Draw("apl");
        cEfficiency->cd(6);
        effDetectionL1CM->Draw("apl");
    }
}
#endif
