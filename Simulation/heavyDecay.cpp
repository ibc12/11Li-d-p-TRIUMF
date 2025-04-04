#ifndef heavyDecay_cpp
#define heavyDecay_cpp

#include "ActDecayGenerator.h"
#include "ActSRIM.h"
#include "ActParticle.h"
#include "ActSilSpecs.h"
#include "ActTPCParameters.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "Math/Point3D.h"
#include "TH1D.h"
#include "TH2D.h"

#include "../Histos.h"

void heavyDecay ()
{
    // Get tree from file
    TString fileName ("../Outputs/7.5MeV/transfer_TRIUMF_Eex_0.000_nPS_0_pPS_0.root");
    auto* file {new TFile(fileName, "read")};
    if (!file || file->IsZombie()) 
    {
        std::cerr << "Error: No se pudo abrir el archivo " << fileName << std::endl;
        return;
    }
    auto* tree {static_cast<TTree*>(file->Get("SimulationTTreeHeavy"))};
    if (!tree) 
    {
        std::cerr << "Error: No se encontró 'SimulationTTreeHeavy' en el archivo " << fileName << std::endl;
        file->ls(); // Listar contenido del archivo para depuración
        return;
    }

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


        
    std::cout << "Sils Z centred at : " << tpc.Z() / 2 << " mm" << '\n';
    sils->DrawGeo();

    // SRIM
    auto* srim {new ActPhysics::SRIM};
    srim->ReadTable("11LiGas", "../SRIM files/11Li_900mb_CF4_90-10.txt");
    srim->ReadTable("12LiGas", "../SRIM files/12Li_900mb_CF4_90-10.txt");
    srim->ReadTable("12LiInSil", "../SRIM files/12Li_silicon.txt");

    // Particles for decays
    ActPhysics::Particle Li12 = ActPhysics::Particle(3, 12);
    ActPhysics::Particle Li11 = ActPhysics::Particle(3, 11);
    ActPhysics::Particle Li9 = ActPhysics::Particle(3, 9);
    ActPhysics::Particle neutron = ActPhysics::Particle(0, 1);
    ActPhysics::Particle deuterium = ActPhysics::Particle(1, 2);

    // Decay generator
    auto decayGen11Li {new ActSim::DecayGenerator(Li12, Li11, neutron)};
    auto decayGen9Li {new ActSim::DecayGenerator(Li11, Li9, deuterium)};

    // Read TTree in a for loop and see decay for each event
    double theta4Lab {};
    double phi4Lab {};
    double T4Lab {};
    ROOT::Math::XYZPoint RP {};
    tree->SetBranchAddress("theta4Lab", &theta4Lab);
    tree->SetBranchAddress("phi4Lab", &phi4Lab);
    tree->SetBranchAddress("T4Lab", &T4Lab);
    tree->SetBranchAddress("RP", &RP);

    // Histograms
    auto hSPf0 {Histos::SP.GetHistogram()};
    hSPf0->SetTitle("SP for f0");
    auto hSPf2 {Histos::SP.GetHistogram()};
    hSPf2->SetTitle("SP for f2");
    auto hTheta12Li {Histos::ThetaLabHeavy.GetHistogram()};
    hTheta12Li->SetTitle("ThetaLab for 12Li");
    auto hTheta11Li {Histos::ThetaLabHeavy.GetHistogram()};
    hTheta11Li->SetTitle("ThetaLab for 11Li");
    
    // Loop over the tree
    for(int i = 0; i < tree->GetEntries(); i++)
    {
        tree->GetEntry(i);
        // Decay
        decayGen11Li->SetDecay(T4Lab, theta4Lab, phi4Lab);
        decayGen11Li->Generate();
        auto LorentzVectorLi11 {decayGen11Li->GetLorentzVector(0)};
        auto thetaLi11 {LorentzVectorLi11->Theta()};
        auto phiLi11 {LorentzVectorLi11->Phi()};
        auto TLi11 {LorentzVectorLi11->E() - LorentzVectorLi11->M()};
        
        ROOT::Math::XYZVector directionLi11 {TMath::Cos(thetaLi11), TMath::Sin(thetaLi11) * TMath::Sin(phiLi11),
                             TMath::Sin(thetaLi11) * TMath::Cos(phiLi11)};

        int silIndex = -1;
        ROOT::Math::XYZPoint silPoint;
        std::string layer;
        for(auto layer : silLayers)
        {
            std::tie(silIndex, silPoint) = sils->FindSPInLayer(layer, RP, directionLi11);
            if(silIndex != -1)
            {
                layer = layer;
                break;
            }                
        }

        if(layer == "f0")
        {
            hSPf0->Fill(silPoint.Y(), silPoint.Z());
        }
        if(layer == "f2")
        {
            hSPf2->Fill(silPoint.Y(), silPoint.Z());
        }


        // Fill histos
        hTheta12Li->Fill(theta4Lab * TMath::RadToDeg());
        hTheta11Li->Fill(thetaLi11 * TMath::RadToDeg());
    }

    // Plot
    auto c1 {new TCanvas("c1", "c1", 800, 600)};
    c1->DivideSquare(4);
    c1->cd(1);
    hTheta12Li->DrawClone();
    c1->cd(2);
    hTheta11Li->DrawClone();

    auto* cSP {new TCanvas {"cSP", "Sil Points"}};
    cSP->DivideSquare(3);
    cSP->cd(1);
    hSPf0->DrawClone("colz");
    sils->GetLayer("f0").GetSilMatrix()->Draw();
    cSP->cd(2);
    hSPf2->DrawClone("colz");
    sils->GetLayer("f2").GetSilMatrix()->Draw();


}

#endif