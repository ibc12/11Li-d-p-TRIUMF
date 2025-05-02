#ifndef heavyDecay_cpp
#define heavyDecay_cpp

#include "ActDecayGenerator.h"
#include "ActSRIM.h"
#include "ActParticle.h"
#include "ActSilSpecs.h"
#include "ActTPCParameters.h"
#include "ActKinematics.h"
#include "ActLine.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "Math/Point3D.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TPolyLine3D.h"


#include "../Histos.h"

ROOT::Math::XYZPointF ComputeLimitPoint(ROOT::Math::XYZVector directionHeavy, ROOT::Math::XYZPoint RP)
{
    // Compute the distance from RP to end of pad plane 0<X<256 0<Y<256 0<Z<256

    // Define the direction given by theta and phi
    ROOT::Math::XYZPointF firstPoint {RP};
    ROOT::Math::XYZPointF secondPoint {RP + 100 * directionHeavy};

    // Let's create a line and then compute the point where x=256 (low angle in heavy, so x will be always the limit)
    auto* line {new ActRoot::Line(firstPoint, secondPoint)};

    auto pointLimit {line->MoveToX(256)};

    double distance {(RP - pointLimit).R()}; // distance in mm
    delete line;
    return pointLimit;
}

void heavyDecay ()
{
    // Get tree from file of transfer dp
    TString fileNameTransfer_dp ("../Outputs/7.5MeV/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0.root");
    auto* filedp {new TFile(fileNameTransfer_dp, "read")};
    if (!filedp || filedp->IsZombie()) 
    {
        std::cerr << "Error: No se pudo abrir el archivo " << fileNameTransfer_dp << std::endl;
        return;
    }
    auto* treeTransfer_dp {static_cast<TTree*>(filedp->Get("SimulationTTreeHeavy"))};
    if (!treeTransfer_dp) 
    {
        std::cerr << "Error: No se encontró 'SimulationTTreeHeavy' en el archivo " << fileNameTransfer_dp << std::endl;
        filedp->ls(); // Listar contenido del archivo para depuración
        return;
    }
    // Get tree from file of transfer dp for light particle
    TString fileNameTransfer_dp_light ("../Outputs/7.5MeV/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0.root");
    auto* filedp_light {new TFile(fileNameTransfer_dp, "read")};
    if (!filedp_light || filedp_light->IsZombie()) 
    {
        std::cerr << "Error: No se pudo abrir el archivo " << fileNameTransfer_dp_light << std::endl;
        return;
    }
    auto* treeTransfer_dp_light {static_cast<TTree*>(filedp_light->Get("SimulationTTree"))};
    if (!treeTransfer_dp_light) 
    {
        std::cerr << "Error: No se encontró 'SimulationTTreeHeavy' en el archivo " << fileNameTransfer_dp_light << std::endl;
        filedp_light->ls(); // Listar contenido del archivo para depuración
        return;
    }
    // Get tree from file of transfer dt
    TString fileNameTransfer_dt ("../Outputs/7.5MeV/2H_3H_TRIUMF_Eex_0.000_nPS_0_pPS_0.root");
    auto* filedt {new TFile(fileNameTransfer_dp, "read")};
    if (!filedt || filedt->IsZombie()) 
    {
        std::cerr << "Error: No se pudo abrir el archivo " << fileNameTransfer_dt << std::endl;
        return;
    }
    auto* treeTransfer_dt {static_cast<TTree*>(filedt->Get("SimulationTTreeHeavy"))};
    if (!treeTransfer_dt) 
    {
        std::cerr << "Error: No se encontró 'SimulationTTreeHeavy' en el archivo " << fileNameTransfer_dt << std::endl;
        filedt->ls(); // Listar contenido del archivo para depuración
        return;
    }
    // Get tree from file of inelastic
    TString fileNameInelastic ("../Outputs/7.5MeV/2H_2H_TRIUMF_Eex_0.818_nPS_2_pPS_0.root");
    auto* fileInelastic {new TFile(fileNameInelastic, "read")};
    if (!fileInelastic || fileInelastic->IsZombie()) 
    {
        std::cerr << "Error: No se pudo abrir el archivo " << fileNameInelastic << std::endl;
        return;
    }
    auto* treeInelastic {static_cast<TTree*>(fileInelastic->Get("SimulationTTreeHeavy"))};
    if (!treeInelastic) 
    {
        std::cerr << "Error: No se encontró 'SimulationTTreeHeavy' en el archivo " << fileNameInelastic << std::endl;
        fileInelastic->ls(); // Listar contenido del archivo para depuración
        return;
    }
    // Get tree from file of elastic
    TString fileNameElastic ("../Outputs/7.5MeV/2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0.root");
    auto* fileElastic {new TFile(fileNameElastic, "read")};
    if (!fileElastic || fileElastic->IsZombie()) 
    {
        std::cerr << "Error: No se pudo abrir el archivo " << fileNameInelastic << std::endl;
        return;
    }
    auto* treeElastic {static_cast<TTree*>(fileElastic->Get("SimulationTTreeHeavy"))};
    if (!treeElastic) 
    {
        std::cerr << "Error: No se encontró 'SimulationTTreeHeavy' en el archivo " << fileNameInelastic << std::endl;
        fileElastic->ls(); // Listar contenido del archivo para depuración
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
    srim->ReadTable("9LiGas", "../SRIM files/9Li_900mb_CF4_90-10.txt");
    srim->ReadTable("11LiGas", "../SRIM files/11Li_900mb_CF4_90-10.txt");
    srim->ReadTable("12LiGas", "../SRIM files/12Li_900mb_CF4_90-10.txt");
    srim->ReadTable("9LiInSil", "../SRIM files/9Li_silicon.txt");
    srim->ReadTable("11LiInSil", "../SRIM files/11Li_silicon.txt");
    srim->ReadTable("12LiInSil", "../SRIM files/12Li_silicon.txt");

    // Particles for decays
    ActPhysics::Particle Li12 = ActPhysics::Particle(3, 12);
    // Li12.SetEx(0.43525); // it says that not enough energy to decay to that state
    ActPhysics::Particle Li11 = ActPhysics::Particle(3, 11);
    // Li11.SetEx(1.26642); // it says thet not enough energy to decay to that state
    ActPhysics::Particle Li11excited = ActPhysics::Particle(3, 11);
    Li11excited.SetEx(1.26642);
    ActPhysics::Particle Li9 = ActPhysics::Particle(3, 9);
    ActPhysics::Particle Li10 = ActPhysics::Particle(3, 10);
    ActPhysics::Particle neutron = ActPhysics::Particle(0, 1);
    ActPhysics::Particle deuterium = ActPhysics::Particle(1, 2);

    // Decay generator
    auto decayGen11Li {new ActSim::DecayGenerator(Li12, Li11, neutron)};
    auto decayGen9Li_from11Li {new ActSim::DecayGenerator(Li11excited, Li9, neutron, neutron)};
    auto decayGen9Li_from10Li {new ActSim::DecayGenerator(Li10, Li9, neutron)};

    // kinematics
    auto kinElastic {new ActPhysics::Kinematics {"11Li(d,d)@82.5|0"}};
    auto kinTransfer_dp {new ActPhysics::Kinematics {"11Li(d,p)@82.5|0"}};
    auto kinTransfer_dt {new ActPhysics::Kinematics {"11Li(d,t)@82.5|0"}};
    auto kinInelastic {new ActPhysics::Kinematics {"11Li(d,d)@82.5|1"}}; // Ex centered at (1.26642 + 0.36928) / 2, and is a gaussian distribution of sigma 0.2, so will not be accurate

    // Read TTree in a for loop and see decay for each event for TRANSFER
    double theta4Lab {};
    double phi4Lab {};
    double T4Lab {};
    ROOT::Math::XYZPoint* RP {};
    double theta3Lab {};
    double phi3Lab {};
    double T3Lab {};

    treeTransfer_dp_light->SetBranchAddress("theta3Lab", &theta3Lab);
    treeTransfer_dp_light->SetBranchAddress("phi3CM", &phi3Lab);
    treeTransfer_dp_light->SetBranchAddress("EVertex", &T3Lab);

    treeTransfer_dp->SetBranchAddress("theta4Lab", &theta4Lab);
    treeTransfer_dp->SetBranchAddress("phi4Lab", &phi4Lab);
    treeTransfer_dp->SetBranchAddress("T4Lab", &T4Lab);
    treeTransfer_dp->SetBranchAddress("RP", &RP);

    treeTransfer_dt->SetBranchAddress("theta4Lab", &theta4Lab);
    treeTransfer_dt->SetBranchAddress("phi4Lab", &phi4Lab);
    treeTransfer_dt->SetBranchAddress("T4Lab", &T4Lab);
    treeTransfer_dt->SetBranchAddress("RP", &RP);

    treeInelastic->SetBranchAddress("theta4Lab", &theta4Lab);
    treeInelastic->SetBranchAddress("phi4Lab", &phi4Lab);
    treeInelastic->SetBranchAddress("T4Lab", &T4Lab);
    treeInelastic->SetBranchAddress("RP", &RP);

    treeElastic->SetBranchAddress("theta4Lab", &theta4Lab);
    treeElastic->SetBranchAddress("phi4Lab", &phi4Lab);
    treeElastic->SetBranchAddress("T4Lab", &T4Lab);
    treeElastic->SetBranchAddress("RP", &RP);

    // Read TTree in a for loop and see decay for each event for INELASTIC 

    // Histograms
    auto hSPf0 {Histos::SP.GetHistogram()};
    hSPf0->SetTitle("SP for f0");
    auto hSPf2 {Histos::SP.GetHistogram()};
    hSPf2->SetTitle("SP for f2");
    auto hTheta12Li {Histos::ThetaLabHeavy.GetHistogram()};
    hTheta12Li->SetTitle("ThetaLab for 12Li");
    auto hTheta11Li {Histos::ThetaLabHeavy.GetHistogram()};
    hTheta11Li->SetTitle("ThetaLab for 11Li");
    auto hkin11Li {Histos::KinHeavy.GetHistogram()};
    hkin11Li->SetTitle("Heavy kinematics for11Li");
    auto hkin9Li {Histos::KinHeavy.GetHistogram()};
    hkin9Li->SetTitle("Heavy kinematics for 9Li");
    auto hkin {Histos::KinHeavy.GetHistogram()};
    hkin->SetTitle("Heavy kinematics for all cases");
    auto hkin12Li {Histos::KinHeavy.GetHistogram()};
    hkin12Li->SetTitle("Heavy kinematics for 12Li");
    auto hkin10Li {Histos::KinHeavy.GetHistogram()};
    hkin10Li->SetTitle("Heavy kinematics for 10Li");
    auto hkin11LiInelastic {Histos::KinHeavy.GetHistogram()};
    hkin11LiInelastic->SetTitle("Heavy kinematics for 11Li inelastic");
    auto hkin11LiElastic {Histos::KinHeavy.GetHistogram()};
    hkin11LiElastic->SetTitle("Heavy kinematics for 11Li elastic");
    // PID's
    auto hPIDtransfer_dp {Histos::PID.GetHistogram()};
    hPIDtransfer_dp->SetTitle("PID for transfer dp");
    auto hPIDtransfer_dt {Histos::PID.GetHistogram()};
    hPIDtransfer_dt->SetTitle("PID for transfer dt");
    auto hPIDinelastic {Histos::PID.GetHistogram()};
    hPIDinelastic->SetTitle("PID for inelastic");
    auto hPIDelastic {Histos::PID.GetHistogram()};
    hPIDelastic->SetTitle("PID for elastic");
    auto hPIDtransferLength_dp {Histos::PIDlength.GetHistogram()};
    hPIDtransferLength_dp->SetTitle("PID for transfer length dp");
    auto hPIDtransferLength_dt {Histos::PIDlength.GetHistogram()};
    hPIDtransferLength_dt->SetTitle("PID for transfer length dt");
    auto hPIDinelasticLength {Histos::PIDlength.GetHistogram()};
    hPIDinelasticLength->SetTitle("PID for inelastic length");
    auto hPIDelasticLength {Histos::PIDlength.GetHistogram()};
    hPIDelasticLength->SetTitle("PID for elastic length");
    auto hPIDAll {Histos::PID.GetHistogram()};
    hPIDAll->SetTitle("PID for all");
    auto hPIDAllLength {Histos::PIDlength.GetHistogram()};
    hPIDAllLength->SetTitle("PID for all length");

    // Histos for debuging
    auto hTheta11LiDebug {Histos::ThetaLabHeavy.GetHistogram()};
    hTheta11LiDebug->SetTitle("ThetaLab for 11Li");
    auto hkin11LiDebug {Histos::KinHeavy.GetHistogram()};
    hkin11LiDebug->SetTitle("Heavy kinematics for 11Li");
    auto hDistances11Li {Histos::DistanceSP.GetHistogram()};
    hDistances11Li->SetTitle("Distances");
    // Debug broad blob
    auto hTheta4LabBroad {Histos::ThetaLabHeavy.GetHistogram()};
    hTheta4LabBroad->SetTitle("ThetaLab for 11Li");
    auto hT4LabBroad {Histos::T4Lab.GetHistogram()};
    hT4LabBroad->SetTitle("T4 for 11Li");
    auto hkin11LiBroad {Histos::KinHeavy.GetHistogram()};
    hkin11LiBroad->SetTitle("Heavy kinematics for 11Li");
    auto hVertexX {Histos::DistanceSP.GetHistogram()};
    hVertexX->SetTitle("Vertex X");
    auto hLossGasvsT4Lab {Histos::PID.GetHistogram()};
    hLossGasvsT4Lab->SetTitle("Loss in gas vs T4");
    // Debug broad blob with light particle information
    auto hTheta3LabBroad {Histos::ThetaLab.GetHistogram()};
    hTheta3LabBroad->SetTitle("ThetaLab for 1H");
    auto hT3LabBroad {Histos::T4Lab.GetHistogram()};
    hT3LabBroad->SetTitle("T3 for 1H");
    auto hkin1HBroad {Histos::Kin.GetHistogram()};
    hkin1HBroad->SetTitle("Kinematics for 1H");
    
    // Loop over the tree
    // for(int r = 0; r < 10; r++)
    // {
    for(int i = 0; i < treeTransfer_dp->GetEntries(); i++)
    {
        treeTransfer_dp->GetEntry(i);
        treeTransfer_dp_light->GetEntry(i);
        if(T4Lab > 60.)
        {
            continue;
        }
        // Decay
        hkin12Li->Fill(theta4Lab * TMath::RadToDeg(), T4Lab);
        decayGen11Li->SetDecay(T4Lab, theta4Lab, phi4Lab);
        decayGen11Li->Generate();
        // Get the decay products
        auto LorentzVectorLi11 {decayGen11Li->GetLorentzVector(0)};
        auto thetaLi11 {LorentzVectorLi11->Theta()};
        auto phiLi11 {LorentzVectorLi11->Phi()};
        auto TLi11 {LorentzVectorLi11->E() - LorentzVectorLi11->M()};
        
        ROOT::Math::XYZVector directionLi11 {TMath::Cos(thetaLi11), TMath::Sin(thetaLi11) * TMath::Sin(phiLi11),
                              TMath::Sin(thetaLi11) * TMath::Cos(phiLi11)};

        // Check hit for the 11Li 
        int silIndex = -1;
        ROOT::Math::XYZPoint silPoint;
        std::string layerHit;
        for(auto layer : silLayers)
        {
            std::tie(silIndex, silPoint) = sils->FindSPInLayer(layer, *RP, directionLi11);
            if(silIndex != -1)
            {
                layerHit = layer;
                break;
            }                
        }
        if(silIndex == -1)
        {
            continue; // If a silicon is not reached, don't continue with punchthough calculation
        }
        // Calculation of DeltaE-E
        // auto rangeInitial {srim->EvalRange("11LiGas", TLi11)};
        // First, slow inside detector volume
        auto limitPointGas {ComputeLimitPoint(directionLi11, *RP)};
        double distanceInside {(*RP - limitPointGas).R()};
        auto energyAfterInside {srim->SlowWithStraggling("11LiGas", TLi11, distanceInside)};
        // auto DeltaE {TLi11 - energyAfterInside};
        // Second, slow in gas before silicon
        auto distanceInterGas {(limitPointGas - silPoint).R()};
        auto energyAfterInterGas {srim->SlowWithStraggling("11LiGas", energyAfterInside, distanceInterGas)};
        auto DeltaE {TLi11 - energyAfterInside};
        // Finally, slow in silicon
        auto energyAfterSil {srim->SlowWithStraggling("11LiInSil", energyAfterInterGas, sils->GetLayer(layerHit).GetUnit().GetThickness(), thetaLi11)};
        auto eLoss {energyAfterInterGas - energyAfterSil};
        if(energyAfterSil > 0)
        {
            std::cout <<" Energy in silicon: " << energyAfterSil << " Angle: " << thetaLi11 * TMath::RadToDeg() << " Layer: " << layerHit.size() <<'\n';
            hTheta11LiDebug->Fill(thetaLi11 * TMath::RadToDeg());
            hkin11LiDebug->Fill(thetaLi11 * TMath::RadToDeg(), TLi11);
        }
        else
        {
            hPIDtransfer_dp->Fill(eLoss, DeltaE);
            hPIDtransferLength_dp->Fill(eLoss, DeltaE / distanceInside);
            hDistances11Li->Fill(RP->X());
            hPIDAll->Fill(eLoss, DeltaE);
            hPIDAllLength->Fill(eLoss, DeltaE / distanceInside);
        }
        

        // Fill SP histos
        if(layerHit == "f0")
        {
            hSPf0->Fill(silPoint.Y(), silPoint.Z());
        }
        if(layerHit == "f2")
        {
            hSPf2->Fill(silPoint.Y(), silPoint.Z());
        }
        // Fill rest histos
        hTheta12Li->Fill(theta4Lab * TMath::RadToDeg());
        hTheta11Li->Fill(thetaLi11 * TMath::RadToDeg());
        if(layerHit == "f2")
        {
            hkin11Li->Fill(thetaLi11 * TMath::RadToDeg(), TLi11);
            hkin->Fill(thetaLi11 * TMath::RadToDeg(), TLi11);
        }
        // fill broad histos
        hTheta4LabBroad->Fill(theta4Lab * TMath::RadToDeg());
        hT4LabBroad->Fill(T4Lab);
        hkin11LiBroad->Fill(thetaLi11 * TMath::RadToDeg(), TLi11);
        hVertexX->Fill(RP->X());
        hLossGasvsT4Lab->Fill(T4Lab, DeltaE);
        // now for light particle
        hTheta3LabBroad->Fill(theta3Lab);
        hT3LabBroad->Fill(T3Lab);
        hkin1HBroad->Fill(theta3Lab, T3Lab);
        

    }
// }
    // Loop over the tree  dt
    for(int i = 0; i < treeTransfer_dt->GetEntries(); i++)
    {
        treeTransfer_dt->GetEntry(i);
        // Decay
        hkin10Li->Fill(theta4Lab * TMath::RadToDeg(), T4Lab);
        decayGen9Li_from10Li->SetDecay(T4Lab, theta4Lab, phi4Lab);
        decayGen9Li_from10Li->Generate();
        // Get the decay products
        auto LorentzVectorLi9 {decayGen9Li_from10Li->GetLorentzVector(0)};
        auto thetaLi9 {LorentzVectorLi9->Theta()};
        auto phiLi9 {LorentzVectorLi9->Phi()};
        auto TLi9 {LorentzVectorLi9->E() - LorentzVectorLi9->M()};
        
        ROOT::Math::XYZVector directionLi11 {TMath::Cos(thetaLi9), TMath::Sin(thetaLi9) * TMath::Sin(phiLi9),
                              TMath::Sin(thetaLi9) * TMath::Cos(phiLi9)};

        // Check hit for the 11Li 
        int silIndex = -1;
        ROOT::Math::XYZPoint silPoint;
        std::string layerHit;
        for(auto layer : silLayers)
        {
            std::tie(silIndex, silPoint) = sils->FindSPInLayer(layer, *RP, directionLi11);
            if(silIndex != -1)
            {
                layerHit = layer;
                break;
            }                
        }
        if(silIndex == -1)
        {
            continue; // If a silicon is not reached, don't continue with punchthough calculation
        }
        // Calculation of DeltaE-E
        // auto rangeInitial {srim->EvalRange("11LiGas", TLi11)};
        // First, slow inside detector volume
        auto limitPointGas {ComputeLimitPoint(directionLi11, *RP)};
        double distanceInside {(*RP - limitPointGas).R()};
        auto energyAfterInside {srim->SlowWithStraggling("11LiGas", TLi9, distanceInside)};
        // auto DeltaE {TLi11 - energyAfterInside};
        // Second, slow in gas before silicon
        auto distanceInterGas {(limitPointGas - silPoint).R()};
        auto energyAfterInterGas {srim->SlowWithStraggling("11LiGas", energyAfterInside, distanceInterGas)};
        auto DeltaE {TLi9 - energyAfterInside};
        // Finally, slow in silicon
        auto energyAfterSil {srim->SlowWithStraggling("11LiInSil", energyAfterInterGas, sils->GetLayer(layerHit).GetUnit().GetThickness(), thetaLi9)};
        auto eLoss {energyAfterInterGas - energyAfterSil};
        if(energyAfterSil > 0)
        {
            std::cout <<" Energy in silicon: " << energyAfterSil << " Angle: " << thetaLi9 * TMath::RadToDeg() << " Layer: " << layerHit.size() <<'\n';
        }
        else
        {
            hPIDtransfer_dt->Fill(eLoss, DeltaE);
            hPIDtransferLength_dt->Fill(eLoss, DeltaE / distanceInside);
            hPIDAll->Fill(eLoss, DeltaE);
            hPIDAllLength->Fill(eLoss, DeltaE / distanceInside);
        }
        

        // Fill SP histos
        if(layerHit == "f0")
        {
            hSPf0->Fill(silPoint.Y(), silPoint.Z());
        }
        if(layerHit == "f2")
        {
            hSPf2->Fill(silPoint.Y(), silPoint.Z());
        }
        // Fill rest histos
        // hTheta10Li->Fill(theta4Lab * TMath::RadToDeg());
        // hTheta9Li->Fill(thetaLi11 * TMath::RadToDeg());
        if(layerHit == "f2")
        {
            hkin9Li->Fill(thetaLi9 * TMath::RadToDeg(), TLi9);
            hkin->Fill(thetaLi9 * TMath::RadToDeg(), TLi9);
        }
        // PID
        

    }
    // Loop over the tree of inelastic
    for(int i = 0; i < treeInelastic->GetEntries(); i++)
    {
        treeInelastic->GetEntry(i);
        // There isno decay from 11Li to 9Li, phase space directly
        // hkin11LiInelastic->Fill(theta4Lab * TMath::RadToDeg(), T4Lab);
        ROOT::Math::XYZVector directionLi9 {TMath::Cos(theta4Lab), TMath::Sin(theta4Lab) * TMath::Sin(phi4Lab),
                              TMath::Sin(theta4Lab) * TMath::Cos(phi4Lab)};

        // Check hit for the 11Li 
        int silIndex = -1;
        ROOT::Math::XYZPoint silPoint;
        std::string layerHit;
        for(auto layer : silLayers)
        {
            std::tie(silIndex, silPoint) = sils->FindSPInLayer(layer, *RP, directionLi9);
            if(silIndex != -1)
            {
                layerHit = layer;
                break;
            }                
        }
        // Calculation of DeltaE-E
        auto rangeInitial {srim->EvalRange("9LiGas", T4Lab)};
        // First, slow inside detector volume
        auto limitPointGas {ComputeLimitPoint(directionLi9, *RP)};
        double distanceInside {(*RP - limitPointGas).R()};
        auto energyAfterInside {srim->SlowWithStraggling("9LiGas", T4Lab, distanceInside)};
        auto DeltaE {T4Lab - energyAfterInside};
        // Second, slow in gas before silicon
        auto distanceInterGas {(limitPointGas - silPoint).R()};
        auto energyAfterInterGas {srim->SlowWithStraggling("9LiGas", energyAfterInside, distanceInterGas)};
        // Finally, slow in silicon
        auto energyAfterSil {srim->SlowWithStraggling("9LiInSil", energyAfterInterGas, sils->GetLayer(layerHit).GetUnit().GetThickness())};
        auto eLoss {energyAfterInterGas - energyAfterSil};


        // Fill all the histos
        if(layerHit == "f0")
        {
            hSPf0->Fill(silPoint.Y(), silPoint.Z());
        }
        if(layerHit == "f2")
        {
            hSPf2->Fill(silPoint.Y(), silPoint.Z());
        }

        if(layerHit == "f2")
        {
            hkin9Li->Fill(theta4Lab * TMath::RadToDeg(), T4Lab);
            hkin->Fill(theta4Lab * TMath::RadToDeg(), T4Lab);
        }
        // PID
        hPIDinelastic->Fill(eLoss, DeltaE);
        hPIDinelasticLength->Fill(eLoss, DeltaE / distanceInside);
        hPIDAll->Fill(eLoss, DeltaE);
        hPIDAllLength->Fill(eLoss, DeltaE / distanceInside);
    }
    // Loop over the elastic events
    for(int i = 0; i < treeElastic->GetEntries(); i++)
    {
        treeElastic->GetEntry(i);
        // For elastic there is no decay
        ROOT::Math::XYZVector directionLi11 {TMath::Cos(theta4Lab), TMath::Sin(theta4Lab) * TMath::Sin(phi4Lab),
                              TMath::Sin(theta4Lab) * TMath::Cos(phi4Lab)};

        // Check hit for the 11Li 
        int silIndex = -1;
        ROOT::Math::XYZPoint silPoint;
        std::string layerHit;
        for(auto layer : silLayers)
        {
            std::tie(silIndex, silPoint) = sils->FindSPInLayer(layer, *RP, directionLi11);
            if(silIndex != -1)
            {
                layerHit = layer;
                break;
            }                
        }
        // Calculation of DeltaE-E
        auto rangeInitial {srim->EvalRange("11LiGas", T4Lab)};
        // First, slow inside detector volume
        auto limitPointGas {ComputeLimitPoint(directionLi11, *RP)};
        double distanceInside {(*RP - limitPointGas).R()};
        auto energyAfterInside {srim->SlowWithStraggling("11LiGas", T4Lab, distanceInside)};
        auto DeltaE {T4Lab - energyAfterInside};
        // Second, slow in gas before silicon
        auto distanceInterGas {(limitPointGas - silPoint).R()};
        auto energyAfterInterGas {srim->SlowWithStraggling("11LiGas", energyAfterInside, distanceInterGas)};
        // Finally, slow in silicon
        auto energyAfterSil {srim->SlowWithStraggling("11LiInSil", energyAfterInterGas, sils->GetLayer(layerHit).GetUnit().GetThickness())};
        auto eLoss {energyAfterInterGas - energyAfterSil};


        // Fill all the histos
        if(layerHit == "f0")
        {
            hSPf0->Fill(silPoint.Y(), silPoint.Z());
        }
        if(layerHit == "f2")
        {
            hSPf2->Fill(silPoint.Y(), silPoint.Z());
        }

        if(layerHit == "f2")
        {
            hkin11LiElastic->Fill(theta4Lab * TMath::RadToDeg(), T4Lab);
            hkin->Fill(theta4Lab * TMath::RadToDeg(), T4Lab);
        }
        // PID
        hPIDelastic->Fill(eLoss, DeltaE);
        hPIDelasticLength->Fill(eLoss, DeltaE / distanceInside);
        hPIDAll->Fill(eLoss, DeltaE);
        hPIDAllLength->Fill(eLoss, DeltaE / distanceInside);
    }


    // Plot
    auto c1 {new TCanvas("c1", "c1", 800, 600)};
    c1->DivideSquare(6);
    c1->cd(1);
    hkin12Li->DrawClone("colz");
    //hTheta12Li->DrawClone();
    c1->cd(2);
    hTheta11Li->DrawClone();
    c1->cd(3);
    hkin->DrawClone("colz");
    c1->cd(4);
    hkin11Li->DrawClone("colz");
    // kinElastic->GetKinematicLine4()->Draw("l");
    c1->cd(5);
    hkin9Li->DrawClone("colz");
    c1->cd(6);
    hDistances11Li->DrawClone();


    auto* cSP {new TCanvas {"cSP", "Sil Points"}};
    cSP->DivideSquare(3);
    cSP->cd(1);
    hSPf0->DrawClone("colz");
    sils->GetLayer("f0").GetSilMatrix()->Draw();
    cSP->cd(2);
    hSPf2->DrawClone("colz");
    sils->GetLayer("f2").GetSilMatrix()->Draw();

    auto* ckins {new TCanvas {"ckins", "Kinematics"}};
    ckins->DivideSquare(6);
    ckins->cd(1);
    hkin12Li->DrawClone("colz");
    kinTransfer_dp->GetKinematicLine4()->Draw("l");
    ckins->cd(4);
    hkin11Li->DrawClone("colz");
    ckins->cd(2);
    hkin11LiInelastic->DrawClone("colz");
    kinInelastic->GetKinematicLine4()->Draw("l");
    ckins->cd(5);
    hkin9Li->DrawClone("colz");
    ckins->cd(3);
    hkin11LiElastic->DrawClone("colz");
    kinElastic->GetKinematicLine4()->Draw("l");
    ckins->cd(6);
    hkin->DrawClone("colz");

    auto* cPIDs {new TCanvas {"cPIDs", "PID"}};
    cPIDs->DivideSquare(8);
    cPIDs->cd(1);
    hPIDtransfer_dp->DrawClone("colz");
    cPIDs->cd(2);
    hPIDinelastic->DrawClone("colz");
    cPIDs->cd(3);
    hPIDelastic->DrawClone("colz");
    cPIDs->cd(4);
    hPIDtransfer_dt->DrawClone("colz");
    cPIDs->cd(5);
    hPIDtransferLength_dp->DrawClone("colz");
    cPIDs->cd(6);
    hPIDinelasticLength->DrawClone("colz");
    cPIDs->cd(7);
    hPIDelasticLength->DrawClone("colz");
    cPIDs->cd(8);
    hPIDtransferLength_dt->DrawClone("colz");

    auto* cPIDAll {new TCanvas {"cPIDAll", "PID All"}};
    cPIDAll->DivideSquare(2);
    cPIDAll->cd(1);
    hPIDAll->DrawClone("colz");
    cPIDAll->cd(2);
    hPIDAllLength->DrawClone("colz");
    

    auto* cBroadBlob {new TCanvas {"cBroadBlob", "Broad Blob"}};
    cBroadBlob->DivideSquare(9);
    cBroadBlob->cd(1);
    hTheta4LabBroad->DrawClone();
    cBroadBlob->cd(2);
    hT4LabBroad->DrawClone();
    cBroadBlob->cd(3);
    hkin11LiBroad->DrawClone("colz");
    cBroadBlob->cd(4);
    hVertexX->DrawClone();
    cBroadBlob->cd(5);
    hLossGasvsT4Lab->DrawClone("colz");
    cBroadBlob->cd(6);
    hPIDtransferLength_dp->DrawClone("colz");
    cBroadBlob->cd(7);
    hTheta3LabBroad->DrawClone();
    cBroadBlob->cd(8);
    hT3LabBroad->DrawClone();
    cBroadBlob->cd(9);
    hkin1HBroad->DrawClone("colz");

    auto* cBroadLight {new TCanvas {"cBroadLight", "Broad Light"}};
    cBroadLight->DivideSquare(3);
    cBroadLight->cd(1);
    hTheta3LabBroad->DrawClone();
    cBroadLight->cd(2);
    hT3LabBroad->DrawClone();   
    cBroadLight->cd(3);
    hkin1HBroad->DrawClone("colz");
}
#endif