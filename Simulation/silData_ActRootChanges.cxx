#include "ROOT/RDataFrame.hxx"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"


#include "ActColors.h"
#include "ActSilData.h"
#include "ActTPCData.h"
#include "ActSilData.h"
#include "ActModularData.h"
#include "ActLine.h"
#include "ActCluster.h"

void silData_ActRootChanges()
{
    // Read SilData from Output Ttrees and work with them
    auto file = TFile::Open("../Outputs/7.5MeV/2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_silicons_reverse.root", "READ");
    if (!file || !file->IsOpen())
    {
        std::cout << "Error opening file!" << std::endl;
        return;
    }
    auto tree = static_cast<TTree *>(file->Get("SimulationTTreeNoCuts"));
    if (!tree)
    {
        std::cout << "Error retrieving tree!" << std::endl;
        return;
    }

    // Prepare Input
    ActRoot::SilData *silData {};
    double theta3Lab{};
    double theta4Lab{};
    double phi3Lab{};
    double phi4Lab{};
    double T3Lab{};
    double T4Lab{};
    ROOT::Math::XYZPoint *rp{};
    tree->SetBranchAddress("silData", &silData);
    tree->SetBranchAddress("theta3Lab", &theta3Lab);
    tree->SetBranchAddress("theta4Lab", &theta4Lab);
    tree->SetBranchAddress("phi3CM", &phi3Lab);
    tree->SetBranchAddress("phi4CM", &phi4Lab);
    tree->SetBranchAddress("T3Lab", &T3Lab);
    tree->SetBranchAddress("T4Lab", &T4Lab);
    tree->SetBranchAddress("RP", &rp);

    // Prepare output variables
    auto TPCDataOut{new ActRoot::TPCData()};
    auto silDataOut{new ActRoot::SilData()};
    auto modularDataOut{new ActRoot::ModularData()};
    // Prepare outpur file
    auto outFileTPC{TFile::Open("../DebugOutputs/TPCData.root", "RECREATE")};
    auto outFileSils{TFile::Open("../DebugOutputs/test_merger_0001.root", "RECREATE")};
    auto outTreeTPC{new TTree("GETTree", "TPC Data Tree")};
    auto outTreeSils{new TTree("VXITree", "Silicon Data Tree")};
    outTreeTPC->Branch("TPCData", &TPCDataOut);
    outTreeSils->Branch("SilData", &silDataOut);
    outTreeSils->Branch("ModularData", &modularDataOut);
    // Fill the trees
    int nEntries = tree->GetEntries();
    std::cout << BOLDGREEN;
    const int percentPrint{5};
    int step{nEntries / (100 / percentPrint)};
    int nextPrint{step};
    int percent{};
    for (int i = 0; i < nEntries; i++)
    {
        if (i >= nextPrint)
        {
            percent = 100 * i / nEntries;
            std::cout << "\r" << std::string(percent / percentPrint, '|') << percent << "%";
            std::cout.flush();
            nextPrint += step;
        }
        tree->GetEntry(i);
        // Fill SilData
        silDataOut = silData;
        // Fill TPCData
        ROOT::Math::XYZVectorF dirLight{
                        static_cast<float>(TMath::Cos(theta3Lab)),
                        static_cast<float>(TMath::Sin(theta3Lab) * TMath::Sin(phi3Lab)),
                        static_cast<float>(TMath::Sin(theta3Lab) * TMath::Cos(phi3Lab))
                    };
        ROOT::Math::XYZPointF pointLight{rp->X() + dirLight.X() * 2, 
                                         rp->Y() + dirLight.Y() * 2, 
                                         rp->Z() + dirLight.Z() * 2};
        ROOT::Math::XYZPointF rpF{rp->X(), rp->Y(), rp->Z()};
        ActRoot::Line lineLight(rpF, pointLight);
        ActRoot::Cluster clusterLight;
        clusterLight.SetLine(lineLight);
        TPCDataOut->fClusters.push_back(clusterLight);
        TPCDataOut->fRPs.push_back(rpF);
        ROOT::Math::XYZVectorF dirHeavy{
                        static_cast<float>(TMath::Cos(theta4Lab)),
                        static_cast<float>(TMath::Sin(theta4Lab) * TMath::Sin(phi4Lab)),
                        static_cast<float>(TMath::Sin(theta4Lab) * TMath::Cos(phi4Lab))
                    };
        ROOT::Math::XYZPointF pointHeavy{rp->X() + dirHeavy.X() * 2, 
                                         rp->Y() + dirHeavy.Y() * 2, 
                                         rp->Z() + dirHeavy.Z() * 2};
        ActRoot::Line lineHeavy(rpF, pointHeavy);
        ActRoot::Cluster clusterHeavy;
        clusterHeavy.SetLine(lineHeavy);
        TPCDataOut->fClusters.push_back(clusterHeavy);

        ROOT::Math::XYZVectorF dirBeam{
                        static_cast<float>(1),
                        static_cast<float>(0),
                        static_cast<float>(0)
                    };
        ROOT::Math::XYZPointF pointBeam{rp->X() + dirBeam.X() * 2, 
                             rp->Y() + dirBeam.Y() * 2, 
                             rp->Z() + dirBeam.Z() * 2};
        ActRoot::Line lineBeam(rpF, pointBeam);
        ActRoot::Cluster clusterBeam;
        clusterBeam.SetLine(lineBeam);
        clusterBeam.SetBeamLike(true);
        TPCDataOut->fClusters.push_back(clusterBeam);
        // Fill ModularData
        bool LateralAndHeavy {silData && (silData->fSiN["l0"].size() != 0 || silData->fSiN["r0"].size() != 0)
                                     && (silData->fSiN["f2"].size() != 0 || silData->fSiN["f3"].size() != 0)};
        bool FrontAndHeavy {silData && (silData->fSiN["f0"].size() != 0 || silData->fSiN["f1"].size())
                                     && (silData->fSiN["f2"].size() != 0 || silData->fSiN["f3"].size() != 0)};
        bool BothFront {silData && (silData->fSiN["l0"].size() == 0 || silData->fSiN["r0"].size() == 0)
                                     && (silData->fSiN["f0"].size() != 0 || silData->fSiN["f0"].size() != 0)
                                     && (silData->fSiN["f2"].size() == 0 || silData->fSiN["f3"].size() == 0)};
        if(LateralAndHeavy)
        {
            modularDataOut->fLeaves["GATCONF"] = 4;
        }
        else if(FrontAndHeavy)
        {
            modularDataOut->fLeaves["GATCONF"] = 8;
        }
        else if(BothFront)
        {
            modularDataOut->fLeaves["GATCONF"] = 12;
        }
        else
        {
            // No heavy particle or no silicon hits
        }
        outTreeTPC->Fill();
        outTreeSils->Fill();

        TPCDataOut->Clear();
        silDataOut->Clear();
        modularDataOut->Clear();
    }
    // Write the output files
    outFileTPC->Write();
    outFileSils->Write();
    outFileTPC->Close();
    outFileSils->Close();
}