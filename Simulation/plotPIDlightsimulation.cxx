#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RResultPtr.hxx"
#include "Rtypes.h"

#include "TAttLine.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TMath.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TVirtualPad.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TMultiGraph.h"

#include <string>
#include <vector>

void plotPIDlightsimulation()
{
    double T1{7.5};
    std::ostringstream oss;
    oss << "../Outputs/" << std::fixed << std::setprecision(1) << T1 << "MeV/";
    std::string inputPath = oss.str();

    std::string silConfig{"silspecs_spacer"};
    std::vector<std::string> files{
        //inputPath + "2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_" + silConfig + ".root",
        //inputPath + "2H_1H_TRIUMF_Eex_0.130_nPS_0_pPS_0_" + silConfig + ".root",
        //inputPath + "2H_1H_TRIUMF_Eex_0.435_nPS_0_pPS_0_" + silConfig + ".root",
        //inputPath + "2H_3H_TRIUMF_Eex_0.000_nPS_0_pPS_0_" + silConfig + ".root",
        //inputPath + "2H_2H_TRIUMF_Eex_0.000_nPS_0_pPS_0_" + silConfig + ".root",
        inputPath + "2H_3He_TRIUMF_Eex_0.000_nPS_0_pPS_0_" + silConfig + ".root",
        inputPath + "2H_4He_TRIUMF_Eex_0.000_nPS_0_pPS_0_" + silConfig + ".root",
    };
    const ROOT::RDF::TH2DModel PIDLightlength {"hPIDlength", "PID length;E_{Sil} [MeV];#DeltaE_{gas} [MeV]", 350, 0, 90, 350, 0,
                                90};
    auto chainCuts{new TChain("SimulationTTree")};
    auto chainHeavy{new TChain("SimulationTTreeHeavy")};
    for (const auto &file : files)
    {
        chainCuts->Add(file.c_str());
        chainHeavy->Add(file.c_str());
    }
    chainCuts->AddFriend(chainHeavy);
    // chainCuts->Print();

    auto df = ROOT::RDataFrame("SimulationTTree", files);
    auto histoPIDgas{df.Histo2D(PIDLightlength, "DeltaESil", "DeltaEgas")};
    histoPIDgas->GetXaxis()->SetTitle("E_{Sil} [MeV]");
    histoPIDgas->GetYaxis()->SetTitle("#DeltaE_{gas} [keV]");
    auto histoPIDSil{df.Histo2D(PIDLightlength, "DeltaESil1", "DeltaESil")};
    histoPIDSil->GetXaxis()->SetTitle("E_{Sil1} [MeV]");
    histoPIDSil->GetYaxis()->SetTitle("E_{Sil1} [MeV]");

    //auto hT3T4{new TH2D("hT3T4", "Kinetic Energy of light vs heavy;T3_{Lab} [MeV];T4_{Lab} [MeV]", 100, 0, 60, 100, 0, 100)};
    //df.Foreach([&hT3T4](double T3Lab, double T4Lab)
    //           { hT3T4->Fill(T3Lab, T4Lab); }, {"T3Lab", "T4Lab"});
    //
    //auto dfcut = ROOT::RDataFrame(*chainCuts);
    //auto histoCut{dfcut.Histo2D({"hT3T4", "Kinetic correlations;T_{3} [MeV];T_{4} [MeV]", 100, 0, 60, 100, 0, 100}, "EVertex", "T4Lab")};

    auto c1 = new TCanvas("c1", "PID Heavy", 800, 600);
    c1->DivideSquare(4);
    c1->cd(1);
    histoPIDgas->DrawClone("colz");
    c1->cd(2);
    histoPIDSil->DrawClone("colz");
    // c1->cd(3);
    // hT3T4->DrawClone("colz");
    // c1->cd(4);
    // histoCut->DrawClone("colz");
}