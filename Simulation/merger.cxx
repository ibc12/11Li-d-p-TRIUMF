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

#include "ActCrossSection.h"


#include <string>
#include <vector>
void merger()
{
    // ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(10000);


    ROOT::EnableImplicitMT();

    double T1 {7.5};
    std::vector<double> Exs {0., 0.130, 0.435};

 // Construct the output folder path based on T1 with fixed precision
    std::ostringstream oss;
    oss << "../Outputs/" << std::fixed << std::setprecision(1) << T1 << "MeV/";
    std::string outputPath = oss.str();

    std::vector<std::string> files {
        outputPath + "2H_1H_TRIUMF_Eex_0.000_nPS_0_pPS_0_silicons_reverse.root",
        outputPath + "2H_1H_TRIUMF_Eex_0.130_nPS_0_pPS_0_silicons_reverse.root",
        outputPath + "2H_1H_TRIUMF_Eex_0.435_nPS_0_pPS_0_silicons_reverse.root"
    };

    // Read dfs
    std::vector<ROOT::RDF::RNode> dfs;
    for(const auto& file : files)
    {
        dfs.push_back(ROOT::RDataFrame {"SimulationTTree", file});
    }
    // Compute scaling factors
    double gasDensity {2.428e-4}; // g/cm3
    double gasMolarDensity {0.9 * 4.0282 + 0.1 * 58.12}; // g/mol 
    //double Nt {(gasDensity/gasMolarDensity) * 6.022e23 * 25.6 * (95 * 2 /(95*2 +5*4+5*10))}; // particles/cm3 * ACTAR length
    double Nt {4.6688e19 * 25.6 * 0.8877};
    double Np {3000 * 6 * 24 * 3600}; // 3e5 pps 6 days
    double Nit {1.e6};

    // Set histogram
    int nbins {200};
    double xmin1 {-1};
    double xmax1 {2};
    auto* hEx {new TH1D {
        "hEx", TString::Format("Ex para todos os picos;E_{x} [MeV];Contas / %.0f keV", (xmax1 - xmin1) / nbins * 1000), nbins,
        xmin1, xmax1}};
    hEx->Sumw2();
    double ymin {0};
    double ymax {40};
    double xmin2 {0};
    double xmax2 {180};
    auto* hKin {new TH2D {
        "hKin", "Kinematics for all E_{x};Theta_{Lab} [degree];E_{Lab} [MeV]", nbins,
        xmin2, xmax2, nbins, ymin, ymax}};
    hKin->Sumw2();
    std::vector<TH1D*> hs1;
    std::vector<TH2D*> hs2;

    int contador {0};
    for(auto& df : dfs)
    {
        double Ex = Exs[contador];
        auto* xs {new ActSim::CrossSection()};
        if(Ex == 0)
        {
            TString data_to_read {TString::Format("../Inputs/TheoXS/%.1fMeV/dp/angs12nospin.dat", T1)};
            xs->ReadFile(data_to_read.Data());
            std::cout << xs->GetTotalXSmbarn() << std::endl;
        }
        else if(Ex == 0.130)
        {
            TString data_to_read {TString::Format("../Inputs/TheoXS/%.1fMeV/dp/angp12nospin.dat", T1)};
            xs->ReadFile(data_to_read.Data());
            std::cout << xs->GetTotalXSmbarn() << std::endl;
        }
        else if(Ex == 0.435)
        {
            TString data_to_read {TString::Format("../Inputs/TheoXS/%.1fMeV/dp/angp32nospin.dat", T1)};
            xs->ReadFile(data_to_read.Data());
            std::cout << xs->GetTotalXSmbarn() << std::endl;
        }

        double totalXS {xs->GetTotalXScm2()};
        double scaling {(Nt * Np * totalXS) / Nit};

        // FILTRO por theta3CM en el rango deseado
        auto dfFiltered = df.Filter("theta3CM >= 0 && theta3CM <= 180", "Angular cut on theta3CM");

        auto h1 {dfFiltered.Histo1D(
            {"h1", "Ex in for loop", hEx->GetNbinsX(), hEx->GetXaxis()->GetXmin(), hEx->GetXaxis()->GetXmax()}, "Eex")};
        h1->Scale(scaling);
        hEx->Add(h1.GetPtr());
        hs1.push_back((TH1D*)h1->Clone());

        auto h2 {dfFiltered.Histo2D(
            {"h2", "Kinematics in for loop", hKin->GetNbinsX(), hKin->GetXaxis()->GetXmin(), hKin->GetXaxis()->GetXmax(), 
                hKin->GetNbinsY(), hKin->GetYaxis()->GetXmin(), hKin->GetYaxis()->GetXmax()}, "theta3Lab", "EVertex")};
        h2->Scale(scaling);
        hKin->Add(h2.GetPtr());
        hs2.push_back((TH2D*)h2->Clone());

        contador += 1;
        delete xs;
    }

    auto* f {new TF1{"f", "[0] * TMath::Voigt(x - [1], [2], [3]) + [4] * TMath::Voigt(x - [5], [6], [7])  + [8] * TMath::Voigt(x - [9], [10], [11]) ", -2, 2}};
    Double_t params[12] = {150, 0, 0.1018, 0.1, 250, 0.13, 0.08895, 0.02, 140, 0.4, 0.09646, 0.08};
    f->SetParameters(params);
    // f->FixParameter(2, 0.102165);
    // f->FixParameter(6, 0.0892859);
    // f->FixParameter(10, 0.0959508);
    f->FixParameter(1,0);
    //f->FixParameter(7, 0.015);
    f->SetParLimits(0, 0, 3500);
    f->SetParLimits(1, -0.01, 0.01);
    f->SetParLimits(3, 0.05, 0.5);
    f->SetParLimits(4, 10, 3500);
    f->SetParLimits(5, 0.121, 0.140);
    f->SetParLimits(7, 0.01, 0.025);
    f->SetParLimits(8, 10, 3500);
    f->SetParLimits(9, 0.4, 0.5);
    f->SetParLimits(11, 0.03, 0.2);

    // Set the minimizer to Minuit2
    //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Simplex");

    hEx->Fit(f, "0M+I10000");
    
    auto* fGS {new TF1{"fGS", "[0] * TMath::Voigt(x - [1], [2], [3])", -2, 2}};
    double* paramsGS = f->GetParameters();
    fGS->SetParameters(paramsGS);
    auto* f1st {new TF1{"fGS", "[0] * TMath::Voigt(x - [1], [2], [3])", -2, 2}};
    double params1st[4];
    params1st[0] = f->GetParameter(4);
    params1st[1] = f->GetParameter(5);
    params1st[2] = f->GetParameter(6);
    params1st[3] = f->GetParameter(7);
    f1st->SetParameters(params1st);
    auto* f2nd {new TF1{"fGS", "[0] * TMath::Voigt(x - [1], [2], [3])", -2, 2}};
    double params2nd[4];
    params2nd[0] = f->GetParameter(8);
    params2nd[1] = f->GetParameter(9);
    params2nd[2] = f->GetParameter(10);
    params2nd[3] = f->GetParameter(11);
    f2nd->SetParameters(params2nd);

    // plot
    std::vector<int> colors {6, 8, 46};
    auto* c0 {new TCanvas {"c0", "Merger canvas 1D"}};
    gStyle->SetOptStat(0);
    hEx->SetLineWidth(2);
    hEx->Draw("hist");
    hEx->GetXaxis()->SetRangeUser(-0.8, 1);
    f->Draw("same");
    fGS->SetLineColor(colors[0]);
    fGS->Draw("same");
    f1st->SetLineColor(colors[1]);
    f1st->Draw("same");
    f2nd->SetLineColor(colors[2]);
    f2nd->Draw("same");
    // hEx->SaveAs("hit_merger_Exs.root");
    std::vector<std::string> labels {"0", "0.130", "0.435"};
    auto* leg1 {new TLegend {0.2, 0.2}};
    leg1->SetHeader("E_{x} [MeV]");
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    for(int i = 0; i < hs1.size(); i++)
    {
        hs1[i]->SetLineColor(colors[i]);
        hs1[i]->SetLineStyle(kDashed);
        hs1[i]->SetLineWidth(2);
        leg1->AddEntry(hs1[i], labels[i].c_str());
        hs1[i]->Draw("hist same");
    }
    leg1->Draw();
}