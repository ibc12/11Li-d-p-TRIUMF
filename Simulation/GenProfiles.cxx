#include "ActKinematics.h"
#include "ActParticle.h"
#include "ActSRIM.h"

#include "Rtypes.h"

#include "TBranch.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"
#include "TMath.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TSpline.h"
#include "TString.h"
#include "TTree.h"
#include "TVirtualPad.h"

#include <algorithm>
#include <cstdlib>
#include <ios>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "./GetFuncs.cxx"


TGraph* FuncDerivative(TH1D* h, double& xmin, double& xmax)
{
    // Create function
    auto func {std::make_unique<TF1>(
        "f", [=](double* x, double* p) { return h->Interpolate(x[0]); }, h->GetXaxis()->GetXmin(),
        h->GetXaxis()->GetXmax(), 1)};
    func->SetLineColor(kBlue);
    func->SetLineWidth(2);
    func->SetTitle("Derivative;X [mm];Hist derivative");
    func->SetNpx((int)h->GetXaxis()->GetXmax() * 4);
    // Create TGraph
    auto* gd {new TGraph {func.get(), "d"}};
    // Get maximum and minimum and their positions along X
    auto y {gd->GetY()};
    auto x {gd->GetX()};
    auto n {gd->GetN()};
    auto imax {std::max_element(y, y + n) - y};
    auto imin {std::min_element(y, y + n) - y};
    xmax = x[imax];
    xmin = x[imin];
    return gd;
}

void HistDerivative(TH1D* h, double& xmin, double& xmax, TGraph*& g, TSpline3*& spe)
{
    // Manually build derivative
    auto* hd {new TH1D {"hd", "Derivative;X [mm];Derivative", h->GetNbinsX(), h->GetXaxis()->GetXmin(),
                        h->GetXaxis()->GetXmax()}};
    for(auto bin = 1; bin <= h->GetNbinsX(); bin++)
    {
        double val {};
        if(bin > 1)
        {
            auto y1 {h->GetBinContent(bin)};
            auto y0 {h->GetBinContent(bin - 1)};
            val = (y1 - y0);
        }
        hd->SetBinContent(bin, val);
    }
    // Build spline
    spe = new TSpline3 {hd, "b2,e2", 0, 0};
    spe->SetLineColor(kBlue);
    spe->SetLineStyle(2);
    spe->SetLineWidth(2);
    // Build graph
    g = new TGraph {hd};
    g->SetLineColor(kTeal);
    g->SetLineWidth(2);
    // Get min and max
    auto y {g->GetY()};
    auto x {g->GetX()};
    auto n {g->GetN()};
    auto imax {std::max_element(y, y + n) - y};
    auto imin {std::min_element(y, y + n) - y};
    xmax = x[imax];
    xmin = x[imin];
    // Delete temporary histo
    delete hd;
}

std::pair<double, double> GetMinMaxFrame(TH1D* h, TGraph* g)
{
    // Max is given by histogram
    auto binMax {h->GetMaximumBin()};
    auto ymax {h->GetBinContent(binMax)};
    // Minimum by min of graph
    auto ymin {TMath::MinElement(g->GetN(), g->GetY())};
    return {ymin * 1.1, ymax * 1.1};
}

void GenProfiles()
{
    // Set reaction
    std::string beam {"8Li"};
    std::string target {"a"};
    std::string light {"n"};

    bool draw {true};

    // Set energy of the beam
    double Ebeam {4.8}; // MeV

    // Set parameters of ACTAR
    double lActarMM {256}; // mm
    int npads {128};

    // Step in X to compute eloss
    double xstep {2}; // mm

    // Number of iterations
    int iterations {1000};

    // Init kinematics
    ActPhysics::Particle pb {beam};
    // Ebeam *= pb.GetAMU();
    ActPhysics::Kinematics kin {beam, target, light, Ebeam};

    // Init SRIM
    ActPhysics::SRIM srim;
    srim.ReadTable("beam", "../SRIM/8Li_He_CO2.txt");
    if(light != "n")
        srim.ReadTable("light", TString::Format("../SRIM/%s_He_CO2.txt", light.c_str()).Data());
    srim.ReadTable("heavy", TString::Format("../SRIM/%s_He_CO2.txt", kin.GetParticle(4).GetName().c_str()).Data());

    // Histograms
    auto* hX {new TProfile {"hX", "Eloss per X;X [mm];ELoss [MeV/mm]", npads, 0, lActarMM}};
    hX->SetLineWidth(2);
    auto* hDiff1 {new TH1D {"hDiff1", "Difference 1;Diff [mm]", 100, -5, 5}};
    auto* hDiff2 {(TH1D*)hDiff1->Clone("hDiff2")};
    hDiff2->SetTitle("Difference 2");
    // Draw histograms to check kinematics
    auto* hKinLight {
        new TH2D {"hKinLight", "Kinematics of light;#theta_{Lab} [#circ];E_{Light} [MeV]", 200, 0, 180, 200, 0, 15}};
    auto* hKinHeavy {
        new TH2D {"hKinHeavy", "Kinematics of heavy;#theta_{Lab} [#circ];E_{Heavy} [MeV]", 200, 0, 180, 200, 0, 15}};

    // Canvas
    auto* c0 {new TCanvas {"c0", "Eloss canvas"}};

    // Save
    TFile* file {};
    TTree* tree {};
    auto* data {new std::vector<double>};
    if(!draw)
    {
        file = new TFile {TString::Format("./Outputs/with_%s.root", light.c_str()), "recreate"};
        tree = new TTree {"tree", "Good reactions tree"};
        tree->Branch("data", &data);
    }

    // Run and print
    const int percentPrint {5};
    int step {iterations / (100 / percentPrint)};
    int nextPrint {step};
    int percent {};
    for(auto it = 0; it < iterations; it++)
    {
        // Print progress
        if(it >= nextPrint)
        {
            percent = 100 * it / iterations;
            std::cout << "\r" << std::string(percent / percentPrint, '|') << percent << "%";
            std::cout.flush();
            nextPrint += step;
        }

        // 1 -> Sample vertex
        auto vertexX {gRandom->Uniform() * lActarMM};
        // Assume point at the centre of actar
        XYZPoint vertex {vertexX, lActarMM / 2, lActarMM / 2};

        // 2 -> Eloss until vertex (aka beam)
        FillHist(&srim, "beam", Ebeam, {0, lActarMM / 2, lActarMM / 2}, {1, 0, 0}, vertexX, hX);

        // 3-> Reaction
        auto Ereac {srim.Slow("beam", Ebeam, vertexX)};
        kin.SetBeamEnergy(Ereac);
        // Flat distribution
        // auto thetaCM {TMath::ACos(gRandom->Uniform(-1, 1))};
        auto thetaCM {gRandom->Uniform(TMath::PiOver2(), TMath::Pi())};
        auto phi {gRandom->Uniform(0, TMath::TwoPi())};
        kin.ComputeRecoilKinematics(thetaCM, phi);

        // 4-> Light particle
        auto T3 {kin.GetT3Lab()};
        auto theta3 {kin.GetTheta3Lab()};
        auto phi3 {kin.GetPhi3Lab()};
        XYZVector dirLight {TMath::Cos(theta3), TMath::Sin(theta3) * TMath::Sin(phi3),
                            TMath::Sin(theta3) * TMath::Cos(phi3)};
        hKinLight->Fill(theta3 * TMath::RadToDeg(), T3);

        if(light != "n")
        {
            auto rangeLight {srim.EvalDirect("light", T3)};
            FillHist(&srim, "light", T3, vertex, dirLight, rangeLight, hX);
        }

        // 5-> Heavy particle
        auto T4 {kin.GetT4Lab()};
        auto theta4 {kin.GetTheta4Lab()};
        auto phi4 {kin.GetPhi4Lab()};
        // Build direction
        XYZVector dirHeavy {TMath::Cos(theta4), TMath::Sin(theta4) * TMath::Sin(phi4),
                            TMath::Sin(theta4) * TMath::Cos(phi4)};
        // Range of heavy recoil
        auto rangeHeavy {srim.EvalDirect("heavy", kin.GetT4Lab())};
        FillHist(&srim, "heavy", T4, vertex, dirHeavy, rangeHeavy, hX);
        hKinHeavy->Fill(theta4 * TMath::RadToDeg(), T4);

        // Draw
        if(draw)
        {
            // Compute RP
            TGraph* gd2 {};
            TSpline3* spe {};
            double xmin2 {};
            double xmax2 {};
            HistDerivative(hX, xmin2, xmax2, gd2, spe);
            // Compute difference
            auto diff2 {xmax2 - vertexX};
            hDiff2->Fill(diff2);


            c0->DivideSquare(4);
            c0->cd(1);
            hX->SetTitle(TString::Format(
                "Evt : %d at X = %.1f mm and #theta_{Lab} = %.1f #circ T_{Lab} = %.1f MeV T_{Beam} = %.1f MeV", it,
                vertexX, theta4 * TMath::RadToDeg(), T4, Ereac));
            hX->Draw("hist");
            gPad->Update();

            auto* lv {new TLine {vertexX, gPad->GetUymin(), vertexX, gPad->GetUymax()}};
            lv->SetLineWidth(2);
            lv->SetLineColor(kMagenta);
            lv->Draw();

            c0->cd(2);
            double pmin, pmax;
            std::tie(pmin, pmax) = GetMinMaxFrame(hX, gd2);
            auto* f2 {gPad->DrawFrame(0, pmin, 256 + 10, pmax)};
            f2->GetXaxis()->SetTitle("X [mm]");
            hX->Draw("hist same");
            gd2->Draw("l");
            spe->Draw("l same");
            // New vertex
            auto* lmax2 {new TLine {xmax2, pmin, xmax2, pmax}};
            lmax2->SetLineWidth(2);
            lmax2->SetLineColor(kOrange);
            lmax2->Draw("same");
            lv->DrawLine(vertexX, pmin, vertexX, pmax);

            c0->Update();
            c0->WaitPrimitive("what", "");
            hX->GetXaxis()->UnZoom();
            hX->GetYaxis()->UnZoom();
            c0->Clear();
            c0->Update();
            gROOT->SetSelectedPad(nullptr);
            // Delete lines
            for(auto* l : {lv, lmax2})
                delete l;
            delete gd2;
            delete spe;
        }

        // Save
        if(file)
        {
            data->clear();
            // WARNING: We must use projectionX of the histogram, otherwise the returned array does not contain the
            // mean! Idk what it contains
            auto* proj {hX->ProjectionX()};
            auto* p {proj->GetArray()};
            *data = std::vector<double>(p, p + hX->GetNbinsX());
            tree->Fill();
            delete proj;
        }

        // Mandatory
        hX->Reset();
    }
    std::cout << '\n';

    c0->Close();
    auto* c1 {new TCanvas {"c1", "Final canvas"}};
    c1->DivideSquare(6);
    c1->cd(1);
    hDiff1->Draw();
    c1->cd(2);
    hDiff2->Draw();
    c1->cd(3);
    hKinLight->Draw("colz");
    c1->cd(4);
    hKinHeavy->Draw("colz");

    if(file)
    {
        tree->Write();
        file->Close();
    }
}
