#ifndef do_profiles_cxx
#define do_profiles_cxx

#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"
#include "TRandom.h"
#include "TLine.h"
#include "TROOT.h"
#include "Math/Point3D.h"
#include "TProfile.h"

#include "ActSRIM.h"
#include "ActKinematics.h"

using XYZPoint = ROOT::Math::XYZPoint;
using XYZVector = ROOT::Math::XYZVector;

bool IsInsideACTAR(XYZPoint point)
{
    bool a {0 <= point.X() && point.X() <= 256};
    bool b {0 <= point.Y() && point.Y() <= 256};
    bool c {0 <= point.Z() && point.Z() <= 256};
    return a && b && c;
}

void FillHistogram (ActPhysics::SRIM* srim, const std::string& particle, double T3Lab, XYZPoint vertex, XYZVector direction, double step, TGraph* graph, std::vector<double>& finalPoints)
{
    double initialRange {srim->EvalRange(particle, T3Lab)};
    double Eiter = T3Lab;
    bool isOutside {false};
    for(double r = step; r < initialRange; r += step)
    {
        double EpostSlow {srim->Slow(particle, Eiter, step)};
        double eLoss {Eiter - EpostSlow};

        XYZPoint stepPoint {vertex + r * direction.Unit()};
        if(!IsInsideACTAR(stepPoint) && !isOutside)
        {
            finalPoints.push_back(r);
            isOutside = true;
            std::cout << "Out of ACTAR at r = " << r << std::endl;
        }

        // Fill the graph
        graph->SetPoint(graph->GetN(), r, eLoss);
        // std::cout << "r: " << r << ", eLoss: " << eLoss << std::endl;

        Eiter = EpostSlow;
    }
    if (!isOutside)
        finalPoints.push_back(-100);
}

void FillHistogramWithUncertanty (ActPhysics::SRIM* srim, const std::string& particle, double T3Lab, XYZPoint vertex, XYZVector direction, double step, TProfile* profile, std::vector<double>& finalPoints, int nIterations, double sigma_r)
{
    double initialRange {srim->EvalRange(particle, T3Lab)};
    bool isOutside {false};
    std::vector<double> finalPointAllIterations {};

    for (int iter = 0; iter < nIterations; iter++) // Repit the process to introduce uncertainty
    {
        double Eiter = T3Lab;
        for(double r = step; r < initialRange; r += step)
        {
            double EpostSlow {srim->Slow(particle, Eiter, step)};
            double eLoss {Eiter - EpostSlow};

            // Introducir incertidumbre en r
            double r_measured = gRandom->Gaus(r, sigma_r);

            XYZPoint stepPoint {vertex + r_measured * direction.Unit()};
            if(!IsInsideACTAR(stepPoint) && !isOutside)
            {
                finalPointAllIterations.push_back(r_measured);
                isOutside = true;
                std::cout << "Out of ACTAR at r = " << r_measured << std::endl;
            }

            // Llenar el perfil con el valor de energía perdida
            profile->Fill(r_measured, eLoss);

            Eiter = EpostSlow;
        }
    }

    if (!isOutside)
    {
        finalPoints.push_back(-100);
    }
    else
    {
        double mean {0};
        for (auto& r : finalPointAllIterations)
        {
            mean += r;
        }
        mean /= finalPointAllIterations.size();
        finalPoints.push_back(mean);
    }
}

void do_profiles() 
{
    // SRIM
    auto* srim {new ActPhysics::SRIM};
    srim->SetUseSpline(false);
    srim->ReadTable("light", "../SRIM files/proton_900mb_CF4_95-5.txt");
    srim->ReadTable("beam", "../SRIM files/11Li_900mb_CF4_95-5.txt");
    srim->ReadTable("heavy", "../SRIM files/12Li_900mb_CF4_95-5.txt");

    // Kinematics
    const std::string& beam {"11Li"};
    const std::string& target {"2H"};
    const std::string& light {"1H"};
    const std::string& heavy {"12Li"};
    const double Tbeam {11 * 7.5}; // MeV
    const double Ex {0.};

    auto* kin {new ActPhysics::Kinematics {beam, target, light, heavy, Tbeam, Ex}};

    // We will do the profiles for a limited number of angles
    //const std::vector<double> thetas {0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170.};
    const std::vector<double> thetas {10., 20., 30, 40};
    std::vector<TGraph*> graphs;

    std::vector<double> finalPoints {}; 

    for(double thetaCM : thetas)
    {
        // Set the angle
        double phi3CM {gRandom -> Uniform(0, 2 * TMath::Pi())};
        kin->ComputeRecoilKinematics(thetaCM * TMath::DegToRad(), phi3CM);
        
        // Get the lab angle and energy
        double theta3Lab {kin->GetTheta3Lab()};
        double T3Lab {kin->GetT3Lab()};

        XYZPoint vertex {128, 128, 128};
        XYZVector direction {TMath::Cos(theta3Lab), TMath::Sin(theta3Lab) * TMath::Sin(phi3CM),
                             TMath::Sin(theta3Lab) * TMath::Cos(phi3CM)};
        
        // Create the graph and fill it
        double step {0.5}; // mm
        auto gProfile {new TGraph};
        gProfile->SetTitle(TString::Format("Eloss profile for #theta_{Lab} = %.1f #circ;distance [mm];dE/dx [MeV/%.2fmm]", theta3Lab * TMath::RadToDeg(), step));

        FillHistogram(srim, "light", T3Lab, vertex, direction, step, gProfile, finalPoints);
        graphs.push_back(gProfile);
    }

    // Draw
    auto* c0 {new TCanvas {"c0", "Eloss profiles"}};
    c0->DivideSquare(graphs.size());
    for(int i = 0; i < graphs.size(); i++)
    {
        c0->cd(i + 1);
        graphs[i]->DrawClone("APL");
        gPad->Update();
        std::cout<<finalPoints[i]<<std::endl;   
        auto* line {new TLine {finalPoints[i], gPad->GetUymin(), finalPoints[i], gPad->GetUymax()}};
        line->SetLineWidth(2);
        line->Draw();
        gROOT->SetSelectedPad(nullptr);

    }

}
#endif