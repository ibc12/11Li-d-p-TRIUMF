#ifndef do_profiles_cxx
#define do_profiles_cxx

#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"
#include "TRandom.h"


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

void do_profiles() 
{
    // SRIM
    auto* srim {new ActPhysics::SRIM};
    srim->ReadTable("light", "../SRIM files/proton_900mb_CF4_90-10.txt");
    srim->ReadTable("beam", "../SRIM files/11Li_900mb_CF4_90-10.txt");
    srim->ReadTable("heavy", "../SRIM files/12Li_900mb_CF4_90-10.txt");
    srim->ReadTable("lightInSil", "../SRIM files/protons_silicon.txt");
    srim->ReadTable("heavyInSil", "../SRIM files/12Li_silicon.txt");

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
    const std::vector<double> thetas {10., 20.};
    std::vector<TGraph*> graphs;

    for(double thetaCM : thetas)
    {
        // Set the angle
        double phi3CM {gRandom -> Uniform(0, 2 * TMath::Pi())};
        kin->ComputeRecoilKinematics(thetaCM * TMath::DegToRad(), phi3CM);
        
        // Get the lab angle and energy
        double theta3Lab {kin->GetTheta3Lab()};
        double T3Lab {kin->GetT3Lab()};

        double initialRange {srim->EvalRange("light", T3Lab)};

        double step {initialRange / 1000.};
        std::cout << "Initial range: " << initialRange << ", step: " << step << std::endl;

        XYZPoint vertex {128, 128, 128};
        XYZVector direction {TMath::Cos(theta3Lab), TMath::Sin(theta3Lab) * TMath::Sin(phi3CM),
                             TMath::Sin(theta3Lab) * TMath::Cos(phi3CM)};
        

        double Eiter = T3Lab;

        auto gProfile {new TGraph};
        gProfile->SetTitle(TString::Format("Eloss profile for #theta_{CM} = %.1f #circ;distance [mm];dE/dx [MeV/%.2fmm]", thetaCM, step));

        for(double r = step; r < initialRange; r += step)
        {
            double EpostSlow {srim->Slow("light", Eiter, step)};
            double eLoss {Eiter - EpostSlow};

            XYZPoint stepPoint {vertex + r * direction.Unit()};
            if(!IsInsideACTAR(stepPoint))
            {

                break;
            }

            // Fill the graph
            gProfile->SetPoint(gProfile->GetN(), r, eLoss);
            std::cout << "r: " << r << ", eLoss: " << eLoss << std::endl;

            Eiter = EpostSlow;
        }

        graphs.push_back(gProfile);
    }

    // Draw
    auto* c0 {new TCanvas {"c0", "Eloss profiles"}};
    c0->DivideSquare(graphs.size());
    for(int i = 0; i < graphs.size(); i++)
    {
        c0->cd(i + 1);
        graphs[i]->DrawClone("AL");
    }

}
#endif