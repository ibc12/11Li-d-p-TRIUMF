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
#include "TSpline.h"
#include "TLegend.h"
#include "TLatex.h"

#include "ActSRIM.h"
#include "ActKinematics.h"

using XYZPoint = ROOT::Math::XYZPoint;
using XYZVector = ROOT::Math::XYZVector;

bool IsInsideACTAR(XYZPoint point)
{
    bool a{0 <= point.X() && point.X() <= 256};
    bool b{0 <= point.Y() && point.Y() <= 256};
    bool c{0 <= point.Z() && point.Z() <= 256};
    return a && b && c;
}

void FillHistogram(ActPhysics::SRIM *srim, const std::string &particle, double initialRange, XYZPoint vertex, XYZVector direction, double step, TGraph *graph, std::vector<double> &finalPoints)
{
    double T3Lab{srim->EvalEnergy(particle, initialRange)};
    double Eiter = T3Lab;
    bool isOutside{false};
    for (double r = step; r < initialRange; r += step)
    {
        double EpostSlow{srim->Slow(particle, Eiter, step)};
        double eLoss{Eiter - EpostSlow};

        XYZPoint stepPoint{vertex + r * direction.Unit()};
        if (!IsInsideACTAR(stepPoint) && !isOutside)
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

void FillHistogramWithUncertantyInRange(ActPhysics::SRIM *srim, const std::string &particle, double initialRange, XYZPoint vertex, XYZVector direction, double step, TProfile *profile, std::vector<double> &finalPoints, int nIterations, double sigma_r)
{
    double T3Lab{srim->EvalEnergy(particle, initialRange)};
    bool isOutside{false};
    std::vector<double> finalPointAllIterations{};

    for (int iter = 0; iter < nIterations; iter++) // Repit the process to introduce uncertainty
    {
        double Eiter = T3Lab;
        for (double r = step; r < initialRange; r += step)
        {
            double EpostSlow{srim->Slow(particle, Eiter, step)};
            double eLoss{Eiter - EpostSlow};

            // Introducir incertidumbre en r
            double r_measured = gRandom->Gaus(r, sigma_r);

            XYZPoint stepPoint{vertex + r_measured * direction.Unit()};
            if (!IsInsideACTAR(stepPoint) && !isOutside)
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
        double mean{0};
        for (auto &r : finalPointAllIterations)
        {
            mean += r;
        }
        mean /= finalPointAllIterations.size();
        finalPoints.push_back(mean);
    }
}

void FillHistogramWithUncertantyInRangeAndELoss(ActPhysics::SRIM *srim, const std::string &particle, double initialRange, XYZPoint vertex, XYZVector direction, double step, TProfile *profile, std::vector<double> &finalPoints, int nIterations, double sigma_r, double sigma_eLoss_percent)
{
    double T3Lab{srim->EvalEnergy(particle, initialRange)};
    bool isOutside{false};
    std::vector<double> finalPointAllIterations{};

    for (int iter = 0; iter < nIterations; iter++) // Repit the process to introduce uncertainty
    {
        double Eiter = T3Lab;
        for (double r = step; r < initialRange; r += step)
        {
            double EpostSlow{srim->Slow(particle, Eiter, step)};
            double eLoss{Eiter - EpostSlow};

            // Introduction of uncertainty in r and eLoss
            double r_measured = gRandom->Gaus(r, sigma_r);
            double eLoss_measured = gRandom->Gaus(eLoss, eLoss * sigma_eLoss_percent);

            XYZPoint stepPoint{vertex + r_measured * direction.Unit()};
            if (!IsInsideACTAR(stepPoint) && !isOutside)
            {
                finalPointAllIterations.push_back(r_measured);
                isOutside = true;
                std::cout << "Out of ACTAR at r = " << r_measured << std::endl;
            }

            // Llenar el perfil con el valor de energía perdida
            profile->Fill(r_measured, eLoss_measured);

            Eiter = EpostSlow;
        }
    }

    if (!isOutside)
    {
        finalPoints.push_back(-100);
    }
    else
    {
        double mean{0};
        for (auto &r : finalPointAllIterations)
        {
            mean += r;
        }
        mean /= finalPointAllIterations.size();
        finalPoints.push_back(mean);
    }
}

double plotSplinesAndRangeValue(TProfile *profile, bool draw = true)
{
    if (!profile)
    {
        std::cerr << "Error: TProfile no válido." << std::endl;
        return -1;
    }

    int nPoints = profile->GetNbinsX();
    std::vector<double> x(nPoints), y(nPoints);
    for (int i = 1; i <= nPoints; ++i)
    {
        x[i - 1] = profile->GetBinCenter(i);
        y[i - 1] = profile->GetBinContent(i);
    }

    TGraph graph(nPoints, x.data(), y.data());
    TSpline3 *spline = new TSpline3("spline", &graph);

    // Get max value
    double maxVal = profile->GetMaximum();
    // std::cout << "Máximo: " << maxVal << std::endl;
    // Get max position
    double maxPos = profile->GetBinCenter(profile->GetMaximumBin());

    // Find point where the spline is half the maximum only for x > maxPos
    double halfMax = maxVal / 2;
    double halfMaxPos = -1;
    for (int i = 1; i <= nPoints; ++i)
    {
        if (x[i - 1] > maxPos)
        {
            double ySpline = spline->Eval(x[i - 1]);
            if (ySpline < halfMax)
            {
                halfMaxPos = x[i - 1];
                break;
            }
        }
    }

    if (draw)
    {
        profile->SetLineColor(kBlue);
        profile->Draw();

        spline->SetLineColor(kRed);
        spline->Draw("same");

        if (halfMaxPos != -1)
        {
            TLine *line = new TLine(halfMaxPos, profile->GetMinimum(), halfMaxPos, profile->GetMaximum());
            line->SetLineColor(kGreen);
            line->SetLineStyle(2);
            line->SetLineWidth(2);
            line->Draw("same");
        }
    }

    return halfMaxPos;
}

std::vector<std::string> GetBestMatchingParticles(const std::vector<TGraph *> &graphsProton,
                                                  const std::vector<TGraph *> &graphsDeuteron,
                                                  const std::vector<TGraph *> &graphsTriton,
                                                  const std::vector<TProfile *> &profiles)
{
    std::vector<std::string> bestMatches;

    for (size_t i = 0; i < profiles.size(); ++i)
    {
        auto *prof = profiles[i];

        // Función para comparar una TProfile con un TGraph
        auto computeDifference = [prof](TGraph *gr) -> double
        {
            double totalDiff = 0.0;
            int nMatched = 0;

            for (int j = 1; j <= prof->GetNbinsX(); ++j)
            {
                double xProf = prof->GetBinCenter(j);
                double yProf = prof->GetBinContent(j);

                if (yProf == 0)
                    continue;

                double *xGr = gr->GetX();
                double *yGr = gr->GetY();
                int n = gr->GetN();

                // Buscar el punto más cercano en X
                double minDiffX = 1e9;
                double yClosest = 0;
                for (int k = 0; k < n; ++k)
                {
                    double dx = std::abs(xGr[k] - xProf);
                    if (dx < minDiffX)
                    {
                        minDiffX = dx;
                        yClosest = yGr[k];
                    }
                }

                totalDiff += std::abs(yClosest - yProf);
                ++nMatched;
            }

            return (nMatched > 0) ? totalDiff / nMatched : 1e9;
        };

        double dP = computeDifference(graphsProton[i]);
        double dD = computeDifference(graphsDeuteron[i]);
        double dT = computeDifference(graphsTriton[i]);

        double minDiff = std::min({dP, dD, dT});

        if (minDiff == dP)
            bestMatches.push_back("Se corresponde con p");
        else if (minDiff == dD)
            bestMatches.push_back("Se corresponde con d");
        else
            bestMatches.push_back("Se corresponde con t");
    }

    return bestMatches;
}

void do_profiles()
{
    // Kinematics
    const std::string &beam{"11Li"};
    const std::string &target{"2H"};
    const std::string &light{"1H"};
    const std::string &heavy{"12Li"};
    const double Tbeam{11 * 7.5}; // MeV
    const double Ex{0.};

    // SRIM
    auto *srim{new ActPhysics::SRIM};
    srim->SetUseSpline(true);
    srim->ReadTable("light", "../SRIM files/" + light + "_900mb_CF4_95-5.txt");
    srim->ReadTable("beam", "../SRIM files/" + beam + "_900mb_CF4_95-5.txt");
    srim->ReadTable("heavy", "../SRIM files/" + heavy + "_900mb_CF4_95-5.txt");
    srim->ReadTable("p", "../SRIM files/1H_900mb_CF4_95-5.txt");
    srim->ReadTable("d", "../SRIM files/2H_900mb_CF4_95-5.txt");
    srim->ReadTable("t", "../SRIM files/3H_900mb_CF4_95-5.txt");
    std::string lightLISEgas{"../LISE files/" + light + "_gas_95-5.txt"};
    std::string pLISEgas{"../LISE files/1H_gas_95-5.txt"};
    std::string dLISEgas{"../LISE files/2H_gas_95-5.txt"};
    std::string tLISEgas{"../LISE files/3H_gas_95-5.txt"};
    srim->SetStragglingLISE("light", lightLISEgas);
    srim->SetStragglingLISE("p", pLISEgas);
    srim->SetStragglingLISE("d", dLISEgas);
    srim->SetStragglingLISE("t", tLISEgas);

    auto *kin{new ActPhysics::Kinematics{beam, target, light, heavy, Tbeam, Ex}};

    // We will do the profiles for a limited number of angles
    // const std::vector<double> thetas {0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170.};
    const std::vector<double> thetas{10., 20., 30, 40};
    std::vector<TGraph *> graphs;
    std::vector<TGraph *> graphsProton;
    std::vector<TGraph *> graphsDeuteron;
    std::vector<TGraph *> graphsTriton;
    std::vector<TProfile *> profilesRangeUncertanty;
    std::vector<TProfile *> profilesRangeAndELossUncertanty;

    std::vector<double> finalPointsGraph{};
    std::vector<double> finalPointsProfileRangeUncertanty{};
    std::vector<double> finalPointsProfileRangeElossUncertanty{};

    std::vector<double> interpolatedRangeValues{};

    for (int i = 1; double thetaCM : thetas)
    {
        // Set the angle
        double phi3CM{gRandom->Uniform(0, 2 * TMath::Pi())};
        kin->ComputeRecoilKinematics(thetaCM * TMath::DegToRad(), phi3CM);

        // Get the lab angle and energy
        double theta3Lab{kin->GetTheta3Lab()};
        double T3Lab{kin->GetT3Lab()};
        double range{srim->EvalRange("light", T3Lab)};

        std::cout << "Theta3Lab: " << theta3Lab << ", T3Lab: " << T3Lab << ", Range: " << srim->EvalRange("light", T3Lab) << std::endl;

        XYZPoint vertex{128, 128, 128};
        XYZVector direction{TMath::Cos(theta3Lab), TMath::Sin(theta3Lab) * TMath::Sin(phi3CM),
                            TMath::Sin(theta3Lab) * TMath::Cos(phi3CM)};

        // Create the graph and fill it
        double step{0.5}; // mm
        auto gProfile{new TGraph};
        gProfile->SetTitle(TString::Format("Eloss profile for #theta_{Lab} = %.1f #circ;distance [mm];dE/dx [MeV/%.2fmm]", theta3Lab * TMath::RadToDeg(), step));

        FillHistogram(srim, "light", range, vertex, direction, step, gProfile, finalPointsGraph);
        graphs.push_back(gProfile);

        // Ceate the profile for range uncertanty and fill it
        auto pProfile{new TProfile(TString::Format("profile_uRange_%d", i), "Eloss profile", 1000, 0, srim->EvalRange("light", T3Lab) + 50)};
        pProfile->SetTitle(TString::Format("Eloss profile for #theta_{Lab} = %.1f #circ;distance [mm];dE/dx [MeV/%.2fmm]", theta3Lab * TMath::RadToDeg(), step));
        FillHistogramWithUncertantyInRange(srim, "light", range, vertex, direction, step, pProfile, finalPointsProfileRangeUncertanty, 100, 2);
        profilesRangeUncertanty.push_back(pProfile);

        // Create the profile for range and energy loss uncertanty and fill it
        auto pProfile2{new TProfile(TString::Format("profile_uRange_uEloss_%d", i), "Eloss profile", 1000, 0, srim->EvalRange("light", T3Lab) + 50)};
        pProfile2->SetTitle(TString::Format("Eloss profile for #theta_{Lab} = %.1f #circ;distance [mm];dE/dx [MeV/%.2fmm]", theta3Lab * TMath::RadToDeg(), step));
        FillHistogramWithUncertantyInRangeAndELoss(srim, "light", range, vertex, direction, step, pProfile2, finalPointsProfileRangeElossUncertanty, 100, 2, 0.1);
        profilesRangeAndELossUncertanty.push_back(pProfile2);

        // Compute the range values
        interpolatedRangeValues.push_back(plotSplinesAndRangeValue(pProfile2, false));

        // Compute the profiles for protons, deuterons and tritons
        auto gProfileProton{new TGraph};
        gProfileProton->SetTitle(TString::Format("Eloss profile for #theta_{Lab} = %.1f #circ;distance [mm];dE/dx [MeV/%.2fmm]", theta3Lab * TMath::RadToDeg(), step));
        FillHistogram(srim, "p", range, vertex, direction, step, gProfileProton, finalPointsGraph);
        graphsProton.push_back(gProfileProton);
        auto gProfileDeuteron{new TGraph};
        gProfileDeuteron->SetTitle(TString::Format("Eloss profile for #theta_{Lab} = %.1f #circ;distance [mm];dE/dx [MeV/%.2fmm]", theta3Lab * TMath::RadToDeg(), step));
        FillHistogram(srim, "d", range, vertex, direction, step, gProfileDeuteron, finalPointsGraph);
        graphsDeuteron.push_back(gProfileDeuteron);
        auto gProfileTriton{new TGraph};
        gProfileTriton->SetTitle(TString::Format("Eloss profile for #theta_{Lab} = %.1f #circ;distance [mm];dE/dx [MeV/%.2fmm]", theta3Lab * TMath::RadToDeg(), step));
        FillHistogram(srim, "t", range, vertex, direction, step, gProfileTriton, finalPointsGraph);
        graphsTriton.push_back(gProfileTriton);

        i++;
    }

    // Draw the graphs
    auto *c0{new TCanvas{"c0", "Eloss profiles"}};
    c0->DivideSquare(graphs.size());
    for (int i = 0; i < graphs.size(); i++)
    {
        c0->cd(i + 1);
        graphs[i]->DrawClone("APL");
        gPad->Update();
        // std::cout<<finalPointsGraph[i]<<std::endl;
        auto *line{new TLine{finalPointsGraph[i], gPad->GetUymin(), finalPointsGraph[i], gPad->GetUymax()}};
        line->SetLineWidth(2);
        line->Draw();
        gROOT->SetSelectedPad(nullptr);
    }
    // Draw the profiles for range uncertanty
    auto *c1{new TCanvas{"c1", "Eloss profiles with uncertainty"}};
    c1->DivideSquare(profilesRangeUncertanty.size());
    for (int i = 0; i < profilesRangeUncertanty.size(); i++)
    {
        c1->cd(i + 1);
        profilesRangeUncertanty[i]->DrawClone();
        gPad->Update();
        auto *line{new TLine{finalPointsProfileRangeUncertanty[i], gPad->GetUymin(), finalPointsProfileRangeUncertanty[i], gPad->GetUymax()}};
        line->SetLineWidth(2);
        line->Draw();
        gROOT->SetSelectedPad(nullptr);
    }
    // Draw the profiles for range and energy loss uncertanty
    auto *c2{new TCanvas{"c2", "Eloss profiles with uncertainty"}};
    c2->DivideSquare(profilesRangeAndELossUncertanty.size());
    for (int i = 0; i < profilesRangeAndELossUncertanty.size(); i++)
    {
        c2->cd(i + 1);
        profilesRangeAndELossUncertanty[i]->DrawClone();
        gPad->Update();
        auto *line{new TLine{finalPointsProfileRangeElossUncertanty[i], gPad->GetUymin(), finalPointsProfileRangeElossUncertanty[i], gPad->GetUymax()}};
        line->SetLineWidth(2);
        line->Draw();

        double halfMaxPosition = plotSplinesAndRangeValue(profilesRangeAndELossUncertanty[i], false);
        std::cout << "Half max position: " << halfMaxPosition << ", T3Fit: " << srim->EvalEnergy("light", halfMaxPosition) << std::endl;

        gROOT->SetSelectedPad(nullptr);
    }
    // Draw Splines
    auto c3{new TCanvas{"c3", "Splines"}};
    c3->DivideSquare(profilesRangeAndELossUncertanty.size());
    for (int i = 0; i < profilesRangeAndELossUncertanty.size(); i++)
    {
        c3->cd(i + 1);
        double halfMaxPosition = plotSplinesAndRangeValue(profilesRangeAndELossUncertanty[i]);
    }
    // Compare the profiles for diferent particles
    auto *c4{new TCanvas{"c4", "Eloss profiles for different particles"}};
    c4->DivideSquare(graphsProton.size());
    for (int i = 0; i < graphsProton.size(); i++)
    {
        c4->cd(i + 1);
        graphsProton[i]->SetLineColor(kRed);
        graphsProton[i]->DrawClone("same APL");
        graphsDeuteron[i]->SetLineColor(kBlue);
        graphsDeuteron[i]->DrawClone("same");
        graphsTriton[i]->SetLineColor(kGreen);
        graphsTriton[i]->DrawClone("same");

        // Draw profile with uncertainty
        if (i < profilesRangeAndELossUncertanty.size())
        {
            profilesRangeAndELossUncertanty[i]->SetLineColor(kBlack);
            profilesRangeAndELossUncertanty[i]->SetLineStyle(2); // Dashed line
            profilesRangeAndELossUncertanty[i]->DrawClone("same");
        }

        // Add legend
        auto *leg = new TLegend(0.12, 0.72, 0.40, 0.88);
        leg->AddEntry(graphsProton[i], "Proton", "l");
        leg->AddEntry(graphsDeuteron[i], "Deuteron", "l");
        leg->AddEntry(graphsTriton[i], "Triton", "l");
        if (i < profilesRangeAndELossUncertanty.size())
            leg->AddEntry(profilesRangeAndELossUncertanty[i], "Profile w/ Uncertainty", "l");

        leg->SetBorderSize(0);
        leg->SetTextSize(0.03);
        leg->Draw();
    }
    // Now, check what particle has the best fit
    auto matches = GetBestMatchingParticles(graphsProton, graphsDeuteron, graphsTriton, profilesRangeAndELossUncertanty);

    for (size_t i = 0; i < matches.size(); ++i)
    {
        c4->cd(i + 1);
        TLatex *label = new TLatex(0.15, 0.65, Form("Se corresponde con %s", matches[i].c_str()));
        label->SetNDC(true);
        label->SetTextSize(0.04);
        label->Draw();
    }
}
#endif