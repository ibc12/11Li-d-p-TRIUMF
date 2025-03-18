#ifndef GenFuncs_cxx
#define GenFuncs_cxx
#include "ActSRIM.h"

#include "TProfile.h"
#include "TRandom.h"

#include "Math/Point3D.h"
#include "Math/Vector3D.h"

using XYZPoint = ROOT::Math::XYZPoint;
using XYZVector = ROOT::Math::XYZVector;

bool IsInChamber(XYZPoint point, bool inMM = true)
{
    if(inMM) // convert to pad units
        point /= 2;
    bool a {0 <= point.X() && point.X() <= 128};
    bool b {0 <= point.Y() && point.Y() <= 128};
    bool c {0 <= point.Z() && point.Z() <= 128};
    return a && b && c;
}

void FillHist(ActPhysics::SRIM* srim, const std::string& which, double Tini, const XYZPoint& start,
              const XYZVector& dir, double range, TProfile* h)
{
    // Number of iterations to implement noise and uncertainties
    int iter {10};
    // Sigma in position
    double sx {2.}; // mm
    // Sigma in ELoss, as percent
    double spel {0.1};
    // Step along range
    double rstep {0.5}; // mm
    for(int i = 0; i < iter; i++)
    {
        double Eit {Tini};
        for(double r = 0; r <= range; r += rstep)
        {
            // Position X
            auto point {start + r * dir.Unit()};
            if(!IsInChamber(point))
                break;
            // Energy loss
            auto aux {srim->Slow(which, Eit, rstep)};
            // auto eloss {srim->EvalStoppingPower(which, Eit)};
            auto eloss {Eit - aux};
            // Also randomize
            eloss = gRandom->Gaus(eloss, eloss * spel);
            // Get X of point but randomize
            auto x {point.X()};
            x = gRandom->Gaus(x, sx);
            // Fill histogram
            h->Fill(x, eloss);
            // Prepare for next iteration
            // Eit = srim->Slow(which, Eit, rstep);
            Eit = aux;
        }
    }
}
#endif
