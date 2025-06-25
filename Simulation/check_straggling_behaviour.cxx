#include "ActSRIM.h"

#include "TMath.h"

#include <iostream>

void check_straggling_behaviour()
{
    auto srim = ActPhysics::SRIM();
    srim.ReadTable("11Li", "../SRIM files/11Li_900mb_CF4_95-5.txt");
    // srim.SetStragglingLISE("11Li", "../LISE files/11Li_gas_95-5.txt");

    double Eini = 80.; // MeV
    double RangeIni = srim.EvalRange("11Li", Eini); // mm

    double step = 0.5; // mm

    double Eafter = srim.Slow("11Li", Eini, step);
    double RangeAfter = srim.EvalRange("11Li", Eafter); // mm
    double EafterStraggling = srim.SlowWithStraggling("11Li", Eini, step);

    double StragglingIni = srim.EvalLongStraggling("11Li", RangeIni);
    double StragglingAfter = srim.EvalLongStraggling("11Li", RangeAfter);

    double stepUncertanty = TMath::Sqrt(StragglingIni * StragglingIni - StragglingAfter * StragglingAfter);

    // couts
    std::cout << "Initial energy: " << Eini << " MeV" << "After energy: " << Eafter << " MeV" << std::endl;
    std::cout << "Initial range: " << RangeIni << " mm" << " After range: " << RangeAfter << " mm" << std::endl;
    std::cout << "Initial straggling: " << StragglingIni << " mm" << " After straggling: " << StragglingAfter << " mm" << std::endl;
    std::cout << "Step: " << step << " mm" << std::endl;
    std::cout << "Step uncertanty: " << stepUncertanty << " mm" << std::endl;

}