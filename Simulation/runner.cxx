#include "TString.h"

#include <string>
#include <unordered_map>
#include <vector>

#include "./do_all_simus.cxx"

void runner(TString what = "plot", bool inspect = true)
{
    // Beam energy
    double Tbeam {11 * 7.5}; // MeV
    // Neutron and Proton phase space
    int neutronPS {0}; // number of neutrons in final state
    int protonPS {0};  // number of protons in final state
    // Particles
    std::string beam {"11Li"};
    std::string target {"2H"};
    std::string light {"1H"};
    std::string heavy {"12Li"};
    // Vector with Exs
    std::vector<double> Exs;
    if(neutronPS == 0 && protonPS == 0 && target == "2H" && light == "1H") 
        Exs = {0., 0.130, 0.435};
    else if(neutronPS == 0 && protonPS == 0)
        Exs = {0};
    else if(neutronPS == 2 && protonPS == 0)
        Exs = {(1.26642 + 0.36928) / 2}; // half value between first excited state and the S_2n
    else if(neutronPS > 0 && protonPS == 0)
        Exs = {0}; // only gs for n phase space
    else if(neutronPS == 0 && protonPS > 0)
        Exs = {0};
    else
        throw std::runtime_error("No confs with neutronPS and protonPS enabled at the same time");

    // Run simu or plot
    if(what.Contains("simu"))
    {
        for(const auto& ex : Exs)
        {
            do_all_simus(beam, target, light, heavy, neutronPS, protonPS ,Tbeam, ex, inspect);
            // auto str {TString::Format("root -l -b -x -q \'triumf.cxx(\"%s\",\"%s\",\"%s\",%f,%f,%d)\'", beam.c_str(),
            //                           target.c_str(), light.c_str(), Tbeam, ex, inspect)};
            if(inspect)
                break; // inspect: to debug simulation
        }
    }
    else
    {
        std::cout << "No plot method implemented yet" << '\n';
    }
}
