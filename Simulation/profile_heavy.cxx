#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "ActSRIM.h"

void profile_heavy()
{
    // Energías iniciales en MeV
    std::vector<double> initialEnergies{10, 20, 30, 40, 50, 60, 70, 80};

    // Colors for graph
    std::vector<int> colors = {
        kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1,
        kOrange + 7, kCyan + 1, kViolet + 1, kTeal + 2
    };

    // Partículas
    std::vector<std::string> particles{"11Li", "9Li", "8Li", "7Li"};

    // Paso en mm
    double step{1.};

    // Vector de TGraph: [partícula][energía]
    std::vector<std::vector<TGraph*>> eLossGraphs(particles.size(), std::vector<TGraph*>(initialEnergies.size(), nullptr));

    for (size_t i = 0; i < particles.size(); ++i)
    {
        const std::string& p = particles[i];

        // Cargar SRIM
        auto *srim = new ActPhysics::SRIM;
        std::string path = "../SRIM files/";
        std::string gas = "900mb_CF4_95-5";
        srim->ReadTable(p, path + p + "_" + gas + ".txt");

        for (size_t j = 0; j < initialEnergies.size(); ++j)
        {
            double initialEnergy = initialEnergies[j];
            double initialRange = srim->EvalRange(p, initialEnergy);

            // Crear el TGraph
            std::string gName = "g_" + p + "_" + std::to_string((int)initialEnergy) + "MeV";
            auto *graph = new TGraph();
            graph->SetName(gName.c_str());
            graph->SetTitle((p + ", " + "All posible energies").c_str());
            graph->GetXaxis()->SetTitle("Distance (mm)");
            graph->GetYaxis()->SetTitle("dE per step (keV/mm)");

            double Eiter = initialEnergy;
            int point = 0;

            for (double r = 0; r < initialRange; r += step)
            {
                double EpostSlow = srim->Slow(p, Eiter, step);
                double eLoss = Eiter - EpostSlow;

                graph->SetPoint(point++, r, eLoss * 1000);
                Eiter = EpostSlow;
            }

            eLossGraphs[i][j] = graph;
        }

        delete srim;
    }

     // DIBUJAR: un canvas por partícula con todas sus curvas
    for (size_t i = 0; i < particles.size(); ++i)
    {
        std::string canvasName = "c_" + particles[i];
        auto *c = new TCanvas(canvasName.c_str(), particles[i].c_str(), 800, 600);
        auto *leg = new TLegend(0.7, 0.7, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.03);

        for (size_t j = 0; j < initialEnergies.size(); ++j)
        {
            TGraph* g = eLossGraphs[i][j];
            if (!g) continue;

            if (j == 0)
            {
                g->SetLineColor(colors[j % colors.size()]);
                g->SetLineWidth(2);
                g->Draw("APL");
            }
            else
            {
                g->SetLineColor(colors[j % colors.size()]);
                g->SetLineWidth(2);
                g->Draw("LP SAME");
            }

            leg->AddEntry(g, (std::to_string((int)initialEnergies[j]) + " MeV").c_str(), "l");
        }

        c->SetTitle((particles[i] + " energy loss").c_str());
        leg->Draw();
    }
}