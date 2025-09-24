#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TRandom3.h"
#include <iostream>

// Me aparece la reconstrucción de Ex menor a Sn, voy a comprobar si es por efecto de la resolcuión aplicada en la simulación

// Masses in MeV/c^2 (from AME2020, approx.)
const double M7Li  = 7.01600455 * 931.494;   // ≈ 6533.833 MeV
const double Md    = 2.01410178 * 931.494;   // ≈ 1875.613 MeV
const double M8Li  = 8.02248736 * 931.494;   // ≈ 7471.424 MeV
const double Mp    = 938.27208816;           // proton mass
const double Mn    = 939.56542052;           // neutron mass

void debugTGenPhaseSpace()
{
    // --- Estado inicial: haz + target ---
    double Tbeam = 7 * 7.5; // MeV de energía cinética del haz 7Li
    double Ebeam = M7Li + Tbeam;
    double pz    = sqrt(Ebeam*Ebeam - M7Li*M7Li);

    TLorentzVector Pbeam(0.,0.,pz,Ebeam);
    TLorentzVector Ptarget(0.,0.,0.,Md);
    TLorentzVector Pini = Pbeam + Ptarget;

    // --- Masas finales ---
    double masses[3] = {M7Li, Mp, Mn};

    TGenPhaseSpace event;
    if(!event.SetDecay(Pini,3,masses)) {
        std::cerr << "Decay not allowed" << std::endl;
        return;
    }

    // --- Histograma Ex ---
    TH1D *hEx = new TH1D("hEx","Ex reconstruida de ^{8}Li;E_{x} [MeV];Cuentas",200,0,10);

    // --- Scatter Ep vs Theta ---
    TGraph *gEpTheta = new TGraph();
    gEpTheta->SetTitle("Energia proton vs angulo lab;Theta_p [deg];E_p [MeV]");

    int N = 50000;
    for(int i=0;i<N;++i){
        event.Generate();

        TLorentzVector *p7 = event.GetDecay(0);   // 7Li
        TLorentzVector *pp = event.GetDecay(1);   // proton
        TLorentzVector *pn = event.GetDecay(2);   // neutron

        // --- Ex reconstruida ---
        double Minv = (*p7 + *pn).M();
        double Ex   = Minv - M8Li;
        hEx->Fill(Ex);

        // --- Proton lab ---
        double Ep = pp->E() - Mp;
        double theta = pp->Vect().Angle(TVector3(0,0,1))*180./M_PI;
        gEpTheta->SetPoint(gEpTheta->GetN(),theta,Ep);
    }

    // --- Linea cinematica teorica (aproximacion 2-body CM) ---
    TGraph *gTheory = new TGraph();
    TLorentzVector Ptot_CM = Pini;
    TVector3 beta_CM = Pini.BoostVector(); // boost a CM
    int Nline = 180;
    for(int i=0;i<Nline;++i){
        double theta_deg = i; // 0-179 grados
        double theta_rad = theta_deg*M_PI/180.;

        // Energia max proton en CM (proton vs par 7Li+n)
        double M_pair = M7Li + Mn;
        double Emax_CM = (Pini.M2() + Mp*Mp - M_pair*M_pair)/(2*Pini.M());
        double p_CM = sqrt(Emax_CM*Emax_CM - Mp*Mp);

        TLorentzVector pProton_CM(0.,0.,p_CM,Emax_CM);
        // rotar proton CM segun theta
        double phi = 0.0;
        double px = p_CM * sin(theta_rad) * cos(phi);
        double py = p_CM * sin(theta_rad) * sin(phi);
        double pz = p_CM * cos(theta_rad);
        pProton_CM.SetPxPyPzE(px, py, pz, Emax_CM);

        // boost al laboratorio
        pProton_CM.Boost(beta_CM);
        gTheory->SetPoint(i,theta_deg,pProton_CM.E()-Mp);
    }

    // --- Dibujar ---
    TCanvas *c1 = new TCanvas("c1","Ex + EpTheta",1200,600);
    c1->Divide(2,1);

    c1->cd(1);
    hEx->Draw();

    c1->cd(2);
    gEpTheta->SetMarkerStyle(20);
    gEpTheta->SetMarkerSize(0.5);
    gEpTheta->Draw("AP");
    gTheory->SetLineColor(kRed);
    gTheory->SetLineWidth(2);
    gTheory->Draw("L SAME");
}