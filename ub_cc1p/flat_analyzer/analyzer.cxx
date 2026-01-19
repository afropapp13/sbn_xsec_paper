#define analyzer_cxx
#include "analyzer.h"

#include <TH1D.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TLorentzVector.h>

#include <iostream>
#include <vector>
#include <cmath>

#include "Constants.h"
#include "STV_Tools.h"
#include "Tools.h"
#include "../../helper_functions.h"

using namespace std;
using namespace Constants;

//--------------------------------------------------
void Reweight(TH1D* h, double SF = 1.)
{
    const int nb = h->GetNbinsX();
    for (int i = 1; i <= nb; ++i) {
        const double w = h->GetBinWidth(i);
        h->SetBinContent(i, SF * h->GetBinContent(i) / w);
        h->SetBinError  (i, SF * h->GetBinError(i)   / w);
    }
}
//--------------------------------------------------

void analyzer::Loop()
{
    if (!fChain) return;

    const Long64_t nentries = fChain->GetEntriesFast();

    const double Units = 1E38;
    const double A = 40.;

    Tools tools;

    //--------------------------------------
    // Output
    TString outname = "../mc_files/analyzerOutput_" + fOutputFile + ".root";

    if (fOutputFile == "AR23" && fweights != "" && findex != -1)
        outname = "../mc_files/" + fweights + "_" +
                  int_to_string(findex) +
                  "_analyzerOutput_" + fOutputFile + ".root";

    TFile* outfile = TFile::Open(outname, "RECREATE");

    //--------------------------------------
  
    TFile* syst_file = nullptr;
    TTree* syst_tree = nullptr;

    double tweak_responses[7] = {1,1,1,1,1,1,1};
    int ntweaks = 0;

    if (fOutputFile == "AR23" && fweights != "" && findex != -1) {

        syst_file = TFile::Open(
            "../mc_files/syst_14_1000180400_CC_v3_6_2_AR23_20i_00_000.root",
            "READonly");

        syst_tree = (TTree*)syst_file->Get("events");

        syst_tree->SetBranchStatus("*", 0);
        syst_tree->SetBranchStatus(("tweak_responses_" + fweights).Data(), 1);
        syst_tree->SetBranchStatus(("ntweaks_" + fweights).Data(), 1);

        syst_tree->SetBranchAddress(("ntweaks_" + fweights).Data(), &ntweaks);
        syst_tree->SetBranchAddress(("tweak_responses_" + fweights).Data(), &tweak_responses);

    }

    //--------------------------------------

    outfile->cd();
    
    TH1D* h_mu_ct = new TH1D(
        "TrueMuonCosThetaPlot",
        ";cos(#theta_{#mu})",
        NBinsMuonCosTheta,
        ArrayNBinsMuonCosTheta);

    TH1D* h_singlebin = new TH1D(
        "TrueSingleBinPlot",
        "",
        1,0,1);

    TH1D* h_dpt = new TH1D(
        "TrueDeltaPTPlot",
        LabelXAxisDeltaPT,
        NBinsDeltaPT,
        ArrayNBinsDeltaPT);

    TH1D* h_dat = new TH1D(
        "TrueDeltaAlphaTPlot",
        LabelXAxisDeltaAlphaT,
        NBinsDeltaAlphaT,
        ArrayNBinsDeltaAlphaT);

    TH1D* h_dphit = new TH1D(
        "TrueDeltaPhiTPlot",
        LabelXAxisDeltaPhiT,
        NBinsDeltaPhiT,
        ArrayNBinsDeltaPhiT);

    h_mu_ct->Sumw2();
    h_singlebin->Sumw2();
    h_dpt->Sumw2();
    h_dat->Sumw2();
    h_dphit->Sumw2();

    //--------------------------------------

    int passed = 0;

    //--------------------------------------
    // EVENT LOOP
    for (Long64_t i = 0; i < nentries; ++i) {

        fChain->GetEntry(i);

        if (PDGLep != 13) continue;
        if (fOutputFile != "ACHILLES" && cc != 1) continue;

        int mu = -1, pr = -1;
        int nmu = 0, npr = 0, nchpi = 0, npi0 = 0, nheavy = 0;

        for (int j = 0; j < nfsp; ++j) {
            const double p = std::sqrt(px[j]*px[j] + py[j]*py[j] + pz[j]*pz[j]);

            if (pdg[j] == 13 && p > 0.1 && p < 1.2) { mu = j; ++nmu; }
            else if (pdg[j] == 2212 && p > 0.3 && p < 1.0) { pr = j; ++npr; }
            else if (abs(pdg[j]) == 211 && p > 0.07) ++nchpi;
            else if (pdg[j] == 111) ++npi0;
            else if (pdg[j] != NeutralPionPdg &&
                     abs(pdg[j]) != AbsChargedPionPdg &&
                     tools.is_meson_or_antimeson(pdg[j])) ++nheavy;
        }

        if (nmu != 1 || npr != 1 || nchpi || npi0 || nheavy) continue;

        ++passed;

        //----------------------------------        

        double syst_weight = 1.0;

        if (findex >= 0) {
            
            syst_tree->GetEntry(i);    
            //cout << "findex: " << findex << ", ntweaks: " << ntweaks << endl;            
            syst_weight = tweak_responses[findex];
            //cout << "syst_weight: " << syst_weight << endl;

        }

        const double weight =
            fScaleFactor * Units * A * Weight * syst_weight;

        TLorentzVector mu4(px[mu], py[mu], pz[mu], E[mu]);
        TLorentzVector pr4(px[pr], py[pr], pz[pr], E[pr]);

        STV_Tools stv(
            mu4.Vect(),
            pr4.Vect(),
            mu4.E(),
            sqrt(pr4.Rho()*pr4.Rho() +
                 ProtonMass_GeV*ProtonMass_GeV));

        double dpt = stv.ReturnPt();
        if (dpt > ArrayNBinsDeltaPT[NBinsDeltaPT])
            dpt = 0.9999 * ArrayNBinsDeltaPT[NBinsDeltaPT];

        h_mu_ct->Fill(mu4.CosTheta(), weight);
        h_singlebin->Fill(0.5, weight);
        h_dpt->Fill(dpt, weight);
        h_dat->Fill(stv.ReturnDeltaAlphaT(), weight);
        h_dphit->Fill(stv.ReturnDeltaPhiT(), weight);
        
    }

    //--------------------------------------

    // no bin width division here bc we need to multiply by Ac first
    // Reweight(h_mu_ct);
    // Reweight(h_dpt);
    // Reweight(h_dat);
    // Reweight(h_dphit);

    outfile->Write();
    outfile->Close();
    if (syst_file) syst_file->Close();

    cout << outname << " processed" << endl;

}