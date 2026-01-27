#define FlatTreeAnalyzer_cxx
#include "FlatTreeAnalyzer.h"

#include <TH1D.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include <iomanip>
#include <sstream>
#include <iostream>
#include <vector>
#include <iterator>

using namespace std;

//Function to divide by the bin width and to get xsecs
void Reweight(TH1D* h);

//----------------------------------------//

void FlatTreeAnalyzer::Loop() {

    if (fChain == 0) return; 

    const Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    Long64_t syst_nbytes = 0, syst_nb = 0;

    // Enable ROOT read cache
    fChain->SetCacheSize(128 * 1024 * 1024);   // 128 MB cache
    fChain->AddBranchToCache("*"); 
    
    fChain->SetBranchStatus("*", 0);

    // Enable only what you actually use
    fChain->SetBranchStatus("PDGLep", 1);
    fChain->SetBranchStatus("cc", 1);
    fChain->SetBranchStatus("nfsp", 1);
    fChain->SetBranchStatus("px", 1);
    fChain->SetBranchStatus("py", 1);
    fChain->SetBranchStatus("pz", 1);
    fChain->SetBranchStatus("pdg", 1);
    fChain->SetBranchStatus("Weight", 1);
    fChain->SetBranchStatus("CosLep", 1);
    fChain->SetBranchStatus("fScaleFactor", 1);    

    double Units = 1E39;
    double A = 1.;

    int NInte = 1;
    std::vector<TString> InteractionLabels = {"","QE","MEC","RES","DIS","COH"};

    TString FileNameAndPath = "../mc_files/analyzerOutput_" + fOutputFile + ".root";
    if (fOutputFile == "AR23" && fweights != "" && findex != -1) {
        FileNameAndPath = "../mc_files/" + fweights + "_" + int_to_string(findex)
                        + "_analyzerOutput_" + fOutputFile + ".root";
    }

    TFile* file = new TFile(FileNameAndPath,"recreate");

    // Systematics setup
    TFile* syst_file = nullptr;
    TTree* syst_ttree = nullptr;
    int ntweaks;
    double tweak_responses[7];
    double paramCVWeight;

    const bool doSyst = (fOutputFile == "AR23" && fweights != "" && findex != -1);

    if (doSyst) {
        TString syst_file_path =
          "/pnfs/sbnd/persistent/users/apapadop/GENIETweakedSamples/"
          "t2k_v3_6_2_AR23_20i_00_000_sbn1/"
          "syst_14_1000060120_CC_v3_6_2_AR23_20i_00_000.ghep.root";

        syst_file = TFile::Open(syst_file_path,"readonly");
        syst_ttree = (TTree*)syst_file->Get("events");

        if (doSyst && syst_ttree) {
            syst_ttree->SetCacheSize(64 * 1024 * 1024);
            syst_ttree->AddBranchToCache("*");
        }           

        syst_ttree->SetBranchAddress("ntweaks_" + fweights, &ntweaks);
        syst_ttree->SetBranchAddress("tweak_responses_" + fweights, &tweak_responses);
        syst_ttree->SetBranchAddress("paramCVWeight_" + fweights, &paramCVWeight);

        syst_ttree->SetBranchStatus("*", 0);
        syst_ttree->SetBranchStatus("ntweaks_" + fweights, 1);
        syst_ttree->SetBranchStatus("tweak_responses_" + fweights, 1);
        syst_ttree->SetBranchStatus("paramCVWeight_" + fweights, 1);

    }    

    // Histogram setup
    TH1D* TrueSingleBinPlot[NInte];
    TH1D* TrueDeltaPTPlot[NInte];
    TH1D* TrueDeltaAlphaTPlot[NInte];
    TH1D* TrueDeltaPhiTPlot[NInte];

    for (int inte = 0; inte < NInte; ++inte) {
        TrueSingleBinPlot[inte] =
          new TH1D(InteractionLabels[inte]+"TrueSingleBinPlot",";",1,0,1);

        double DeltaPTbins[9] = {0,0.0800,0.120,0.155,0.200,0.260,0.360,0.510,1.1};
        TrueDeltaPTPlot[inte] =
          new TH1D(InteractionLabels[inte]+"TrueDeltaPTPlot",";DeltaPT",8,DeltaPTbins);

        double DeltaAlphaTbins[9] = {0,0.47000,1.0200,1.5400,1.9800,2.3400,2.6400,2.8900,3.1416};
        TrueDeltaAlphaTPlot[inte] =
          new TH1D(InteractionLabels[inte]+"TrueDeltaAlphaTPlot",";DeltaAlphaT",8,DeltaAlphaTbins);

        double DeltaPhiTbins[9] = {0,0.067,0.140,0.225,0.340,0.520,0.850,1.500,3.1416};
        TrueDeltaPhiTPlot[inte] =
          new TH1D(InteractionLabels[inte]+"TrueDeltaPhiTPlot",";DeltaPhiT",8,DeltaPhiTbins);
    }

    // Event loop
    for (Long64_t jentry = 0; jentry < nentries; ++jentry) {

        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;

        nb = fChain->GetEntry(jentry);
        nbytes += nb;

        if (doSyst) {
            syst_nb = syst_ttree->GetEntry(jentry);
            syst_nbytes += syst_nb;
        }

        if (PDGLep != 13) continue;
        if (cc != 1) continue;
        if (nfsp < 2) continue;

        const double global_weight = fScaleFactor * Units * A;
        const double syst_weight = doSyst ? tweak_responses[findex] : 1.0;
        const double weight = global_weight * Weight * syst_weight;

        double Pp = 0.0, CosP = -999.;
        double Pmu = 0.0;
        double MuonX = 0.0, MuonY = 0.0, MuonZ = 0.0;
        double ProtonX = 0.0, ProtonY = 0.0, ProtonZ = 0.0;

        int ChargedPionTagging = 0;
        int NeutralPionTagging = 0;
        int MuonTagging = 0;

        // Loop over final state particles
        for (int i = 0; i < nfsp; ++i) {

            const double p2 = px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i];
            const double pf = std::sqrt(p2);

            const int abspdg = std::abs(pdg[i]);

            if (pdg[i] == 13) {
                MuonTagging++;
                Pmu = pf;
                MuonX = px[i];
                MuonY = py[i];
                MuonZ = pz[i];
            }

            if (pdg[i] == 2212) {
                if (pf > Pp) {
                    Pp = pf;
                    CosP = pz[i] / pf;
                    ProtonX = px[i];
                    ProtonY = py[i];
                    ProtonZ = pz[i];
                }
            }

            if (abspdg == 211) ChargedPionTagging++;
            if (pdg[i] == 111) NeutralPionTagging++;

            // Early exit if pion already found
            if (ChargedPionTagging > 0 || NeutralPionTagging > 0) break;
        }

        const bool flag_0pi =
          (ChargedPionTagging == 0 && NeutralPionTagging == 0 && MuonTagging == 1);

        if (flag_0pi) {
            TrueSingleBinPlot[0]->Fill(0.5);
        }

        const bool flag_tki =
          flag_0pi && Pmu > 0.25 && CosLep > -0.6 && Pp > 0.45 && Pp < 1 && CosP > 0.4;

        if (!flag_tki) continue;

        const double muTx = MuonX;
        const double muTy = MuonY;
        const double muTmag = std::sqrt(muTx*muTx + muTy*muTy);

        const double pTx = ProtonX;
        const double pTy = ProtonY;
        const double pTmag = std::sqrt(pTx*pTx + pTy*pTy);

        const double ptX = muTx + pTx;
        const double ptY = muTy + pTy;
        const double Pt  = std::sqrt(ptX*ptX + ptY*ptY);

        if (muTmag <= 0 || Pt <= 0 || pTmag <= 0) continue;

        double dot_mu_pt = -(muTx*ptX + muTy*ptY);
        double dot_mu_p  = -(muTx*pTx + muTy*pTy);

        double cosAlpha = dot_mu_pt / (muTmag * Pt);
        if (cosAlpha >  1.0) cosAlpha =  1.0;
        if (cosAlpha < -1.0) cosAlpha = -1.0;

        double cosPhi = dot_mu_p / (muTmag * pTmag);
        if (cosPhi >  1.0) cosPhi =  1.0;
        if (cosPhi < -1.0) cosPhi = -1.0;

        double DeltaAlphaT = std::acos(cosAlpha);
        double DeltaPhiT   = std::acos(cosPhi);

        if (DeltaPhiT > 3.1416) DeltaPhiT -= 3.1416;
        if (DeltaPhiT < 0.0)    DeltaPhiT += 3.1416;

        TrueDeltaPTPlot[0]->Fill(Pt, weight);
        TrueDeltaAlphaTPlot[0]->Fill(DeltaAlphaT, weight);
        TrueDeltaPhiTPlot[0]->Fill(DeltaPhiT, weight);
    }

    // Reweight histograms
    for (int inte = 0; inte < NInte; ++inte) {
        Reweight(TrueDeltaPhiTPlot[inte]);
        Reweight(TrueDeltaPTPlot[inte]);
        Reweight(TrueDeltaAlphaTPlot[inte]);
    }

    file->cd();

    TrueSingleBinPlot[0]->Write("TrueSingleBinPlot");
    TrueDeltaPTPlot[0]->Write("TrueDeltaPTPlot");
    TrueDeltaAlphaTPlot[0]->Write("TrueDeltaAlphaTPlot");
    TrueDeltaPhiTPlot[0]->Write("TrueDeltaPhiTPlot");

    file->Close();
    fFile->Close();

    std::cout << FileNameAndPath << " processed" << std::endl;
}

//----------------------------------------//		

void Reweight(TH1D* h) {

  int NBins = h->GetXaxis()->GetNbins();

  for (int i = 0; i < NBins; i++) {

    double CurrentEntry = h->GetBinContent(i+1);
    double NewEntry = CurrentEntry / h->GetBinWidth(i+1);

    double CurrentError = h->GetBinError(i+1);
    double NewError = CurrentError / h->GetBinWidth(i+1);

    h->SetBinContent(i+1,NewEntry); 
    h->SetBinError(i+1,NewError); 
    //h->SetBinError(i+1,0.000001); 

  }

}

//----------------------------------------//	