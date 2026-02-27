#define analyzer_cxx
#include "analyzer.h"

#include <TH1D.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TGraph.h>

#include <iomanip>
#include <sstream>
#include <iostream>
#include <vector>
#include <iterator>
#include <fstream>

#include "../../helper_functions.h"

using namespace std;

vector<double> ArrayNBinsDeltaPT = {0,0.1,0.2,0.3,0.4,0.5,0.6,1.0};
int NBinsDeltaPT = ArrayNBinsDeltaPT.size()-1;

vector<double> ArrayNBinsCosThetaMuSumP = {-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1};
int NBinsCosThetaMuSumP = ArrayNBinsCosThetaMuSumP.size()-1;

vector<double> ArrayNBinsCosThetaPLPR = {-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1};
int NBinsCosThetaPLPR = ArrayNBinsCosThetaPLPR.size()-1;

//----------------------------------------//

void analyzer::Loop() {

	//----------------------------------------//	

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	Long64_t syst_nbytes = 0, syst_nb = 0;	

	double Units = 1E38; // so that the extracted cross-section is in 10^{-38} cm^{2}
	double A = 40.; // so that we can have xsecs per nucleus

 	int NInte = 1; // Interaction processes: All, QE, MEC, RES, DIS, COH   
	//int NInte = 6; // Interaction processes: All, QE, MEC, RES, DIS, COH
	std::vector<TString> InteractionLabels = {"","QE","MEC","RES","DIS","COH"};

	//----------------------------------------//	

	TString FileNameAndPath = "../mc_files/analyzerOutput_"+fOutputFile+".root";

	if (fOutputFile == "AR23" && fweights != "" && findex != -1) { 

		FileNameAndPath = "../mc_files/" + fweights + "_" + int_to_string(findex) + "_analyzerOutput_"+fOutputFile+".root";

	}

	TFile* file = new TFile(FileNameAndPath,"recreate");

	//----------------------------------------//

	TFile* syst_file = nullptr;
	TTree* syst_ttree = nullptr;	

	int	ntweaks;
   	double tweak_responses[7]; // default to 7, but some knobs might have 6
   	double paramCVWeight;	

	// syst unc only for AR23
	if (fOutputFile == "AR23" && fweights != "" && findex != -1) {

		TString syst_file_path = "/pnfs/sbnd/persistent/users/apapadop/GENIETweakedSamples/v3_6_2_AR23_20i_00_000/syst_14_1000180400_CC_v3_6_2_AR23_20i_00_000.root";
		syst_file = TFile::Open(syst_file_path,"readonly");      

	  	syst_ttree = (TTree*)syst_file->Get("events");
        fChain->AddFriend(syst_ttree);

	  	fChain->SetBranchAddress("ntweaks_" + fweights, &ntweaks);	
	  	fChain->SetBranchAddress("tweak_responses_" + fweights, &tweak_responses);	
	  	fChain->SetBranchAddress("paramCVWeight_" + fweights, &paramCVWeight);						

	}

	//----------------------------------------//

	// Plot declaration

	TH1D* TrueSingleBinPlot[NInte];	
	TH1D* TrueDeltaPTPlot[NInte];
	//TH1D* TrueDeltaAlphaTPlot[NInte];
	//TH1D* TrueDeltaPhiTPlot[NInte];
	TH1D* TrueCosThetaMuSumPPlot[NInte];
	TH1D* TrueCosThetaPLPRPlot[NInte];		

	//----------------------------------------//

	// Loop over the interaction processes

	for (int inte = 0; inte < NInte; inte++) {

	  //--------------------------------------------------//

	  TrueSingleBinPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSingleBinPlot",";",1,0,1);	  
	  TrueDeltaPTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPTPlot",";",NBinsDeltaPT,&ArrayNBinsDeltaPT[0]);	
	  TrueCosThetaMuSumPPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueCosThetaMuSumPPlot",";",NBinsCosThetaMuSumP,&ArrayNBinsCosThetaMuSumP[0]);
	  TrueCosThetaPLPRPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueCosThetaPLPRPlot",";",NBinsCosThetaPLPR,&ArrayNBinsCosThetaPLPR[0]);

	  //--------------------------------------------------//

	} // End of the loop over the interaction processes							

	//----------------------------------------//

	// Counters

	int CounterEventsPassedSelection = 0;

	//----------------------------------------//
	
    const bool use_syst =
        (fOutputFile == "AR23" && fweights != "" && findex != -1);

    vector<int> ProtonID;
    vector<int> MuonID;
    ProtonID.reserve(4);
    MuonID.reserve(2);

    for (Long64_t jentry = 0; jentry < nentries; ++jentry) {

        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;

        fChain->GetEntry(jentry);

        double syst_weight = 1.0;
        if (use_syst) syst_weight = tweak_responses[findex];
        //divide by the MvA CV
        if (use_syst && fweights.Contains("MvA")) { 
            // tweak_responses[6] has the new CV
            syst_weight = syst_weight/tweak_responses[6]; 
        
        }

        double weight = fScaleFactor * Units * A * Weight * syst_weight;
        if (std::isnan(weight)) { continue; }

        // ---------------- Signal selection ----------------

        if (PDGLep != 13) continue;

        int ProtonTagging = 0, ChargedPionTagging = 0;
        int NeutralPionTagging = 0, MuonTagging = 0;

        ProtonID.clear();
        MuonID.clear();

        for (int i = 0; i < nfsp; ++i) {

            const double px_i = px[i];
            const double py_i = py[i];
            const double pz_i = pz[i];

            const double pf2 = px_i*px_i + py_i*py_i + pz_i*pz_i;
            const double pf  = std::sqrt(pf2);

            const int pdg_i = pdg[i];

            if (pdg_i == 13 && (pf > 0.1 && pf < 1.2)) {
                ++MuonTagging;
                MuonID.push_back(i);
            }

            if (pdg_i == 2212 && (pf > 0.3 && pf < 1.0)) {
                ++ProtonTagging;
                ProtonID.push_back(i);
            }

            if (std::abs(pdg_i) == 211 && pf > 0.065) {
                ++ChargedPionTagging;
            }

            if (pdg_i == 111) {
                ++NeutralPionTagging;
            }
        }

        if (!(ProtonTagging == 2 &&
            ChargedPionTagging == 0 &&
            NeutralPionTagging == 0 &&
            MuonTagging == 1)) continue;

        // ---------------- Interaction classification ----------------

        int genie_mode = -1;
        if (fOutputFile == "ACHILLES") {
            genie_mode = 1;
        } else {
            const int m = std::abs(Mode);
            if      (m == 1) genie_mode = 1;
            else if (m == 2) genie_mode = 2;
            else if (m == 10 || m == 11 || m == 12 ||
                    m == 13 || m == 17 || m == 22 || m == 23)
                genie_mode = 3;
            else if (m == 21 || m == 26)
                genie_mode = 4;
            else if (m == 16)
                genie_mode = 5;
            else
                continue;
        }

        ++CounterEventsPassedSelection;

        // ---------------- Kinematics (no ROOT objects) ----------------

        const int imu = MuonID[0];
        const int ip1 = ProtonID[0];
        const int ip2 = ProtonID[1];

        const double px_mu = px[imu], py_mu = py[imu], pz_mu = pz[imu];
        const double px_p1 = px[ip1], py_p1 = py[ip1], pz_p1 = pz[ip1];
        const double px_p2 = px[ip2], py_p2 = py[ip2], pz_p2 = pz[ip2];

        const double px_sum = px_p1 + px_p2;
        const double py_sum = py_p1 + py_p2;
        const double pz_sum = pz_p1 + pz_p2;

        const double protomn_sum_mag2 =
            px_sum*px_sum + py_sum*py_sum + pz_sum*pz_sum;
        const double muon_mag2 =
            px_mu*px_mu + py_mu*py_mu + pz_mu*pz_mu;

        if (protomn_sum_mag2 == 0 || muon_mag2 == 0) continue;

        const double DeltaPT = std::sqrt(px_sum*px_sum + py_sum*py_sum);

        const double mu_dot_sum =
            px_mu*px_sum + py_mu*py_sum + pz_mu*pz_sum;

        const double costheta_mu_sump =
            mu_dot_sum / std::sqrt(muon_mag2 * protomn_sum_mag2);

        const double pl_dot_pr =
            px_p1*px_p2 + py_p1*py_p2 + pz_p1*pz_p2;

        const double pl_mag2 =
            px_p1*px_p1 + py_p1*py_p1 + pz_p1*pz_p1;
        const double pr_mag2 =
            px_p2*px_p2 + py_p2*py_p2 + pz_p2*pz_p2;

        if (pl_mag2 == 0 || pr_mag2 == 0) continue;

        const double costheta_pl_pr =
            pl_dot_pr / std::sqrt(pl_mag2 * pr_mag2);

        double DeltaPT_fill = DeltaPT;
        if (DeltaPT_fill > ArrayNBinsDeltaPT[NBinsDeltaPT]) {
            DeltaPT_fill =
                0.5 * (ArrayNBinsDeltaPT[NBinsDeltaPT] +
                    ArrayNBinsDeltaPT[NBinsDeltaPT - 1]);
        }

        // ---------------- Fill histograms ----------------

        TrueSingleBinPlot[0]->Fill(0.5, weight);
        TrueDeltaPTPlot[0]->Fill(DeltaPT_fill, weight);
        TrueCosThetaMuSumPPlot[0]->Fill(costheta_mu_sump, weight);
        TrueCosThetaPLPRPlot[0]->Fill(costheta_pl_pr, weight);

    }

	//----------------------------------------//		

	file->cd();
	//file->Write();

	TrueSingleBinPlot[0]->Write();	
	TrueDeltaPTPlot[0]->Write();
	// TrueDeltaAlphaTPlot[0]->Write();
	// TrueDeltaPhiTPlot[0]->Write();
	TrueCosThetaMuSumPPlot[0]->Write();
	TrueCosThetaPLPRPlot[0]->Write();		

    for (int inte = 0; inte < NInte; ++inte) {
        delete TrueSingleBinPlot[inte];
        delete TrueDeltaPTPlot[inte];
        delete TrueCosThetaMuSumPPlot[inte];
        delete TrueCosThetaPLPRPlot[inte];
    }    

	file->Close();	
	delete file;

    if (syst_file) {

        fChain->RemoveFriend(syst_ttree);
        syst_file->Close();
        delete syst_file;
        syst_file = nullptr;
    
    }  

	std::cout << endl << FileNameAndPath +" processed" << std::endl; 

	//----------------------------------------//		

} // End of the program