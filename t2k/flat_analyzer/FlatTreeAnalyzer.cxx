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

	//----------------------------------------//	

	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
 	Long64_t syst_nbytes = 0, syst_nb = 0;	 

	double Units = 1E39; // so that the extracted cross-section is in 10^{-38} cm^{2}
	double A = 1.; // so that we can have xsecs per nucleon

	int NInte = 6; // Interaction processes: All, QE, MEC, RES, DIS, COH
	std::vector<TString> InteractionLabels = {"","QE","MEC","RES","DIS","COH"};

	//----------------------------------------//	

  // Output file

	TString FileNameAndPath = "../mc_files/analyzerOutput_"+fOutputFile+".root";

	if (fOutputFile == "AR23" && fweights != "" && findex != -1) { 

		FileNameAndPath = "../mc_files/" + fweights + "_" + int_to_string(findex) + "_analyzerOutput_"+fOutputFile+".root";

	}

	TFile* file = new TFile(FileNameAndPath,"recreate");

	std::cout << std::endl << "------------------------------------------------" << std::endl << std::endl;
	std::cout << "File " << FileNameAndPath << " to be created" << std::endl << std::endl;
	
	//----------------------------------------//

	TFile* syst_file = nullptr;
	TTree* syst_ttree = nullptr;	

	int	ntweaks;
   	double tweak_responses[7]; // default to 7, but some knobs might have 6
   	double paramCVWeight;	

	if (fOutputFile == "AR23" && fweights != "" && findex != -1) {

		TString syst_file_path = "/pnfs/sbnd/persistent/users/apapadop/GENIETweakedSamples/t2k_v3_6_2_AR23_20i_00_000_sbn1/syst_14_1000060120_CC_v3_6_2_AR23_20i_00_000.ghep.root";
		syst_file = TFile::Open(syst_file_path,"readonly");

	  	syst_ttree = (TTree*)syst_file->Get("events");

	  	syst_ttree->SetBranchAddress("ntweaks_" + fweights, &ntweaks);	
	  	syst_ttree->SetBranchAddress("tweak_responses_" + fweights, &tweak_responses);	
	  	syst_ttree->SetBranchAddress("paramCVWeight_" + fweights, &paramCVWeight);						

	}

	//----------------------------------------//

	// Plot declaration

  TH1D* TrueSingleBinPlot[NInte];  
  TH1D* TrueDeltaPTPlot[NInte];
  TH1D* TrueDeltaAlphaTPlot[NInte];
  TH1D* TrueDeltaPhiTPlot[NInte];

	// Loop over the interaction processes

	for (int inte = 0; inte < NInte; inte++) {

    TrueSingleBinPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueSingleBinPlot",";",1,0,1);

    double DeltaPTbins[9] = {0, 0.0800,0.120, 0.155, 0.200, 0.260, 0.360, 0.510, 1.1};
    TrueDeltaPTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPTPlot",";DeltaPT",8,DeltaPTbins);

    double DeltaAlphaTbins[9] = {0, 0.47000, 1.0200, 1.5400, 1.9800, 2.3400, 2.6400, 2.8900, 3.1416};
    TrueDeltaAlphaTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaAlphaTPlot",";DeltaAlphaT",8,DeltaAlphaTbins);

    double DeltaPhiTbins[9] = {0, 0.067000, 0.14000, 0.22500, 0.34000, 0.52000, 0.85000, 1.5000, 3.1416};
    TrueDeltaPhiTPlot[inte] = new TH1D(InteractionLabels[inte]+"TrueDeltaPhiTPlot",";DeltaPhiT",8,DeltaPhiTbins);

 	} // End of the loop over the interaction processes							

	//----------------------------------------//

	// Counters

	int CounterEventsPassedSelection = 0;
	int CounterQEEventsPassedSelection = 0;
	int CounterMECEventsPassedSelection = 0;
	int CounterRESEventsPassedSelection = 0;
	int CounterDISEventsPassedSelection = 0;
	int CounterCOHEventsPassedSelection = 0;	

	//----------------------------------------//
	
	// Loop over the events
	for (Long64_t jentry=0; jentry<nentries;jentry++) {

	  //----------------------------------------//	
	
	  Long64_t ientry = LoadTree(jentry);

	  if (fOutputFile == "AR23" && fweights != "" && findex != -1) {

	  	Long64_t syst_ientry = syst_ttree->LoadTree(jentry);
	  	syst_nb = syst_ttree->GetEntry(jentry); syst_nbytes += syst_nb; 

	  }    

	  if (ientry < 0) break; nb = fChain->GetEntry(jentry); nbytes += nb;
	  if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/nentries*100. << " %"<< std::endl;

	  //----------------------------------------//	
		
	  double syst_weight = 1.;

	  if (fOutputFile == "AR23" && fweights != "" && findex != -1) {

		syst_weight = tweak_responses[findex];

	  }		

	  double weight = fScaleFactor*Units*A*Weight*syst_weight;	

	  //----------------------------------------//	

	  // Signal definition

	  if (PDGLep != 13) { continue; } // make sure that we have only a muon in the final state
	  if (cc != 1) { continue; } // make sure that we have only CC interactions		

   // TotalXs += weight;

	  int ProtonTagging = 0, ChargedPionTagging = 0, NeutralPionTagging = 0, MuonTagging = 0, ProtonTaggingTKI = 0;
    int ElectronTagging = 0, PhotonTagging = 0;

    TVector3 P3Vector(0,0,0);
    double Pp = 0;
    double CosP = -999;

    TVector3 Muon3Vector(0,0,0);
    double Pmu = 0;

	  vector <int> ProtonID; ProtonID.clear();
	  vector <int> MuonID; MuonID.clear();		

	  // Example selection with CC1p0pi (units in GeV/c)
	  // Loop over final state particles

	  for (int i = 0; i < nfsp; i++) {
		
	    double pf = TMath::Sqrt( px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);

	    if (pdg[i] == 13) {

	      MuonTagging ++;
	      MuonID.push_back(i);
        Pmu = pf;
        Muon3Vector.SetXYZ(px[i],py[i],pz[i]);

	    }

	    if (pdg[i] == 2212 && pf > 0.5 ) {

	      ProtonTagging ++;
	      ProtonID.push_back(i);

	    }
      
      // for TKI
      if (pdg[i] == 2212) {
      
        ProtonTaggingTKI++;
        if (pf > Pp) { Pp=pf; CosP=pz[i]/pf; P3Vector.SetXYZ(px[i],py[i],pz[i]);}
      
      }

	    if (fabs(pdg[i]) == 211)  {

	      ChargedPionTagging ++;

	    }

	    if (pdg[i] == 111)  {

	      NeutralPionTagging ++;

	    }

	    if (fabs(pdg[i]) == 11)  {

	      ElectronTagging ++;

	    }

	    if (fabs(pdg[i]) == 22)  {

	      PhotonTagging ++;

	    }

	  } // End of the loop over the final state particles

	  // If the signal definition is not satisfied, continue

    bool flag_0pi = false;
    bool flag_tki = false;
	  if (ChargedPionTagging == 0 && NeutralPionTagging == 0 && MuonTagging ==1) { flag_0pi=true; }
    if (ChargedPionTagging == 0 && NeutralPionTagging == 0 && MuonTagging == 1 && Pmu>0.25 && CosLep>-0.6 && Pp>0.45 && Pp<1 && CosP>0.4) { flag_tki=true; }

	  //----------------------------------------//	

	  // Classify the events based on the interaction type

	  int genie_mode = -1.;
          
    if(flag_0pi){
            
      //TotalXs0pi+=weight;
      CounterEventsPassedSelection++;

	    if (TMath::Abs(Mode) == 1) { CounterQEEventsPassedSelection++; genie_mode = 1; } // QE
	    else if (TMath::Abs(Mode) == 2) { CounterMECEventsPassedSelection++; genie_mode = 2; } // MEC
	    else if (
		   TMath::Abs(Mode) == 10 ||
		   TMath::Abs(Mode) == 11 || TMath::Abs(Mode) == 12 || TMath::Abs(Mode) == 13 ||
		   TMath::Abs(Mode) == 17 || TMath::Abs(Mode) == 22 || TMath::Abs(Mode) == 23
	  	   ) { CounterRESEventsPassedSelection++; genie_mode = 3; } // RES
	    else if (TMath::Abs(Mode) == 21 || TMath::Abs(Mode) == 26) { CounterDISEventsPassedSelection++; genie_mode = 4; } // DIS
	    else if (TMath::Abs(Mode) == 16) { CounterCOHEventsPassedSelection++; genie_mode = 5;} // COH
	    else { genie_mode = 3; }  

      if (genie_mode < 0 || genie_mode >= NInte) continue;
          
    }
	  //----------------------------------------//
    
    if(flag_tki){

      TVector3 MuonVectorTrans;
      MuonVectorTrans.SetXYZ(Muon3Vector.X(),Muon3Vector.Y(),0.);
      double MuonVectorTransMag = MuonVectorTrans.Mag();

      TVector3 MuonVectorLong;
      MuonVectorLong.SetXYZ(0.,0.,Muon3Vector.Z());
      double MuonVectorLongMag = MuonVectorLong.Mag();

      TVector3 ProtonVectorTrans;
      ProtonVectorTrans.SetXYZ(P3Vector.X(),P3Vector.Y(),0.);
      double ProtonVectorTransMag = ProtonVectorTrans.Mag();

      TVector3 ProtonVectorLong;
      ProtonVectorLong.SetXYZ(0.,0.,P3Vector.Z());
      double ProtonVectorLongMag = ProtonVectorLong.Mag();

      TVector3 PtVector = MuonVectorTrans + ProtonVectorTrans;

      double Pt = PtVector.Mag();

	  double DeltaAlphaT = -1.;
	  double DeltaPhiT = -1.;

      if (MuonVectorTransMag > 0 && Pt > 0) {

        double arg = (-MuonVectorTrans * PtVector) /
                    (MuonVectorTransMag * Pt);
        arg = std::clamp(arg, -1.0, 1.0);
        DeltaAlphaT = TMath::ACos(arg);

      } else {

          continue;

      }

      if (MuonVectorTransMag > 0 && ProtonVectorTransMag > 0) {

        double arg = (-MuonVectorTrans * ProtonVectorTrans) /
                    (MuonVectorTransMag * ProtonVectorTransMag);
        arg = std::clamp(arg, -1.0, 1.0);
        DeltaPhiT = TMath::ACos(arg);

      } else {

        continue;

      }
	  
	  if (DeltaPhiT > 3.1416) { DeltaPhiT -= 3.1416; }
      if (DeltaPhiT < 0.) { DeltaPhiT += 3.1416; }

	  TrueSingleBinPlot[0]->Fill(0.5);
	  TrueSingleBinPlot[genie_mode]->Fill(0.5);	  
      TrueDeltaPTPlot[0]->Fill(Pt,weight);
      TrueDeltaPTPlot[genie_mode]->Fill(Pt,weight);
      TrueDeltaAlphaTPlot[0]->Fill(DeltaAlphaT,weight);
      TrueDeltaAlphaTPlot[genie_mode]->Fill(DeltaAlphaT,weight);
      TrueDeltaPhiTPlot[0]->Fill(DeltaPhiT,weight);
      TrueDeltaPhiTPlot[genie_mode]->Fill(DeltaPhiT,weight);
      
    }

	  //----------------------------------------//
	
	} // End of the loop over the events

	//----------------------------------------//	

	std::cout << "Percetage of events passing the selection cuts = " << 
	double(CounterEventsPassedSelection)/ double(nentries)*100. << " %" << std::endl; std::cout << std::endl;

	std::cout << "Success percetage in selecting QE events = " << 
	double(CounterQEEventsPassedSelection)/ double(CounterEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;

	std::cout << "Success percetage in selecting MEC events = " << 
	double(CounterMECEventsPassedSelection)/ double(CounterEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;

	std::cout << "Success percetage in selecting RES events = " << 
	double(CounterRESEventsPassedSelection)/ double(CounterEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;

	std::cout << "Success percetage in selecting DIS events = " << 
	double(CounterDISEventsPassedSelection)/ double(CounterEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;

	std::cout << "Success percetage in selecting COH events = " << 
	double(CounterCOHEventsPassedSelection)/ double(CounterEventsPassedSelection)*100. << " %" << std::endl; std::cout << std::endl;	

	//----------------------------------------//	

	// Division by bin width to get the cross sections	
	// Loop over the interaction processes

	for (int inte = 0; inte < NInte; inte++) {

    Reweight(TrueDeltaPhiTPlot[inte]);
    Reweight(TrueDeltaPTPlot[inte]);
    Reweight(TrueDeltaAlphaTPlot[inte]);

	} // End of the loop over the interaction processes		

	//----------------------------------------//		

	file->cd();
	//file->Write();
	TrueSingleBinPlot[0]->Write("TrueSingleBinPlot");
	TrueDeltaPTPlot[0]->Write("TrueDeltaPTPlot");
	TrueDeltaAlphaTPlot[0]->Write("TrueDeltaAlphaTPlot");
	TrueDeltaPhiTPlot[0]->Write("TrueDeltaPhiTPlot");

  	file->Close();
	fFile->Close();

	std::cout << std::endl;
	std::cout << "File " << FileNameAndPath +" has been created created " << std::endl; 
	std::cout << std::endl;

	std::cout << std::endl << "------------------------------------------------" << std::endl << std::endl;

	//----------------------------------------//		

} // End of the program

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