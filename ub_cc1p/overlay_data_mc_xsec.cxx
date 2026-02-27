#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TMath.h>
#include <TLatex.h>
#include <TMatrixD.h>
#include <TVectorD.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "flat_analyzer/Constants.h"
#include "../helper_functions.h"

using namespace std;
using namespace Constants;

//----------------------------------------//

void overlay_data_mc_xsec() {

	//----------------------------------------//

	int DecimalAccuracy = 2;

	TH1D::SetDefaultSumw2();
	gStyle->SetEndErrorSize(6);		
	gStyle->SetPalette(55); 
	const Int_t NCont = 999; 
	gStyle->SetNumberContours(NCont); 
	gStyle->SetTitleSize(0.07,"t");
	gStyle->SetOptStat(0);

	//----------------------------------------//

	// AR23 uncertainty file

	TFile* ar23_unc_file = TFile::Open("mc_files/covariances.root","readonly");

	//----------------------------------------//

	vector<TString> PlotNames;
	//PlotNames.push_back("SingleBinPlot");	
	PlotNames.push_back("MuonCosThetaPlot");
	PlotNames.push_back("DeltaPTPlot");
	PlotNames.push_back("DeltaAlphaTPlot");	 

	const int N1DPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << N1DPlots << endl;

	//----------------------------------------//

	vector<vector<TH1D*> > PlotsFullUncReco; PlotsFullUncReco.clear();
	vector<vector<TH1D*> > PlotsTrue; PlotsTrue.clear();						

	vector<TString> NameOfSamples; NameOfSamples.clear();
	vector<int> Colors; Colors.clear();		
	vector<TString> Labels; Labels.clear();
	vector<int> LineStyle; LineStyle.clear();
	vector<TString> weights; weights.clear();
	
	// CV

	NameOfSamples.push_back("Overlay9"); Colors.push_back(OverlayColor); Labels.push_back("G18T "); LineStyle.push_back(G18LineStyle); weights.push_back("");

	//----------------------------------------//	

	NameOfSamples.push_back("AR23"); Colors.push_back(kOrange+10); Labels.push_back("AR23 ");LineStyle.push_back(kSolid); weights.push_back("");
	NameOfSamples.push_back("AR25"); Colors.push_back(kGreen+2); Labels.push_back("AR25 "); LineStyle.push_back(kSolid); weights.push_back("");
	//NameOfSamples.push_back("AR23_Minerva"); Colors.push_back(kBlue+1); Labels.push_back("AR23 Minerva");LineStyle.push_back(kSolid); weights.push_back("");
	//NameOfSamples.push_back("AR23_LQCD"); Colors.push_back(kMagenta+1); Labels.push_back("AR23 LQCD");LineStyle.push_back(kSolid); weights.push_back("");
	//NameOfSamples.push_back("G18"); Colors.push_back(kGreen+2); Labels.push_back("G18");LineStyle.push_back(kSolid); weights.push_back("");			

	//----------------------------------------//

	const int NSamples = NameOfSamples.size();
	vector<TFile*> FileSample; FileSample.clear();

	//----------------------------------------//

	// Open the files and grap the relevant plots

	for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

		vector<TH1D*> CurrentPlotsFullUncReco; CurrentPlotsFullUncReco.clear();
		vector<TH1D*> CurrentPlotsTrue; CurrentPlotsTrue.clear();												

		// CV With Statistical Uncertainties

		if (NameOfSamples[WhichSample] == "Overlay9") { // CV with statistical uncertainties only for now

			TString FileSampleName = "data_files/All_XSecs_Combined_v08_00_00_52.root"; 
			FileSample.push_back(TFile::Open(FileSampleName,"readonly")); 

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

				TH1D* histFullUncReco = (TH1D*)(FileSample[WhichSample]->Get("FullUnc_"+PlotNames[WhichPlot]));
				CurrentPlotsFullUncReco.push_back(histFullUncReco);										

				TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get("OverlayGENIE_"+PlotNames[WhichPlot]));
				CurrentPlotsTrue.push_back(histTrue);																						     
		
			}

		}

		else {
		
		  FileSample.push_back(TFile::Open("mc_files/analyzerOutput_"+NameOfSamples[WhichSample]+".root")); 

			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

				TH1D* histFullUncReco = nullptr;
				CurrentPlotsFullUncReco.push_back(histFullUncReco);										

				TH1D* histTrue = (TH1D*)(FileSample[WhichSample]->Get("True"+PlotNames[WhichPlot]));
				CurrentPlotsTrue.push_back(histTrue);				
		
			}

		} // end of opening the files

		PlotsFullUncReco.push_back(CurrentPlotsFullUncReco);
		PlotsTrue.push_back(CurrentPlotsTrue);	

	}

	//----------------------------------------//

	// Loop over the plots

	for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {	

		//----------------------------------------//

		TH2D* Ac = (TH2D*)FileSample[0]->Get("Ac_"+PlotNames[WhichPlot]);

		TString CovString = "UnfCov_"+PlotNames[WhichPlot];

		// already includes the bin area division
		TH2D* Cov = (TH2D*)FileSample[0]->Get(CovString);	

		// -----------------------------------------------------------------------------------------------------------------------------			

		TCanvas* PlotCanvas = new TCanvas(PlotNames[WhichPlot]+"_Combined",PlotNames[WhichPlot]+"_Combined",205,34,1024,768);
		PlotCanvas->cd();
		PlotCanvas->SetBottomMargin(0.14);
		PlotCanvas->SetTopMargin(0.09);
		PlotCanvas->SetLeftMargin(0.18);
		PlotCanvas->SetRightMargin(0.05);				

		TLegend* leg = new TLegend(0.2,0.55,0.55,0.85);

		if (PlotNames[WhichPlot] == "DeltaPTPlot") {

			leg = new TLegend(0.4,0.55,0.75,0.85);

		}

		leg->SetBorderSize(0);
		leg->SetTextSize(TextSize);
		leg->SetTextFont(FontStyle);
		leg->SetNColumns(1);
		leg->SetMargin(0.15);				

		// ------------------------------------------------------------------------------------------------------------------

		// BeamOn Total Uncertainty

		double MaxValue = PlotsFullUncReco[0][WhichPlot]->GetMaximum();
		double MinValue = PlotsFullUncReco[0][WhichPlot]->GetMinimum();

		PlotsFullUncReco[0][WhichPlot]->GetXaxis()->SetTitleFont(FontStyle);			
		PlotsFullUncReco[0][WhichPlot]->GetXaxis()->SetLabelFont(FontStyle);					
		PlotsFullUncReco[0][WhichPlot]->GetXaxis()->SetLabelSize(TextSize);				
		PlotsFullUncReco[0][WhichPlot]->GetXaxis()->SetTitleSize(TextSize);	
		PlotsFullUncReco[0][WhichPlot]->GetXaxis()->SetTitleOffset(1.1);	
		PlotsFullUncReco[0][WhichPlot]->GetXaxis()->CenterTitle();	
		PlotsFullUncReco[0][WhichPlot]->GetXaxis()->SetNdivisions(8);								

		PlotsFullUncReco[0][WhichPlot]->GetYaxis()->SetTitleFont(FontStyle);			
		PlotsFullUncReco[0][WhichPlot]->GetYaxis()->SetLabelFont(FontStyle);					
		PlotsFullUncReco[0][WhichPlot]->GetYaxis()->SetLabelSize(TextSize);				
		PlotsFullUncReco[0][WhichPlot]->GetYaxis()->SetTitleSize(TextSize);			
		PlotsFullUncReco[0][WhichPlot]->GetYaxis()->SetRangeUser(XSecRange[PlotNames[WhichPlot]].first,XSecRange[PlotNames[WhichPlot]].second);
		PlotsFullUncReco[0][WhichPlot]->GetYaxis()->SetTitleOffset(1.3);
		PlotsFullUncReco[0][WhichPlot]->GetYaxis()->CenterTitle();	
		PlotsFullUncReco[0][WhichPlot]->GetYaxis()->SetNdivisions(8);					

		PlotsFullUncReco[0][WhichPlot]->SetLineColor(BeamOnColor);
		PlotsFullUncReco[0][WhichPlot]->SetMarkerColor(BeamOnColor);
		PlotsFullUncReco[0][WhichPlot]->SetMarkerSize(1.);
		PlotsFullUncReco[0][WhichPlot]->SetMarkerStyle(20);
		PlotsFullUncReco[0][WhichPlot]->SetLineWidth(1);	
		PlotsFullUncReco[0][WhichPlot]->GetYaxis()->SetTitle(VarLabel[PlotNames[WhichPlot]]);					

		leg->AddEntry(PlotsFullUncReco[0][WhichPlot],"MicroBooNE Data","ep");					
		PlotCanvas->cd();
			
		//------------------------------//

		PlotsFullUncReco[0][WhichPlot]->Draw("e1x0 same");

		// -----------------------------------------------------------------------------------------------------------------

		// Overlay GENIE v3 + Tune

		PlotsTrue[0][WhichPlot]->SetLineColor(Colors[0]);
		PlotsTrue[0][WhichPlot]->SetMarkerColor(Colors[0]);
		PlotsTrue[0][WhichPlot]->SetLineStyle(LineStyle[0]);

		// -----------------------------------------------------------------------------------------------------------------

		// arrays for NSamples

		double Chi2[NSamples];
		int Ndof[NSamples];
		double pval[NSamples];

		// Clones for the NSamples-1 model predictions
		// index 0 corresponds to nominal overlay / CV

		TH1D* Clone[NSamples-1];			

		for (int WhichSample = 1; WhichSample < NSamples; WhichSample++) {

			// Apply the additional smearing matrix Ac
			Clone[WhichSample-1] = Multiply(PlotsTrue[WhichSample][WhichPlot],Ac);											
			// Divide by the bin width
			divide_bin_width(Clone[WhichSample-1],1.);																

			//Clone[WhichSample-1] = PlotsTrue[WhichSample][WhichPlot];				
			Clone[WhichSample-1]->SetLineColor(Colors[WhichSample]);
			Clone[WhichSample-1]->SetLineStyle(LineStyle[WhichSample]);
			Clone[WhichSample-1]->SetMarkerColor(Colors[WhichSample]);

			Clone[WhichSample-1]->SetLineWidth(3);		
			Clone[WhichSample-1]->Draw("hist same");		
			
			//----------------------------------------//

			TString Chi2NdofAlt = "";

			// only for AR23, add systematics

			if (NameOfSamples[WhichSample] == "AR23") {

				int n = Clone[WhichSample-1]->GetNbinsX();
				TString cov_name = "tot_covariance_True"+PlotNames[WhichPlot];
				TH2D* cov_ar23 = (TH2D*)ar23_unc_file->Get(cov_name);			

				TMatrixD mat_cov_ar23(n,n);
				//diagonal runs lower-left → upper-right
				H2M(cov_ar23, mat_cov_ar23,kFALSE);		
				//PrintMatrix(mat_cov_ar23);

				TMatrixD mat_cov_ar23_fixed(n, n);

				for (int i = 0; i < n; ++i) {
  
					for (int j = 0; j < n; ++j) {
    
						mat_cov_ar23_fixed(i, j) = mat_cov_ar23(n - 1 - i, j);
  
					}

				}

				//diagonal runs upper-left → top-right
				mat_cov_ar23 = mat_cov_ar23_fixed;

				// TCanvas c_mat_cov_ar23("c_mat_cov_ar23", "c_mat_cov_ar23", 800, 800);
				// mat_cov_ar23.Draw("COLZ TEXT");   // COLZ = color map, TEXT = numbers
				// c_mat_cov_ar23.SaveAs("debug_pdf/c_mat_cov_ar23_"+PlotNames[WhichPlot]+".pdf");

				TMatrixD trans_mat_cov_ar23(n,n);
				//trans_mat_cov_ar23.Transpose(mat_cov_ar23);
				trans_mat_cov_ar23 = mat_cov_ar23;						

				//TCanvas c_trans_mat_cov_ar23("c_trans_mat_cov_ar23", "c_trans_mat_cov_ar23", 800, 800);
				//trans_mat_cov_ar23.Draw("COLZ TEXT");   // COLZ = color map, TEXT = numbers
				//c_trans_mat_cov_ar23.SaveAs("debug_pdf/c_trans_mat_cov_ar23_"+PlotNames[WhichPlot]+".pdf");				

				TMatrixD mat_Ac(n,n);
				H2M(Ac, mat_Ac,kFALSE);		
				
				TMatrixD mat_Ac_fixed(n, n);

				for (int i = 0; i < n; ++i) {
  
					for (int j = 0; j < n; ++j) {
    
						mat_Ac_fixed(i, j) = mat_Ac(n - 1 - i, j);
  
					}

				}

				//diagonal runs upper-left → top-right
				mat_Ac = mat_Ac_fixed;				

				//TCanvas c_mat_ac_ar23("c_mat_ac_ar23", "c_mat_ac_ar23", 800, 800);
				//mat_Ac.Draw("COLZ TEXT");   // COLZ = color map, TEXT = numbers
				//c_mat_ac_ar23.SaveAs("debug_pdf/c_mat_ac_ar23_"+PlotNames[WhichPlot]+".pdf");		

				TMatrixD trans_mat_Ac(n,n);

				for (int i = 0; i < n; ++i) {
				
					for (int j = 0; j < n; ++j) {
					
						trans_mat_Ac(i, j) = mat_Ac(n - 1 - j, n - 1 - i);
					
					}
				
				}		

				//TCanvas c_trans_mat_ac_ar23("c_trans_mat_ac_ar23", "c_trans_mat_ac_ar23", 800, 800);
				//trans_mat_Ac.Draw("COLZ TEXT");   // COLZ = color map, TEXT = numbers
				//c_trans_mat_ac_ar23.SaveAs("debug_pdf/c_trans_mat_ac_ar23_"+PlotNames[WhichPlot]+".pdf");						

				// Applying the smearing: Cov' = A * Cov^T * A^T
				TMatrixD smeared_mat_cov_ar23 = mat_Ac * trans_mat_cov_ar23 * trans_mat_Ac;
				//TMatrixD smeared_mat_cov_ar23 = trans_mat_cov_ar23;	// wrong !!! fix it			

				//TCanvas c_smeared_mat_cov_ar23("c_smeared_mat_cov_ar23", "c_smeared_mat_cov_ar23", 800, 800);
				//smeared_mat_cov_ar23.Draw("COLZ TEXT");   // COLZ = color map, TEXT = numbers
				//c_smeared_mat_cov_ar23.SaveAs("debug_pdf/c_smeared_mat_cov_ar23_"+PlotNames[WhichPlot]+".pdf");						

				TH2D* cov_ar23_smeared = (TH2D*)cov_ar23->Clone();
				flip_rows(smeared_mat_cov_ar23, cov_ar23_smeared);
				//TCanvas c_smeared_th2d_cov_ar23("c_smeared_th2d_cov_ar23", "c_smeared_th2d_cov_ar23", 800, 800);
				//cov_ar23_smeared->Draw("COLZ TEXT");   // COLZ = color map, TEXT = numbers
				//c_smeared_th2d_cov_ar23.SaveAs("debug_pdf/c_smeared_th2d_cov_ar23_"+PlotNames[WhichPlot]+".pdf");						

				divide_bin_area(cov_ar23_smeared,Clone[WhichSample-1]);
				set_unc_from_cov(cov_ar23_smeared,Clone[WhichSample-1]);

				//----------------------//

				// error band on AR23 prediction

				TH1D* h = Clone[WhichSample-1];

				// Build graph with symmetric errors
				TGraphAsymmErrors* gerr = new TGraphAsymmErrors(h);

				for (int i = 0; i < gerr->GetN(); ++i) {
					double x, y;
					gerr->GetPoint(i, x, y);

					// True symmetric error
					double eyup  = gerr->GetErrorYhigh(i);
					double eydn  = gerr->GetErrorYlow(i);
					double ey    = 0.5 * (eyup + eydn);  // enforce symmetry

					// Keep point at central value
					gerr->SetPoint(i, x, y);

					// Make symmetric around y
					gerr->SetPointEYhigh(i, ey);
					gerr->SetPointEYlow(i, ey);
				}

				// Style only the symmetric band
				gerr->SetFillColor(Colors[WhichSample]);
				gerr->SetFillStyle(3004);
				gerr->SetLineWidth(0);

				// Draw histogram line
				h->SetLineColor(Colors[WhichSample]);
				h->SetFillStyle(0);
				h->Draw("hist same");

				// Draw symmetric error band ONLY
				gerr->Draw("E2 same");

				TH2D* cov_xsec = (TH2D*)cov_ar23_smeared->Clone();

				// Clone Cov ONLY to inherit binning
				TH2D* cov_xsec_rebinned = (TH2D*)Cov->Clone("cov_xsec_rebinned");
				cov_xsec_rebinned->Reset();  // remove contents

				// Copy bin contents by bin number
				for (int i = 1; i <= Cov->GetNbinsX(); ++i) {
					for (int j = 1; j <= Cov->GetNbinsY(); ++j) {
						cov_xsec_rebinned->SetBinContent(i, j,
							cov_xsec->GetBinContent(i, j));
					}
				}				

				cov_xsec_rebinned->Add(Cov);
				// nominal covariance
				calc_chi2(Clone[WhichSample-1],PlotsFullUncReco[0][WhichPlot],Cov,Chi2[WhichSample],Ndof[WhichSample],pval[WhichSample]);
				double chi2_nom = Chi2[WhichSample];
				// nominal + xsec covariance
				calc_chi2(Clone[WhichSample-1],PlotsFullUncReco[0][WhichPlot],cov_xsec_rebinned,Chi2[WhichSample],Ndof[WhichSample],pval[WhichSample]);
				double chi2_full = Chi2[WhichSample];
				Chi2NdofAlt = "(" + to_string_with_precision(chi2_nom,1) + "/" + TString(std::to_string(Ndof[WhichSample])) +")"\
				 + " [" + to_string_with_precision(chi2_full,1)+ + "/" + TString(std::to_string(Ndof[WhichSample])) +"]";

			}

			else {  

				calc_chi2(Clone[WhichSample-1],PlotsFullUncReco[0][WhichPlot],Cov,Chi2[WhichSample],Ndof[WhichSample],pval[WhichSample]);
				Chi2NdofAlt = "(" + to_string_with_precision(Chi2[WhichSample],1) + "/" + TString(std::to_string(Ndof[WhichSample])) +")";

			}

			TLegendEntry* lGenie = leg->AddEntry(Clone[WhichSample-1],Labels[WhichSample] + Chi2NdofAlt,"l");
			lGenie->SetTextColor(Colors[WhichSample]); 			

			//----------------------------------------//

		}

		//----------------------------------------//

		// chi2
		calc_chi2(PlotsTrue[0][WhichPlot],PlotsFullUncReco[0][WhichPlot],Cov,Chi2[0],Ndof[0],pval[0]);
		TString Chi2NdofNom = "(" + to_string_with_precision(Chi2[0],1) + "/" + TString(std::to_string(Ndof[0])) +")";

		TLegendEntry* lGenie_GenieOverlay = leg->AddEntry(PlotsTrue[0][WhichPlot],Labels[0]+Chi2NdofNom,"l");
		PlotsTrue[0][WhichPlot]->SetLineWidth(3); 
		PlotsTrue[0][WhichPlot]->Draw("hist same"); 
		lGenie_GenieOverlay->SetTextColor(Colors[0]); 

		// ---------------------------------------------------------------------------------------------------------
		// ---------------------------------------------------------------------------------------------------------

		PlotsFullUncReco[0][WhichPlot]->Draw("e1x0 same"); // BeamOn Stat Total
		leg->Draw();					

		TLatex *textSlice = new TLatex();
		textSlice->SetTextFont(FontStyle);
		textSlice->SetTextSize(TextSize);
		TString PlotNameDuplicate = PlotNames[WhichPlot];
		TString ReducedPlotName = PlotNameDuplicate.ReplaceAll("Reco","") ;
		textSlice->DrawLatexNDC(0.24, 0.92, LatexLabel[ ReducedPlotName ].ReplaceAll("All events","") );

		//----------------------------------------//

		// Saving the canvas with the data (total uncertainties) vs overlay & generator predictions

		PlotCanvas->SaveAs("pdf/ub_cc1p_"+PlotNames[WhichPlot]+".pdf");
		delete PlotCanvas;

		//----------------------------------------//

	} // End of the loop over the plots

} // End of the program 
