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

#include "../ub_cc1p/flat_analyzer/Constants.h"
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
	PlotNames.push_back("TrueDeltaPTPlot");
	PlotNames.push_back("TrueDeltaAlphaTPlot");	
	PlotNames.push_back("TrueDeltaPhiTPlot");	 

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
	vector<TH2D*> cov_matrix; cov_matrix.clear();	

	//----------------------------------------//

	// Open the files and grap the relevant plots

	for (int WhichSample = 0; WhichSample < NSamples; WhichSample ++) {

		vector<TH1D*> CurrentPlotsFullUncReco; CurrentPlotsFullUncReco.clear();
		vector<TH1D*> CurrentPlotsTrue; CurrentPlotsTrue.clear();												

		// CV With Statistical Uncertainties

		if (NameOfSamples[WhichSample] == "Overlay9") { // CV with statistical uncertainties only for now


			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

				TString FileSampleName = "data_files/dptResults.root"; 
				if (PlotNames[WhichPlot] == "TrueDeltaAlphaTPlot") { FileSampleName = "data_files/datResults.root"; }
				if (PlotNames[WhichPlot] == "TrueDeltaPhiTPlot") { FileSampleName = "data_files/dphitResults.root"; }	
				TFile* f = TFile::Open(FileSampleName,"readonly"); 					

				TH1D* histFullUncReco = (TH1D*)(f->Get("Result"));
				CurrentPlotsFullUncReco.push_back(histFullUncReco);		
				
				TH2D* cov = (TH2D*)f->Get("Covariance_Matrix");
				if (!cov) {
					std::cerr << "ERROR: Could not find Covariance in file "
							<< FileSampleName << " for plot "
							<< PlotNames[WhichPlot] << std::endl;
				}
				cov_matrix.push_back(cov);

				TH1D* histTrue = nullptr;
				CurrentPlotsTrue.push_back(histTrue);																						     
		
			}

		}

		else {
		
			for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {

				TFile* f = TFile::Open("mc_files/analyzerOutput_"+NameOfSamples[WhichSample]+".root");

				TH1D* histFullUncReco = nullptr;
				CurrentPlotsFullUncReco.push_back(histFullUncReco);										

				TH1D* histTrue = (TH1D*)(f->Get(PlotNames[WhichPlot]));
				CurrentPlotsTrue.push_back(histTrue);				
		
			}

		} // end of opening the files

		PlotsFullUncReco.push_back(CurrentPlotsFullUncReco);
		PlotsTrue.push_back(CurrentPlotsTrue);	

	}

	//----------------------------------------//

	// Loop over the plots

	for (int WhichPlot = 0; WhichPlot < N1DPlots; WhichPlot ++) {	

		// -----------------------------------------------------------------------------------------------------------------------------			

		TCanvas* PlotCanvas = new TCanvas(PlotNames[WhichPlot]+"_Combined",PlotNames[WhichPlot]+"_Combined",205,34,1024,768);
		PlotCanvas->cd();
		PlotCanvas->SetBottomMargin(0.16);
		PlotCanvas->SetTopMargin(0.09);
		PlotCanvas->SetLeftMargin(0.18);
		PlotCanvas->SetRightMargin(0.05);				

		TLegend* leg = new TLegend(0.5,0.55,0.85,0.85);

		if (PlotNames[WhichPlot] == "TrueDeltaAlphaTPlot") {

			leg = new TLegend(0.2,0.6,0.55,0.9);

		}

		leg->SetBorderSize(0);
		leg->SetTextSize(TextSize);
		leg->SetTextFont(FontStyle);
		leg->SetNColumns(1);
		leg->SetMargin(0.15);			
		leg->SetFillStyle(0);   // transparent background
		leg->SetFillColor(0);   // just in case			

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
		PlotsFullUncReco[0][WhichPlot]->GetYaxis()->SetRangeUser(0.,1.6*PlotsFullUncReco[0][WhichPlot]->GetMaximum());
		PlotsFullUncReco[0][WhichPlot]->GetYaxis()->SetTitleOffset(1.3);
		PlotsFullUncReco[0][WhichPlot]->GetYaxis()->CenterTitle();	
		PlotsFullUncReco[0][WhichPlot]->GetYaxis()->SetNdivisions(8);	
		TString reduced_name = PlotNames[WhichPlot];
		reduced_name.ReplaceAll("True","");
		PlotsFullUncReco[0][WhichPlot]->GetYaxis()->SetTitle(VarLabel[reduced_name]);							

		PlotsFullUncReco[0][WhichPlot]->SetLineColor(BeamOnColor);
		PlotsFullUncReco[0][WhichPlot]->SetMarkerColor(BeamOnColor);
		PlotsFullUncReco[0][WhichPlot]->SetMarkerSize(1.);
		PlotsFullUncReco[0][WhichPlot]->SetMarkerStyle(20);
		PlotsFullUncReco[0][WhichPlot]->SetLineWidth(1);	
		PlotsFullUncReco[0][WhichPlot]->SetTitle("");							

		leg->AddEntry(PlotsFullUncReco[0][WhichPlot],"T2K data","ep");					
		PlotCanvas->cd();
			
		//------------------------------//

		PlotsFullUncReco[0][WhichPlot]->Draw("e1x0 same");

		// -----------------------------------------------------------------------------------------------------------------

		// arrays for NSamples

		double Chi2[NSamples];
		int Ndof[NSamples];
		double pval[NSamples];

		// Clones for the NSamples-1 model predictions
		// index 0 corresponds to nominal overlay / CV

		TH1D* Clone[NSamples-1];			

		for (int WhichSample = 1; WhichSample < NSamples; WhichSample++) {
		
			Clone[WhichSample-1] = PlotsTrue[WhichSample][WhichPlot];											
			Clone[WhichSample-1]->Scale(1e-39);															

			//Clone[WhichSample-1] = PlotsTrue[WhichSample][WhichPlot];				
			Clone[WhichSample-1]->SetLineColor(Colors[WhichSample]);
			Clone[WhichSample-1]->SetLineStyle(LineStyle[WhichSample]);
			Clone[WhichSample-1]->SetMarkerColor(Colors[WhichSample]);

			Clone[WhichSample-1]->SetLineWidth(3);		
			Clone[WhichSample-1]->Draw("hist same");	
			
			TH1D* mc_clone = (TH1D*)(Clone[WhichSample-1]->Clone());
			mc_clone->Scale(1e38);

			TH1D* data_clone = (TH1D*)(PlotsFullUncReco[0][WhichPlot]->Clone());
			data_clone->Scale(1e38);			

			calc_chi2(mc_clone,data_clone,cov_matrix[WhichPlot],Chi2[WhichSample],Ndof[WhichSample],pval[WhichSample]);
			TString Chi2NdofAlt = "(" + to_string_with_precision(Chi2[WhichSample],1) + "/" + TString(std::to_string(Ndof[WhichSample])) +")";

			TLegendEntry* lGenie = leg->AddEntry(Clone[WhichSample-1],Labels[WhichSample] + Chi2NdofAlt,"l");
			lGenie->SetTextColor(Colors[WhichSample]); 
			
			//----------------------------------------//

			// only for AR23, add systematics

			if (NameOfSamples[WhichSample] == "AR23") {

				int n = Clone[WhichSample-1]->GetNbinsX();
				TString cov_name = "tot_covariance_"+PlotNames[WhichPlot];
				TH2D* cov_ar23 = (TH2D*)ar23_unc_file->Get(cov_name);	
				cov_ar23->Scale(1e-80);		
			
				divide_bin_area(cov_ar23,Clone[WhichSample-1]);
				set_unc_from_cov(cov_ar23,Clone[WhichSample-1]);

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

			} // end of AR23 

			//----------------------------------------//

		}

		//----------------------------------------//

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

		PlotCanvas->SaveAs("pdf/t2k_ccNp_"+PlotNames[WhichPlot]+".pdf");
		delete PlotCanvas;

		//----------------------------------------//

	} // End of the loop over the plots

} // End of the program 
