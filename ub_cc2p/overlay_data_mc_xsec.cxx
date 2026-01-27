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

#include "../helper_functions.h"

using namespace std;

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

	//----------------------------------------//

	// AR23 uncertainty file

	TFile* ar23_unc_file = TFile::Open("mc_files/covariances.root","readonly");

	//----------------------------------------//

	// data file

	TString fdata_name = "data_files/NuisanceInput_WSVD_2ndDeriv.root"; 
	TFile* fdata = TFile::Open(fdata_name,"readonly"); 	

	//----------------------------------------//

	vector<TString> plot_name;
	plot_name.push_back("TrueDeltaPTPlot");
	plot_name.push_back("TrueCosThetaMuSumPPlot");
	plot_name.push_back("TrueCosThetaPLPRPlot");	 

	const int nplots = plot_name.size();
	cout << "Number of 1D Plots = " << nplots << endl;

	//----------------------------------------//

	//1st index: plot 
	vector<TH1D*> data_plot; data_plot.clear(); data_plot.resize(nplots);
	vector<TH2D*> unfcov; unfcov.clear(); unfcov.resize(nplots);
	vector<TH2D*> ac_mat; ac_mat.clear(); ac_mat.resize(nplots);
	vector<TCanvas*> can; can.clear(); can.resize(nplots);
	vector<TLegend*> leg; leg.clear(); leg.resize(nplots);	

	//1st index: plot, 2nd index: sample
	vector<vector<TH1D*> > true_plot; true_plot.clear(); true_plot.resize(nplots);						

	//----------------------------------------//	

	vector<TString> sample_names; sample_names.clear();
	vector<int> colors; colors.clear();		
	vector<TString> labels; labels.clear();

	sample_names.push_back("AR23"); colors.push_back(kOrange+10); labels.push_back("AR23 ");
	sample_names.push_back("AR25"); colors.push_back(kGreen+2); labels.push_back("AR25 ");

	//----------------------------------------//

	const int nmcsamples = sample_names.size();
	vector<TFile*> fmc; fmc.clear();
	fmc.resize(nmcsamples);

	//----------------------------------------//

	// Loop over the plots

	for (int iplot = 0; iplot < nplots; iplot ++) {	

		//----------------------------------------//
		
		// export the variable name for data plots

		TString var = "delta_p_t";
		TString xlabel = "#deltap_{T} [GeV/c]";
		TString ylabel = "#frac{d#sigma}{d#deltap_{T}} [10^{-38} cm^{2}/(GeV/c)]";		

		if (plot_name.at(iplot) == "TrueCosThetaMuSumPPlot") { 
			
			var = "cos(Mu,P_sum)"; 
		
			xlabel = "cos(#theta_{#mu,P_{sum}})";
			ylabel = "#frac{d#sigma}{dcos(#theta_{#mu,P_{sum}})} [10^{-38} cm^{2}]";

		}	

		if (plot_name.at(iplot) == "TrueCosThetaPLPRPlot") { 
			
			var = "cos(P_L,P_R)"; 
			xlabel = "cos(#theta_{P_{L},P_{R}})";
			ylabel = "#frac{d#sigma}{dcos(#theta_{P_{L},P_{R}})} [10^{-38} cm^{2}]";

		}	

		//----------------------------------------//

		// grab data plots

		data_plot.at(iplot) = (TH1D*)(fdata->Get("reco_"+var+"/reco_"+var+"_DataHist"));
		unfcov.at(iplot) = (TH2D*)(fdata->Get("reco_"+var+"/reco_"+var+"_CovMat"));
		ac_mat.at(iplot) = (TH2D*)(fdata->Get("reco_"+var+"/reco_"+var+"_AC"));

		//----------------------------------------//

		// mc plots

		true_plot.at(iplot).resize(nmcsamples);

		for (int imcsample = 0; imcsample < nmcsamples; imcsample++) {

			TString fmc_name = "mc_files/analyzerOutput_" + sample_names.at(imcsample) + ".root"; 
			fmc.at(imcsample) = TFile::Open(fmc_name,"readonly"); 	

			true_plot.at(iplot).at(imcsample) = (TH1D*)(fmc.at(imcsample)->Get(plot_name.at(iplot)));			

		}

		//----------------------------------------//

		can.at(iplot) = new TCanvas(plot_name.at(iplot),plot_name.at(iplot),205,34,1024,768);
		can.at(iplot)->cd();
		can.at(iplot)->SetBottomMargin(0.16);
		can.at(iplot)->SetTopMargin(0.1);
		can.at(iplot)->SetLeftMargin(0.18);
		can.at(iplot)->SetRightMargin(0.05);	
		
		//----------------------------------------//

		leg.at(iplot) = new TLegend(0.05,0.91,0.95,0.98);
		leg.at(iplot)->SetBorderSize(0);
		leg.at(iplot)->SetTextSize(text_size);
		leg.at(iplot)->SetTextFont(text_font);
		leg.at(iplot)->SetNColumns(3);
		leg.at(iplot)->SetMargin(0.15);				

		//----------------------------------------//

		// data

		double MaxValue = data_plot.at(iplot)->GetMaximum();
		double MinValue = data_plot.at(iplot)->GetMinimum();
		data_plot.at(iplot)->SetMinimum(0.);

		data_plot.at(iplot)->GetXaxis()->SetTitle(xlabel);		
		data_plot.at(iplot)->GetXaxis()->SetTitleFont(text_font);			
		data_plot.at(iplot)->GetXaxis()->SetLabelFont(text_font);					
		data_plot.at(iplot)->GetXaxis()->SetLabelSize(text_size);				
		data_plot.at(iplot)->GetXaxis()->SetTitleSize(text_size);	
		data_plot.at(iplot)->GetXaxis()->SetTitleOffset(1.1);	
		data_plot.at(iplot)->GetXaxis()->CenterTitle();	
		data_plot.at(iplot)->GetXaxis()->SetNdivisions(8);								

		data_plot.at(iplot)->GetYaxis()->SetTitle(ylabel);		
		data_plot.at(iplot)->GetYaxis()->SetTitleFont(text_font);			
		data_plot.at(iplot)->GetYaxis()->SetLabelFont(text_font);					
		data_plot.at(iplot)->GetYaxis()->SetLabelSize(text_size);				
		data_plot.at(iplot)->GetYaxis()->SetTitleSize(text_size);			
		data_plot.at(iplot)->GetYaxis()->SetTitleOffset(1.3);
		data_plot.at(iplot)->GetYaxis()->CenterTitle();	
		data_plot.at(iplot)->GetYaxis()->SetNdivisions(8);					

		data_plot.at(iplot)->SetTitle("");		
		data_plot.at(iplot)->SetLineColor(kBlack);
		data_plot.at(iplot)->SetMarkerColor(kBlack);
		data_plot.at(iplot)->SetMarkerSize(1.);
		data_plot.at(iplot)->SetMarkerStyle(20);
		data_plot.at(iplot)->SetLineWidth(1);					

		leg.at(iplot)->AddEntry(data_plot.at(iplot),"MicroBooNE Data","ep");					
		can.at(iplot)->cd();

		data_plot.at(iplot)->Draw("e1x0 same");

		//------------------------------//		

		// arrays for nmcsamples

		double Chi2[nmcsamples];
		int Ndof[nmcsamples];
		double pval[nmcsamples];

		// Clones for the nmcsamples model predictions

		TH1D* Clone[nmcsamples];			

		for (int imcsample = 0; imcsample < nmcsamples; imcsample++) {

			// Apply the additional smearing matrix Ac
			Clone[imcsample] = Multiply(true_plot.at(iplot).at(imcsample),ac_mat.at(iplot));											
			// Divide by the bin width
			divide_bin_width(Clone[imcsample],1.);																

			Clone[imcsample]->SetLineColor(colors[imcsample]);
			Clone[imcsample]->SetMarkerColor(colors[imcsample]);

			Clone[imcsample]->SetLineWidth(3);		
			Clone[imcsample]->Draw("hist same");		

			calc_chi2(Clone[imcsample],data_plot.at(iplot),unfcov.at(iplot),Chi2[imcsample],Ndof[imcsample],pval[imcsample]);
			TString Chi2NdofAlt = "(" + to_string_with_precision(Chi2[imcsample],1) + "/" + TString(std::to_string(Ndof[imcsample])) +")";

			TLegendEntry* lGenie = leg.at(iplot)->AddEntry(Clone[imcsample],labels[imcsample] + Chi2NdofAlt,"l");
			lGenie->SetTextColor(colors[imcsample]); 
			
			//----------------------------------------//

			// only for AR23, add systematics

			if (sample_names[imcsample] == "AR23") {

				int n = Clone[imcsample]->GetNbinsX();
				TH2D* cov_ar23 = (TH2D*)ar23_unc_file->Get("tot_covariance_"+plot_name[iplot]);			

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
				// c_mat_cov_ar23.SaveAs("debug_pdf/c_mat_cov_ar23_"+PlotNames[iplot]+".pdf");

				TMatrixD trans_mat_cov_ar23(n,n);
				//trans_mat_cov_ar23.Transpose(mat_cov_ar23);
				trans_mat_cov_ar23 = mat_cov_ar23;						

				//TCanvas c_trans_mat_cov_ar23("c_trans_mat_cov_ar23", "c_trans_mat_cov_ar23", 800, 800);
				//trans_mat_cov_ar23.Draw("COLZ TEXT");   // COLZ = color map, TEXT = numbers
				//c_trans_mat_cov_ar23.SaveAs("debug_pdf/c_trans_mat_cov_ar23_"+PlotNames[iplot]+".pdf");				

				TMatrixD mat_Ac(n,n);
				H2M(ac_mat.at(iplot), mat_Ac,kFALSE);		
				
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
				//c_mat_ac_ar23.SaveAs("debug_pdf/c_mat_ac_ar23_"+PlotNames[iplot]+".pdf");		

				TMatrixD trans_mat_Ac(n,n);

				for (int i = 0; i < n; ++i) {
				
					for (int j = 0; j < n; ++j) {
					
						trans_mat_Ac(i, j) = mat_Ac(n - 1 - j, n - 1 - i);
					
					}
				
				}		

				//TCanvas c_trans_mat_ac_ar23("c_trans_mat_ac_ar23", "c_trans_mat_ac_ar23", 800, 800);
				//trans_mat_Ac.Draw("COLZ TEXT");   // COLZ = color map, TEXT = numbers
				//c_trans_mat_ac_ar23.SaveAs("debug_pdf/c_trans_mat_ac_ar23_"+PlotNames[iplot]+".pdf");						

				// Applying the smearing: Cov' = A * Cov^T * A^T
				//TMatrixD smeared_mat_cov_ar23 = mat_Ac * trans_mat_cov_ar23 * trans_mat_Ac;
				TMatrixD smeared_mat_cov_ar23 = trans_mat_cov_ar23;	// wrong !!! fix it			

				//TCanvas c_smeared_mat_cov_ar23("c_smeared_mat_cov_ar23", "c_smeared_mat_cov_ar23", 800, 800);
				//smeared_mat_cov_ar23.Draw("COLZ TEXT");   // COLZ = color map, TEXT = numbers
				//c_smeared_mat_cov_ar23.SaveAs("debug_pdf/c_smeared_mat_cov_ar23_"+PlotNames[iplot]+".pdf");						

				TH2D* cov_ar23_smeared = (TH2D*)cov_ar23->Clone();
				flip_rows(smeared_mat_cov_ar23, cov_ar23_smeared);
				//TCanvas c_smeared_th2d_cov_ar23("c_smeared_th2d_cov_ar23", "c_smeared_th2d_cov_ar23", 800, 800);
				//cov_ar23_smeared->Draw("COLZ TEXT");   // COLZ = color map, TEXT = numbers
				//c_smeared_th2d_cov_ar23.SaveAs("debug_pdf/c_smeared_th2d_cov_ar23_"+PlotNames[iplot]+".pdf");						

				divide_bin_area(cov_ar23_smeared,Clone[imcsample]);
				set_unc_from_cov(cov_ar23_smeared,Clone[imcsample]);

				//----------------------//

				// error band on AR23 prediction

				TH1D* h = Clone[imcsample];

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
				gerr->SetFillColor(colors[imcsample]);
				gerr->SetFillStyle(3004);
				gerr->SetLineWidth(0);

				// Draw histogram line
				h->SetLineColor(colors[imcsample]);
				h->SetFillStyle(0);
				h->Draw("hist same");

				// Draw symmetric error band ONLY
				gerr->Draw("E2 same");

			}

			//----------------------------------------//

		}

		data_plot.at(iplot)->Draw("e1x0 same");
		leg.at(iplot)->Draw();					

		can.at(iplot)->SaveAs("pdf/ub_cc2p_"+plot_name[iplot]+".pdf");
		delete can.at(iplot);

		//----------------------------------------//

	} // End of the loop over the plots

} // End of the program 
