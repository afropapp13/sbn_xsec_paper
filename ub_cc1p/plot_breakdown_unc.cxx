#include <TFile.h>
#include <TF1.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TMatrixD.h>
#include <TVectorD.h>

#include "../helper_functions.h"

using namespace std;

//--------------------//

void plot_breakdown_unc() {

	//--------------------//

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();
	TGaxis::SetMaxDigits(3);
	gStyle->SetOptStat(0);

	//--------------------//

	vector<TString> plot_names;
	plot_names.push_back("TrueSingleBinPlot");	
	plot_names.push_back("TrueMuonCosThetaPlot");	
	plot_names.push_back("TrueDeltaPTPlot"); 
	plot_names.push_back("TrueDeltaAlphaTPlot"); 
	//plot_names.push_back("SerialTrueDeltaPT_DeltaAlphaTPlot");
	//plot_names.push_back("SerialTrueDeltaAlphaT_DeltaPTPlot");

	const int nplots = plot_names.size();

	//--------------------//

	std::vector<TString> knob;

	// knob.push_back("b1");
	knob.push_back("b2");
	knob.push_back("b3");
	knob.push_back("b4");
	// knob.push_back("q0bin1");
	// knob.push_back("q0bin2");
	knob.push_back("q0bin3");
	// ////knob.push_back("q0bin4");
	// ////knob.push_back("q0bin5");
	knob.push_back("QEIntf_dial_0");
	// knob.push_back("QEIntf_dial_1");
	// knob.push_back("QEIntf_dial_2");
	// knob.push_back("QEIntf_dial_3");
	// knob.push_back("QEIntf_dial_4");
	// knob.push_back("QEIntf_dial_5");
	// knob.push_back("VecFFCCQEshape");
	// knob.push_back("CoulombCCQE");
	knob.push_back("NormCCMEC");
	// knob.push_back("NormNCMEC");
	knob.push_back("DecayAngMEC");
	// knob.push_back("MFP_pi");
	// knob.push_back("FrCEx_pi");
	// knob.push_back("FrInel_pi");
	// knob.push_back("FrAbs_pi");
	// knob.push_back("FrPiProd_pi"); 
	// knob.push_back("FrG4_N");
	// knob.push_back("FrINCL_N");
	// knob.push_back("FrG4LoE_N");
	// knob.push_back("FrINCLLoE_N");
	knob.push_back("FrG4M1E_N");
	// knob.push_back("FrINCLM1E_N");
	// knob.push_back("FrG4M2E_N");
	knob.push_back("FrINCLM2E_N");
	// knob.push_back("FrG4HiE_N");
	// knob.push_back("FrINCLHiE_N");
	// knob.push_back("MFPLoE_N");
	knob.push_back("MFPM1E_N");
	// knob.push_back("MFPM2E_N");
	// knob.push_back("MFPHiE_N");
	// knob.push_back("FrKin_PiProFix_N");
	// knob.push_back("FrKin_PiProBias_N");
	knob.push_back("MECResponse");

	int nknobs = knob.size();

	//--------------------//

	// file containing the covariances

	TFile* cov_file = TFile::Open("mc_files/covariances.root","readonly");

	//--------------------//
	
	// declaration of matrices & canvases

	//1st index = plot, 2nd index = knob
	vector<TH2D*> tot_covariance; tot_covariance.resize(nplots);	
	vector<TH1D*> tot_unc_plot; tot_unc_plot.resize(nplots);
	vector<TH1D*> cv_plot; cv_plot.resize(nplots);			
	vector< vector<TH2D*> > knob_cov; knob_cov.resize(nplots);
	vector< vector<TH1D*> > knob_unc_plot; knob_unc_plot.resize(nplots);	

	//1st index = plot
	vector< TCanvas*> canvas; canvas.resize(nplots);	
	vector< TLegend*> leg; leg.resize(nplots);			

	//--------------------//

	// color palette

	std::vector<Color_t> colors = {
		kRed,
		kBlue+2,     // bluish
		kOrange+7,   // orange variants
		kGreen+2,    // greenish
		kRed+2,      // reddish
		kAzure+7,    // azure variants
		kMagenta+1,  // magenta variants
		kGray+1,     // grayish
		kCyan+2,     // cyanish
		kYellow+3,   // yellow variants
		kViolet+4    // violet variants
	};

	//--------------------//	

	// looping over the plots

	for (int iplot = 0; iplot < nplots; iplot++) {

		// create canvas
		canvas.at(iplot) = new TCanvas("c_break_unc_"+plot_names[iplot],"c_break_unc_"+plot_names[iplot],205,34,1024,768);
		canvas.at(iplot)->cd();
		canvas.at(iplot)->SetBottomMargin(0.15);
		canvas.at(iplot)->SetLeftMargin(0.15);
		canvas.at(iplot)->SetRightMargin(0.1);	
		canvas.at(iplot)->SetTopMargin(0.17);

		// legend
		leg.at(iplot) = new TLegend(0.0,0.84,0.9,0.98);
		leg.at(iplot)->SetBorderSize(0);
		leg.at(iplot)->SetTextFont(text_font);
		leg.at(iplot)->SetTextSize(text_size-0.01);
		leg.at(iplot)->SetNColumns(4);	

		// 2D vectors for the covariances
		knob_cov.at(iplot).resize(nknobs);	
		knob_unc_plot.at(iplot).resize(nknobs);			
		
		tot_covariance.at(iplot) = (TH2D*)(cov_file->Get("tot_covariance_"+plot_names[iplot]));	
		cv_plot.at(iplot) = (TH1D*)(cov_file->Get("cv_plot_"+plot_names[iplot]));			
		tot_unc_plot.at(iplot) = get_frac_unc_from_cov(tot_covariance.at(iplot),cv_plot.at(iplot));
		
		canvas.at(iplot)->cd();

		tot_unc_plot.at(iplot)->GetXaxis()->CenterTitle();
		tot_unc_plot.at(iplot)->GetXaxis()->SetTitleFont(text_font);
		tot_unc_plot.at(iplot)->GetXaxis()->SetTitleSize(text_size);
		tot_unc_plot.at(iplot)->GetXaxis()->SetLabelFont(text_font);
		tot_unc_plot.at(iplot)->GetXaxis()->SetLabelSize(text_size);
		tot_unc_plot.at(iplot)->GetXaxis()->SetNdivisions(6);	
		
		tot_unc_plot.at(iplot)->GetXaxis()->CenterTitle();
		tot_unc_plot.at(iplot)->GetXaxis()->SetTitleFont(text_font);
		tot_unc_plot.at(iplot)->GetXaxis()->SetTitleSize(text_size);
		tot_unc_plot.at(iplot)->GetXaxis()->SetLabelFont(text_font);
		tot_unc_plot.at(iplot)->GetXaxis()->SetLabelSize(text_size);
		tot_unc_plot.at(iplot)->GetXaxis()->SetNdivisions(6);

		tot_unc_plot.at(iplot)->GetYaxis()->CenterTitle();
		tot_unc_plot.at(iplot)->GetYaxis()->SetTitleFont(text_font);
		tot_unc_plot.at(iplot)->GetYaxis()->SetTitleSize(text_size);
		tot_unc_plot.at(iplot)->GetYaxis()->SetLabelFont(text_font);
		tot_unc_plot.at(iplot)->GetYaxis()->SetLabelSize(text_size);
		tot_unc_plot.at(iplot)->GetYaxis()->SetNdivisions(6);
		tot_unc_plot.at(iplot)->GetYaxis()->SetTitle("fractional uncertainty (%)");		
			
		tot_unc_plot.at(iplot)->SetMaximum( 59.);
		tot_unc_plot.at(iplot)->SetMinimum( 0 );	
		
		if (plot_names[iplot] == "TrueSingleBinPlot") {

			tot_unc_plot.at(iplot)->SetMaximum( 29.);

		}

		tot_unc_plot.at(iplot)->SetLineWidth(3);			
		tot_unc_plot.at(iplot)->SetLineColor(kBlack);
		tot_unc_plot.at(iplot)->Draw("hist");
		TString total_unc_str = "total (" + to_string_with_precision(get_one_bin_unc(tot_unc_plot.at(iplot)), 1) + "%)";		
		leg.at(iplot)->AddEntry(tot_unc_plot.at(iplot), "total", "");		
		//cout << total_unc_str << endl;

		tot_unc_plot.at(iplot)->Draw();

		// looping over the knobs

		for (int iknob = 0; iknob < nknobs; iknob++) {

			knob_cov.at(iplot).at(iknob) = (TH2D*)(cov_file->Get("covariance_"+knob[iknob]+"_"+plot_names[iplot]));	
			knob_unc_plot.at(iplot).at(iknob) = get_frac_unc_from_cov(knob_cov.at(iplot).at(iknob),cv_plot.at(iplot));	
			
			knob_unc_plot.at(iplot).at(iknob)->SetLineColor(colors.at(iknob)); // add me
			knob_unc_plot.at(iplot).at(iknob)->Draw("same");
			//leg.at(iplot)->AddEntry(knob_unc_plot.at(iplot).at(iknob), knob[iknob], "l");
			
			TLegendEntry* lGenie = leg.at(iplot)->AddEntry(knob_unc_plot.at(iplot).at(iknob), knob[iknob],"");
			lGenie->SetTextColor(colors.at(iknob)); // add me		
			TString unc_str = knob[iknob] + " (" + to_string_with_precision(get_one_bin_unc(knob_unc_plot.at(iplot).at(iknob)), 1) + "%)";		
			//cout << unc_str << endl;							

		} // end of the loop over the knobs

		leg.at(iplot)->Draw("same");
		canvas.at(iplot)->SaveAs("pdf/frav_unc_"+plot_names[iplot]+".pdf");				

	} // end of the loop over the plots

	//--------------------//	

	cov_file->Close();	

} // end of the program
