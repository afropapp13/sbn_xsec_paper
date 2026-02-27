#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include "Math/DistFunc.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <fstream>
#include <limits>
#include <algorithm> // Required for std::max_element

using namespace std;

//---------------------------//

const int text_font = 132;
const double text_size = 0.06;
const int ndivs = 6; 

//---------------------------//

double get_one_bin_unc(TH1D* hErr) {

	double sigma2 = 0.0;

	for (int i = 1; i <= hErr->GetNbinsX(); ++i) {
		double si = hErr->GetBinContent(i);
		double width = hErr->GetBinWidth(i);   // handles irregular binning
		sigma2 += (si * width) * (si * width);
	}

	double sigma_tot = std::sqrt(sigma2);

	return sigma_tot;

}

//---------------------------//

void PrintMatrix(const TMatrixD& M) {

  for (int i = 0; i < M.GetNrows(); ++i) {

    for (int j = 0; j < M.GetNcols(); ++j) {

      std::cout << std::setw(14)
                << std::setprecision(6)
                << std::fixed
                << M(i, j);

    }

    std::cout << '\n';

  }
  
}

//---------------------------//

void set_unc_from_cov(TH2D* cov, TH1D* h) {

	if (!cov || !h) return;

	int nbins = h->GetNbinsX();

	for (int i = 1; i <= nbins; ++i) {
		double variance = cov->GetBinContent(i, i);
		double error = (variance >= 0) ? sqrt(variance) : 0.0;
		h->SetBinError(i, error);
	}

}

//---------------------------//

TH1D* get_frac_unc_from_cov(TH2D* cov, TH1D* h) {

	TH1D* clone = (TH1D*)h->Clone();

	if (!cov || !h) return clone;

	int nbins = h->GetNbinsX();

	for (int i = 1; i <= nbins; ++i) {

		double value = h->GetBinContent(i);
		double variance = cov->GetBinContent(i, i);
		double error = (variance >= 0) ? sqrt(variance) : 0.0;
		double frac_unc = (value != 0) ? 100 * error / value : 0.0;

		clone->SetBinContent(i, frac_unc);
		clone->SetBinError(i, 0.);		

	}

	return clone;

}

//---------------------------//

void divide_bin_area(TH2D* h, TH1D* pred) {

    if (!h) return;

    TAxis* xax = pred->GetXaxis();
    TAxis* yax = pred->GetXaxis();

    int nx = h->GetNbinsX();
    int ny = h->GetNbinsY();

    for (int ix = 1; ix <= nx; ++ix) {
        double dx = xax->GetBinWidth(ix);
        for (int iy = 1; iy <= ny; ++iy) {
            double dy = yax->GetBinWidth(iy);
            double area = dx * dy;

            double val = h->GetBinContent(ix, iy);
            double err = h->GetBinError(ix, iy);

            if (area > 0) {
                h->SetBinContent(ix, iy, val / area);
                h->SetBinError(ix, iy, err / area);
            }
        }
    }
}

//---------------------------//

TString int_to_string(int num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;
  
}

//---------------------------//

TString double_to_string(double num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}

//----------------------------------------//

TString to_string_with_precision(double a_value, const int n = 1) {
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return TString(out.str());

}

//----------------------------------------//	

TString ToString(int num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}

//----------------------------------------//	

TString ToString(double num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}

//----------------------------------------//		

void divide_bin_width(TH1D* h, double SF=1.) {

	int NBins = h->GetXaxis()->GetNbins();
  
	for (int i = 0; i < NBins; i++) {
  
	  double CurrentEntry = h->GetBinContent(i+1);
	  double NewEntry = SF * CurrentEntry / h->GetBinWidth(i+1);
  
	  double CurrentError = h->GetBinError(i+1);
	  double NewError = SF * CurrentError / h->GetBinWidth(i+1);
  
	  h->SetBinContent(i+1,NewEntry); 
	  h->SetBinError(i+1,NewError); 
	  //h->SetBinError(i+1,0.000001); 
  
	}
  
  }

  //----------------------------------------//		

void rm_bin_width(TH1D* h, double SF=1.) {

	int NBins = h->GetXaxis()->GetNbins();
  
	for (int i = 0; i < NBins; i++) {
  
	  double CurrentEntry = h->GetBinContent(i+1);
	  double NewEntry = SF * CurrentEntry * h->GetBinWidth(i+1);
  
	  double CurrentError = h->GetBinError(i+1);
	  double NewError = SF * CurrentError * h->GetBinWidth(i+1);
  
	  h->SetBinContent(i+1,NewEntry); 
	  h->SetBinError(i+1,NewError); 
	  //h->SetBinError(i+1,0.000001); 
  
	}
  
  }

   //----------------------------------------//	 

void H2V(const TH1D* histo, TVectorD& vec)
{
    // Fill 1D histogram into matrix
    for(Int_t i=0; i<histo->GetNbinsX(); i++)
    {
        vec(i) = histo->GetBinContent(i+1);
    }
}

  //----------------------------------------//	

void flip_rows(const TMatrixD mat, TH2D* histo)
{

	const Int_t nrows = mat.GetNrows();
	const Int_t ncols = mat.GetNcols();

	for (Int_t i = 0; i < nrows; ++i)
	{
		for (Int_t j = 0; j < ncols; ++j)
		{
			histo->SetBinContent(
                i + 1,        // keep rows the same
                ncols - j,    // flip column index
                mat(i, j)
			);
		}
	}

}

//----------------------------------------//	

void M2H(const TMatrixD mat, TH2D* histo)
{

	// Fill matrix to histogram
	for(Int_t i=0; i<mat.GetNrows(); i++)
	{
		for(Int_t j=0; j<mat.GetNcols(); j++)
		{
			histo->SetBinContent(i+1, j+1, mat(i, j));
		}
	}

}

//----------------------------------------//

void V2H(const TVectorD vec, TH1D* histo)
{
    // Fill vector to histogram,
    for(Int_t i=0; i<vec.GetNrows(); i++)
    {
        histo->SetBinContent(i+1, vec(i));
    }
}

  //----------------------------------------//	

void H2M(const TH2D* histo, TMatrixD& mat, bool rowcolumn) {

    // Fill 2D histogram into matrix
    // If TH2D(i, j) = Matrix(i, j), rowcolumn = kTRUE, else rowcolumn = kFALSE

    for (Int_t i=0; i<histo->GetNbinsX(); i++) {

        for (Int_t j=0; j<histo->GetNbinsY(); j++) {

            if (rowcolumn) { mat(i, j) = histo->GetBinContent(i+1, j+1); }
            else { mat(j, i) = histo->GetBinContent(i+1, j+1); }

        }

    }

}

  //----------------------------------------//	

  TH1D* Multiply(TH1D* True, TH2D* SmearMatrix) {

	TH1D* TrueClone = (TH1D*)(True->Clone());

	int XBins = SmearMatrix->GetXaxis()->GetNbins();
	int YBins = SmearMatrix->GetYaxis()->GetNbins();

	if (XBins != YBins) { std::cout << "Not symmetric matrix" << std::endl; }

	TVectorD signal(XBins);
	TMatrixD response(XBins,XBins);

	H2V(TrueClone, signal);
	H2M(SmearMatrix, response, kTRUE);

	TVectorD RecoSpace = response * signal;
	V2H(RecoSpace, TrueClone);	

	return TrueClone;

}

//----------------------------------------//	

double Chi2Prob(double chi2, int ndof)
{
    return ROOT::Math::chisquared_cdf_c(chi2, ndof);
}

//----------------------------------------//	

void calc_chi2(TH1D* h_model, TH1D* h_data, TH2D* cov, double &chi, int &ndof, double &pval) {

	// Clone them so we can scale them 

	TH1D* h_model_clone = (TH1D*)h_model->Clone();
	TH1D* h_data_clone  = (TH1D*)h_data->Clone();
	TH2D* h_cov_clone   = (TH2D*)cov->Clone();
	int NBins = h_cov_clone->GetNbinsX();

	// Getting covariance matrix in TMatrix form

	TMatrixD cov_m;
	cov_m.Clear();
	cov_m.ResizeTo(NBins,NBins);

	// loop over rows

	for (int i = 0; i < NBins; i++) {			

		// loop over columns

		for (int j = 0; j < NBins; j++) {

			cov_m[i][j] = h_cov_clone->GetBinContent(i+1, j+1);
 
		}
	
	}

	TMatrixD copy_cov_m = cov_m;

	// Inverting the covariance matrix
	TMatrixD inverse_cov_m = cov_m.Invert();

	// Calculating the chi2 = Summation_ij{ (x_i - mu_j)*E_ij^(-1)*(x_j - mu_j)  }
	// x = data, mu = model, E^(-1) = inverted covariance matrix 

	chi = 0.;
	
	for (int i = 0; i < NBins; i++) {

		//double XWidth = h_data_clone->GetBinWidth(i+1);

		for (int j = 0; j < NBins; j++) {

			//double YWidth = h_data_clone->GetBinWidth(i+1);

			double diffi = h_data_clone->GetBinContent(i+1) - h_model_clone->GetBinContent(i+1);
			double diffj = h_data_clone->GetBinContent(j+1) - h_model_clone->GetBinContent(j+1);
			double LocalChi = diffi * inverse_cov_m[i][j] * diffj; 
			chi += LocalChi;

		}

	}

	ndof = h_data_clone->GetNbinsX();
	pval = Chi2Prob(chi, ndof);

	delete h_model_clone;
	delete h_data_clone;
	delete h_cov_clone;

}

//----------------------------------------//	