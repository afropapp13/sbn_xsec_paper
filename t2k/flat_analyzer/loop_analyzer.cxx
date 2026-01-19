{
	vector<TString> WhichSample; vector<TString> WhichName;

	//---------------------------//

	vector< tuple<TString,int> > tweaks = {
		//{"b1", 1}
		{"b1", 7},
		{"b2", 7},
		{"b3", 7},
		{"b4", 7},
		{"q0bin1", 6},
		{"q0bin2", 6},
		{"q0bin3", 6},
		//{"q0bin4", 6},
		//{"q0bin5", 6},
		{"QEIntf_dial_0", 6},
		{"QEIntf_dial_1", 6},
		{"QEIntf_dial_2", 6},
		{"QEIntf_dial_3", 6},
		{"QEIntf_dial_4", 6},
		{"QEIntf_dial_5", 6},
		{"VecFFCCQEshape", 6},
		{"CoulombCCQE", 6},
		{"NormCCMEC", 6},
		{"NormNCMEC", 6},
		{"DecayAngMEC", 6},
		{"MFP_pi", 6},
		{"FrCEx_pi", 6},
		{"FrInel_pi", 6},
		{"FrAbs_pi", 6},
		{"FrPiProd_pi", 6},
		{"FrG4_N", 6},
		{"FrINCL_N", 6},
		{"FrG4LoE_N", 6},
		{"FrINCLLoE_N", 6},
		{"FrG4M1E_N", 6},
		{"FrINCLM1E_N", 6},
		{"FrG4M2E_N", 6},
		{"FrINCLM2E_N", 6},
		{"FrG4HiE_N", 6},
		{"FrINCLHiE_N", 6},
		{"MFPLoE_N", 6},
		{"MFPM1E_N", 6},
		{"MFPM2E_N", 6},
		{"MFPHiE_N", 6},
		{"FrKin_PiProFix_N", 6},
		{"FrKin_PiProBias_N", 6},
		{"MECResponse", 6}

	};

	//----------------------------------------//

	WhichSample.push_back("/pnfs/sbnd/persistent/users/apapadop/GENIETweakedSamples/t2k_v3_6_2_AR23_20i_00_000_sbn1/14_1000060120_CC_v3_6_2_AR23_20i_00_000.flat.root"); WhichName.push_back("AR23");
	WhichSample.push_back("/pnfs/sbnd/persistent/users/apapadop/GENIETweakedSamples/t2k_v3_6_2_AR25_20i_00_000_sbn1/14_1000060120_CC_v3_6_2_AR25_20i_00_000.flat.root"); WhichName.push_back("AR25");

	//----------------------------------------//

	gROOT->ProcessLine(".L FlatTreeAnalyzer.cxx++");

	for (int i =0;i < (int)(WhichSample.size()); i++) {

		gROOT->ProcessLine("FlatTreeAnalyzer(\""+WhichSample[i]+"\",\""+WhichName[i]+"\").Loop()");

		if (WhichName[i] == "AR23") {

			int n_tweaks = (int)(tweaks.size());

			for (int itweak=0; itweak < n_tweaks; itweak++) {

				int n_indices= std::get<1>(tweaks[itweak]);

				for (int i_index=0; i_index < n_indices; i_index++) {

					gROOT->ProcessLine("FlatTreeAnalyzer(\""+WhichSample[i]+"\",\""+WhichName[i]+"\",\""+std::get<0>(tweaks[itweak])+"\","+ i_index +").Loop()");

				}

			}

		}

	}

	//gROOT->ProcessLine(".q");

};
