{
	vector<TString> WhichSample; vector<TString> WhichName;

	//---------------------------//

	vector< tuple<TString,int> > tweaks = {

		// ---------------- ZExp PCA (7 universes) ----------------
		/*{"ZExpPCAWeighter_SBNNuSyst_multisigma_D_ZExp_b1", 7},
		{"ZExpPCAWeighter_SBNNuSyst_multisigma_D_ZExp_b2", 7},
		{"ZExpPCAWeighter_SBNNuSyst_multisigma_D_ZExp_b3", 7},
		{"ZExpPCAWeighter_SBNNuSyst_multisigma_D_ZExp_b4", 7},*/

		{"ZExpPCAWeighter_SBNNuSyst_multisigma_MvA_ZExp_b1", 7},
		{"ZExpPCAWeighter_SBNNuSyst_multisigma_MvA_ZExp_b2", 7},
		{"ZExpPCAWeighter_SBNNuSyst_multisigma_MvA_ZExp_b3", 7},
		{"ZExpPCAWeighter_SBNNuSyst_multisigma_MvA_ZExp_b4", 7},

		// ---------------- CCQE Template (SF) ----------------
		{"CCQETemplateReweight_SBNNuSyst_multisigma_SF_q0bin1", 6},
		{"CCQETemplateReweight_SBNNuSyst_multisigma_SF_q0bin2", 6},
		{"CCQETemplateReweight_SBNNuSyst_multisigma_SF_q0bin3", 6},
		{"CCQETemplateReweight_SBNNuSyst_multisigma_SF_q0bin4", 6},
		{"CCQETemplateReweight_SBNNuSyst_multisigma_SF_q0bin5", 6},

		// ---------------- CCQE Template (HF) ----------------
		{"CCQETemplateReweight_SBNNuSyst_multisigma_HF_q0bin1", 6},
		{"CCQETemplateReweight_SBNNuSyst_multisigma_HF_q0bin2", 6},
		{"CCQETemplateReweight_SBNNuSyst_multisigma_HF_q0bin3", 6},
		{"CCQETemplateReweight_SBNNuSyst_multisigma_HF_q0bin4", 6},
		{"CCQETemplateReweight_SBNNuSyst_multisigma_HF_q0bin5", 6},

		// ---------------- CCQE Template (CRPA) ----------------
		{"CCQETemplateReweight_SBNNuSyst_multisigma_CRPA_q0bin1", 6},
		{"CCQETemplateReweight_SBNNuSyst_multisigma_CRPA_q0bin2", 6},
		{"CCQETemplateReweight_SBNNuSyst_multisigma_CRPA_q0bin3", 6},
		{"CCQETemplateReweight_SBNNuSyst_multisigma_CRPA_q0bin4", 6},
		{"CCQETemplateReweight_SBNNuSyst_multisigma_CRPA_q0bin5", 6},

		// ---------------- QE Interference ----------------
		{"QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_0", 6},
		{"QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_1", 6},
		{"QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_2", 6},
		{"QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_3", 6},
		{"QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_4", 6},
		{"QEInterference_SBNNuSyst_multisigma_INT_QEIntf_dial_5", 6},

		// ---------------- GENIE EDepFSI knobs ----------------
		{"GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_VecFFCCQEshape", 6},
		{"GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_CoulombCCQE", 6},
		{"GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_NormCCMEC", 6},
		{"GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_NormNCMEC", 6},
		{"GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_DecayAngMEC", 6},
		{"GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_MFP_pi", 6},
		{"GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrCEx_pi", 6},
		{"GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrInel_pi", 6},
		{"GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrAbs_pi", 6},
		{"GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrPiProd_pi", 6},
		{"GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrG4_N", 6},
		{"GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrINCL_N", 6},
		{"GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrG4LoE_N", 6},
		{"GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrINCLLoE_N", 6},
		{"GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrG4M1E_N", 6},
		{"GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrINCLM1E_N", 6},
		{"GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrG4M2E_N", 6},
		{"GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrINCLM2E_N", 6},
		{"GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrG4HiE_N", 6},
		{"GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrINCLHiE_N", 6},
		{"GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_MFPLoE_N", 6},
		{"GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_MFPM1E_N", 6},
		{"GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_MFPM2E_N", 6},
		{"GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_MFPHiE_N", 6},
		{"GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrKin_PiProFix_N", 6},
		{"GENIEReWeight_SBNNuSyst_multisigma_EDepFSI_FrKin_PiProBias_N", 6},

		// ---------------- MEC Valencia ----------------
		{"MECq0q3InterpWeighting_SuSAv2ToValenica_q0binned_MECResponse_q0bin0", 6},
		{"MECq0q3InterpWeighting_SuSAv2ToValenica_q0binned_MECResponse_q0bin1", 6},
		{"MECq0q3InterpWeighting_SuSAv2ToValenica_q0binned_MECResponse_q0bin2", 6},
		{"MECq0q3InterpWeighting_SuSAv2ToValenica_q0binned_MECResponse_q0bin3", 6},

		// ---------------- MEC Martini ----------------
		{"MECq0q3InterpWeighting_SuSAv2ToMartini_q0binned_MECResponse_q0bin0", 6},
		{"MECq0q3InterpWeighting_SuSAv2ToMartini_q0binned_MECResponse_q0bin1", 6},
		{"MECq0q3InterpWeighting_SuSAv2ToMartini_q0binned_MECResponse_q0bin2", 6},
		{"MECq0q3InterpWeighting_SuSAv2ToMartini_q0binned_MECResponse_q0bin3", 6}
	};

	//----------------------------------------//

	WhichSample.push_back("/pnfs/sbnd/persistent/users/apapadop/GENIETweakedSamples/v3_6_2_AR23_20i_00_000/14_1000180400_CC_v3_6_2_AR23_20i_00_000.flat.root"); WhichName.push_back("AR23");
	WhichSample.push_back("/pnfs/sbnd/persistent/users/apapadop/GENIETweakedSamples/v3_6_2_AR25_20i_00_000_sbn1/14_1000180400_CC_v3_6_2_AR25_20i_00_000.flat.root"); WhichName.push_back("AR25");

	//WhichSample.push_back("/pnfs/sbnd/persistent/users/apapadop/GENIETweakedSamples/v3_6_2_AR23_32i_00_000_liang_bugfix_microboone/14_1000180400_CC_v3_6_2_AR23_32i_00_000.flat.root"); WhichName.push_back("AR23_LQCD");
	//WhichSample.push_back("/pnfs/sbnd/persistent/users/apapadop/GENIETweakedSamples/v3_6_2_AR23_31i_00_000_liang_bugfix_microboone/14_1000180400_CC_v3_6_2_AR23_31i_00_000.flat.root"); WhichName.push_back("AR23_Minerva");		

	//----------------------------------------//

	gROOT->ProcessLine(".L analyzer.cxx++");

	for (int i =0;i < (int)(WhichSample.size()); i++) {

		gROOT->ProcessLine("analyzer(\""+WhichSample[i]+"\",\""+WhichName[i]+"\").Loop()");

		if (WhichName[i] == "AR23") {

			int n_tweaks = (int)(tweaks.size());

			for (int itweak=0; itweak < n_tweaks; itweak++) {

				int n_indices= std::get<1>(tweaks[itweak]);

				for (int i_index=0; i_index < n_indices; i_index++) {

					gROOT->ProcessLine("analyzer(\""+WhichSample[i]+"\",\""+WhichName[i]+"\",\""+std::get<0>(tweaks[itweak])+"\","+ i_index +").Loop()");

				}

			}

		}

	}

	//gROOT->ProcessLine(".q");

};
