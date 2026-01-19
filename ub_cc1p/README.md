# run the event selection on the mc files
cd flat_analyzer
root -b loop_analyzer.cxx

# run the covariances
cd ..
root -b make_covariance.cxx

# uncertainty breakdown
root -b plot_breakdown_unc.cxx

# overlay data/mc
root -b overlay_data_mc_xsec.cxx