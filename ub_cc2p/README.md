# run over the mc files
cd flat_analyzer
root -b loop_analyzer.cxx

#go back to the main cc2p firectory
cd ..

# make the covariances
root -b make_covariance.cxx

# overlay data and mc
root -b overlay_data_mc_xsec.cxx