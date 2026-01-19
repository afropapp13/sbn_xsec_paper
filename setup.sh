source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
#setup root v6_28_12 -q e26:p3915:prof
setup sbndcode v10_14_02 -q e26:prof

#htgettoken -a htvaultprod.fnal.gov -i sbnd
httokensh -a htvaultprod.fnal.gov -i sbnd -- /bin/bash

export TERM=screen