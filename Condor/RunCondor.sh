#!/bin/bash

mkdir kAnalysis
cd kAnalysis

####
OUTDIR="root://cmseos.fnal.gov//store/user/keanet/tchannel/SVJP_08272020"

# transfering macro files
xrdcp -s $OUTDIR/tchannel_coffea.tgz .
tar -xf tchannel_coffea.tgz # untarring the tar ball
cd tchannel_coffea

source init.sh

python test.py ${1}

echo "xrdcp -f ${1}.root ${OUTDIR}/Histos/${1}.root"
xrdcp -f ${1}.root ${OUTDIR}/Histos/${1}.root

# Remove everything
pwd
cd ../..
rm -rf kAnalysis
