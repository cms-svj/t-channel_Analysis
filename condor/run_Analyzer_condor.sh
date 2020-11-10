#!/bin/bash

dataset=$1
nfiles=$2
startfile=$3
base_dir=`pwd`

printf "\n\n ls output\n"
ls -l

printf "\n\n unpacking CMSSW tar file.\n\n"
tar -xzf tchannel.tar.gz
source init.sh

printf "\n\n ls output\n"
ls -l

printf "\n\n output of uname -s : "
uname -s
printf "\n\n"

cp ${base_dir}/exestuff.tar.gz .
tar xzf exestuff.tar.gz
cd exestuff/

printf "\n\n Attempting to run MyAnalysis executable.\n\n"
python analyze.py --condor -d ${dataset} -N ${nfiles} -M ${startfile}

printf "\n\n ls output\n"
ls -l

mv MyAnalysis*.root ${base_dir}
cd ${base_dir}

printf "\n\n ls output\n"
ls -l

