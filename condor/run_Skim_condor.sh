#!/usr/bin/env bash

dataset_longname=$1
nfiles=$2
startfile=$3
analyzeFile=$4
year=$5
base_dir=`pwd`
eosStore="/eos/uscms/store/user/lpcdarkqcd/tchannel_UL/skims/"


# echo "singularity container:"
# echo $SINGULARITY_CONTAINER

echo "ls output"
ls -l

echo "unpacking tar file"
mkdir tchannel_skims
mv Skim.tar.gz tchannel_skims/.
cd tchannel_skims
tar -xzf Skim.tar.gz
ls -l

# Setup the activation script for the virtual environment
# echo "\nSetting up the activation script for the virtual environment ... "
# source init.sh ""

# echo "\nNeed an environment similar to initCondor.sh otherwise ROOT module will be not found..."
# echo "source initROOT.sh """
# source initROOT.sh ""
echo "source root"
storage_dir=$(readlink -f $PWD)
export TCHANNEL_BASE=${storage_dir}
source /cvmfs/sft.cern.ch/lcg/views/LCG_101cuda/x86_64-centos7-gcc8-opt/setup.sh


echo "ls output"
ls -l

# echo "output of uname -s : "
# uname -s

# echo "unpacking exestuff"
# cp ${base_dir}/exestuff.tar.gz .
# tar xzf exestuff.tar.gz
# mv exestuff/* .
# ls -l

echo "\n\n Attempting to run Skim executable.\n\n"
echo ${dataset_longname}
echo ${nfiles}
echo ${startfile}
python3.8 ${analyzeFile} --condor -d ${dataset_longname} -n ${nfiles} -M ${startfile} 

echo "\n\n ls output\n"
ls -l

echo ${eosStore}/${year}/${dataset_longname}/.
xrdcp -f Skim_*.root ${eosStore}/${year}/${dataset_longname}/.
rm Skim_*.root


cd ${base_dir}
rm docker_stderror

echo "\n\n ls output\n"
ls -l

