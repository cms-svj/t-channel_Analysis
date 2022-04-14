#!/usr/bin/env bash

dataset_longname=$1
nfiles=$2
startfile=$3
workers=$4
chunksize=$5
analyzeFile=$6
NNTrainingOut=$7
base_dir=`pwd`

echo "singularity container:"
echo $SINGULARITY_CONTAINER

echo "ls output"
ls -l

echo "unpacking tar file"
mkdir tchannel
mv tchannel.tar.gz tchannel/.
cd tchannel
tar -xzf tchannel.tar.gz
ls -l

# Setup the activation script for the virtual environment
echo "\nSetting up the activation script for the virtual environment ... "
source init.sh

echo "ls output"
ls -l

echo "output of uname -s : "
uname -s

echo "unpacking exestuff"
cp ${base_dir}/exestuff.tar.gz .
tar xzf exestuff.tar.gz
mv exestuff/* .
ls -l

echo "\n\n Attempting to run MyAnalysis executable.\n\n"
echo ${dataset_longname}
echo ${nfiles}
echo ${startfile}
echo ${workers}
python ${analyzeFile} --condor -d ${dataset_longname} -N ${nfiles} -M ${startfile} -w ${workers} -s ${chunksize}

echo "\n\n ls output\n"
ls -l

if [[ ${analyzeFile} == analyze.py ]]
then
  mv MyAnalysis*.root ${base_dir}
else
  xrdcp -f MyAnalysis*.root ${NNTrainingOut}/.
  rm MyAnalysis*.root
fi

cd ${base_dir}
rm docker_stderror

echo "\n\n ls output\n"
ls -l
