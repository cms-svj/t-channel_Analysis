#!/bin/bash

dataset_longname=$1
nfiles=$2
startfile=$3
base_dir=`pwd`

echo "ls output"
ls -l

echo "unpacking tar file"
mkdir tchannel
mv tchannel.tar.gz tchannel/.
cd tchannel
tar -xzf tchannel.tar.gz
source init.sh

echo "ls output"
ls -l

echo "output of uname -s : "
uname -s

#printf "\n\n"
#cp ${base_dir}/exestuff.tar.gz .
#tar xzf exestuff.tar.gz

echo "\n\n Attempting to run MyAnalysis executable.\n\n"
echo ${dataset_longname}
echo ${nfiles}
echo ${startfile}
python analyze.py --condor -d ${dataset_longname} -N ${nfiles} -M ${startfile}

echo "\n\n ls output\n"
ls -l

mv MyAnalysis*.root ${base_dir}
cd ${base_dir}
rm docker_stderror

echo "\n\n ls output\n"
ls -l

