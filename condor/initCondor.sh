source /cvmfs/sft.cern.ch/lcg/views/LCG_101cuda/x86_64-centos7-gcc8-opt/setup.sh
storage_dir=$(readlink -f $PWD/../)
export TCHANNEL_BASE=${storage_dir}
