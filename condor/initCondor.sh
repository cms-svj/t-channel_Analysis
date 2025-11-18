# source ../local.sh
# source ${TCHANNEL_LCG}/setup.sh
# source /cvmfs/sft.cern.ch/lcg/views/LCG_106_cuda/x86_64-el9-gcc11-opt/setup.sh
source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-el9-gcc12-opt/setup.sh
storage_dir=$(readlink -f $PWD/../)
export TCHANNEL_BASE=${storage_dir}
