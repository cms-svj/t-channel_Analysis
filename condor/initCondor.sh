source ../local.sh
source ${TCHANNEL_LCG}/setup.sh
storage_dir=$(readlink -f $PWD/../)
export TCHANNEL_BASE=${storage_dir}
