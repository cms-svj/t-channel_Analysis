#!/bin/bash

source local.sh
NAME=$1
if [ -z "$NAME" ]; then
	NAME=coffeaenv
fi

if [[ "$BASH_SOURCE" == "initLCG.sh" ]]; then
	NAME=${NAME}LCG
elif [[ -z "$SINGULARITY_CONTAINER" ]]; then
	./launchSingularity.sh "$BASH_SOURCE"
	return
fi

# vars for jupyter
storage_dir=$(readlink -f $PWD)
export TCHANNEL_BASE=/srv
export JUPYTER_PATH=/srv/.jupyter
export JUPYTER_RUNTIME_DIR=/srv/.local/share/jupyter/runtime
export JUPYTER_DATA_DIR=/srv/.local/share/jupyter
export IPYTHONDIR=/srv/.ipython
export CONDOR_CONFIG=/srv/.condor_config
export COFFEA_IMAGE=$TCHANNEL_SC
export COFFEA_IMAGE_FULL=/cvmfs/unpacked.cern.ch/registry.hub.docker.com/${TCHANNEL_SC}
echo "Sourcing virtual env from $NAME ..."
source $NAME/bin/activate

# keep terminal open
if [[ -n "$SINGULARITY_CONTAINER" ]]; then
	/bin/bash
	
fi
