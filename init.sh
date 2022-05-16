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
export TCHANNEL_BASE=${storage_dir}
export JUPYTER_PATH=${storage_dir}/.jupyter
export JUPYTER_RUNTIME_DIR=${storage_dir}/.local/share/jupyter/runtime
export JUPYTER_DATA_DIR=${storage_dir}/.local/share/jupyter
export IPYTHONDIR=${storage_dir}/.ipython

echo "Sourcing virtual env from $NAME ..."
source $NAME/bin/activate

# keep terminal open
if [[ -n "$SINGULARITY_CONTAINER" ]]; then
	/bin/bash
fi
