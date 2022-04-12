#!/bin/bash

NAME=coffeaenv

if [[ "$1" != "" ]]; then
    NAME=$1
fi

# vars for jupyter
storage_dir=$(readlink -f $PWD)
export TCHANNEL_BASE=${storage_dir}
export JUPYTER_PATH=${storage_dir}/.jupyter
export JUPYTER_RUNTIME_DIR=${storage_dir}/.local/share/jupyter/runtime
export JUPYTER_DATA_DIR=${storage_dir}/.local/share/jupyter
export IPYTHONDIR=${storage_dir}/.ipython

# source virtual environment
source $NAME/bin/activate

