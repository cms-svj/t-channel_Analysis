#!/bin/bash

if [ -z "$SINGULARITY_BIND" ]; then
	if [[ "$(uname -a)" == *cms*.fnal.gov* ]]; then
		export SINGULARITY_BIND="$(readlink -f $HOME),/uscms_data,$(readlink -f ${HOME}/nobackup/)"
	fi
fi

singularity run --nv --bind /cvmfs /cvmfs/unpacked.cern.ch/registry.hub.docker.com/fnallpc/fnallpc-docker:pytorch-1.9.0-cuda11.1-cudnn8-runtime-singularity
