#!/bin/bash

source local.sh

if [ -z "$SINGULARITY_BIND" ]; then
	if [[ "$(uname -a)" == *cms*.fnal.gov* ]]; then
		export SINGULARITY_BIND="$(readlink -f $HOME),/uscms_data,$(readlink -f ${HOME}/nobackup/)"
	fi
fi

singularity run --nv --bind /cvmfs $TCHANNEL_SC $@
