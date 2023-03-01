#!/bin/bash

source local.sh

if [ -z "$SINGULARITY_BIND" ]; then
	if [[ "$(uname -a)" == *cms*.fnal.gov* ]]; then
		export SINGULARITY_BIND="$(readlink -f $HOME),/uscms_data,$(readlink -f ${HOME}/nobackup/)"
	fi
fi

singularity exec -p -B ${PWD}:/srv -B /uscmst1b_scratch --pwd /srv \
    /cvmfs/unpacked.cern.ch/registry.hub.docker.com/${TCHANNEL_SC} \
    /bin/bash $@