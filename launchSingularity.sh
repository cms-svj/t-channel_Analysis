#!/bin/bash

source local.sh

grep -v '^include' /etc/condor/config.d/01_cmslpc_interactive > .condor_config

if [ -z "$SINGULARITY_BIND" ]; then
	if [[ "$(uname -a)" == *cms*.fnal.gov* ]]; then
		export SINGULARITY_BIND="$(readlink -f $HOME),/uscms_data,$(readlink -f ${HOME}/nobackup/)"
                export APPTAINER_BIND="$(readlink -f $HOME),/uscms_data,$(readlink -f ${HOME}/nobackup/)"
	fi
fi

singularity exec -B ${PWD}:/srv -B /cvmfs -B /uscmst1b_scratch --pwd /srv \
    /cvmfs/unpacked.cern.ch/registry.hub.docker.com/${TCHANNEL_SC} \
    /bin/bash $@
