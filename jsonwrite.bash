#!/bin/bash
year=$1

# python writeJSON.py -p ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8 -y $year -k ST_s-channel_hadronicDecays
# python writeJSON.py -p ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8 -y $year -k ST_s-channel_leptonDecays
# python writeJSON.py -p ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8 -y $year -k ST_t-channel_antitop_Incl
# python writeJSON.py -p ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8 -y $year -k ST_tW_top_Incl
# python writeJSON.py -p ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8 -y $year -k ST_tW_antitop_Incl 
# python writeJSON.py -p ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8 -y $year -k ST_tW_top_Incl 
# python writeJSON.py -p tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-amcatnlo-pythia8 -y $year -k ST_tZq_ll_ckm

python writeJSON.py -p METData --dest input/sampleJSONs/data/ -y $year -k METData --isSkim