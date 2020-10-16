#!/bin/bash 

cfgDir=${CMSSW_BASE}/src/t-channel_Analysis/EventLoopFramework/cfg

#Default options
sampleSCfgFile="sampleSets_v6.cfg"
#sampleSCfgFile="sampleSets_v6_2.cfg"
sampleCCfgFile="sampleCollections_v6.cfg"

sampleSCfgFileSkim='sampleSets_v2_skim.cfg'
sampleCCfgFileSkim='sampleCollections_v2_skim.cfg'

function print_help {
    echo "Usage:"
    echo "getSamplesCfg.sh                                                 | Makes default softlinks"
    echo "getSamplesCfg.sh skims                                           | Makes default skim softlinks"
    echo "getSamplesCfg.sh -s [SampleSetConfig] -c [SampleCollectionConfg] | Makes softlinks to given SampleSet and SampleCollection (path already included)"
    echo ""
}

while getopts "h?s:c:" opt; do
    case "$opt" in
    h|\?)
        print_help
        exit 0
        ;;
    s)  sampleSCfgFile=$OPTARG
        ;;
    c)  sampleCCfgFile=$OPTARG
        ;;
    esac
done

subcommand=$1
case "$subcommand" in
    skims)
        sampleSCfgFile=$sampleSCfgFileSkim
        sampleCCfgFile=$sampleCCfgFileSkim
esac

if [ ! -f sampleSets.cfg ] && [ ! -f sampleCollections.cfg ] 
then
    ln -s $cfgDir/$sampleSCfgFile sampleSets.cfg
    ln -s $cfgDir/$sampleCCfgFile sampleCollections.cfg
    echo "Made soft link: "$sampleSCfgFile
    echo "Made soft link: "$sampleCCfgFile
else
    echo "Soft links for sampleSets.cfg and sampleCollections.cfg found"
    echo "Remove soft links for sampleSets.cfg and sampleCollections.cfg and try again"
fi
