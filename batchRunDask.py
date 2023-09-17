from utils import samples as s
import os
import sys
import numpy as np


# 2016
# reRunDetail='''2016_QCD_Pt_1000to1400
# 2016_QCD_Pt_120to170
# 2016_QCD_Pt_1400to1800
# 2016_QCD_Pt_170to300
# 2016_QCD_Pt_1800to2400
# 2016_QCD_Pt_2400to3200
# 2016_QCD_Pt_300to470
# 2016_QCD_Pt_3200toInf
# 2016_QCD_Pt_470to600
# 2016_QCD_Pt_600to800
# 2016_QCD_Pt_800to1000
# 2016_QCD_Pt_80to120
# 2016_ST_s-channel_4f_hadronicDecays
# 2016_ST_s-channel_4f_leptonDecays
# 2016_ST_t-channel_antitop
# 2016_ST_t-channel_top
# 2016_ST_tW_antitop
# 2016_ST_tW_top
# 2016_ST_tZq_ll_4f_ckm
# 2016_TTJets_DiLept_genMET-150
# 2016_TTJets_DiLept
# 2016_TTJets_HT-1200to2500
# 2016_TTJets_HT-2500toInf
# 2016_TTJets_HT-600to800
# 2016_TTJets_HT-800to1200
# 2016_TTJets_SingleLeptFromTbar_genMET-150
# 2016_TTJets_SingleLeptFromTbar
# 2016_TTJets_SingleLeptFromT_genMET-150
# 2016_TTJets_SingleLeptFromT
# 2016_TTJets_Incl
# 2016_WJetsToLNu_HT-100To200
# 2016_WJetsToLNu_HT-1200To2500
# 2016_WJetsToLNu_HT-200To400
# 2016_WJetsToLNu_HT-2500ToInf
# 2016_WJetsToLNu_HT-400To600
# 2016_WJetsToLNu_HT-600To800
# 2016_WJetsToLNu_HT-800To1200
# 2016_ZJetsToNuNu_HT-100To200
# 2016_ZJetsToNuNu_HT-1200To2500
# 2016_ZJetsToNuNu_HT-200To400
# 2016_ZJetsToNuNu_HT-2500ToInf
# 2016_ZJetsToNuNu_HT-400To600
# 2016_ZJetsToNuNu_HT-600To800
# 2016_ZJetsToNuNu_HT-800To1200
# 2016_HTMHTData
# 2016_JetHTData
# 2016_METData'''.splitlines()

# 2018
# reRunDetail='''2018_Skim_METData
# 2018_JetHTData'''.splitlines()



# 2017
# reRunDetail='''2017_QCD_Pt_1000to1400'''.splitlines()
# 2017_QCD_Pt_120to170
# 2017_QCD_Pt_1400to1800
# 2017_QCD_Pt_170to300
# 2017_QCD_Pt_1800to2400
# 2017_QCD_Pt_2400to3200
# 2017_QCD_Pt_300to470
# 2017_QCD_Pt_3200toInf
# 2017_QCD_Pt_470to600
# 2017_QCD_Pt_600to800
# 2017_QCD_Pt_800to1000
# 2017_QCD_Pt_80to120
# 2017_ST_s-channel_4f_hadronicDecays
# 2017_ST_s-channel_4f_leptonDecays
# 2017_ST_t-channel_antitop
# 2017_ST_t-channel_top
# 2017_ST_tW_antitop
# 2017_ST_tW_top
# 2017_ST_tZq_ll_4f_ckm
reRunDetail='''2017_TTJets_DiLept_genMET-150_TuneCP5_13TeV-madgraphMLM-pythia8
2017_TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8
2017_TTJets_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8
2017_TTJets_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8
2017_TTJets_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8
2017_TTJets_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8
2017_TTJets_SingleLeptFromTbar_genMET-150_TuneCP5_13TeV-madgraphMLM-pythia8
2017_TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8
2017_TTJets_SingleLeptFromT_genMET-150_TuneCP5_13TeV-madgraphMLM-pythia8
2017_TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8
2017_TTJets_Incl'''.splitlines()
# 2017_WJetsToLNu_HT-100To200
# 2017_WJetsToLNu_HT-1200To2500
# 2017_WJetsToLNu_HT-200To400
# 2017_WJetsToLNu_HT-2500ToInf
# 2017_WJetsToLNu_HT-400To600
# 2017_WJetsToLNu_HT-600To800
# 2017_WJetsToLNu_HT-800To1200
# 2017_ZJetsToNuNu_HT-100To200
# 2017_ZJetsToNuNu_HT-1200To2500
# 2017_ZJetsToNuNu_HT-200To400
# 2017_ZJetsToNuNu_HT-2500ToInf
# 2017_ZJetsToNuNu_HT-400To600
# 2017_ZJetsToNuNu_HT-600To800
# 2017_ZJetsToNuNu_HT-800To1200
# 2017_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1'''.splitlines()
# 2017_HTMHTData
# 2017_JetHTData
# 2017_METData'''.splitlines()


# reRunDetail = '''2018_QCD_Pt_300to470
# 2018_TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8'''.splitlines()
# reRunDetail = '''2018_QCD_Pt_1000to1400
# 2018_QCD_Pt_1400to1800
# 2018_QCD_Pt_1800to2400
# 2018_QCD_Pt_2400to3200
# 2018_QCD_Pt_300to470
# 2018_QCD_Pt_3200toInf
# 2018_QCD_Pt_470to600
# 2018_QCD_Pt_600to800
# 2018_QCD_Pt_800to1000
# 2018_TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8
# 2018_TTJets_DiLept_genMET-150_TuneCP5_13TeV-madgraphMLM-pythia8
# 2018_TTJets_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8
# 2018_TTJets_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8
# 2018_TTJets_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8
# 2018_TTJets_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8
# 2018_TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8
# 2018_TTJets_SingleLeptFromT_genMET-150_TuneCP5_13TeV-madgraphMLM-pythia8
# 2018_TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8
# 2018_TTJets_SingleLeptFromTbar_genMET-150_TuneCP5_13TeV-madgraphMLM-pythia8
# 2018_TTJets_Incl_TuneCP5_13TeV-madgraphMLM-pythia8'''.splitlines()

nFilesPerJob = 4
maxJobs = 50
maxFiles = nFilesPerJob*maxJobs

for rerun in reRunDetail:
    sample = rerun.replace("_TuneCP5_13TeV-madgraphMLM-pythia8","")
    print(sample)
    fileset = s.getFileset(sample, True)
    totalNumberOfFiles = len(list(fileset.items())[0][1])
    numberOfJobsRequired = int(np.ceil(totalNumberOfFiles/nFilesPerJob))
    print("totalNumberOfFiles",totalNumberOfFiles)
    print("numberOfJobsRequired",numberOfJobsRequired)
    if numberOfJobsRequired > maxJobs:
        for mVal in np.arange(0,totalNumberOfFiles,maxFiles):
            if numberOfJobsRequired >= maxJobs:
                command = "python analyze.py -d {} -N 250 -M {} -b {} -s 1000 --condor --dask -j".format(sample,mVal,maxJobs)
                print(command)
                os.system(command)
                numberOfJobsRequired -= maxJobs
            else:
                command = "python analyze.py -d {} -M {} -b {} -s 1000 --condor --dask -j".format(sample,mVal,numberOfJobsRequired)
                print(command)
                os.system(command)
    else:
        command = "python analyze.py -d {} -b {} -s 1000 --condor --dask -j".format(sample,numberOfJobsRequired)
        print(command)
        os.system(command)
        
