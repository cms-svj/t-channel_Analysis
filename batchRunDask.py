from utils import samples as s
import os
import sys
import numpy as np

# 2018
reRunDetail = '''2018_QCD_Pt_300to470
2018_TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8'''.splitlines()
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

nFilesPerJob = 5
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
                command = "python analyze.py -d {} -N 250 -M {} -b {} -s 1000 --condor --dask".format(sample,mVal,maxJobs)
                print(command)
                os.system(command)
                numberOfJobsRequired -= maxJobs
            else:
                command = "python analyze.py -d {} -M {} -b {} -s 1000 --condor --dask".format(sample,mVal,numberOfJobsRequired)
                print(command)
                os.system(command)
    else:
        command = "python analyze.py -d {} -b {} -s 1000 --condor --dask".format(sample,numberOfJobsRequired)
        print(command)
        os.system(command)
        
