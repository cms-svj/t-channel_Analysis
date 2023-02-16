from utils import samples as s
import os
import sys
import numpy as np

# run this on the lpg-gpu since they have more RAM

# 2018
reRunDetail = {
"2018_QCD_Pt_1000to1400": 1,
"2018_QCD_Pt_1400to1800": 1,
"2018_QCD_Pt_1800to2400": 1,
"2018_QCD_Pt_2400to3200": 1,
"2018_QCD_Pt_300to470": 200,
"2018_QCD_Pt_3200toInf": 1,
"2018_QCD_Pt_470to600": 10,
"2018_QCD_Pt_600to800": 2,
"2018_QCD_Pt_800to1000": 1,
"2018_TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8": 200,
"2018_TTJets_DiLept_genMET-150_TuneCP5_13TeV-madgraphMLM-pythia8": 10,
"2018_TTJets_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8": 2,
"2018_TTJets_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8": 1,
"2018_TTJets_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8": 5,
"2018_TTJets_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8": 3,
"2018_TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8": 60,
"2018_TTJets_SingleLeptFromT_genMET-150_TuneCP5_13TeV-madgraphMLM-pythia8": 20,
"2018_TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8": 80,
"2018_TTJets_SingleLeptFromTbar_genMET-150_TuneCP5_13TeV-madgraphMLM-pythia8": 20,
"2018_TTJets_Incl_TuneCP5_13TeV-madgraphMLM-pythia8": 290
}

batch = 10

for rerun,num in reRunDetail.items():
    sample = rerun.replace("_TuneCP5_13TeV-madgraphMLM-pythia8","")
    if num <= batch:
        command = "python analyze_root_varModule.py -d {} -N {} -w 4 -s 1000 --training NN".format(sample,num)
        print(command)
        os.system(command)
        command = "mv trainFile.root tree_{}_NN.root".format(sample.replace("2018_",""))
        os.system(command)
    else:
        for mVal in np.arange(0,num,batch):
            command = "python analyze_root_varModule.py -d {} -N 10 -M {} -w 8 -s 1000 --training NN".format(sample,mVal)
            print(command)
            os.system(command)
            command = "mv trainFile.root tree_{}_NN_M{}.root".format(sample.replace("2018_",""),mVal)
            os.system(command)
        
