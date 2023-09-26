from utils import samples as s
import os
import sys
import numpy as np
from glob import glob
# run this on the lpg-gpu since they have more RAM

def fileExists(existingFiles,fileName):
    fileExists = False
    for exFile in existingFiles:
        if fileName in exFile:
            fileExists = True
            break
    if fileExists == False:
        print(fileName)
    return fileExists


tcut = "_pre" # _pre or _pre_1PSVJ
outHistF = "output/eventTaggerOutput/trainingFiles/wpt0p9{}/".format(tcut)
resume = True

# 2018
# pre
if tcut == "_pre":
    reRunDetail = {
    "2018_QCD_Pt_300to470": 200,
    "2018_QCD_Pt_470to600": 10,
    "2018_QCD_Pt_600to800": 2,
    "2018_QCD_Pt_800to1000": 1,
    "2018_QCD_Pt_1000to1400": 1,
    "2018_QCD_Pt_1400to1800": 1,
    "2018_QCD_Pt_1800to2400": 1,
    "2018_QCD_Pt_2400to3200": 1,
    "2018_QCD_Pt_3200toInf": 1,
    "2018_TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8": 400, # 200
    "2018_TTJets_DiLept_genMET-150_TuneCP5_13TeV-madgraphMLM-pythia8": 20, # 10
    "2018_TTJets_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8": 10, # 5
    "2018_TTJets_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8": 6, # 3
    "2018_TTJets_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8": 4, # 2
    "2018_TTJets_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8": 2, # 1
    "2018_TTJets_Incl_TuneCP5_13TeV-madgraphMLM-pythia8": 290, # 290
    "2018_TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8": 120, # 60
    "2018_TTJets_SingleLeptFromT_genMET-150_TuneCP5_13TeV-madgraphMLM-pythia8": 40, # 20
    "2018_TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8": 160, # 80
    "2018_TTJets_SingleLeptFromTbar_genMET-150_TuneCP5_13TeV-madgraphMLM-pythia8": 40, # 20
    }
    batch = 10
elif tcut == "_pre_1PSVJ":
    reRunDetail = {
    '2018_QCD_Pt_300to470': 200, 
    '2018_QCD_Pt_470to600': 20, 
    '2018_QCD_Pt_600to800': 10, 
    '2018_QCD_Pt_800to1000': 5, 
    '2018_QCD_Pt_1000to1400': 5, 
    '2018_QCD_Pt_1400to1800': 5, 
    '2018_QCD_Pt_2400to3200': 5, 
    '2018_QCD_Pt_1800to2400': 5,
    '2018_TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8': 300,  
    '2018_TTJets_DiLept_genMET-150_TuneCP5_13TeV-madgraphMLM-pythia8': 100, 
    '2018_TTJets_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8': 50, 
    '2018_TTJets_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8': 30, 
    '2018_TTJets_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8': 10, 
    '2018_TTJets_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8': 20, 
    '2018_TTJets_Incl_TuneCP5_13TeV-madgraphMLM-pythia8': 290, 
    '2018_TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8': 500, 
    '2018_TTJets_SingleLeptFromT_genMET-150_TuneCP5_13TeV-madgraphMLM-pythia8': 200, 
    '2018_TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8': 500, 
    '2018_TTJets_SingleLeptFromTbar_genMET-150_TuneCP5_13TeV-madgraphMLM-pythia8': 200, 
    }
    batch = 20

existingFiles = glob("{}/*.root".format(outHistF))

for rerun,num in reRunDetail.items():
    sample = rerun.replace("_TuneCP5_13TeV-madgraphMLM-pythia8","")
    sampleName = sample.replace("_TuneCP5_13TeV-madgraphMLM-pythia8","")
    if num <= batch:
        fileName = "tree_{}_M0".format(sampleName)
        if resume:
            if (fileExists(existingFiles,fileName) == False):
                command = "python analyze_root_varModule.py -d {} -N {} -w 16 -s 1000 --training NN --outHistF {} --tcut {}".format(sample,num,outHistF,tcut)
                print(command)
                os.system(command)
        else:
            command = "python analyze_root_varModule.py -d {} -N {} -w 16 -s 1000 --training NN --outHistF {} --tcut {}".format(sample,num,outHistF,tcut)
            print(command)
            os.system(command)

    else:
        for mVal in np.arange(0,num,batch):
            fileName = "tree_{}_M{}".format(sampleName,mVal)
            NVal = batch
            if num - mVal < batch:
                NVal = num-mVal
            command = "python analyze_root_varModule.py -d {} -N {} -M {} -w 16 -s 1000 --training NN --outHistF {} --tcut {}".format(sample,NVal,mVal,outHistF,tcut)
            if resume:
                if (fileExists(existingFiles,fileName) == False):
                    print(command)
                    os.system(command)
        
            else:
                print(command)
                os.system(command)