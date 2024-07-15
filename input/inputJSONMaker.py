# use regular LPC environment for this, it will break in the coffea env because of the eosls command.
# create json files that store dictionary of input files of different samples
import json
import os
from inputDictionary import yearDict,backgrounds,signals,data,postpendToRemove

skimSource = False
if skimSource:
    tmNtpLocation = ""
    tmSigLocation = ""
    ntupleKind = "skims"
else:
    tmNtpLocation = "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20"
    tmSigLocation = "/store/user/lpcdarkqcd/tchannel_UL"
    ntupleKind = "treeMakerNtuples"


def getEosContent(eosSource):
    os.system(f"eosls {eosSource} > out.txt")
    with open(f"out.txt","r") as f:
        backgroundList = [sample[:-1] for sample in f.readlines()]
    os.system("rm out.txt")
    return backgroundList

def removePostpend(sampleName,sampleGroup):
    newSampleName = sampleName
    for postpend in postpendToRemove:
        newSampleName = newSampleName.replace(postpend,"")
    if sampleName == "TTJets_TuneCP5_13TeV-madgraphMLM-pythia8": # to make sure [year]_TTJets means all TTJets samples instead of the TTJets inclusive sample only
        newSampleName += "_Incl"
    if "Data" in sampleGroup:
        newSampleName += "Data"
    return newSampleName

def makeJSONFiles(sampleDict, yearDict):
    for sampleGroup, sampleList in sampleDict.items():
        year = sampleGroup[:sampleGroup.find("_")]
        print(year)
        if "Data" in sampleGroup:
            yearAliasList = yearDict[year]["data"]
        else:
            yearAliasList = yearDict[year]["background"]
        outputDir = f"sampleJSONs/{ntupleKind}/{sampleGroup}"
        os.makedirs(outputDir,exist_ok=True)
        for sampleName in sampleList:
            print(sampleName)
            sampleFileName = f"{year}_{removePostpend(sampleName,sampleGroup)}"
            allSampleLocationList = []
            for yearAlias in yearAliasList:
                sampleNameList = getEosContent(f"{tmNtpLocation}/{yearAlias}")
                if sampleName in sampleNameList:
                    sampleFileList = getEosContent(f"{tmNtpLocation}/{yearAlias}/{sampleName}")
                    print(f"{len(sampleFileList)} files")
                    sampleLocationList = [f"root://cmseos.fnal.gov/{tmNtpLocation}/{yearAlias}/{sampleName}/{sampleFile}" for sampleFile in sampleFileList]
                    allSampleLocationList += sampleLocationList
                    print()

            with open(f"{outputDir}/{sampleFileName}.json","w") as fp:
                json.dump({sampleFileName:allSampleLocationList}, fp, indent=4)

makeJSONFiles(backgrounds, yearDict)
makeJSONFiles(data, yearDict)

