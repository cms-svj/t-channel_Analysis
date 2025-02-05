# use regular LPC environment for this, it will break in the coffea env because of the eosls command.
# create json files that store dictionary of input files of different samples
import json
import os

skimSource = True
if skimSource:
    from inputDictionary_skim import backgrounds,signals,data
else:
    from inputDictionary import yearDict,backgrounds,signals,data,postpendToRemove

# options are: t_channel_pre_selection, t_channel_lost_lepton_control_region (see root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims/{year}/})
skimModule = "t_channel_pre_selection" 
# skimModule = "t_channel_lost_lepton_control_region" 
specialNote = "" # use "" for normal circumstances

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

def cleanUpEmptyFolders(outputDir):
    content = os.listdir(outputDir)
    if len(content) == 0:
        os.system(f"rm -r {outputDir}")


def makeJSONFiles(sampleDict):
    sampleGroupLabelList = []
    sampleLabelList = []
    for sampleGroup, sampleList in sampleDict.items():
        year = sampleGroup[:sampleGroup.find("_")]
        if skimSource:
            skimLocation = f"/store/user/lpcdarkqcd/tchannel_UL/skims/{year}/{skimModule}/nominal/"
            ntupleKind = "skims"
        else:
            dataBkgLocation = "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20"
            sigLocation = "/store/user/lpcdarkqcd/tchannel_UL"
            ntupleKind = "treeMakerNtuples"
        print(year)
        sampleGroupLabelList.append(sampleGroup)
        if skimSource:
            outputDir = f"sampleJSONs/{ntupleKind}/{skimModule}{specialNote}/{sampleGroup}"
        else:
            if "Data" in sampleGroup:
                yearAliasList = yearDict[year]["data"]
            elif "SVJ" in sampleGroup: # signal
                yearAliasList = [year]
            else:
                yearAliasList = yearDict[year]["background"]
            outputDir = f"sampleJSONs/{ntupleKind}/{sampleGroup}"
        os.makedirs(outputDir,exist_ok=True)
        for sampleName in sampleList:
            print(sampleName)
            if skimSource:
                sampleFileName = f"{year}_{sampleName}"
                if sampleName == "TTJets":
                    sampleFileName += "_Incl"
                elif "SVJ" in sampleGroup:
                    sampleFileName = sampleFileName.replace("t-channel_","")
                elif "Data" in sampleGroup:
                    sampleFileName += "Data"
                sampleFileList = getEosContent(f"{skimLocation}/{sampleName}")
                if len(sampleFileList) > 0:
                    sampleLocationList = [f"root://cmseos.fnal.gov/{skimLocation}/{sampleName}/{sampleFile}" for sampleFile in sampleFileList]
                    print(f"{outputDir}/{sampleFileName}.json")
                    with open(f"{outputDir}/{sampleFileName}.json","w") as fp:
                        json.dump({sampleFileName:sampleLocationList}, fp, indent=4)
            else:
                sampleFileName = f"{year}_{removePostpend(sampleName,sampleGroup)}"
                allSampleLocationList = []
                for yearAlias in yearAliasList:
                    if "SVJ" in sampleGroup: # signal
                        sampleGroupLocation = f"{sigLocation}/{yearAlias}/Full/PrivateSamples"
                    else:
                        sampleGroupLocation = f"{dataBkgLocation}/{yearAlias}"
                    sampleNameList = getEosContent(sampleGroupLocation)
                    if sampleName in sampleNameList:
                        sampleFileList = getEosContent(f"{sampleGroupLocation}/{sampleName}")
                        print(f"{len(sampleFileList)} files")
                        sampleLocationList = [f"root://cmseos.fnal.gov/{sampleGroupLocation}/{sampleName}/{sampleFile}" for sampleFile in sampleFileList]
                        allSampleLocationList += sampleLocationList
                print(f"{outputDir}/{sampleFileName}.json")
                sampleLabelList.append(sampleFileName)
                print()

                if len(allSampleLocationList) > 0:
                    with open(f"{outputDir}/{sampleFileName}.json","w") as fp:
                        json.dump({sampleFileName:allSampleLocationList}, fp, indent=4)
        cleanUpEmptyFolders(outputDir)

makeJSONFiles(backgrounds)
makeJSONFiles(data)
makeJSONFiles(signals)
