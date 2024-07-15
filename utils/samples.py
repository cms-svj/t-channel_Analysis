import os
import json
from glob import glob

def getSamplesFromGroup(sampleGroup,skimSource=False):
    ntupleKind = "treeMakerNtuples"
    if skimSource:
        ntupleKind = "skims"
    sampleInputFolder = f"{os.getcwd()}/input/sampleJSONs/{ntupleKind}/{sampleGroup}"
    return [sample.replace(".json","") for sample in os.listdir(sampleInputFolder)]    

def getGroupFromSample(sample,skimSource=False):
    ntupleKind = "treeMakerNtuples"
    if skimSource:
        ntupleKind = "skims"
    inputFolder = f"{os.getcwd()}/input/sampleJSONs/{ntupleKind}/"
    allSampleGroups = os.listdir(inputFolder)
    for sampleGroup in allSampleGroups:
        allSamples = os.listdir(f"{inputFolder}/{sampleGroup}")
        if f"{sample}.json" in allSamples:
            return sampleGroup
    raise Exception(f"{sample} does not belong to any existing sample group.")

def getFileset(sample,startFile=0,nFiles=-1,skimSource=False,verbose=False):
    # find all json files for the sample or sample collection
    ntupleKind = "treeMakerNtuples"
    if skimSource:
        ntupleKind = "skims"
    sampleInputFolder = f"{os.getcwd()}/input/sampleJSONs/{ntupleKind}/"
    allSampleGroups = os.listdir(sampleInputFolder)
    inputSamples = []
    if sample in allSampleGroups:
        inputSamples = glob(f"{sampleInputFolder}/{sample}/*.json")
    else:
        allFiles = glob(f"{sampleInputFolder}/**/*.json",recursive=True)
        for f in allFiles:
            if f"{sample}.json" in f:
                inputSamples.append(f)
                break
    if len(inputSamples) == 0:
        raise Exception(f"{sample} does not exist.")

    # open all json files and dump them into a dictionary
    fileset_all = {}
    for s in inputSamples:
        fileset_all.update(json.load(open(s ,'r')))

    # process a subset of a sample
    fileset = {}
    for n, rFiles in fileset_all.items():
        nf = nFiles
        fileset[n] = []

        if nf < 0:
            nf = len(rFiles)
        fn = startFile

        while fn < startFile+nf and fn < len(rFiles):
            # print("Name of sample: %-40s" % (n), "start file number: " + str(fn), "Number of total files: " + str(len(rFiles)))
            fileset[n].append(rFiles[fn])
            fn+=1

    # print the sample names if verbose
    if verbose:
        for key, _ in fileset.items():
            print(key)

    return fileset

def getFilesetFromList(sampleList,options,verbose=False):
    startFileList = options.startFile 
    nFilesList = options.nFiles
    skimSource = options.skimSource
    allFileset = {}
    for i in range(len(sampleList)):
        sample = sampleList[i]
        startFile = startFileList[i]
        nFiles = nFilesList[i]
        fileset = getFileset(sample,startFile,nFiles,skimSource,verbose)
        allFileset.update(fileset)
    return allFileset

def getAllFilesets():
    fBG = getFileset("*_*", False)
    fSG1 = getFileset("*_mMed", False)
    fSG2 = getFileset("*_mZprime", False)
    fBG.update(fSG1)
    fBG.update(fSG2)
    return fBG

def sfGetter(sample,scaleOn=False,tcut=""):
    if scaleOn:
        # scaleFactor for particleNet/event tagger training files (pre)
        # if tcut == "_pre":
            # scaleFactor = {
            #     "2018_QCD_Pt_1000to1400": 19730000/48000.,
            #     "2018_QCD_Pt_1400to1800": 10982000/21000.,
            #     "2018_QCD_Pt_1800to2400": 5491000/48000.,
            #     "2018_QCD_Pt_2400to3200": 2997000/30000.,
            #     "2018_QCD_Pt_300to470": 57910000/9423000.,
            #     "2018_QCD_Pt_3200toInf": 1000000/3000.,
            #     "2018_QCD_Pt_470to600": 52448000/447000.,
            #     "2018_QCD_Pt_600to800": 67508000/96000.,
            #     "2018_QCD_Pt_800to1000": 37160000/48000.,
            #     "2018_TTJets_DiLept": 29290487/16993527.,
            #     "2018_TTJets_DiLept_genMET-150": 10592111/195331.,
            #     "2018_TTJets_HT-600to800": 15258099/179188.,
            #     "2018_TTJets_HT-800to1200": 9201990/135827.,
            #     "2018_TTJets_HT-1200to2500": 2009331/70237.,
            #     "2018_TTJets_HT-2500toInf": 1001084/40694.,
            #     "2018_TTJets_Incl": 10494353/10339641.,
            #     "2018_TTJets_SingleLeptFromT": 58237254/5080853.,
            #     "2018_TTJets_SingleLeptFromT_genMET-150": 13337428/845757.,
            #     "2018_TTJets_SingleLeptFromTbar": 58510607/6941755.,
            #     "2018_TTJets_SingleLeptFromTbar_genMET-150": 12597052/968711.,
            #     }
        if tcut == "_pre": # this version increases the number of high MET qcd samples
            scaleFactor = {
                "2018_QCD_Pt_300to470": 57910000/9423000.,
                "2018_QCD_Pt_470to600": 52448000/927000.,
                "2018_QCD_Pt_600to800": 67508000/960000.,
                "2018_QCD_Pt_800to1000": 37160000/906000.,
                "2018_QCD_Pt_1000to1400": 19730000/936000.,
                "2018_QCD_Pt_1400to1800": 10982000/843000.,
                "2018_QCD_Pt_1800to2400": 5491000/654000.,
                "2018_QCD_Pt_2400to3200": 2997000/786000.,
                "2018_QCD_Pt_3200toInf": 1000000/484000.,
                }
        # scaleFactor for particleNet/event tagger training files (pre)
        elif tcut == "_pre_1PSVJ":
            scaleFactor = {
                "2018_QCD_Pt_1000to1400": 19730000/225000.,
                "2018_QCD_Pt_1400to1800": 10982000/213000.,
                "2018_QCD_Pt_1800to2400": 5491000/183000.,
                "2018_QCD_Pt_2400to3200": 2997000/207000.,
                "2018_QCD_Pt_300to470": 57910000/9423000.,
                "2018_QCD_Pt_3200toInf": 1000000/195000.,
                "2018_QCD_Pt_470to600": 52448000/927000.,
                "2018_QCD_Pt_600to800": 67508000/480000.,
                "2018_QCD_Pt_800to1000": 37160000/207000.,
                "2018_TTJets_DiLept": 29290487/12663839.,
                "2018_TTJets_DiLept_genMET-150": 10592111/1870456.,
                "2018_TTJets_HT-1200to2500": 2009331/591459.,
                "2018_TTJets_HT-2500toInf": 1001084/129330.,
                "2018_TTJets_HT-600to800": 15258099/974946.,
                "2018_TTJets_HT-800to1200": 9201990/611853.,
                "2018_TTJets_SingleLeptFromT": 58237254/19505097.,
                "2018_TTJets_SingleLeptFromT_genMET-150": 13337428/3905981.,
                "2018_TTJets_SingleLeptFromTbar": 58510607/18798625.,
                "2018_TTJets_SingleLeptFromTbar_genMET-150": 12597052/4489650.,
                "2018_TTJets_Incl": 10494353/10339641.,
                }
        elif tcut == "_pre_MET50":
            scaleFactor = {
                "2018_QCD_Pt_300to470": 57910000/7023000.,
                "2018_QCD_Pt_470to600": 52448000/591000.,
                "2018_QCD_Pt_600to800": 67508000/192000.,
                "2018_QCD_Pt_800to1000": 37160000/48000.,
                "2018_QCD_Pt_1000to1400": 19730000/48000.,
                "2018_QCD_Pt_1400to1800": 10982000/21000.,
                "2018_QCD_Pt_1800to2400": 5491000/48000.,
                "2018_QCD_Pt_2400to3200": 2997000/30000.,
                }
        elif tcut == "_pre_MET75":
            scaleFactor = {
                "2018_QCD_Pt_300to470": 57910000/11749000.,
                "2018_QCD_Pt_470to600": 52448000/1023000.,
                "2018_QCD_Pt_600to800": 67508000/288000.,
                "2018_QCD_Pt_800to1000": 37160000/48000.,
                "2018_QCD_Pt_1000to1400": 19730000/48000.,
                "2018_QCD_Pt_1400to1800": 10982000/21000.,
                "2018_QCD_Pt_1800to2400": 5491000/48000.,
                "2018_QCD_Pt_2400to3200": 2997000/30000.,
                }
        elif tcut == "_pre_MET100":
            scaleFactor = {
                "2018_QCD_Pt_300to470": 57910000/21457000.,
                "2018_QCD_Pt_470to600": 52448000/1860000.,
                "2018_QCD_Pt_600to800": 67508000/528000.,
                "2018_QCD_Pt_800to1000": 37160000/96000.,
                "2018_QCD_Pt_1000to1400": 19730000/48000.,
                "2018_QCD_Pt_1400to1800": 10982000/21000.,
                "2018_QCD_Pt_1800to2400": 5491000/48000.,
                "2018_QCD_Pt_2400to3200": 2997000/30000.,
                }
        elif tcut == "_pre_MET125":
            scaleFactor = {
                "2018_QCD_Pt_300to470": 57910000/39730000.,
                "2018_QCD_Pt_470to600": 52448000/3444000.,
                "2018_QCD_Pt_600to800": 67508000/912000.,
                "2018_QCD_Pt_800to1000": 37160000/159000.,
                "2018_QCD_Pt_1000to1400": 19730000/48000.,
                "2018_QCD_Pt_1400to1800": 10982000/21000.,
                "2018_QCD_Pt_1800to2400": 5491000/48000.,
                "2018_QCD_Pt_2400to3200": 2997000/30000.,
                }
        elif tcut == "_pre_MET200":
            scaleFactor = {
                "2018_QCD_Pt_300to470": 57910000/57910000.,
                "2018_QCD_Pt_470to600": 52448000/12885000.,
                "2018_QCD_Pt_600to800": 67508000/3624000.,
                "2018_QCD_Pt_800to1000": 37160000/591000.,
                "2018_QCD_Pt_1000to1400": 19730000/177000.,
                "2018_QCD_Pt_1400to1800": 10982000/21000.,
                "2018_QCD_Pt_1800to2400": 5491000/48000.,
                "2018_QCD_Pt_2400to3200": 2997000/30000.,
                }

        # for 20 files per sample, used for small subset of samples to make histograms
        elif tcut == "":
            scaleFactor = {
                "2018_QCD_Pt_80to120": 29685000/1899000.,
                "2018_QCD_Pt_120to170": 29949000/1824000.,
                "2018_QCD_Pt_170to300": 29676000/960000.,
                "2018_QCD_Pt_300to470": 57910000/870000.,
                "2018_QCD_Pt_470to600": 52448000/927000.,
                "2018_QCD_Pt_600to800": 67508000/960000.,
                "2018_QCD_Pt_800to1000": 37160000/906000.,
                "2018_QCD_Pt_1000to1400": 19730000/936000.,
                "2018_QCD_Pt_1400to1800": 10982000/843000.,
                "2018_QCD_Pt_1800to2400": 5491000/654000.,
                "2018_QCD_Pt_2400to3200": 2997000/786000.,
                "2018_QCD_Pt_3200toInf": 1000000/484000.,
                "2018_TTJets_DiLept": 29290487/871719.,
                "2018_TTJets_DiLept_genMET-150": 10592111/195331.,
                "2018_TTJets_HT-1200to2500": 2009331/591459.,
                "2018_TTJets_HT-2500toInf": 1001084/184038.,
                "2018_TTJets_HT-600to800": 15258099/345290.,
                "2018_TTJets_HT-800to1200": 9201990/392841.,
                "2018_TTJets_SingleLeptFromT": 58237254/835004.,
                "2018_TTJets_SingleLeptFromT_genMET-150": 13337428/414361.,
                "2018_TTJets_SingleLeptFromTbar": 58510607/863322.,
                "2018_TTJets_SingleLeptFromTbar_genMET-150": 12597052/503985.,
                "2018_TTJets_Incl": 10494353/672958.,
                "2017_QCD_Pt_1000to1400": 19781000/1374000,
                "2017_QCD_Pt_120to170": 28896000/1356000,
                "2017_QCD_Pt_1400to1800": 10994000/1566000,
                "2017_QCD_Pt_170to300": 29811000/1680000,
                "2017_QCD_Pt_1800to2400": 5488000/1095000,
                "2017_QCD_Pt_2400to3200": 2997000/828000,
                "2017_QCD_Pt_300to470": 55690000/1449000,
                "2017_QCD_Pt_3200toInf": 1000000/889000,
                "2017_QCD_Pt_470to600": 50885000/1569000,
                "2017_QCD_Pt_600to800": 67379000/1662000,
                "2017_QCD_Pt_800to1000": 36890000/1635000,
                "2017_QCD_Pt_80to120": 29403000/1614000,
                "2017_ST_s-channel_4f_hadronicDecays": 11696999/1125000,
                "2017_ST_s-channel_4f_leptonDecays": 13882000/1414000,
                "2017_ST_t-channel_antitop": 70203000/1133000,
                "2017_ST_t-channel_top": 141910000/1605000,
                "2017_ST_tW_antitop": 5674000/645000,
                "2017_ST_tW_top": 5649000/562000,
                "2017_ST_tZq_ll_4f_ckm": 9530000/1490000.,
                "2017_TTJets_DiLept_genMET-150": 9565099/408290.,
                "2017_TTJets_DiLept": 30187083/645469,
                "2017_TTJets_HT-1200to2500": 2039938/301687,
                "2017_TTJets_HT-2500toInf": 981339/351243,
                "2017_TTJets_HT-600to800": 15594718/388177,
                "2017_TTJets_HT-800to1200": 10410425/412379,
                "2017_TTJets_SingleLeptFromTbar_genMET-150": 12181705/264579,
                "2017_TTJets_SingleLeptFromTbar": 58683345/584633,
                "2017_TTJets_SingleLeptFromT_genMET-150": 14094979/426582,
                "2017_TTJets_SingleLeptFromT": 58846842/239983,
                "2017_TTJets_Incl": 10052151/537399,
                "2017_WJetsToLNu_HT-100To200": 47424468/1552964,
                "2017_WJetsToLNu_HT-1200To2500": 4955636/934696,
                "2017_WJetsToLNu_HT-200To400": 42602407/1678203,
                "2017_WJetsToLNu_HT-2500ToInf": 1185699/694391,
                "2017_WJetsToLNu_HT-400To600": 5468473/929762,
                "2017_WJetsToLNu_HT-600To800": 5545298/1017485,
                "2017_WJetsToLNu_HT-800To1200": 5088483/570667,
                "2017_ZJetsToNuNu_HT-100To200": 18983897/1182572,
                "2017_ZJetsToNuNu_HT-1200To2500": 267125/267125,
                "2017_ZJetsToNuNu_HT-200To400": 17349597/1380671,
                "2017_ZJetsToNuNu_HT-2500ToInf": 176201/176201,
                "2017_ZJetsToNuNu_HT-400To600": 13963690/1367374,
                "2017_ZJetsToNuNu_HT-600To800": 4418971/947483,
                "2017_ZJetsToNuNu_HT-800To1200": 1513585/1029984,
                "2017_JetHTData": 388179354/45608920., ## For data the scales are for 500 files 
                "2017_METData": 411664094/57158579.,
                "2017_HTMHTData": 41371874/34117296.,
                }
        # for 20 files per sample, used for small subset of samples to make histograms
        else:
            scaleFactor = {
                "2018_QCD_Pt_80to120": 29685000/1899000.,
                "2018_QCD_Pt_120to170": 29949000/1824000.,
                "2018_QCD_Pt_170to300": 29676000/960000.,
                "2018_QCD_Pt_1000to1400": 19730000/936000.,
                "2018_QCD_Pt_1400to1800": 10982000/843000.,
                "2018_QCD_Pt_1800to2400": 5491000/654000.,
                "2018_QCD_Pt_2400to3200": 2997000/786000.,
                "2018_QCD_Pt_300to470": 57910000/870000.,
                "2018_QCD_Pt_3200toInf": 1000000/484000.,
                "2018_QCD_Pt_470to600": 52448000/927000.,
                "2018_QCD_Pt_600to800": 67508000/960000.,
                "2018_QCD_Pt_800to1000": 37160000/906000.,
                "2018_TTJets_DiLept": 29290487/871719.,
                "2018_TTJets_DiLept_genMET-150": 10592111/195331.,
                "2018_TTJets_HT-1200to2500": 2009331/591459.,
                "2018_TTJets_HT-2500toInf": 1001084/184038.,
                "2018_TTJets_HT-600to800": 15258099/345290.,
                "2018_TTJets_HT-800to1200": 9201990/392841.,
                "2018_TTJets_SingleLeptFromT": 58237254/835004.,
                "2018_TTJets_SingleLeptFromT_genMET-150": 13337428/414361.,
                "2018_TTJets_SingleLeptFromTbar": 58510607/863322.,
                "2018_TTJets_SingleLeptFromTbar_genMET-150": 12597052/503985.,
                "2018_TTJets_Incl": 10494353/672958.,
                
                }
        sf = 1
        for key,item in scaleFactor.items():
            if key == sample:
                sf = item
                break
    else:
        sf = 1
    return sf