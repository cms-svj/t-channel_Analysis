from os import system, environ
import json
import glob

def getFileset(sample,verbose=False,startFile=0,nFiles=-1,mlTraining=False):
    # find all json files for the sample or sample collection
    f_ = sample.find("_")
    year = sample[:f_]
    detailKey = sample[f_+1:]
    kind = "backgrounds"
    if ("mMed" in detailKey) or ("mZprime" in detailKey):
        kind = "signals"
    dataKeys = ["HTMHTData","JetHTData","METData","SingleElectronData","SingleMuonData","SinglePhotonData","EGammaData"]
    for dKey in dataKeys:
        if "{}".format(dKey) in detailKey:
            kind = "data"
            break
    if "Incl" in detailKey:
        ii = detailKey.find("Incl")
        detailKey = detailKey[:ii] + "Tune"
    elif "TTJets_SingleLeptFromT" in sample or "TTJets_DiLept" in detailKey:
        detailKey += "_Tune"

    JSONDir = environ['TCHANNEL_BASE'] + '/input/sampleJSONs/' + kind + "/" + year + "/"

    # inputSamples = glob.glob(JSONDir+"*"+detailKey+"*.json")
    inputSamples = glob.glob(JSONDir+year+"_"+detailKey+"*.json")
    if len(inputSamples) == 0:
        print("Error: no json file found with name:", JSONDir)
    # else:
    #     print(inputSamples)

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

def getAllFilesets():
    fBG = getFileset("*_*", False)
    fSG1 = getFileset("*_mMed", False)
    fSG2 = getFileset("*_mZprime", False)
    fBG.update(fSG1)
    fBG.update(fSG2)
    return fBG

def sfGetter(sample,on=False,tcut=""):
    if on:
        # scaleFactor for particleNet/event tagger training files (pre)
        if tcut == "_pre":
            scaleFactor = {
                "2018_QCD_Pt_1000to1400": 19730000/48000.,
                "2018_QCD_Pt_1400to1800": 10982000/21000.,
                "2018_QCD_Pt_1800to2400": 5491000/48000.,
                "2018_QCD_Pt_2400to3200": 2997000/30000.,
                "2018_QCD_Pt_300to470": 57910000/9423000.,
                "2018_QCD_Pt_3200toInf": 1000000/3000.,
                "2018_QCD_Pt_470to600": 52448000/447000.,
                "2018_QCD_Pt_600to800": 67508000/96000.,
                "2018_QCD_Pt_800to1000": 37160000/48000.,
                "2018_TTJets_DiLept": 29290487/16993527.,
                "2018_TTJets_DiLept_genMET-150": 10592111/195331.,
                "2018_TTJets_HT-600to800": 15258099/179188.,
                "2018_TTJets_HT-800to1200": 9201990/135827.,
                "2018_TTJets_HT-1200to2500": 2009331/70237.,
                "2018_TTJets_HT-2500toInf": 1001084/40694.,
                "2018_TTJets_Incl": 10494353/10339641.,
                "2018_TTJets_SingleLeptFromT": 58237254/5080853.,
                "2018_TTJets_SingleLeptFromT_genMET-150": 13337428/845757.,
                "2018_TTJets_SingleLeptFromTbar": 58510607/6941755.,
                "2018_TTJets_SingleLeptFromTbar_genMET-150": 12597052/968711.,
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
                "2018_TTJets_Incl": 10494353/672958.
                }
        sf = 1
        for key,item in scaleFactor.items():
            if key == sample:
                sf = item
                break
    else:
        sf = 1
    return sf
