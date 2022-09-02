from os import system, environ
import json
import glob

def getFileset(sample,verbose=True,startFile=0,nFiles=-1,mlTraining=False):
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

    if mlTraining:
        JSONDir = environ['TCHANNEL_BASE'] + '/input/sampleJSONs/jetConstTrainingNtuples/' + kind + "/" + year + "/"
    else:
        JSONDir = environ['TCHANNEL_BASE'] + '/input/sampleJSONs/' + kind + "/" + year + "/"
    inputSamples = glob.glob(JSONDir+"*"+detailKey+"*.json")
    if len(inputSamples) == 0:
        print("Error: no json file found with name:", JSONDir)

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

def sfGetter(sample):
    scaleFactor = {
        "2018_mQCDmini_Pt_80to120":      29535000./3565000., # 29535000./575000.
        "2018_mQCDmini_Pt_120to170":     25255000./3940000.,
        "2018_mQCDmini_Pt_170to300":     29710000./3175000.,
        "2018_mQCDmini_Pt_300to470":     41744000./3600000.,
        "2018_mQCDmini_Pt_470to600":     42459973./3712000., # only the second half of the mini samples
        "2018_mQCDmini_Pt_600to800":     64061000./1848000.,
        "2018_mQCDmini_Pt_800to1000":    37598000./1902000.,
        "2018_mQCDmini_Pt_1000to1400":   18485000./2037000.,
        "2018_mQCDmini_Pt_1400to1800":   6928000./1972000.,
        "2018_mQCDmini_Pt_1800to2400":   4017800./1445800.,
        "2018_mQCDmini_Pt_2400to3200":   2394000./1440000.,
        "2018_mQCDmini_Pt_3200toInf":    800000./800000.,
        "2018_mTTJetsmini_Incl":         10244307./841560.
    }
    sf = 1
    for key,item in scaleFactor.items():
        if key == sample:
            sf = item
            break
    return sf



