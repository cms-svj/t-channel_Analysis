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

    JSONDir = environ['TCHANNEL_BASE'] + '/input/sampleJSONs/' + kind + "/" + year + "/"

    # inputSamples = glob.glob(JSONDir+"*"+detailKey+"*.json")
    inputSamples = glob.glob(JSONDir+year+"_"+detailKey+"*.json")
    if len(inputSamples) == 0:
        print("Error: no json file found with name:", JSONDir)
    else:
        print(inputSamples)

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

def sfGetter(sample,on=False):
    # scaleFactor for particleNet training files
    if on:
        print("Yes")
        scaleFactor = {
                "2018_QCD_Pt_1000to1400": 19730000.0/48000,
                "2018_QCD_Pt_1400to1800": 10982000.0/21000,
                "2018_QCD_Pt_1800to2400": 5491000.0/48000,
                "2018_QCD_Pt_2400to3200": 2997000.0/30000,
                "2018_QCD_Pt_300to470": 57910000.0/9423000,
                "2018_QCD_Pt_3200toInf": 1000000.0/51000,
                "2018_QCD_Pt_470to600": 52448000.0/447000,
                "2018_QCD_Pt_600to800": 67508000.0/96000,
                "2018_QCD_Pt_800to1000": 37160000.0/48000,
                "2018_TTJets_DiLept": 29290487.0/8532037,
                "2018_TTJets_DiLept_genMET-150": 10592111.0/103384,
                "2018_TTJets_HT-1200to2500": 2009331.0/29630,
                "2018_TTJets_HT-2500toInf": 1001084.0/78961,
                "2018_TTJets_HT-600to800": 15258099.0/83062,
                "2018_TTJets_HT-800to1200": 9201990.0/69470,
                "2018_TTJets_SingleLeptFromT": 58237254.0/2572246,
                "2018_TTJets_SingleLeptFromT_genMET-150": 13337428.0/414361,
                "2018_TTJets_SingleLeptFromTbar": 58510607.0/3589663,
                "2018_TTJets_SingleLeptFromTbar_genMET-150": 12597052.0/503985,
                "2018_TTJets_Incl": 10494353.0/10339641,
        }
        sf = 1
        for key,item in scaleFactor.items():
            if key == sample:
                print("Yes")
                sf = item
                break
    else:
        sf = 1
    return sf



