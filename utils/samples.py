from os import system, environ
import json
import glob

def getFileset(sample,verbose=True,startFile=0,nFiles=-1):
    # find all json files for the sample or sample collection
    f_ = sample.find("_")
    year = sample[:f_]
    detailKey = sample[f_+1:]
    kind = "signals" if ("mMed" in detailKey or "mZprime" in detailKey) else "backgrounds"
    if "Incl" in detailKey:
        ii = detailKey.find("Incl")
        detailKey = detailKey[:ii] + "Tune"

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
        fileset[n] = []
        if nFiles < 0:
            nFiles = len(rFiles)
        fn = startFile
        while fn < startFile+nFiles and fn < len(rFiles):
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

