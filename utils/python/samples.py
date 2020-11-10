from os import system, environ
import json
import glob

class Samples:

    def getFileset(self,sample,verbose=True,startFile=0,nFiles=-1):
        # find all json files for the sample or sample collection
        f_ = sample.find("_")
        year = sample[:f_]    
        detailKey = sample[f_:]
        kind = "signals" if ("mMed" in detailKey or "mZprime" in detailKey) else "backgrounds"
        
        if "Incl" in detailKey:
            ii = detailKey.find("Incl")
            detailKey = detailKey[:ii] + "Tune"

        JSONDir = environ['TCHANNEL_BASE'] + '/input/sampleJSONs/' + kind + "/" + year + "/"
        allfiles = glob.glob(JSONDir+"*.json")
        if len(allfiles) == 0:
            print("Error: no json file found with name:", JSONDir)
        inputSamples = list(f for f in allfiles if detailKey in f)

        # open all json files and dump them into a dictionary
        fileset = {}
        for s in inputSamples: 
            fileset.update(json.load(open(s ,'r')))    

        # condor specific code to process a subset of a sample
        if len(fileset.keys()) == 1:
            fs = {}
            for n, rFiles in fileset.items():
                fs[n] = []
                if nFiles < 0:
                    nFiles = len(rFiles)
                fn = startFile
                while fn < startFile+nFiles and fn < len(rFiles):
                    #print(fn, len(rFiles))
                    fs[n].append(rFiles[fn])
                    fn+=1                    
            fileset = fs

        # print the sample names if verbose
        if verbose:
            for key, _ in fileset.items(): 
                print(key)

        return fileset
    
    def getAllFilesets(self):
        fBG = self.getFileset("*_", False)
        fSG1 = self.getFileset("*_mMed", False)
        fSG2 = self.getFileset("*_mZprime", False)
        fBG.update(fSG1)
        return fBG.update(fSG2)

