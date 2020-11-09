#!/usr/bin/env python
from coffea import hist, processor
from processors.mainProcessor import MainProcessor  
import uproot
import sys
import json
import glob
import time
from optparse import OptionParser

def getFileset(sample):
    f_ = sample.find("_")
    year = sample[:f_]    
    detailKey = sample[f_:]
    kind = "signals" if "mMed" in detailKey else "backgrounds"
    
    if "Incl" in detailKey:
        ii = detailKey.find("Incl")
        detailKey = detailKey[:ii] + "Tune"
    
    JSONDir = 'input/sampleJSONs/' + kind + "/" + year + "/"
    allfiles = glob.glob(JSONDir+"*.json")
    inputSamples = list(f for f in allfiles if detailKey in f)

    fileset = {}
    for s in inputSamples: 
        fileset.update(json.load(open(s ,'r')))    
    for key, _ in fileset.items(): 
        print(key)
    return fileset

def main():
    # start run time clock
    tstart = time.time()

    # get options from command line
    parser = OptionParser()
    parser.add_option('-d', '--dataset',   help='dataset', dest='dataset')
    parser.add_option('-w', '--workers',   help='Number of workers to use for multi-worker executors (e.g. futures or condor)', dest='workers', type=int, default=8)
    options, args = parser.parse_args()
    
    # e.g. sample collection = "2016_QCD, 2016_WJets, 2016_TTJets, 2016_ZJets"
    sample = options.dataset
    
    # getting dictionary of files
    fileset = getFileset(sample)

    # run processor
    output = processor.run_uproot_job(
        fileset,
        treename='TreeMaker2/PreSelection',
        processor_instance=MainProcessor(),
        executor=processor.futures_executor,
        executor_args={'workers': options.workers, 'flatten': False},
        #chunksize=10000,
    )
        
    # export the histograms to root files
    hnames = list(output.keys())[:-1]
    fout = uproot.recreate('test.root')
    for h in hnames:
        fout[h] = hist.export1d(output[h])
    fout.close()

    # print run time in seconds
    dt = time.time() - tstart
    print("run time: %.2f [sec]" % (dt))

if __name__ == "__main__":
    main()
