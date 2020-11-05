import sys, os
from os import system, environ
import subprocess

from samples import SampleCollection
import optparse
from glob import glob
import datetime
import shutil
import ROOT

def red(string):
     CRED = "\033[91m"
     CEND = "\033[0m"
     return CRED + str(string) + CEND

def checkNumEvents(nEvents, rootFile, sampleCollection, log):
    try:
         f = ROOT.TFile.Open(rootFile)
         f.cd()
         try:
              h = f.Get("EventCounter")
              nNeg = h.GetBinContent(1)
              nPos = h.GetBinContent(2)
              diff = nEvents-(nPos-nNeg)
              if abs(diff) > 5.0:
                   message = "Error: Sample: "+sampleCollection+" Expected nEvents:  "+str(nEvents)+" EventCounter nEvents: "+str(nPos-nNeg)+" = "+str(nPos)+" "+str(-nNeg)
                   log.append(message)
                   print red("----------------------------------------------------------------------------------------------------------")
                   print red("Num events in \"EventCounter\" doesn't match the number in \"sampleSet.cfg\"")
                   print red(message)
                   print red("----------------------------------------------------------------------------------------------------------")
         except:
              print red("Error: Problem opening and reading from histogram \"EventCounter\"")
              pass
         f.Close()
    except:
         print red("Error: Can't open rootFile: %s" % rootFile)
         pass

    return log

def getDataSets(inPath):
    l = glob(inPath+"/*")
    print "-----------------------------------------------------------------------------" 
    print red("Warning: No dataset specified: using all directory names in input path")
    print "-----------------------------------------------------------------------------\n" 
    return list(s[len(inPath)+1:] for s in l)

def main():
    # Parse command line arguments
    parser = optparse.OptionParser("usage: %prog [options]\n")
    parser.add_option('-d', dest='datasets', type='string', default='',             help="Lists of datasets, comma separated")
    parser.add_option('-H', dest='outDir',   type='string', default='rootfiles',    help="Can pass in the output directory name")
    parser.add_option('-p', dest='inPath',   type='string', default='output-files', help="Can pass in the input directory name")
    parser.add_option('-y', dest='year',     type='string', default='',             help="Can pass in the year for this data")
    parser.add_option('-o', action='store_true',                                    help="Overwrite output directory")
    parser.add_option('--noHadd', action='store_true',                              help="Dont hadd the the root files")
    parser.add_option('--haddOther', action='store_true',                           help="Do the hack to make BG_OTHER.root")
    parser.add_option('--haddData', action='store_true',                            help="Do the hack to make Data.root")
    options, args = parser.parse_args()

    # Get input directory path
    inPath = options.inPath
        
    # Checks if user specified a dataset(s)
    datasets = []
    if options.datasets:
        datasets = options.datasets.split(',')
    else:
        datasets = getDataSets(inPath)
    
    # Check if output directory exits and makes it if not
    outDir = options.outDir
    overwrite = options.o
    if os.path.exists(outDir):
        if overwrite: 
            print red("Warning: Overwriting output directory")
            shutil.rmtree(outDir)
            os.makedirs(outDir)
        else:
            print red("Error: Output directory %s already exits" % ('"'+outDir+'"'))
            exit(0)    
    else:
        os.makedirs(outDir) 
    
    # Loop over all sample options to find files to hadd
    log = []
    sc = SampleCollection("../sampleSets.cfg", "../sampleCollections.cfg")
    scl = sc.sampleCollectionList()
    for sampleCollection in scl:
        sl = sc.sampleList(sampleCollection)
        if sampleCollection in datasets:
            directory = sampleCollection
            files = ""
            print "-----------------------------------------------------------"
            print sampleCollection
            print "-----------------------------------------------------------"
            
            # hadd signal root files
            sampleSetsToHadd = ["2016_AllSignal", "2017_AllSignal", "2017_AllSignal_CP5", "2018pre_AllSignal", 
                                "2018post_AllSignal", "2018_t-channel_AllSignal", "2018_s-channel_AllSignal"]
            if sampleCollection in sampleSetsToHadd:
                for sample in sl:
                    files = " " + " ".join(glob("%s/%s/MyAnalysis_%s_*.root" % (inPath, directory, sample[1])))
                    outfile = "%s/%s.root" % (outDir,sample[1])
                    command = "hadd %s/%s.root %s" % (outDir, sample[1], files)
                    print command
                    if not options.noHadd: system(command)
                    log = checkNumEvents(nEvents=float(sample[2]), rootFile=outfile, sampleCollection=sample[1], log=log)
    
            # hadd other condor jobs
            else:
                nEvents=0.0
                for sample in sl:
                    files += " " + " ".join(glob("%s/%s/MyAnalysis_%s_*.root" % (inPath, directory, sample[1])))
                    nEvents+=float(sample[2])
    
                outfile = "%s/%s.root" % (outDir,sampleCollection)
                command = "hadd %s %s" % (outfile, files)
                try:
                    if not options.noHadd: 
                        process = subprocess.Popen(command, shell=True)
                        process.wait()
                except:
                    print red("Warning: Too many files to hadd, using the exception setup")
                    command = "hadd %s/%s.root %s/%s/*" % (outDir, sampleCollection, inPath, sampleCollection)
                    if not options.noHadd: system(command)
                    pass
    
                log = checkNumEvents(nEvents=nEvents, rootFile=outfile, sampleCollection=sampleCollection, log=log)
    
    #Print log of hadd at the end
    if len(log) > 0:
         print red("------------------------------------------------------------------------------------------------")
         print red("There was some jobs that didn't match the epected number of events, see summary below")
         for l in log:
              print red(l)
         print red("------------------------------------------------------------------------------------------------")

    if options.haddOther:
        # Hack to make the BG_OTHER.root file
        sigNttbar_old = ["AllSignal", "TT", "TTJets", "Data_SingleMuon", "Data_SingleElectron"]
        sigNttbar_2016 = ["2016_AllSignal", "2016_TT", "2016_TTJets", "2016_Data_SingleMuon", "2016_Data_SingleElectron","2016_TT_isrUp", "2016_TT_isrDown", "2016_TT_fsrUp", "2016_TT_fsrDown", "2016_TTX" , "2016_QCD", "2016_TT_erdOn", "2016_TT_hdampUp", "2016_TT_hdampDown", "2016_TT_underlyingEvtUp", "2016_TT_underlyingEvtDown"]
        sigNttbar_2017 = ["2017_AllSignal", "2017_AllSignal_CP5", "2017_TT", "2017_TTJets", "2017_Data_SingleMuon", "2017_Data_SingleElectron", "2017_TTX", "2017_QCD", "2017_TT_erdOn", "2017_TT_hdampUp", "2017_TT_hdampDown", "2017_TT_underlyingEvtUp", "2017_TT_underlyingEvtDown"]
        sigNttbar_2018pre  = ["2018pre_AllSignal", "2018pre_TT", "2018pre_TTJets", "2018pre_Data_SingleMuon", "2018pre_Data_SingleElectron", "2018pre_TTX", "2018pre_QCD", "2018pre_TT_erdOn", "2018pre_TT_hdampUp", "2018pre_TT_hdampDown", "2018pre_TT_underlyingEvtUp", "2018pre_TT_underlyingEvtDown"]
        sigNttbar_2018post = ["2018post_AllSignal", "2018post_TT", "2018post_TTJets", "2018post_Data_SingleMuon", "2018post_Data_SingleElectron", "2018post_TTX", "2018post_QCD",  "2018post_TT_erdOn", "2018post_TT_hdampUp", "2018post_TT_hdampDown", "2018post_TT_underlyingEvtUp", "2018post_TT_underlyingEvtDown"]
        sigNttbar = sigNttbar_old+sigNttbar_2016+sigNttbar_2017+sigNttbar_2018pre+sigNttbar_2018post
        files = ""
        for sampleCollection in scl:
            sl = sc.sampleList(sampleCollection)
            if sampleCollection in datasets:
                if sampleCollection not in sigNttbar: 
                    directory = sampleCollection
                    files += " %s/%s.root " % (outDir, directory)
        if options.year:
            command = "hadd %s/%s_BG_OTHER.root %s" % (outDir, options.year, files)
        else:
            command = "hadd %s/BG_OTHER.root %s" % (outDir, files)
        print "-----------------------------------------------------------"
        print command
        print "-----------------------------------------------------------"
        system(command)
        
    if options.haddData:
        # Hack to make the Data.root file (hadd all the data together)
        dataFiles = ["Data_SingleMuon.root", "Data_SingleElectron.root", 
                     "2016_Data_SingleMuon.root", "2016_Data_SingleElectron.root", 
                     "2017_Data_SingleMuon.root", "2017_Data_SingleElectron.root",
                     "2018pre_Data_SingleMuon.root", "2018pre_Data_SingleElectron.root",
                     "2018post_Data_SingleMuon.root", "2018post_Data_SingleElectron.root"]
        if options.year:
            command = "hadd %s/%s_Data.root " % (outDir,options.year)
        else:
            command = "hadd %s/Data.root " % outDir
        for f in dataFiles:
            if os.path.exists(outDir+"/"+f):
                command += " %s/%s" % (outDir, f)
        print "-----------------------------------------------------------"
        print command
        print "-----------------------------------------------------------"
        system(command)

if __name__ == "__main__":
    main()
