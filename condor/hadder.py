import sys, os
from os import system, environ
sys.path = [environ["TCHANNEL_BASE"],] + sys.path
from utils.python.samples import Samples as s
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
                   print(red("----------------------------------------------------------------------------------------------------------"))
                   print(red("Num events in \"EventCounter\" doesn't match the number in \"sampleSet.cfg\""))
                   print(red(message))
                   print(red("----------------------------------------------------------------------------------------------------------"))
         except:
              print(red("Error: Problem opening and reading from histogram \"EventCounter\""))
              pass
         f.Close()
    except:
         print(red("Error: Can't open rootFile: %s" % rootFile))
         pass

    return log

def getDataSets(inPath):
    l = glob(inPath+"/*")
    print("-----------------------------------------------------------------------------")
    print(red("Warning: No dataset specified: using all directory names in input path"))
    print("-----------------------------------------------------------------------------\n")
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
            print(red("Warning: Overwriting output directory"))
            shutil.rmtree(outDir)
            os.makedirs(outDir)
        else:
            print(red("Error: Output directory %s already exits" % ('"'+outDir+'"')))
            exit(0)
    else:
        os.makedirs(outDir)

    # Loop over all sample options to find files to hadd
    log = []
    scl = s().getAllFilesets()
    for sampleCollection in datasets:
        sl = s().getFileset(sampleCollection, False)
        directory = sampleCollection
        files = ""
        print("-----------------------------------------------------------")
        print(sampleCollection)
        print("-----------------------------------------------------------")

        # hadd signal root files
        sampleSetsToHadd = ["2016_mMed", "2017_mMed", "2018_mMed"]
        if sampleCollection in sampleSetsToHadd:
            for sample in sl.keys():
                files = " " + " ".join(glob("%s/%s/MyAnalysis_%s_*.root" % (inPath, directory, sample)))
                outfile = "%s/%s.root" % (outDir,sample)
                command = "hadd %s/%s.root %s" % (outDir, sample, files)
                if not options.noHadd: system(command)
                #log = checkNumEvents(nEvents=float(sample[2]), rootFile=outfile, sampleCollection=sample[1], log=log)
        # hadd other condor jobs
        else:
            nEvents=0.0
            for sample in sl.keys():
                files += " " + " ".join(glob("%s/%s/MyAnalysis_%s_*.root" % (inPath, directory, sample)))
                #nEvents+=float(sample[2])

            outfile = "%s/%s.root" % (outDir,sampleCollection)
            command = "hadd %s %s" % (outfile, files)
            try:
                if not options.noHadd:
                    process = subprocess.Popen(command, shell=True)
                    process.wait()
            except:
                print(red("Warning: Too many files to hadd, using the exception setup"))
                command = "hadd %s/%s.root %s/%s/*" % (outDir, sampleCollection, inPath, sampleCollection)
                if not options.noHadd: system(command)
                pass

            #log = checkNumEvents(nEvents=nEvents, rootFile=outfile, sampleCollection=sampleCollection, log=log)

    #Print log of hadd at the end
    if len(log) > 0:
         print(red("------------------------------------------------------------------------------------------------"))
         print(red("There was some jobs that didn't match the epected number of events, see summary below"))
         for l in log:
              print(red(l))
         print(red("------------------------------------------------------------------------------------------------"))


if __name__ == "__main__":
    main()
