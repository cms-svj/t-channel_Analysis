from utils import samples as s
from utils.python import lpc_dask as ld
import os
import sys
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--nFilesPerJob",       type=int, default=5, help="Number of files to run per job. Higher number can cause memory issue.")
parser.add_argument("--maxFilesPerSample",  type=int, default=-1, help="The number of files to run over per sample group. -1 = all available files.")
parser.add_argument("--maxJobs",            type=int, default=50, help="Maximum number of jobs to run on condor at a time. If inferencing on particleNet (not the case for skims), higher number than 100 can cause instability.")
parser.add_argument("--submissionMode",     type=int, default=1, help="0 = submit jobs without mixing different samples together, 1 = submit jobs while mixing samples together (faster since lower number of submissions)")
parser.add_argument("--outHistF",           type=str, default="output/QCD_cl0p02_net64", help="Location to save all the output histograms")
parser.add_argument("--eTagName",           type=str, default="sdt_QCD_disco_0p001_closure_0p02_damp_1_net_64_32_16_8_1Evt_pn", help="Name of the event classifier.")
parser.add_argument("--skimSource",         action="store_true", help="Use skims instead of TreeMaker ntuples. (under development)")
parser.add_argument("--runSignalLocal",     action="store_true", help="Run over signals locally. Slow, but more stable than running on condor.")
parser.add_argument("--printOnly",          action="store_true", help="Print the commands without running them.")
parser.add_argument("--rerunMissingFiles",  action="store_true", help="Rerun the jobs based on missing files.")
parser.add_argument("--haddAll",            action="store_true", help="Hadd all the output files by their sample group. (This has to be run outside of the coffeaenv. Try running in `source /cvmfs/sft.cern.ch/lcg/views/LCG_106_cuda/x86_64-el9-gcc11-opt/setup.sh`)")
parser.add_argument('--runJetTag',          action='store_true', help='Run jet tagger.', )
parser.add_argument('--runEvtClass',        action='store_true', help='Run event classifier.')
parser.add_argument("--skimCut",            type=str, default="t_channel_pre_selection", help='The selection of cuts that have been applied to the TM ntuples when making the skims: t_channel_pre_selection, t_channel_lost_lepton_control_region, etc.')

args = parser.parse_args()

nFilesPerJob = args.nFilesPerJob 
maxFilesPerSample = args.maxFilesPerSample 
maxJobs = args.maxJobs
submissionMode = args.submissionMode 
skimSource = args.skimSource
printOnly = args.printOnly 
rerunMissingFiles = args.rerunMissingFiles
outHistF = args.outHistF
eTagName = args.eTagName
haddAll = args.haddAll
listOfSampleGroupsToRun = [
                        "2018_QCD",
                        "2018_ST",
                        "2018_TTJets",
                        "2018_WJets",
                        "2018_ZJets",
                        "2018_Data",
                        "2018_SVJ_t",
]
runSignalLocal = args.runSignalLocal
evtTaggerLoc = f"utils/data/DNNEventClassifier/{eTagName}"
preCommand = "python analyze.py"
# adding extra flags
if maxFilesPerSample != -1: # if not all the files are used, then turn the scale factor flag on
    preCommand += " -f" 
if skimSource:
    preCommand += f" --skimSource --skimCut {args.skimCut}" 
if args.runJetTag:
    preCommand += " --runJetTag" 
if args.runEvtClass:
    preCommand += " --runEvtClass" 
maxFiles = nFilesPerJob*maxJobs
expectedFilesDict = {}

def runOrPrintCommand(command,haddAll=False,printOnly=False):
    if printOnly:
        print(command)
    elif haddAll:
        pass
    else:
        print(command)
        os.system(command)

def runMissingFile(command,outHistF,dataset,nFiles,startFile,condor,dask,hemPeriod):
    outFile = ld.out_file_name_creator(outHistF,dataset.split(" "),nFiles.split(" "),startFile.split(" "),condor,dask,hemPeriod).replace(f"{outHistF}/","")
    if outFile not in os.listdir(outHistF):
        print("Rerunning:")
        runOrPrintCommand(command,haddAll,printOnly)

def addExpectedFile(outHistF,dataset,nFiles,startFile,condor,dask,hemPeriod,expectedFilesDict,sampleGroupToRun):
    outFile = ld.out_file_name_creator(outHistF,dataset.split(" "),nFiles.split(" "),startFile.split(" "),condor,dask,hemPeriod).replace(f"{outHistF}/","")
    if "SVJ" in sampleGroupToRun:
        expectedFilesDict[sampleGroupToRun].append([outFile,dataset])
    else:    
        expectedFilesDict[sampleGroupToRun].append(outFile)

def runHadd(expectedFilesDict,outHistF):
    for sampleGroup, sampleList in expectedFilesDict.items():
        if "SVJ" in sampleGroup:
            for sample in sampleList:
                if sample[0] in os.listdir(outHistF):
                    print(f"mv {outHistF}/{sample[0]} {outHistF}/{sample[1]}.root")
                else:
                    print(f"{sample[0]} is missing.")
                    os.system(f"mv {sample[0]} {sample[1]}.root")
        else:
            command = f"hadd -f {sampleGroup}.root "
            missingFile = False
            for sample in sampleList:
                if sample in os.listdir(outHistF):
                    command += f"{sample} "
                else:
                    missingFile = True
                    break
            if missingFile:
                print(f"{sampleGroup} has missing files. hadd was not performed.")
            else:
                print(command)
                os.system(command)

for sampleGroupToRun in listOfSampleGroupsToRun:
    expectedFilesDict[sampleGroupToRun] = []
    runDetail = s.getSamplesFromGroup(sampleGroupToRun,args.skimCut,skimSource=skimSource)
    # figuring out how many files each subsample has for the given constraints
    jobDetails = []
    # since signal samples are much smaller than the background and data samples, we can always process all of them
    if "SVJ" in sampleGroupToRun:
        for sample in runDetail:
            jobDetails.append([sample,"0","-1"])
        for job in jobDetails:
            sample, mVal, nVal = job
            # since the signal samples are smaller, it is possible to process them offline, or process them with less workers on condor.
            if runSignalLocal:
                command = f"{preCommand} -d {sample} -N {nVal} -M {mVal} -w 4 --outHistF {outHistF} -t {evtTaggerLoc} -j -s 1000"
                addExpectedFile(outHistF,sample,nVal,mVal,False,False,"",expectedFilesDict,sampleGroupToRun)
                if rerunMissingFiles:
                    runMissingFile(command,outHistF,sample,nVal,mVal,False,False,"")
            else:
                command = f"{preCommand} -d {sample} -N {nVal} -M {mVal} -b 20 --outHistF {outHistF} -t {evtTaggerLoc} -j -s 1000 --condor --dask"
                addExpectedFile(outHistF,sample,nVal,mVal,True,True,"",expectedFilesDict,sampleGroupToRun)
                if rerunMissingFiles:
                    runMissingFile(command,outHistF,sample,nVal,mVal,True,True,"")
            if rerunMissingFiles == False:
                runOrPrintCommand(command,haddAll,printOnly)
    else:
        for sample in runDetail:
            fileset = s.getFileset(sample,skimCut=args.skimCut,startFile=0,nFiles=-1,skimSource=args.skimSource,verbose=False)
            if maxFilesPerSample == -1:
                totalNumberOfFiles = len(list(fileset.items())[0][1]) # run over all the MC samples
            else:
                totalNumberOfFiles = maxFilesPerSample
            # totalNumberOfFiles = 20  # run over 20 of files
            numberOfJobsRequired = int(np.ceil(totalNumberOfFiles/nFilesPerJob)) 
            remainingFiles = totalNumberOfFiles
            if numberOfJobsRequired > maxJobs:
                for mVal in np.arange(0,totalNumberOfFiles,maxFiles):
                    if numberOfJobsRequired >= maxJobs:
                        jobDetails.append([sample,mVal,maxFiles])
                        numberOfJobsRequired -= maxJobs
                        remainingFiles -= maxFiles
                    else:
                        jobDetails.append([sample,mVal,remainingFiles])
            else:
                jobDetails.append([sample,0,totalNumberOfFiles])
        if (submissionMode == 0):
            for job in jobDetails:
                sample, mVal, nVal = job
                mVal = str(mVal)
                nVal = str(nVal)
                command = f"{preCommand} -d {sample} -N {nVal} -M {mVal} -b {maxJobs} --outHistF {outHistF} -t {evtTaggerLoc} -j -s 1000 --condor --dask"
                addExpectedFile(outHistF,sample,nVal,mVal,True,True,"",expectedFilesDict,sampleGroupToRun)
                if rerunMissingFiles:
                    runMissingFile(command,outHistF,sample,nVal,mVal,True,True,"")
                else:
                    runOrPrintCommand(command,haddAll,printOnly)
        else:
            # group different runs together such that the total number of files are always maxFiles
            sampleGroupList = []
            mValGroupList = []
            nValGroupList = []
            sampleGroup = []
            mValGroup = []
            nValGroup = []
            numberOfFilesInGroup = 0
            for i in range(len(jobDetails)):
                job = jobDetails[i]
                sample, mVal, nVal = job
                sampleGroup.append(sample)
                mValGroup.append(mVal)
                nValGroup.append(nVal)
                numberOfFilesInGroup += nVal
                if numberOfFilesInGroup == maxFiles:
                    sampleGroupList.append(sampleGroup)
                    mValGroupList.append(mValGroup)
                    nValGroupList.append(nValGroup)
                    sampleGroup = []
                    mValGroup = []
                    nValGroup = []
                    numberOfFilesInGroup = 0
                elif numberOfFilesInGroup > maxFiles:
                    sampleGroupList.append(sampleGroup)
                    mValGroupList.append(mValGroup)
                    allowedNumOfFiles = maxFiles - (numberOfFilesInGroup - nVal)
                    nValGroup[-1] = allowedNumOfFiles
                    nValGroupList.append(nValGroup)
                    sampleGroup = [sample]
                    mValGroup = [mVal + allowedNumOfFiles]
                    nValGroup = [nVal - allowedNumOfFiles]
                    numberOfFilesInGroup = nVal - allowedNumOfFiles
                if i == len(jobDetails)-1:
                    sampleGroupList.append(sampleGroup)
                    mValGroupList.append(mValGroup)
                    nValGroupList.append(nValGroup)

            for i in range(len(sampleGroupList)):
                sampleGroup = sampleGroupList[i]
                mValGroup = mValGroupList[i]
                nValGroup = nValGroupList[i]
            if not haddAll:
                print("Number of batches of submission without combining:",len(jobDetails))
                print("Number of batches of submission with combining:",len(sampleGroupList))
            for i in range(len(sampleGroupList)):
                sampleGroup = sampleGroupList[i]
                mValGroup = mValGroupList[i]
                nValGroup = nValGroupList[i]
                sampleList = ""
                mValList = ""
                nValList = ""
                for i in range(len(sampleGroup)):
                    sampleList += f"{sampleGroup[i]} "
                    mValList += f"{mValGroup[i]} "
                    nValList += f"{nValGroup[i]} "
                command = f"{preCommand} -d {sampleList}-N {nValList}-M {mValList}-b {maxJobs} --outHistF {outHistF} -t {evtTaggerLoc} -j -s 1000 --condor --dask"
                addExpectedFile(outHistF,sampleList[:-1],nValList[:-1],mValList[:-1],True,True,"",expectedFilesDict,sampleGroupToRun)
                if rerunMissingFiles:
                    runMissingFile(command,outHistF,sampleList[:-1],nValList[:-1],mValList[:-1],True,True,"")
                else:
                    runOrPrintCommand(command,haddAll,printOnly)
if haddAll:
    runHadd(expectedFilesDict,outHistF)
