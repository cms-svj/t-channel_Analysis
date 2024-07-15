from utils import samples as s
import os
import sys
import numpy as np

nFilesPerJob = 5 # usually 5
maxFilesPerSample = -1 # use -1 to run over all the samples
maxJobs = 50
submissionMode = 1 # 0 = submit jobs without mixing different samples together (the old way), 1 = submit jobs while mixing samples together (lower the number of submissions)
skimSource = False
outHistF = "output/QCD_cl0p02_net64"
eTagName = "sdt_QCD_disco_0p001_closure_0p02_damp_1_net_64_32_16_8_1Evt_pn"
listOfSamplesToRun = [
                        "2018_QCD",
                        "2018_ST",
                        "2018_TTJets",
                        "2018_WJets",
                        "2018_ZJets",
                        "2018_Data",
                        "2018_SVJ_t",
]
runSignalLocal = False # process signal samples locally, more reliable, but slow
evtTaggerLoc = f"utils/data/DNNEventClassifier/{eTagName}"
preCommand = "python analyze.py"
# adding extra flags
if maxFilesPerSample != -1: # if not all the files are used, then turn the scale factor flag on
    preCommand += " -f" 
if skimSource:
    preCommand += " --skimSource" 
maxFiles = nFilesPerJob*maxJobs

for sampleToRun in listOfSamplesToRun:
    reRunDetail = s.getSamplesFromGroup(sampleToRun,skimSource=skimSource)
    # figuring out how many files each subsample has for the given constraints
    jobDetails = []
    for sample in reRunDetail:
        fileset = s.getFileset(sample, False)
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

    # since signal samples are much smaller than the background and data samples, we can always process all of them
    if "SVJ" in sampleToRun:
        for job in jobDetails:
            sample, mVal, nVal = job
            if runSignalLocal:
                command = f"{preCommand} -d {sample} -w 4 --outHistF {outHistF} -t {evtTaggerLoc} -j -s 1000"
            else:
                command = f"{preCommand} -d {sample} -b 20 --outHistF {outHistF} -t {evtTaggerLoc} -j -s 1000 --condor --dask"
            print(command)
            os.system(command)
    else:
        if (submissionMode == 0):
            for job in jobDetails:
                sample, mVal, nVal = job
                command = f"{preCommand} -d {sample} -N {nVal} -M {mVal} -b {maxJobs} --outHistF {outHistF} -t {evtTaggerLoc} -j -s 1000 --condor --dask"
                print(command)
                os.system(command)
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
                print(command)
                os.system(command)