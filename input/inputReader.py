import json


def grabInput(sample):
    baseloc = "root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV17/"
    sigloc = "root://cmseos.fnal.gov//store/user/keanet/tchannel/SVJP_08272020/NTuples/"

    f_ = sample.find("_")
    year = sample[:f_]


    if year == "2016":
        yrkey = "Summer16v3"
    elif year == "2017":
        yrkey = "Fall17"
    elif year == "2018":
        yrkey = "Autumn18"

    detailKey = sample[f_+1:]

    if "mMed" in detailKey:
        readFrom = "tchannel_Files.txt"
        loc = sigloc
        yrkey = ""
    else:
        readFrom = "Run2ProductionV17_Files.txt"
        loc = baseloc

    if "Incl" in detailKey:
        ii = detailKey.find("Incl")
        detailKey = detailKey[:ii] + "Tune"

    fileReadFrom = open(readFrom,'r')
    Lines = fileReadFrom.readlines()

    fileList = []

    for line in Lines:
        if yrkey in line and detailKey in line:
            fileList.append(loc+line[:-1])
    return fileList

sampleLabels = open('sampleLabels.txt','r')
sampleLists = sampleLabels.readlines()

for sample in sampleLists:
    f_ = sample.find("_")
    year = sample[:f_]
    detailKey = sample[f_+1:-1]

    if "Incl" in detailKey:
        ii = detailKey.find("Incl")
        detailKey = detailKey[:ii] + "Tune"

    slist = {sample[:-1]:grabInput(sample[:-1])}

    firstSample = slist[sample[:-1]][0]
    si = firstSample.find("_0_")
    fi = firstSample.find(detailKey)

    if "Tune" in detailKey:
        sampleLab = year + "_" + detailKey + firstSample[fi+len(detailKey):si]
    else:
        sampleLab = year + "_" + detailKey + "_" + firstSample[fi+len(detailKey)+1:si]

    if "_ext" in sampleLab:
        sampleLab = sampleLab[:sampleLab.find("_ext")]

    print (sampleLab)
    if "mMed" in detailKey:
        saveFile = "sampleJSONs/signals/" + year + "/" + sampleLab + ".json"
    else:
        saveFile = "sampleJSONs/backgrounds/" + year + "/" + sampleLab + ".json"

    with open(saveFile,"w") as fp:
        json.dump(slist, fp, indent=4)
