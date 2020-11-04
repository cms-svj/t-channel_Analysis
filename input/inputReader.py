def grabInput(bkgLab):
    # bkgLab = sys.argv[1] # background label, e.g. QCD16_Pt_80to120
    years = ["16","17","18"]
    yrV17s = ["Summer16v3","Fall17","Autumn18"]
    bkgs = ["QCD","TTJets","WJets","ZJets"]

    baseloc = "root://cmseos.fnal.gov//store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV17/"
    sigloc = "root://cmseos.fnal.gov//store/user/keanet/tchannel/SVJP_08272020/NTuples/"
    dict_ = {}

    fui = bkgLab.find("_") # index of the first underscore
    comKey = [] # complete key string for identifying the correct files

    if fui != -1: # background label identification; signal labels do not contain underscore, so .find would have returned -1.
        # inFile = "Run2ProductionV17_Files.txt" # use this for test in the input directory
        inFile = "input/Run2ProductionV17_Files.txt" # use this for running test.py
        yr_key = 0
        bkg_key = ""

        # first identify the year from the label
        for iy in range(len(years)):
            year = years[iy]
            if year in bkgLab[:fui]:
                yr_key = yrV17s[iy]

        # identify the background
        for bkg in bkgs:
            if bkg in bkgLab[:fui]:
                if bkg == "WJets":
                    bkg_key = "WJetsToLNu"
                elif bkg == "ZJets":
                    bkg_key = "ZJetsToNuNu"
                else:
                    bkg_key = bkg

        # identify the kind of background
        bkinds = ["Pt","HT","SingleLeptFromT_genMET-150",
        "SingleLeptFromTbar_genMET-150","DiLept_genMET-150",
        "SingleLeptFromTbar","SingleLeptFromT","DiLept"]

        bkind = ""

        for bk in bkinds:
            if bk in bkgLab[fui+1:]:
                bkind = bk
                # the following condition is to avoid _genMET version of the background gets included
                if bk == "SingleLeptFromTbar" or bk == "SingleLeptFromT" or bk == "DiLept":
                    bkind = bk + "_Tune"
                break
            else:
                bkind = "Tune"

        ran = ""
        if bkind == "Pt" or bkind == "HT":
            sui = bkgLab[fui+1:].find("_") + fui
            ran = bkgLab[sui+2:]

        lastSep = "_"
        if bkind == "HT":
            lastSep = "-"
        elif "Tune" in bkind:
            lastSep = ""

        # the complete key for identifying the ntuple files
        comKey.append(yr_key + "." + bkg_key + "_" + bkind + lastSep + ran)

    else: # signal label identification
        # inFile = "tchannel_Files.txt" # use this for test in the input directory
        inFile = "input/tchannel_Files.txt" # use this for running test.py
        pAlias = ["M","d","r","a"]
        pCom = ["mMed","mDark","rinv","alpha"]

        pvar = ""

        # figure out parameter that is different from baseline
        if bkgLab == "base":
            comKey.append("PrivateSamples.SVJ_2018_t-channel_mMed-3000_mDark-20_rinv-0p3_alpha-peak")
        else:
            for pi in range(len(pAlias)):
                pA = pAlias[pi]
                if pA in bkgLab:
                    pvar = pCom[pi]

        yi = bkgLab.find("y") # "y" is a key character for yukawa
        if yi != -1:
            pval = bkgLab[1:yi]
            comKey.append("yukawa-" + bkgLab[yi+1:] + "_")
        else:
            pval = bkgLab[1:] # if we don't specify it in the label, just assume yukawa-1
            comKey.append("yukawa-1_")

        if len(comKey) < 2:
            comKey.append(pvar + "-" + pval + "_")

    fileList = []

    file1 = open(inFile, 'r')
    Lines = file1.readlines()

    for line in Lines:
        cMet = 1
        for ck in comKey:
            if ck not in line:
                cMet = 0
        if cMet == 1:
            if fui != -1:
                fileList.append(baseloc+line[:-1])
            else:
                fileList.append(sigloc+line[:-1])

    return fileList
