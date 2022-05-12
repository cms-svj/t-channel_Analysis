import numpy as np
import awkward as ak
from . import triggerDict as tD

def TTStitch(dataset,events):
    # # TT Stiching mask
    ttStitchMask = None
    if "TTJets" in dataset:
        madHT = events.madHT
        nEle = ak.num(events.GenElectrons)
        nMu = ak.num(events.GenMuons)
        nTau = ak.num(events.GenTaus)
        genMET = events.GenMET
        if "TTJets_Inc" in dataset or "mTTJetsmini_Incl" in dataset:
            ttStitchMask = (madHT < 600) & (nEle==0) & (nMu==0) & (nTau==0)
        elif "TTJets_HT" in dataset:
            ttStitchMask = (madHT >= 600)
        elif ("TTJets_DiLept" in dataset and "genMET" not in dataset) or ("TTJets_SingleLeptFromT" in dataset and "genMET" not in dataset):
            ttStitchMask = (madHT < 600) & (genMET < 150)
        elif ("TTJets_DiLept" in dataset and "genMET" in dataset) or ("TTJets_SingleLeptFromT" in dataset and "genMET" in dataset):
            ttStitchMask = (madHT < 600) & (genMET >= 150)
    else:
        ttStitchMask = np.ones(len(events),dtype=bool)

    return ttStitchMask

def vetoPhiSpike(etaLead,phiLead,etaSub,phiSub,rad,eta,phi):
    veto = True
    if len(eta) >= 2:
        for iep in range(len(etaLead)):
            if (etaLead[iep] - eta[0])**2 + (phiLead[iep] - phi[0])**2 < rad:
                veto = False
                break
        for iep in range(len(etaSub)):
            if (etaSub[iep] - eta[1])**2 + (phiSub[iep] - phi[1])**2 < rad:
                veto = False
                break
    return veto

def PhiSpikeFilter(dataset,Jets):
    # phi spike filters
    rad = 0.028816 # half the length of the diagonal of the eta-phi rectangular cell
    rad *= 0.35 # the factor of 0.35 was optimized from the signal vs. background sensitivity study

    etalead = None
    etasub = None
    philead = None
    phisub = None
    if "16" in dataset:
        etalead = [0.048,0.24,1.488,1.584,-1.008]
        philead = [-0.35,-0.35,-0.77,-0.77,-1.61]
        etasub = [-1.2,-0.912,-0.912,-0.816,-0.72,-0.72,-0.528,-0.432,-0.336,-0.24,-0.24,-0.144,-0.144,-0.048,0.144,0.912,0.912,1.008,1.296,-1.584,-0.816,-0.72,-0.144,-0.048,-0.048,0.048,1.104,1.488]
        phisub = [-1.19,2.03,3.01,-1.75,-2.17,-0.77,2.73,2.73,0.21,0.07,0.21,-2.59,0.77,0.91,1.75,1.75,2.87,0.63,-0.49,0.63,1.47,-2.31,0.07,-2.59,0.77,0.91,-3.15,2.73]
    elif "17" in dataset:
        etalead = [0.144,1.488,1.488,1.584,-0.624]
        philead = [-0.35,-0.77,-0.63,-0.77,0.91]
        etasub = [-0.912,-0.912,-0.816,-0.72,-0.528,-0.336,-0.24,-0.24,-0.144,-0.144,-0.048,0.144,0.912,0.912,1.008,-1.2,-0.72,-0.72,-0.432,0.336,0.624,1.104,1.296]
        phisub = [2.03,3.01,-1.75,-0.77,2.73,0.21,0.07,0.21,-2.59,0.77,0.91,1.75,1.75,2.87,0.63,-1.19,-2.31,-2.17,2.73,-0.77,-0.77,-3.15,-0.49]
    elif "18" in dataset:
        etalead = [1.488,1.488,1.584]
        philead = [-0.77,-0.63,-0.77]
        etasub = [-1.584,-1.2,-0.912,-0.912,-0.816,-0.816,-0.72,-0.72,-0.528,-0.432,-0.336,-0.24,-0.24,-0.144,-0.144,-0.144,-0.048,-0.048,0.144,0.912,0.912,1.008,1.296,-0.72,1.104,1.488,1.776]
        phisub = [0.63,-1.19,2.03,3.01,-1.75,-0.77,-2.17,-0.77,2.73,2.73,0.21,0.07,0.21,-2.59,0.07,0.77,0.77,0.91,1.75,1.75,2.87,0.63,-0.49,-2.31,-3.15,-0.21,0.77]
    else:
        raise ValueError("No year label in dataset's name ")

    vfunc = np.vectorize(lambda eta,phi: vetoPhiSpike(etalead,philead,etasub,phisub,rad,eta,phi))
    return vfunc(Jets.eta,Jets.phi)

def METFilters(events):
    # MET filters
    BPFM = events.BadPFMuonFilter
    EBCR = events.ecalBadCalibFilter
    eeBS = events.eeBadScFilter
    gSTH = events.globalSuperTightHalo2016Filter
    HBHEIN = events.HBHEIsoNoiseFilter
    HBHEN = events.HBHENoiseFilter
    nV = events.NVtx

    return ((gSTH == 1) & (HBHEN == 1) & (HBHEIN == 1) & (BPFM == 1) &
    (EBCR == 1) & (eeBS == 1) & (nV > 0))

def Preselection(qualityCuts,nl):
    return (qualityCuts & (nl == 0))

def trgListtoInd(trigDict,trgList):
    return [trigDict.get("HLT_{}_v".format(trg)) for trg in trgList]

def PassTrigger(triggerPass,indices):
    triggerPass = ak.to_numpy(triggerPass)
    nTrigs= len(triggerPass[0])
    trigReq = []
    trigReq = np.zeros(nTrigs, dtype=int)
    trigReq[indices] = 1
    mult = triggerPass*trigReq
    return np.any(mult==1,axis=1)

def cutList(dataset,events,vars_noCut,SVJCut=True):
    evtw = vars_noCut["evtw"]
    nl = vars_noCut["nl"]
    nnim = vars_noCut["nnim"]
    njets = vars_noCut["njets"]
    njetsAK8 = vars_noCut["njetsAK8"]
    nsvjJetsAK8 = vars_noCut["nsvjJetsAK8"]
    nb = vars_noCut["nb"]
    met = vars_noCut["met"]
    ht = vars_noCut["ht"]
    dPhiMinj = vars_noCut["dPhiMinjMET"]
    dPhiMinjAK8 = vars_noCut["dPhiMinjMETAK8"]
    triggerPass = events.TriggerPass
    jetID = events.JetID
    jetIDAK8 = events.JetIDAK8
    ttStitch = TTStitch(dataset,events)
    metFilters = METFilters(events)
    # psFilter = PhiSpikeFilter(dataset,vars_noCut['jets'])
    qualityCuts = metFilters & (nl == 0)
    # qualityCuts = metFilters & psFilter # NN training files
    # preselection = Preselection(qualityCuts,nl)
    # cuts to get over trigger plateau
    metCut = met > 266
    htCut = ht > 1280
    stCut = (ht + met) > 1350
    trgPlat = metCut & htCut

    # trigger choices
    years = ["2016","2017","2018"]
    yr = 0
    for year in years:
        if year in dataset:
            yr = year
    trigDict = tD.trigDicts[yr]
    trgSelection = tD.trgSelections[yr]
    trgSelectionsQCDCR = tD.trgSelectionsQCDCR[yr]
    tch_trgs =  trgListtoInd(trigDict,trgSelection)
    tch_trgs_QCDCR =  trgListtoInd(trigDict,trgSelectionsQCDCR)
    passTrigger = PassTrigger(triggerPass,tch_trgs)

    # Define all cuts for histo making
    cuts = {
            ""                  : np.ones(len(evtw),dtype=bool),
            #"_qual"             : qualityCuts,
            #"_qual_met"         : qualityCuts & metCut,
            #"_qual_ht"          : qualityCuts & htCut,
            #"_qual_st"          : qualityCuts & stCut,
            #"_qual_trg"         : qualityCuts & passTrigger,
            #"_qual_trg_met"     : qualityCuts & passTrigger & metCut,
            #"_qual_trg_ht"      : qualityCuts & passTrigger & htCut,
            "_qual_trg_st"             : qualityCuts & passTrigger & stCut,        
            "_qual_trg_st_0nim"        : qualityCuts & passTrigger & stCut & (nnim == 0),
            "_qual_trg_st_0nim_2J"     : qualityCuts & passTrigger & stCut & (nnim == 0) & (nsvjJetsAK8 == 2),
            "_qual_trg_st_0nim_ge2J"   : qualityCuts & passTrigger & stCut & (nnim == 0) & (nsvjJetsAK8 >= 2),

            "_qual_trg_st_ge1nim"      : qualityCuts & passTrigger & stCut & (nnim >= 1),
            "_qual_trg_st_ge1nim_2J"   : qualityCuts & passTrigger & stCut & (nnim >= 1) & (nsvjJetsAK8 == 2),
            "_qual_trg_st_ge1nim_ge2J" : qualityCuts & passTrigger & stCut & (nnim >= 1) & (nsvjJetsAK8 >= 2),
            #"_metfilter_0l_1nim_trgQCDCR" : metFilters & (nl == 0) & (nnim == 1) & PassTrigger(triggerPass,tch_trgs_QCDCR),
    }

    # cuts with svj
    #if SVJCut == True:
    #     cuts["_pre_trg_1PSVJ"] =  preselection & triggerCut & (nsvjJetsAK8 >= 1)
    #     cuts["_pre_trg_2PSVJ"] =  preselection & triggerCut & (nsvjJetsAK8 >= 2)
    #     # nsvj characterization
    #     cuts["_pre_trg_0SVJ"] =   preselection & triggerCut & (nsvjJetsAK8 == 0)
    #     cuts["_pre_trg_1SVJ"] =   preselection & triggerCut & (nsvjJetsAK8 == 1)
    #     cuts["_pre_trg_2SVJ"] =   preselection & triggerCut & (nsvjJetsAK8 == 2)
    #     cuts["_pre_trg_3SVJ"] =   preselection & triggerCut & (nsvjJetsAK8 == 3)
    #     cuts["_pre_trg_4SVJ"] =   preselection & triggerCut & (nsvjJetsAK8 == 4)
    #     cuts["_pre_trg_5PSVJ"] =  preselection & triggerCut & (nsvjJetsAK8 >= 5)
    return cuts
