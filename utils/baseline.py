import numpy as np
import awkward as ak

def TTStitch(dataset,events):
    # # TT Stiching mask
    madHT = events.madHT
    nEle = ak.num(events.GenElectrons)
    nMu = ak.num(events.GenMuons)
    nTau = ak.num(events.GenTaus)
    genMET = events.GenMET

    ttStitchMask = None
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

def PassTrigger(triggerPass):
    indicesOfHighEffTrig = [11,12,13,14,67,107,108,131,8,90,98,116]
    triggerPass = ak.to_numpy(triggerPass)
    nTrigs= len(triggerPass[0])
    trigReq = []
    for i in range(nTrigs):
        if i in indicesOfHighEffTrig:
            trigReq.append(1)
        else:
            trigReq.append(0)
    return np.sum(np.multiply(trigReq,triggerPass),axis=1) > 0

def cutList(dataset,events,vars_noCut,SVJCut=True):
    evtw = vars_noCut["evtw"][0]
    nl = vars_noCut["nl"][0]
    njets = vars_noCut["njets"][0]
    njetsAK8 = vars_noCut["njetsAK8"][0]
    nb = vars_noCut["nb"][0]
    met = vars_noCut["met"][0]
    ht = vars_noCut["ht"][0]
    dPhiMinj = vars_noCut["dPhiMinjMET"][0]
    dPhiMinjAK8 = vars_noCut["dPhiMinjMETAK8"][0]
    deltaR12jAK8 = vars_noCut["dRJ12AK8"][0]
    triggerPass = events.TriggerPass
    jetID = events.JetID
    jetIDAK8 = events.JetIDAK8
    ttStitch = TTStitch(dataset,events)
    metFilters = METFilters(events)
    triggerCut = PassTrigger(triggerPass)
    # psFilter = PhiSpikeFilter(dataset,vars_noCut['jets'])
    # qualityCuts = metFilters & psFilter & ttStitch
    # qualityCuts = metFilters & psFilter # NN training files
    # preselection = Preselection(qualityCuts,nl)
    # cuts to get over trigger plateau
    metCut = met > 300
    htCut = 1250
    trgPlat = metCut & htCut

    cuts = {
            ""                          : np.ones(len(evtw),dtype=bool),
            "_trg"                         : triggerCut,
            "_trg_metfilter"             : triggerCut & metFilters,
            "_metfilter_0l"              : metFilters & (nl == 0),
            "_trg_metfilter_0l"          : triggerCut & metFilters & (nl == 0),
            # "_trg_metfilter_0l_trgPlat"  : triggerCut & metFilters & (nl == 0) & trgPlat,
    }
    # cuts with svj
    if SVJCut == True:
        nsvjJetsAK8 = vars_noCut["nsvjJetsAK8"][0]
        # cuts["_pre_trg_1PSVJ"] =  preselection & triggerCut & (nsvjJetsAK8 >= 1)
        # cuts["_pre_trg_2PSVJ"] =  preselection & triggerCut & (nsvjJetsAK8 >= 2)
        # # nsvj characterization
        # cuts["_pre_trg_0SVJ"] =   preselection & triggerCut & (nsvjJetsAK8 == 0)
        # cuts["_pre_trg_1SVJ"] =   preselection & triggerCut & (nsvjJetsAK8 == 1)
        # cuts["_pre_trg_2SVJ"] =   preselection & triggerCut & (nsvjJetsAK8 == 2)
        # cuts["_pre_trg_3SVJ"] =   preselection & triggerCut & (nsvjJetsAK8 == 3)
        # cuts["_pre_trg_4SVJ"] =   preselection & triggerCut & (nsvjJetsAK8 == 4)
        # cuts["_pre_trg_5PSVJ"] =  preselection & triggerCut & (nsvjJetsAK8 >= 5)
    return cuts
