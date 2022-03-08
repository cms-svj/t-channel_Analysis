import numpy as np
import awkward as ak

# UL2018 Ntuples Good triggers
trigDict = {
'HLT_AK8PFHT800_TrimMass50_v': 4,
'HLT_AK8PFHT850_TrimMass50_v': 5,
'HLT_AK8PFHT900_TrimMass50_v': 6,
'HLT_AK8PFJet400_TrimMass30_v': 8,
'HLT_AK8PFJet420_TrimMass30_v': 9,
'HLT_AK8PFJet500_v': 11,
'HLT_AK8PFJet550_v': 12,
'HLT_AK8PFJetFwd400_v': 13,
'HLT_CaloJet500_NoJetID_v': 14,
'HLT_CaloJet550_NoJetID_v': 15,
'HLT_CaloMET350_HBHECleaned_v': 16,
'HLT_DiPFJetAve300_HFJEC_v': 21,
'HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_v': 22,
'HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350_v': 25,
'HLT_Ele115_CaloIdVT_GsfTrkIdT_v': 31,
'HLT_Ele135_CaloIdVT_GsfTrkIdT_v': 32,
'HLT_Ele145_CaloIdVT_GsfTrkIdT_v': 33,
'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v': 43,
'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v': 44,
'HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v': 48,
'HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v': 49,
'HLT_Ele32_WPTight_Gsf_v': 50,
'HLT_Ele35_WPTight_Gsf_v': 51,
'HLT_IsoMu27_v': 71,
'HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_v': 76,
'HLT_MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight_v': 77,
'HLT_MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight_v': 78,
'HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v': 79,
'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v': 88,
'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v': 89,
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v': 94,
'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v': 95,
'HLT_Mu50_v': 101,
'HLT_Mu55_v': 102,
'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v': 107,
'HLT_PFHT1050_v': 108,
'HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5_v': 115,
'HLT_PFHT500_PFMET100_PFMHT100_IDTight_v': 139,
'HLT_PFHT500_PFMET110_PFMHT110_IDTight_v': 140,
'HLT_PFHT700_PFMET85_PFMHT85_IDTight_v': 148,
'HLT_PFHT700_PFMET95_PFMHT95_IDTight_v': 149,
'HLT_PFHT800_PFMET75_PFMHT75_IDTight_v': 154,
'HLT_PFHT800_PFMET85_PFMHT85_IDTight_v': 155,
'HLT_PFJet500_v': 160,
'HLT_PFJet550_v': 161,
'HLT_PFJetFwd400_v': 162,
'HLT_PFMET120_PFMHT120_IDTight_v': 170,
'HLT_PFMET130_PFMHT130_IDTight_v': 172,
'HLT_PFMET140_PFMHT140_IDTight_v': 174,
'HLT_PFMET200_HBHE_BeamHaloCleaned_v': 176,
'HLT_PFMET250_HBHECleaned_v': 177,
'HLT_PFMET300_HBHECleaned_v': 178,
'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v': 190,
'HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v': 192,
'HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v': 194,
'HLT_PFMETTypeOne140_PFMHT140_IDTight_v': 198,
'HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v': 199,
'HLT_Photon200_v': 204,
'HLT_Photon300_NoHE_v': 205,
'HLT_TkMu100_v': 209
}

schTriggers = {
'HLT_AK8PFJet500_v':11,
'HLT_AK8PFJet550_v':12,
'HLT_CaloJet500_NoJetID_v':14,
'HLT_CaloJet550_NoJetID_v':15,
'HLT_PFHT1050_v':108,
'HLT_PFJet500_v':160,
'HLT_PFJet550_v':161,
}

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

def trgListtoInd(trgList):
    indList = []
    for trg in trgList:
        indList.append(trigDict["HLT_{}_v".format(trg)])
    return indList

def PassTrigger(triggerPass,indices):
    triggerPass = ak.to_numpy(triggerPass)
    nTrigs= len(triggerPass[0])
    trigReq = []
    for i in range(nTrigs):
        if i in indices:
            trigReq.append(1)
        else:
            trigReq.append(0)
    mult = triggerPass*trigReq
    return np.any(mult==1,axis=1)

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
            "_metfilter_0l"              : metFilters & (nl == 0),
    }

    # trigger choices
    HETrg_noSch = [
        "AK8PFHT800_TrimMass50",
        "AK8PFJet400_TrimMass30",
        "AK8PFJet500",
        "AK8PFJetFwd400",
        "CaloJet500_NoJetID",
        "DiPFJetAve300_HFJEC",
        "PFHT1050",
        "PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5",
        "PFHT500_PFMET100_PFMHT100_IDTight",
        "PFHT700_PFMET85_PFMHT85_IDTight",
        "PFHT800_PFMET75_PFMHT75_IDTight",
        "PFMET120_PFMHT120_IDTight",
        "PFMETNoMu120_PFMHTNoMu120_IDTight",
        "PFMETTypeOne140_PFMHT140_IDTight",
        "PFMETTypeOne200_HBHE_BeamHaloCleaned"
    ]
    HETrg_noSch_ind =  trgListtoInd(HETrg_noSch)
    cuts["_HETrg_noSch"] = PassTrigger(triggerPass,HETrg_noSch_ind)
    HETrg_wSch = [
        "AK8PFHT800_TrimMass50",
        "AK8PFJet400_TrimMass30",
        "PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5",
        "PFHT500_PFMET100_PFMHT100_IDTight",
        "PFHT700_PFMET85_PFMHT85_IDTight",
        "PFHT800_PFMET75_PFMHT75_IDTight",
        "PFMET120_PFMHT120_IDTight",
        "PFMETNoMu120_PFMHTNoMu120_IDTight",
        "PFMETTypeOne200_HBHE_BeamHaloCleaned"
    ]
    HETrg_wSch_ind =  trgListtoInd(HETrg_wSch)
    schTrigIns = list(schTriggers.values())
    cuts["_HETrg_wSch"] = PassTrigger(triggerPass,HETrg_wSch_ind + schTrigIns)
    oldHEsch = [
        "PFMETNoMu120_PFMHTNoMu120_IDTight",
        "AK8PFJet400_TrimMass30",
        "PFHT500_PFMET100_PFMHT100_IDTight",
        "PFHT700_PFMET85_PFMHT85_IDTight",
        "PFMET120_PFMHT120_IDTight",
    ]
    oldHEsch_ind =  trgListtoInd(oldHEsch)
    cuts["_oldHEsch"] = PassTrigger(triggerPass,oldHEsch_ind + schTrigIns)
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
