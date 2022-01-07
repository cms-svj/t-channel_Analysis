import numpy as np
import awkward1 as ak

def TTStitch(df):
    dataset = df['dataset']

    # # TT Stiching mask
    madHT = df['madHT']
    nEle = df['GenElectrons'].counts
    nMu = df['GenMuons'].counts
    nTau = df['GenTaus'].counts
    genMET = df['GenMET']

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
        ttStitchMask = np.ones(len(df['MET']),dtype=bool)

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

def PhiSpikeFilter(df,Jets):
    dataset = df['dataset']
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

def METFilters(df):
    # MET filters
    gSTH = df['globalSuperTightHalo2016Filter']
    HBHEN = df['HBHENoiseFilter']
    HBHEIN = df['HBHEIsoNoiseFilter']
    BPFM = df['BadPFMuonFilter']
    BCC = df['BadChargedCandidateFilter']
    EDCTP = df['ecalBadCalibReducedFilter']
    eeBS = df['eeBadScFilter']
    nV = df['NVtx']

    return ((gSTH == 1) & (HBHEN == 1) & (HBHEIN == 1) & (BPFM == 1) & (BCC == 1) &
    (EDCTP == 1) & (eeBS == 1) & (nV > 0))

def Preselection(qualityCuts,electrons,muons):
    return (qualityCuts & (electrons.counts == 0) & (muons.counts == 0))

def PassTrigger(triggerPass):
    # indicesOfHighEffTrig = [4, 5, 6, 8, 9, 11, 12, 13, 14, 17, 22, 23, 24, 34, 38, 39, 40, 50, 65, 66, 67, 90, 91, 98, 99, 101, 102, 107, 108, 116, 118, 120, 131, 133, 135, 141, 142, 146] # all good triggers
    indicesOfHighEffTrig = [11,12,13,14,67,107,108,131,8,90,98,116] # all s-channel + 5 highest signal efficiency (no muon trigger)
    # indicesOfHighEffTrig = [11,12,13,14,67,107,108,131,116] # all s-channel + 2 of the 5 highest signal efficiency (no cross MET, HT triggers, or PFJet trigger); for sanity check
    # indicesOfHighEffTrig = [11,12,13,14,67,107,108] # all s-channel
    tPassedHEList = []
    tPassedList = []
    for evt in triggerPass:
        tPassed = []
        tPassedHE = []
        for tp in range(len(evt)):
            if evt[tp] == 1:
                if tp in indicesOfHighEffTrig:
                    tPassedHE.append(tp)
        tPassedList.append(tPassed)
        tPassedHEList.append(tPassedHE)
    tPassedList = ak.Array(tPassedList)
    tPassedHEList = ak.Array(tPassedHEList)
    return ak.count(tPassedHEList,axis=-1) > 0

def cutList(df,inpObj):
    evtw = inpObj["evtw"]
    electrons = inpObj["electrons"]
    muons = inpObj["muons"]
    jets = inpObj["jets"]
    fjets = inpObj["fjets"]
    bjets = inpObj["bjets"]
    met = inpObj["met"]
    triggerPass = inpObj["triggerPass"]
    jetID = inpObj["jetID"]
    jetIDAK8 = inpObj["jetIDAK8"]
    ht = inpObj["ht"]
    dPhiMinj = inpObj["dPhiMinj"]
    dPhiMinjAK8 = inpObj["dPhiMinjAK8"]
    deltaR12jAK8 = inpObj["deltaR12jAK8"]
    ttStitch = TTStitch(df)
    metFilters = METFilters(df)
    triggerCut = PassTrigger(triggerPass)
    psFilter = PhiSpikeFilter(df,jets)
    qualityCuts = metFilters & psFilter & ttStitch
    preselection = Preselection(qualityCuts,electrons,muons)

    cuts = {
            ""                          : np.ones(len(evtw),dtype=bool),
            # "_trg"                      : triggerCut,
            # "_qc"                       : ttStitch & metFilters & psFilter,
            # "_qc_trg"                   : ttStitch & metFilters & psFilter & triggerCut,
            # "_qc_trg_0l"                : ttStitch & metFilters & psFilter & triggerCut & (electrons.counts == 0) & (muons.counts == 0),
            "_pre"                      : preselection,
            # "_pre_ge2AK4j"              : preselection & (jets.counts >= 2) & (jetID == True),
            # "_pre_ge2AK8j"              : preselection & (fjets.counts >= 2) & (jetIDAK8 == True),
            "_pre_trg"                  : preselection & triggerCut,
            ##  cuts for optimizing significance for full t-channel production
            # "_pre_ge4AK8j_trg"          : preselection & (fjets.counts >= 4) & (jetIDAK8 == True) & triggerCut,
            # "_pre_ge4AK8j_trg_met100"   : preselection & (fjets.counts >= 4) & (jetIDAK8 == True) & triggerCut & (met > 100),
            # "_pre_ge4AK8j_trg_nb2"      : preselection & (fjets.counts >= 4) & (jetIDAK8 == True) & triggerCut & (bjets.counts >= 2),
            # "_pre_ge4AK8j_trg_dR3p3"    : preselection & (fjets.counts >= 4) & (jetIDAK8 == True) & triggerCut & (deltaR12jAK8 < 3.3),
            # "_pre_ge4AK8j_trg_met100_dR3p3"    : preselection & (fjets.counts >= 4) & (jetIDAK8 == True) & triggerCut & (met > 100) & (deltaR12jAK8 < 3.3),
            # "_pre_ge4AK8j_trg_nb2_dR3p3"       : preselection & (fjets.counts >= 4) & (jetIDAK8 == True) & triggerCut & (bjets.counts >= 2) & (deltaR12jAK8 < 3.3),
            # "_pre_ge4AK8j_trg_met100_nb2"      : preselection & (fjets.counts >= 4) & (jetIDAK8 == True) & triggerCut & (met > 100) & (bjets.counts >= 2),
            # "_pre_ge4AK8j_trg_met100_nb2_dR3p3"      : preselection & (fjets.counts >= 4) & (jetIDAK8 == True) & triggerCut & (met > 100) & (bjets.counts >= 2) & (deltaR12jAK8 < 3.3),
            ## optimum cuts for pairProduction
            "_pre_trg_nb5"              : preselection & triggerCut & (bjets.counts >= 5),
            "_pre_trg_met300"           : preselection & triggerCut & (met > 300),
            "_pre_trg_nb5_met300"       : preselection & triggerCut & (bjets.counts >= 5) & (met > 300),
            ## cuts for NN training files
            "_npz"                      : metFilters & psFilter & triggerCut & (electrons.counts == 0) & (muons.counts == 0),
    }
    return cuts
