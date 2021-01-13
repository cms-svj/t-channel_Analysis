import numpy as np

def ttStich(df):
    dataset = df['dataset']

    # # TT Stiching mask
    madHT = df['madHT']
    nEle = df['GenElectrons'].counts
    nMu = df['GenMuons'].counts
    nTau = df['GenTaus'].counts
    genMET = df['GenMET']

    ttStitchMask = None
    if "TTJets_Inc" in dataset:
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

def vetoPhiSpike(etaHSL,phiHSL,rad,eta,phi):
    veto = np.ones(len(eta),dtype=bool)
    for iep in range(len(etaHSL)):
        veto = veto & ( (etaHSL[iep] - eta)**2 + (phiHSL[iep] - phi)**2 > rad )

    return veto

def phiSpikeFilter(df,Jets):
    dataset = df['dataset']
    # phi spike filters
    rad = 0.028816 # half the length of the diagonal of the eta-phi rectangular cell
    rad = rad * 0.35 # the factor of 0.35 was optimized from the signal vs. background sensitivity study

    # hot spots for leading jets
    eta16lead = [0.048,0.24,1.488,1.584,-1.008]
    phi16lead = [-0.35,-0.35,-0.77,-0.77,-1.61]

    eta17lead = [0.144,1.488,1.488,1.584,-0.624]
    phi17lead = [-0.35,-0.77,-0.63,-0.77,0.91]

    eta18lead = [1.488,1.488,1.584]
    phi18lead = [-0.77,-0.63,-0.77]

    # hot spots for subleading jets
    eta16sub = [-1.2,-0.912,-0.912,-0.816,-0.72,-0.72,-0.528,-0.432,-0.336,-0.24,-0.24,-0.144,-0.144,-0.048,0.144,
    0.912,0.912,1.008,1.296,-1.584,-0.816,-0.72,-0.144,-0.048,-0.048,0.048,1.104,1.488]
    phi16sub = [-1.19,2.03,3.01,-1.75,-2.17,-0.77,2.73,2.73,0.21,0.07,0.21,-2.59,0.77,0.91,1.75,1.75,2.87,0.63,
    -0.49,0.63,1.47,-2.31,0.07,-2.59,0.77,0.91,-3.15,2.73]

    eta17sub = [-0.912,-0.912,-0.816,-0.72,-0.528,-0.336,-0.24,-0.24,-0.144,-0.144,-0.048,0.144,0.912,0.912,1.008,
    -1.2,-0.72,-0.72,-0.432,0.336,0.624,1.104,1.296]
    phi17sub = [2.03,3.01,-1.75,-0.77,2.73,0.21,0.07,0.21,-2.59,0.77,0.91,1.75,1.75,2.87,0.63,-1.19,-2.31,-2.17,
    2.73,-0.77,-0.77,-3.15,-0.49]

    eta18sub = [-1.584,-1.2,-0.912,-0.912,-0.816,-0.816,-0.72,-0.72,-0.528,-0.432,-0.336,-0.24,-0.24,-0.144,-0.144,
    -0.144,-0.048,-0.048,0.144,0.912,0.912,1.008,1.296,-0.72,1.104,1.488,1.776]
    phi18sub = [0.63,-1.19,2.03,3.01,-1.75,-0.77,-2.17,-0.77,2.73,2.73,0.21,0.07,0.21,-2.59,0.07,0.77,0.77,0.91,
    1.75,1.75,2.87,0.63,-0.49,-2.31,-3.15,-0.21,0.77]

    psFilter = None
    if "16" in dataset:
        psFilter = vetoPhiSpike(eta16sub,phi16sub,rad,Jets.eta[:,1],Jets.phi[:,1]) & vetoPhiSpike(eta16lead,phi16lead,rad,Jets.eta[:,0],Jets.phi[:,0])
    elif "17" in dataset:
        psFilter = vetoPhiSpike(eta17sub,phi17sub,rad,Jets.eta[:,1],Jets.phi[:,1]) & vetoPhiSpike(eta17lead,phi17lead,rad,Jets.eta[:,0],Jets.phi[:,0])
    elif "18" in dataset:
        psFilter = vetoPhiSpike(eta18sub,phi18sub,rad,Jets.eta[:,1],Jets.phi[:,1]) & vetoPhiSpike(eta18lead,phi18lead,rad,Jets.eta[:,0],Jets.phi[:,0])
    else:
        psFilter = np.ones(len(Jets.eta[:,1]),dtype=bool)

    return psFilter

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

def preselection(electrons,muons,met):
    return ((electrons.counts == 0) & (muons.counts == 0))
    # return ((electrons.counts == 0) & (muons.counts == 0) & (met > 100))
