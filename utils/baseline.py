import numpy as np
import awkward1 as ak

def ttStitch(df):
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

def phiSpikeFilter(df,Jets):
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

def preselection(qualityCuts,electrons,muons,met):
    return (qualityCuts & (electrons.counts == 0) & (muons.counts == 0) & (met > 100))
    # return ((electrons.counts == 0) & (muons.counts == 0) & (met > 100))


def passTrigger(triggerPass):
    indicesOfHighEffTrig = [11,12,13,14,67,107,108,131,8,90,98,116] # all s-channel + 5 highest signal efficiency

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

    # oneTrigger = ak.count(tPassedList,axis=-1) > 0
    return ak.count(tPassedHEList,axis=-1) > 0
