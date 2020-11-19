import numpy as np

def getBaselineMask(df):
    dataset = df['dataset']
    # # TT Stiching mask
    ttStitchMask = np.ones(len(df['MET']),dtype=bool)

    madHT = df['madHT']
    nEle = df['GenElectrons'].counts
    nMu = df['GenMuons'].counts
    nTau = df['GenTaus'].counts
    genMET = df['GenMET']

    # put the stitching in utility
    if "TTJets_Inc" in dataset:
        ttStitchMask = (madHT < 600) & (nEle==0) & (nMu==0) & (nTau==0)
    elif "TTJets_HT" in dataset:
        ttStitchMask = (madHT >= 600)
    elif ("TTJets_DiLept" in dataset and "genMET" not in dataset) or ("TTJets_SingleLeptFromT" in dataset and "genMET" not in dataset):
        ttStitchMask = (madHT < 600) & (genMET < 150)
    elif ("TTJets_DiLept" in dataset and "genMET" in dataset) or ("TTJets_SingleLeptFromT" in dataset and "genMET" in dataset):
        ttStitchMask = (madHT < 600) & (genMET >= 150)

    return ttStitchMask
