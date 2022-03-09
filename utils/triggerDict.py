# UL2018 Ntuples Good triggers
trigDict = {
'AK8PFHT800_TrimMass50': 4,
'AK8PFHT850_TrimMass50': 5,
'AK8PFHT900_TrimMass50': 6,
'AK8PFJet400_TrimMass30': 8,
'AK8PFJet420_TrimMass30': 9,
'AK8PFJet500': 11,
'AK8PFJet550': 12,
'AK8PFJetFwd400': 13,
'CaloJet500_NoJetID': 14,
'CaloJet550_NoJetID': 15,
'CaloMET350_HBHECleaned': 16,
'DiPFJetAve300_HFJEC': 21,
'DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350': 22,
'DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350': 25,
'Ele115_CaloIdVT_GsfTrkIdT': 31,
'Ele135_CaloIdVT_GsfTrkIdT': 32,
'Ele145_CaloIdVT_GsfTrkIdT': 33,
'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ': 43,
'Ele23_Ele12_CaloIdL_TrackIdL_IsoVL': 44,
'Ele28_eta2p1_WPTight_Gsf_HT150': 48,
'Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned': 49,
'Ele32_WPTight_Gsf': 50,
'Ele35_WPTight_Gsf': 51,
'IsoMu27': 71,
'MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight': 76,
'MonoCentralPFJet80_PFMETNoMu130_PFMHTNoMu130_IDTight': 77,
'MonoCentralPFJet80_PFMETNoMu140_PFMHTNoMu140_IDTight': 78,
'Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ': 79,
'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8': 88,
'Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8': 89,
'Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ': 94,
'Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL': 95,
'Mu50': 101,
'Mu55': 102,
'Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ': 107,
'PFHT1050': 108,
'PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5': 115,
'PFHT500_PFMET100_PFMHT100_IDTight': 139,
'PFHT500_PFMET110_PFMHT110_IDTight': 140,
'PFHT700_PFMET85_PFMHT85_IDTight': 148,
'PFHT700_PFMET95_PFMHT95_IDTight': 149,
'PFHT800_PFMET75_PFMHT75_IDTight': 154,
'PFHT800_PFMET85_PFMHT85_IDTight': 155,
'PFJet500': 160,
'PFJet550': 161,
'PFJetFwd400': 162,
'PFMET120_PFMHT120_IDTight': 170,
'PFMET130_PFMHT130_IDTight': 172,
'PFMET140_PFMHT140_IDTight': 174,
'PFMET200_HBHE_BeamHaloCleaned': 176,
'PFMET250_HBHECleaned': 177,
'PFMET300_HBHECleaned': 178,
'PFMETNoMu120_PFMHTNoMu120_IDTight': 190,
'PFMETNoMu130_PFMHTNoMu130_IDTight': 192,
'PFMETNoMu140_PFMHTNoMu140_IDTight': 194,
'PFMETTypeOne140_PFMHT140_IDTight': 198,
'PFMETTypeOne200_HBHE_BeamHaloCleaned': 199,
'Photon200': 204,
'Photon300_NoHE': 205,
'TkMu100': 209
}

# s-channel 2018 triggers
schTriggers = [
    'AK8PFJet500',
    'AK8PFJet550',
    'CaloJet500_NoJetID',
    'CaloJet550_NoJetID',
    'PFHT1050',
    'PFJet500',
    'PFJet550'
]

# UL2018 Highest Efficiency Triggers (full t-channel signal samples)
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

# UL2018 Highest Efficiency Triggers (full t-channel signal samples; started  with the s-channel triggers)
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

# 2018 Highest Efficiency Triggers (full t-channel baseline signal sample only; started  with the s-channel triggers)
oldHEsch = [
    "PFMETNoMu120_PFMHTNoMu120_IDTight",
    "AK8PFJet400_TrimMass30",
    "PFHT500_PFMET100_PFMHT100_IDTight",
    "PFHT700_PFMET85_PFMHT85_IDTight",
    "PFMET120_PFMHT120_IDTight",
]
