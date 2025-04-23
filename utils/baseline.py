import numpy as np
import awkward as ak
import pickle
from . import triggerDict as tD
from . import utility as util
from . import samples as s
# dataset here is just a string that is the input of -d
# events are the tree in the ntuple

def hemVeto(ak4Jets,electrons,muons,etaMin,etaMax,phiMin,phiMax):
    ak4jHemCond = (ak4Jets.eta > etaMin) & (ak4Jets.eta < etaMax) & (ak4Jets.phi > phiMin) & (ak4Jets.phi < phiMax)
    elecHemCond = (electrons.eta > etaMin) & (electrons.eta < etaMax) & (electrons.phi > phiMin) & (electrons.phi < phiMax)
    muonHemCond = (muons.eta > etaMin) & (muons.eta < etaMax) & (muons.phi > phiMin) & (muons.phi < phiMax)
    # hemEvents = ((ak.num(ak4Jets) > 0) & ak.any(ak4jHemCond,axis=1)) | ((ak.num(muons) > 0) & ak.any(muonHemCond,axis=1)) | ((ak.num(electrons) > 0) & ak.any(elecHemCond,axis=1))
    hemEvents = ak.any(ak4jHemCond,axis=1) | ak.any(muonHemCond,axis=1) | ak.any(elecHemCond,axis=1)
    return ~hemEvents

def hemPeriodMask(dataset,events,ak4Jets,electrons,muons,hemPeriod,etaMin,etaMax,phiMin,phiMax):
    runNum = events.RunNum
    if ("2018" in dataset):
        if (hemPeriod == "PreHEM") and ("Data" in dataset):
            return runNum < 319077
        elif (hemPeriod == "PostHEM") and ("Data" in dataset):
            return (runNum >= 319077) & hemVeto(ak4Jets,electrons,muons,etaMin,etaMax,phiMin,phiMax)
        elif (hemPeriod == "") and ("Data" in dataset):
            return (runNum < 319077) | ((runNum >= 319077) & hemVeto(ak4Jets,electrons,muons,etaMin,etaMax,phiMin,phiMax))
        elif hemPeriod == "PostHEM":
            return hemVeto(ak4Jets,electrons,muons,etaMin,etaMax,phiMin,phiMax)
        else:
            return np.ones(len(events),dtype=bool)
    else:
        return np.ones(len(events),dtype=bool)

def ttStitchMask(dataset,events):
    # # TT Stiching mask
    ttStitchMask = None
    if "TTJets" in dataset:
        madHT = events.madHT
        nEle = ak.num(events.GenElectrons)
        nMu = ak.num(events.GenMuons)
        nTau = ak.num(events.GenTaus)
        genMET = events.GenMET
        if "TTJets_Incl" in dataset or "mTTJetsmini_Incl" in dataset:
            ttStitchMask = (madHT < 600) & (nEle==0) & (nMu==0) & (nTau==0)
        elif "TTJets_HT" in dataset:
            ttStitchMask = (madHT >= 600)
        elif ("TTJets_DiLept" in dataset and "genMET" not in dataset) or ("TTJets_SingleLeptFromT" in dataset and "genMET" not in dataset):
            ttStitchMask = (madHT < 600) & (genMET < 150)
        elif ("TTJets_DiLept" in dataset and "genMET" in dataset) or ("TTJets_SingleLeptFromT" in dataset and "genMET" in dataset):
            ttStitchMask = (madHT < 600) & (genMET >= 150)
    else:
        ttStitchMask = np.ones(len(events),dtype=bool) # passes all true

    return ttStitchMask

def phiSpikeVeto(hotSpotsDict,varName,jEtai,jPhii,rad):
    hotEtasi, hotPhisi = hotSpotsDict[varName]
    jEtaiReshaped = np.broadcast_to(list(jEtai),(len(hotEtasi),len(jEtai)))
    hotEtasiReshaped = np.reshape(hotEtasi,(len(hotEtasi),1))
    jPhiiReshaped = np.broadcast_to(list(jPhii),(len(hotPhisi),len(jPhii)))
    hotPhisiReshaped = np.reshape(hotPhisi,(len(hotPhisi),1))
    return np.prod((jEtaiReshaped - hotEtasiReshaped)**2 + (jPhiiReshaped - hotPhisiReshaped)**2 > rad, axis=0, dtype=bool)

def phiSpikeFilter(yr, jets, hotSpotsPkl):
    with open(hotSpotsPkl,"rb") as infile:
        phiSpikeHotSpots = pickle.load(infile)
    rad = 0.028816*0.35 # the factor of 0.35 was optimized from the signal vs. background sensitivity study for s-channel
    hotSpotsDict = phiSpikeHotSpots[yr]
    eta = jets.eta
    phi = jets.phi
    pass1 = phiSpikeVeto(hotSpotsDict,"j1Phivsj1Eta",util.jetVar_i(eta,0),util.jetVar_i(phi,0),rad)
    pass2 = phiSpikeVeto(hotSpotsDict,"j2Phivsj2Eta",util.jetVar_i(eta,1),util.jetVar_i(phi,1),rad)
    pass3 = phiSpikeVeto(hotSpotsDict,"j3Phivsj3Eta",util.jetVar_i(eta,2),util.jetVar_i(phi,2),rad)
    pass4 = phiSpikeVeto(hotSpotsDict,"j4Phivsj4Eta",util.jetVar_i(eta,3),util.jetVar_i(phi,3),rad)
    return pass1 & pass2 & pass3 & pass4


def metFilterMask(events):
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

def removeOverlap(dataset, events, yr):
    ### coffea code ###
    """Function to remove overlap between different dataset by removing the events with overlapping triggers. Currently the preferred dataset is JetHT.""" 
    trgDict = tD.trigDicts[yr]
    trgListMet = trgListtoInd(trgDict,tD.trgMET[yr])
    trgListJetHT = trgListtoInd(trgDict,tD.trgJetHT[yr])
    passMetTrg = passTriggerMask(events.TriggerPass, trgListMet)
    passJetHTTrg = passTriggerMask(events.TriggerPass, trgListJetHT)
    
    OverlapDataMask = None
    if "METData" in dataset:
        OverlapDataMask = np.bitwise_not(passJetHTTrg)
    elif "HTMHTData" in dataset:
        OverlapDataMask = np.bitwise_not(passJetHTTrg + passMetTrg)
    else:
        OverlapDataMask = np.ones(len(events),dtype=bool)
    #Preliminary checks
    # for i in range(0,len(events)):
    #     if passMetTrg[i] and passJetHTTrg[i]:
    #         print("Event no = {} | Lumi no = {} | Run no. = {} | PassMetTrg = {} | PassJetHTTrg = {} | PassHTMHT = {} | Mask = {} --- ".format(i,events.EvtNum[i], events.LumiBlockNum[i], events.RunNum[i],passMetTrg[i],passJetHTTrg[i],passHTMHTTrg[i],OverlapDataMask[i]))
    #########

    return OverlapDataMask

def removeOverlapCR(dataset, events, yr):
    """Function to remove overlap between different dataset by removing the events with overlapping triggers."""     
    # For single lepton control region - Main dataset is Single Muon, removing overlap from the Single Electron
    passSingleM = passTriggerMask(events.TriggerPass,trgListtoInd(tD.trigDicts[yr],tD.trgSingleM[yr]))
    passSingleE = passTriggerMask(events.TriggerPass,trgListtoInd(tD.trigDicts[yr],tD.trgSingleE[yr])) 

    OverlapDataMask = None
    # Adding single lepton control region overlap removal
    if yr == 2018:
        if "EGammaData" in dataset:
            OverlapDataMask = np.bitwise_not(passSingleM)
        else:
            OverlapDataMask = np.ones(len(events),dtype=bool)
    else:
        if "SingleElectronData" in dataset:
            OverlapDataMask = np.bitwise_not(passSingleM)
            
        elif "SinglePhotonData" in dataset:
            OverlapDataMask = np.bitwise_not(passSingleM + passSingleE)
        else:
            OverlapDataMask = np.ones(len(events),dtype=bool)
  
    return OverlapDataMask

def trgListtoInd(trigDict,trgList):
    return [trigDict.get("HLT_{}_v".format(trg)) for trg in trgList]

def passTriggerMask(triggerPass,indices):
    triggerPass = ak.to_numpy(triggerPass)
    nTrigs= len(triggerPass[0])
    trigReq = []
    trigReq = np.zeros(nTrigs, dtype=int)
    trigReq[indices] = 1
    mult = triggerPass*trigReq
    return np.any(mult==1,axis=1) # correct version

def conditionMask(var, condition):
    FullMask = []
    for val in var:
        Mask = None
        Mask = (val == 0)
        FullMask.append(Mask)
        print(val)
    print("Mask = ",FullMask)

def cutList(dataset,events,vars_noCut,hemStudy,trgEffStudy,hemPeriod,skimCut,skimSource,runJetTag=True):
    evtw = vars_noCut["evtw"]
    nl = vars_noCut["nl"]
    # nnim = vars_noCut["nnim"] # no of Isolated Muons
    # njets = vars_noCut["njets"]
    njetsAK8 = vars_noCut["njetsAK8"]
    nb = vars_noCut["nb"]
    met = vars_noCut["MET"]
    ht = vars_noCut["HT"]
    st = vars_noCut["ST"]
    JAK8Pt = vars_noCut["J1AK8Pt"]
    # dPhiMinj = vars_noCut["dPhiMinjMET"]
    dPhiMinjAK8 = vars_noCut["dPhiMinjMETAK8"]
    nelectron = ak.num(vars_noCut["electrons"])
    nmuon = ak.num(vars_noCut["muons"])
    ncrMuons = ak.num(vars_noCut["crMuons"])
    ncrElectrons = ak.num(vars_noCut["crElectrons"])
    mtMetCRMuon = vars_noCut["mtMETCRMuon"]
    mtMetCRElectron = vars_noCut["mtMETCRElectron"]

    triggerPass = events.TriggerPass
    jetID = events.JetID
    jetIDAK8 = vars_noCut["fJetsID"]
    ttStitch = ttStitchMask(dataset,events)

    years = ["2016","2017","2018"]
    yr = 0
    for year in years:
        if year in dataset:
            yr = year

    psFilterFolder = "utils"
    psFilter =   phiSpikeFilter(yr, vars_noCut["jets"], f"{psFilterFolder}/phiSpikeHotSpots.pkl")

    if f"{yr}_Data" == s.getGroupFromSample(dataset,skimCut,skimSource=skimSource):
        DataMask = removeOverlap(dataset, events, yr)
    else: 
        DataMask = np.ones(len(events),dtype=bool)
    metFilters = metFilterMask(events)
    # values taken from HEM veto optimization study
    etaMin = -3.05
    etaMax = -1.35
    phiMin = -1.62
    phiMax = -0.82
    hemMask = hemPeriodMask(dataset,events,vars_noCut["jets"],vars_noCut["electrons"],vars_noCut["muons"],hemPeriod,etaMin,etaMax,phiMin,phiMax)
    # psFilter = phiSpikeFilter(dataset,vars_noCut['jets'])
    qualityCuts = metFilters & (nl == 0) & ttStitch & DataMask & hemMask & psFilter
    qualityWithLepton = metFilters & ttStitch & DataMask & hemMask & psFilter
    # qualityCuts = metFilters & psFilter # NN training files
    # cuts to get over trigger plateau
    # metCut = met > 266
    # htCut = ht > 1280
    htCut = ht > 500
    stCut = st > 1300
    metcut = met > 200

    # trgPlat = metCut & htCut

    # trigger choices
    trigDict = tD.trigDicts[yr]
    trgSelection = tD.trgSelections[yr]
    trgSelectionsCR = tD.trgSelectionsCR[yr]
    trgSelectionsQCDCR = tD.trgSelectionsQCDCR[yr]
    tch_trgs =  trgListtoInd(trigDict,trgSelection)
    tch_trgs_CR =  trgListtoInd(trigDict,trgSelectionsCR)
    tch_trgs_QCDCR =  trgListtoInd(trigDict,trgSelectionsQCDCR)
    passTrigger = passTriggerMask(triggerPass,tch_trgs)
    # trigger study for MCs
    # nOffMuons = vars_noCut['nOffMuons']
    passTrigger_muon = passTriggerMask(triggerPass,tch_trgs_CR)
    # preselection_offLineMuons = metFilters & ttStitch & DataMask & hemMask & passTrigger_muon & (njetsAK8 >= 2) & (nOffMuons >= 1)
    # preselection_offLineMuons_tchTrg = preselection_offLineMuons & passTrigger
    cr_muon_cut                  = qualityWithLepton & passTrigger & (ncrMuons == 1)     & (nelectron == 0) & stCut & (njetsAK8 >= 2) & (dPhiMinjAK8 <= 1.5) & metcut
    cr_electron_cut              = qualityWithLepton & passTrigger & (ncrElectrons == 1) & (nmuon == 0) & stCut & (njetsAK8 >= 2) & (dPhiMinjAK8 <= 1.5) & metcut
    preselection =                 qualityCuts       & passTrigger & (njetsAK8 >=2) & stCut & (dPhiMinjAK8 <= 1.5) & metcut & jetIDAK8    
    preselection_noMETCut =        qualityCuts       & passTrigger & (njetsAK8 >=2) & stCut & (dPhiMinjAK8 <= 1.5) & jetIDAK8    
    lcr_preselection_loose =       qualityWithLepton & passTrigger & (njetsAK8 >=2) & stCut & (dPhiMinjAK8 <= 1.5) & metcut & jetIDAK8 & (nl >= 1) 
    lcr_preselection =             qualityWithLepton & passTrigger & (njetsAK8 >=2) & stCut & (dPhiMinjAK8 <= 1.5) & metcut & jetIDAK8 & (nl == 1)
    lcr_preselection_noMETCut =    qualityWithLepton & passTrigger & (njetsAK8 >=2) & stCut & (dPhiMinjAK8 <= 1.5) & jetIDAK8 & (nl == 1)
    cuts = {
            "_nocut":                       np.ones(len(evtw),dtype=bool), 
            # "_pdOverlap":                   DataMask,
            "_trigger":                     DataMask & passTrigger,
            "_st":                          DataMask & passTrigger & stCut,
            "_metFilters":                  DataMask & passTrigger & stCut & metFilters,
            "_hemVeto":                     DataMask & passTrigger & stCut & metFilters & hemMask,
            "_goodJetFilters":              DataMask & passTrigger & stCut & metFilters & hemMask & jetIDAK8,
            "_2JetsAK8":                    DataMask & passTrigger & stCut & metFilters & hemMask & jetIDAK8 & (njetsAK8 >=2),
            "_leptonVeto":                  DataMask & passTrigger & stCut & metFilters & hemMask & jetIDAK8 & (njetsAK8 >=2) & (nl == 0),
            "_dPhiMin":                     DataMask & passTrigger & stCut & metFilters & hemMask & jetIDAK8 & (njetsAK8 >=2) & (nl == 0) & (dPhiMinjAK8 <= 1.5),
            "_met":                         DataMask & passTrigger & stCut & metFilters & hemMask & jetIDAK8 & (njetsAK8 >=2) & (nl == 0) & (dPhiMinjAK8 <= 1.5) & metcut,
            "_psFilter":                    DataMask & passTrigger & stCut & metFilters & hemMask & jetIDAK8 & (njetsAK8 >=2) & (nl == 0) & (dPhiMinjAK8 <= 1.5) & metcut & psFilter,
            "_metOnly":                     metcut,
            "_pre_noHEM":                   DataMask & passTrigger & stCut & metFilters & jetIDAK8 & (njetsAK8 >=2) & (nl == 0) & (dPhiMinjAK8 <= 1.5) & metcut & ttStitch,
            "_pre":                         preselection,
            # "_pre_psFilterSig2":            preselection & psFilterSig2,
            # "_pre_psFilterSig2p5":          preselection & psFilterSig2p5,
            # "_pre_psFilterSig3":            preselection & psFilterSig3,
            # "_pre_psFilterSig3p5":          preselection & psFilterSig3p5,
            # lost lepton control region
            "_lcr_pre":                       lcr_preselection,
            # "_lcr_pre_loose":               lcr_preselection_loose,
            # "_lcr_pre_noMet":               lcr_preselection_noMETCut,
            # "_cr_muon_":                    cr_muon_cut, 
            # "_cr_electron_":                cr_electron_cut,
    }

    if hemStudy:
        cutsHemStudy = {}
        etaLows = np.arange(-3.0,-3.51,-0.05)
        etaHighs = np.arange(-1.40,-0.89,0.05)
        phiLows = np.arange(-1.57,-2.08,-0.05)
        phiHighs = np.arange(-0.87,-0.36,0.05)
        for i in range(len(etaLows)):
            hemMask = hemPeriodMask(dataset,events,vars_noCut["jets"],vars_noCut["electrons"],vars_noCut["muons"],hemPeriod,etaLows[i],etaHighs[i],phiLows[i],phiHighs[i])
            etaCutLabel = str(np.round(etaLows[i],2)).replace(".","p").replace("-","")
            cutsHemStudy[f"_pre_eta_{etaCutLabel}"] = DataMask & passTrigger & stCut & metFilters & jetIDAK8 & (njetsAK8 >=2) & (nl == 0) & (dPhiMinjAK8 <= 1.5) & metcut & hemMask
            cutsHemStudy[f"_lcr_eta_{etaCutLabel}"] = DataMask & passTrigger & stCut & metFilters & jetIDAK8 & (njetsAK8 >=2) & (nl == 1) & (dPhiMinjAK8 <= 1.5) & metcut & hemMask
            cutsHemStudy[f"_lcr_loose_eta_{etaCutLabel}"] = DataMask & passTrigger & stCut & metFilters & jetIDAK8 & (njetsAK8 >=2) & (nl >= 1) & (dPhiMinjAK8 <= 1.5) & metcut & hemMask
        cuts.update(cutsHemStudy)

    if trgEffStudy:
        cuts = {
                    "_nocut":                       np.ones(len(evtw),dtype=bool),
                    "_trigger":                     DataMask & passTrigger,
                } 
        for trgName in trgSelection:
            trgList = []
            for trg in trgSelection:
                if trg != trgName:
                    trgList.append(trg)
            tch_trgs =  trgListtoInd(trigDict,trgList)
            passTriggerInd = passTriggerMask(triggerPass,tch_trgs)
            cuts["_all_except_"+trgName] = passTriggerInd

    # orthogonal dataset: SingleMuon
    if "Muon" in dataset:
        nHLTMatchedMuons = vars_noCut['nHLTMatchedMuons'] 
        preselection_matchedHLTMuons = metFilters & ttStitch & DataMask & hemMask & passTrigger_muon & (njetsAK8 >= 2) & (nHLTMatchedMuons >= 1)
        # preselection_offLineMuons = metFilters & ttStitch & DataMask & hemMask & passTrigger_muon & (njetsAK8 >= 2) & (nOffMuons >= 1)
        preselection_matchedHLTMuons_tchTrg = preselection_matchedHLTMuons & passTrigger
        # preselection_offLineMuons_tchTrg = preselection_offLineMuons & passTrigger
        # cuts = {
        #             ""                                          : np.ones(len(evtw),dtype=bool),
        #             "_preselec_matchedHLTMuons"                 : preselection_matchedHLTMuons,
        #             "_preselec_offLineMuons"                    : preselection_offLineMuons,  
        #             "_preselection_matchedHLTMuons_tchTrg"      : preselection_matchedHLTMuons_tchTrg,
        #             "_preselection_offLineMuons_tchTrg"         : preselection_offLineMuons_tchTrg,
        #             "_preselec_matchedHLTMuons_st"              : preselection_matchedHLTMuons & stCut,
        #             "_preselec_offLineMuons_st"                 : preselection_offLineMuons & stCut,  
        #             "_preselection_matchedHLTMuons_tchTrg_st"   : preselection_matchedHLTMuons_tchTrg & stCut,
        #             "_preselection_offLineMuons_tchTrg_st"      : preselection_offLineMuons_tchTrg & stCut,
        # }
    if runJetTag:
        nsvjJetsAK8 = vars_noCut["nsvjJetsAK8"]
        nsvjWNAE = vars_noCut["nsvjWNAE"]
        cutsWithNSVJ = {
            "_pre_0SVJ":                        preselection & (nsvjJetsAK8 == 0),
            "_pre_1SVJ":                        preselection & (nsvjJetsAK8 == 1),
            "_pre_2SVJ":                        preselection & (nsvjJetsAK8 == 2),
            "_pre_3SVJ":                        preselection & (nsvjJetsAK8 == 3),
            "_pre_4PSVJ":                       preselection & (nsvjJetsAK8 >= 4),


            "_pre_WNAE_0SVJ":                        preselection & (nsvjWNAE == 0),
            "_pre_WNAE_1SVJ":                        preselection & (nsvjWNAE == 1),
            "_pre_WNAE_2SVJ":                        preselection & (nsvjWNAE == 2),
            "_pre_WNAE_3SVJ":                        preselection & (nsvjWNAE == 3),
            "_pre_WNAE_4PSVJ":                       preselection & (nsvjWNAE >= 4),

            "_lcr_pre_0SVJ":                    lcr_preselection & (nsvjJetsAK8 == 0),
            "_lcr_pre_1SVJ":                    lcr_preselection & (nsvjJetsAK8 == 1),
            "_lcr_pre_2SVJ":                    lcr_preselection & (nsvjJetsAK8 == 2),
            "_lcr_pre_2PSVJ":                   lcr_preselection & (nsvjJetsAK8 >= 2),
            "_lcr_pre_3SVJ":                    lcr_preselection & (nsvjJetsAK8 == 3),
            "_lcr_pre_4PSVJ":                   lcr_preselection & (nsvjJetsAK8 >= 4),
            "_lcr_pre_noMet_0SVJ":              lcr_preselection_noMETCut & (nsvjJetsAK8 == 0),
            "_lcr_pre_noMet_1SVJ":              lcr_preselection_noMETCut & (nsvjJetsAK8 == 1),
            "_lcr_pre_noMet_2SVJ":              lcr_preselection_noMETCut & (nsvjJetsAK8 == 2),
            "_lcr_pre_noMet_2PSVJ":             lcr_preselection_noMETCut & (nsvjJetsAK8 >= 2),
            "_lcr_pre_noMet_3SVJ":              lcr_preselection_noMETCut & (nsvjJetsAK8 == 3),
            "_lcr_pre_noMet_4PSVJ":             lcr_preselection_noMETCut & (nsvjJetsAK8 >= 4),
            "_cr_muon_0SVJ":                    cr_muon_cut & (nsvjJetsAK8 == 0), 
            "_cr_muon_1SVJ":                    cr_muon_cut & (nsvjJetsAK8 == 1), 
            "_cr_muon_2SVJ":                    cr_muon_cut & (nsvjJetsAK8 == 2), 
            "_cr_muon_2PSVJ":                   cr_muon_cut & (nsvjJetsAK8 >= 2), 
            "_cr_muon_3SVJ":                    cr_muon_cut & (nsvjJetsAK8 == 3), 
            "_cr_muon_4PSVJ":                   cr_muon_cut & (nsvjJetsAK8 >= 4), 
            "_cr_electron_0SVJ":                cr_electron_cut & (nsvjJetsAK8 == 0),
            "_cr_electron_1SVJ":                cr_electron_cut & (nsvjJetsAK8 == 1),
            "_cr_electron_2SVJ":                cr_electron_cut & (nsvjJetsAK8 == 2),
            "_cr_electron_2PSVJ":               cr_electron_cut & (nsvjJetsAK8 >= 2),
            "_cr_electron_3SVJ":                cr_electron_cut & (nsvjJetsAK8 == 3),
            "_cr_electron_4PSVJ":               cr_electron_cut & (nsvjJetsAK8 >= 4),
            "_pre_1PSVJ":                       preselection & (nsvjJetsAK8 >= 1),
            "_pre_0SVJ":                        preselection & (nsvjJetsAK8 == 0),
            "_pre_1SVJ":                        preselection & (nsvjJetsAK8 == 1),
            "_pre_2SVJ":                        preselection & (nsvjJetsAK8 == 2),
            "_pre_3SVJ":                        preselection & (nsvjJetsAK8 == 3),
            "_pre_4PSVJ":                       preselection & (nsvjJetsAK8 >= 4),
            "_pre_2J_1PSVJ":                    preselection & (njetsAK8 == 2) & (nsvjJetsAK8 >= 1),
            "_pre_3J_1PSVJ":                    preselection & (njetsAK8 == 3) & (nsvjJetsAK8 >= 1),
            "_pre_4J_1PSVJ":                    preselection & (njetsAK8 == 4) & (nsvjJetsAK8 >= 1),
            "_pre_5PJ_1PSVJ":                   preselection & (njetsAK8 >= 5) & (nsvjJetsAK8 >= 1), 
            # "_qual_trg_st_0nim_0SVJ":           qualityCuts & passTrigger & stCut & (nnim == 0) & (nsvjJetsAK8 == 0),
            # "_qual_trg_st_0nim_1SVJ":           qualityCuts & passTrigger & stCut & (nnim == 0) & (nsvjJetsAK8 == 1),
            # "_qual_trg_st_0nim_2SVJ":           qualityCuts & passTrigger & stCut & (nnim == 0) & (nsvjJetsAK8 == 2),
            # "_qual_trg_st_0nim_ge1SVJ" :        qualityCuts & passTrigger & stCut & (nnim == 0) & (nsvjJetsAK8 >= 1),
            # "_qual_trg_st_0nim_ge2SVJ" :        qualityCuts & passTrigger & stCut & (nnim == 0) & (nsvjJetsAK8 >= 2),
            # "_qual_trg_st_0nim_4J_0SVJ":        qualityCuts & passTrigger & stCut & (nnim == 0) & (njetsAK8 == 4) & (nsvjJetsAK8 == 0),
            # "_qual_trg_st_0nim_4J_1SVJ":        qualityCuts & passTrigger & stCut & (nnim == 0) & (njetsAK8 == 4) & (nsvjJetsAK8 == 1),
            # "_qual_trg_st_0nim_4J_2SVJ":        qualityCuts & passTrigger & stCut & (nnim == 0) & (njetsAK8 == 4) & (nsvjJetsAK8 == 2),
            # "_qual_trg_st_0nim_4J_3SVJ":        qualityCuts & passTrigger & stCut & (nnim == 0) & (njetsAK8 == 4) & (nsvjJetsAK8 == 3),
            # "_qual_trg_st_0nim_4J_4SVJ":        qualityCuts & passTrigger & stCut & (nnim == 0) & (njetsAK8 == 4) & (nsvjJetsAK8 == 4),
            # "_qual_trg_st_ge1nim_0SVJ":         qualityCuts & passTrigger & stCut & (nnim >= 1) & (nsvjJetsAK8 == 0),
            # "_qual_trg_st_ge1nim_1SVJ":         qualityCuts & passTrigger & stCut & (nnim >= 1) & (nsvjJetsAK8 == 1),
            # "_qual_trg_st_ge1nim_2SVJ":         qualityCuts & passTrigger & stCut & (nnim >= 1) & (nsvjJetsAK8 == 2),
            # "_qual_trg_st_ge1nim_ge1SVJ":       qualityCuts & passTrigger & stCut & (nnim >= 1) & (nsvjJetsAK8 >= 1),
            # "_qual_trg_st_ge1nim_ge2SVJ":       qualityCuts & passTrigger & stCut & (nnim >= 1) & (nsvjJetsAK8 >= 2),
        }
        cuts.update(cutsWithNSVJ)
    return cuts
