import numpy as np
import awkward as ak
from . import triggerDict as tD
from . import utility as util
# dataset here is just a string that is the input of -d
# events are the tree in the ntuple

def hemVeto(ak4Jets,electrons,muons):
    ak4jHemCond = (ak4Jets.eta > -3.05) & (ak4Jets.eta < -1.35) & (ak4Jets.phi > -1.62) & (ak4Jets.phi < -0.82)
    elecHemCond = (electrons.eta > -3.05) & (electrons.eta < -1.35) & (electrons.phi > -1.62) & (electrons.phi < -0.82)
    muonHemCond = (muons.eta > -3.05) & (muons.eta < -1.35) & (muons.phi > -1.62) & (muons.phi < -0.82)
    hemEvents = ((ak.num(ak4Jets) > 0) & ak.any(ak4jHemCond,axis=1)) | ((ak.num(muons) > 0) & ak.any(muonHemCond,axis=1)) | ((ak.num(electrons) > 0) & ak.any(elecHemCond,axis=1))
    eta = ak4Jets.eta
    phi = ak4Jets.phi
    return ~hemEvents

def hemPeriodMask(dataset,events,ak4Jets,electrons,muons,hemPeriod):
    runNum = events.RunNum
    if ("2018" in dataset):
        if (hemPeriod == "PreHEM") and ("Data" in dataset):
            return runNum < 319077
        elif (hemPeriod == "PostHEM") and ("Data" in dataset):
            return (runNum >= 319077) & hemVeto(ak4Jets)
        elif hemPeriod == "PostHEM":
            return hemVeto(ak4Jets,electrons,muons)
        else:
            return np.ones(len(events),dtype=bool)
    else:
        return np.ones(len(events),dtype=bool)

def TTStitch(dataset,events):
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

def RemoveOverlap(dataset, events, yr):
    """Function to remove overlap between different dataset by removing the events with overlapping triggers. Currently the preferred dataset is JetHT.""" 
    trgDict = tD.trigDicts[yr]
    trgListMet = trgListtoInd(trgDict,tD.trgMET[yr])
    trgListJetHT = trgListtoInd(trgDict,tD.trgJetHT[yr])
    # trgListHTMHT = trgListtoInd(trgDict,tD.trgHTMHT[yr])
    passMetTrg = PassTrigger(events.TriggerPass, trgListMet)
    # passHTMHTTrg = PassTrigger(events.TriggerPass, trgListHTMHT)
    passJetHTTrg = PassTrigger(events.TriggerPass, trgListJetHT)
    
    # For single lepton control region - Main dataset is Single Muon, removing overlap from the Single Electron
    # passSingleM = PassTrigger(events.TriggerPass,trgListtoInd(tD.trigDicts[yr],tD.trgSingleM[yr]))
    # passSingleE = PassTrigger(events.TriggerPass,trgListtoInd(tD.trigDicts[yr],tD.trgSingleE[yr])) 

    OverlapDataMask = None
    # if "JetHTData" in dataset:
    #     OverlapDataMask = np.bitwise_not(passMetTrg)
    if "METData" in dataset:
        OverlapDataMask = np.bitwise_not(passJetHTTrg)
    elif "HTMHTData" in dataset:
        OverlapDataMask = np.bitwise_not(passJetHTTrg + passMetTrg)
    # Adding single lepton control region overlap removal
    # elif "SingleElectronData" in dataset:
    #     OverlapDataMask = np.bitwise_not(passSingleM)
    # elif "EGammaData" in dataset:
    #     OverlapDataMask = np.bitwise_not(passSingleM)    
    # elif "SinglePhotonData" in dataset:
    #     OverlapDataMask = np.bitwise_not(passSingleM + passSingleE)
    else:
        OverlapDataMask = np.ones(len(events),dtype=bool)
    #Preliminary checks
    # for i in range(0,len(events)):
    #     if passMetTrg[i] and passJetHTTrg[i]:
    #         print("Event no = {} | Lumi no = {} | Run no. = {} | PassMetTrg = {} | PassJetHTTrg = {} | PassHTMHT = {} | Mask = {} --- ".format(i,events.EvtNum[i], events.LumiBlockNum[i], events.RunNum[i],passMetTrg[i],passJetHTTrg[i],passHTMHTTrg[i],OverlapDataMask[i]))
  
    return OverlapDataMask

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

def ConditionMask(var, condition):
    FullMask = []
    for val in var:
        Mask = None
        Mask = (val == 0)
        FullMask.append(Mask)
        print(val)
    print("Mask = ",FullMask)

def cutList(dataset,events,vars_noCut,hemPeriod,SVJCut=True):
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
    jetIDAK8 = events.JetIDAK8
    ttStitch = TTStitch(dataset,events)

    # ConditionMask(ptdAk8,0)

    years = ["2016","2017","2018"]
    yr = 0
    for year in years:
        if year in dataset:
            yr = year

    if "Data" in dataset:
        DataMask = RemoveOverlap(dataset, events,yr)
    else: 
        DataMask = np.ones(len(events),dtype=bool)
    metFilters = METFilters(events)
    hemMask = hemPeriodMask(dataset,events,vars_noCut["jets"],vars_noCut["electrons"],vars_noCut["muons"],hemPeriod)
    # psFilter = PhiSpikeFilter(dataset,vars_noCut['jets'])
    qualityCuts = metFilters & (nl == 0) & ttStitch & DataMask & hemMask
    qualityWithLepton = metFilters & ttStitch & DataMask & hemMask
    # qualityCuts = metFilters & psFilter # NN training files
    # preselection = Preselection(qualityCuts,nl)
    # cuts to get over trigger plateau
    # metCut = met > 266
    # htCut = ht > 1280
    htCut = ht > 500
    stCut = st > 1300
    # metcut = met > 1200

    # trgPlat = metCut & htCut

    # trigger choices
    trigDict = tD.trigDicts[yr]
    trgSelection = tD.trgSelections[yr]
    trgSelectionsCR = tD.trgSelectionsCR[yr]
    trgSelectionsQCDCR = tD.trgSelectionsQCDCR[yr]
    tch_trgs =  trgListtoInd(trigDict,trgSelection)
    tch_trgs_CR =  trgListtoInd(trigDict,trgSelectionsCR)
    tch_trgs_QCDCR =  trgListtoInd(trigDict,trgSelectionsQCDCR)
    passTrigger = PassTrigger(triggerPass,tch_trgs)
    preselection = qualityCuts & passTrigger & stCut & (njetsAK8 >= 2) & (dPhiMinjAK8 <= 1.5)
    # trigger study for MCs
    # nOffMuons = vars_noCut['nOffMuons']
    passTrigger_muon = PassTrigger(triggerPass,tch_trgs_CR)
    # preselection_offLineMuons = metFilters & ttStitch & DataMask & hemMask & passTrigger_muon & (njetsAK8 >= 2) & (nOffMuons >= 1)
    # preselection_offLineMuons_tchTrg = preselection_offLineMuons & passTrigger
    cr_muon_cut                  = qualityWithLepton & passTrigger & (ncrMuons == 1)     & (nelectron == 0) & stCut & (njetsAK8 >= 2) & (dPhiMinjAK8 <= 1.5)
    cr_electron_cut              = qualityWithLepton & passTrigger & (ncrElectrons == 1) & (nmuon == 0) & stCut & (njetsAK8 >= 2) & (dPhiMinjAK8 <= 1.5)
    # cr_dphimin                   = preselection & (dPhiMinjAK8 > 1.5)
    # cuts = {
    #             ""                                          : np.ones(len(evtw),dtype=bool),
    #             "_preselec_offLineMuons"                    : preselection_offLineMuons,  
    #             "_preselection_offLineMuons_tchTrg"         : preselection_offLineMuons_tchTrg,
    #             "_preselec_offLineMuons_st"                 : preselection_offLineMuons & stCut,  
    #             "_preselection_offLineMuons_tchTrg_st"      : preselection_offLineMuons_tchTrg & stCut,
    # }

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
    else:
        # print("No problem before this")
        # Define all cuts for histo making
        cuts = {
                    # ""                          : np.ones(len(evtw),dtype=bool),
                    # "_2PJ"                      : (njetsAK8 >= 2),
                    # "_2PJ_nl"            : (njetsAK8 >= 2) & (nl == 0),
                    # # "_data_mask"            : DataMask,
                    # "_st"                : stCut,
                    # "_ht"               : htCut,
                    # "_trg"               : passTrigger,
                    # "_qual"             : qualityCuts,
                    # "_qual_ht"              : qualityCuts & htCut,
                    # "_qual_trg"             : qualityCuts & passTrigger,
                    # "_qual_st"             : qualityCuts & stCut,
                    
                    # "_qual_trg_st"              : qualityCuts & passTrigger & stCut,
                    # "_qual_trg_st_dphimin"      : qualityCuts & passTrigger & stCut & dPhiMinjAK8Cut,
                    # "_qual_trg_st_ht"           : qualityCuts & passTrigger & stCut & htCut,
                    # "_qual_trg_st_1PJ"          : qualityCuts & passTrigger & stCut & (njetsAK8 >= 1),
                    # "_qual_trg_st_2PJ"          : qualityCuts & passTrigger & stCut & (njetsAK8 >= 2),
                    # "_qual_trg_st_ht_2PJ"          : qualityCuts & passTrigger & stCut & htCut & (njetsAK8 >= 2),
                    # "_qual_trg_st_ht_2PJ_dphimin"          : qualityCuts & passTrigger & stCut & htCut & (njetsAK8 >= 2) & dPhiMinjAK8Cut,
                    # "_all_cuts_ptd"    : qualityCuts & passTrigger & stCut & htCut & (njetsAK8 >= 2) & dPhiMinjAK8Cut & (ptdAk8 == 0),
                    # "_all_cuts_girth"  : qualityCuts & passTrigger & stCut & htCut & (njetsAK8 >= 2) & dPhiMinjAK8Cut & (girthAK8 > 0.5),
                    # "_all_cuts_met"    : qualityCuts & passTrigger & stCut & htCut & (njetsAK8 >= 2) & dPhiMinjAK8Cut & (met > 1200)
                    # "_qual_st"             : qualityCuts & stCut,
                }

    
    if SVJCut == True:
        nsvjJetsAK8 = vars_noCut["nsvjJetsAK8"]
        cuts = {
                "_nocut" : np.ones(len(evtw),dtype=bool), 
                "_qual" : qualityWithLepton,
                "_qual_passTrig" : qualityWithLepton & passTrigger,
                "_qual_2PJ" : qualityWithLepton & passTrigger & (njetsAK8 >=2),
                "_qual_2PJ_st" : qualityWithLepton & passTrigger & (njetsAK8 >=2) & stCut,
                "_qual_2PJ_st_nl" : qualityWithLepton & passTrigger & (njetsAK8 >=2) & stCut & (nl == 0),
                
                "_pre_dphimin_":      preselection,
                # "_pre_dphimin_0SVJ":  preselection & (nsvjJetsAK8 == 0),
                # "_pre_dphimin_1SVJ":  preselection & (nsvjJetsAK8 == 1),
                # "_pre_dphimin_2SVJ":  preselection & (nsvjJetsAK8 == 2),
                # "_pre_dphimin_3SVJ":  preselection & (nsvjJetsAK8 == 3),
                # "_pre_dphimin_4PSVJ": preselection & (nsvjJetsAK8 >= 4),
            
                "_cr_muon_"                  : cr_muon_cut, 
                "_cr_electron_"              : cr_electron_cut,
                "_cr_muon_0SVJ"                  : cr_muon_cut & (nsvjJetsAK8 == 0), 
                "_cr_electron_0SVJ"              : cr_electron_cut & (nsvjJetsAK8 == 0),
                "_cr_muon_1SVJ"                  : cr_muon_cut & (nsvjJetsAK8 == 1), 
                "_cr_electron_1SVJ"              : cr_electron_cut & (nsvjJetsAK8 == 1),
                "_cr_muon_2SVJ"                  : cr_muon_cut & (nsvjJetsAK8 == 2), 
                "_cr_electron_2SVJ"              : cr_electron_cut & (nsvjJetsAK8 == 2),
                "_cr_muon_3SVJ"                  : cr_muon_cut & (nsvjJetsAK8 == 3), 
                "_cr_electron_3SVJ"              : cr_electron_cut & (nsvjJetsAK8 == 3),
                "_cr_muon_4PSVJ"                 : cr_muon_cut & (nsvjJetsAK8 >= 4), 
                "_cr_electron_4PSVJ"             : cr_electron_cut & (nsvjJetsAK8 >= 4),

                # ""            : np.ones(len(evtw),dtype=bool),
                # "_pre"        : preselection,
                # "_pre_1PSVJ"  : preselection & (nsvjJetsAK8 >= 1),
                # "_pre_2J"     : preselection & (njetsAK8 == 2),
                # "_pre_3J"     : preselection & (njetsAK8 == 3),
                # "_pre_4J"     : preselection & (njetsAK8 == 4),
                # "_pre_5PJ"    : preselection & (njetsAK8 >= 5),
                # "_pre_0SVJ"   : preselection & (nsvjJetsAK8 == 0),
                # "_pre_1SVJ"   : preselection & (nsvjJetsAK8 == 1),
                # "_pre_2SVJ"   : preselection & (nsvjJetsAK8 == 2),
                # "_pre_3SVJ"   : preselection & (nsvjJetsAK8 == 3),
                # "_pre_4PSVJ"   : preselection & (nsvjJetsAK8 >= 4),
                # "_pre_2J_1PSVJ"     : preselection & (njetsAK8 == 2) & (nsvjJetsAK8 >= 1),
                # "_pre_3J_1PSVJ"     : preselection & (njetsAK8 == 3) & (nsvjJetsAK8 >= 1),
                # "_pre_4J_1PSVJ"     : preselection & (njetsAK8 == 4) & (nsvjJetsAK8 >= 1),
                # "_pre_5PJ_1PSVJ"    : preselection & (njetsAK8 >= 5) & (nsvjJetsAK8 >= 1),   
                # #"_qual"             : qualityCuts,
                #"_qual_met"         : qualityCuts & metCut,
                #"_qual_ht"          : qualityCuts & htCut,
                #"_qual_st"          : qualityCuts & stCut,
                #"_qual_trg"         : qualityCuts & passTrigger,
                #"_qual_trg_met"     : qualityCuts & passTrigger & metCut,
                #"_qual_trg_ht"      : qualityCuts & passTrigger & htCut,
                # "_qual_trg_st"             : qualityCuts & passTrigger & stCut,
                # "_qual_trg_st_0nim"        : qualityCuts & passTrigger & stCut & (nnim == 0),
                # "_qual_trg_st_0nim_0SVJ"   : qualityCuts & passTrigger & stCut & (nnim == 0) & (nsvjJetsAK8 == 0),
                # "_qual_trg_st_0nim_1SVJ"   : qualityCuts & passTrigger & stCut & (nnim == 0) & (nsvjJetsAK8 == 1),
                # "_qual_trg_st_0nim_2SVJ"   : qualityCuts & passTrigger & stCut & (nnim == 0) & (nsvjJetsAK8 == 2),
                # "_qual_trg_st_0nim_ge1SVJ" : qualityCuts & passTrigger & stCut & (nnim == 0) & (nsvjJetsAK8 >= 1),
                # "_qual_trg_st_0nim_ge2SVJ" : qualityCuts & passTrigger & stCut & (nnim == 0) & (nsvjJetsAK8 >= 2),
                # "_qual_trg_st_0nim_0J"     : qualityCuts & passTrigger & stCut & (nnim == 0) & (njetsAK8 == 0),
                # "_qual_trg_st_0nim_1J"     : qualityCuts & passTrigger & stCut & (nnim == 0) & (njetsAK8 == 1),
                # "_qual_trg_st_0nim_2J"     : qualityCuts & passTrigger & stCut & (nnim == 0) & (njetsAK8 == 2),

                # "_qual_trg_st_0nim_4J"     : qualityCuts & passTrigger & stCut & (nnim == 0) & (njetsAK8 == 4),
                # "_qual_trg_st_0nim_4J_0SVJ": qualityCuts & passTrigger & stCut & (nnim == 0) & (njetsAK8 == 4) & (nsvjJetsAK8 == 0),
                # "_qual_trg_st_0nim_4J_1SVJ": qualityCuts & passTrigger & stCut & (nnim == 0) & (njetsAK8 == 4) & (nsvjJetsAK8 == 1),
                # "_qual_trg_st_0nim_4J_2SVJ": qualityCuts & passTrigger & stCut & (nnim == 0) & (njetsAK8 == 4) & (nsvjJetsAK8 == 2),
                # "_qual_trg_st_0nim_4J_3SVJ": qualityCuts & passTrigger & stCut & (nnim == 0) & (njetsAK8 == 4) & (nsvjJetsAK8 == 3),
                # "_qual_trg_st_0nim_4J_4SVJ": qualityCuts & passTrigger & stCut & (nnim == 0) & (njetsAK8 == 4) & (nsvjJetsAK8 == 4),

                # "_qual_trg_st_0nim_ge1J"   : qualityCuts & passTrigger & stCut & (nnim == 0) & (njetsAK8 >= 1),
                # "_qual_trg_st_0nim_ge2J"   : qualityCuts & passTrigger & stCut & (nnim == 0) & (njetsAK8 >= 2),

                # "_qual_trg_st_ge1nim"        : qualityCuts & passTrigger & stCut & (nnim >= 1),
                # "_qual_trg_st_ge1nim_0SVJ"   : qualityCuts & passTrigger & stCut & (nnim >= 1) & (nsvjJetsAK8 == 0),
                # "_qual_trg_st_ge1nim_1SVJ"   : qualityCuts & passTrigger & stCut & (nnim >= 1) & (nsvjJetsAK8 == 1),
                # "_qual_trg_st_ge1nim_2SVJ"   : qualityCuts & passTrigger & stCut & (nnim >= 1) & (nsvjJetsAK8 == 2),
                # "_qual_trg_st_ge1nim_ge1SVJ" : qualityCuts & passTrigger & stCut & (nnim >= 1) & (nsvjJetsAK8 >= 1),
                # "_qual_trg_st_ge1nim_ge2SVJ" : qualityCuts & passTrigger & stCut & (nnim >= 1) & (nsvjJetsAK8 >= 2),
                # "_qual_trg_st_ge1nim_0J"     : qualityCuts & passTrigger & stCut & (nnim >= 1) & (njetsAK8 == 0),
                # "_qual_trg_st_ge1nim_1J"     : qualityCuts & passTrigger & stCut & (nnim >= 1) & (njetsAK8 == 1),
                # "_qual_trg_st_ge1nim_2J"     : qualityCuts & passTrigger & stCut & (nnim >= 1) & (njetsAK8 == 2),
                # "_qual_trg_st_ge1nim_ge1J"   : qualityCuts & passTrigger & stCut & (nnim >= 1) & (njetsAK8 >= 1),
                # "_qual_trg_st_ge1nim_ge2J"   : qualityCuts & passTrigger & stCut & (nnim >= 1) & (njetsAK8 >= 2),
                #"_metfilter_0l_1nim_trgQCDCR" : metFilters & (nl == 0) & (nnim == 1) & PassTrigger(triggerPass,tch_trgs_QCDCR),
        }

    return cuts
