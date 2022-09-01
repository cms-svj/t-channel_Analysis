import numpy as np
import awkward as ak
# from mt2 import mt2
from . import objects as ob
from itertools import combinations

def awkwardReshape(akArray,npArray):
    if len(akArray) == 0:
        return ak.Array([])
    else:
        return ak.broadcast_arrays(akArray.pt,1.0)[1] * npArray

def arrayConcatenate(array1,array2):
    if ak.any(array1) != True:
        return [array2]
    else:
        return ak.concatenate((array1, [array2]),axis=0)

# the two functions below return infinity when the event doesn't have the required
# number of jets or return the correct value for the jet variable
def jetVar_i(var,i):
    paddedVar = ak.fill_none(ak.pad_none(var,i+1),np.Inf)
    return paddedVar[:,i]

# convert phi values into spherical coordinate?
def phi_(x,y):
    phi = np.arctan2(y,x)
    return ak.where(phi < 0, phi + 2*np.pi, phi)

def deltaPhi(jetphiL,metPhiL):
    phi1 = phi_( np.cos(jetphiL), np.sin(jetphiL) )
    phi2 = phi_( np.cos(metPhiL), np.sin(metPhiL) )
    dphi = phi1 - phi2
    dphi_edited = ak.where(dphi < -np.pi, dphi + 2*np.pi, dphi)
    dphi_edited = ak.where(dphi_edited > np.pi, dphi_edited - 2*np.pi, dphi_edited)
    return abs(dphi_edited)

def deltaEta(eta0,eta1):
    return abs(eta0 - eta1)

def delta_R(eta0,eta1,phi0,phi1):
    dp = deltaPhi(phi0,phi1)
    deta = deltaEta(eta0,eta1)
    deltaR2 = deta * deta + dp * dp
    return np.sqrt(deltaR2)

def tauRatio(tau_a,tau_b,i):
    Ji_tau_a = jetVar_i(tau_a,i)
    Ji_tau_b = jetVar_i(tau_b,i)
    Ji_tau_ab = Ji_tau_a/Ji_tau_b
    return Ji_tau_ab

## Need to update mt2 code.
# def lorentzVector(pt,eta,phi,mass,i):
#     px = pt[i]*np.cos(phi[i])
#     py = pt[i]*np.sin(phi[i])
#     pz = pt[i]*np.sinh(eta[i])
#     p2 = px**2 + py**2 + pz**2
#     m2 = mass[i]**2
#     energy = np.sqrt(m2+p2)
#     return np.array([energy,px,py,pz])
#
# def mass4vec(m4vec):
#     E = m4vec[0]
#     px = m4vec[1]
#     py = m4vec[2]
#     pz = m4vec[3]
#     return np.sqrt(E**2 - (px**2 + py**2 + pz**2))
#
# def eta4vec(m4vec):
#     px = m4vec[1]
#     py = m4vec[2]
#     pz = m4vec[3]
#     pt = np.sqrt(px**2 + py**2)
#     return np.arcsinh(pz/pt)
#
# def phi4vec(m4vec):
#     px = m4vec[1]
#     py = m4vec[2]
#     pt = np.sqrt(px**2 + py**2)
#     return np.arcsin(py/pt)
#
# def M_2J(j1,j2):
#     totJets = j1+j2
#     return mass4vec(totJets)
#
# def MT2Cal(FDjet0,FSMjet0,FDjet1,FSMjet1,met,metPhi):
#     Fjet0 = FDjet0 + FSMjet0
#     Fjet1 = FDjet1 + FSMjet1
#     METx = met*np.cos(metPhi)
#     METy = met*np.sin(metPhi)
#     MT2v = mt2(
#     mass4vec(Fjet0), Fjet0[1], Fjet0[2],
#     mass4vec(Fjet1), Fjet1[1], Fjet1[2],
#     METx, METy, 0.0, 0.0, 0
#     )
#     return MT2v
#
# def f4msmCom(pt,eta,phi,mass,met,metPhi,cut):
#     MT2 = 0
#     if len(pt) >= 4:
#         List4jets_3Com = [[0,1,2,3],[0,2,1,3],[0,3,1,2]]
#         diffList = []
#         jetList = []
#         for c in List4jets_3Com:
#             jc1 = lorentzVector(pt,eta,phi,mass,c[0])
#             jc2 = lorentzVector(pt,eta,phi,mass,c[1])
#             jc3 = lorentzVector(pt,eta,phi,mass,c[2])
#             jc4 = lorentzVector(pt,eta,phi,mass,c[3])
#             jetList.append([jc1,jc2,jc3,jc4])
#             diffList.append(abs(M_2J(jc1,jc2) - M_2J(jc3,jc4)))
#         comIndex = np.argmin(diffList)
#         msmJet = jetList[comIndex]
#         # applying angular cut
#         dEtaCut = 1.8
#         dPhiCut = 0.9
#         dRCutLow = 1.5
#         dRCutHigh = 4.0
#         dPhiMETCut = 1.5
#
#         if cut == "dEta":
#             if deltaEta(eta4vec(msmJet[0]),eta4vec(msmJet[1])) < dEtaCut and deltaEta(eta4vec(msmJet[2]),eta4vec(msmJet[3])) < dEtaCut:
#                 MT2 = MT2Cal(msmJet[0],msmJet[1],msmJet[2],msmJet[3],met,metPhi)
#         elif cut == "dPhi":
#             if deltaPhiji(phi4vec(msmJet[0]),phi4vec(msmJet[1])) < dPhiCut and deltaPhiji(phi4vec(msmJet[2]),phi4vec(msmJet[3])) > dPhiCut:
#                 MT2 = MT2Cal(msmJet[0],msmJet[1],msmJet[2],msmJet[3],met,metPhi)
#         elif cut == "dR":
#             if (dRCutLow < delta_R(eta4vec(msmJet[0]),eta4vec(msmJet[1]),phi4vec(msmJet[0]),phi4vec(msmJet[1])) < dRCutHigh) and (dRCutLow < delta_R(eta4vec(msmJet[2]),eta4vec(msmJet[3]),phi4vec(msmJet[2]),phi4vec(msmJet[3])) < dRCutHigh):
#                 MT2 = MT2Cal(msmJet[0],msmJet[1],msmJet[2],msmJet[3],met,metPhi)
#         elif cut == "":
#             MT2 = MT2Cal(msmJet[0],msmJet[1],msmJet[2],msmJet[3],met,metPhi)
#     return MT2
#
# def f4msmCom_vec(pt,eta,phi,mass,met,metPhi,cut):
#     vfunc = np.vectorize(lambda pt,eta,phi,mass,met,metPhi,cut: f4msmCom(pt,eta,phi,mass,met,metPhi,cut))
#     return vfunc(pt,eta,phi,mass,met,metPhi,cut)

def decode(hvCat,nlabel,clabel,catList):
    if hvCat >= nlabel:
        hvCat -= nlabel
        catList.append(clabel)
    return hvCat

def tch_hvCat_decode(hvCat):
    catList = []
    hvCat = decode(hvCat,16,"QsM",catList)
    hvCat = decode(hvCat,8,"QdM",catList)
    hvCat = decode(hvCat,4,"Gd",catList)
    hvCat = decode(hvCat,2,"Qd",catList)
    hvCat = decode(hvCat,1,"stableD",catList)
    return catList

def baselineVar(dataset,events,scaleFactor):
    varVal = {}
    dataKeys = ["HTMHT","JetHT","MET","SingleElectron","SingleMuon","SinglePhoton","EGamma"]
    isData = False
    for dKey in dataKeys:
        if dKey in dataset:
            isData = True
            break
    evtw = np.ones(len(events))
    if not isData:
        luminosity = 59692.692 # 2018 lumi
        if "2016" in dataset:
            luminosity = 35921.036
        elif "2017" in dataset:
            luminosity = 41521.331
        evtw = luminosity*events.Weight*scaleFactor
    eCounter = np.where(evtw >= 0, 1, -1)
    obj = ob.Objects(events)
    jets = obj.goodJets()
    bjets = obj.goodBJets(dataset,jets)
    fjets = obj.goodFatJets()
    electrons = obj.goodElectrons()
    muons = obj.goodMuons()
    nonIsoMuons = obj.nonIsoMuons()
    met = events.MET
    metPhi = events.METPhi
    mtAK8 = events.MT_AK8
    ht = ak.sum(jets.pt,axis=1)
    st = ht + met
    # AK4 Jet Variables
    jetPhi = jets.phi
    dPhij = deltaPhi(jetPhi,metPhi)
    dPhiMinj = ak.min(dPhij,axis=1,mask_identity=False)
    # AK8 Jet Variables
    jetAK8Phi = fjets.phi
    dPhijAK8 = deltaPhi(jetAK8Phi,metPhi)
    dPhiMinjAK8 = ak.min(dPhijAK8,axis=1,mask_identity=False)

    if len(bjets) > 0:
        nBJets = ak.num(bjets)
    else:
        nBJets = np.zeros(len(evtw))

    varVal['jets'] = jets
    varVal['bjets'] = bjets
    varVal['fjets'] = fjets
    varVal['electrons'] = electrons
    varVal['muons'] = muons
    varVal['nonIsoMuons'] = nonIsoMuons
    varVal['evtw'] = evtw
    varVal['nl'] = (ak.num(electrons) + ak.num(muons))
    varVal['nnim'] = ak.num(nonIsoMuons)
    varVal['njets'] = ak.num(jets)
    varVal['njetsAK8'] = ak.num(fjets)
    varVal['nb'] = nBJets
    varVal['met'] = met
    varVal['metPhi'] = metPhi
    varVal['mT'] = mtAK8
    varVal['ht'] = ht
    varVal['st'] = st
    varVal['METrHT_pt30'] = met/ht
    varVal['METrST_pt30'] = met/st
    varVal['dPhiMinjMET'] = dPhiMinj
    varVal['dPhiMinjMETAK8'] = dPhiMinjAK8
    return varVal

def jConstVarGetter(dataset,events,varVal,cut):
    evtw = varVal["evtw"][cut]
    fjets = varVal["fjets"][cut]
    ## GenJetsAK8_hvCategory is only present in the signal samples, not the V17 background
    jetCats = []
    fjw = awkwardReshape(fjets,evtw)
    evtNum = events.EvtNum
    fjEvtNum = awkwardReshape(fjets,evtNum)
    bkgKeys = ["QCD","TTJets","WJets","ZJets"]
    isSignal = False
    if "mMed" in dataset:
        isSignal = True
    if isSignal:
        jetsAK8GenInd = fjets.genIndex
        for gji in range(len(jetsAK8GenInd)):
            genInd = jetsAK8GenInd[gji]
            GenJetsAK8 = events.GenJetsAK8
            gfjets = GenJetsAK8[GenJetsAK8.pt > 170 & (abs(GenJetsAK8.eta) < 5.0)]
            genCat = gfjets.hvCategory[gji]
            if (len(genCat) > 0) and (len(genInd) > 0):
                if np.max(genInd) < len(genCat):
                    jetCats.append(list(genCat[genInd]))
                else:
                    jetCats.append([-1]*len(genInd))
            else:
                jetCats.append([-1]*len(genInd))
        jetCats = ak.Array(jetCats)
    else:
        jetCats = awkwardReshape(fjets,np.ones(len(events))*-1)

    jetConstituents = events.JetsConstituents
    JetsAK8_constituentsIndex = fjets.constituentsIndex
    jCstPt = jetConstituents.pt
    jCst4vec = {}
    jCstVar = {}
    jCst4vec["jCstPt"] = jCstPt
    jCst4vec["jCstEta"] = jetConstituents.eta
    jCst4vec["jCstPhi"] = jetConstituents.phi
    jCst4vec["jCstEnergy"] = jetConstituents.energy
    jCst4vec["jCstPdgId"] = jetConstituents.PdgId
    jCstVar["jCstPtAK8"] = [fjets.pt,[False],np.array([])]
    jCstVar["jCstAxismajorAK8"] = [fjets.axismajor,[False],np.array([])]
    jCstVar["jCstAxisminorAK8"] = [fjets.axisminor,[False],np.array([])]
    # jCstVar["jCstChEMEFractAK8"] = [fjets.chargedEmEnergyFraction,[False],np.array([])]
    # jCstVar["jCstChHadEFractAK8"] = [fjets.chargedHadronEnergyFraction,[False],np.array([])]
    # jCstVar["jCstChHadMultAK8"] = [fjets.chargedHadronMultiplicity,[False],np.array([])]
    # jCstVar["jCstChMultAK8"] = [fjets.chargedMultiplicity,[False],np.array([])]
    jCstVar["jCstdoubleBDiscriminatorAK8"] = [fjets.doubleBDiscriminator,[False],np.array([])]
    # jCstVar["jCstecfN2b1AK8"] = [fjets.ecfN2b1,[False],np.array([])]
    # jCstVar["jCstecfN2b2AK8"] = [fjets.ecfN2b2,[False],np.array([])]
    # jCstVar["jCstecfN3b1AK8"] = [fjets.ecfN3b1,[False],np.array([])]
    # jCstVar["jCstecfN3b2AK8"] = [fjets.ecfN3b2,[False],np.array([])]
    # jCstVar["jCstEleEFractAK8"] = [fjets.electronEnergyFraction,[False],np.array([])]
    # jCstVar["jCstEleMultAK8"] = [fjets.electronMultiplicity,[False],np.array([])]
    # jCstVar["jCstGirthAK8"] = [fjets.girth,[False],np.array([])]
    # jCstVar["jCstHfEMEFractAK8"] = [fjets.hfEMEnergyFraction,[False],np.array([])]
    # jCstVar["jCstHfHadEFractAK8"] = [fjets.hfHadronEnergyFraction,[False],np.array([])]
    # jCstVar["jCstMultAK8"] = [fjets.multiplicity,[False],np.array([])]
    # jCstVar["jCstMuEFractAK8"] = [fjets.muonEnergyFraction,[False],np.array([])]
    # jCstVar["jCstMuMultAK8"] = [fjets.muonMultiplicity,[False],np.array([])]
    # jCstVar["jCstNeuEmEFractAK8"] = [fjets.neutralEmEnergyFraction,[False],np.array([])]
    # jCstVar["jCstNeuHadEFractAK8"] = [fjets.neutralHadronEnergyFraction,[False],np.array([])]
    # jCstVar["jCstNeuHadMultAK8"] = [fjets.neutralHadronMultiplicity,[False],np.array([])]
    # jCstVar["jCstNeuMultAK8"] = [fjets.neutralMultiplicity,[False],np.array([])]
    jCstVar["jCstTau1AK8"] = [fjets.NsubjettinessTau1,[False],np.array([])]
    jCstVar["jCstTau2AK8"] = [fjets.NsubjettinessTau2,[False],np.array([])]
    jCstVar["jCstTau3AK8"] = [fjets.NsubjettinessTau3,[False],np.array([])]
    jCstVar["jCstNumBhadronsAK8"] = [fjets.NumBhadrons,[False],np.array([])]
    jCstVar["jCstNumChadronsAK8"] = [fjets.NumChadrons,[False],np.array([])]
    # jCstVar["jCstPhoEFractAK8"] = [fjets.photonEnergyFraction,[False],np.array([])]
    # jCstVar["jCstPhoMultAK8"] = [fjets.photonMultiplicity,[False],np.array([])]
    jCstVar["jCstPtDAK8"] = [fjets.ptD,[False],np.array([])]
    jCstVar["jCstSoftDropMassAK8"] = [fjets.softDropMass,[False],np.array([])]
    jCstVar["jCsthvCategory"] = [jetCats,[False],np.array([])]
    jCstVar["jCstWeightAK8"] = [fjw,[False],np.array([])]
    jCstVar["jCstEvtNum"] = [fjEvtNum,[False],np.array([])]
    jCstVar["jCstJNum"] = [ak.local_index(fjw),[False],np.array([])]

    # looping over events
    for i in range(len(JetsAK8_constituentsIndex)):
        jcPt = jCstPt[i]
        JetsAK8_constituentsIndexPerEvent = JetsAK8_constituentsIndex[i]
        # if an event doesn't have any jet
        if len(JetsAK8_constituentsIndexPerEvent) == 0:
            for key,details in jCstVar.items():
                jCstVar[key][1] = arrayConcatenate(jCstVar[key][1],[np.inf])
            continue
        # looping over jets in each event
        for j in range(len(JetsAK8_constituentsIndexPerEvent)):
            consInd = ak.to_numpy(JetsAK8_constituentsIndexPerEvent[j])
            for key,details in jCstVar.items():
                if j == 0:
                    jetConstIndFlatten = np.zeros(len(jcPt)) + np.inf
                    jetConstIndFlatten.flat[consInd] = details[0][i][j]
                    details[2] = jetConstIndFlatten
                else:
                    details[2].flat[consInd] = details[0][i][j]
        for key in jCstVar.keys():
            jCstVar[key][1] = arrayConcatenate(jCstVar[key][1],jCstVar[key][2])
            jCstVar[key][2] = np.array([])
    return jCst4vec,jCstVar

def varGetter(dataset,events,varVal,jNVar=False):
    jets = varVal['jets'] 
    bjets = varVal['bjets'] 
    fjets = varVal['fjets'] 
    electrons = varVal['electrons'] 
    muons = varVal['muons'] 
    nonIsoMuons = varVal['nonIsoMuons'] 
    evtw = varVal['evtw'] 
    nBJets = varVal['nb']
    met = varVal['met']
    metPhi = varVal['metPhi']
    mtAK8 = varVal['mT']
    ht = varVal['ht']
    dPhiMinj = varVal['dPhiMinjMET']
    dPhiMinjAK8 = varVal['dPhiMinjMETAK8']

    eCounter = np.where(evtw >= 0, 1, -1)
    jetAK8Eta = fjets.eta
    jetAK8Phi = fjets.phi
    j1_etaAK8 = jetVar_i(jetAK8Eta,0)
    j2_etaAK8 = jetVar_i(jetAK8Eta,1)
    j1_phiAK8 = jetVar_i(jetAK8Phi,0)
    j2_phiAK8 = jetVar_i(jetAK8Phi,1)

    ## GenJetsAK8_hvCategory is only present in the signal samples, not the background
    jetCats = []
    bkgKeys = ["QCD","TTJets","WJets","ZJets"]
    isSignal = False
    if "mMed" in dataset:
        isSignal = True
    if isSignal:
        jetsAK8GenInd = fjets.genIndex
        for gji in range(len(jetsAK8GenInd)):
            genInd = jetsAK8GenInd[gji]
            GenJetsAK8 = events.GenJetsAK8
            gfjets = GenJetsAK8[GenJetsAK8.pt > 170 & (abs(GenJetsAK8.eta) < 5.0)]
            genCat = gfjets.hvCategory[gji]
            if (len(genCat) > 0) and (len(genInd) > 0):
                if np.max(genInd) < len(genCat):
                    jetCats.append(list(genCat[genInd]))
                else:
                    jetCats.append([-1]*len(genInd))
            else:
                jetCats.append([-1]*len(genInd))
        jetCats = ak.Array(jetCats)
        varVal['JetsAK8_hvCategory'] = jetCats
    else:
        varVal['JetsAK8_hvCategory'] = awkwardReshape(fjets,np.ones(len(evtw))*-1)
        
    ew = awkwardReshape(electrons,evtw)
    mw = awkwardReshape(muons,evtw)
    nimw = awkwardReshape(nonIsoMuons,evtw)
    jw = awkwardReshape(jets,evtw)
    fjw = awkwardReshape(fjets,evtw)

    # AK4 Jet Variables
    jetPhi = jets.phi
    jetEta = jets.eta
    j1_eta = jetVar_i(jetEta,0)
    j2_eta = jetVar_i(jetEta,1)
    j1_phi = jetVar_i(jetPhi,0)
    j2_phi = jetVar_i(jetPhi,1)
    dPhij1 = deltaPhi(j1_phi,metPhi)
    dPhij2 = deltaPhi(j2_phi,metPhi)
    dPhij1rdPhij2 = dPhij1/dPhij2
    dPhiMinj = ak.min(deltaPhi(jetPhi,metPhi),axis=1,mask_identity=False)
    dEtaj12 = deltaEta(j1_eta,j2_eta)
    deltaR12j = delta_R(j1_eta,j2_eta,j1_phi,j2_phi)

    # AK8 Jet Variables
    jetAK8pT = fjets.pt
    jetAK8Phi = fjets.phi
    jetAK8Eta = fjets.eta
    jetAK8M = fjets.mass
    j1_etaAK8 = jetVar_i(jetAK8Eta,0)
    j2_etaAK8 = jetVar_i(jetAK8Eta,1)
    j1_phiAK8 = jetVar_i(jetAK8Phi,0)
    j2_phiAK8 = jetVar_i(jetAK8Phi,1)
    dPhij1AK8 = deltaPhi(j1_phiAK8,metPhi)
    dPhij2AK8 = deltaPhi(j2_phiAK8,metPhi)
    dPhij1rdPhij2AK8 = dPhij1AK8/dPhij2AK8
    dPhijAK8 = deltaPhi(jetAK8Phi,metPhi)
    dPhiMinjAK8 = ak.min(dPhijAK8,axis=1,mask_identity=False)
    dEtaj12AK8 = deltaEta(j1_etaAK8,j2_etaAK8)
    deltaR12jAK8 = delta_R(j1_etaAK8,j2_etaAK8,j1_phiAK8,j2_phiAK8)
    tau1 = fjets.NsubjettinessTau1
    tau2 = fjets.NsubjettinessTau2
    tau3 = fjets.NsubjettinessTau3
    J_tau21 = tau2/tau1
    J_tau32 = tau3/tau2
    J1_tau21 = tauRatio(tau2,tau1,0)
    J1_tau32 = tauRatio(tau3,tau2,0)
    J2_tau21 = tauRatio(tau2,tau1,1)
    J2_tau32 = tauRatio(tau3,tau2,1)

    varVal['eCounter'] = eCounter
    varVal['jw'] = jw
    varVal['fjw'] = fjw
    varVal['ew'] = ew
    varVal['mw'] = mw
    varVal['nimw'] = nimw
    # varVal['madHT'] = madHT
    varVal['jPt'] = jets.pt
    varVal['jEta'] = jetEta
    varVal['jPhi'] = jetPhi
    varVal['jAxismajor'] = jets.axismajor
    varVal['jAxisminor'] = jets.axisminor
    varVal['jPtD'] = jets.ptD
    varVal['dPhiMinjMET'] = dPhiMinj
    varVal['jPtAK8'] = fjets.pt
    varVal['jEtaAK8'] = jetAK8Eta
    varVal['jPhiAK8'] = jetAK8Phi
    varVal['jAxismajorAK8'] = fjets.axismajor
    varVal['jAxisminorAK8'] = fjets.axisminor
    varVal['jChEMEFractAK8'] = fjets.chargedEmEnergyFraction
    varVal['jChHadEFractAK8'] = fjets.chargedHadronEnergyFraction
    varVal['jChHadMultAK8'] = fjets.chargedHadronMultiplicity
    varVal['jChMultAK8'] = fjets.chargedMultiplicity
    varVal['jdoubleBDiscriminatorAK8'] = fjets.doubleBDiscriminator
    varVal['jecfN2b1AK8'] = fjets.ecfN2b1
    varVal['jecfN2b2AK8'] = fjets.ecfN2b2
    varVal['jecfN3b1AK8'] = fjets.ecfN3b1
    varVal['jecfN3b2AK8'] = fjets.ecfN3b2
    varVal['jEleEFractAK8'] = fjets.electronEnergyFraction
    varVal['jEleMultAK8'] = fjets.electronMultiplicity
    varVal['jGirthAK8'] = fjets.girth
    varVal['jHfEMEFractAK8'] = fjets.hfEMEnergyFraction
    varVal['jHfHadEFractAK8'] = fjets.hfHadronEnergyFraction
    varVal['jMultAK8'] = fjets.multiplicity
    varVal['jMuEFractAK8'] = fjets.muonEnergyFraction
    varVal['jMuMultAK8'] = fjets.muonMultiplicity
    varVal['jNeuEmEFractAK8'] = fjets.neutralEmEnergyFraction
    varVal['jNeuHadEFractAK8'] = fjets.neutralHadronEnergyFraction
    varVal['jNeuHadMultAK8'] = fjets.neutralHadronMultiplicity
    varVal['jNeuMultAK8'] = fjets.neutralMultiplicity
    varVal['jTau1AK8'] = tau1
    varVal['jTau2AK8'] = tau2
    varVal['jTau3AK8'] = tau3
    varVal['jTau21AK8'] = J_tau21
    varVal['jTau32AK8'] = J_tau32
    varVal['jNumBhadronsAK8'] = fjets.NumBhadrons
    varVal['jNumChadronsAK8'] = fjets.NumChadrons
    varVal['jPhoEFractAK8'] = fjets.photonEnergyFraction
    varVal['jPhoMultAK8'] = fjets.photonMultiplicity
    varVal['jPtDAK8'] = fjets.ptD
    varVal['jSoftDropMassAK8'] = fjets.softDropMass
    varVal['dPhijMETAK8'] = dPhijAK8
    varVal['dPhiMinjMETAK8'] = dPhiMinjAK8
    varVal['dEtaj12AK8'] = dEtaj12AK8
    varVal['dRJ12AK8'] = deltaR12jAK8
    varVal['dPhij1rdPhij2AK8'] = dPhij1rdPhij2AK8
    varVal['electronsIso'] = electrons.iso
    varVal['muonsIso'] = muons.iso
    varVal['nonIsoMuonsPt'] = nonIsoMuons.pt
    varVal['nonIsoMuonsIso'] = nonIsoMuons.iso
    if jNVar:
        # preparing histograms for jN variables
        maxN = 4
        for i in range(maxN):
            varVal['j{}Pt'.format(i+1)] = jetVar_i(jets.pt,i)
            varVal['j{}Eta'.format(i+1)] = jetVar_i(jetEta,i)
            varVal['j{}Phi'.format(i+1)] = jetVar_i(jetPhi,i)
            varVal['j{}Axismajor'.format(i+1)] = jetVar_i(jets.axismajor,i)
            varVal['j{}Axisminor'.format(i+1)] = jetVar_i(jets.axisminor,i)
            varVal['j{}PtD'.format(i+1)] = jetVar_i(jets.ptD,i)
            varVal['dPhij{}MET'.format(i+1)] = deltaPhi(jetVar_i(jetPhi,i),metPhi)
            varVal['j{}PtAK8'.format(i+1)] = jetVar_i(fjets.pt,i)
            varVal['j{}EtaAK8'.format(i+1)] = jetVar_i(jetAK8Eta,i)
            varVal['j{}PhiAK8'.format(i+1)] = jetVar_i(jetAK8Phi,i)
            varVal['j{}AxismajorAK8'.format(i+1)] = jetVar_i(fjets.axismajor,i)
            varVal['j{}AxisminorAK8'.format(i+1)] = jetVar_i(fjets.axisminor,i)
            varVal['j{}GirthAK8'.format(i+1)] = jetVar_i(fjets.girth,i)
            varVal['j{}PtDAK8'.format(i+1)] = jetVar_i(fjets.ptD,i)
            varVal['j{}Tau1AK8'.format(i+1)] = jetVar_i(tau1,i)
            varVal['j{}Tau2AK8'.format(i+1)] = jetVar_i(tau2,i)
            varVal['j{}Tau3AK8'.format(i+1)] = jetVar_i(tau3,i)
            varVal['j{}Tau21AK8'.format(i+1)] = tauRatio(tau2,tau1,i)
            varVal['j{}Tau32AK8'.format(i+1)] = tauRatio(tau3,tau2,i)
            varVal['j{}SoftDropMassAK8'.format(i+1)] = jetVar_i(fjets.softDropMass,i)
            varVal['dPhij{}METAK8'.format(i+1)] = deltaPhi(jetVar_i(jetAK8Phi,i),metPhi)
        allComs = list(combinations(range(maxN),2))
        for com in allComs:
            j1 = com[0]
            j2 = com[1]
            j1_eta = jetVar_i(jetEta,j1)
            j2_eta = jetVar_i(jetEta,j2)
            j1_phi = jetVar_i(jetPhi,j1)
            j2_phi = jetVar_i(jetPhi,j2)
            dPhij1 = deltaPhi(j1_phi,metPhi)
            dPhij2 = deltaPhi(j2_phi,metPhi)
            j1_etaAK8 = jetVar_i(jetAK8Eta,j1)
            j2_etaAK8 = jetVar_i(jetAK8Eta,j2)
            j1_phiAK8 = jetVar_i(jetAK8Phi,j1)
            j2_phiAK8 = jetVar_i(jetAK8Phi,j2)
            dPhij1AK8 = deltaPhi(j1_phiAK8,metPhi)
            dPhij2AK8 = deltaPhi(j2_phiAK8,metPhi)
            varVal['dEtaj{}{}'.format(j1+1,j2+1)] = deltaEta(j1_eta,j2_eta)
            varVal['dPhij{}{}'.format(j1+1,j2+1)] = deltaPhi(j1_phi,j2_phi)
            varVal['dRj{}{}'.format(j1+1,j2+1)] = delta_R(j1_eta,j2_eta,j1_phi,j2_phi)
            varVal['dPhij{}rdPhij{}'.format(j1+1,j2+1)] = dPhij1/dPhij2
            varVal['dEtaj{}{}AK8'.format(j1+1,j2+1)] = deltaEta(j1_etaAK8,j2_etaAK8)
            varVal['dPhij{}{}AK8'.format(j1+1,j2+1)] = deltaPhi(j1_phiAK8,j2_phiAK8)
            varVal['dRj{}{}AK8'.format(j1+1,j2+1)] = delta_R(j1_etaAK8,j2_etaAK8,j1_phiAK8,j2_phiAK8)
            varVal['dPhij{}rdPhij{}AK8'.format(j1+1,j2+1)] = dPhij1AK8/dPhij2AK8
    # varVal['GenJetsAK8_hvCategory'] = GenJetsAK8_hvCategory.flatten()
    # varVal['mT2_f4_msm'] = f4msmCom_vec(jetAK8pT,jetAK8Eta,jetAK8Phi,jetAK8M,met,metPhi,"")
    # varVal['mT2_f4_msm_dEta'] = f4msmCom_vec(jetAK8pT,jetAK8Eta,jetAK8Phi,jetAK8M,met,metPhi,"dEta")
    # varVal['mT2_f4_msm_dPhi'] = f4msmCom_vec(jetAK8pT,jetAK8Eta,jetAK8Phi,jetAK8M,met,metPhi,"dPhi")
    # varVal['mT2_f4_msm_dR'] = f4msmCom_vec(jetAK8pT,jetAK8Eta,jetAK8Phi,jetAK8M,met,metPhi,"dR")
    # varVal['GenJetsAK8_darkPtFrac'] = GenJetsAK8_darkPtFrac.flatten()
    # varVal['GenMT2_AK8'] = GenMT2_AK8
