import numpy as np
import awkward as ak
from . import objects as ob

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
    return np.nan_to_num(abs(dphi_edited),nan=np.Inf)

def varGetter(dataset,events,scaleFactor):
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
    # gfjets = obj.goodGenFatJets()
    electrons = obj.goodElectrons()
    muons = obj.goodMuons()
    nonIsoMuons = obj.nonIsoMuons()

    ## GenJetsAK8_hvCategory is only present in the signal samples, not the V17 background
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
    else:
        jetCats = awkwardReshape(fjets,np.ones(len(evtw))*-1)

    met = events.MET
    metPhi = events.METPhi
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
    jCstVar["jCstAxismajorAK8"] = [fjets.axismajor,[False],np.array([])]
    jCstVar["jCstAxisminorAK8"] = [fjets.axisminor,[False],np.array([])]
    jCstVar["jCstChEMEFractAK8"] = [fjets.chargedEmEnergyFraction,[False],np.array([])]
    jCstVar["jCstChHadEFractAK8"] = [fjets.chargedHadronEnergyFraction,[False],np.array([])]
    jCstVar["jCstChHadMultAK8"] = [fjets.chargedHadronMultiplicity,[False],np.array([])]
    jCstVar["jCstChMultAK8"] = [fjets.chargedMultiplicity,[False],np.array([])]
    jCstVar["jCstdoubleBDiscriminatorAK8"] = [fjets.doubleBDiscriminator,[False],np.array([])]
    jCstVar["jCstecfN2b1AK8"] = [fjets.ecfN2b1,[False],np.array([])]
    jCstVar["jCstecfN2b2AK8"] = [fjets.ecfN2b2,[False],np.array([])]
    jCstVar["jCstecfN3b1AK8"] = [fjets.ecfN3b1,[False],np.array([])]
    jCstVar["jCstecfN3b2AK8"] = [fjets.ecfN3b2,[False],np.array([])]
    jCstVar["jCstEleEFractAK8"] = [fjets.electronEnergyFraction,[False],np.array([])]
    jCstVar["jCstEleMultAK8"] = [fjets.electronMultiplicity,[False],np.array([])]
    jCstVar["jCstGirthAK8"] = [fjets.girth,[False],np.array([])]
    jCstVar["jCstHfEMEFractAK8"] = [fjets.hfEMEnergyFraction,[False],np.array([])]
    jCstVar["jCstHfHadEFractAK8"] = [fjets.hfHadronEnergyFraction,[False],np.array([])]
    jCstVar["jCstMultAK8"] = [fjets.multiplicity,[False],np.array([])]
    jCstVar["jCstMuEFractAK8"] = [fjets.muonEnergyFraction,[False],np.array([])]
    jCstVar["jCstMuMultAK8"] = [fjets.muonMultiplicity,[False],np.array([])]
    jCstVar["jCstNeuEmEFractAK8"] = [fjets.neutralEmEnergyFraction,[False],np.array([])]
    jCstVar["jCstNeuHadEFractAK8"] = [fjets.neutralHadronEnergyFraction,[False],np.array([])]
    jCstVar["jCstNeuHadMultAK8"] = [fjets.neutralHadronMultiplicity,[False],np.array([])]
    jCstVar["jCstNeuMultAK8"] = [fjets.neutralMultiplicity,[False],np.array([])]
    jCstVar["jCstTau1AK8"] = [fjets.NsubjettinessTau1,[False],np.array([])]
    jCstVar["jCstTau2AK8"] = [fjets.NsubjettinessTau2,[False],np.array([])]
    jCstVar["jCstTau3AK8"] = [fjets.NsubjettinessTau3,[False],np.array([])]
    jCstVar["jCstNumBhadronsAK8"] = [fjets.NumBhadrons,[False],np.array([])]
    jCstVar["jCstNumChadronsAK8"] = [fjets.NumChadrons,[False],np.array([])]
    jCstVar["jCstPhoEFractAK8"] = [fjets.photonEnergyFraction,[False],np.array([])]
    jCstVar["jCstPhoMultAK8"] = [fjets.photonMultiplicity,[False],np.array([])]
    jCstVar["jCstPtDAK8"] = [fjets.ptD,[False],np.array([])]
    jCstVar["jCstSoftDropMassAK8"] = [fjets.softDropMass,[False],np.array([])]
    jCstVar["jCsthvCategory"] = [jetCats,[False],np.array([])]
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

    ht = ak.sum(jets.pt,axis=1)
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

    varVal['fjets'] = fjets
    varVal['evtw'] = evtw
    varVal['nl'] = (ak.num(electrons) + ak.num(muons))
    varVal['nnim'] = ak.num(nonIsoMuons)
    varVal['njets'] = ak.num(jets)
    varVal['njetsAK8'] = ak.num(fjets)
    varVal['nb'] = nBJets
    varVal['met'] = met
    varVal['ht'] = ht
    varVal['dPhiMinjMET'] = dPhiMinj
    varVal['dPhiMinjMETAK8'] = dPhiMinjAK8

    return varVal,jCst4vec,jCstVar
