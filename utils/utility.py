import numpy as np
import awkward as ak
from mt2 import mt2
from . import objects as ob
from itertools import combinations

def awkwardReshape(akArray,npArray):
    if len(akArray) == 0:
        return ak.Array([])
    else:
        return ak.broadcast_arrays(akArray.pt,1.0)[1] * npArray

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
        if "2016" in dataset:
            luminosity = 35921.036
        elif "2017" in dataset:
            luminosity = 41521.331
        elif "2018" in dataset:
            luminosity = 59692.692
        evtw = luminosity*events.Weight*scaleFactor
    eCounter = np.where(evtw >= 0, 1, -1)
    obj = ob.Objects(events)
    jets = obj.goodJets()
    bjets = obj.goodBJets(dataset,jets)
    fjets = obj.goodFatJets()
    # gfjets = obj.goodGenFatJets()
    electrons = obj.goodElectrons()
    muons = obj.goodMuons()
    met = events.MET
    metPhi = events.METPhi
    mtAK8 = events.MT_AK8
    jetAK8Eta = fjets.eta
    jetAK8Phi = fjets.phi
    j1_etaAK8 = jetVar_i(jetAK8Eta,0)
    j2_etaAK8 = jetVar_i(jetAK8Eta,1)
    j1_phiAK8 = jetVar_i(jetAK8Phi,0)
    j2_phiAK8 = jetVar_i(jetAK8Phi,1)

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
        varVal['JetsAK8_hvCategory'] = [jetCats,'fjw']
    else:
        varVal['JetsAK8_hvCategory'] = [awkwardReshape(fjets,np.ones(len(evtw))*-1),'fjw']

    ew = awkwardReshape(electrons,evtw)
    mw = awkwardReshape(muons,evtw)
    jw = awkwardReshape(jets,evtw)
    fjw = awkwardReshape(fjets,evtw)
    ht = ak.sum(jets.pt,axis=1)
    st = ht + met
    metrht = met/ht
    metrst = met/st

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

    if len(bjets) > 0:
        nBJets = ak.num(bjets)
    else:
        nBJets = np.zeros(len(evtw))

    varVal['jets'] = jets
    varVal['fjets'] = fjets
    varVal['eCounter'] = [eCounter,'w1']
    varVal['evtw'] = [evtw,'evtw']
    varVal['jw'] = [jw,'jw']
    varVal['fjw'] = [fjw,'fjw']
    varVal['ew'] = [ew,'ew']
    varVal['mw'] = [mw,'mw']
    varVal['njets'] = [ak.num(jets),'evtw']
    varVal['njetsAK8'] = [ak.num(fjets),'evtw']
    varVal['nb'] = [nBJets,'evtw']
    varVal['nl'] = [(ak.num(electrons) + ak.num(muons)),'evtw']
    varVal['ht'] = [ht,'evtw']
    varVal['st'] = [st,'evtw']
    varVal['met'] = [met,'evtw']
    varVal['metPhi'] = [metPhi, 'evtw']
    # varVal['madHT'] = [madHT,'evtw']
    varVal['jPt'] = [jets.pt,'jw']
    varVal['jEta'] = [jetEta,'jw']
    varVal['jPhi'] = [jetPhi,'jw']
    varVal['jAxismajor'] = [jets.axismajor,'jw']
    varVal['jAxisminor'] = [jets.axisminor,'jw']
    varVal['jPtD'] = [jets.ptD,'jw']
    varVal['dPhiMinjMET'] = [dPhiMinj,'evtw']
    varVal['jPtAK8'] = [fjets.pt,'fjw']
    varVal['jEtaAK8'] = [jetAK8Eta,'fjw']
    varVal['jPhiAK8'] = [jetAK8Phi,'fjw']
    varVal['jAxismajorAK8'] = [fjets.axismajor,'fjw']
    varVal['jAxisminorAK8'] = [fjets.axisminor,'fjw']
    varVal['jChEMEFractAK8'] = [fjets.chargedEmEnergyFraction,'fjw']
    varVal['jChHadEFractAK8'] = [fjets.chargedHadronEnergyFraction,'fjw']
    varVal['jChHadMultAK8'] = [fjets.chargedHadronMultiplicity,'fjw']
    varVal['jChMultAK8'] = [fjets.chargedMultiplicity,'fjw']
    varVal['jdoubleBDiscriminatorAK8'] = [fjets.doubleBDiscriminator,'fjw']
    varVal['jecfN2b1AK8'] = [fjets.ecfN2b1,'fjw']
    varVal['jecfN2b2AK8'] = [fjets.ecfN2b2,'fjw']
    varVal['jecfN3b1AK8'] = [fjets.ecfN3b1,'fjw']
    varVal['jecfN3b2AK8'] = [fjets.ecfN3b2,'fjw']
    varVal['jEleEFractAK8'] = [fjets.electronEnergyFraction,'fjw']
    varVal['jEleMultAK8'] = [fjets.electronMultiplicity,'fjw']
    varVal['jGirthAK8'] = [fjets.girth,'fjw']
    varVal['jHfEMEFractAK8'] = [fjets.hfEMEnergyFraction,'fjw']
    varVal['jHfHadEFractAK8'] = [fjets.hfHadronEnergyFraction,'fjw']
    varVal['jMultAK8'] = [fjets.multiplicity,'fjw']
    varVal['jMuEFractAK8'] = [fjets.muonEnergyFraction,'fjw']
    varVal['jMuMultAK8'] = [fjets.muonMultiplicity,'fjw']
    varVal['jNeuEmEFractAK8'] = [fjets.neutralEmEnergyFraction,'fjw']
    varVal['jNeuHadEFractAK8'] = [fjets.neutralHadronEnergyFraction,'fjw']
    varVal['jNeuHadMultAK8'] = [fjets.neutralHadronMultiplicity,'fjw']
    varVal['jNeuMultAK8'] = [fjets.neutralMultiplicity,'fjw']
    varVal['jTau1AK8'] = [tau1,'fjw']
    varVal['jTau2AK8'] = [tau2,'fjw']
    varVal['jTau3AK8'] = [tau3,'fjw']
    varVal['jTau21AK8'] = [J_tau21,'fjw']
    varVal['jTau32AK8'] = [J_tau32,'fjw']
    varVal['jNumBhadronsAK8'] = [fjets.NumBhadrons,'fjw']
    varVal['jNumChadronsAK8'] = [fjets.NumChadrons,'fjw']
    varVal['jPhoEFractAK8'] = [fjets.photonEnergyFraction,'fjw']
    varVal['jPhoMultAK8'] = [fjets.photonMultiplicity,'fjw']
    varVal['jPtDAK8'] = [fjets.ptD,'fjw']
    varVal['jSoftDropMassAK8'] = [fjets.softDropMass,'fjw']
    varVal['dPhijMETAK8'] = [dPhijAK8,'fjw']
    varVal['dPhiMinjMETAK8'] = [dPhiMinjAK8,'evtw']
    varVal['mT'] = [mtAK8,'evtw']
    varVal['METrHT_pt30'] = [metrht,'evtw']
    varVal['METrST_pt30'] = [metrst,'evtw']
    varVal['electronsIso'] = [electrons.iso,'ew']
    varVal['muonsIso'] = [muons.iso,'mw']
    # preparing histograms for jN variables
    maxN = 4
    for i in range(maxN):
        varVal['j{}Pt'.format(i+1)] = [jetVar_i(jets.pt,i),'evtw']
        varVal['j{}Eta'.format(i+1)] = [jetVar_i(jetEta,i),'evtw']
        varVal['j{}Phi'.format(i+1)] = [jetVar_i(jetPhi,i),'evtw']
        varVal['j{}Axismajor'.format(i+1)] = [jetVar_i(jets.axismajor,i),'evtw']
        varVal['j{}Axisminor'.format(i+1)] = [jetVar_i(jets.axisminor,i),'evtw']
        varVal['j{}PtD'.format(i+1)] = [jetVar_i(jets.ptD,i),'evtw']
        varVal['dPhij{}MET'.format(i+1)] = [deltaPhi(jetVar_i(jetPhi,i),metPhi),'evtw']
        varVal['j{}PtAK8'.format(i+1)] = [jetVar_i(fjets.pt,i),'evtw']
        varVal['j{}EtaAK8'.format(i+1)] = [jetVar_i(jetAK8Eta,i),'evtw']
        varVal['j{}PhiAK8'.format(i+1)] = [jetVar_i(jetAK8Phi,i),'evtw']
        varVal['j{}AxismajorAK8'.format(i+1)] = [jetVar_i(fjets.axismajor,i),'evtw']
        varVal['j{}AxisminorAK8'.format(i+1)] = [jetVar_i(fjets.axisminor,i),'evtw']
        varVal['j{}GirthAK8'.format(i+1)] = [jetVar_i(fjets.girth,i),'evtw']
        varVal['j{}PtDAK8'.format(i+1)] = [jetVar_i(fjets.ptD,i),'evtw']
        varVal['j{}Tau1AK8'.format(i+1)] = [jetVar_i(tau1,i),'evtw']
        varVal['j{}Tau2AK8'.format(i+1)] = [jetVar_i(tau2,i),'evtw']
        varVal['j{}Tau3AK8'.format(i+1)] = [jetVar_i(tau3,i),'evtw']
        varVal['j{}Tau21AK8'.format(i+1)] = [tauRatio(tau2,tau1,i),'evtw']
        varVal['j{}Tau32AK8'.format(i+1)] = [tauRatio(tau3,tau2,i),'evtw']
        varVal['j{}SoftDropMassAK8'.format(i+1)] = [jetVar_i(fjets.softDropMass,i),'evtw']
        varVal['dPhij{}METAK8'.format(i+1)] = [deltaPhi(jetVar_i(jetAK8Phi,i),metPhi),'evtw']
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
        varVal['dEtaj{}{}'.format(j1+1,j2+1)] = [deltaEta(j1_eta,j2_eta),'evtw']
        varVal['dPhij{}{}'.format(j1+1,j2+1)] = [deltaPhi(j1_phi,j2_phi),'evtw']
        varVal['dRj{}{}'.format(j1+1,j2+1)] = [delta_R(j1_eta,j2_eta,j1_phi,j2_phi),'evtw']
        varVal['dPhij{}rdPhij{}'.format(j1+1,j2+1)] = [dPhij1/dPhij2,'evtw']
        varVal['dEtaj{}{}AK8'.format(j1+1,j2+1)] = [deltaEta(j1_etaAK8,j2_etaAK8),'evtw']
        varVal['dPhij{}{}AK8'.format(j1+1,j2+1)] = [deltaPhi(j1_phiAK8,j2_phiAK8),'evtw']
        varVal['dRj{}{}AK8'.format(j1+1,j2+1)] = [delta_R(j1_etaAK8,j2_etaAK8,j1_phiAK8,j2_phiAK8),'evtw']
        varVal['dPhij{}rdPhij{}AK8'.format(j1+1,j2+1)] = [dPhij1AK8/dPhij2AK8,'evtw']
    # varVal['GenJetsAK8_hvCategory'] = [GenJetsAK8_hvCategory.flatten(),gfjweight]
    # varVal['mT2_f4_msm'] = [f4msmCom_vec(jetAK8pT,jetAK8Eta,jetAK8Phi,jetAK8M,met,metPhi,""),'evtw']
    # varVal['mT2_f4_msm_dEta'] = [f4msmCom_vec(jetAK8pT,jetAK8Eta,jetAK8Phi,jetAK8M,met,metPhi,"dEta"),'evtw']
    # varVal['mT2_f4_msm_dPhi'] = [f4msmCom_vec(jetAK8pT,jetAK8Eta,jetAK8Phi,jetAK8M,met,metPhi,"dPhi"),'evtw']
    # varVal['mT2_f4_msm_dR'] = [f4msmCom_vec(jetAK8pT,jetAK8Eta,jetAK8Phi,jetAK8M,met,metPhi,"dR"),'evtw']
    # varVal['GenJetsAK8_darkPtFrac'] = [GenJetsAK8_darkPtFrac.flatten(),gfjweight]
    # varVal['GenMT2_AK8'] = [GenMT2_AK8,'evtw']
    return varVal
