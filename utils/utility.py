import numpy as np
import awkward1 as ak
import awkward
from mt2 import mt2

def awkwardReshape(akArray,npArray):
    if len(akArray) == 0:
        return ak.Array([])
    else:
        return ak.broadcast_arrays(akArray.pt,1.0)[1] * npArray

# the two functions below return infinity when the event doesn't have the required
# number of jets or return the correct value for the jet variable
def jetVar(var,i):
    val = np.Inf
    if len(var) > i:
        val = var[i]
    return val

def jetVar_vec(var,i):
    vfunc = np.vectorize(lambda var,i: jetVar(var,i))
    return vfunc(var,i)

# convert phi values into spherical coordinate?
def phi_(x,y):
    phi = np.arctan2(y,x)
    if type(x) == awkward.array.jagged.JaggedArray:
        for i in range(len(phi)):
            for j in range(len(phi[i])):
                if phi[i][j] < 0:
                    phi[i][j] += 2*np.pi
        return phi
    elif type(x) == np.ndarray:
        return np.where(phi < 0, phi + 2*np.pi, phi)
    else:
        if phi < 0:
            phi += 2*np.pi
        return phi

def deltaPhi(jetphiL,metPhiL):
    phi1 = phi_( np.cos(jetphiL), np.sin(jetphiL) )
    phi2 = phi_( np.cos(metPhiL), np.sin(metPhiL) )
    dphi = phi1 - phi2

    if type(dphi) == awkward.array.jagged.JaggedArray:
        # this loops through the events
        for i in range(len(dphi)):
            # this loops through the jets in each event
            for j in range(len(dphi[i])):
                if dphi[i][j] < -np.pi:
                    dphi[i][j] += 2*np.pi
                elif dphi[i][j] > np.pi:
                    dphi[i][j] -= 2*np.pi
        return abs(dphi)
    elif type(dphi) == np.ndarray:
        dphi_edited = np.where(dphi < -np.pi, dphi + 2*np.pi, dphi)
        dphi_edited = np.where(dphi_edited > np.pi, dphi_edited - 2*np.pi, dphi_edited)
        return abs(dphi_edited)

def deltaPhiji(jetphi,metPhi):
    val = np.Inf
    if np.isfinite(jetphi) and np.isfinite(metPhi):
        phi1 = phi_( np.cos(jetphi), np.sin(jetphi) )
        phi2 = phi_( np.cos(metPhi), np.sin(metPhi) )
        dphi = phi1 - phi2
        if dphi < -np.pi:
            dphi += 2*np.pi
        elif dphi > np.pi:
            dphi -= 2*np.pi
        val = dphi
    return abs(val)

def deltaPhiji_vec(jetphi,metPhi):
    vfunc = np.vectorize(lambda jetphi,metPhi: deltaPhiji(jetphi,metPhi))
    return vfunc(jetphi,metPhi)

def deltaEta(eta0,eta1):
    val = np.Inf
    if np.isfinite(eta0) and np.isfinite(eta1):
        val = abs(eta0 - eta1)
    return val

def deltaEta_vec(eta0,eta1):
    vfunc = np.vectorize(lambda eta0,eta1: deltaEta(eta0,eta1))
    return vfunc(eta0,eta1)

def delta_R(eta0,eta1,phi0,phi1):
    dp = deltaPhiji_vec(phi0,phi1)
    deta = deltaEta_vec(eta0,eta1)
    deltaR2 = deta * deta + dp * dp
    return np.sqrt(deltaR2)

# the two functions below make sure that "divide" won't give us warning messages about invalid value
def divide(a,b):
    val = np.Inf
    if np.isfinite(a) and np.isfinite(b) and b > 0:
        val = a/b
    return val

def divide_vec(a,b):
    if len(a) > 0 and len(b) > 0:
        vfunc = np.vectorize(lambda a,b: divide(a,b))
        return vfunc(a,b)
    else:
        return np.array([])

def tauRatio(tau_a,tau_b,i):
    Ji_tau_a = jetVar_vec(tau_a,i)
    Ji_tau_b = jetVar_vec(tau_b,i)
    Ji_tau_ab = divide_vec(Ji_tau_a,Ji_tau_b)
    return Ji_tau_ab

def lorentzVector(pt,eta,phi,mass,i):
    px = pt[i]*np.cos(phi[i])
    py = pt[i]*np.sin(phi[i])
    pz = pt[i]*np.sinh(eta[i])
    p2 = px**2 + py**2 + pz**2
    m2 = mass[i]**2
    energy = np.sqrt(m2+p2)
    return np.array([energy,px,py,pz])

def mass4vec(m4vec):
    E = m4vec[0]
    px = m4vec[1]
    py = m4vec[2]
    pz = m4vec[3]
    return np.sqrt(E**2 - (px**2 + py**2 + pz**2))

def eta4vec(m4vec):
    px = m4vec[1]
    py = m4vec[2]
    pz = m4vec[3]
    pt = np.sqrt(px**2 + py**2)
    return np.arcsinh(pz/pt)

def phi4vec(m4vec):
    px = m4vec[1]
    py = m4vec[2]
    pt = np.sqrt(px**2 + py**2)
    return np.arcsin(py/pt)

def M_2J(j1,j2):
    totJets = j1+j2
    return mass4vec(totJets)

def MT2Cal(FDjet0,FSMjet0,FDjet1,FSMjet1,met,metPhi):
    Fjet0 = FDjet0 + FSMjet0
    Fjet1 = FDjet1 + FSMjet1
    METx = met*np.cos(metPhi)
    METy = met*np.sin(metPhi)
    MT2v = mt2(
    mass4vec(Fjet0), Fjet0[1], Fjet0[2],
    mass4vec(Fjet1), Fjet1[1], Fjet1[2],
    METx, METy, 0.0, 0.0, 0
    )
    return MT2v

def f4msmCom(pt,eta,phi,mass,met,metPhi,cut):
    MT2 = 0
    if len(pt) >= 4:
        List4jets_3Com = [[0,1,2,3],[0,2,1,3],[0,3,1,2]]
        diffList = []
        jetList = []
        for c in List4jets_3Com:
            jc1 = lorentzVector(pt,eta,phi,mass,c[0])
            jc2 = lorentzVector(pt,eta,phi,mass,c[1])
            jc3 = lorentzVector(pt,eta,phi,mass,c[2])
            jc4 = lorentzVector(pt,eta,phi,mass,c[3])
            jetList.append([jc1,jc2,jc3,jc4])
            diffList.append(abs(M_2J(jc1,jc2) - M_2J(jc3,jc4)))
        comIndex = np.argmin(diffList)
        msmJet = jetList[comIndex]
        # applying angular cut
        dEtaCut = 1.8
        dPhiCut = 0.9
        dRCutLow = 1.5
        dRCutHigh = 4.0
        dPhiMETCut = 1.5

        if cut == "dEta":
            if deltaEta(eta4vec(msmJet[0]),eta4vec(msmJet[1])) < dEtaCut and deltaEta(eta4vec(msmJet[2]),eta4vec(msmJet[3])) < dEtaCut:
                MT2 = MT2Cal(msmJet[0],msmJet[1],msmJet[2],msmJet[3],met,metPhi)
        elif cut == "dPhi":
            if deltaPhiji(phi4vec(msmJet[0]),phi4vec(msmJet[1])) < dPhiCut and deltaPhiji(phi4vec(msmJet[2]),phi4vec(msmJet[3])) > dPhiCut:
                MT2 = MT2Cal(msmJet[0],msmJet[1],msmJet[2],msmJet[3],met,metPhi)
        elif cut == "dR":
            if (dRCutLow < delta_R(eta4vec(msmJet[0]),eta4vec(msmJet[1]),phi4vec(msmJet[0]),phi4vec(msmJet[1])) < dRCutHigh) and (dRCutLow < delta_R(eta4vec(msmJet[2]),eta4vec(msmJet[3]),phi4vec(msmJet[2]),phi4vec(msmJet[3])) < dRCutHigh):
                MT2 = MT2Cal(msmJet[0],msmJet[1],msmJet[2],msmJet[3],met,metPhi)
        elif cut == "":
            MT2 = MT2Cal(msmJet[0],msmJet[1],msmJet[2],msmJet[3],met,metPhi)
    return MT2

def f4msmCom_vec(pt,eta,phi,mass,met,metPhi,cut):
    vfunc = np.vectorize(lambda pt,eta,phi,mass,met,metPhi,cut: f4msmCom(pt,eta,phi,mass,met,metPhi,cut))
    return vfunc(pt,eta,phi,mass,met,metPhi,cut)

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

def varGetter(inpObj):
    varVal = {}
    jets = inpObj["jets"]
    bjets = inpObj["bjets"]
    fjets = inpObj["fjets"]
    gfjets = inpObj["gfjets"]
    electrons = inpObj["electrons"]
    muons = inpObj["muons"]
    met = inpObj["met"]
    metPhi = inpObj["metPhi"]
    mtAK8 = inpObj["mtAK8"]
    madHT = inpObj["madHT"]
    # GenMT2_AK8 = inpObj["GenMT2_AK8"]
    # GenJetsAK8_darkPtFrac = inpObj["GenJetsAK8_darkPtFrac"]
    # GenJetsAK8_hvCategory = inpObj["GenJetsAK8_hvCategory"]
    # JetsAK8_hvCategory = inpObj["JetsAK8_hvCategory"]
    eCounter = inpObj["eCounter"]
    evtw = inpObj["evtw"]
    ew = awkwardReshape(electrons,evtw)
    mw = awkwardReshape(muons,evtw)
    jw = awkwardReshape(jets,evtw)
    fjw = awkwardReshape(fjets,evtw)
    gfjw = awkwardReshape(gfjets,evtw)
    ht = ak.sum(jets.pt,axis=1)
    st = ht + met
    metrht = divide_vec(met,ht)
    metrst = divide_vec(met,st)

    # AK4 Jet Variables
    jetPhi = jets.phi
    jetEta = jets.eta
    j1_eta = jetVar_vec(jetEta,0)
    j2_eta = jetVar_vec(jetEta,1)
    j1_phi = jetVar_vec(jetPhi,0)
    j2_phi = jetVar_vec(jetPhi,1)
    dPhij1 = deltaPhiji_vec(j1_phi,metPhi)
    dPhij2 = deltaPhiji_vec(j2_phi,metPhi)
    dPhij1rdPhij2 = divide_vec(dPhij1,dPhij2)
    dPhiMinj = deltaPhi(jetPhi,metPhi).min()
    dEtaj12 = deltaEta_vec(j1_eta,j2_eta)
    deltaR12j = delta_R(j1_eta,j2_eta,j1_phi,j2_phi)
    jweight = ak.flatten(jw)

    # AK8 Jet Variables
    jetAK8pT = fjets.pt
    jetAK8Phi = fjets.phi
    jetAK8Eta = fjets.eta
    jetAK8M = fjets.mass
    j1_etaAK8 = jetVar_vec(jetAK8Eta,0)
    j2_etaAK8 = jetVar_vec(jetAK8Eta,1)
    j1_phiAK8 = jetVar_vec(jetAK8Phi,0)
    j2_phiAK8 = jetVar_vec(jetAK8Phi,1)
    dPhij1AK8 = deltaPhiji_vec(j1_phiAK8,metPhi)
    dPhij2AK8 = deltaPhiji_vec(j2_phiAK8,metPhi)
    dPhij1rdPhij2AK8 = divide_vec(dPhij1AK8,dPhij2AK8)
    dPhijAK8 = deltaPhi(jetAK8Phi,metPhi)
    dPhiMinjAK8 = dPhijAK8.min()
    dEtaj12AK8 = deltaEta_vec(j1_etaAK8,j2_etaAK8)
    deltaR12jAK8 = delta_R(j1_etaAK8,j2_etaAK8,j1_phiAK8,j2_phiAK8)
    fjweight = ak.flatten(fjw)
    gfjweight = ak.flatten(gfjw)
    tau1 = fjets.tau1
    tau2 = fjets.tau2
    tau3 = fjets.tau3
    ecfN2b1 = fjets.ecfN2b1
    ecfN2b2 = fjets.ecfN2b2
    ecfN3b1 = fjets.ecfN3b1
    ecfN3b2 = fjets.ecfN3b2
    fEle = fjets.fEle
    fMu = fjets.fMu
    fNeuHad = fjets.fNeuHad
    fPho = fjets.fPho
    fNeuEM = fjets.fNeuEM
    fHFHad = fjets.fHFHad
    fHFEM = fjets.fHFEM
    fChEM = fjets.fChEM
    nPho = fjets.nPho
    nNeu = fjets.nNeu
    nNeuHad = fjets.nNeuHad
    nMu = fjets.nMu
    nEle = fjets.nEle
    nChHad = fjets.nChHad
    nCh = fjets.nCh
    mult = fjets.mult
    J_tau21 = divide_vec(tau2.flatten(),tau1.flatten())
    J_tau32 = divide_vec(tau3.flatten(),tau2.flatten())
    J1_tau21 = tauRatio(tau2,tau1,0)
    J1_tau32 = tauRatio(tau3,tau2,0)
    J2_tau21 = tauRatio(tau2,tau1,1)
    J2_tau32 = tauRatio(tau3,tau2,1)

    if len(bjets) > 0:
        nBJets = bjets.counts
    else:
        nBJets = np.zeros(len(evtw))

    varVal['eCounter'] = [eCounter,np.ones(len(eCounter))]
    varVal['evtw'] = [evtw,evtw]
    varVal['jw'] = [jweight,jweight]
    varVal['fjw'] = [fjweight,fjweight]
    varVal['njets'] = [jets.counts,evtw]
    varVal['njetsAK8'] = [fjets.counts,evtw]
    varVal['nb'] = [nBJets,evtw]
    varVal['nl'] = [(electrons.counts + muons.counts),evtw]
    varVal['ht'] = [ht,evtw]
    varVal['st'] = [st,evtw]
    varVal['met'] = [met,evtw]
    varVal['madHT'] = [madHT,evtw]
    varVal['jPt'] = [jets.pt.flatten(),jweight]
    varVal['jEta'] = [jetEta.flatten(),jweight]
    varVal['jPhi'] = [jetPhi.flatten(),jweight]
    varVal['jAxismajor'] = [jets.axismajor.flatten(),jweight]
    varVal['jAxisminor'] = [jets.axisminor.flatten(),jweight]
    varVal['jPtD'] = [jets.ptD.flatten(),jweight]
    varVal['dPhiMinjMET'] = [dPhiMinj,evtw]
    varVal['dEtaj12'] = [dEtaj12,evtw]
    varVal['dRJ12'] = [deltaR12j,evtw]
    varVal['jPtAK8'] = [fjets.pt.flatten(),fjweight]
    varVal['jEtaAK8'] = [jetAK8Eta.flatten(),fjweight]
    varVal['jPhiAK8'] = [jetAK8Phi.flatten(),fjweight]
    varVal['jAxismajorAK8'] = [fjets.axismajor.flatten(),fjweight]
    varVal['jAxisminorAK8'] = [fjets.axisminor.flatten(),fjweight]
    varVal['jGirthAK8'] = [fjets.girth.flatten(),fjweight]
    varVal['jPtDAK8'] = [fjets.ptD.flatten(),fjweight]
    varVal['jTau1AK8'] = [tau1.flatten(),fjweight]
    varVal['jTau2AK8'] = [tau2.flatten(),fjweight]
    varVal['jTau3AK8'] = [tau3.flatten(),fjweight]
    varVal['jTau21AK8'] = [J_tau21,fjweight]
    varVal['jTau32AK8'] = [J_tau32,fjweight]
    varVal['jSoftDropMassAK8'] = [fjets.softDropMass.flatten(),fjweight]
    varVal['jecfN2b1AK8'] = [ecfN2b1.flatten(),fjweight]
    varVal['jecfN2b2AK8'] = [ecfN2b2.flatten(),fjweight]
    varVal['jecfN3b1AK8'] = [ecfN3b1.flatten(),fjweight]
    varVal['jecfN3b2AK8'] = [ecfN3b2.flatten(),fjweight]
    varVal['jEleEFractAK8'] = [fEle.flatten(),fjweight]
    varVal['jMuEFractAK8'] = [fMu.flatten(),fjweight]
    varVal['jNeuHadEFractAK8'] = [fNeuHad.flatten(),fjweight]
    varVal['jPhoEFractAK8'] = [fPho.flatten(),fjweight]
    varVal['jNeuEmEFractAK8'] = [fNeuEM.flatten(),fjweight]
    varVal['jHfHadEFractAK8'] = [fHFHad.flatten(),fjweight]
    varVal['jHfEMEFractAK8'] = [fHFEM.flatten(),fjweight]
    varVal['jChEMEFractAK8'] = [fChEM.flatten(),fjweight]
    varVal['jPhoMultAK8'] = [nPho.flatten(),fjweight]
    varVal['jNeuMultAK8'] = [nNeu.flatten(),fjweight]
    varVal['jNeuHadMultAK8'] = [nNeuHad.flatten(),fjweight]
    varVal['jMuMultAK8'] = [nMu.flatten(),fjweight]
    varVal['jEleMultAK8'] = [nEle.flatten(),fjweight]
    varVal['jChHadMultAK8'] = [nChHad.flatten(),fjweight]
    varVal['jChMultAK8'] = [nCh.flatten(),fjweight]
    varVal['jMultAK8'] = [mult.flatten(),fjweight]
    varVal['dPhijMETAK8'] = [dPhijAK8.flatten(),fjweight]
    varVal['dPhiMinjMETAK8'] = [dPhiMinjAK8,evtw]
    varVal['dEtaj12AK8'] = [dEtaj12AK8,evtw]
    varVal['dRJ12AK8'] = [deltaR12jAK8,evtw]
    varVal['mT'] = [mtAK8,evtw]
    varVal['METrHT_pt30'] = [metrht,evtw]
    varVal['METrST_pt30'] = [metrst,evtw]
    varVal['j1Pt'] = [jetVar_vec(jets.pt,0),evtw]
    varVal['j1Eta'] = [j1_eta,evtw]
    varVal['j1Phi'] = [j1_phi,evtw]
    varVal['j1Axismajor'] = [jetVar_vec(jets.axismajor,0),evtw]
    varVal['j1Axisminor'] = [jetVar_vec(jets.axisminor,0),evtw]
    varVal['j1PtD'] = [jetVar_vec(jets.ptD,0),evtw]
    varVal['dPhij1MET'] = [dPhij1,evtw]
    varVal['j2Pt'] = [jetVar_vec(jets.pt,1),evtw]
    varVal['j2Eta'] = [j2_eta,evtw]
    varVal['j2Phi'] = [j2_phi,evtw]
    varVal['j2Axismajor'] = [jetVar_vec(jets.axismajor,1),evtw]
    varVal['j2Axisminor'] = [jetVar_vec(jets.axisminor,1),evtw]
    varVal['j2PtD'] = [jetVar_vec(jets.ptD,1),evtw]
    varVal['dPhij2MET'] = [dPhij2,evtw]
    varVal['dPhij1rdPhij2'] = [dPhij1rdPhij2,evtw]
    varVal['j1PtAK8'] = [jetVar_vec(fjets.pt,0),evtw]
    varVal['j1EtaAK8'] = [j1_etaAK8,evtw]
    varVal['j1PhiAK8'] = [j1_phiAK8,evtw]
    varVal['j1AxismajorAK8'] = [jetVar_vec(fjets.axismajor,0),evtw]
    varVal['j1AxisminorAK8'] = [jetVar_vec(fjets.axisminor,0),evtw]
    varVal['j1GirthAK8'] = [jetVar_vec(fjets.girth,0),evtw]
    varVal['j1PtDAK8'] = [jetVar_vec(fjets.ptD,0),evtw]
    varVal['j1Tau1AK8'] = [jetVar_vec(fjets.tau1,0),evtw]
    varVal['j1Tau2AK8'] = [jetVar_vec(fjets.tau2,0),evtw]
    varVal['j1Tau3AK8'] = [jetVar_vec(fjets.tau3,0),evtw]
    varVal['j1Tau21AK8'] = [J1_tau21,evtw]
    varVal['j1Tau32AK8'] = [J1_tau32,evtw]
    varVal['j1SoftDropMassAK8'] = [jetVar_vec(fjets.softDropMass,0),evtw]
    varVal['dPhij1METAK8'] = [dPhij1AK8,evtw]
    varVal['j2PtAK8'] = [jetVar_vec(fjets.pt,1),evtw]
    varVal['j2EtaAK8'] = [j2_etaAK8,evtw]
    varVal['j2PhiAK8'] = [j2_phiAK8,evtw]
    varVal['j2AxismajorAK8'] = [jetVar_vec(fjets.axismajor,1),evtw]
    varVal['j2AxisminorAK8'] = [jetVar_vec(fjets.axisminor,1),evtw]
    varVal['j2GirthAK8'] = [jetVar_vec(fjets.girth,1),evtw]
    varVal['j2PtDAK8'] = [jetVar_vec(fjets.ptD,1),evtw]
    varVal['j2Tau1AK8'] = [jetVar_vec(fjets.tau1,1),evtw]
    varVal['j2Tau2AK8'] = [jetVar_vec(fjets.tau2,1),evtw]
    varVal['j2Tau3AK8'] = [jetVar_vec(fjets.tau3,1),evtw]
    varVal['j2Tau21AK8'] = [J2_tau21,evtw]
    varVal['j2Tau32AK8'] = [J2_tau32,evtw]
    varVal['j2SoftDropMassAK8'] = [jetVar_vec(fjets.softDropMass,1),evtw]
    varVal['dPhij2METAK8'] = [dPhij2AK8,evtw]
    varVal['dPhij1rdPhij2AK8'] = [dPhij1rdPhij2AK8,evtw]
    # varVal['mT2_f4_msm'] = [f4msmCom_vec(jetAK8pT,jetAK8Eta,jetAK8Phi,jetAK8M,met,metPhi,""),evtw]
    # varVal['mT2_f4_msm_dEta'] = [f4msmCom_vec(jetAK8pT,jetAK8Eta,jetAK8Phi,jetAK8M,met,metPhi,"dEta"),evtw]
    # varVal['mT2_f4_msm_dPhi'] = [f4msmCom_vec(jetAK8pT,jetAK8Eta,jetAK8Phi,jetAK8M,met,metPhi,"dPhi"),evtw]
    # varVal['mT2_f4_msm_dR'] = [f4msmCom_vec(jetAK8pT,jetAK8Eta,jetAK8Phi,jetAK8M,met,metPhi,"dR"),evtw]
    # varVal['GenMT2_AK8'] = [GenMT2_AK8,evtw]
    # varVal['GenJetsAK8_hvCategory'] = [GenJetsAK8_hvCategory.flatten(),gfjweight]
    # varVal['GenJetsAK8_darkPtFrac'] = [GenJetsAK8_darkPtFrac.flatten(),gfjweight]
    # varVal['JetsAK8_hvCategory'] = [ak.flatten(JetsAK8_hvCategory),fjweight]
    return varVal
