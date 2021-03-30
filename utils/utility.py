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
    vec4 = ak.zip(
        {
            "pt": [pt[i]],
            "eta": [eta[i]],
            "phi": [phi[i]],
            "mass": [mass[i]]
        },
        with_name="PtEtaPhiMLorentzVector",
    )
    return vec4

def M_2J(j1,j2):
    totJets = j1+j2
    return totJets.mass

def MT2Cal(FDjet0,FSMjet0,FDjet1,FSMjet1,met,metPhi):
    Fjet0 = FDjet0 + FSMjet0
    Fjet1 = FDjet1 + FSMjet1
    METx = met*np.cos(metPhi)
    METy = met*np.sin(metPhi)
    MT2v = mt2(
    Fjet0.mass[0], Fjet0.pt[0] * np.cos(Fjet0.phi[0]), Fjet0.pt[0] * np.sin(Fjet0.phi[0]),
    Fjet1.mass[0], Fjet1.pt[0] * np.cos(Fjet1.phi[0]), Fjet1.pt[0] * np.sin(Fjet1.phi[0]),
    METx, METy, 0.0, 0.0, 0
    )
    return MT2v

def f4msmCom(pt,eta,phi,mass,met,metPhi,cut):
    MT2 = np.Inf
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
            if deltaEta(msmJet[0].eta[0],msmJet[1].eta[0]) < dEtaCut and deltaEta(msmJet[2].eta[0],msmJet[3].eta[0]) < dEtaCut:
                MT2 = MT2Cal(msmJet[0],msmJet[1],msmJet[2],msmJet[3],met,metPhi)
        elif cut == "dPhi":
            if deltaPhiji(msmJet[0].phi[0],msmJet[1].phi[0]) > dPhiCut and deltaPhiji(msmJet[2].phi[0],msmJet[3].phi[0]) > dPhiCut:
                MT2 = MT2Cal(msmJet[0],msmJet[1],msmJet[2],msmJet[3],met,metPhi)
        elif cut == "dR":
            if (dRCutLow < delta_R(msmJet[0].eta[0],msmJet[1].eta[0],msmJet[0].phi[0],msmJet[1].phi[0]) < dRCutHigh) and (dRCutLow < delta_R(msmJet[2].eta[0],msmJet[3].eta[0],msmJet[2].phi[0],msmJet[3].phi[0]) < dRCutHigh):
                MT2 = MT2Cal(msmJet[0],msmJet[1],msmJet[2],msmJet[3],met,metPhi)
        elif cut == "":
            MT2 = MT2Cal(msmJet[0],msmJet[1],msmJet[2],msmJet[3],met,metPhi)
    return MT2

def f4msmCom_vec(pt,eta,phi,mass,met,metPhi,cut):
    vfunc = np.vectorize(lambda pt,eta,phi,mass,met,metPhi,cut: f4msmCom(pt,eta,phi,mass,met,metPhi,cut))
    return vfunc(pt,eta,phi,mass,met,metPhi,cut)
