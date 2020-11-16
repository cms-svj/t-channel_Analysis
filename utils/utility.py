import numpy as np
import awkward1 as ak
import awkward

def awkwardReshape(akArray,npArray):
    return ak.broadcast_arrays(akArray.pt,1.0)[1] * npArray

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

def deltaPhi(jetphiL,metPhiL):
    phi1 = phi_( np.cos(jetphiL), np.sin(jetphiL) )
    phi2 = phi_( np.cos(metPhiL), np.sin(metPhiL) )
    dphi = phi1 - phi2

    if type(dphi) == awkward.array.jagged.JaggedArray:
        for i in range(len(dphi)):
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

def delta_R(jet1_phi,jet2_phi,jet1_eta,jet2_eta):
    dp = deltaPhi(jet1_phi,jet2_phi)
    deltaR2 = (jet1_eta - jet2_eta) * (jet1_eta - jet2_eta) + dp * dp
    return np.sqrt(deltaR2)

