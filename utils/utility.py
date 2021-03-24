import numpy as np
import awkward1 as ak
import awkward

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
