import numpy as np 
import awkward as ak

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
    return dphi_edited

def deltaEta(eta0,eta1):
  return (eta0 - eta1)

def delta_R(eta0,eta1,phi0,phi1):
    dp = deltaPhi(phi0,phi1)
    deta = deltaEta(eta0,eta1)
    deltaR2 = deta * deta + dp * dp
    return np.sqrt(deltaR2)

def run_jet_constituent_matching(events, orig_jets):
    jets = orig_jets[:]
    jetConstituents = events.JetsConstituents[:]
    jetsAK8_constituentsIndex = jets.constituentsIndex
    indices = ak.flatten(jetsAK8_constituentsIndex,axis=-1)
    jetConstituentsForJets = jetConstituents[indices]
    jetConstJetArray = ak.unflatten(jetConstituentsForJets,ak.flatten(ak.num(jetsAK8_constituentsIndex,axis=-1)),axis=1)
    jets["Constituents"] = jetConstJetArray
    return jets
