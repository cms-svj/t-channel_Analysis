from coffea.analysis_objects import JaggedCandidateArray
import numpy as np
import awkward1 as ak
from . import utility as utl

class Objects:
    def __init__(self, df):
        #self.df = df
        self.electrons = JaggedCandidateArray.candidatesfromcounts(
            df['Electrons'].counts,
            px=df['Electrons'].fP.fX.flatten(),
            py=df['Electrons'].fP.fY.flatten(),
            pz=df['Electrons'].fP.fZ.flatten(),
            energy=df['Electrons'].fE.flatten(),
            passIso=df['Electrons_passIso'].flatten(),
            Electrons_MiniIso = df['Electrons_MiniIso'].flatten(),
            charge=df['Electrons_charge'].flatten(),
            iD=df['Electrons_mediumID'].flatten()
        )
        self.muons = JaggedCandidateArray.candidatesfromcounts(
            df['Muons'].counts,
            px=df['Muons'].fP.fX.flatten(),
            py=df['Muons'].fP.fY.flatten(),
            pz=df['Muons'].fP.fZ.flatten(),
            energy=df['Muons'].fE.flatten(),
            passIso=df['Muons_passIso'].flatten(),
            Muons_MiniIso = df['Muons_MiniIso'].flatten(),
            charge=df['Muons_charge'].flatten(),
            iD=df['Muons_mediumID'].flatten()
        )
        self.jets = JaggedCandidateArray.candidatesfromcounts(
            df['Jets'].counts,
            px=df['Jets'].fP.fX.flatten(),
            py=df['Jets'].fP.fY.flatten(),
            pz=df['Jets'].fP.fZ.flatten(),
            energy=df['Jets'].fE.flatten(),
            axismajor=df['Jets_axismajor'].flatten(),
            axisminor=df['Jets_axisminor'].flatten(),
            ptD=df['Jets_ptD'].flatten(),
            bDeepCSVprobb=df['Jets_bJetTagDeepCSVprobb'].flatten(),
            bDeepCSVprobbb=df['Jets_bJetTagDeepCSVprobbb'].flatten()
        )
        self.fjets = JaggedCandidateArray.candidatesfromcounts(
            df['JetsAK8'].counts,
            px=df['JetsAK8'].fP.fX.flatten(),
            py=df['JetsAK8'].fP.fY.flatten(),
            pz=df['JetsAK8'].fP.fZ.flatten(),
            energy=df['JetsAK8'].fE.flatten(),
            axismajor=df['JetsAK8_axismajor'].flatten(),
            axisminor=df['JetsAK8_axisminor'].flatten(),
            girth=df['JetsAK8_girth'].flatten(),
            ptD=df['JetsAK8_ptD'].flatten(),
            tau1=df['JetsAK8_NsubjettinessTau1'].flatten(),
            tau2=df['JetsAK8_NsubjettinessTau2'].flatten(),
            tau3=df['JetsAK8_NsubjettinessTau3'].flatten(),
            ecfN2b1=df['JetsAK8_ecfN2b1'].flatten(),
            ecfN2b2=df['JetsAK8_ecfN2b2'].flatten(),
            ecfN3b1=df['JetsAK8_ecfN3b1'].flatten(),
            ecfN3b2=df['JetsAK8_ecfN3b2'].flatten(),
            fEle=df['JetsAK8_electronEnergyFraction'].flatten(),
            fMu=df['JetsAK8_muonEnergyFraction'].flatten(),
            fNeuHad=df['JetsAK8_neutralHadronEnergyFraction'].flatten(),
            fPho=df['JetsAK8_photonEnergyFraction'].flatten(),
            fNeuEM=df['JetsAK8_neutralEmEnergyFraction'].flatten(),
            fHFHad=df['JetsAK8_hfHadronEnergyFraction'].flatten(),
            fHFEM=df['JetsAK8_hfEMEnergyFraction'].flatten(),
            fChEM=df['JetsAK8_chargedEmEnergyFraction'].flatten(),
            nPho=df['JetsAK8_photonMultiplicity'].flatten(),
            nNeu=df['JetsAK8_neutralMultiplicity'].flatten(),
            nNeuHad=df['JetsAK8_neutralHadronMultiplicity'].flatten(),
            nMu=df['JetsAK8_muonMultiplicity'].flatten(),
            nEle=df['JetsAK8_electronMultiplicity'].flatten(),
            nChHad=df['JetsAK8_chargedHadronMultiplicity'].flatten(),
            nCh=df['JetsAK8_chargedMultiplicity'].flatten(),
            mult=df['JetsAK8_multiplicity'].flatten(),
            softDropMass=df['JetsAK8_softDropMass'].flatten(),
        )
        self.gfjets = JaggedCandidateArray.candidatesfromcounts(
            df['GenJetsAK8'].counts,
            px=df['GenJetsAK8'].fP.fX.flatten(),
            py=df['GenJetsAK8'].fP.fY.flatten(),
            pz=df['GenJetsAK8'].fP.fZ.flatten(),
            energy=df['GenJetsAK8'].fE.flatten(),
        )
        # Quality cut
        self.etaCut = 2.4
        self.leptonPt = 10

    def goodElectrons(self):
        # # Good Electrons
        electronQualityCut = (self.electrons.pt > self.leptonPt) & (abs(self.electrons.eta) < self.etaCut) & (self.electrons.Electrons_MiniIso < 0.1)
        return self.electrons[electronQualityCut]

    def goodMuons(self):
        # # Good Muons
        muonQualityCut = (self.muons.pt > self.leptonPt) & (abs(self.muons.eta) < self.etaCut) & (self.muons.Muons_MiniIso < 0.4)
        return self.muons[muonQualityCut]

    def goodJets(self):
        # # Good AK4 Jets Cut
        ak4QualityCut = (self.jets.pt > 30) & (abs(self.jets.eta) < 5.0)
        return self.jets[ak4QualityCut]

    def goodFatJets(self):
        # # Good AK8 Jets Cut
        ak8QualityCut = (self.fjets.pt > 170) & (abs(self.fjets.eta) < 5.0)
        return self.fjets[ak8QualityCut]

    def goodFatJetCut(self):
        # # Good AK8 Jets Cut
        ak8QualityCut = (self.fjets.pt > 170) & (abs(self.fjets.eta) < 5.0)
        return ak8QualityCut

    def goodGenFatJets(self):
        # # Good AK8 Jets Cut
        ak8QualityCut = self.gfjets.pt > 170 & (abs(self.gfjets.eta) < 5.0)
        return self.gfjets[ak8QualityCut]

    def goodBJets(self,df,jets):
        if len(jets.flatten()) > 0:
            jets = jets[~np.isnan(jets.bDeepCSVprobb)]
            deepCSV = jets.bDeepCSVprobb + jets.bDeepCSVprobbb
            dataset = df['dataset']
            # use medium working point
            if "2016" in dataset:
                return jets[deepCSV >= 0.6321]
            elif "2017" in dataset:
                return jets[deepCSV >= 0.4941]
            elif "2018" in dataset:
                return jets[deepCSV >= 0.4184]
        else:
            return ak.Array([])

def inpObj(df,scaleFactor):
    inpObj = {}
    luminosity = 21071.0+38654.0
    inpObj["evtw"] = luminosity*df['Weight']*scaleFactor
    inpObj["eCounter"] = np.where(inpObj["evtw"] >= 0, 1, -1)

    obj = Objects(df)
    inpObj["electrons"] = obj.goodElectrons()
    inpObj["muons"] = obj.goodMuons()
    inpObj["jets"] = obj.goodJets()
    inpObj["fjets"] = obj.goodFatJets()
    inpObj["gfjets"] = obj.goodGenFatJets()
    inpObj["bjets"] = obj.goodBJets(df,inpObj["jets"])
    inpObj["metPhi"] = df['METPhi']
    inpObj["madHT"] = df['madHT']
    inpObj["genMET"] = df['GenMET']
    inpObj["met"] = df['MET']
    inpObj["mtAK8"] = df['MT_AK8']
    inpObj["triggerPass"] = df['TriggerPass']
    inpObj["jetID"] = df['JetID']
    inpObj["jetIDAK8"] = df['JetIDAK8']
    inpObj["ht"] = ak.sum(inpObj["jets"].pt,axis=1)
    inpObj["dPhiMinj"] = utl.deltaPhi(inpObj["jets"].phi,inpObj["metPhi"]).min()
    inpObj["dPhiMinjAK8"] = utl.deltaPhi(inpObj["fjets"].phi,inpObj["metPhi"]).min()
    jetAK8Eta = inpObj["fjets"].eta
    jetAK8Phi = inpObj["fjets"].phi
    j1_etaAK8 = utl.jetVar_vec(jetAK8Eta,0)
    j2_etaAK8 = utl.jetVar_vec(jetAK8Eta,1)
    j1_phiAK8 = utl.jetVar_vec(jetAK8Phi,0)
    j2_phiAK8 = utl.jetVar_vec(jetAK8Phi,1)
    inpObj["deltaR12jAK8"] = utl.delta_R(j1_etaAK8,j2_etaAK8,j1_phiAK8,j2_phiAK8)
    ## matching AK8 Jet to AK8 GenJet to get jet category for AK8 Jet
    ## still need to fix this, not working fully yet
    ## GenJetsAK8_hvCategory is only present in the signal samples, not the V17 background
    # jetCats = []
    # jetsAK8GenInd = df["JetsAK8_genIndex"][obj.goodFatJetCut()]
    # for gji in range(len(jetsAK8GenInd)):
    #     genInd = jetsAK8GenInd[gji]
    #     genCat = df["GenJetsAK8_hvCategory"][gji]
    #     if len(genCat) > 0:
    #         jetCats.append(list(genCat[genInd]))
    # jetCats = ak.Array(jetCats)
    # inpObj["JetsAK8_hvCategory"] = jetCats
    # inpObj["GenJetsAK8_hvCategory"] = df['GenJetsAK8_hvCategory']
    # inpObj["GenMT2_AK8"] = df['GenMT2_AK8']
    # inpObj["GenJetsAK8_darkPtFrac"] = df['GenJetsAK8_darkPtFrac']

    return inpObj
