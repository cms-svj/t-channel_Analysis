import numpy as np
import awkward as ak
from . import utility as utl

class Objects:
    def __init__(self, events):
        self.electrons = events.Electrons
        self.muons = events.Muons
        self.jets = events.Jets
        self.fjets = events.JetsAK8
        self.gfjets = events.GenJetsAK8
        # Quality cut
        self.etaCut = 2.4
        self.leptonPt = 10.0

    def goodElectrons(self):
        # # Good Electrons
        electronQualityCut = (self.electrons.pt > self.leptonPt) & (abs(self.electrons.eta) < self.etaCut) & (self.electrons.iso < 0.1)
        return self.electrons[electronQualityCut]

    def goodMuons(self):
        # # Good Muons
        muonQualityCut = (self.muons.pt > self.leptonPt) & (abs(self.muons.eta) < self.etaCut) & (self.muons.iso < 0.4)
        return self.muons[muonQualityCut]

    def nonIsoMuons(self):
        # # Good Muons
        nonIsomuonQualityCut = (self.muons.pt > 55.0) & (abs(self.muons.eta) < self.etaCut) & (self.muons.passIso == False) & (self.muons.mediumID == True)
        return self.muons[nonIsomuonQualityCut]

    def goodJets(self):
        # # Good AK4 Jets Cut
        ak4QualityCut = (self.jets.pt > 30) & (abs(self.jets.eta) < 2.4) & (self.jets.ID == True)
        return self.jets[ak4QualityCut]

    def goodFatJets(self):
        # # Good AK8 Jets Cut
        ak8QualityCut = (self.fjets.pt > 170) & (abs(self.fjets.eta) < 5.0) & (self.fjets.ID == True)
        return self.fjets[ak8QualityCut]

    def goodGenFatJets(self):
        # # Good AK8 Jets Cut
        ak8QualityCut = self.gfjets.pt > 170 & (abs(self.gfjets.eta) < 5.0)
        return self.gfjets[ak8QualityCut]

    def goodBJets(self,dataset,jets):
        if len(jets) > 0:
            jets = jets[~np.isnan(jets.bJetTagDeepCSVprobb)]
            deepCSV = jets.bJetTagDeepCSVprobb + jets.bJetTagDeepCSVprobbb
            # use medium working point
            if "2016" in dataset:
                return jets[deepCSV >= 0.6321]
            elif "2017" in dataset:
                return jets[deepCSV >= 0.4941]
            elif "2018" in dataset:
                return jets[deepCSV >= 0.4184]
        else:
            return ak.Array([])
