from coffea.analysis_objects import JaggedCandidateArray
import numpy as np
import awkward1 as ak

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
            softDropMass=df['JetsAK8_softDropMass'].flatten()
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
        ak4QualityCut = (self.jets.pt > 30) & (abs(self.jets.eta) < self.etaCut)
        return self.jets[ak4QualityCut]

    def goodFatJets(self):
        # # Good AK8 Jets Cut
        ak8QualityCut = (self.fjets.pt > 200) & (abs(self.fjets.eta) < self.etaCut)
        return self.fjets[ak8QualityCut]

    def goodgfjets(self):
        # # Good AK8 Jets Cut
        ak8QualityCut = (self.gfjets.pt > 200) & (abs(self.gfjets.eta) < self.etaCut)
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
