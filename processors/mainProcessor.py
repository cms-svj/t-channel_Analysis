from coffea import hist, processor
import numpy as np
import awkward1 as ak
import awkward
from utils import utility as utl
from utils.objects import Objects
from utils import baseline as bl
from itertools import combinations

class MainProcessor(processor.ProcessorABC):
        def __init__(self):
                self._accumulator = processor.dict_accumulator({})
                self.setupHistos = None
                #         'cutflow':                          processor.defaultdict_accumulator(int),
        @property
        def accumulator(self):
                return self._accumulator

        def setupHistogram(self,cuts):
                histograms = {}
                histograms['EventCounter']  = hist.Hist("EventCounter",     hist.Bin("EventCounter",          "EventCounter",                                     2,      -1.1,   1.1 ))
                histograms['h_weight']      = hist.Hist("h_weight",         hist.Bin("Weight",                "Weight",                                           2,      -1.1,   1.1))
                for name,cut in cuts.items():
                    histograms['h_ht'+name]                 = hist.Hist("h_ht"+name,                hist.Bin("ht",                    r"$H_{T}$ (GeV)",                                   500,    0.0,    5000.0))
                    histograms['h_st'+name]                 = hist.Hist("h_st"+name,                hist.Bin("st",                    r"$S_{T}$ (GeV)",                                   500,    0.0,    5000.0))
                    histograms['h_njets'+name]              = hist.Hist('h_njets'+name,             hist.Bin("njets",                 "Number of Jets",                                   20,     0.0,    20.0))
                    histograms['h_njetsAK8'+name]           = hist.Hist('h_njetsAK8'+name,          hist.Bin("njets",                 "Number of Jets",                                   20,     0.0,    20.0))
                    histograms['h_nb'+name]                 = hist.Hist('h_nb'+name,                hist.Bin("nb",                    "Number of b",                                      10,     0.0,    10.0))
                    histograms['h_nl'+name]                 = hist.Hist('h_nl'+name,                hist.Bin("nl",                    "Number of Leptons",                                10,     0.0,    10.0))
                    histograms['h_met'+name]                = hist.Hist('h_met'+name,               hist.Bin("MET",                   "MET [GeV]",                                        500,    0.0,    2000.0))
                    histograms['h_jPt'+name]                = hist.Hist('h_jPt'+name,               hist.Bin("pt",                    r"$p_{T}$ [GeV]",                                   200,    0.0,    2000.0))
                    histograms['h_jEta'+name]               = hist.Hist('h_jEta'+name,              hist.Bin("eta",                   r"$\eta$",                                          200,   -6.0,    6.0))
                    histograms['h_jPhi'+name]               = hist.Hist('h_jPhi'+name,              hist.Bin("phi",                   r"$\phi$",                                          200,   -4.0,    4.0))
                    histograms['h_jAxismajor'+name]         = hist.Hist('h_jAxismajor'+name,        hist.Bin("axismajor",             r"$\sigma_{major}(j)$",                             40,     0.0,    0.5))
                    histograms['h_jAxisminor'+name]         = hist.Hist('h_jAxisminor'+name,        hist.Bin("axisminor",             r"$\sigma_{minor}(j)$",                             40,     0.0,    0.3))
                    histograms['h_jPtD'+name]               = hist.Hist('h_jPtD'+name,              hist.Bin("ptD",                   "ptD",                                              40,     0.0,    1.2))
                    histograms['h_dPhiMinjMET'+name]        = hist.Hist('h_dPhiMinjMET'+name,       hist.Bin("dPhiJMET",              r"$\Delta\phi_{min}(j,MET)$",                       100,    0.0,    4.0))
                    histograms['h_dEtaj12'+name]            = hist.Hist('h_dEtaj12'+name,           hist.Bin("dEtaJ12",               r"$\Delta\eta(J_{1},J_{2})$",                       200,    0.0,    10.0))
                    histograms['h_dRJ12'+name]              = hist.Hist('h_dRJ12'+name,             hist.Bin("dRJ12",                 r"$\Delta R(J_{1},J_{2})$",                         100,    0.0,   10.0))
                    histograms['h_jPtAK8'+name]             = hist.Hist('h_jPtAK8'+name,            hist.Bin("pt",                    r"$p_{T}$ [GeV]",                                   200,    0.0,    2000.0))
                    histograms['h_jEtaAK8'+name]            = hist.Hist('h_jEtaAK8'+name,           hist.Bin("eta",                   r"$\eta$",                                          200,   -6.0,    6.0))
                    histograms['h_jPhiAK8'+name]            = hist.Hist('h_jPhiAK8'+name,           hist.Bin("phi",                   r"$\phi$",                                          200,   -4.0,    4.0))
                    histograms['h_jAxismajorAK8'+name]      = hist.Hist('h_jAxismajorAK8'+name,     hist.Bin("axismajor",             r"$\sigma_{major}(j)$",                             40,     0.0,    0.5))
                    histograms['h_jAxisminorAK8'+name]      = hist.Hist('h_jAxisminorAK8'+name,     hist.Bin("axisminor",             r"$\sigma_{minor}(j)$",                             40,     0.0,    0.3))
                    histograms['h_jGirthAK8'+name]          = hist.Hist('h_jGirthAK8'+name,         hist.Bin("girth",                 "girth(j)",                                         40,     0.0,    0.5))
                    histograms['h_jPtDAK8'+name]            = hist.Hist('h_jPtDAK8'+name,           hist.Bin("ptD",                   "ptD",                                              40,     0.0,    1.2))
                    histograms['h_jTau1AK8'+name]           = hist.Hist('h_jTau1AK8'+name,          hist.Bin("tau1",                  r"$\tau_{1}(j)$",                                   40,     0.0,    0.8))
                    histograms['h_jTau2AK8'+name]           = hist.Hist('h_jTau2AK8'+name,          hist.Bin("tau2",                  r"$\tau_{2}(j)$",                                   40,     0.0,    0.65))
                    histograms['h_jTau3AK8'+name]           = hist.Hist('h_jTau3AK8'+name,          hist.Bin("tau3",                  r"$\tau_{3}(j)$",                                   40,     0.0,    0.35))
                    histograms['h_jTau21AK8'+name]          = hist.Hist('h_jTau21AK8'+name,         hist.Bin("tau21",                 r"$	au_{21}(j)$",                                 40,     0.0,    1.3))
                    histograms['h_jTau32AK8'+name]          = hist.Hist('h_jTau32AK8'+name,         hist.Bin("tau32",                 r"$	au_{32}(j)$",                                 40,     0.0,    1.3))
                    histograms['h_jSoftDropMassAK8'+name]   = hist.Hist('h_jSoftDropMassAK8'+name,  hist.Bin("softDropMass",          r"$m_{SD}(j)$",                                     40,     0.0,    200))
                    histograms['h_dPhiMinjMETAK8'+name]     = hist.Hist('h_dPhiMinjMETAK8'+name,    hist.Bin("dPhiJMET",              r"$\Delta\phi_{min}(j,MET)$",                       100,    0.0,    4.0))
                    histograms['h_dEtaj12AK8'+name]         = hist.Hist('h_dEtaj12AK8'+name,        hist.Bin("dEtaJ12",               r"$\Delta\eta(J_{1},J_{2})$",                       200,    0.0,    10.0))
                    histograms['h_dRJ12AK8'+name]           = hist.Hist('h_dRJ12AK8'+name,          hist.Bin("dRJ12",                 r"$\Delta R(J_{1},J_{2})$",                         100,    0.0,   10.0))
                    histograms['h_mT'+name]                 = hist.Hist('h_mT'+name,                hist.Bin("mT",                    r"$m_{T} (GeV)$",                                   500,    0.0,    5000.0))
                    histograms['h_METrHT_pt30'+name]        = hist.Hist('h_METrHT_pt30'+name,       hist.Bin("METrHT_pt30",           r"$MET/H_{T}$",                                     100,    0.0,    3.0))
                    histograms['h_METrST_pt30'+name]        = hist.Hist('h_METrST_pt30'+name,       hist.Bin("METrST_pt30",           r"$MET/S_{T}",                                      100,    0.0,    1.0))
                    histograms['h_j1Pt'+name]               = hist.Hist('h_j1Pt'+name,              hist.Bin("pt",                    r"$p_{T}$ [GeV]",                                   200,    0.0,    2000.0))
                    histograms['h_j1Eta'+name]              = hist.Hist('h_j1Eta'+name,             hist.Bin("eta",                   r"$\eta$",                                          200,   -6.0,    6.0))
                    histograms['h_j1Phi'+name]              = hist.Hist('h_j1Phi'+name,             hist.Bin("phi",                   r"$\phi$",                                          200,   -4.0,    4.0))
                    histograms['h_j1Axismajor'+name]        = hist.Hist('h_j1Axismajor'+name,       hist.Bin("axismajor",             r"$\sigma_{major}(j)$",                             40,     0.0,    0.5))
                    histograms['h_j1Axisminor'+name]        = hist.Hist('h_j1Axisminor'+name,       hist.Bin("axisminor",             r"$\sigma_{minor}(j)$",                             40,     0.0,    0.3))
                    histograms['h_j1PtD'+name]              = hist.Hist('h_j1PtD'+name,             hist.Bin("ptD",                   "ptD",                                              40,     0.0,    1.2))
                    histograms['h_dPhij1MET'+name]          = hist.Hist('h_dPhij1MET'+name,         hist.Bin("dPhiJMET",              r"$\Delta\phi(J_{1},MET)$",                         100,    0.0,    4.0))
                    histograms['h_j2Pt'+name]               = hist.Hist('h_j2Pt'+name,              hist.Bin("pt",                    r"$p_{T}$ [GeV]",                                   200,    0.0,    2000.0))
                    histograms['h_j2Eta'+name]              = hist.Hist('h_j2Eta'+name,             hist.Bin("eta",                   r"$\eta$",                                          200,   -6.0,    6.0))
                    histograms['h_j2Phi'+name]              = hist.Hist('h_j2Phi'+name,             hist.Bin("phi",                   r"$\phi$",                                          200,   -4.0,    4.0))
                    histograms['h_j2Axismajor'+name]        = hist.Hist('h_j2Axismajor'+name,       hist.Bin("axismajor",             r"$\sigma_{major}(j)$",                             40,     0.0,    0.5))
                    histograms['h_j2Axisminor'+name]        = hist.Hist('h_j2Axisminor'+name,       hist.Bin("axisminor",             r"$\sigma_{minor}(j)$",                             40,     0.0,    0.3))
                    histograms['h_j2PtD'+name]              = hist.Hist('h_j2PtD'+name,             hist.Bin("ptD",                   "ptD",                                              40,     0.0,    1.2))
                    histograms['h_dPhij2MET'+name]          = hist.Hist('h_dPhij2MET'+name,         hist.Bin("dPhiJMET",              r"$\Delta\phi(J_{1},MET)$",                         100,    0.0,    4.0))
                    histograms['h_dPhij1rdPhij2'+name]      = hist.Hist('h_dPhij1rdPhij2'+name,     hist.Bin("dPhiJ1METrdPhiJ2MET",   r"$\Delta\phi(J_{1},MET)/\Delta\phi(J_{2},MET)$",   100,    0.0,    100.0))
                    histograms['h_j1PtAK8'+name]            = hist.Hist('h_j1PtAK8'+name,           hist.Bin("pt",                    r"$p_{T}$ [GeV]",                                   200,    0.0,    2000.0))
                    histograms['h_j1EtaAK8'+name]           = hist.Hist('h_j1EtaAK8'+name,          hist.Bin("eta",                   r"$\eta$",                                          200,   -6.0,    6.0))
                    histograms['h_j1PhiAK8'+name]           = hist.Hist('h_j1PhiAK8'+name,          hist.Bin("phi",                   r"$\phi$",                                          200,   -4.0,    4.0))
                    histograms['h_j1AxismajorAK8'+name]     = hist.Hist('h_j1AxismajorAK8'+name,    hist.Bin("axismajor",             r"$\sigma_{major}(j)$",                             40,     0.0,    0.5))
                    histograms['h_j1AxisminorAK8'+name]     = hist.Hist('h_j1AxisminorAK8'+name,    hist.Bin("axisminor",             r"$\sigma_{minor}(j)$",                             40,     0.0,    0.3))
                    histograms['h_j1GirthAK8'+name]         = hist.Hist('h_j1GirthAK8'+name,        hist.Bin("girth",                 "girth(j)",                                         40,     0.0,    0.5))
                    histograms['h_j1PtDAK8'+name]           = hist.Hist('h_j1PtDAK8'+name,          hist.Bin("ptD",                   "ptD",                                              40,     0.0,    1.2))
                    histograms['h_j1Tau1AK8'+name]          = hist.Hist('h_j1Tau1AK8'+name,         hist.Bin("tau1",                  r"$	au_{1}(j)$",                                  40,     0.0,    0.8))
                    histograms['h_j1Tau2AK8'+name]          = hist.Hist('h_j1Tau2AK8'+name,         hist.Bin("tau2",                  r"$	au_{2}(j)$",                                  40,     0.0,    0.65))
                    histograms['h_j1Tau3AK8'+name]          = hist.Hist('h_j1Tau3AK8'+name,         hist.Bin("tau3",                  r"$	au_{3}(j)$",                                  40,     0.0,    0.35))
                    histograms['h_j1Tau21AK8'+name]         = hist.Hist('h_j1Tau21AK8'+name,        hist.Bin("tau21",                 r"$	au_{21}(j)$",                                 40,     0.0,    1.3))
                    histograms['h_j1Tau32AK8'+name]         = hist.Hist('h_j1Tau32AK8'+name,        hist.Bin("tau32",                 r"$	au_{32}(j)$",                                 40,     0.0,    1.3))
                    histograms['h_j1SoftDropMassAK8'+name]  = hist.Hist('h_j1SoftDropMassAK8'+name, hist.Bin("softDropMass",          r"$m_{SD}(j)$",                                     40,     0.0,    200))
                    histograms['h_dPhij1METAK8'+name]       = hist.Hist('h_dPhij1METAK8'+name,      hist.Bin("dPhiJMET",              r"$\Delta\phi(J_{1},MET)$",                         100,    0.0,    4.0))
                    histograms['h_j2PtAK8'+name]            = hist.Hist('h_j2PtAK8'+name,           hist.Bin("pt",                    r"$p_{T}$ [GeV]",                                   200,    0.0,    2000.0))
                    histograms['h_j2EtaAK8'+name]           = hist.Hist('h_j2EtaAK8'+name,          hist.Bin("eta",                   r"$\eta$",                                          200,   -6.0,    6.0))
                    histograms['h_j2PhiAK8'+name]           = hist.Hist('h_j2PhiAK8'+name,          hist.Bin("phi",                   r"$\phi$",                                          200,   -4.0,    4.0))
                    histograms['h_j2AxismajorAK8'+name]     = hist.Hist('h_j2AxismajorAK8'+name,    hist.Bin("axismajor",             r"$\sigma_{major}(j)$",                             40,     0.0,    0.5))
                    histograms['h_j2AxisminorAK8'+name]     = hist.Hist('h_j2AxisminorAK8'+name,    hist.Bin("axisminor",             r"$\sigma_{minor}(j)$",                             40,     0.0,    0.3))
                    histograms['h_j2GirthAK8'+name]         = hist.Hist('h_j2GirthAK8'+name,        hist.Bin("girth",                 "girth(j)",                                         40,     0.0,    0.5))
                    histograms['h_j2PtDAK8'+name]           = hist.Hist('h_j2PtDAK8'+name,          hist.Bin("ptD",                   "ptD",                                              40,     0.0,    1.2))
                    histograms['h_j2Tau1AK8'+name]          = hist.Hist('h_j2Tau1AK8'+name,         hist.Bin("tau1",                  r"$	au_{1}(j)$",                                  40,     0.0,    0.8))
                    histograms['h_j2Tau2AK8'+name]          = hist.Hist('h_j2Tau2AK8'+name,         hist.Bin("tau2",                  r"$	au_{2}(j)$",                                  40,     0.0,    0.65))
                    histograms['h_j2Tau3AK8'+name]          = hist.Hist('h_j2Tau3AK8'+name,         hist.Bin("tau3",                  r"$	au_{3}(j)$",                                  40,     0.0,    0.35))
                    histograms['h_j2Tau21AK8'+name]         = hist.Hist('h_j2Tau21AK8'+name,        hist.Bin("tau21",                 r"$	au_{21}(j)$",                                 40,     0.0,    1.3))
                    histograms['h_j2Tau32AK8'+name]         = hist.Hist('h_j2Tau32AK8'+name,        hist.Bin("tau32",                 r"$	au_{32}(j)$",                                 40,     0.0,    1.3))
                    histograms['h_j2SoftDropMassAK8'+name]  = hist.Hist('h_j2SoftDropMassAK8'+name, hist.Bin("softDropMass",          r"$m_{SD}(j)$",                                     40,     0.0,    200))
                    histograms['h_dPhij2METAK8'+name]       = hist.Hist('h_dPhij2METAK8'+name,      hist.Bin("dPhiJMET",              r"$\Delta\phi(J_{1},MET)$",                         100,    0.0,    4.0))
                    histograms['h_dPhij1rdPhij2AK8'+name]   = hist.Hist('h_dPhij1rdPhij2AK8'+name,  hist.Bin("dPhiJ1METrdPhiJ2MET",   r"$\Delta\phi(J_{1},MET)/\Delta\phi(J_{2},MET)$",   100,    0.0,    100.0))
                    histograms['h_madHT'+name]              = hist.Hist("h_madHT"+name,             hist.Bin("ht",                    r"$H_{T}$ (GeV)",                                   500,    0.0,    5000.0))
                self._accumulator = processor.dict_accumulator(histograms)
                self.setupHistos = True

        def process(self, df):
                # set up event counter: useful for checking that we ran over the correct samples
                luminosity = 21071.0+38654.0
                eventWeight = luminosity*df['Weight']
                eCounter = np.where(eventWeight >= 0, 1, -1)

                # cut loop
                ## objects used for cuts
                obj = Objects(df)
                electrons_noCut = obj.goodElectrons()
                muons_noCut = obj.goodMuons()
                jets_noCut = obj.goodJets()
                bjets_noCut = obj.goodBJets(df,jets_noCut)
                if len(bjets_noCut) > 0:
                    nBJets_noCut = bjets_noCut.counts
                else:
                    nBJets_noCut = np.zeros(len(eventWeight))
                fjets_noCut = obj.goodFatJets()
                ## variables used for cuts
                metPhi_noCut = df['METPhi']
                ht_noCut = ak.sum(jets_noCut.pt,axis=1)
                dPhiMinj_noCut = utl.deltaPhi(jets_noCut.phi,metPhi_noCut).min()
                dPhiMinjAK8_noCut = utl.deltaPhi(fjets_noCut.phi,metPhi_noCut).min()
                ## defining cuts
                ttStitch = bl.ttStitch(df)
                metFilters = bl.METFilters(df)
                triggerCut = bl.passTrigger(df['TriggerPass'])
                psFilter = bl.phiSpikeFilter(df,jets_noCut)
                qualityCuts = ttStitch & metFilters & psFilter & triggerCut
                preselection = bl.preselection(qualityCuts,electrons_noCut,muons_noCut,df['MET'])
                ht_Cut = ht_noCut > 400
                dPhiMinj_Cut = dPhiMinj_noCut < 0.65
                dPhiMinjAK8_Cut = dPhiMinjAK8_noCut < 0.65
                # Cuts from theory paper
                # cuts = {
                #         ""                          : np.ones(len(df["Weight"]),dtype=bool),
                #         "_trigPresAK4"              : (utl.jetVar_vec(jets_noCut.pt,0) > 250) & (df['MET'] > 200),
                #         "_trigPresAK4_MET"          : (utl.jetVar_vec(jets_noCut.pt,0) > 250) & (df['MET'] > 200) & (df['MET'] > 800),
                #         "_trigPresAK4_MET_dPhiG"    : (utl.jetVar_vec(jets_noCut.pt,0) > 250) & (df['MET'] > 200) & (df['MET'] > 800) & (dPhiMinj_noCut > 0.4),
                #         "_trigPresAK4_MET_dPhiL"    : (utl.jetVar_vec(jets_noCut.pt,0) > 250) & (df['MET'] > 200) & (df['MET'] > 800) & (dPhiMinj_noCut < 0.4),
                #         "_trigPresAK4_dPhiL"        : (utl.jetVar_vec(jets_noCut.pt,0) > 250) & (df['MET'] > 200) & (dPhiMinj_noCut < 0.4),
                #         "_trigPresAK8"              : (utl.jetVar_vec(fjets_noCut.pt,0) > 250) & (df['MET'] > 200),
                #         "_trigPresAK8_MET"          : (utl.jetVar_vec(fjets_noCut.pt,0) > 250) & (df['MET'] > 200) & (df['MET'] > 800),
                #         "_trigPresAK8_MET_dPhiG"    : (utl.jetVar_vec(fjets_noCut.pt,0) > 250) & (df['MET'] > 200) & (df['MET'] > 800) & (dPhiMinjAK8_noCut > 0.4),
                #         "_trigPresAK8_MET_dPhiL"    : (utl.jetVar_vec(fjets_noCut.pt,0) > 250) & (df['MET'] > 200) & (df['MET'] > 800) & (dPhiMinjAK8_noCut < 0.4),
                #         "_trigPresAK8_dPhiL"        : (utl.jetVar_vec(fjets_noCut.pt,0) > 250) & (df['MET'] > 200) & (dPhiMinjAK8_noCut < 0.4),
                # }
                # Our preselection
                cuts = {
                        ""                          : np.ones(len(df["Weight"]),dtype=bool),
                        "_qc"                       : ttStitch & metFilters & psFilter,
                        "_qc_trg"                   : ttStitch & metFilters & psFilter & triggerCut,
                        "_qc_trg_MET300"            : ttStitch & metFilters & psFilter & triggerCut & (df["MET"] > 300), # trigger turn on at MET > 250
                        "_qc_trg_ht1300"            : ttStitch & metFilters & psFilter & triggerCut & (ht_noCut > 1300), # trigger turn on at hT > 1300
                        "_qc_trg_ht1300_MET300"     : ttStitch & metFilters & psFilter & triggerCut & (df["MET"] > 300) & (ht_noCut > 1300), # trigger turn on at hT > 1300
                        "_qc_trg_0l"                : ttStitch & metFilters & psFilter & triggerCut & (electrons_noCut.counts == 0) & (muons_noCut.counts == 0),
                        "_pre"                      : preselection,
                        "_pre_ge2AK4j"              : preselection & (jets_noCut.counts >= 2) & (df["JetID"] == True),
                        "_pre_ge2AK4j_MET220"       : preselection & (jets_noCut.counts >= 2) & (df["JetID"] == True) & (df["MET"] > 220), # 50% relative signal efficiency
                        "_pre_ge2AK4j_MET800"       : preselection & (jets_noCut.counts >= 2) & (df["JetID"] == True) & (df["MET"] > 800), # theory paper MET cut
                        "_pre_ge2AK4j_ht400"        : preselection & (jets_noCut.counts >= 2) & (df["JetID"] == True) & ht_Cut,
                        "_pre_ge2AK4j_ht750"        : preselection & (jets_noCut.counts >= 2) & (df["JetID"] == True) & (ht_noCut > 750), # 50% relative signal efficiency
                        "_pre_ge2AK4j_ht1300"       : preselection & (jets_noCut.counts >= 2) & (df["JetID"] == True) & (ht_noCut > 1300), # trigger turn on at hT > 1300
                        "_pre_ge2AK4j_nb2"          : preselection & (jets_noCut.counts >= 2) & (df["JetID"] == True) & (nBJets_noCut >= 2),
                        "_pre_ge2AK4j_nb3"          : preselection & (jets_noCut.counts >= 2) & (df["JetID"] == True) & (nBJets_noCut >= 3),
                        "_pre_ge2AK4j_nb4"          : preselection & (jets_noCut.counts >= 2) & (df["JetID"] == True) & (nBJets_noCut >= 4),
                        "_pre_ge2AK4j_dpjp65"       : preselection & (jets_noCut.counts >= 2) & (df["JetID"] == True) & dPhiMinj_Cut,
                        "_pre_ge2AK4j_ht400_dpjp65" : preselection & (jets_noCut.counts >= 2) & (df["JetID"] == True) & ht_Cut & dPhiMinj_Cut,
                        "_pre_ge2AK8j"              : preselection & (fjets_noCut.counts >= 2) & (df["JetIDAK8"] == True),
                        "_pre_ge2AK8j_MET220"       : preselection & (fjets_noCut.counts >= 2) & (df["JetIDAK8"] == True) & (df["MET"] > 220), # 50% relative signal efficiency
                        "_pre_ge2AK8j_MET800"       : preselection & (fjets_noCut.counts >= 2) & (df["JetIDAK8"] == True) & (df["MET"] > 800), # theory paper MET cut
                        "_pre_ge2AK8j_ht400"        : preselection & (fjets_noCut.counts >= 2) & (df["JetIDAK8"] == True) & ht_Cut,
                        "_pre_ge2AK8j_ht750"        : preselection & (fjets_noCut.counts >= 2) & (df["JetIDAK8"] == True) & (ht_noCut > 750), # 50% relative signal efficiency
                        "_pre_ge2AK8j_ht1300"       : preselection & (fjets_noCut.counts >= 2) & (df["JetIDAK8"] == True) & (ht_noCut > 1300), # trigger turn on at hT > 1300
                        "_pre_ge2AK8j_nb2"          : preselection & (fjets_noCut.counts >= 2) & (df["JetIDAK8"] == True) & (nBJets_noCut >= 2),
                        "_pre_ge2AK8j_nb3"          : preselection & (fjets_noCut.counts >= 2) & (df["JetIDAK8"] == True) & (nBJets_noCut >= 3),
                        "_pre_ge2AK8j_nb4"          : preselection & (fjets_noCut.counts >= 2) & (df["JetIDAK8"] == True) & (nBJets_noCut >= 4),
                        "_pre_ge2AK8j_dpJp65"       : preselection & (fjets_noCut.counts >= 2) & (df["JetIDAK8"] == True) & dPhiMinjAK8_Cut,
                        "_pre_ge2AK8j_ht400__dpJp65": preselection & (fjets_noCut.counts >= 2) & (df["JetIDAK8"] == True) & ht_Cut & dPhiMinjAK8_Cut
                }
# add cutflow (see coffea cutflow)
# go through the paper to get cutflows, get the cuts for t-channel from the paper
# try to generate tchannel production signals: don't use madgraph flag, use t-channel, don't need LHE

                # setup histograms
                if self.setupHistos is None:
                    self.setupHistogram(cuts)
                output = self.accumulator.identity()

                # fill initial histograms
                if len(eCounter) > 0:
                    output["EventCounter"].fill(EventCounter=eCounter,weight=np.ones(len(eCounter)))
                    output["h_weight"].fill(Weight=eCounter,weight=eventWeight)

                # run cut loop
                for name,cut in cuts.items():
                    # defining objects
                    electrons = electrons_noCut[cut]
                    muons = muons_noCut[cut]
                    jets = jets_noCut[cut]
                    fjets = fjets_noCut[cut]
                    bjets = obj.goodBJets(df,jets)
                    madHT = df['madHT'][cut]
                    genMET = df['GenMET'][cut]
                    met = df['MET'][cut]
                    metPhi = metPhi_noCut[cut]
                    mtAK8 = df['MT_AK8'][cut]
                    triggerPass = df['TriggerPass'][cut]
                    evtw = eventWeight[cut]
                    if len(bjets) > 0:
                        nBJets = bjets.counts
                    else:
                        nBJets = np.zeros(len(evtw))
                    # defining weights for awkward arrays
                    ew = utl.awkwardReshape(electrons,evtw)
                    mw = utl.awkwardReshape(muons,evtw)
                    jw = utl.awkwardReshape(jets,evtw)
                    fjw = utl.awkwardReshape(fjets,evtw)
                    # bjw = utl.awkwardReshape(bjets,evtw)

                    if len(evtw) > 0:
                        # Getting subset of variables based on number of AK8 jets
                        # calculating variables

                        ht = ak.sum(jets.pt,axis=1)
                        st = ht + met
                        metrht = utl.divide_vec(met,ht)
                        metrst = utl.divide_vec(met,st)

                        # if name == "_pre_ge2AK4j_ht400":
                        #     htTruth = ht>400
                        #     allTrue = True
                        #     for t in htTruth:
                        #         if t == False:
                        #             allTrue = 0
                        #             break
                        #     if allTrue == False:
                        #         print(ht)
                        #         print(htTruth)
                        # AK4 Jet Variables
                        jetPhi = jets.phi
                        jetEta = jets.eta
                        j1_eta = utl.jetVar_vec(jetEta,0)
                        j2_eta = utl.jetVar_vec(jetEta,1)
                        j1_phi = utl.jetVar_vec(jetPhi,0)
                        j2_phi = utl.jetVar_vec(jetPhi,1)
                        dPhij1 = utl.deltaPhiji_vec(j1_phi,metPhi)
                        dPhij2 = utl.deltaPhiji_vec(j2_phi,metPhi)
                        dPhij1rdPhij2 = utl.divide_vec(dPhij1,dPhij2)
                        dPhiMinj = utl.deltaPhi(jetPhi,metPhi).min()
                        dEtaj12 = utl.deltaEta_vec(j1_eta,j2_eta)
                        deltaR12j = utl.delta_R(j1_eta,j2_eta,j1_phi,j2_phi)

                        # AK8 Jet Variables
                        jetAK8Phi = fjets.phi
                        jetAK8Eta = fjets.eta
                        j1_etaAK8 = utl.jetVar_vec(jetAK8Eta,0)
                        j2_etaAK8 = utl.jetVar_vec(jetAK8Eta,1)
                        j1_phiAK8 = utl.jetVar_vec(jetAK8Phi,0)
                        j2_phiAK8 = utl.jetVar_vec(jetAK8Phi,1)
                        dPhij1AK8 = utl.deltaPhiji_vec(j1_phiAK8,metPhi)
                        dPhij2AK8 = utl.deltaPhiji_vec(j2_phiAK8,metPhi)
                        dPhij1rdPhij2AK8 = utl.divide_vec(dPhij1AK8,dPhij2AK8)
                        dPhiMinjAK8 = utl.deltaPhi(jetAK8Phi,metPhi).min()
                        dEtaj12AK8 = utl.deltaEta_vec(j1_etaAK8,j2_etaAK8)
                        deltaR12jAK8 = utl.delta_R(j1_etaAK8,j2_etaAK8,j1_phiAK8,j2_phiAK8)

                        tau1 = fjets.tau1
                        tau2 = fjets.tau2
                        tau3 = fjets.tau3
                        J_tau21 = utl.divide_vec(tau2.flatten(),tau1.flatten())
                        J_tau32 = utl.divide_vec(tau3.flatten(),tau2.flatten())
                        J1_tau21 = utl.tauRatio(tau2,tau1,0)
                        J1_tau32 = utl.tauRatio(tau3,tau2,0)
                        J2_tau21 = utl.tauRatio(tau2,tau1,1)
                        J2_tau32 = utl.tauRatio(tau3,tau2,1)

                        # filling histograms
                        output['h_njets'+name].fill(njets=jets.counts,weight=evtw)
                        output['h_njetsAK8'+name].fill(njets=fjets.counts,weight=evtw)
                        output['h_nb'+name].fill(nb=nBJets,weight=evtw)
                        output['h_nl'+name].fill(nl=(electrons.counts + muons.counts),weight=evtw)
                        output['h_ht'+name].fill(ht=ht,weight=evtw)
                        output['h_st'+name].fill(st=st,weight=evtw)
                        output['h_met'+name].fill(MET=met,weight=evtw)
                        output['h_madHT'+name].fill(ht=madHT,weight=evtw)
                        output['h_jPt'+name].fill(pt=jets.pt.flatten(),weight=ak.flatten(jw))
                        output['h_jEta'+name].fill(eta=jetEta.flatten(),weight=ak.flatten(jw))
                        output['h_jPhi'+name].fill(phi=jetPhi.flatten(),weight=ak.flatten(jw))
                        output['h_jAxismajor'+name].fill(axismajor=jets.axismajor.flatten(),weight=ak.flatten(jw))
                        output['h_jAxisminor'+name].fill(axisminor=jets.axisminor.flatten(),weight=ak.flatten(jw))
                        output['h_jPtD'+name].fill(ptD=jets.ptD.flatten(),weight=ak.flatten(jw))
                        output['h_dPhiMinjMET'+name].fill(dPhiJMET=dPhiMinj,weight=evtw)
                        output['h_dEtaj12'+name].fill(dEtaJ12=dEtaj12,weight=evtw)
                        output['h_dRJ12'+name].fill(dRJ12=deltaR12j,weight=evtw)
                        output['h_jPtAK8'+name].fill(pt=fjets.pt.flatten(),weight=ak.flatten(fjw))
                        output['h_jEtaAK8'+name].fill(eta=jetAK8Eta.flatten(),weight=ak.flatten(fjw))
                        output['h_jPhiAK8'+name].fill(phi=jetAK8Phi.flatten(),weight=ak.flatten(fjw))
                        output['h_jAxismajorAK8'+name].fill(axismajor=fjets.axismajor.flatten(),weight=ak.flatten(fjw))
                        output['h_jAxisminorAK8'+name].fill(axisminor=fjets.axisminor.flatten(),weight=ak.flatten(fjw))
                        output['h_jGirthAK8'+name].fill(girth=fjets.girth.flatten(),weight=ak.flatten(fjw))
                        output['h_jPtDAK8'+name].fill(ptD=fjets.ptD.flatten(),weight=ak.flatten(fjw))
                        output['h_jTau1AK8'+name].fill(tau1=tau1.flatten(),weight=ak.flatten(fjw))
                        output['h_jTau2AK8'+name].fill(tau2=tau2.flatten(),weight=ak.flatten(fjw))
                        output['h_jTau3AK8'+name].fill(tau3=tau3.flatten(),weight=ak.flatten(fjw))
                        output['h_jTau21AK8'+name].fill(tau21=J_tau21,weight=ak.flatten(fjw))
                        output['h_jTau32AK8'+name].fill(tau32=J_tau32,weight=ak.flatten(fjw))
                        output['h_jSoftDropMassAK8'+name].fill(softDropMass=fjets.softDropMass.flatten(),weight=ak.flatten(fjw))
                        output['h_dPhiMinjMETAK8'+name].fill(dPhiJMET=dPhiMinjAK8,weight=evtw)
                        output['h_dEtaj12AK8'+name].fill(dEtaJ12=dEtaj12AK8,weight=evtw)
                        output['h_dRJ12AK8'+name].fill(dRJ12=deltaR12jAK8,weight=evtw)
                        output['h_mT'+name].fill(mT=mtAK8,weight=evtw)
                        output['h_METrHT_pt30'+name].fill(METrHT_pt30=metrht,weight=evtw)
                        output['h_METrST_pt30'+name].fill(METrST_pt30=metrst,weight=evtw)
                        output['h_j1Pt'+name].fill(pt=utl.jetVar_vec(jets.pt,0),weight=evtw)
                        output['h_j1Eta'+name].fill(eta=j1_eta,weight=evtw)
                        output['h_j1Phi'+name].fill(phi=j1_phi,weight=evtw)
                        output['h_j1Axismajor'+name].fill(axismajor=utl.jetVar_vec(jets.axismajor,0),weight=evtw)
                        output['h_j1Axisminor'+name].fill(axisminor=utl.jetVar_vec(jets.axisminor,0),weight=evtw)
                        output['h_j1PtD'+name].fill(ptD=utl.jetVar_vec(jets.ptD,0),weight=evtw)
                        output['h_dPhij1MET'+name].fill(dPhiJMET=dPhij1,weight=evtw)
                        output['h_j2Pt'+name].fill(pt=utl.jetVar_vec(jets.pt,1),weight=evtw)
                        output['h_j2Eta'+name].fill(eta=j2_eta,weight=evtw)
                        output['h_j2Phi'+name].fill(phi=j2_phi,weight=evtw)
                        output['h_j2Axismajor'+name].fill(axismajor=utl.jetVar_vec(jets.axismajor,1),weight=evtw)
                        output['h_j2Axisminor'+name].fill(axisminor=utl.jetVar_vec(jets.axisminor,1),weight=evtw)
                        output['h_j2PtD'+name].fill(ptD=utl.jetVar_vec(jets.ptD,1),weight=evtw)
                        output['h_dPhij2MET'+name].fill(dPhiJMET=dPhij2,weight=evtw)
                        output['h_dPhij1rdPhij2'+name].fill(dPhiJ1METrdPhiJ2MET=dPhij1rdPhij2,weight=evtw)
                        output['h_j1PtAK8'+name].fill(pt=utl.jetVar_vec(fjets.pt,0),weight=evtw)
                        output['h_j1EtaAK8'+name].fill(eta=j1_etaAK8,weight=evtw)
                        output['h_j1PhiAK8'+name].fill(phi=j1_phiAK8,weight=evtw)
                        output['h_j1AxismajorAK8'+name].fill(axismajor=utl.jetVar_vec(fjets.axismajor,0),weight=evtw)
                        output['h_j1AxisminorAK8'+name].fill(axisminor=utl.jetVar_vec(fjets.axisminor,0),weight=evtw)
                        output['h_j1GirthAK8'+name].fill(girth=utl.jetVar_vec(fjets.girth,0),weight=evtw)
                        output['h_j1PtDAK8'+name].fill(ptD=utl.jetVar_vec(fjets.ptD,0),weight=evtw)
                        output['h_j1Tau1AK8'+name].fill(tau1=utl.jetVar_vec(fjets.tau1,0),weight=evtw)
                        output['h_j1Tau2AK8'+name].fill(tau2=utl.jetVar_vec(fjets.tau2,0),weight=evtw)
                        output['h_j1Tau3AK8'+name].fill(tau3=utl.jetVar_vec(fjets.tau3,0),weight=evtw)
                        output['h_j1Tau21AK8'+name].fill(tau21=J1_tau21,weight=evtw)
                        output['h_j1Tau32AK8'+name].fill(tau32=J1_tau32,weight=evtw)
                        output['h_j1SoftDropMassAK8'+name].fill(softDropMass=utl.jetVar_vec(fjets.softDropMass,0),weight=evtw)
                        output['h_dPhij1METAK8'+name].fill(dPhiJMET=dPhij1AK8,weight=evtw)
                        output['h_j2PtAK8'+name].fill(pt=utl.jetVar_vec(fjets.pt,1),weight=evtw)
                        output['h_j2EtaAK8'+name].fill(eta=j2_etaAK8,weight=evtw)
                        output['h_j2PhiAK8'+name].fill(phi=j2_phiAK8,weight=evtw)
                        output['h_j2AxismajorAK8'+name].fill(axismajor=utl.jetVar_vec(fjets.axismajor,1),weight=evtw)
                        output['h_j2AxisminorAK8'+name].fill(axisminor=utl.jetVar_vec(fjets.axisminor,1),weight=evtw)
                        output['h_j2GirthAK8'+name].fill(girth=utl.jetVar_vec(fjets.girth,1),weight=evtw)
                        output['h_j2PtDAK8'+name].fill(ptD=utl.jetVar_vec(fjets.ptD,1),weight=evtw)
                        output['h_j2Tau1AK8'+name].fill(tau1=utl.jetVar_vec(fjets.tau1,1),weight=evtw)
                        output['h_j2Tau2AK8'+name].fill(tau2=utl.jetVar_vec(fjets.tau2,1),weight=evtw)
                        output['h_j2Tau3AK8'+name].fill(tau3=utl.jetVar_vec(fjets.tau3,1),weight=evtw)
                        output['h_j2Tau21AK8'+name].fill(tau21=J2_tau21,weight=evtw)
                        output['h_j2Tau32AK8'+name].fill(tau32=J2_tau32,weight=evtw)
                        output['h_j2SoftDropMassAK8'+name].fill(softDropMass=utl.jetVar_vec(fjets.softDropMass,1),weight=evtw)
                        output['h_dPhij2METAK8'+name].fill(dPhiJMET=dPhij2AK8,weight=evtw)
                        output['h_dPhij1rdPhij2AK8'+name].fill(dPhiJ1METrdPhiJ2MET=dPhij1rdPhij2AK8,weight=evtw)
                #
                    # if len(ak.flatten(tPassedList)) > 0:
                    #     output['h_trigger'].fill(trigger=ak.flatten(tPassedList),weight=np.ones(len(ak.flatten(tPassedList))))
                return output

        def postprocess(self, accumulator):
                return accumulator
