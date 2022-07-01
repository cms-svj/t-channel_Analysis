from itertools import combinations
from coffea import hist

class SVJHistoInfo:
    def __init__(self, varXName, xbins, npzInfo, flattenInfo, weightName, varYName=None, ybins=None):
        self.varXName = varXName
        self.xbins = xbins
        self.npzInfo = npzInfo
        self.flattenInfo = flattenInfo
        self.weightName = weightName
        self.varYName = varYName
        self.ybins = ybins
        self.dim = 2 if varYName else 1

def variables(jNVar=False):
    # [xlabel,number of bins,xmin,xmax,whether to keep it in npz for training (0=do not keep, 1=keep as is, 2=keep but make sure the length is the same as AK8 variables), whether to flatten the array or not when filling histogram (2 = ak.flatten(), 1 = .flatten(), 0 = do not flatten)]
    allVars = {
        'eCounter':                SVJHistoInfo("eCounter",                 hist.Bin("val", "h_eCounter",                                      2,    -1.1,    1.1   ),     0,     0,       'w1'  ),
        'evtw':                    SVJHistoInfo("evtw",                     hist.Bin("val", "h_evtw",                                          2,    -1.1,    1.1   ),     0,     0,       'evtw'),
        'jw':                      SVJHistoInfo("jw",                       hist.Bin("val", "h_jw",                                            2,    -1.1,    1.1   ),     0,     2,       'jw'  ),
        'fjw':                     SVJHistoInfo("fjw",                      hist.Bin("val", "h_fjw",                                           2,    -1.1,    1.1   ),     1,     2,       'fjw' ),
        'njets':                   SVJHistoInfo("njets",                    hist.Bin("val", "Number of Jets",                                 20,     0.0,    20.0  ),     2,     0,       'evtw'),
        'njetsAK8':                SVJHistoInfo("njetsAK8",                 hist.Bin("val", "Number of AK8Jets",                              20,     0.0,    20.0  ),     2,     0,       'evtw'),
        'nb':                      SVJHistoInfo("nb",                       hist.Bin("val", "Number of b",                                    10,     0.0,    10.0  ),     2,     0,       'evtw'),
        'nl':                      SVJHistoInfo("nl",                       hist.Bin("val", "Number of Leptons",                              10,     0.0,    10.0  ),     2,     0,       'evtw'),
        'nnim':                    SVJHistoInfo("nnim",                     hist.Bin("val", "Number of NonIsoMuons",                          10,     0.0,    10.0  ),     2,     0,       'evtw'),
        'ht':                      SVJHistoInfo("ht",                       hist.Bin("val", r"$H_{T}$ (GeV)",                                 500,    0.0,    5000.0),     2,     0,       'evtw'),
        'st':                      SVJHistoInfo("st",                       hist.Bin("val", r"$S_{T}$ (GeV)",                                 500,    0.0,    5000.0),     2,     0,       'evtw'),
        'met':                     SVJHistoInfo("met",                      hist.Bin("val", "MET [GeV]",                                      500,    0.0,    2000.0),     2,     0,       'evtw'),
        'metPhi':                  SVJHistoInfo("metPhi",                   hist.Bin("val", r"MET $\phi$ [GeV]",                              40,    -4.0,    4.0   ),     0,     0,       'evtw'),
        # 'madHT':                 SVJHistoInfo("madHT",                    hist.Bin("val", r"$H_{T}$ (GeV)",                                 500,    0.0,    5000.0),     2,     0,       'evtw'),
        'jPt':                     SVJHistoInfo("jPt",                      hist.Bin("val", r"$p_{T}$ [GeV]",                                 200,    0.0,    2500.0),     0,     1,       'jw'  ),
        'jEta':                    SVJHistoInfo("jEta",                     hist.Bin("val", r"$\eta$",                                        200,   -6.0,    6.0   ),     0,     1,       'jw'  ),
        'jPhi':                    SVJHistoInfo("jPhi",                     hist.Bin("val", r"$\phi$",                                        200,   -4.0,    4.0   ),     0,     1,       'jw'  ),
        'jAxismajor':              SVJHistoInfo("jAxismajor",               hist.Bin("val", r"$\sigma_{major}(j)$",                           40,     0.0,    0.5   ),     0,     1,       'jw'  ),
        'jAxisminor':              SVJHistoInfo("jAxisminor",               hist.Bin("val", r"$\sigma_{minor}(j)$",                           40,     0.0,    0.3   ),     0,     1,       'jw'  ),
        'jPtD':                    SVJHistoInfo("jPtD",                     hist.Bin("val", "ptD",                                            40,     0.0,    1.2   ),     0,     1,       'jw'  ),
        'dPhiMinjMET':             SVJHistoInfo("dPhiMinjMET",              hist.Bin("val", r"$\Delta\phi_{min}(j,MET)$",                     100,    0.0,    4.0   ),     0,     0,       'evtw'),
        'jPtAK8':                  SVJHistoInfo("jPtAK8",                   hist.Bin("val", r"$p_{T}$ [GeV]",                                 280,    0.0,    2800.0),     1,     1,       'fjw' ),
        'jEtaAK8':                 SVJHistoInfo("jEtaAK8",                  hist.Bin("val", r"$\eta$",                                        200,   -6.0,    6.0   ),     1,     1,       'fjw' ),
        'jPhiAK8':                 SVJHistoInfo("jPhiAK8",                  hist.Bin("val", r"$\phi$",                                        200,   -4.0,    4.0   ),     1,     1,       'fjw' ),
        'jAxismajorAK8':           SVJHistoInfo("jAxismajorAK8",            hist.Bin("val", r"$\sigma_{major}(J)$",                           40,     0.0,    0.6   ),     1,     1,       'fjw' ),
        'jAxisminorAK8':           SVJHistoInfo("jAxisminorAK8",            hist.Bin("val", r"$\sigma_{minor}(J)$",                           40,     0.0,    0.4   ),     1,     1,       'fjw' ),
        'jChEMEFractAK8':          SVJHistoInfo("jChEMEFractAK8",           hist.Bin("val", "fChEM(J)",                                       50,     0.0,    1.0   ),     1,     1,       'fjw' ),
        'jChHadEFractAK8':         SVJHistoInfo("jChHadEFractAK8",          hist.Bin("val", "fChHad(J)",                                      50,     0.0,    1.0   ),     1,     1,       'fjw' ),
        'jChHadMultAK8':           SVJHistoInfo("jChHadMultAK8",            hist.Bin("val", "nChHad(J)",                                      145,    0.0,    145.0 ),     1,     1,       'fjw' ),
        'jChMultAK8':              SVJHistoInfo("jChMultAK8",               hist.Bin("val", "nCh(J)",                                         145,    0.0,    145.0 ),     1,     1,       'fjw' ),
        'jdoubleBDiscriminatorAK8':SVJHistoInfo("jdoubleBDiscriminatorAK8", hist.Bin("val", "doubleBDiscriminator(J)",                        100,    -1.0,   1.0   ),     1,     1,       'fjw' ),
        'jecfN2b1AK8':             SVJHistoInfo("jecfN2b1AK8",              hist.Bin("val", "ecfN2b1(J)",                                     50,     0.0,    0.6   ),     1,     1,       'fjw' ),
        'jecfN2b2AK8':             SVJHistoInfo("jecfN2b2AK8",              hist.Bin("val", "ecfN2b2(J)",                                     50,     0.0,    0.4   ),     1,     1,       'fjw' ),
        'jecfN3b1AK8':             SVJHistoInfo("jecfN3b1AK8",              hist.Bin("val", "ecfN3b1(J)",                                     50,     0.0,    6.0   ),     1,     1,       'fjw' ),
        'jecfN3b2AK8':             SVJHistoInfo("jecfN3b2AK8",              hist.Bin("val", "ecfN3b2(J)",                                     50,     0.0,    5.0   ),     1,     1,       'fjw' ),
        'jEleEFractAK8':           SVJHistoInfo("jEleEFractAK8",            hist.Bin("val", "fEle(J)",                                        50,     0.0,    1.0   ),     1,     1,       'fjw' ),
        'jEleMultAK8':             SVJHistoInfo("jEleMultAK8",              hist.Bin("val", "nEle(J)",                                        8,      0.0,    8.0   ),     1,     1,       'fjw' ),
        'jGirthAK8':               SVJHistoInfo("jGirthAK8",                hist.Bin("val", "girth(J)",                                       40,     0.0,    0.7   ),     1,     1,       'fjw' ),
        'jHfEMEFractAK8':          SVJHistoInfo("jHfEMEFractAK8",           hist.Bin("val", "fHFEM(J)",                                       50,     0.0,    1.0   ),     1,     1,       'fjw' ),
        'jHfHadEFractAK8':         SVJHistoInfo("jHfHadEFractAK8",          hist.Bin("val", "fHFHad(J)",                                      50,     0.0,    1.0   ),     1,     1,       'fjw' ),
        'jMultAK8':                SVJHistoInfo("jMultAK8",                 hist.Bin("val", "mult(J)",                                        250,    0.0,    250.0 ),     1,     1,       'fjw' ),
        'jMuEFractAK8':            SVJHistoInfo("jMuEFractAK8",             hist.Bin("val", "fMu(J)",                                         50,     0.0,    1.0   ),     1,     1,       'fjw' ),
        'jMuMultAK8':              SVJHistoInfo("jMuMultAK8",               hist.Bin("val", "nMu(J)",                                         8,      0.0,    10.0  ),     1,     1,       'fjw' ),
        'jNeuEmEFractAK8':         SVJHistoInfo("jNeuEmEFractAK8",          hist.Bin("val", "fNeuEM(J)",                                      50,     0.0,    1.0   ),     1,     1,       'fjw' ),
        'jNeuHadEFractAK8':        SVJHistoInfo("jNeuHadEFractAK8",         hist.Bin("val", "fNeuHad(J)",                                     50,     0.0,    1.0   ),     1,     1,       'fjw' ),
        'jNeuHadMultAK8':          SVJHistoInfo("jNeuHadMultAK8",           hist.Bin("val", "nNeuHad(J)",                                     25,     0.0,    25.0  ),     1,     1,       'fjw' ),
        'jNeuMultAK8':             SVJHistoInfo("jNeuMultAK8",              hist.Bin("val", "nNeu(J)",                                        120,    0.0,    120.0 ),     1,     1,       'fjw' ),
        'jTau1AK8':                SVJHistoInfo("jTau1AK8",                 hist.Bin("val", r"$\tau_{1}(J)$",                                 40,     0.0,    0.8   ),     1,     1,       'fjw' ),
        'jTau2AK8':                SVJHistoInfo("jTau2AK8",                 hist.Bin("val", r"$\tau_{2}(J)$",                                 40,     0.0,    0.65  ),     1,     1,       'fjw' ),
        'jTau3AK8':                SVJHistoInfo("jTau3AK8",                 hist.Bin("val", r"$\tau_{3}(J)$",                                 40,     0.0,    0.35  ),     1,     1,       'fjw' ),
        'jTau21AK8':               SVJHistoInfo("jTau21AK8",                hist.Bin("val", r"$\tau_{21}(J)$",                                40,     0.0,    1.3   ),     1,     2,       'fjw' ),
        'jTau32AK8':               SVJHistoInfo("jTau32AK8",                hist.Bin("val", r"$\tau_{32}(J)$",                                40,     0.0,    1.3   ),     1,     2,       'fjw' ),
        'jNumBhadronsAK8':         SVJHistoInfo("jNumBhadronsAK8",          hist.Bin("val", "nBHad(J)",                                       30,     0.0,    30.0  ),     1,     1,       'fjw' ),
        'jNumChadronsAK8':         SVJHistoInfo("jNumChadronsAK8",          hist.Bin("val", "nCHad(J)",                                       30,     0.0,    30.0  ),     1,     1,       'fjw' ),
        'jPhoEFractAK8':           SVJHistoInfo("jPhoEFractAK8",            hist.Bin("val", "fPho(J)",                                        50,     0.0,    1.0   ),     1,     1,       'fjw' ),
        'jPhoMultAK8':             SVJHistoInfo("jPhoMultAK8",              hist.Bin("val", "nPho(J)",                                        110,    0.0,    110.0 ),     1,     1,       'fjw' ),
        'jPtDAK8':                 SVJHistoInfo("jPtDAK8",                  hist.Bin("val", "ptD",                                            40,     0.0,    1.2   ),     1,     1,       'fjw' ),
        'jSoftDropMassAK8':        SVJHistoInfo("jSoftDropMassAK8",         hist.Bin("val", r"$m_{SD}(J)$",                                   200,    0.0,    900   ),     1,     1,       'fjw' ),
        'dPhijMETAK8':             SVJHistoInfo("dPhijMETAK8",              hist.Bin("val", r"$\Delta\phi(J,MET)$",                           100,    0.0,    4.0   ),     1,     1,       'fjw' ),
        'dEtaj12AK8':              SVJHistoInfo("dEtaj12AK8",               hist.Bin("val", r"$\Delta\eta(J_{1},J_{2})$",                     200,    0.0,    10.0  ),     2,     0,       'evtw'),
        'dRJ12AK8':                SVJHistoInfo("dRJ12AK8",                 hist.Bin("val", r"$\Delta R(J_{1},J_{2})$",                       100,    0.0,    10.0  ),     2,     0,       'evtw'),
        'dPhiMinjMETAK8':          SVJHistoInfo("dPhiMinjMETAK8",           hist.Bin("val", r"$\Delta\phi_{min}(j,MET)$",                     100,    0.0,    4.0   ),     2,     0,       'evtw'),
        'mT':                      SVJHistoInfo("mT",                       hist.Bin("val", r"$m_{T} (GeV)$",                                 500,    0.0,    6000.0),     2,     0,       'evtw'),
        'METrHT_pt30':             SVJHistoInfo("METrHT_pt30",              hist.Bin("val", r"$MET/H_{T}$",                                   100,    0.0,    3.0   ),     2,     0,       'evtw'),
        'METrST_pt30':             SVJHistoInfo("METrST_pt30",              hist.Bin("val", r"$MET/S_{T}",                                    100,    0.0,    1.0   ),     2,     0,       'evtw'),
        'dPhij1rdPhij2AK8':        SVJHistoInfo("dPhij1rdPhij2AK8",         hist.Bin("val", r"$\Delta\phi(J_{1},MET)/\Delta\phi(J_{2},MET)$", 100,    0.0,    100.0 ),     2,     0,       'evtw'),
        'electronsIso':            SVJHistoInfo("electronsIso",             hist.Bin("val", "electrons iso",                                  100,    0.0,    1.0   ),     0,     1,         'ew'),
        'muonsIso':                SVJHistoInfo("muonsIso",                 hist.Bin("val", "muons iso",                                      100,    0.0,    1.0   ),     0,     1,         'mw'),
        'nonIsoMuonsIso':          SVJHistoInfo("nonIsoMuonsIso",           hist.Bin("val", "NonIsoMuons iso",                                200,    0.0,    10.0  ),     0,     1,       'nimw'),
        'nonIsoMuonsPt':           SVJHistoInfo("nonIsoMuonsPt",            hist.Bin("val", "NonIsoMuons $p_{T}$ [GeV]",                      500,    0.0,    2500.0),     0,     1,       'nimw'),
        # 'mT2_f4_msm':            SVJHistoInfo("mT2_f4_msm",               hist.Bin("val", r"$m_{T2} (GeV)$",                                500,    0.0,    5000.0),     2,     0,       'evtw'),
        # 'mT2_f4_msm_dEta':       SVJHistoInfo("mT2_f4_msm_dEta",          hist.Bin("val", r"$m_{T2} (GeV)$",                                500,    0.0,    5000.0),     2,     0,       'evtw'),
        # 'mT2_f4_msm_dPhi':       SVJHistoInfo("mT2_f4_msm_dPhi",          hist.Bin("val", r"$m_{T2} (GeV)$",                                500,    0.0,    5000.0),     2,     0,       'evtw'),
        # 'mT2_f4_msm_dR':         SVJHistoInfo("mT2_f4_msm_dR",            hist.Bin("val", r"$m_{T2} (GeV)$",                                500,    0.0,    5000.0),     2,     0,       'evtw'),
        # "GenJetsAK8_hvCategory": SVJHistoInfo("GenJetsAK8_hvCategory",    hist.Bin("val", "GenJetAK8 hvCategory",                           32,     0.0,    32.0  ),     0,     1,        'fjw'),
        # "JetsAK8_hvCategory":    SVJHistoInfo("JetsAK8_hvCategory",       hist.Bin("val", "JetAK8 hvCategory",                              32,     0.0,    32.0  ),     1,     2,        'fjw'),
        # "GenMT2_AK8":            SVJHistoInfo("GenMT2_AK8",               hist.Bin("val", r"$m_{T2} (GeV)$",                                500,    0.0,    5000.0),     0,     0,,      'evtw'),
        # "GenJetsAK8_darkPtFrac": SVJHistoInfo("GenJetsAK8_darkPtFrac",    hist.Bin("val", "GenJetAK8 Dark pT Fraction",                     100,    0.0,    1.0   ),     0,     1,        'fjw'),
        'nsvjJetsAK8':             SVJHistoInfo("nsvjJetsAK8",              hist.Bin("val", "Number of SVJ AK8Jets",                          20,     0.0,    20.0  ),     0,     0, 'evtw'      ),
        'nsvjJetsAK8Plus1':        SVJHistoInfo("nsvjJetsAK8Plus1",         hist.Bin("val", "Number of SVJ AK8Jets pred1Jets",                20,     0.0,    20.0  ),     0,     0, 'pred1_evtw'),
        'nsvjJetsAK8Plus2':        SVJHistoInfo("nsvjJetsAK8Plus2",         hist.Bin("val", "Number of SVJ AK8Jets pred2Jets",                20,     0.0,    20.0  ),     0,     0, 'pred2_evtw'),
        'nsvjJetsAK8Plus3':        SVJHistoInfo("nsvjJetsAK8Plus3",         hist.Bin("val", "Number of SVJ AK8Jets pred3Jets",                20,     0.0,    20.0  ),     0,     0, 'pred3_evtw'),
        'nsvjJetsAK8Plus4':        SVJHistoInfo("nsvjJetsAK8Plus4",         hist.Bin("val", "Number of SVJ AK8Jets pred4Jets",                20,     0.0,    20.0  ),     0,     0, 'pred4_evtw'),
        'nnOutput':                SVJHistoInfo("nnOutput",                 hist.Bin("val", "nnOutput",                                       100,    0.0,    1.0   ),     0,     1, 'fjw'       ),
        'svjPtAK8':                SVJHistoInfo("svjPtAK8",                 hist.Bin("val", r"$p_{T}$ [GeV]",                                 280,    0.0,    2800.0),     1,     1, 'svfjw'     ),
        'svjEtaAK8':               SVJHistoInfo("svjEtaAK8",                hist.Bin("val", r"$\eta$",                                        200,   -6.0,    6.0   ),     1,     1, 'svfjw'     ),
    }
    if jNVar:
        # preparing histograms for jN variables
        maxN = 4
        for i in range(maxN):
            jNList = {
                'j{}Pt'.format(i+1):              SVJHistoInfo('j{}Pt'.format(i+1),               hist.Bin("val", r"p_{T}(j_"+str(i+1)+") [GeV]",           200,    0.0,    2500.0),       0,     0,       'evtw'),
                'j{}Eta'.format(i+1):             SVJHistoInfo('j{}Eta'.format(i+1),              hist.Bin("val", r"$\eta$(j_"+str(i+1)+")",                200,   -6.0,    6.0   ),       0,     0,       'evtw'),
                'j{}Phi'.format(i+1):             SVJHistoInfo('j{}Phi'.format(i+1),              hist.Bin("val", r"$\phi$(j_"+str(i+1)+")",                200,   -4.0,    4.0   ),       0,     0,       'evtw'),
                'j{}Axismajor'.format(i+1):       SVJHistoInfo('j{}Axismajor'.format(i+1),        hist.Bin("val", r"$\sigma_{major}(j_"+str(i+1)+")$",      40,     0.0,    0.5   ),       0,     0,       'evtw'),
                'j{}Axisminor'.format(i+1):       SVJHistoInfo('j{}Axisminor'.format(i+1),        hist.Bin("val", r"$\sigma_{minor}(j_"+str(i+1)+")$",      40,     0.0,    0.3   ),       0,     0,       'evtw'),
                'j{}PtD'.format(i+1):             SVJHistoInfo('j{}PtD'.format(i+1),              hist.Bin("val", "ptD(j_"+str(i+1)+")",                    40,     0.0,    1.2   ),       0,     0,       'evtw'),
                'dPhij{}MET'.format(i+1):         SVJHistoInfo('dPhij{}MET'.format(i+1),          hist.Bin("val", r"$\Delta\phi(j_{"+str(i+1)+"},MET)$",    100,    0.0,    4.0   ),       0,     0,       'evtw'),
                'j{}PtAK8'.format(i+1):           SVJHistoInfo('j{}PtAK8'.format(i+1),            hist.Bin("val", r"p_{T}(J_"+str(i+1)+") [GeV]",           200,    0.0,    2500.0),       2,     0,       'evtw'),
                'j{}EtaAK8'.format(i+1):          SVJHistoInfo('j{}EtaAK8'.format(i+1),           hist.Bin("val", r"$\eta$(J_"+str(i+1)+")",                200,   -6.0,    6.0   ),       2,     0,       'evtw'),
                'j{}PhiAK8'.format(i+1):          SVJHistoInfo('j{}PhiAK8'.format(i+1),           hist.Bin("val", r"$\phi$(J_"+str(i+1)+")",                200,   -4.0,    4.0   ),       2,     0,       'evtw'),
                'j{}AxismajorAK8'.format(i+1):    SVJHistoInfo('j{}AxismajorAK8'.format(i+1),     hist.Bin("val", r"$\sigma_{major}(J_"+str(i+1)+")$",      40,     0.0,    0.6   ),       2,     0,       'evtw'),
                'j{}AxisminorAK8'.format(i+1):    SVJHistoInfo('j{}AxisminorAK8'.format(i+1),     hist.Bin("val", r"$\sigma_{minor}(J_"+str(i+1)+")$",      40,     0.0,    0.4   ),       2,     0,       'evtw'),
                'j{}GirthAK8'.format(i+1):        SVJHistoInfo('j{}GirthAK8'.format(i+1),         hist.Bin("val", "girth(J_"+str(i+1)+")",                  40,     0.0,    0.7   ),       2,     0,       'evtw'),
                'j{}PtDAK8'.format(i+1):          SVJHistoInfo('j{}PtDAK8'.format(i+1),           hist.Bin("val", "ptD(J_"+str(i+1)+")",                    40,     0.0,    1.2   ),       2,     0,       'evtw'),
                'j{}Tau1AK8'.format(i+1):         SVJHistoInfo('j{}Tau1AK8'.format(i+1),          hist.Bin("val", r"$\tau_{1}(J_"+str(i+1)+")$",            40,     0.0,    0.8   ),       2,     0,       'evtw'),
                'j{}Tau2AK8'.format(i+1):         SVJHistoInfo('j{}Tau2AK8'.format(i+1),          hist.Bin("val", r"$\tau_{2}(J_"+str(i+1)+")$",            40,     0.0,    0.65  ),       2,     0,       'evtw'),
                'j{}Tau3AK8'.format(i+1):         SVJHistoInfo('j{}Tau3AK8'.format(i+1),          hist.Bin("val", r"$\tau_{3}(J_"+str(i+1)+")$",            40,     0.0,    0.35  ),       2,     0,       'evtw'),
                'j{}Tau21AK8'.format(i+1):        SVJHistoInfo('j{}Tau21AK8'.format(i+1),         hist.Bin("val", r"$\tau_{21}(J_"+str(i+1)+")$",           40,     0.0,    1.3   ),       2,     0,       'evtw'),
                'j{}Tau32AK8'.format(i+1):        SVJHistoInfo('j{}Tau32AK8'.format(i+1),         hist.Bin("val", r"$\tau_{32}(J_"+str(i+1)+")$",           40,     0.0,    1.3   ),       2,     0,       'evtw'),
                'j{}SoftDropMassAK8'.format(i+1): SVJHistoInfo('j{}SoftDropMassAK8'.format(i+1),  hist.Bin("val", r"$m_{SD}(J_"+str(i+1)+")$",              200,    0.0,    900   ),       2,     0,       'evtw'),
                'dPhij{}METAK8'.format(i+1):      SVJHistoInfo('dPhij{}METAK8'.format(i+1),       hist.Bin("val", r"$\Delta\phi(J_{"+str(i+1)+"},MET)$",    100,    0.0,    4.0   ),       2,     0,       'evtw'),
            }
            allVars.update(jNList)

        allComs = list(combinations(range(maxN),2))
        for com in allComs:
            j1 = com[0]+1
            j2 = com[1]+1
            jNjMList = {
                'dEtaj{}{}'.format(j1,j2):         SVJHistoInfo('dEtaj{}{}'.format(j1,j2),         hist.Bin("val",r"$\Delta\eta(j_{"+str(j1)+"},j_{"+str(j2)+"})$",                    200,0.0,5.3  ), 0, 0, 'evtw'),
                'dPhij{}{}'.format(j1,j2):         SVJHistoInfo('dPhij{}{}'.format(j1,j2),         hist.Bin("val",r"$\Delta\phi(j_{"+str(j1)+"},j_{"+str(j2)+"})$",                    100,0.0,4.0  ), 2, 0, 'evtw'),
                'dRj{}{}'.format(j1,j2):           SVJHistoInfo('dRj{}{}'.format(j1,j2),           hist.Bin("val",r"$\Delta R(j_{"+str(j1)+"},j_{"+str(j2)+"})$",                      100,0.0,6.0  ), 0, 0, 'evtw'),
                'dPhij{}rdPhij{}'.format(j1,j2):   SVJHistoInfo('dPhij{}rdPhij{}'.format(j1,j2),   hist.Bin("val",r"$\Delta\phi(j_{"+str(j1)+"},MET)/\Delta\phi(j_{"+str(j2)+"},MET)$",100,0.0,120.0), 0, 0, 'evtw'),
                'dEtaj{}{}AK8'.format(j1,j2):      SVJHistoInfo('dEtaj{}{}AK8'.format(j1,j2),      hist.Bin("val",r"$\Delta\eta(J_{"+str(j1)+"},J_{"+str(j2)+"})$",                    200,0.0,5.3  ), 2, 0, 'evtw'),
                'dPhij{}{}AK8'.format(j1,j2):      SVJHistoInfo('dPhij{}{}AK8'.format(j1,j2),      hist.Bin("val",r"$\Delta\phi(J_{"+str(j1)+"},J_{"+str(j2)+"})$",                    100,0.0,4.0  ), 2, 0, 'evtw'),
                'dRj{}{}AK8'.format(j1,j2):        SVJHistoInfo('dRj{}{}AK8'.format(j1,j2),        hist.Bin("val",r"$\Delta R(J_{"+str(j1)+"},J_{"+str(j2)+"})$",                      100,0.0,6.0  ), 2, 0, 'evtw'),
                'dPhij{}rdPhij{}AK8'.format(j1,j2):SVJHistoInfo('dPhij{}rdPhij{}AK8'.format(j1,j2),hist.Bin("val",r"$\Delta\phi(J_{"+str(j1)+"},MET)/\Delta\phi(J_{"+str(j2)+"},MET)$",100,0.0,120.0), 2, 0, 'evtw'),
            }
            allVars.update(jNjMList)

    #varXName, xbins, npzInfo, flattenInfo, weightName, varYName=None, ybins=None
    histos2D = {
        'jPtAK8vsjEtaAK8': SVJHistoInfo(varXName='jPtAK8', varYName='jEtaAK8',npzInfo=0, flattenInfo=1, weightName='fjw', 
                                        xbins=hist.Bin("val1", "Jet $p_{T}$", 280,  0.0, 2800.0), ybins=hist.Bin("val2", "Jet $\eta$",  200, -6.0,    6.0)),
    }
    allVars.update(histos2D)

    return allVars

