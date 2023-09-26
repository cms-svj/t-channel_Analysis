from itertools import combinations
import hist as h

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
        'eCounter':                SVJHistoInfo("eCounter",                 h.axis.Regular(name="x", label="h_eCounter",                                     bins= 2,    start=-1.1,    stop=1.1   ),     0,     0,       'w1'  ),
        'evtw':                    SVJHistoInfo("evtw",                     h.axis.Regular(name="x", label="h_evtw",                                         bins= 2,    start=-1.1,    stop=1.1   ),     2,     0,       'evtw'),
        'jw':                      SVJHistoInfo("jw",                       h.axis.Regular(name="x", label="h_jw",                                           bins= 2,    start=-1.1,    stop=1.1   ),     0,     2,       'jw'  ),
        'fjw':                     SVJHistoInfo("fjw",                      h.axis.Regular(name="x", label="h_fjw",                                          bins= 2,    start=-1.1,    stop=1.1   ),     1,     2,       'fjw' ),
        'njets':                   SVJHistoInfo("njets",                    h.axis.Regular(name="x", label="Number of Jets",                                 bins=20,    start= 0.0,    stop=20.0  ),     2,     0,       'evtw'),
        'njetsAK8':                SVJHistoInfo("njetsAK8",                 h.axis.Regular(name="x", label="Number of AK8Jets",                              bins=20,    start= 0.0,    stop=20.0  ),     2,     0,       'evtw'),
        'nTruthSVJ':               SVJHistoInfo("nTruthSVJ",                h.axis.Regular(name="x", label="Number of SVJ (Truth)",                          bins=20,    start= 0.0,    stop=20.0  ),     2,     0,       'evtw'),
        'nb':                      SVJHistoInfo("nb",                       h.axis.Regular(name="x", label="Number of b",                                    bins=10,    start= 0.0,    stop=10.0  ),     2,     0,       'evtw'),
        'nl':                      SVJHistoInfo("nl",                       h.axis.Regular(name="x", label="Number of Leptons",                              bins=10,    start= 0.0,    stop=10.0  ),     2,     0,       'evtw'),
        'nnim':                    SVJHistoInfo("nnim",                     h.axis.Regular(name="x", label="Number of NonIsoMuons",                          bins=10,    start= 0.0,    stop=10.0  ),     2,     0,       'evtw'),
        'ht':                      SVJHistoInfo("ht",                       h.axis.Regular(name="x", label=r"$H_{T}$ (GeV)",                                 bins=500,   start= 0.0,    stop=5000.0),     2,     0,       'evtw'),
        'st':                      SVJHistoInfo("st",                       h.axis.Regular(name="x", label=r"$S_{T}$ (GeV)",                                 bins=500,   start= 0.0,    stop=5000.0),     2,     0,       'evtw'),
        'met':                     SVJHistoInfo("met",                      h.axis.Regular(name="x", label="MET [GeV]",                                      bins=500,   start= 0.0,    stop=2000.0),     2,     0,       'evtw'),
        'logMET':                  SVJHistoInfo("logMET",                   h.axis.Regular(name="x", label="log MET",                                        bins=100,   start= 0.0,    stop=10.0),       2,     0,       'evtw'),
        'metSig':                  SVJHistoInfo("metSig",                   h.axis.Regular(name="x", label="MET Significance",                               bins=200,   start= 0.0,    stop=200.0),      2,     0,       'evtw'),
        'logMETSig':               SVJHistoInfo("logMETSig",                h.axis.Regular(name="x", label="log MET Significance",                           bins=100,   start= 0.0,    stop=10.0),       2,     0,       'evtw'),
        'metPhi':                  SVJHistoInfo("metPhi",                   h.axis.Regular(name="x", label=r"MET $\phi$ [GeV]",                              bins=40,    start=-4.0,    stop=4.0   ),     2,     0,       'evtw'),
        # 'madHT':                 SVJHistoInfo("madHT",                    h.axis.Regular(name="x", label=r"$H_{T}$ (GeV)",                                 bins=500,   start= 0.0,    stop=5000.0),     2,     0,       'evtw'),
        'jPt':                     SVJHistoInfo("jPt",                      h.axis.Regular(name="x", label=r"$p_{T}$(j) [GeV]",                              bins=200,   start= 0.0,    stop=2500.0),     0,     1,       'jw'  ),
        'jEta':                    SVJHistoInfo("jEta",                     h.axis.Regular(name="x", label=r"$\eta$(j)",                                     bins=200,   start=-6.0,    stop=6.0   ),     0,     1,       'jw'  ),
        'jPhi':                    SVJHistoInfo("jPhi",                     h.axis.Regular(name="x", label=r"$\phi$(j)",                                     bins=200,   start=-4.0,    stop=4.0   ),     0,     1,       'jw'  ),
        'jE':                      SVJHistoInfo("jE",                       h.axis.Regular(name="x", label=r"$E$(j) [GeV]",                                  bins=200,   start= 0.0,    stop=4000.0),     0,     1,       'jw'  ),
        'jAxismajor':              SVJHistoInfo("jAxismajor",               h.axis.Regular(name="x", label=r"$\sigma_{major}(j)$",                           bins=40,    start= 0.0,    stop=0.5   ),     0,     1,       'jw'  ),
        'jAxisminor':              SVJHistoInfo("jAxisminor",               h.axis.Regular(name="x", label=r"$\sigma_{minor}(j)$",                           bins=40,    start= 0.0,    stop=0.3   ),     0,     1,       'jw'  ),
        'jPtD':                    SVJHistoInfo("jPtD",                     h.axis.Regular(name="x", label="ptD",                                            bins=40,    start= 0.0,    stop=1.2   ),     0,     1,       'jw'  ),
        'dPhiMinjMET':             SVJHistoInfo("dPhiMinjMET",              h.axis.Regular(name="x", label=r"$\Delta\phi_{min}(j,MET)$",                     bins=100,   start= 0.0,    stop=4.0   ),     0,     0,       'evtw'),
        'jPtAK8':                  SVJHistoInfo("jPtAK8",                   h.axis.Regular(name="x", label=r"$p_{T}$(J) [GeV]",                              bins=280,   start= 0.0,    stop=2800.0),     1,     1,       'fjw' ),
        'jEtaAK8':                 SVJHistoInfo("jEtaAK8",                  h.axis.Regular(name="x", label=r"$\eta$(J)",                                     bins=200,   start=-6.0,    stop=6.0   ),     1,     1,       'fjw' ),
        'jPhiAK8':                 SVJHistoInfo("jPhiAK8",                  h.axis.Regular(name="x", label=r"$\phi$(J)",                                     bins=200,   start=-4.0,    stop=4.0   ),     1,     1,       'fjw' ),
        'jEAK8':                   SVJHistoInfo("jEAK8",                    h.axis.Regular(name="x", label=r"$E$(J) [GeV]",                                  bins=300,   start= 0.0,    stop=4500.0),     1,     1,       'fjw' ),
        'jAxismajorAK8':           SVJHistoInfo("jAxismajorAK8",            h.axis.Regular(name="x", label=r"$\sigma_{major}(J)$",                           bins=40,    start= 0.0,    stop=0.6   ),     1,     1,       'fjw' ),
        'jAxisminorAK8':           SVJHistoInfo("jAxisminorAK8",            h.axis.Regular(name="x", label=r"$\sigma_{minor}(J)$",                           bins=40,    start= 0.0,    stop=0.4   ),     1,     1,       'fjw' ),
        'jChEMEFractAK8':          SVJHistoInfo("jChEMEFractAK8",           h.axis.Regular(name="x", label="fChEM(J)",                                       bins=50,    start= 0.0,    stop=1.0   ),     1,     1,       'fjw' ),
        'jChHadEFractAK8':         SVJHistoInfo("jChHadEFractAK8",          h.axis.Regular(name="x", label="fChHad(J)",                                      bins=50,    start= 0.0,    stop=1.0   ),     1,     1,       'fjw' ),
        'jChHadMultAK8':           SVJHistoInfo("jChHadMultAK8",            h.axis.Regular(name="x", label="nChHad(J)",                                      bins=145,   start= 0.0,    stop=145.0 ),     1,     1,       'fjw' ),
        'jChMultAK8':              SVJHistoInfo("jChMultAK8",               h.axis.Regular(name="x", label="nCh(J)",                                         bins=145,   start= 0.0,    stop=145.0 ),     1,     1,       'fjw' ),
        'jecfN2b1AK8':             SVJHistoInfo("jecfN2b1AK8",              h.axis.Regular(name="x", label="ecfN2b1(J)",                                     bins=50,    start= 0.0,    stop=0.6   ),     1,     1,       'fjw' ),
        'jecfN2b2AK8':             SVJHistoInfo("jecfN2b2AK8",              h.axis.Regular(name="x", label="ecfN2b2(J)",                                     bins=50,    start= 0.0,    stop=0.4   ),     1,     1,       'fjw' ),
        'jecfN3b1AK8':             SVJHistoInfo("jecfN3b1AK8",              h.axis.Regular(name="x", label="ecfN3b1(J)",                                     bins=50,    start= 0.0,    stop=6.0   ),     1,     1,       'fjw' ),
        'jecfN3b2AK8':             SVJHistoInfo("jecfN3b2AK8",              h.axis.Regular(name="x", label="ecfN3b2(J)",                                     bins=50,    start= 0.0,    stop=5.0   ),     1,     1,       'fjw' ),
        'jEleEFractAK8':           SVJHistoInfo("jEleEFractAK8",            h.axis.Regular(name="x", label="fEle(J)",                                        bins=50,    start= 0.0,    stop=1.0   ),     1,     1,       'fjw' ),
        'jEleMultAK8':             SVJHistoInfo("jEleMultAK8",              h.axis.Regular(name="x", label="nEle(J)",                                        bins=8,     start= 0.0,    stop=8.0   ),     1,     1,       'fjw' ),
        'jGirthAK8':               SVJHistoInfo("jGirthAK8",                h.axis.Regular(name="x", label="girth(J)",                                       bins=40,    start= 0.0,    stop=0.7   ),     1,     1,       'fjw' ),
        'jHfEMEFractAK8':          SVJHistoInfo("jHfEMEFractAK8",           h.axis.Regular(name="x", label="fHFEM(J)",                                       bins=50,    start= 0.0,    stop=1.0   ),     1,     1,       'fjw' ),
        'jHfHadEFractAK8':         SVJHistoInfo("jHfHadEFractAK8",          h.axis.Regular(name="x", label="fHFHad(J)",                                      bins=50,    start= 0.0,    stop=1.0   ),     1,     1,       'fjw' ),
        'jMultAK8':                SVJHistoInfo("jMultAK8",                 h.axis.Regular(name="x", label="mult(J)",                                        bins=250,   start= 0.0,    stop=250.0 ),     1,     1,       'fjw' ),
        'jMuEFractAK8':            SVJHistoInfo("jMuEFractAK8",             h.axis.Regular(name="x", label="fMu(J)",                                         bins=50,    start= 0.0,    stop=1.0   ),     1,     1,       'fjw' ),
        'jMuMultAK8':              SVJHistoInfo("jMuMultAK8",               h.axis.Regular(name="x", label="nMu(J)",                                         bins=8,     start= 0.0,    stop=10.0  ),     1,     1,       'fjw' ),
        'jNeuEmEFractAK8':         SVJHistoInfo("jNeuEmEFractAK8",          h.axis.Regular(name="x", label="fNeuEM(J)",                                      bins=50,    start= 0.0,    stop=1.0   ),     1,     1,       'fjw' ),
        'jNeuHadEFractAK8':        SVJHistoInfo("jNeuHadEFractAK8",         h.axis.Regular(name="x", label="fNeuHad(J)",                                     bins=50,    start= 0.0,    stop=1.0   ),     1,     1,       'fjw' ),
        'jNeuHadMultAK8':          SVJHistoInfo("jNeuHadMultAK8",           h.axis.Regular(name="x", label="nNeuHad(J)",                                     bins=25,    start= 0.0,    stop=25.0  ),     1,     1,       'fjw' ),
        'jNeuMultAK8':             SVJHistoInfo("jNeuMultAK8",              h.axis.Regular(name="x", label="nNeu(J)",                                        bins=120,   start= 0.0,    stop=120.0 ),     1,     1,       'fjw' ),
        'jTau1AK8':                SVJHistoInfo("jTau1AK8",                 h.axis.Regular(name="x", label=r"$\tau_{1}(J)$",                                 bins=40,    start= 0.0,    stop=0.8   ),     1,     1,       'fjw' ),
        'jTau2AK8':                SVJHistoInfo("jTau2AK8",                 h.axis.Regular(name="x", label=r"$\tau_{2}(J)$",                                 bins=40,    start= 0.0,    stop=0.65  ),     1,     1,       'fjw' ),
        'jTau3AK8':                SVJHistoInfo("jTau3AK8",                 h.axis.Regular(name="x", label=r"$\tau_{3}(J)$",                                 bins=40,    start= 0.0,    stop=0.35  ),     1,     1,       'fjw' ),
        'jTau21AK8':               SVJHistoInfo("jTau21AK8",                h.axis.Regular(name="x", label=r"$\tau_{21}(J)$",                                bins=40,    start= 0.0,    stop=1.3   ),     1,     2,       'fjw' ),
        'jTau32AK8':               SVJHistoInfo("jTau32AK8",                h.axis.Regular(name="x", label=r"$\tau_{32}(J)$",                                bins=40,    start= 0.0,    stop=1.3   ),     1,     2,       'fjw' ),
        'jPhoEFractAK8':           SVJHistoInfo("jPhoEFractAK8",            h.axis.Regular(name="x", label="fPho(J)",                                        bins=50,    start= 0.0,    stop=1.0   ),     1,     1,       'fjw' ),
        'jPhoMultAK8':             SVJHistoInfo("jPhoMultAK8",              h.axis.Regular(name="x", label="nPho(J)",                                        bins=110,   start= 0.0,    stop=110.0 ),     1,     1,       'fjw' ),
        'jPtDAK8':                 SVJHistoInfo("jPtDAK8",                  h.axis.Regular(name="x", label="ptD",                                            bins=40,    start= 0.0,    stop=1.2   ),     1,     1,       'fjw' ),
        'jSoftDropMassAK8':        SVJHistoInfo("jSoftDropMassAK8",         h.axis.Regular(name="x", label=r"$m_{SD}(J)$",                                   bins=200,   start= 0.0,    stop=900   ),     1,     1,       'fjw' ),
        'dPhijMETAK8':             SVJHistoInfo("dPhijMETAK8",              h.axis.Regular(name="x", label=r"$\Delta\phi(J,MET)$",                           bins=100,   start= 0.0,    stop=4.0   ),     1,     1,       'fjw' ),
        'dEtaj12AK8':              SVJHistoInfo("dEtaj12AK8",               h.axis.Regular(name="x", label=r"$\Delta\eta(J_{1},J_{2})$",                     bins=200,   start= 0.0,    stop=10.0  ),     2,     0,       'evtw'),
        'dRJ12AK8':                SVJHistoInfo("dRJ12AK8",                 h.axis.Regular(name="x", label=r"$\Delta R(J_{1},J_{2})$",                       bins=100,   start= 0.0,    stop=10.0  ),     2,     0,       'evtw'),
        'dPhiMinjMETAK8':          SVJHistoInfo("dPhiMinjMETAK8",           h.axis.Regular(name="x", label=r"$\Delta\phi_{min}(j,MET)$",                     bins=100,   start= 0.0,    stop=4.0   ),     2,     0,       'evtw'),
        'mT':                      SVJHistoInfo("mT",                       h.axis.Regular(name="x", label=r"$m_{T} (GeV)$",                                 bins=500,   start= 0.0,    stop=6000.0),     2,     0,       'evtw'),
        'METrHT_pt30':             SVJHistoInfo("METrHT_pt30",              h.axis.Regular(name="x", label=r"$MET/H_{T}$",                                   bins=100,   start= 0.0,    stop=3.0   ),     2,     0,       'evtw'),
        'METrST_pt30':             SVJHistoInfo("METrST_pt30",              h.axis.Regular(name="x", label=r"$MET/S_{T}",                                    bins=100,   start= 0.0,    stop=1.0   ),     2,     0,       'evtw'),
        'dPhij1rdPhij2AK8':        SVJHistoInfo("dPhij1rdPhij2AK8",         h.axis.Regular(name="x", label=r"$\Delta\phi(J_{1},MET)/\Delta\phi(J_{2},MET)$", bins=100,   start= 0.0,    stop=100.0 ),     2,     0,       'evtw'),
        'electronsIso':            SVJHistoInfo("electronsIso",             h.axis.Regular(name="x", label="electrons iso",                                  bins=100,   start= 0.0,    stop=1.0   ),     0,     1,         'ew'),
        #'muonsIso':                SVJHistoInfo("muonsIso",                 h.axis.Regular(name="x", label="muons iso",                                      bins=100,   start= 0.0,    stop=1.0   ),     0,     1,         'mw'),
        #'nonIsoMuonsIso':          SVJHistoInfo("nonIsoMuonsIso",           h.axis.Regular(name="x", label="NonIsoMuons iso",                                bins=200,   start= 0.0,    stop=10.0  ),     0,     1,       'nimw'),
        #'nonIsoMuonsPt':           SVJHistoInfo("nonIsoMuonsPt",            h.axis.Regular(name="x", label="NonIsoMuons $p_{T}$ [GeV]",                      bins=500,   start= 0.0,    stop=2500.0),     0,     1,       'nimw'),
        # 'mT2_f4_msm':            SVJHistoInfo("mT2_f4_msm",               h.axis.Regular(name="x", label=r"$m_{T2} (GeV)$",                                bins=500,   start= 0.0,    stop=5000.0),     2,     0,       'evtw'),
        # 'mT2_f4_msm_dEta':       SVJHistoInfo("mT2_f4_msm_dEta",          h.axis.Regular(name="x", label=r"$m_{T2} (GeV)$",                                bins=500,   start= 0.0,    stop=5000.0),     2,     0,       'evtw'),
        # 'mT2_f4_msm_dPhi':       SVJHistoInfo("mT2_f4_msm_dPhi",          h.axis.Regular(name="x", label=r"$m_{T2} (GeV)$",                                bins=500,   start= 0.0,    stop=5000.0),     2,     0,       'evtw'),
        # 'mT2_f4_msm_dR':         SVJHistoInfo("mT2_f4_msm_dR",            h.axis.Regular(name="x", label=r"$m_{T2} (GeV)$",                                bins=500,   start= 0.0,    stop=5000.0),     2,     0,       'evtw'),
        #"GenJetsAK8_hvCategory":   SVJHistoInfo("GenJetsAK8_hvCategory",    h.axis.Regular(name="x", label="GenJetAK8 hvCategory",                           bins=32,    start= 0.0,    stop=32.0  ),     0,     1,        'gfjw'),
        "nNMedEvent":              SVJHistoInfo("nNMedEvent",               h.axis.Regular(name="x", label="N-Med Event",                                    bins=4,    start= 0.0,     stop=4.0  ),      0,     0,        'evtw'),
        "JetsAK8_hvCategory":      SVJHistoInfo("JetsAK8_hvCategory",       h.axis.Regular(name="x", label="JetAK8 hvCategory",                              bins=32,    start= 0.0,    stop=32.0  ),     1,     2,        'fjw'),
        # "GenMT2_AK8":            SVJHistoInfo("GenMT2_AK8",               h.axis.Regular(name="x", label=r"$m_{T2} (GeV)$",                                bins=500,   start= 0.0,    stop=5000.0),     0,     0,,      'evtw'),
        "JetsAK8_darkPtFrac":      SVJHistoInfo("JetsAK8_darkPtFrac",       h.axis.Regular(name="x", label="JetAK8 Dark pT Fraction",                        bins=100,   start= 0.0,    stop=1.0   ),     1,     2,        'fjw'),
        'nsvjJetsAK8':             SVJHistoInfo("nsvjJetsAK8",              h.axis.Regular(name="x", label="Number of SVJ AK8Jets",                          bins=20,    start= 0.0,    stop=20.0  ),     0,     0, 'evtw'      ),
        'nsvjJetsAK8Plus1':        SVJHistoInfo("nsvjJetsAK8Plus1",         h.axis.Regular(name="x", label="Number of SVJ AK8Jets pred1Jets",                bins=20,    start= 0.0,    stop=20.0  ),     0,     0, 'pred1_evtw'),
        'nsvjJetsAK8Plus2':        SVJHistoInfo("nsvjJetsAK8Plus2",         h.axis.Regular(name="x", label="Number of SVJ AK8Jets pred2Jets",                bins=20,    start= 0.0,    stop=20.0  ),     0,     0, 'pred2_evtw'),
        'nsvjJetsAK8Plus3':        SVJHistoInfo("nsvjJetsAK8Plus3",         h.axis.Regular(name="x", label="Number of SVJ AK8Jets pred3Jets",                bins=20,    start= 0.0,    stop=20.0  ),     0,     0, 'pred3_evtw'),
        'nsvjJetsAK8Plus4':        SVJHistoInfo("nsvjJetsAK8Plus4",         h.axis.Regular(name="x", label="Number of SVJ AK8Jets pred4Jets",                bins=20,    start= 0.0,    stop=20.0  ),     0,     0, 'pred4_evtw'),
        'nnOutput':                SVJHistoInfo("nnOutput",                 h.axis.Regular(name="x", label="nnOutput",                                       bins=100,   start= 0.0,    stop=1.0   ),     1,     1, 'fjw'       ),
        'nnEventOutput':           SVJHistoInfo("nnEventOutput",            h.axis.Regular(name="x", label="nnEventOutput",                                  bins=100,   start= 0.0,    stop=1.0   ),     0,     0, 'evtw'       ),
        'nnEventOutputrMET':       SVJHistoInfo("nnEventOutputrMET",        h.axis.Regular(name="x", label="nnEventOutputrMET",                              bins=100,   start= 0.0,    stop=1.0   ),     0,     0, 'evtw'       ),
        'nnEventOutputrST':        SVJHistoInfo("nnEventOutputrST",         h.axis.Regular(name="x", label="nnEventOutputrST",                               bins=100,   start= 0.0,    stop=0.02  ),     0,     0, 'evtw'       ),
        #'svjPtAK8':                SVJHistoInfo("svjPtAK8",                 h.axis.Regular(name="x", label=r"$p_{T}$ [GeV]",                                 bins=280,   start= 0.0,    stop=2800.0),     0,     1, 'svfjw'     ),
        #'svjEtaAK8':               SVJHistoInfo("svjEtaAK8",                h.axis.Regular(name="x", label=r"$\eta$",                                        bins=200,   start=-6.0,    stop=6.0   ),     0,     1, 'svfjw'     ),
    }
    if jNVar:
        # preparing histograms for jN variables
        maxN = 4
        for i in range(maxN):
            jNList = {
                'j{}Pt'.format(i+1):              SVJHistoInfo('j{}Pt'.format(i+1),               h.axis.Regular(name="x", label=r"p_{T}(j_"+str(i+1)+") [GeV]",           bins=200,   start= 0.0,    stop=2500.0),       0,     0,       'evtw'),
                'j{}Eta'.format(i+1):             SVJHistoInfo('j{}Eta'.format(i+1),              h.axis.Regular(name="x", label=r"$\eta$(j_"+str(i+1)+")",                bins=200,   start=-6.0,    stop=6.0   ),       0,     0,       'evtw'),
                'j{}Phi'.format(i+1):             SVJHistoInfo('j{}Phi'.format(i+1),              h.axis.Regular(name="x", label=r"$\phi$(j_"+str(i+1)+")",                bins=200,   start=-4.0,    stop=4.0   ),       0,     0,       'evtw'),
                'j{}E'.format(i+1):               SVJHistoInfo('j{}E'.format(i+1),                h.axis.Regular(name="x", label=r"E(j_"+str(i+1)+") [GeV]",               bins=200,   start= 0.0,    stop=4000.0),       0,     0,       'evtw'),
                'j{}Axismajor'.format(i+1):       SVJHistoInfo('j{}Axismajor'.format(i+1),        h.axis.Regular(name="x", label=r"$\sigma_{major}(j_"+str(i+1)+")$",      bins=40,    start= 0.0,    stop=0.5   ),       0,     0,       'evtw'),
                'j{}Axisminor'.format(i+1):       SVJHistoInfo('j{}Axisminor'.format(i+1),        h.axis.Regular(name="x", label=r"$\sigma_{minor}(j_"+str(i+1)+")$",      bins=40,    start= 0.0,    stop=0.3   ),       0,     0,       'evtw'),
                'j{}PtD'.format(i+1):             SVJHistoInfo('j{}PtD'.format(i+1),              h.axis.Regular(name="x", label="ptD(j_"+str(i+1)+")",                    bins=40,    start= 0.0,    stop=1.2   ),       0,     0,       'evtw'),
                'dPhij{}MET'.format(i+1):         SVJHistoInfo('dPhij{}MET'.format(i+1),          h.axis.Regular(name="x", label=r"$\Delta\phi(j_{"+str(i+1)+"},MET)$",    bins=100,   start= 0.0,    stop=4.0   ),       0,     0,       'evtw'),
                'j{}PtAK8'.format(i+1):           SVJHistoInfo('j{}PtAK8'.format(i+1),            h.axis.Regular(name="x", label=r"p_{T}(J_"+str(i+1)+") [GeV]",           bins=200,   start= 0.0,    stop=2500.0),       2,     0,       'evtw'),
                'j{}EtaAK8'.format(i+1):          SVJHistoInfo('j{}EtaAK8'.format(i+1),           h.axis.Regular(name="x", label=r"$\eta$(J_"+str(i+1)+")",                bins=200,   start=-6.0,    stop=6.0   ),       2,     0,       'evtw'),
                'j{}PhiAK8'.format(i+1):          SVJHistoInfo('j{}PhiAK8'.format(i+1),           h.axis.Regular(name="x", label=r"$\phi$(J_"+str(i+1)+")",                bins=200,   start=-4.0,    stop=4.0   ),       2,     0,       'evtw'),
                'j{}EAK8'.format(i+1):            SVJHistoInfo('j{}EAK8'.format(i+1),             h.axis.Regular(name="x", label=r"E(J_"+str(i+1)+") [GeV]",               bins=200,   start= 0.0,    stop=4500.0),       2,     0,       'evtw'),
                'j{}AxismajorAK8'.format(i+1):    SVJHistoInfo('j{}AxismajorAK8'.format(i+1),     h.axis.Regular(name="x", label=r"$\sigma_{major}(J_"+str(i+1)+")$",      bins=40,    start= 0.0,    stop=0.6   ),       2,     0,       'evtw'),
                'j{}AxisminorAK8'.format(i+1):    SVJHistoInfo('j{}AxisminorAK8'.format(i+1),     h.axis.Regular(name="x", label=r"$\sigma_{minor}(J_"+str(i+1)+")$",      bins=40,    start= 0.0,    stop=0.4   ),       2,     0,       'evtw'),
                'j{}GirthAK8'.format(i+1):        SVJHistoInfo('j{}GirthAK8'.format(i+1),         h.axis.Regular(name="x", label="girth(J_"+str(i+1)+")",                  bins=40,    start= 0.0,    stop=0.7   ),       2,     0,       'evtw'),
                'j{}PtDAK8'.format(i+1):          SVJHistoInfo('j{}PtDAK8'.format(i+1),           h.axis.Regular(name="x", label="ptD(J_"+str(i+1)+")",                    bins=40,    start= 0.0,    stop=1.2   ),       2,     0,       'evtw'),
                'j{}Tau1AK8'.format(i+1):         SVJHistoInfo('j{}Tau1AK8'.format(i+1),          h.axis.Regular(name="x", label=r"$\tau_{1}(J_"+str(i+1)+")$",            bins=40,    start= 0.0,    stop=0.8   ),       2,     0,       'evtw'),
                'j{}Tau2AK8'.format(i+1):         SVJHistoInfo('j{}Tau2AK8'.format(i+1),          h.axis.Regular(name="x", label=r"$\tau_{2}(J_"+str(i+1)+")$",            bins=40,    start= 0.0,    stop=0.65  ),       2,     0,       'evtw'),
                'j{}Tau3AK8'.format(i+1):         SVJHistoInfo('j{}Tau3AK8'.format(i+1),          h.axis.Regular(name="x", label=r"$\tau_{3}(J_"+str(i+1)+")$",            bins=40,    start= 0.0,    stop=0.35  ),       2,     0,       'evtw'),
                'j{}Tau21AK8'.format(i+1):        SVJHistoInfo('j{}Tau21AK8'.format(i+1),         h.axis.Regular(name="x", label=r"$\tau_{21}(J_"+str(i+1)+")$",           bins=40,    start= 0.0,    stop=1.3   ),       2,     0,       'evtw'),
                'j{}Tau32AK8'.format(i+1):        SVJHistoInfo('j{}Tau32AK8'.format(i+1),         h.axis.Regular(name="x", label=r"$\tau_{32}(J_"+str(i+1)+")$",           bins=40,    start= 0.0,    stop=1.3   ),       2,     0,       'evtw'),
                'j{}SoftDropMassAK8'.format(i+1): SVJHistoInfo('j{}SoftDropMassAK8'.format(i+1),  h.axis.Regular(name="x", label=r"$m_{SD}(J_"+str(i+1)+")$",              bins=200,   start= 0.0,    stop=900   ),       2,     0,       'evtw'),
                'dPhij{}METAK8'.format(i+1):      SVJHistoInfo('dPhij{}METAK8'.format(i+1),       h.axis.Regular(name="x", label=r"$\Delta\phi(J_{"+str(i+1)+"},MET)$",    bins=100,   start= 0.0,    stop=4.0   ),       2,     0,       'evtw'),
                "J{}_hvCategory".format(i+1):     SVJHistoInfo("J{}_hvCategory".format(i+1),      h.axis.Regular(name="x", label="J_"+str(i+1)+" hvCategory",              bins=32,    start= 0.0,    stop=32.0  ),       0,     0,       'evtw'),
                "J{}_darkPtFrac".format(i+1):     SVJHistoInfo("J{}_darkPtFrac".format(i+1),      h.axis.Regular(name="x", label="J_"+str(i+1)+" Dark pT Fraction",        bins=100,   start= 0.0,    stop=1.0   ),       0,     1,       'evtw'),
            }
            allVars.update(jNList)

        allComs = list(combinations(range(maxN),2))
        for com in allComs:
            j1 = com[0]+1
            j2 = com[1]+1
            jNjMList = {
                'dEtaj{}{}'.format(j1,j2):         SVJHistoInfo('dEtaj{}{}'.format(j1,j2),         h.axis.Regular(name="x",label=r"$\Delta\eta(j_{"+str(j1)+"},j_{"+str(j2)+"})$",                    bins=200,start=0.0,stop=5.3  ), 0, 0, 'evtw'),
                'dPhij{}{}'.format(j1,j2):         SVJHistoInfo('dPhij{}{}'.format(j1,j2),         h.axis.Regular(name="x",label=r"$\Delta\phi(j_{"+str(j1)+"},j_{"+str(j2)+"})$",                    bins=100,start=0.0,stop=4.0  ), 2, 0, 'evtw'),
                'dRj{}{}'.format(j1,j2):           SVJHistoInfo('dRj{}{}'.format(j1,j2),           h.axis.Regular(name="x",label=r"$\Delta R(j_{"+str(j1)+"},j_{"+str(j2)+"})$",                      bins=100,start=0.0,stop=6.0  ), 0, 0, 'evtw'),
                'dPhij{}rdPhij{}'.format(j1,j2):   SVJHistoInfo('dPhij{}rdPhij{}'.format(j1,j2),   h.axis.Regular(name="x",label=r"$\Delta\phi(j_{"+str(j1)+"},MET)/\Delta\phi(j_{"+str(j2)+"},MET)$",bins=100,start=0.0,stop=120.0), 0, 0, 'evtw'),
                'dEtaj{}{}AK8'.format(j1,j2):      SVJHistoInfo('dEtaj{}{}AK8'.format(j1,j2),      h.axis.Regular(name="x",label=r"$\Delta\eta(J_{"+str(j1)+"},J_{"+str(j2)+"})$",                    bins=200,start=0.0,stop=5.3  ), 2, 0, 'evtw'),
                'dPhij{}{}AK8'.format(j1,j2):      SVJHistoInfo('dPhij{}{}AK8'.format(j1,j2),      h.axis.Regular(name="x",label=r"$\Delta\phi(J_{"+str(j1)+"},J_{"+str(j2)+"})$",                    bins=100,start=0.0,stop=4.0  ), 2, 0, 'evtw'),
                'dRj{}{}AK8'.format(j1,j2):        SVJHistoInfo('dRj{}{}AK8'.format(j1,j2),        h.axis.Regular(name="x",label=r"$\Delta R(J_{"+str(j1)+"},J_{"+str(j2)+"})$",                      bins=100,start=0.0,stop=6.0  ), 2, 0, 'evtw'),
                'dPhij{}rdPhij{}AK8'.format(j1,j2):SVJHistoInfo('dPhij{}rdPhij{}AK8'.format(j1,j2),h.axis.Regular(name="x",label=r"$\Delta\phi(J_{"+str(j1)+"},MET)/\Delta\phi(J_{"+str(j2)+"},MET)$",bins=100,start=0.0,stop=120.0), 2, 0, 'evtw'),
            }
            allVars.update(jNjMList)

    #varXName, xbins, npzInfo, flattenInfo, weightName, varYName=None, ybins=None
    histos2D = {
        'jPtAK8vsjEtaAK8': SVJHistoInfo(varXName='jPtAK8', varYName='jEtaAK8',npzInfo=0, flattenInfo=1, weightName='fjw', 
                                        xbins=h.axis.Regular(name="x", label="Jet $p_{T}$", bins=280, start=0.0, stop=2800.0), ybins=h.axis.Regular(name="y", label="Jet $\eta$",  bins=200, start=-6.0, stop=6.0)),
        'njetsAK8vsDNN': SVJHistoInfo(varXName='njetsAK8', varYName='nnEventOutput',npzInfo=0, flattenInfo=0, weightName='evtw', 
                                        xbins=h.axis.Regular(name="x", label="Number of AK8Jets",      bins=20, start=0.0, stop=20.0), ybins=h.axis.Regular(name="y", label="nnEventOutput",  bins=100, start=0.0, stop=1.0 )),
        'nTagSVJvsnTruSVJ': SVJHistoInfo(varXName='nsvjJetsAK8', varYName='nTruthSVJ',npzInfo=0, flattenInfo=0, weightName='evtw', 
                                        xbins=h.axis.Regular(name="x", label="Number of Tagged SVJ",      bins=20, start=0.0, stop=20.0), ybins=h.axis.Regular(name="y", label="Number of Truth SVJ",  bins=20, start=0.0, stop=20.0 )),
        'METvsDNN': SVJHistoInfo(varXName='met', varYName='nnEventOutput',npzInfo=0, flattenInfo=0, weightName='evtw', 
                                        xbins=h.axis.Regular(name="x", label="MET [GeV]",      bins=500, start=0.0, stop=2000.0), ybins=h.axis.Regular(name="y", label="nnEventOutput",  bins=100, start=0.0, stop=1.0 )),
        'logMETvsDNN': SVJHistoInfo(varXName='logMET', varYName='nnEventOutput',npzInfo=0, flattenInfo=0, weightName='evtw', 
                                        xbins=h.axis.Regular(name="x", label="log MET",      bins=100, start=0.0, stop=10.0), ybins=h.axis.Regular(name="y", label="nnEventOutput",  bins=100, start=0.0, stop=1.0 )),
        'metSigvsDNN': SVJHistoInfo(varXName='metSig', varYName='nnEventOutput',npzInfo=0, flattenInfo=0, weightName='evtw', 
                                        xbins=h.axis.Regular(name="x", label="MET Significance",      bins=200, start=0.0, stop=200.0), ybins=h.axis.Regular(name="y", label="nnEventOutput",  bins=100, start=0.0, stop=1.0 )),
        'logMETSigvsDNN': SVJHistoInfo(varXName='logMETSig', varYName='nnEventOutput',npzInfo=0, flattenInfo=0, weightName='evtw', 
                                        xbins=h.axis.Regular(name="x", label="log MET Significance",      bins=100, start=0.0, stop=10.0), ybins=h.axis.Regular(name="y", label="nnEventOutput",  bins=100, start=0.0, stop=1.0 )),
        'STvsDNN':  SVJHistoInfo(varXName='st',  varYName='nnEventOutput',npzInfo=0, flattenInfo=0, weightName='evtw', 
                                        xbins=h.axis.Regular(name="x", label=r"$S_{T}$ (GeV)", bins=500, start=0.0, stop=5000.0), ybins=h.axis.Regular(name="y", label="nnEventOutput",  bins=100, start=0.0, stop=1.0 )),
        'METvsDNNrMET': SVJHistoInfo(varXName='met', varYName='nnEventOutputrMET',npzInfo=0, flattenInfo=0, weightName='evtw', 
                                        xbins=h.axis.Regular(name="x", label="MET [GeV]",      bins=500, start=0.0, stop=2000.0), ybins=h.axis.Regular(name="y", label="nnEventOutput",  bins=100, start=0.0, stop=1.0 )),
        'STvsDNNrST':  SVJHistoInfo(varXName='st',  varYName='nnEventOutputrST',npzInfo=0, flattenInfo=0, weightName='evtw', 
                                        xbins=h.axis.Regular(name="x", label=r"$S_{T}$ (GeV)", bins=500, start=0.0, stop=5000.0), ybins=h.axis.Regular(name="y", label="nnEventOutput",  bins=100, start=0.0, stop=0.02 )),
    }
    allVars.update(histos2D)

    return allVars
