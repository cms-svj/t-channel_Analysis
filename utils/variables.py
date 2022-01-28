# [xlabel,number of bins,xmin,xmax,whether to keep it in npz for training (0=do not keep, 1=keep as is, 2=keep but make sure the length is the same as AK8 variables), whether to flatten the array or not when filling histogram (2 = ak.flatten(), 1 = .flatten(), 0 = do not flatten)]
variables = {
'eCounter':               ["h_eCounter",                                        2,    -1.1,    1.1,         0,     0],
'evtw':                   ["h_evtw",                                            2,    -1.1,    1.1,         0,     0],
'jw':                     ["h_jw",                                              2,    -1.1,    1.1,         0,     2],
'fjw':                    ["h_fjw",                                             2,    -1.1,    1.1,         1,     2],
'njets':                  ["Number of Jets",                                   20,     0.0,    20.0,        2,     0],
'njetsAK8':               ["Number of AK8Jets",                                20,     0.0,    20.0,        2,     0],
'nb':                     ["Number of b",                                      10,     0.0,    10.0,        2,     0],
'nl':                     ["Number of Leptons",                                10,     0.0,    10.0,        2,     0],
'ht':                     [r"$H_{T}$ (GeV)",                                   500,    0.0,    5000.0,      2,     0],
'st':                     [r"$S_{T}$ (GeV)",                                   500,    0.0,    5000.0,      2,     0],
'met':                    ["MET [GeV]",                                        500,    0.0,    2000.0,      2,     0],
'metPhi':                 [r"MET $\phi$ [GeV]",                                40,    -4.0,    4.0,         0,     1],
'madHT':                  [r"$H_{T}$ (GeV)",                                   500,    0.0,    5000.0,      2,     0],
'jPt':                    [r"$p_{T}$ [GeV]",                                   200,    0.0,    2000.0,      0,     1],
'jEta':                   [r"$\eta$",                                          200,   -6.0,    6.0,         0,     1],
'jPhi':                   [r"$\phi$",                                          200,   -4.0,    4.0,         0,     1],
'jAxismajor':             [r"$\sigma_{major}(j)$",                             40,     0.0,    0.5,         0,     1],
'jAxisminor':             [r"$\sigma_{minor}(j)$",                             40,     0.0,    0.3,         0,     1],
'jPtD':                   ["ptD",                                              40,     0.0,    1.2,         0,     1],
'dPhiMinjMET':            [r"$\Delta\phi_{min}(j,MET)$",                       100,    0.0,    4.0,         0,     0],
'dEtaj12':                [r"$\Delta\eta(J_{1},J_{2})$",                       200,    0.0,    10.0,        0,     0],
'dRJ12':                  [r"$\Delta R(J_{1},J_{2})$",                         100,    0.0,    10.0,        0,     0],
'jPtAK8':                 [r"$p_{T}$ [GeV]",                                   200,    0.0,    2000.0,      1,     1],
'jEtaAK8':                [r"$\eta$",                                          200,   -6.0,    6.0,         1,     1],
'jPhiAK8':                [r"$\phi$",                                          200,   -4.0,    4.0,         1,     1],
'jAxismajorAK8':          [r"$\sigma_{major}(J)$",                             40,     0.0,    0.5,         1,     1],
'jAxisminorAK8':          [r"$\sigma_{minor}(J)$",                             40,     0.0,    0.3,         1,     1],
'jGirthAK8':              ["girth(J)",                                         40,     0.0,    0.5,         1,     1],
'jPtDAK8':                ["ptD",                                              40,     0.0,    1.2,         1,     1],
'jTau1AK8':               [r"$\tau_{1}(J)$",                                   40,     0.0,    0.8,         1,     1],
'jTau2AK8':               [r"$\tau_{2}(J)$",                                   40,     0.0,    0.65,        1,     1],
'jTau3AK8':               [r"$\tau_{3}(J)$",                                   40,     0.0,    0.35,        1,     1],
'jTau21AK8':              [r"$\tau_{21}(J)$",                                  40,     0.0,    1.3,         1,     2],
'jTau32AK8':              [r"$\tau_{32}(J)$",                                  40,     0.0,    1.3,         1,     2],
'jSoftDropMassAK8':       [r"$m_{SD}(J)$",                                     40,     0.0,    200,         1,     1],
'jecfN2b1AK8':            ["ecfN2b1(J)",                                       50,     0.0,    0.6,         1,     1],
'jecfN2b2AK8':            ["ecfN2b2(J)",                                       50,     0.0,    0.4,         1,     1],
'jecfN3b1AK8':            ["ecfN3b1(J)",                                       50,     0.0,    5.0,         1,     1],
'jecfN3b2AK8':            ["ecfN3b2(J)",                                       50,     0.0,    5.0,         1,     1],
'jEleEFractAK8':          ["fEle(J)",                                          50,     0.0,    1.0,         1,     1],
'jMuEFractAK8':           ["fMu(J)",                                           50,     0.0,    1.0,         1,     1],
'jNeuHadEFractAK8':       ["fNeuHad(J)",                                       50,     0.0,    1.0,         1,     1],
'jPhoEFractAK8':          ["fPho(J)",                                          50,     0.0,    1.0,         1,     1],
'jNeuEmEFractAK8':        ["fNeuEM(J)",                                        50,     0.0,    1.0,         1,     1],
'jHfHadEFractAK8':        ["fHFHad(J)",                                        50,     0.0,    1.0,         1,     1],
'jHfEMEFractAK8':         ["fHFEM(J)",                                         50,     0.0,    1.0,         1,     1],
'jChEMEFractAK8':         ["fChEM(J)",                                         50,     0.0,    1.0,         1,     1],
'jPhoMultAK8':            ["nPho(J)",                                          110,    0.0,    110.0,       1,     1],
'jNeuMultAK8':            ["nNeu(J)",                                          120,    0.0,    120.0,       1,     1],
'jNeuHadMultAK8':         ["nNeuHad(J)",                                       25,     0.0,    25.0,        1,     1],
'jMuMultAK8':             ["nMu(J)",                                           8,      0.0,    8.0,         1,     1],
'jEleMultAK8':            ["nEle(J)",                                          8,      0.0,    8.0,         1,     1],
'jChHadMultAK8':          ["nChHad(J)",                                        145,    0.0,    145.0,       1,     1],
'jChMultAK8':             ["nCh(J)",                                           145,    0.0,    145.0,       1,     1],
'jMultAK8':               ["mult(J)",                                          250,    0.0,    250.0,       1,     1],
'dPhijMETAK8':            [r"$\Delta\phi(J,MET)$",                             100,    0.0,    4.0,         1,     1],
'dPhiMinjMETAK8':         [r"$\Delta\phi_{min}(j,MET)$",                       100,    0.0,    4.0,         2,     0],
'dEtaj12AK8':             [r"$\Delta\eta(J_{1},J_{2})$",                       200,    0.0,    10.0,        2,     0],
'dRJ12AK8':               [r"$\Delta R(J_{1},J_{2})$",                         100,    0.0,    10.0,        2,     0],
'mT':                     [r"$m_{T} (GeV)$",                                   500,    0.0,    5000.0,      2,     0],
'METrHT_pt30':            [r"$MET/H_{T}$",                                     100,    0.0,    3.0,         2,     0],
'METrST_pt30':            [r"$MET/S_{T}",                                      100,    0.0,    1.0,         2,     0],
'j1Pt':                   [r"$p_{T}$ [GeV]",                                   200,    0.0,    2000.0,      0,     0],
'j1Eta':                  [r"$\eta$",                                          200,   -6.0,    6.0,         0,     0],
'j1Phi':                  [r"$\phi$",                                          200,   -4.0,    4.0,         0,     0],
'j1Axismajor':            [r"$\sigma_{major}(j)$",                             40,     0.0,    0.5,         0,     0],
'j1Axisminor':            [r"$\sigma_{minor}(j)$",                             40,     0.0,    0.3,         0,     0],
'j1PtD':                  ["ptD",                                              40,     0.0,    1.2,         0,     0],
'dPhij1MET':              [r"$\Delta\phi(J_{1},MET)$",                         100,    0.0,    4.0,         0,     0],
'j2Pt':                   [r"$p_{T}$ [GeV]",                                   200,    0.0,    2000.0,      0,     0],
'j2Eta':                  [r"$\eta$",                                          200,   -6.0,    6.0,         0,     0],
'j2Phi':                  [r"$\phi$",                                          200,   -4.0,    4.0,         0,     0],
'j2Axismajor':            [r"$\sigma_{major}(j)$",                             40,     0.0,    0.5,         0,     0],
'j2Axisminor':            [r"$\sigma_{minor}(j)$",                             40,     0.0,    0.3,         0,     0],
'j2PtD':                  ["ptD",                                              40,     0.0,    1.2,         0,     0],
'dPhij2MET':              [r"$\Delta\phi(J_{1},MET)$",                         100,    0.0,    4.0,         0,     0],
'dPhij1rdPhij2':          [r"$\Delta\phi(J_{1},MET)/\Delta\phi(J_{2},MET)$",   100,    0.0,    100.0,       0,     0],
'j1PtAK8':                [r"$p_{T}$ [GeV]",                                   200,    0.0,    2000.0,      2,     0],
'j1EtaAK8':               [r"$\eta$",                                          200,   -6.0,    6.0,         2,     0],
'j1PhiAK8':               [r"$\phi$",                                          200,   -4.0,    4.0,         2,     0],
'j1AxismajorAK8':         [r"$\sigma_{major}(j)$",                             40,     0.0,    0.5,         2,     0],
'j1AxisminorAK8':         [r"$\sigma_{minor}(j)$",                             40,     0.0,    0.3,         2,     0],
'j1GirthAK8':             ["girth(j)",                                         40,     0.0,    0.5,         2,     0],
'j1PtDAK8':               ["ptD",                                              40,     0.0,    1.2,         2,     0],
'j1Tau1AK8':              [r"$	au_{1}(j)$",                                   40,     0.0,    0.8,         2,     0],
'j1Tau2AK8':              [r"$	au_{2}(j)$",                                   40,     0.0,    0.65,        2,     0],
'j1Tau3AK8':              [r"$	au_{3}(j)$",                                   40,     0.0,    0.35,        2,     0],
'j1Tau21AK8':             [r"$	au_{21}(j)$",                                  40,     0.0,    1.3,         2,     0],
'j1Tau32AK8':             [r"$	au_{32}(j)$",                                  40,     0.0,    1.3,         2,     0],
'j1SoftDropMassAK8':      [r"$m_{SD}(j)$",                                     40,     0.0,    200,         2,     0],
'dPhij1METAK8':           [r"$\Delta\phi(J_{1},MET)$",                         100,    0.0,    4.0,         2,     0],
'j2PtAK8':                [r"$p_{T}$ [GeV]",                                   200,    0.0,    2000.0,      2,     0],
'j2EtaAK8':               [r"$\eta$",                                          200,   -6.0,    6.0,         2,     0],
'j2PhiAK8':               [r"$\phi$",                                          200,   -4.0,    4.0,         2,     0],
'j2AxismajorAK8':         [r"$\sigma_{major}(j)$",                             40,     0.0,    0.5,         2,     0],
'j2AxisminorAK8':         [r"$\sigma_{minor}(j)$",                             40,     0.0,    0.3,         2,     0],
'j2GirthAK8':             ["girth(j)",                                         40,     0.0,    0.5,         2,     0],
'j2PtDAK8':               ["ptD",                                              40,     0.0,    1.2,         2,     0],
'j2Tau1AK8':              [r"$	au_{1}(j)$",                                   40,     0.0,    0.8,         2,     0],
'j2Tau2AK8':              [r"$	au_{2}(j)$",                                   40,     0.0,    0.65,        2,     0],
'j2Tau3AK8':              [r"$	au_{3}(j)$",                                   40,     0.0,    0.35,        2,     0],
'j2Tau21AK8':             [r"$	au_{21}(j)$",                                  40,     0.0,    1.3,         2,     0],
'j2Tau32AK8':             [r"$	au_{32}(j)$",                                  40,     0.0,    1.3,         2,     0],
'j2SoftDropMassAK8':      [r"$m_{SD}(j)$",                                     40,     0.0,    200,         2,     0],
'dPhij2METAK8':           [r"$\Delta\phi(J_{1},MET)$",                         100,    0.0,    4.0,         2,     0],
'dPhij1rdPhij2AK8':       [r"$\Delta\phi(J_{1},MET)/\Delta\phi(J_{2},MET)$",   100,    0.0,    100.0,       2,     0],
# 'mT2_f4_msm':             [r"$m_{T2} (GeV)$",                                  500,    0.0,    5000.0,      2,     0],
# 'mT2_f4_msm_dEta':        [r"$m_{T2} (GeV)$",                                  500,    0.0,    5000.0,      2,     0],
# 'mT2_f4_msm_dPhi':        [r"$m_{T2} (GeV)$",                                  500,    0.0,    5000.0,      2,     0],
# 'mT2_f4_msm_dR':          [r"$m_{T2} (GeV)$",                                  500,    0.0,    5000.0,      2,     0],
# "GenJetsAK8_hvCategory":  ["GenJetAK8 hvCategory",                             32,     0.0,    32.0,        0,     1],
# "JetsAK8_hvCategory":     ["JetAK8 hvCategory",                                32,     0.0,    32.0,        1,     2],
# "GenMT2_AK8":             [r"$m_{T2} (GeV)$",                                  500,    0.0,    5000.0,      0,     0],
# "GenJetsAK8_darkPtFrac":  ["GenJetAK8 Dark pT Fraction",                       100,    0.0,    1.0,         0,     1],
'nsvjJetsAK8':            ["Number of SVJ AK8Jets",                            20,     0.0,    20.0,        0,     0],
'nnOutput':               ["nnOutput",                                         100,    0.0,    1.0,         0,     1],
}
