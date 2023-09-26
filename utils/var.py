from itertools import combinations

def var(jNVar=False):
    myVars = {
        # 'eCounter':  ["h_eCounter", 2, -1.1, 1.1   ,     0,     0,       'w1'  ],
        # 'evtw':      [ "h_evtw",   2, -1.1, 1.1   ,     0,     0,       'evtw'],
        # 'jw':        [ "h_jw",     2, -1.1, 1.1   ,     0,     2,       'jw'  ],
        # 'fjw':       ["h_fjw",    2, -1.1, 1.1   ,     1,     2,       'fjw' ],
        'njets':          [ "Number of AK4Jets", 20, 0.0, 20.0, 2, 0, 'evtw'],
        'njetsAK8':       [ "Number of AK8Jets", 20, 0.0, 20.0, 2, 0, 'evtw'],
        # 'nb':             [ "Number of b",10, 0.0, 10.0, 2, 0, 'evtw'],
        # 'nl':             [ "Number of Leptons",  10, 0.0, 10.0, 2, 0, 'evtw'],
        # 'nnim':           [ "Number of NonIsoMuons", 10, 0.0, 10.0, 2, 0, 'evtw'],
        'ht':             [ r"H_{T} (GeV)", 500, 0.0, 5000.0, 2, 0,'evtw'],
        'st':             [ r"S_{T} (GeV)", 500, 0.0, 5000.0, 2, 0, 'evtw'],
        'met':            [ "p_{T}^{miss} (GeV)", 500, 0.0, 2000.0, 2, 0, 'evtw'],
        'metPhi':         [ r"#phi_{p_{T}^{miss}} ", 40, -4.0, 4.0, 0, 0, 'evtw'],
        # # 'madHT':                 [ r"H_{T} (GeV)",        500, 0.0,   5000.0,     2,     0,       'evtw'],
        'jPt':                     [r"p_{T} [GeV]",        200, 0.0,   4000.0,     0,     1,       'jw'  ],
        'jEta':                    [r"#eta",               200, -6.0,   6.0   ,     0,     1,       'jw'  ],
        'jPhi':                    [ r"#phi",               200, -3.5,   3.5   ,     0,     1,       'jw'  ],
        'jAxismajor':              [ r"#sigma_{major}(j)",  40,  0.0,   0.8   ,     0,     1,       'jw'  ],
        'jAxisminor':              [ r"#sigma_{minor}(j)",  40,  0.0,   0.3   ,     0,     1,       'jw'  ],
        'jPtD':                    [  "p_{t}^{D}",                   40,  0.0,   1.2   ,     0,     1,       'jw'  ],
        'dPhiMinjMET':             [ r"#Delta#phi_{min}(j,p_{T}^{miss})", 100, 0.0,   4.0   ,     0,     0,       'evtw'],
        'jPtAK8':                  [ r"p_{T}(J) (GeV)",        200, 0.0,   4000.0,     1,     1,       'fjw' ],
        'jEtaAK8':                 [ r"#eta(J)",               200, -6.0,   6.0   ,     1,     1,       'fjw' ],
        'jPhiAK8':                 [ r"#phi(J)",               200,  -4.0,   4.0   ,     1,     1,       'fjw' ],
        'jAxismajorAK8':           [ r"#sigma_{major}(J)",  40,  0.0,   0.8   ,     1,     1,       'fjw' ],
        'jAxisminorAK8':           [ r"#sigma_{minor}(J)",  40,  0.0,   0.4   ,     1,     1,       'fjw' ],

        # # 'jChEMEFractAK8':          [ "fChEM(J)",              50,  0.0,   1.0   ,     1,     1,       'fjw' ],
        # # 'jChHadEFractAK8':         [ "fChHad(J)",             50,  0.0,   1.0   ,     1,     1,       'fjw' ],
        # # 'jChHadMultAK8':           [ "nChHad(J)",             145, 0.0,   145.0 ,     1,     1,       'fjw' ],
        'jChMultAK8':              [ "nCh(J)",                145, 0.0,   145.0 ,     1,     1,       'fjw' ],
        # # 'jecfN2b1AK8':             [ "ecfN2b1(J)",            50,  0.0,   0.6   ,     1,     1,       'fjw' ],
        # # 'jecfN2b2AK8':             [ "ecfN2b2(J)",            50,  0.0,   0.4   ,     1,     1,       'fjw' ],
        # # 'jecfN3b1AK8':             [ "ecfN3b1(J)",            50,  0.0,   6.0   ,     1,     1,       'fjw' ],
        # # 'jecfN3b2AK8':             [ "ecfN3b2(J)",            50,  0.0,   5.0   ,     1,     1,       'fjw' ],
        # # 'jEleEFractAK8':           [ "fEle(J)",               50,  0.0,   1.0  ,     1,     1,       'fjw' ],
        # # 'jEleMultAK8':             [ "nEle(J)",               8,   0.0,   8.0   ,     1,     1,       'fjw' ],
        'jGirthAK8':               [ "girth(J)",              40,  0.0,   0.8  ,     1,     1,       'fjw' ],
        # # 'jHfEMEFractAK8':          [ "fHFEM(J)",              50,  0.0,   1.0   ,     1,     1,       'fjw' ],
        # # 'jHfHadEFractAK8':         [ "fHFHad(J)",             50,  0.0,   1.0   ,     1,     1,       'fjw' ],
        # # 'jMultAK8':                [ "mult(J)",               250, 0.0,   250.0 ,     1,     1,       'fjw' ],
        # # 'jMuEFractAK8':            [ "fMu(J)",                50,  0.0,   1.0   ,     1,     1,       'fjw' ],
        # # 'jMuMultAK8':              [ "nMu(J)",                8,   0.0,   10.0  ,     1,     1,       'fjw' ],
        # # 'jNeuEmEFractAK8':         [ "fNeuEM(J)",             50,  0.0,   1.0   ,     1,     1,       'fjw' ],
        # # 'jNeuHadEFractAK8':        [ "fNeuHad(J)",            50,  0.0,   1.0   ,     1,     1,       'fjw' ],
        # # 'jNeuHadMultAK8':          [ "nNeuHad(J)",            25,  0.0,   25.0  ,     1,     1,       'fjw' ],
        'jNeuMultAK8':             [ "nNeu(J)",               120, 0.0,   120.0 ,     1,     1,       'fjw' ],
        'jTau1AK8':                [ r"#tau_{1}(J)",        40,  0.0,   0.8   ,     1,     1,       'fjw' ],
        'jTau2AK8':                [ r"#tau_{2}(J)",        40,  0.0,   0.65  ,     1,     1,       'fjw' ],
        'jTau3AK8':                [ r"#tau_{3}(J)",        40,  0.0,   0.35  ,     1,     1,       'fjw' ],
        'jTau21AK8':               [ r"#tau_{21}(J)",       40,  0.0,   1.1   ,     1,     2,       'fjw' ],
        'jTau32AK8':               [ r"#tau_{32}(J)",       40,  0.0,   1.3   ,     1,     2,       'fjw' ],
        # # 'jPhoEFractAK8':           [ "fPho(J)",               50,  0.0,   1.0   ,     1,     1,       'fjw' ],
        # # 'jPhoMultAK8':             [ "nPho(J)",               110, 0.0,   110.0 ,     1,     1,       'fjw' ],
        'jPtDAK8':                 [ "p_{t}^{D}(J)",                   40,  0.0,   1.2  ,     1,     1,       'fjw' ],
        # # 'jSoftDropMassAK8':        [ r"m_{SD}(J)",          200, 0.0,   900   ,     1,     1,       'fjw' ],
        'dPhijMETAK8':             [ r"#Delta#phi(J,p_{T}^{miss})",  100, 0.0,   4.0   ,     1,     1,       'fjw' ],
        'dEtaj12AK8':              [ r"#Delta#eta(J_{1},J_{2})", 200, 0.0,   10.0  ,     2,     0,       'evtw'],
        'dRJ12AK8':                [ r"#Delta R(J_{1},J_{2})", 100, 0.0,   10.0  ,     2,     0,       'evtw'],
        'dPhiMinjMETAK8':          [ r"#Delta#phi_{min}(J,p_{T}^{miss})", 100, 0.0,   4.0   ,     2,     0,       'evtw'],
        'mT':                      [ r"m_{T} (GeV)",        500, 0.0,   6000.0,     2,     0,       'evtw'],
        'METrHT_pt30':             [ r"p_{T}^{miss}/H_{T}",          100, 0.0,   3.0   ,     2,     0,       'evtw'],
        'METrST_pt30':             [ r"p_{T}^{miss}/S_{T}",           100, 0.0,   1.0   ,     2,     0,       'evtw'],
        'dPhij1rdPhij2AK8':        [ r"#Delta#phi(J_{1},p_{T}^{miss})/#Delta#phi(J_{2},p_{T}^{miss})", 100, 0.0,   100.0 ,     2,     0,       'evtw'],
        # # 'electronsIso':            [ "electrons iso",         100, 0.0,   1.0   ,     0,     1,         'ew'],
        # # # 'muonsIso':                [ "muons iso",             100, 0.0,   1.0   ,     0,     1,         'mw'],
        # # #'nonIsoMuonsIso':          [ "NonIsoMuons iso",       200, 0.0,   10.0  ,     0,     1,       'nimw'],
        # # #'nonIsoMuonsPt':           [ "NonIsoMuons p_{T} [GeV]", 500, 0.0,   2500.0,     0,     1,       'nimw'],
        # # # 'mT2_f4_msm':            [ r"m_{T2} (GeV)",       500, 0.0,   5000.0,     2,     0,       'evtw'],
        # # # 'mT2_f4_msm_dEta':       [ r"m_{T2} (GeV)",       500, 0.0,   5000.0,     2,     0,       'evtw'],
        # # # 'mT2_f4_msm_dPhi':       [ r"m_{T2} (GeV)",       500, 0.0,   5000.0,     2,     0,       'evtw'],
        # # # 'mT2_f4_msm_dR':         [ r"m_{T2} (GeV)$,       500, 0.0,   5000.0,     2,     0,       'evtw'],
        # # #"GenJetsAK8_hvCategory":   [ "GenJetAK8 hvCategory",  32,  0.0,   32.0  ,     0,     1,        'gfjw'],
        # # "nNMedEvent":              [ "N-Med Event",           4,  0.0,    4.0  ,      0,     0,        'evtw'],
        # # "JetsAK8_hvCategory":      [ "JetAK8 hvCategory",     32,  0.0,   32.0  ,     1,     2,        'fjw'],
        # # # "GenMT2_AK8":            [ r"$m_{T2} (GeV)$",       500, 0.0,   5000.0,     0,     0,,      'evtw'],
        # # # "GenJetsAK8_darkPtFrac": [ "GenJetAK8 Dark pT Fraction", 100, 0.0,   1.0   ,     0,     1,        'fjw'],
        # # # 'nsvjJetsAK8':             [ "Number of SVJ AK8Jets", 20,  0.0,   20.0  ,     0,     0, 'evtw'      ],
        # # # 'nsvjJetsAK8Plus1':        [ "Number of SVJ AK8Jets pred1Jets",  20,  0.0,   20.0  ,     0,     0, 'pred1_evtw'],
        # # # 'nsvjJetsAK8Plus2':        [ "Number of SVJ AK8Jets pred2Jets",  20,  0.0,   20.0  ,     0,     0, 'pred2_evtw'],
        # # # 'nsvjJetsAK8Plus3':        [ "Number of SVJ AK8Jets pred3Jets",  20,  0.0,   20.0  ,     0,     0, 'pred3_evtw'],
        # # # 'nsvjJetsAK8Plus4':        [ "Number of SVJ AK8Jets pred4Jets",  20,  0.0,   20.0  ,     0,     0, 'pred4_evtw'],
        # # # 'nnOutput':                [ "nnOutput",              100, 0.0,   1.0   ,     0,     1, 'fjw'       ],
        # # # #'svjPtAK8':                [ r"p_{T} [GeV]",        280, 0.0,   2800.0,     0,     1, 'svfjw'     ],
        # # # #'svjEtaAK8':               [ r"#eta",               200, -6.0,   6.0   ,     0,     1, 'svfjw'     ],
        'electronPT':              [r"p_{T}(electron) (GeV)",         100,    0.0,    500.0,     0,     1,       'ew'],
        'electronPhi':              [r"\phi (electron) ",             40,    -4.0,    4.0   ,     0,     1,       'ew'],
        'electronEta':              [r"\eta (electron)",                   200,   -6.0,    6.0   ,     0,     1,       'ew' ],
        'muonPT':                   [r"p_{T}(muon) (GeV)",                 100,   0.0,  500.0,     0,     1,       'mw'],
        'muonPhi':                  [r"\phi (muon) ",                     40,   -4.0,    4.0   ,     0,     1,       'mw'],
        'muonEta':                   [r"\eta (muon)",                           200,   -6.0,    6.0   ,     0,     1,       'mw' ],
        

    }
    
    if jNVar:
        # preparing histograms for jN variables
        maxN = 4
        for i in range(maxN):
            jNList = {
                'j{}Pt'.format(i+1):              [  r"p_{T}(j_{"+str(i+1)+"}) [GeV]",          200,     0.0,     4000.0,       0,     0,       'evtw'],
                'j{}Eta'.format(i+1):             [   r"#eta(j_{"+str(i+1)+"})",               200,    -6.0,     6.0   ,       0,     0,       'evtw'],
                'j{}Phi'.format(i+1):             [  r"#phi(j_{"+str(i+1)+"})",               200,    -4.0,     4.0   ,       0,     0,       'evtw'],
                'j{}Axismajor'.format(i+1):       [    r"#sigma_{major}(j_{"+str(i+1)+"})",     40,      0.0,     0.5   ,       0,     0,       'evtw'],
                'j{}Axisminor'.format(i+1):       [    r"#sigma_{minor}(j_{"+str(i+1)+"})",     40,      0.0,     0.3   ,       0,     0,       'evtw'],
                'j{}PtD'.format(i+1):             [     "ptD(j_"+str(i+1)+")",                   40,      0.0,     1.2   ,       0,     0,       'evtw'],
                'dPhij{}MET'.format(i+1):         [ r"#Delta#phi(j_{"+str(i+1)+"},MET)",   100,     0.0,     4.0   ,       0,     0,       'evtw'],
                'j{}PtAK8'.format(i+1):           [ r"p_{T}(J_{"+str(i+1)+"}) [GeV]",          200,     0.0,     4000.0,       2,     0,       'evtw'],
                'j{}EtaAK8'.format(i+1):          [ r"#eta(J_{"+str(i+1)+"})",               200,    -6.0,     6.0   ,       2,     0,       'evtw'],
                'j{}PhiAK8'.format(i+1):          [  r"#phi(J_{"+str(i+1)+"})",               200,    -4.0,     4.0   ,       2,     0,       'evtw'],
                'j{}AxismajorAK8'.format(i+1):    [ r"#sigma_{major}(J_{"+str(i+1)+"})",     40,      0.0,     0.6   ,       2,     0,       'evtw'],
                'j{}AxisminorAK8'.format(i+1):    [ r"#sigma_{minor}(J_{"+str(i+1)+"})",     40,      0.0,     0.4   ,       2,     0,       'evtw'],
                'j{}GirthAK8'.format(i+1):        [ "girth(J_{"+str(i+1)+"})",                 40,      0.0,     0.8   ,       2,     0,       'evtw'],
                'j{}PtDAK8'.format(i+1):          [ "ptD(J_{"+str(i+1)+"})",                   40,      0.0,     1.2   ,       2,     0,       'evtw'],
                # 'j{}Tau1AK8'.format(i+1):         [  r"#tau_{1}(J_{"+str(i+1)+"})",           40,      0.0,     0.8   ,       2,     0,       'evtw'],
                # 'j{}Tau2AK8'.format(i+1):         [  r"#tau_{2}(J_{"+str(i+1)+"})",           40,      0.0,     0.65  ,       2,     0,       'evtw'],
                # 'j{}Tau3AK8'.format(i+1):         [ r"#tau_{3}(J_{"+str(i+1)+"})",           40,      0.0,     0.35  ,       2,     0,       'evtw'],
                'j{}Tau21AK8'.format(i+1):        [  r"#tau_{21}(J_{"+str(i+1)+"})",          40,      0.0,     1.3   ,       2,     0,       'evtw'],
                'j{}Tau32AK8'.format(i+1):        [r"#tau_{32}(J_{"+str(i+1)+"})",          40,      0.0,     1.3   ,       2,     0,       'evtw'],
                'j{}SoftDropMassAK8'.format(i+1): [ r"m_{SD}(J_{"+str(i+1)+"})",             200,     0.0,     900   ,       2,     0,       'evtw'],
                'dPhij{}METAK8'.format(i+1):      [   r"#Delta#phi(J_{"+str(i+1)+"},MET)",   100,     0.0,     4.0   ,       2,     0,       'evtw'],
                # "J{}_hvCategory".format(i+1):     [ "J_"+str(i+1)+" hvCategory",             32,      0.0,     32.0  ,       0,     0,       'evtw'],
            }
            myVars.update(jNList)

        allComs = list(combinations(range(maxN),2))
        for com in allComs:
            j1 = com[0]+1
            j2 = com[1]+1
            jNjMList = {
                # 'dEtaj{}{}'.format(j1,j2):         [ r"#Delta#eta(j_{"+str(j1)+"},j_{"+str(j2)+"})",                   200, 0.0, 5.3  , 0, 0, 'evtw'],
                # 'dPhij{}{}'.format(j1,j2):         [ r"#Delta#phi(j_{"+str(j1)+"},j_{"+str(j2)+"})",                   100, 0.0, 4.0  , 2, 0, 'evtw'],
                # 'dRj{}{}'.format(j1,j2):           [ r"#Delta R(j_{"+str(j1)+"},j_{"+str(j2)+"})",                     100, 0.0, 6.0  , 0, 0, 'evtw'],
                # 'dPhij{}rdPhij{}'.format(j1,j2):   [ r"#Delta#phi(j_{"+str(j1)+"},MET)/#Delta#phi(j_{"+str(j2)+"},MET)",100, 0.0, 120.0, 0, 0, 'evtw'],
                # 'dEtaj{}{}AK8'.format(j1,j2):      [ r"#Delta#eta(J_{"+str(j1)+"},J_{"+str(j2)+"})",                   200, 0.0, 5.3  , 2, 0, 'evtw'],
                # 'dPhij{}{}AK8'.format(j1,j2):      [ r"#Delta#phi(J_{"+str(j1)+"},J_{"+str(j2)+"})",                   100, 0.0, 4.0  , 2, 0, 'evtw'],
                # 'dRj{}{}AK8'.format(j1,j2):        [ r"#Delta R(J_{"+str(j1)+"},J_{"+str(j2)+"})",                     100, 0.0, 6.0  , 2, 0, 'evtw'],
                # 'dPhij{}rdPhij{}AK8'.format(j1,j2):[ r"#Delta#phi(J_{"+str(j1)+"},MET)/#Delta#phi(J_{"+str(j2)+"},MET)",100, 0.0, 120.0, 2, 0, 'evtw'],
            }
            myVars.update(jNjMList)

    #varXName, xbins, npzInfo, flattenInfo, weightName, varYName=None, ybins=None
    # histos2D = {
    #     'jPtAK8vsjEtaAK8': [varXName='jPtAK8', varYName='jEtaAK8',npzInfo=0, flattenInfo=1, weightName='fjw', 
    #                                     xbins= "Jet $p_{T}$", bins=280,  0.0,  2800.0), ybins=h.axis.Regular(name="y", label="Jet $\eta$", 200,  -6.0,  6.0)),
    # }
    # myVars.update(histos2D)

    return myVars
