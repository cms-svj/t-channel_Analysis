from coffea import hist, processor
import numpy as np
import awkward1 as ak
import awkward
from utils import utility as utl
from utils.objects import Objects
from utils import baseline as bl

# add vector of Trigger names and Trigger pass: want to store the events that pass at least one trigger and then use
# those events to plot all the kinematic variables. Can be done interactively on one signal first.

# bugs: if we run over too few events, we may have an empty array for variables like jets.pt,
# there will be a lot of bugs when that happens. For example, ak.broadcast_arrays returns error
# if the array is empty
class MainProcessor(processor.ProcessorABC):
        def __init__(self):
                # dataset_axis = hist.Cat("dataset", "Primary dataset")
                # pt_axis = hist.Bin("pt", r"$p_{T}$ [GeV]", 40, 0, 3500)
                eventCounter_axis           = hist.Bin("EventCounter",          "EventCounter",                                     2,      -1.1,   1.1 )
                njet_axis                   = hist.Bin("njets",                 "Number of Jets",                                   20,     0.0,    20.0)
                nb_axis                     = hist.Bin("nb",                    "Number of b",                                      10,     0.0,    10.0)
                nl_axis                     = hist.Bin("nl",                    "Number of Leptons",                                10,     0.0,    10.0)
                ht_axis                     = hist.Bin("ht",                    r"$H_{T}$ (GeV)",                                   500,    0.0,    5000.0)
                st_axis                     = hist.Bin("st",                    r"$S_{T}$ (GeV)",                                   500,    0.0,    5000.0)
                MET_axis                    = hist.Bin("MET",                   "MET [GeV]",                                        500,    0.0,    2000.0)
                pt_axis                     = hist.Bin("pt",                    r"$p_{T}$ [GeV]",                                   200,    0.0,    2000.0)
                eta_axis                    = hist.Bin("eta",                   r"$\eta$",                                          200,   -6.0,    6.0)
                phi_axis                    = hist.Bin("phi",                   r"$\phi$",                                          200,   -4.0,    4.0)
                axismajor_axis              = hist.Bin("axismajor",             r"$\sigma_{major}(j)$",                             40,     0.0,    0.5)
                axisminor_axis              = hist.Bin("axisminor",             r"$\sigma_{minor}(j)$",                             40,     0.0,    0.3)
                girth_axis                  = hist.Bin("girth",                 "girth(j)",                                         40,     0.0,    0.5)
                ptD_axis                    = hist.Bin("ptD",                   "ptD",                                              40,     0.0,    1.2)
                tau1_axis                   = hist.Bin("tau1",                  r"$\tau_{1}(j)$",                                   40,     0.0,    0.8)
                tau2_axis                   = hist.Bin("tau2",                  r"$\tau_{2}(j)$",                                   40,     0.0,    0.65)
                tau3_axis                   = hist.Bin("tau3",                  r"$\tau_{3}(j)$",                                   40,     0.0,    0.35)
                tau21_axis                  = hist.Bin("tau21",                 r"$\tau_{21}(j)$",                                  40,     0.0,    1.3)
                tau32_axis                  = hist.Bin("tau32",                 r"$\tau_{32}(j)$",                                  40,     0.0,    1.3)
                softDropMass_axis           = hist.Bin("softDropMass",          r"$m_{SD}(j)$",                                     40,     0.0,    200)
                weight_axis                 = hist.Bin("Weight",                "Weight",                                           200,   -5.0,    5.0)
                dEtaJ12_axis                = hist.Bin("dEtaJ12",               r"$\Delta\eta(J_{1},J_{2})$",                       200,    0.0,    10.0)
                dRJ12_axis                  = hist.Bin("dRJ12",                 r"$\Delta R(J_{1},J_{2})$",                         100,    0.0,   10.0)
                dPhiJMET_axis               = hist.Bin("dPhiJMET",              r"$\Delta\phi(J_{1},MET)$",                         100,    0.0,    4.0)
                dPhiJ1METrdPhiJ2MET_axis    = hist.Bin("dPhiJ1METrdPhiJ2MET",   r"$\Delta\phi(J_{1},MET)/\Delta\phi(J_{2},MET)$",   100,    0.0,    100.0)
                mjjM_axis                   = hist.Bin("mjjM",                  r"$M(J_{1},J_{2}) (GeV)$",                          500,    0.0,    5000.0)
                mjjPt_axis                  = hist.Bin("mjjPt",                 r"$p_{T}(J_{1},J_{2}) (GeV)$",                      200,    0.0,    2000.0)
                mjjEta_axis                 = hist.Bin("mjjEta",                r"$\eta(J_{1},J_{2})$",                             200,   -6.0,    6.0)
                mT_axis                     = hist.Bin("mT",                    r"$m_{T} (GeV)$",                                   500,    0.0,    5000.0)
                METrHT_pt30_axis            = hist.Bin("METrHT_pt30",           r"$MET/H_{T}$",                                     100,    0.0,    20.0)
                METrST_pt30_axis            = hist.Bin("METrST_pt30",           r"$MET/S_{T}",                                      100,    0.0,    1.0)

                self._accumulator = processor.dict_accumulator({
                        # 'jtpt':hist.Hist("Counts", dataset_axis, pt_axis),
                        # 'jteta':hist.Hist("Counts",dataset_axis,eta_axis),
                        'EventCounter':                     hist.Hist("EventCounter",                   eventCounter_axis),
                        'h_weight':                         hist.Hist("Weight",                         weight_axis),
                        'h_njets':                          hist.Hist("h_njets",                        njet_axis),
                        'h_njetsAK8':                       hist.Hist("h_njetsAK8",                     njet_axis),
                        'h_nb':                             hist.Hist("h_nb",                           nb_axis),
                        'h_nl':                             hist.Hist("h_nl",                           nl_axis),
                        'h_ht':                             hist.Hist("h_ht",                           ht_axis),
                        'h_st':                             hist.Hist("h_st",                           st_axis),
                        'h_met':                            hist.Hist("h_met",                          MET_axis),
                        'h_jPt':                            hist.Hist("h_jPt",                          pt_axis),
                        'h_jEta':                           hist.Hist("h_jEta",                         eta_axis),
                        'h_jPhi':                           hist.Hist("h_jPhi",                         phi_axis),
                        'h_jAxismajor':                     hist.Hist("h_jAxismajor",                   axismajor_axis),
                        'h_jAxisminor':                     hist.Hist("h_jAxisminor",                   axisminor_axis),
                        'h_jPtD':                           hist.Hist("h_jPtD",                         ptD_axis),
                        'h_jPtAK8':                         hist.Hist("h_jPtAK8",                       pt_axis),
                        'h_jEtaAK8':                        hist.Hist("h_jEtaAK8",                      eta_axis),
                        'h_jPhiAK8':                        hist.Hist("h_jPhiAK8",                      phi_axis),
                        'h_jAxismajorAK8':                  hist.Hist("h_jAxismajorAK8",                axismajor_axis),
                        'h_jAxisminorAK8':                  hist.Hist("h_jAxisminorAK8",                axisminor_axis),
                        'h_jGirthAK8':                      hist.Hist("h_jGirthAK8",                    girth_axis),
                        'h_jPtDAK8':                        hist.Hist("h_jPtDAK8",                      ptD_axis),
                        'h_jTau1AK8':                       hist.Hist("h_jTau1AK8",                     tau1_axis),
                        'h_jTau2AK8':                       hist.Hist("h_jTau2AK8",                     tau2_axis),
                        'h_jTau3AK8':                       hist.Hist("h_jTau3AK8",                     tau3_axis),
                        'h_jTau21AK8':                      hist.Hist("h_jTau21AK8",                    tau21_axis),
                        'h_jTau32AK8':                      hist.Hist("h_jTau32AK8",                    tau32_axis),
                        'h_jSoftDropMassAK8':               hist.Hist("h_jSoftDropMassAK8",             softDropMass_axis),
                        # 'h_mjjM':                   hist.Hist("h_mjjM",                     mjjM_axis),
                        # 'h_mjjPt':                  hist.Hist("h_mjjPt",                    mjjPt_axis),
                        # 'h_mjjEta':                 hist.Hist("h_mjjEta",                   mjjEta_axis),
                        'h_mT':                             hist.Hist("h_mT",                           mT_axis),
                        'h_METrHT_pt30':                    hist.Hist("h_METrHT_pt30",                  METrHT_pt30_axis),
                        'h_METrST_pt30':                    hist.Hist("h_METrST_pt30",                  METrST_pt30_axis),
                        # Cut: at least 2 AK8 jets
                        'h_njets_ge2AK8j':                  hist.Hist("h_njets_ge2AK8j",                njet_axis),
                        'h_njetsAK8_ge2AK8j':               hist.Hist("h_njetsAK8_ge2AK8j",             njet_axis),
                        'h_ht_ge2AK8j':                     hist.Hist("h_ht_ge2AK8j",                   ht_axis),
                        'h_st_ge2AK8j':                     hist.Hist("h_st_ge2AK8j",                   st_axis),
                        'h_met_ge2AK8j':                    hist.Hist("h_met_ge2AK8j",                  MET_axis),
                        'h_dEtaJ12_ge2AK8j':                hist.Hist("h_dEtaJ12_ge2AK8j",              dEtaJ12_axis),
                        'h_dRJ12_ge2AK8j':                  hist.Hist("h_dRJ12_ge2AK8j",                dRJ12_axis),
                        'h_dPhiJ1MET_ge2AK8j':              hist.Hist("h_dPhiJ1MET_ge2AK8j",            dPhiJMET_axis),
                        'h_dPhiJ2MET_ge2AK8j':              hist.Hist("h_dPhiJ2MET_ge2AK8j",            dPhiJMET_axis),
                        'h_dPhiJ1METrdPhiJ2MET_ge2AK8j':    hist.Hist("h_dPhiJ1METrdPhiJ2MET_ge2AK8j",  dPhiJ1METrdPhiJ2MET_axis),
                        'h_dPhiMinJMET_ge2AK8j':            hist.Hist("h_dPhiMinJMET_ge2AK8j",          dPhiJMET_axis),
                        'h_METrHT_pt30_ge2AK8j':            hist.Hist("h_METrHT_pt30_ge2AK8j",          METrHT_pt30_axis),
                        'h_METrST_pt30_ge2AK8j':            hist.Hist("h_METrST_pt30_ge2AK8j",          METrST_pt30_axis),
                        'h_mT_ge2AK8j':                     hist.Hist("h_mT_ge2AK8j",                   mT_axis),
                        'h_j1PtAK8_ge2AK8j':                hist.Hist("h_j1PtAK8_ge2AK8j",              pt_axis),
                        'h_j1EtaAK8_ge2AK8j':               hist.Hist("h_j1EtaAK8_ge2AK8j",             eta_axis),
                        'h_j1PhiAK8_ge2AK8j':               hist.Hist("h_j1PhiAK8_ge2AK8j",             phi_axis),
                        'h_j1AxismajorAK8_ge2AK8j':         hist.Hist("h_j1AxismajorAK8_ge2AK8j",       axismajor_axis),
                        'h_j1AxisminorAK8_ge2AK8j':         hist.Hist("h_j1AxisminorAK8_ge2AK8j",       axisminor_axis),
                        'h_j1GirthAK8_ge2AK8j':             hist.Hist("h_j1GirthAK8_ge2AK8j",           girth_axis),
                        'h_j1PtDAK8_ge2AK8j':               hist.Hist("h_j1PtDAK8_ge2AK8j",             ptD_axis),
                        'h_j1Tau1AK8_ge2AK8j':              hist.Hist("h_j1Tau1AK8_ge2AK8j",            tau1_axis),
                        'h_j1Tau2AK8_ge2AK8j':              hist.Hist("h_j1Tau2AK8_ge2AK8j",            tau2_axis),
                        'h_j1Tau3AK8_ge2AK8j':              hist.Hist("h_j1Tau3AK8_ge2AK8j",            tau3_axis),
                        'h_j1Tau21AK8_ge2AK8j':             hist.Hist("h_j1Tau21AK8_ge2AK8j",           tau21_axis),
                        'h_j1Tau32AK8_ge2AK8j':             hist.Hist("h_j1Tau32AK8_ge2AK8j",           tau32_axis),
                        'h_j1SoftDropMassAK8_ge2AK8j':      hist.Hist("h_j1SoftDropMassAK8_ge2AK8j",    softDropMass_axis),
                        'h_j2PtAK8_ge2AK8j':                hist.Hist("h_j2PtAK8_ge2AK8j",              pt_axis),
                        'h_j2EtaAK8_ge2AK8j':               hist.Hist("h_j2EtaAK8_ge2AK8j",             eta_axis),
                        'h_j2PhiAK8_ge2AK8j':               hist.Hist("h_j2PhiAK8_ge2AK8j",             phi_axis),
                        'h_j2AxismajorAK8_ge2AK8j':         hist.Hist("h_j2AxismajorAK8_ge2AK8j",       axismajor_axis),
                        'h_j2AxisminorAK8_ge2AK8j':         hist.Hist("h_j2AxisminorAK8_ge2AK8j",       axisminor_axis),
                        'h_j2GirthAK8_ge2AK8j':             hist.Hist("h_j2GirthAK8_ge2AK8j",           girth_axis),
                        'h_j2PtDAK8_ge2AK8j':               hist.Hist("h_j2PtDAK8_ge2AK8j",             ptD_axis),
                        'h_j2Tau1AK8_ge2AK8j':              hist.Hist("h_j2Tau1AK8_ge2AK8j",            tau1_axis),
                        'h_j2Tau2AK8_ge2AK8j':              hist.Hist("h_j2Tau2AK8_ge2AK8j",            tau2_axis),
                        'h_j2Tau3AK8_ge2AK8j':              hist.Hist("h_j2Tau3AK8_ge2AK8j",            tau3_axis),
                        'h_j2Tau21AK8_ge2AK8j':             hist.Hist("h_j2Tau21AK8_ge2AK8j",           tau21_axis),
                        'h_j2Tau32AK8_ge2AK8j':             hist.Hist("h_j2Tau32AK8_ge2AK8j",           tau32_axis),
                        'h_j2SoftDropMassAK8_ge2AK8j':      hist.Hist("h_j2SoftDropMassAK8_ge2AK8j",    softDropMass_axis),
                        # At least 2 AK4 jets
                        'h_njets_ge2AK4j':                  hist.Hist("h_njets_ge2AK4j",                njet_axis),
                        'h_njetsAK8_ge2AK4j':               hist.Hist("h_njetsAK8_ge2AK4j",             njet_axis),
                        'h_ht_ge2AK4j':                     hist.Hist("h_ht_ge2AK4j",                   ht_axis),
                        'h_st_ge2AK4j':                     hist.Hist("h_st_ge2AK4j",                   st_axis),
                        'h_met_ge2AK4j':                    hist.Hist("h_met_ge2AK4j",                  MET_axis),
                        'h_dEtaj12_ge2AK4j':                hist.Hist("h_dEtaj12_ge2AK4j",              dEtaJ12_axis),
                        'h_dRj12_ge2AK4j':                  hist.Hist("h_dRj12_ge2AK4j",                dRJ12_axis),
                        'h_dPhij1MET_ge2AK4j':              hist.Hist("h_dPhij1MET_ge2AK4j",            dPhiJMET_axis),
                        'h_dPhij2MET_ge2AK4j':              hist.Hist("h_dPhij2MET",                    dPhiJMET_axis),
                        'h_dPhij1METrdPhij2MET_ge2AK4j':    hist.Hist("h_dPhij1METrdPhij2MET",          dPhiJ1METrdPhiJ2MET_axis),
                        'h_dPhiMinjMET_ge2AK4j':            hist.Hist("h_dPhiMinjMET_ge2AK4j",          dPhiJMET_axis),
                        'h_METrHT_pt30_ge2AK4j':            hist.Hist("h_METrHT_pt30_ge2AK4j",          METrHT_pt30_axis),
                        'h_METrST_pt30_ge2AK4j':            hist.Hist("h_METrST_pt30_ge2AK4j",          METrST_pt30_axis),
                        'h_mT_ge2AK4j':                     hist.Hist("h_mT_ge2AK4j",                   mT_axis),
                        'h_j1Pt_ge2AK4j':                   hist.Hist("h_j1Pt_ge2AK4j",                 pt_axis),
                        'h_j1Eta_ge2AK4j':                  hist.Hist("h_j1Eta_ge2AK4j",                eta_axis),
                        'h_j1Phi_ge2AK4j':                  hist.Hist("h_j1Phi_ge2AK4j",                phi_axis),
                        'h_j1Axismajor_ge2AK4j':            hist.Hist("h_j1Axismajor_ge2AK4j",          axismajor_axis),
                        'h_j1Axisminor_ge2AK4j':            hist.Hist("h_j1Axisminor_ge2AK4j",          axisminor_axis),
                        'h_j1PtD_ge2AK4j':                  hist.Hist("h_j1PtD_ge2AK4j",                ptD_axis),
                        'h_j2Pt_ge2AK4j':                   hist.Hist("h_j2Pt_ge2AK4j",                 pt_axis),
                        'h_j2Eta_ge2AK4j':                  hist.Hist("h_j2Eta_ge2AK4j",                eta_axis),
                        'h_j2Phi_ge2AK4j':                  hist.Hist("h_j2Phi_ge2AK4j",                phi_axis),
                        'h_j2Axismajor_ge2AK4j':            hist.Hist("h_j2Axismajor_ge2AK4j",          axismajor_axis),
                        'h_j2Axisminor_ge2AK4j':            hist.Hist("h_j2Axisminor_ge2AK4j",          axisminor_axis),
                        'h_j2PtD_ge2AK4j':                  hist.Hist("h_j2PtD_ge2AK4j",                ptD_axis),
                        # for verifying TTStitching
                        'h_genMET':                         hist.Hist("h_genMET",                       MET_axis),
                        'h_madHT':                          hist.Hist("madHT_Stitched",                 ht_axis),
                        'cutflow':                          processor.defaultdict_accumulator(int),
                })

        @property
        def accumulator(self):
                return self._accumulator

        def process(self, df):
                output = self.accumulator.identity()

                # important variables
                obj = Objects(df)
                mask = bl.getBaselineMask(df)

                electrons = obj.goodElectrons()[mask]
                muons = obj.goodMuons()[mask]
                jets = obj.goodJets()[mask]
                fjets = obj.goodFatJets()[mask]
                bjets = obj.goodBJets(df,jets)

                madHT_cut = df['madHT'][mask]
                genMET = df['GenMET'][mask]
                met = df['MET'][mask]
                metPhi = df['METPhi'][mask]
                mtAK8 = df['MT_AK8'][mask]
                output['cutflow']['all events'] += fjets.size

                # defining weights
                luminosity = 21071.0+38654.0
                evtw = df['Weight'][mask]*luminosity
                eCounter = np.where(evtw >= 0, 1, -1)
                ew = utl.awkwardReshape(electrons,evtw)
                mw = utl.awkwardReshape(muons,evtw)
                jw = utl.awkwardReshape(jets,evtw)
                fjw = utl.awkwardReshape(fjets,evtw)
                # bjw = utl.awkwardReshape(bjets,evtw)

                # Getting subset of variables based on number of AK8 jets
                # calculating event variables
                ht = ak.sum(jets.pt,axis=1)
                st = ht + met
                dPhiMinJ = utl.deltaPhi(fjets.phi,metPhi).min()
                dPhiMinj = utl.deltaPhi(jets.phi,metPhi).min()

                ## at least 2 AK4 Jets
                cut_2j = jets.counts >= 2
                jets_2j = jets[cut_2j]
                metPhi_2j = metPhi[cut_2j]
                evtw_2j = evtw[cut_2j]
                j1_eta = np.array(jets_2j.eta[:,0])
                j2_eta = np.array(jets_2j.eta[:,1])
                j1_phi = np.array(jets_2j.phi[:,0])
                j2_phi = np.array(jets_2j.phi[:,1])
                dEtaj12 = abs(j1_eta - j2_eta)
                deltaR12j = utl.delta_R(j1_phi,j2_phi,j1_eta,j2_eta)
                dPhij1 = utl.deltaPhi(j1_phi,metPhi_2j)
                dPhij2 = utl.deltaPhi(j2_phi,metPhi_2j)

                ## at least 2 AK8 Jets
                cut_2fj = fjets.counts >= 2
                fjets_2fj = fjets[cut_2fj]
                metPhi_2fj = metPhi[cut_2fj]
                evtw_2fj = evtw[cut_2fj]
                J1_eta = np.array(fjets_2fj.eta[:,0])
                J2_eta = np.array(fjets_2fj.eta[:,1])
                J1_phi = np.array(fjets_2fj.phi[:,0])
                J2_phi = np.array(fjets_2fj.phi[:,1])
                dEtaJ12 = abs(J1_eta - J2_eta)
                deltaR12J = utl.delta_R(J1_phi,J2_phi,J1_eta,J2_eta)
                dPhiJ1 = utl.deltaPhi(J1_phi,metPhi_2fj)
                dPhiJ2 = utl.deltaPhi(J2_phi,metPhi_2fj)
                J1_tau1_2fj = fjets_2fj.tau1[:,0]
                J1_tau2_2fj = fjets_2fj.tau2[:,0]
                J1_tau3_2fj = fjets_2fj.tau3[:,0]
                J1_tau21_2fj = J1_tau2_2fj[J1_tau1_2fj>0]/J1_tau1_2fj[J1_tau1_2fj>0]
                J1_tau32_2fj = J1_tau3_2fj[J1_tau2_2fj>0]/J1_tau2_2fj[J1_tau2_2fj>0]
                J2_tau1_2fj = fjets_2fj.tau1[:,1]
                J2_tau2_2fj = fjets_2fj.tau2[:,1]
                J2_tau3_2fj = fjets_2fj.tau3[:,1]
                J2_tau21_2fj = J2_tau2_2fj[J2_tau1_2fj>0]/J2_tau1_2fj[J2_tau1_2fj>0]
                J2_tau32_2fj = J2_tau3_2fj[J2_tau2_2fj>0]/J2_tau2_2fj[J2_tau2_2fj>0]


                # twofjets = (fjets.counts >= 2)
                # difjets = fjets[twofjets]
                # ptcut = (difjets.pt[:,0] > 200) & (difjets.pt[:,1] > 200)
                # difjets_pt200 = difjets[ptcut]

                twofjets = (fjets.counts >= 2)
                output['cutflow']['two fjets'] += twofjets.sum()
                # difjets = fjets[twofjets]
                # difjets_pt200 = difjets[ptcut]
                # output['jtpt'].fill(dataset=dataset, pt=fjets.pt.flatten())
                # output['jteta'].fill(dataset=dataset, eta=fjets.eta.flatten())
                output["EventCounter"].fill(EventCounter=eCounter,weight=np.ones(len(evtw)))
                output["h_weight"].fill(Weight=evtw,weight=np.ones(len(evtw)))
                output['h_njets'].fill(njets=jets.counts.flatten(),weight=evtw)
                output['h_njetsAK8'].fill(njets=fjets.counts.flatten(),weight=evtw)
                output['h_nb'].fill(nb=bjets.counts,weight=evtw)
                output['h_nl'].fill(nl=(electrons.counts + muons.counts),weight=evtw)
                output['h_ht'].fill(ht=ht,weight=evtw)
                output['h_st'].fill(st=st,weight=evtw)
                output['h_met'].fill(MET=met,weight=evtw)
                output['h_jPt'].fill(pt=jets.pt.flatten(),weight=ak.flatten(jw))
                output['h_jEta'].fill(eta=jets.eta.flatten(),weight=ak.flatten(jw))
                output['h_jPhi'].fill(phi=jets.phi.flatten(),weight=ak.flatten(jw))
                output['h_jAxismajor'].fill(axismajor=jets.axismajor.flatten(),weight=ak.flatten(jw))
                output['h_jAxisminor'].fill(axisminor=jets.axisminor.flatten(),weight=ak.flatten(jw))
                output['h_jPtD'].fill(ptD=jets.ptD.flatten(),weight=ak.flatten(jw))
                output['h_jPtAK8'].fill(pt=fjets.pt.flatten(),weight=ak.flatten(fjw))
                output['h_jEtaAK8'].fill(eta=fjets.eta.flatten(),weight=ak.flatten(fjw))
                output['h_jPhiAK8'].fill(phi=fjets.phi.flatten(),weight=ak.flatten(fjw))
                output['h_jAxismajorAK8'].fill(axismajor=fjets.axismajor.flatten(),weight=ak.flatten(fjw))
                output['h_jAxisminorAK8'].fill(axisminor=fjets.axisminor.flatten(),weight=ak.flatten(fjw))
                output['h_jGirthAK8'].fill(girth=fjets.girth.flatten(),weight=ak.flatten(fjw))
                output['h_jPtDAK8'].fill(ptD=fjets.ptD.flatten(),weight=ak.flatten(fjw))
                output['h_jTau1AK8'].fill(tau1=fjets.tau1.flatten(),weight=ak.flatten(fjw))
                output['h_jTau2AK8'].fill(tau2=fjets.tau2.flatten(),weight=ak.flatten(fjw))
                output['h_jTau3AK8'].fill(tau3=fjets.tau3.flatten(),weight=ak.flatten(fjw))
                output['h_jTau21AK8'].fill(tau21=fjets[fjets.tau1 > 0].tau2.flatten()/fjets[fjets.tau1 > 0].tau1.flatten(),weight=ak.flatten(fjw)[(fjets.tau1 > 0).flatten()])
                output['h_jTau32AK8'].fill(tau32=fjets[fjets.tau2 > 0].tau3.flatten()/fjets[fjets.tau2 > 0].tau2.flatten(),weight=ak.flatten(fjw)[(fjets.tau2 > 0).flatten()])
                output['h_jSoftDropMassAK8'].fill(softDropMass=fjets.softDropMass.flatten(),weight=ak.flatten(fjw))
                output['h_mT'].fill(mT=mtAK8,weight=evtw)
                output['h_METrHT_pt30'].fill(METrHT_pt30=met[ht>0]/ht[ht>0],weight=evtw[ht>0])
                output['h_METrST_pt30'].fill(METrST_pt30=met[st>0]/st[st>0],weight=evtw[st>0])
                # Cut: at least 2 AK8 jets
                output['h_njets_ge2AK8j'].fill(njets=jets[cut_2fj].counts.flatten(),weight=evtw_2fj)
                output['h_njetsAK8_ge2AK8j'].fill(njets=fjets[cut_2fj].counts.flatten(),weight=evtw_2fj)
                output['h_ht_ge2AK8j'].fill(ht=ht[cut_2fj],weight=evtw_2fj)
                output['h_st_ge2AK8j'].fill(st=st[cut_2fj],weight=evtw_2fj)
                output['h_met_ge2AK8j'].fill(MET=met[cut_2fj],weight=evtw_2fj)
                output['h_dEtaJ12_ge2AK8j'].fill(dEtaJ12=dEtaJ12,weight=evtw_2fj)
                output['h_dRJ12_ge2AK8j'].fill(dRJ12=deltaR12J,weight=evtw_2fj)
                output['h_dPhiJ1MET_ge2AK8j'].fill(dPhiJMET=dPhiJ1,weight=evtw_2fj)
                output['h_dPhiJ2MET_ge2AK8j'].fill(dPhiJMET=dPhiJ2,weight=evtw_2fj)
                output['h_dPhiJ1METrdPhiJ2MET_ge2AK8j'].fill(dPhiJ1METrdPhiJ2MET=dPhiJ1[dPhiJ2>0]/dPhiJ2[dPhiJ2>0],weight=evtw_2fj[dPhiJ2>0])
                output['h_dPhiMinJMET_ge2AK8j'].fill(dPhiJMET=dPhiMinJ[(cut_2fj) & (dPhiMinJ<4.0)],weight=evtw[(cut_2fj) & (dPhiMinJ<4.0)])
                output['h_METrHT_pt30_ge2AK8j'].fill(METrHT_pt30=met[(cut_2fj) & (ht>0)]/ht[(cut_2fj) & (ht>0)],weight=evtw[(cut_2fj) & (ht>0)])
                output['h_METrST_pt30_ge2AK8j'].fill(METrST_pt30=met[(cut_2fj) & (st>0)]/st[(cut_2fj) & (st>0)],weight=evtw[(cut_2fj) & (st>0)])
                output['h_mT_ge2AK8j'].fill(mT=mtAK8[cut_2fj],weight=evtw_2fj)
                output['h_j1PtAK8_ge2AK8j'].fill(pt=fjets_2fj.pt[:,0],weight=evtw_2fj)
                output['h_j1EtaAK8_ge2AK8j'].fill(eta=fjets_2fj.eta[:,0],weight=evtw_2fj)
                output['h_j1PhiAK8_ge2AK8j'].fill(phi=fjets_2fj.phi[:,0],weight=evtw_2fj)
                output['h_j1AxismajorAK8_ge2AK8j'].fill(axismajor=fjets_2fj.axismajor[:,0],weight=evtw_2fj)
                output['h_j1AxisminorAK8_ge2AK8j'].fill(axisminor=fjets_2fj.axisminor[:,0],weight=evtw_2fj)
                output['h_j1GirthAK8_ge2AK8j'].fill(girth=fjets_2fj.girth[:,0],weight=evtw_2fj)
                output['h_j1PtDAK8_ge2AK8j'].fill(ptD=fjets_2fj.ptD[:,0],weight=evtw_2fj)
                output['h_j1Tau1AK8_ge2AK8j'].fill(tau1=fjets_2fj.tau1[:,0],weight=evtw_2fj)
                output['h_j1Tau2AK8_ge2AK8j'].fill(tau2=fjets_2fj.tau2[:,0],weight=evtw_2fj)
                output['h_j1Tau3AK8_ge2AK8j'].fill(tau3=fjets_2fj.tau3[:,0],weight=evtw_2fj)
                output['h_j1Tau21AK8_ge2AK8j'].fill(tau21=J1_tau21_2fj,weight=evtw_2fj[J1_tau1_2fj>0])
                output['h_j1Tau32AK8_ge2AK8j'].fill(tau32=J1_tau32_2fj,weight=evtw_2fj[J1_tau2_2fj>0])
                output['h_j1SoftDropMassAK8_ge2AK8j'].fill(softDropMass=fjets_2fj.softDropMass[:,0],weight=evtw_2fj)
                output['h_j2PtAK8_ge2AK8j'].fill(pt=fjets_2fj.pt[:,1],weight=evtw_2fj)
                output['h_j2EtaAK8_ge2AK8j'].fill(eta=fjets_2fj.eta[:,1],weight=evtw_2fj)
                output['h_j2PhiAK8_ge2AK8j'].fill(phi=fjets_2fj.phi[:,1],weight=evtw_2fj)
                output['h_j2AxismajorAK8_ge2AK8j'].fill(axismajor=fjets_2fj.axismajor[:,1],weight=evtw_2fj)
                output['h_j2AxisminorAK8_ge2AK8j'].fill(axisminor=fjets_2fj.axisminor[:,1],weight=evtw_2fj)
                output['h_j2GirthAK8_ge2AK8j'].fill(girth=fjets_2fj.girth[:,1],weight=evtw_2fj)
                output['h_j2PtDAK8_ge2AK8j'].fill(ptD=fjets_2fj.ptD[:,1],weight=evtw_2fj)
                output['h_j2Tau1AK8_ge2AK8j'].fill(tau1=fjets_2fj.tau1[:,1],weight=evtw_2fj)
                output['h_j2Tau2AK8_ge2AK8j'].fill(tau2=fjets_2fj.tau2[:,1],weight=evtw_2fj)
                output['h_j2Tau3AK8_ge2AK8j'].fill(tau3=fjets_2fj.tau3[:,1],weight=evtw_2fj)
                output['h_j2Tau21AK8_ge2AK8j'].fill(tau21=J2_tau21_2fj,weight=evtw_2fj[J2_tau1_2fj>0])
                output['h_j2Tau32AK8_ge2AK8j'].fill(tau32=J2_tau32_2fj,weight=evtw_2fj[J2_tau2_2fj>0])
                output['h_j2SoftDropMassAK8_ge2AK8j'].fill(softDropMass=fjets_2fj.softDropMass[:,1],weight=evtw_2fj)
                # Cut: at least 2 AK4 jets
                output['h_njets_ge2AK4j'].fill(njets=jets[cut_2j].counts.flatten(),weight=evtw[cut_2j])
                output['h_njetsAK8_ge2AK4j'].fill(njets=fjets[cut_2j].counts.flatten(),weight=evtw[cut_2j])
                output['h_ht_ge2AK4j'].fill(ht=ht[cut_2j],weight=evtw[cut_2j])
                output['h_st_ge2AK4j'].fill(st=st[cut_2j],weight=evtw[cut_2j])
                output['h_met_ge2AK4j'].fill(MET=met[cut_2j],weight=evtw[cut_2j])
                output['h_dEtaj12_ge2AK4j'].fill(dEtaJ12=dEtaj12,weight=evtw_2j)
                output['h_dRj12_ge2AK4j'].fill(dRJ12=deltaR12j,weight=evtw_2j)
                output['h_dPhij1MET_ge2AK4j'].fill(dPhiJMET=dPhij1,weight=evtw_2j)
                output['h_dPhij2MET_ge2AK4j'].fill(dPhiJMET=dPhij2,weight=evtw_2j)
                output['h_dPhij1METrdPhij2MET_ge2AK4j'].fill(dPhiJ1METrdPhiJ2MET=dPhij1[dPhij2>0]/dPhij2[dPhij2>0],weight=evtw_2j[dPhij2>0])
                output['h_dPhiMinjMET_ge2AK4j'].fill(dPhiJMET=dPhiMinj[(cut_2j) & (dPhiMinj<4.0)],weight=evtw[(cut_2j) & (dPhiMinj<4.0)])
                output['h_METrHT_pt30_ge2AK4j'].fill(METrHT_pt30=met[(cut_2j) & (ht>0)]/ht[(cut_2j) & (ht>0)],weight=evtw[(cut_2j) & (ht>0)])
                output['h_METrST_pt30_ge2AK4j'].fill(METrST_pt30=met[(cut_2j) & (st>0)]/st[(cut_2j) & (st>0)],weight=evtw[(cut_2j) & (st>0)])
                output['h_mT_ge2AK4j'].fill(mT=mtAK8[cut_2j],weight=evtw[cut_2j])
                output['h_j1Pt_ge2AK4j'].fill(pt=jets_2j.pt[:,0],weight=evtw_2j)
                output['h_j1Eta_ge2AK4j'].fill(eta=jets_2j.eta[:,0],weight=evtw_2j)
                output['h_j1Phi_ge2AK4j'].fill(phi=jets_2j.phi[:,0],weight=evtw_2j)
                output['h_j1Axismajor_ge2AK4j'].fill(axismajor=jets_2j.axismajor[:,0],weight=evtw_2j)
                output['h_j1Axisminor_ge2AK4j'].fill(axisminor=jets_2j.axisminor[:,0],weight=evtw_2j)
                output['h_j1PtD_ge2AK4j'].fill(ptD=jets_2j.ptD[:,0],weight=evtw_2j)
                output['h_j2Pt_ge2AK4j'].fill(pt=jets_2j.pt[:,1],weight=evtw_2j)
                output['h_j2Eta_ge2AK4j'].fill(eta=jets_2j.eta[:,1],weight=evtw_2j)
                output['h_j2Phi_ge2AK4j'].fill(phi=jets_2j.phi[:,1],weight=evtw_2j)
                output['h_j2Axismajor_ge2AK4j'].fill(axismajor=jets_2j.axismajor[:,1],weight=evtw_2j)
                output['h_j2Axisminor_ge2AK4j'].fill(axisminor=jets_2j.axisminor[:,1],weight=evtw_2j)
                output['h_j2PtD_ge2AK4j'].fill(ptD=jets_2j.ptD[:,1],weight=evtw_2j)
                # for verifying TTStitching
                output['h_madHT'].fill(ht=madHT_cut,weight=evtw)
                output['h_genMET'].fill(MET=genMET,weight=evtw)
                return output

        def postprocess(self, accumulator):
                return accumulator
