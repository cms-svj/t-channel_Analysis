from coffea import hist, processor
import numpy as np
import awkward1 as ak
import awkward
from utils import utility as utl
from utils.objects import Objects
from utils import baseline as bl
from itertools import combinations

def col_accumulator(a):
    return processor.column_accumulator(np.array(a))

class MainProcessor(processor.ProcessorABC):
        def __init__(self):
                self._accumulator = processor.dict_accumulator({
                "evtw": processor.column_accumulator(np.zeros(shape=(0))),
                "jw": processor.column_accumulator(np.zeros(shape=(0))),
                "fjw": processor.column_accumulator(np.zeros(shape=(0))),
                "njets": processor.column_accumulator(np.zeros(shape=(0))),
                "njetsAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "nb": processor.column_accumulator(np.zeros(shape=(0))),
                "nl": processor.column_accumulator(np.zeros(shape=(0))),
                "ht": processor.column_accumulator(np.zeros(shape=(0))),
                "st": processor.column_accumulator(np.zeros(shape=(0))),
                "met": processor.column_accumulator(np.zeros(shape=(0))),
                "madHT": processor.column_accumulator(np.zeros(shape=(0))),
                "jPt": processor.column_accumulator(np.zeros(shape=(0))),
                "jEta": processor.column_accumulator(np.zeros(shape=(0))),
                "jPhi": processor.column_accumulator(np.zeros(shape=(0))),
                "jAxismajor": processor.column_accumulator(np.zeros(shape=(0))),
                "jAxisminor": processor.column_accumulator(np.zeros(shape=(0))),
                "jPtD": processor.column_accumulator(np.zeros(shape=(0))),
                "dPhiMinjMET": processor.column_accumulator(np.zeros(shape=(0))),
                "dEtaj12": processor.column_accumulator(np.zeros(shape=(0))),
                "dRJ12": processor.column_accumulator(np.zeros(shape=(0))),
                "jPtAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "jEtaAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "jPhiAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "jAxismajorAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "jAxisminorAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "jGirthAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "jPtDAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "jTau1AK8": processor.column_accumulator(np.zeros(shape=(0))),
                "jTau2AK8": processor.column_accumulator(np.zeros(shape=(0))),
                "jTau3AK8": processor.column_accumulator(np.zeros(shape=(0))),
                "jTau21AK8": processor.column_accumulator(np.zeros(shape=(0))),
                "jTau32AK8": processor.column_accumulator(np.zeros(shape=(0))),
                "jSoftDropMassAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "dPhiMinjMETAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "dEtaj12AK8": processor.column_accumulator(np.zeros(shape=(0))),
                "dRJ12AK8": processor.column_accumulator(np.zeros(shape=(0))),
                "mT": processor.column_accumulator(np.zeros(shape=(0))),
                "METrHT_pt30": processor.column_accumulator(np.zeros(shape=(0))),
                "METrST_pt30": processor.column_accumulator(np.zeros(shape=(0))),
                "j1Pt": processor.column_accumulator(np.zeros(shape=(0))),
                "j1Eta": processor.column_accumulator(np.zeros(shape=(0))),
                "j1Phi": processor.column_accumulator(np.zeros(shape=(0))),
                "j1Axismajor": processor.column_accumulator(np.zeros(shape=(0))),
                "j1Axisminor": processor.column_accumulator(np.zeros(shape=(0))),
                "j1PtD": processor.column_accumulator(np.zeros(shape=(0))),
                "dPhij1MET": processor.column_accumulator(np.zeros(shape=(0))),
                "j2Pt": processor.column_accumulator(np.zeros(shape=(0))),
                "j2Eta": processor.column_accumulator(np.zeros(shape=(0))),
                "j2Phi": processor.column_accumulator(np.zeros(shape=(0))),
                "j2Axismajor": processor.column_accumulator(np.zeros(shape=(0))),
                "j2Axisminor": processor.column_accumulator(np.zeros(shape=(0))),
                "j2PtD": processor.column_accumulator(np.zeros(shape=(0))),
                "dPhij2MET": processor.column_accumulator(np.zeros(shape=(0))),
                "dPhij1rdPhij2": processor.column_accumulator(np.zeros(shape=(0))),
                "j1PtAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "j1EtaAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "j1PhiAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "j1AxismajorAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "j1AxisminorAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "j1GirthAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "j1PtDAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "j1Tau1AK8": processor.column_accumulator(np.zeros(shape=(0))),
                "j1Tau2AK8": processor.column_accumulator(np.zeros(shape=(0))),
                "j1Tau3AK8": processor.column_accumulator(np.zeros(shape=(0))),
                "j1Tau21AK8": processor.column_accumulator(np.zeros(shape=(0))),
                "j1Tau32AK8": processor.column_accumulator(np.zeros(shape=(0))),
                "j1SoftDropMassAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "dPhij1METAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "j2PtAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "j2EtaAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "j2PhiAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "j2AxismajorAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "j2AxisminorAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "j2GirthAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "j2PtDAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "j2Tau1AK8": processor.column_accumulator(np.zeros(shape=(0))),
                "j2Tau2AK8": processor.column_accumulator(np.zeros(shape=(0))),
                "j2Tau3AK8": processor.column_accumulator(np.zeros(shape=(0))),
                "j2Tau21AK8": processor.column_accumulator(np.zeros(shape=(0))),
                "j2Tau32AK8": processor.column_accumulator(np.zeros(shape=(0))),
                "j2SoftDropMassAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "dPhij2METAK8": processor.column_accumulator(np.zeros(shape=(0))),
                "dPhij1rdPhij2AK8": processor.column_accumulator(np.zeros(shape=(0)))
                })
        @property
        def accumulator(self):
                return self._accumulator

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
                metPhi_noCut = df['METPhi']
                cuts = {
                            "" : np.ones(len(df["Weight"]),dtype=bool)
                        }
# add cutflow (see coffea cutflow)
# go through the paper to get cutflows, get the cuts for t-channel from the paper
# try to generate tchannel production signals: don't use madgraph flag, use t-channel, don't need LHE

                output = self.accumulator.identity()

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
                    jw = utl.awkwardReshape(jets,evtw)
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
                        ht = ak.sum(jets.pt,axis=1)
                        st = ht + met
                        metrht = utl.divide_vec(met,ht)
                        metrst = utl.divide_vec(met,st)
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
                        jetAK8pT = fjets.pt
                        jetAK8Phi = fjets.phi
                        jetAK8Eta = fjets.eta
                        jetAK8M = fjets.mass
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
                        # Getting subset of variables based on number of AK8 jets
                        # calculating variables
                        output['evtw'] += col_accumulator(evtw)
                        output['jw'] += col_accumulator(ak.flatten(jw))
                        output['fjw'] += col_accumulator(ak.flatten(fjw))
                        output['njets'] += col_accumulator(jets.counts)
                        output['njetsAK8'] += col_accumulator(fjets.counts)
                        output['nb'] += col_accumulator(nBJets)
                        output['nl'] += col_accumulator((electrons.counts + muons.counts))
                        output['ht'] += col_accumulator(ht)
                        output['st'] += col_accumulator(st)
                        output['met'] += col_accumulator(met)
                        output['madHT'] += col_accumulator(madHT)
                        output['jPt'] += col_accumulator(jets.pt.flatten())
                        output['jEta'] += col_accumulator(jetEta.flatten())
                        output['jPhi'] += col_accumulator(jetPhi.flatten())
                        output['jAxismajor'] += col_accumulator(jets.axismajor.flatten())
                        output['jAxisminor'] += col_accumulator(jets.axisminor.flatten())
                        output['jPtD'] += col_accumulator(jets.ptD.flatten())
                        output['dPhiMinjMET'] += col_accumulator(dPhiMinj)
                        output['dEtaj12'] += col_accumulator(dEtaj12)
                        output['dRJ12'] += col_accumulator(deltaR12j)
                        output['jPtAK8'] += col_accumulator(fjets.pt.flatten())
                        output['jEtaAK8'] += col_accumulator(jetAK8Eta.flatten())
                        output['jPhiAK8'] += col_accumulator(jetAK8Phi.flatten())
                        output['jAxismajorAK8'] += col_accumulator(fjets.axismajor.flatten())
                        output['jAxisminorAK8'] += col_accumulator(fjets.axisminor.flatten())
                        output['jGirthAK8'] += col_accumulator(fjets.girth.flatten())
                        output['jPtDAK8'] += col_accumulator(fjets.ptD.flatten())
                        output['jTau1AK8'] += col_accumulator(tau1.flatten())
                        output['jTau2AK8'] += col_accumulator(tau2.flatten())
                        output['jTau3AK8'] += col_accumulator(tau3.flatten())
                        output['jTau21AK8'] += col_accumulator(J_tau21)
                        output['jTau32AK8'] += col_accumulator(J_tau32)
                        output['jSoftDropMassAK8'] += col_accumulator(fjets.softDropMass.flatten())
                        output['dPhiMinjMETAK8'] += col_accumulator(dPhiMinjAK8)
                        output['dEtaj12AK8'] += col_accumulator(dEtaj12AK8)
                        output['dRJ12AK8'] += col_accumulator(deltaR12jAK8)
                        output['mT'] += col_accumulator(mtAK8)
                        output['METrHT_pt30'] += col_accumulator(metrht)
                        output['METrST_pt30'] += col_accumulator(metrst)
                        output['j1Pt'] += col_accumulator(utl.jetVar_vec(jets.pt,0))
                        output['j1Eta'] += col_accumulator(j1_eta)
                        output['j1Phi'] += col_accumulator(j1_phi)
                        output['j1Axismajor'] += col_accumulator(utl.jetVar_vec(jets.axismajor,0))
                        output['j1Axisminor'] += col_accumulator(utl.jetVar_vec(jets.axisminor,0))
                        output['j1PtD'] += col_accumulator(utl.jetVar_vec(jets.ptD,0))
                        output['dPhij1MET'] += col_accumulator(dPhij1)
                        output['j2Pt'] += col_accumulator(utl.jetVar_vec(jets.pt,1))
                        output['j2Eta'] += col_accumulator(j2_eta)
                        output['j2Phi'] += col_accumulator(j2_phi)
                        output['j2Axismajor'] += col_accumulator(utl.jetVar_vec(jets.axismajor,1))
                        output['j2Axisminor'] += col_accumulator(utl.jetVar_vec(jets.axisminor,1))
                        output['j2PtD'] += col_accumulator(utl.jetVar_vec(jets.ptD,1))
                        output['dPhij2MET'] += col_accumulator(dPhij2)
                        output['dPhij1rdPhij2'] += col_accumulator(dPhij1rdPhij2)
                        output['j1PtAK8'] += col_accumulator(utl.jetVar_vec(fjets.pt,0))
                        output['j1EtaAK8'] += col_accumulator(j1_etaAK8)
                        output['j1PhiAK8'] += col_accumulator(j1_phiAK8)
                        output['j1AxismajorAK8'] += col_accumulator(utl.jetVar_vec(fjets.axismajor,0))
                        output['j1AxisminorAK8'] += col_accumulator(utl.jetVar_vec(fjets.axisminor,0))
                        output['j1GirthAK8'] += col_accumulator(utl.jetVar_vec(fjets.girth,0))
                        output['j1PtDAK8'] += col_accumulator(utl.jetVar_vec(fjets.ptD,0))
                        output['j1Tau1AK8'] += col_accumulator(utl.jetVar_vec(fjets.tau1,0))
                        output['j1Tau2AK8'] += col_accumulator(utl.jetVar_vec(fjets.tau2,0))
                        output['j1Tau3AK8'] += col_accumulator(utl.jetVar_vec(fjets.tau3,0))
                        output['j1Tau21AK8'] += col_accumulator(J1_tau21)
                        output['j1Tau32AK8'] += col_accumulator(J1_tau32)
                        output['j1SoftDropMassAK8'] += col_accumulator(utl.jetVar_vec(fjets.softDropMass,0))
                        output['dPhij1METAK8'] += col_accumulator(dPhij1AK8)
                        output['j2PtAK8'] += col_accumulator(utl.jetVar_vec(fjets.pt,1))
                        output['j2EtaAK8'] += col_accumulator(j2_etaAK8)
                        output['j2PhiAK8'] += col_accumulator(j2_phiAK8)
                        output['j2AxismajorAK8'] += col_accumulator(utl.jetVar_vec(fjets.axismajor,1))
                        output['j2AxisminorAK8'] += col_accumulator(utl.jetVar_vec(fjets.axisminor,1))
                        output['j2GirthAK8'] += col_accumulator(utl.jetVar_vec(fjets.girth,1))
                        output['j2PtDAK8'] += col_accumulator(utl.jetVar_vec(fjets.ptD,1))
                        output['j2Tau1AK8'] += col_accumulator(utl.jetVar_vec(fjets.tau1,1))
                        output['j2Tau2AK8'] += col_accumulator(utl.jetVar_vec(fjets.tau2,1))
                        output['j2Tau3AK8'] += col_accumulator(utl.jetVar_vec(fjets.tau3,1))
                        output['j2Tau21AK8'] += col_accumulator(J2_tau21)
                        output['j2Tau32AK8'] += col_accumulator(J2_tau32)
                        output['j2SoftDropMassAK8'] += col_accumulator(utl.jetVar_vec(fjets.softDropMass,1))
                        output['dPhij2METAK8'] += col_accumulator(dPhij2AK8)
                        output['dPhij1rdPhij2AK8'] += col_accumulator(dPhij1rdPhij2AK8)
                #
                    # if len(ak.flatten(tPassedList)) > 0:
                    #     output['h_trigger'].fill(trigger=ak.flatten(tPassedList),weight=np.ones(len(ak.flatten(tPassedList))))
                return output

        def postprocess(self, accumulator):
                return accumulator
