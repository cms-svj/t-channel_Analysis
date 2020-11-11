from coffea import hist, processor
from coffea.analysis_objects import JaggedCandidateArray
import uproot
import numpy as np
import awkward as ak
import sys
from matplotlib import pyplot as plt
import json
import glob
# e.g. sample collection = "QCD"
# sample = sys.argv[1]

class FancyDimuonProcessor(processor.ProcessorABC):
    def __init__(self):
        # dataset_axis = hist.Cat("dataset", "Primary dataset")
        # pt_axis = hist.Bin("pt", r"$p_{T}$ [GeV]", 40, 0, 3500)
        njet_axis                   = hist.Bin("njets",                 "Number of Jets",                                   20,     0.0,    20.0)
        ht_axis                     = hist.Bin("ht",                    r"$H_{T}$ (GeV)",                                   500,    0.0,    5000.0)
        st_axis                     = hist.Bin("st",                    r"$S_{T}$ (GeV)",                                   500,    0.0,    5000.0)
        MET_axis                    = hist.Bin("MET",                   "MET [GeV]",                                        500,    0.0,    5000.0)
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
        dPhiJMET_axis               = hist.Bin("dPhiJMET",             r"$\Delta\phi(J_{1},MET)$",                         100,    0.0,    4.0)
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
            'h_njets':                  hist.Hist("h_njets",                njet_axis),
            'h_njetsAK8':               hist.Hist("h_njetsAK8",             njet_axis),
            'h_ht':                     hist.Hist("h_ht",                   ht_axis),
            'h_st':                     hist.Hist("h_st",                   st_axis),
            'h_met':                    hist.Hist("h_met",                  MET_axis),
            'h_jPt':                    hist.Hist("h_jPt",                  pt_axis),
            'h_jEta':                   hist.Hist("h_jEta",                 eta_axis),
            'h_jPhi':                   hist.Hist("h_jPhi",                 phi_axis),
            'h_jAxismajor':             hist.Hist("h_jAxismajor",           axismajor_axis),
            'h_jAxisminor':             hist.Hist("h_jAxisminor",           axisminor_axis),
            'h_jPtD':                   hist.Hist("h_jPtD",                 ptD_axis),
            'h_jPtAK8':                 hist.Hist("h_jPtAK8",               pt_axis),
            'h_jEtaAK8':                hist.Hist("h_jEtaAK8",              eta_axis),
            'h_jPhiAK8':                hist.Hist("h_jPhiAK8",              phi_axis),
            'h_jAxismajorAK8':          hist.Hist("h_jAxismajorAK8",        axismajor_axis),
            'h_jAxisminorAK8':          hist.Hist("h_jAxisminorAK8",        axisminor_axis),
            'h_jGirthAK8':              hist.Hist("h_jGirthAK8",            girth_axis),
            'h_jPtDAK8':                hist.Hist("h_jPtDAK8",              ptD_axis),
            'h_jTau1AK8':               hist.Hist("h_jTau1AK8",             tau1_axis),
            'h_jTau2AK8':               hist.Hist("h_jTau2AK8",             tau2_axis),
            'h_jTau3AK8':               hist.Hist("h_jTau3AK8",             tau3_axis),
            'h_jTau21AK8':              hist.Hist("h_jTau21AK8",            tau21_axis),
            'h_jTau32AK8':              hist.Hist("h_jTau32AK8",            tau32_axis),
            'h_jSoftDropMassAK8':       hist.Hist("h_jSoftDropMassAK8",     softDropMass_axis),
            'h_weight':                 hist.Hist("h_weight",               weight_axis),
            'h_dEtaJ12':                hist.Hist("h_dEtaJ12",              dEtaJ12_axis),
            'h_dRJ12':                  hist.Hist("h_dRJ12",                dRJ12_axis),
            'h_dPhiJ1MET':              hist.Hist("h_dPhiJ1MET",            dPhiJMET_axis),
            'h_dPhiJ2MET':              hist.Hist("h_dPhiJ2MET",            dPhiJMET_axis),
            'h_dPhiMinJMET':            hist.Hist("h_dPhiMinJMET",          dPhiJMET_axis),
            'h_dPhiJ1METrdPhiJ2MET':    hist.Hist("h_dPhiJ1METrdPhiJ2MET",  dPhiJ1METrdPhiJ2MET_axis),
            'h_mjjM':                   hist.Hist("h_mjjM",                 mjjM_axis),
            'h_mjjPt':                  hist.Hist("h_mjjPt",                mjjPt_axis),
            'h_mjjEta':                 hist.Hist("h_mjjEta",               mjjEta_axis),
            'h_mT':                     hist.Hist("h_mT",                   mT_axis),
            'h_METrHT_pt30':            hist.Hist("h_METrHT_pt30",          METrHT_pt30_axis),
            'h_METrST_pt30':            hist.Hist("h_METrST_pt30",          METrST_pt30_axis),
            'cutflow':                  processor.defaultdict_accumulator(int),
        })

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, df):
        output = self.accumulator.identity()

        dataset = df['dataset']

        jets = JaggedCandidateArray.candidatesfromcounts(
            df['Jets'].counts,
            px=df['Jets'].fP.fX.flatten(),
            py=df['Jets'].fP.fY.flatten(),
            pz=df['Jets'].fP.fZ.flatten(),
            energy=df['Jets'].fE.flatten(),
            axismajor=df['Jets_axismajor'].flatten(),
            axisminor=df['Jets_axisminor'].flatten(),
            ptD=df['Jets_ptD'].flatten()
            )

        fjets = JaggedCandidateArray.candidatesfromcounts(
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

        output['cutflow']['all events'] += fjets.size

        # important variables
        MET_ = df['MET']
        MT_AK8_ = df['MT_AK8']
        w = df["Weight"]

        # # Good AK8 Jets Cut
        AK8qualityCut = (fjets.pt > 200) & (abs(fjets.eta) < 2.4)
        fjets = fjets[AK8qualityCut]

        # # Good AK4 Jets Cut
        AK4qualityCut = (jets.pt > 30) & (abs(jets.eta) < 2.4)
        jets = jets[AK4qualityCut]


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

        output['h_njets'].fill(njets=jets.counts.flatten(),weight=w.flatten())
        output['h_njetsAK8'].fill(njets=fjets.counts.flatten())
        output['h_jPtAK8'].fill(pt=fjets.pt.flatten())
        output['h_jEtaAK8'].fill(eta=fjets.eta.flatten())
        output['h_jPhiAK8'].fill(phi=fjets.phi.flatten())
        output['h_jAxismajorAK8'].fill(axismajor=fjets.axismajor.flatten())
        output['h_jAxisminorAK8'].fill(axisminor=fjets.axisminor.flatten())
        output['h_jGirthAK8'].fill(girth=fjets.girth.flatten())
        output['h_jPtDAK8'].fill(ptD=fjets.ptD.flatten())
        output['h_jTau1AK8'].fill(tau1=fjets.tau1.flatten())
        output['h_jTau2AK8'].fill(tau2=fjets.tau2.flatten())
        output['h_jTau3AK8'].fill(tau3=fjets.tau3.flatten())
        output['h_jTau21AK8'].fill(tau21=fjets.tau2.flatten()/fjets.tau1.flatten())
        output['h_jTau32AK8'].fill(tau32=fjets.tau3.flatten()/fjets.tau2.flatten())
        output['h_jSoftDropMassAK8'].fill(softDropMass=fjets.softDropMass.flatten())
        output['h_jPt'].fill(pt=jets.pt.flatten())
        output['h_jEta'].fill(eta=jets.eta.flatten())
        output['h_jPhi'].fill(phi=jets.phi.flatten())
        output['h_jAxismajor'].fill(axismajor=jets.axismajor.flatten())
        output['h_jAxisminor'].fill(axisminor=jets.axisminor.flatten())
        output['h_jPtD'].fill(ptD=jets.ptD.flatten())
        output['h_met'].fill(MET=MET_)
        return output

    def postprocess(self, accumulator):
        return accumulator

import time

tstart = time.time()

# getting dictionary of files
# f_ = sample.find("_")
# year = sample[:f_]
#
# detailKey = sample[f_:]
# if "mMed" in detailKey:
#     kind = "signals"
# else:
#     kind = "backgrounds"
#
# if "Incl" in detailKey:
#     ii = detailKey.find("Incl")
#     detailKey = detailKey[:ii] + "Tune"
#
# JSONDir = 'input/sampleJSONs/' + kind + "/" + year + "/"
# allfiles = glob.glob(JSONDir+"*.json")
# print (detailKey)
# for file in allfiles:
#     print (file)
#     if detailKey in file:
#         inputSample = file
#
# fileset = json.load(open(inputSample ,'r'))
fileset = {"sample":["sample.root"]}

output = processor.run_uproot_job(
    fileset,
    treename='TreeMaker2/PreSelection',
    processor_instance=FancyDimuonProcessor(),
    executor=processor.futures_executor,
    executor_args={'workers': 6, 'flatten': False},
    chunksize=10000,
)

elapsed = time.time() - tstart
print(output)

# import matplotlib
# ax = hist.plot1d(output['jtpt'], overlay='dataset')
# ax.set_xlim(70,150)
# ax.set_ylim(0.01, 5000)
# ax.set_yscale("log")

# hnames = ['h_njets','h_njetsAK8','h_met','h_jSoftDropMassAK8']
# hnames = list(output.keys())[:-1]
hnames = []
for key in list(output.keys())[:-1]:
    if len(output[key].values()) > 0:
        hnames.append(key)

# export the histograms to root files
fout = uproot.recreate('test.root')
for h in hnames:
    fout[h] = hist.export1d(output[h])
fout.close()
