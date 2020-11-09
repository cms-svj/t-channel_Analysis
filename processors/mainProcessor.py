from coffea import hist, processor
from coffea.analysis_objects import JaggedCandidateArray

class MainProcessor(processor.ProcessorABC):
    def __init__(self):
        # dataset_axis = hist.Cat("dataset", "Primary dataset")
        # pt_axis = hist.Bin("pt", r"$p_{T}$ [GeV]", 40, 0, 3500)
        pt_axis = hist.Bin("pt", r"$p_{T}$ [GeV]", 40, 0, 3500)
        eta_axis = hist.Bin("eta", r"$\eta$", 15, -7, 7)
        phi_axis = hist.Bin("phi", r"$\phi$", 15, -3.2, 3.2)
        axismajor_axis = hist.Bin("axismajor", r"$\sigma_{major}(j)$", 40, 0, 0.5)
        axisminor_axis = hist.Bin("axisminor", r"$\sigma_{minor}(j)$", 40, 0, 0.3)
        girth_axis = hist.Bin("girth","girth(j)",40, 0, 0.5)
        ptD_axis = hist.Bin("ptD", "ptD", 40, 0, 1.2)
        tau1_axis = hist.Bin("tau1",r"$\tau_{1}(j)$",40, 0, 0.8)
        tau2_axis = hist.Bin("tau2",r"$\tau_{2}(j)$",40, 0, 0.65)
        tau3_axis = hist.Bin("tau3",r"$\tau_{3}(j)$",40, 0, 0.35)
        tau21_axis = hist.Bin("tau21",r"$\tau_{21}(j)$",40, 0, 1.3)
        tau32_axis = hist.Bin("tau32",r"$\tau_{32}(j)$",40, 0, 1.3)
        softDropMass_axis = hist.Bin("softDropMass",r"$m_{SD}(j)$",40,0,200)
        MET_axis = hist.Bin("MET", "MET [GeV]", 20, 0, 2000)

        self._accumulator = processor.dict_accumulator({
            # 'jtpt':hist.Hist("Counts", dataset_axis, pt_axis),
            # 'jteta':hist.Hist("Counts",dataset_axis,eta_axis),
            'AK8Jets_pt':hist.Hist("AK8Jets_pt", pt_axis),
            'AK8Jets_eta':hist.Hist("AK8Jets_eta",eta_axis),
            'AK8Jets_phi':hist.Hist("AK8Jets_phi",phi_axis),
            'AK8Jets_axismajor':hist.Hist("AK8Jets_axismajor", axismajor_axis),
            'AK8Jets_axisminor':hist.Hist("AK8Jets_axisminor",axisminor_axis),
            'AK8Jets_girth':hist.Hist("AK8Jets_girth",girth_axis),
            'AK8Jets_ptD':hist.Hist("AK8Jets_ptD",ptD_axis),
            'AK8Jets_tau1':hist.Hist("AK8Jets_tau1",tau1_axis),
            'AK8Jets_tau2':hist.Hist("AK8Jets_tau2",tau2_axis),
            'AK8Jets_tau3':hist.Hist("AK8Jets_tau3", tau3_axis),
            'AK8Jets_tau21':hist.Hist("AK8Jets_tau21",tau21_axis),
            'AK8Jets_tau32':hist.Hist("AK8Jets_tau32",tau32_axis),
            'AK8Jets_softDropMass':hist.Hist("AK8Jets_softDropMass",softDropMass_axis),
            'AK4Jets_pt':hist.Hist("AK4Jets_pt", pt_axis),
            'AK4Jets_eta':hist.Hist("AK4Jets_eta",eta_axis),
            'AK4Jets_phi':hist.Hist("AK4Jets_phi",phi_axis),
            'AK4Jets_axismajor':hist.Hist("AK4Jets_axismajor", axismajor_axis),
            'AK4Jets_axisminor':hist.Hist("AK4Jets_axisminor",axisminor_axis),
            'MET':hist.Hist("MET",MET_axis),
            'cutflow': processor.defaultdict_accumulator(int),
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

        # # Good AK8 Jets Cut
        # ptAK8cut = (fjets.pt > 200 and abs(fjets.eta) < 2.4 )
        # fjets = fjets[ptAK8cut]
        #
        # # Good AK4 Jets Cut
        # ptAK4cut = (jets.pt > 30 and abs(jets.eta) < 2.4 )
        # jets = jets[ptAK4cut]

        twofjets = (fjets.counts >= 2)
        output['cutflow']['two fjets'] += twofjets.sum()

        # difjets = fjets[twofjets]
        # difjets_pt200 = difjets[ptcut]
        # output['jtpt'].fill(dataset=dataset, pt=fjets.pt.flatten())
        # output['jteta'].fill(dataset=dataset, eta=fjets.eta.flatten())

        output['AK8Jets_pt'].fill(pt=fjets.pt.flatten())
        output['AK8Jets_eta'].fill(eta=fjets.eta.flatten())
        output['AK8Jets_phi'].fill(phi=fjets.phi.flatten())
        output['AK8Jets_axismajor'].fill(axismajor=fjets.axismajor.flatten())
        output['AK8Jets_axisminor'].fill(axisminor=fjets.axisminor.flatten())
        output['AK8Jets_girth'].fill(girth=fjets.girth.flatten())
        output['AK8Jets_ptD'].fill(ptD=fjets.ptD.flatten())
        output['AK8Jets_tau1'].fill(tau1=fjets.tau1.flatten())
        output['AK8Jets_tau2'].fill(tau2=fjets.tau2.flatten())
        output['AK8Jets_tau3'].fill(tau3=fjets.tau3.flatten())
        output['AK8Jets_tau21'].fill(tau21=fjets.tau2.flatten()/fjets.tau1.flatten())
        output['AK8Jets_tau32'].fill(tau32=fjets.tau3.flatten()/fjets.tau2.flatten())
        output['AK8Jets_softDropMass'].fill(softDropMass=fjets.softDropMass.flatten())
        output['AK4Jets_pt'].fill(pt=jets.pt.flatten())
        output['AK4Jets_eta'].fill(eta=jets.eta.flatten())
        output['AK4Jets_phi'].fill(phi=jets.phi.flatten())
        output['AK4Jets_axismajor'].fill(axismajor=jets.axismajor.flatten())
        output['AK4Jets_axisminor'].fill(axisminor=jets.axisminor.flatten())
        output['MET'].fill(MET=df['MET'])
        return output

    def postprocess(self, accumulator):
        return accumulator

