from coffea import hist, processor
import numpy as np
import awkward as ak
from utils import utility as utl
from utils import baseline as bl
from utils.variables import variables

def col_accumulator(a):
    return processor.column_accumulator(np.array(a))

class MainProcessor(processor.ProcessorABC):
        def __init__(self,dataset,sf):
            self._accumulator = processor.dict_accumulator({})
            self.dataset = dataset
            self.setupNPArr = None
            self.scaleFactor = sf
        @property
        def accumulator(self):
                return self._accumulator

        def setupNPArray(self,varNames):
            varDict = {}
            for v in varNames:
                varDict[v] = processor.column_accumulator(np.zeros(shape=(0)))
            self._accumulator = processor.dict_accumulator(varDict)
            self.setupNPArr = True

        def process(self, events):
                ## objects used for cuts
                baselineVars = utl.baselineVar(self.dataset,events,self.scaleFactor)

                # Our preselection
                cuts = bl.cutList(self.dataset,events,baselineVars,SVJCut=False)

                inputVars = [
                                "jCstPt",
                                "jCstEta",
                                "jCstPhi",
                                "jCstEnergy",
                                "jCstPdgId",
                                "jCstPtAK8",
                                "jCstEtaAK8",
                                "jCstPhiAK8",
                                "jCstEnergyAK8",
                                "jCstAxismajorAK8",
                                "jCstAxisminorAK8",
                                "jCstdoubleBDiscriminatorAK8",
                                "jCstTau1AK8",
                                "jCstTau2AK8",
                                "jCstTau3AK8",
                                "jCstNumBhadronsAK8",
                                "jCstNumChadronsAK8",
                                "jCstPtDAK8",
                                "jCstSoftDropMassAK8",
                                "jCsthvCategory",
                                "jCstWeightAK8",
                                "jCstEvtNum",
                                "jCstJNum"
                ]

                if self.setupNPArr is None:
                    self.setupNPArray(inputVars)
                output = self.accumulator.identity()

                # run cut loop
                cut = cuts["_qual_trg_st_1PJ"]
                # for c in cut:
                #     print(c)
                if (len(events) > 0) and (np.any(cut)):
                    eventsCut = events[cut]
                    if (len(eventsCut) > 0):
                        evtw = baselineVars["evtw"]
                        evtwCut = evtw[cut]
                        fjets = baselineVars["fjets"]
                        fjetsCut = fjets[cut]
                        jCst4vec_noCut,jCstVar_noCut = utl.jConstVarGetter(self.dataset,eventsCut,baselineVars,cut)
                        for varName in jCstVar_noCut.keys():
                            hIn = ak.flatten(jCstVar_noCut[varName][1])
                            finiteMask = np.isfinite(hIn)
                            output['{}'.format(varName)] += col_accumulator(hIn[finiteMask])
                        for varName in jCst4vec_noCut.keys():
                            hIn = ak.flatten(jCst4vec_noCut[varName])
                            output['{}'.format(varName)] += col_accumulator(hIn[finiteMask])

                return output

        def postprocess(self, accumulator):
                return accumulator
