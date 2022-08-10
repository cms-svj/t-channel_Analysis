from coffea import processor
import numpy as np
import awkward as ak
from utils import utilityML as utl
import utils.objects as ob
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
                vars_noCut,jCst4vec_noCut,jCstVar_noCut = utl.varGetter(self.dataset,events,self.scaleFactor)
                # Our preselection
                cuts = bl.cutList(self.dataset,events,vars_noCut,SVJCut=False)

                if self.setupNPArr is None:
                    self.setupNPArray(list(jCst4vec_noCut.keys()) + list(jCstVar_noCut.keys()))
                output = self.accumulator.identity()

                # run cut loop
                cut = cuts["_qual_trg_st"]
                if len(events) > 0:
                    for varName in jCstVar_noCut.keys():
                        hIn = ak.flatten(jCstVar_noCut[varName][1][cut])
                        finiteMask = np.isfinite(hIn)
                        output['{}'.format(varName)] += col_accumulator(hIn[finiteMask])
                    for varName in jCst4vec_noCut.keys():
                        hIn = ak.flatten(jCst4vec_noCut[varName][cut])
                        output['{}'.format(varName)] += col_accumulator(hIn[finiteMask])

                return output

        def postprocess(self, accumulator):
                return accumulator
