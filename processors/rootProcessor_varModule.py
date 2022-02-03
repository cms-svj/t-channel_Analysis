from coffea import hist, processor
import numpy as np
import awkward1 as ak
import awkward
from utils import utility as utl
import utils.objects as ob
from utils import baseline as bl
from utils.variables import variables

def col_accumulator(a):
    return processor.column_accumulator(np.array(a))

class MainProcessor(processor.ProcessorABC):
        def __init__(self,sf):
            self._accumulator = processor.dict_accumulator({})
            self.setupNPArr = None
            self.scaleFactor = sf
        @property
        def accumulator(self):
                return self._accumulator

        def setupNPArray(self,cuts,variables):
            varDict = {}
            for v,d in variables.items():
                if d[4] > 0:
                    varDict[v] = processor.column_accumulator(np.zeros(shape=(0)))
            for cutName in cuts.keys():
                varDict[cutName] =    processor.column_accumulator(np.zeros(shape=(0)))
            self._accumulator = processor.dict_accumulator(varDict)
            self.setupNPArr = True

        def process(self, df):
                ## objects used for cuts
                vars_noCut = utl.varGetter(df,self.scaleFactor)
                # Our preselection
                cuts = bl.cutList(df,vars_noCut,SVJCut=False)

                if self.setupNPArr is None:
                    self.setupNPArray(cuts,variables)
                output = self.accumulator.identity()

                # run cut loop
                cut = cuts["_npz"]
                weight = vars_noCut["evtw"][0][cut]
                if len(weight) > 0:
                    for varName,varDetail in variables.items():
                        # only store jetAK8 variables
                        if varDetail[4] == 1:
                            hIn = vars_noCut[varName][0][cut]
                            # properly flatten certain inputs
                            if varDetail[5] == 1:
                                hIn = hIn.flatten()
                            elif varDetail[5] == 2:
                                hIn = ak.flatten(hIn)
                            output['{}'.format(varName)] += col_accumulator(hIn)

                return output

        def postprocess(self, accumulator):
                return accumulator
