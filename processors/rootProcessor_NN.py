from coffea import hist, processor
import numpy as np
import awkward as ak
from utils import utility as utl
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

        def setupNPArray(self,variables):
            varDict = {}
            for v,d in variables.items():
                if d.npzInfo > 0:
                    varDict[v] = processor.column_accumulator(np.zeros(shape=(0)))
            self._accumulator = processor.dict_accumulator(varDict)
            self.setupNPArr = True

        def process(self, events):
                ## objects used for cuts
                vars_noCut = utl.baselineVar(self.dataset,events,self.scaleFactor)
                # Our preselection
                cuts = bl.cutList(self.dataset,events,vars_noCut,SVJCut=False)
                
                if self.setupNPArr is None:
                    self.setupNPArray(variables())
                output = self.accumulator.identity()

                # run cut loop
                cut = cuts["_qual_trg_st_1PJ"]
                if (len(events) > 0) and (np.any(cut)):
                    eventsCut = events[cut]
                    utl.varGetter(self.dataset,eventsCut,vars_noCut,cut,False)
                    if (len(eventsCut) > 0):
                        for varName,varDetail in variables().items():
                            # only store jetAK8 variables
                            if varDetail.npzInfo == 1:
                                hIn = vars_noCut[varName]
                                # properly flatten certain inputs
                                hIn = ak.flatten(hIn)
                                output['{}'.format(varName)] += col_accumulator(hIn)

                return output

        def postprocess(self, accumulator):
                return accumulator
