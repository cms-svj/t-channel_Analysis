from coffea import hist, processor
import numpy as np
import awkward1 as ak
import awkward
from utils import utility as utl
import utils.objects as ob
from utils import baseline as bl
from itertools import combinations
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
            for name in cuts.keys():
                varDict[name] =    processor.column_accumulator(np.zeros(shape=(0)))
            self._accumulator = processor.dict_accumulator(varDict)
            self.setupNPArr = True

        def process(self, df):
                ## objects used for cuts
                inpObj_noCut = ob.inpObj(df,self.scaleFactor)
                # Our preselection
                cuts = bl.cutList(df,inpObj_noCut)

                if self.setupNPArr is None:
                    self.setupNPArray(cuts,variables)
                output = self.accumulator.identity()

                # run cut loop
                for name,cut in cuts.items():
                    # defining objects
                    inpObj = {}
                    for key,item in inpObj_noCut.items():
                        inpObj[key] = item[cut]

                    if len(inpObj['evtw']) > 0:
                        if name == "_npz":
                            varValDict = utl.varGetter(inpObj)
                            for varName,varDetail in varValDict.items():
                                if variables[varName][4] == 1:
                                    output['{}'.format(varName)] += col_accumulator(varDetail[0])
                                elif variables[varName][4] == 2:
                                    output['{}'.format(varName)] += col_accumulator(np.repeat(ak.to_numpy(varDetail[0]),ak.to_numpy(varValDict["njetsAK8"][0])))
                return output

        def postprocess(self, accumulator):
                return accumulator
