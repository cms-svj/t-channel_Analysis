from coffea import hist, processor
import numpy as np
import awkward1 as ak
import awkward
from utils import utility as utl
import utils.objects as ob
from utils import baseline as bl
from utils.variables import variables
from itertools import combinations

class MainProcessor(processor.ProcessorABC):
        def __init__(self,sf):
                self._accumulator = processor.dict_accumulator({})
                self.setupHistos = None
                self.scaleFactor = sf
        @property
        def accumulator(self):
                return self._accumulator

        def setupHistogram(self,cuts):
                histograms = {}
                for cutName,cut in cuts.items():
                    for histName,histDetail in variables.items():
                        histograms['h_{}{}'.format(histName,cutName)] = hist.Hist('h_{}{}'.format(histName,cutName), hist.Bin("val", histDetail[0], histDetail[1], histDetail[2], histDetail[3]))
                self._accumulator = processor.dict_accumulator(histograms)
                self.setupHistos = True

        def process(self, df):
                # cut loop
                ## objects used for cuts
                inpObj_noCut = ob.inpObj(df,self.scaleFactor)
                # Our preselection
                cuts = bl.cutList(df,inpObj_noCut)

                # setup histograms
                if self.setupHistos is None:
                    self.setupHistogram(cuts)
                output = self.accumulator.identity()

                # run cut loop
                for cutName,cut in cuts.items():
                    # defining objects
                    inpObj = {}
                    for key,item in inpObj_noCut.items():
                        inpObj[key] = item[cut]
                    if len(inpObj["evtw"]) > 0:
                        varValDict = utl.varGetter(inpObj)
                        # filling histograms
                        for varName,varDetail in varValDict.items():
                            output['h_{}{}'.format(varName,cutName)].fill(val=varDetail[0],weight=varDetail[1])
                        ## filling histograms by jet category
                        # jetCat = ak.flatten(inpObj["JetsAK8_hvCategory"]) == 17 # 9 = QdM, 17 = QsM
                        # for varName,varDetail in varValDict.items():
                        #     if variables[varName][4] == 1:
                        #         output['h_{}{}'.format(varName,cutName)].fill(val=varDetail[0][jetCat],weight=varDetail[1][jetCat])
                        #     else:
                        #         output['h_{}{}'.format(varName,cutName)].fill(val=varDetail[0],weight=varDetail[1])
                return output

        def postprocess(self, accumulator):
                return accumulator
