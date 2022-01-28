from coffea import hist, processor
import numpy as np
import awkward1 as ak
import awkward
from utils import utility as utl
import utils.objects as ob
from utils import baseline as bl
from utils.variables import variables
from itertools import combinations
from utils.runNeuralNetwork import runNN

class MainProcessor(processor.ProcessorABC):
        def __init__(self,sf,model,varSet,normMean,normStd):
                self._accumulator = processor.dict_accumulator({})
                self.setupHistos = None
                self.scaleFactor = sf
                self.model = model
                self.varSet = varSet
                self.normMean = normMean
                self.normStd = normStd
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
                vars_noCut = utl.varGetter(df,self.scaleFactor)
                runNN(self.model,vars_noCut,self.varSet,self.normMean,self.normStd)
                # Our preselection
                cuts = bl.cutList(df,vars_noCut)

                # setup histograms
                if self.setupHistos is None:
                    self.setupHistogram(cuts)
                output = self.accumulator.identity()

                # run cut loop
                for cutName,cut in cuts.items():
                    # defining objects
                    weight = vars_noCut["evtw"][0][cut]
                    jweight = ak.flatten(vars_noCut["jw"][0][cut])
                    fjweight = ak.flatten(vars_noCut["fjw"][0][cut])
                    if len(weight) > 0:
                        ## filling histograms
                        for varName,varDetail in variables.items():
                            hIn = vars_noCut[varName][0][cut]
                            hW = weight
                            wKey = vars_noCut[varName][1]
                            # properly flatten certain inputs
                            if varDetail[5] == 1:
                                hIn = hIn.flatten()
                            elif varDetail[5] == 2:
                                hIn = ak.flatten(hIn)
                            # make sure the correct weights are applied
                            if wKey == "jw":
                                hW = jweight
                            elif wKey == "fjw":
                                hW = fjweight
                            elif wKey == "w1":
                                hW = np.ones(len(weight))
                            if len(hIn) > 0:
                                output['h_{}{}'.format(varName,cutName)].fill(val=hIn,weight=hW)
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
