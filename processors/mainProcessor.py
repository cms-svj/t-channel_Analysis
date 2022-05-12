from coffea import hist, processor
import numpy as np
import awkward as ak
from utils import utility as utl
from utils import baseline as bl
from utils.variables import variables
from utils.runNeuralNetwork import runNN

class MainProcessor(processor.ProcessorABC):
        def __init__(self,dataset,sf,model,varSet,normMean,normStd,jNVar):
                self._accumulator = processor.dict_accumulator({})
                self.setupHistos = None
                self.dataset = dataset
                self.scaleFactor = sf
                self.model = model
                self.varSet = varSet
                self.normMean = normMean
                self.normStd = normStd
                self.jNVar = jNVar
        @property
        def accumulator(self):
                return self._accumulator

        def setupHistogram(self,cuts):
                histograms = {}
                for cutName,cut in cuts.items():
                    for histName,histDetail in variables(self.jNVar).items():
                        histograms['h_{}{}'.format(histName,cutName)] = hist.Hist('h_{}{}'.format(histName,cutName), hist.Bin("val", histDetail[0], histDetail[1], histDetail[2], histDetail[3]))
                self._accumulator = processor.dict_accumulator(histograms)
                self.setupHistos = True

        def process(self, events):
                # cut loop
                ## objects used for cuts
                vars_noCut = utl.varGetter(self.dataset,events,self.scaleFactor,self.jNVar)
                runNN(self.model,vars_noCut,self.varSet,self.normMean,self.normStd)
                # Our preselection
                cuts = bl.cutList(self.dataset,events,vars_noCut,SVJCut=False)

                # setup histograms
                if self.setupHistos is None:
                    self.setupHistogram(cuts)
                output = self.accumulator.identity()

                # run cut loop
                for cutName,cut in cuts.items():
                    # defining objects
                    weight = vars_noCut["evtw"][cut]
                    jweight = ak.flatten(vars_noCut["jw"][cut])
                    fjweight = ak.flatten(vars_noCut["fjw"][cut])
                    eweight = ak.flatten(vars_noCut["ew"][cut])
                    mweight = ak.flatten(vars_noCut["mw"][cut])
                    nimweight = ak.flatten(vars_noCut["nimw"][cut])
                    pred1_evtw = vars_noCut["pred1_evtw"][cut]
                    pred2_evtw = vars_noCut["pred2_evtw"][cut]
                    if len(events) > 0:
                        ## filling histograms
                        for varName,varDetail in variables(self.jNVar).items():
                            hIn = vars_noCut[varName][cut]
                            hW = weight
                            wKey = varDetail[6]
                            # properly flatten certain inputs
                            if varDetail[5] >= 1:
                                hIn = ak.flatten(hIn)
                            # make sure the correct weights are applied
                            if wKey == "jw":
                                hW = jweight
                            elif wKey == "fjw":
                                hW = fjweight
                            elif wKey == "ew":
                                hW = eweight
                            elif wKey == "mw":
                                hW = mweight
                            elif wKey == "nimw":    
                                hW = nimweight
                            elif wKey == "pred1_evtw":
                                hW = pred1_evtw
                            elif wKey == "pred2_evtw":
                                hW = pred2_evtw
                            elif wKey == "w1":
                                hW = np.ones(len(weight))
                            if len(hIn) > 0:
                                output['h_{}{}'.format(varName,cutName)].fill(val=hIn,weight=hW)
                        ## filling histograms by jet category
                        # jetCat = ak.flatten(inpObj["JetsAK8_hvCategory"]) == 17 # 9 = QdM, 17 = QsM
                        # for varName,varDetail in varValDict.items():
                        #     if variables(self.jNVar)[varName][4] == 1:
                        #         output['h_{}{}'.format(varName,cutName)].fill(val=varDetail[0][jetCat],weight=varDetail[1][jetCat])
                        #     else:
                        #         output['h_{}{}'.format(varName,cutName)].fill(val=varDetail[0],weight=varDetail[1])
                return output

        def postprocess(self, accumulator):
                return accumulator
