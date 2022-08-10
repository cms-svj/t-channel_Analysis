from coffea import processor
import hist as h
import numpy as np
import awkward as ak
from utils import utility as utl
from utils import baseline as bl
from utils.variables import variables
from utils.runNeuralNetwork import runNN
import uproot

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
                self.fakerateHisto = self.getHistoFromFile("fakerate.root", "jPt_Fakerate_SR;1") 

        @property
        def accumulator(self):
                return self._accumulator

        def getHistoFromFile(self, fName, hName):
                try:
                    f = uproot.open(fName)
                    h = f[hName]
                    return h
                except FileNotFoundError:
                    print("\n\n\tError: No such file or directory: '{}'".format(fName))                    
                    print("\tWill use default fakerate of 1.0 for each jet\n\n")
                    return None
                except uproot.exceptions.KeyInFileError:
                    print("\n\n\tError: Histogram '{}' not found in file '{}'".format(hName,fName))
                    print("\tWill use default fakerate of 1.0 for each jet\n\n")
                    return None                    

        def getSFEvaluator(self, rootFileName, histoName):
                ext = lookup_tools.extractor()
                ext.add_weight_sets(["{} {} {}".format(histoName, histoName, rootFileName)])
                ext.finalize()
                evaluator = ext.make_evaluator()
                return evaluator[histoName]

        def setupHistogram(self,cuts):
                histograms = {}
                for cutName,cut in cuts.items():
                    for histName, histDetail in variables(self.jNVar).items():
                        if   histDetail.dim == 1:
                            histograms['h_{}{}'.format(histName,cutName)] = h.Hist(histDetail.xbins, storage="weight")
                        elif histDetail.dim == 2:                        
                            histograms['h_{}{}'.format(histName,cutName)] = h.Hist(histDetail.xbins, histDetail.ybins, storage="weight")

                self._accumulator = histograms
                self.setupHistos = True

        def process(self, events):
                # cut loop
                ## objects used for cuts
                vars_noCut = utl.varGetter(self.dataset,events,self.scaleFactor,self.jNVar)
                runNN(self.model,vars_noCut,self.varSet,self.normMean,self.normStd,self.fakerateHisto)
                # Our preselection
                cuts = bl.cutList(self.dataset,events,vars_noCut,SVJCut=False)

                # setup histograms
                if self.setupHistos is None:
                    self.setupHistogram(cuts)
                output = self.accumulator

                # run cut loop
                for cutName,cut in cuts.items():
                    # defining objects
                    weights = {
                            "evtw" : vars_noCut["evtw"][cut],
                            "jw"   : ak.flatten(vars_noCut["jw"][cut]),
                            "fjw"  : ak.flatten(vars_noCut["fjw"][cut]),
                            "ew"   : ak.flatten(vars_noCut["ew"][cut]),
                            "mw"   : ak.flatten(vars_noCut["mw"][cut]),
                            "nimw" : ak.flatten(vars_noCut["nimw"][cut]),
                            "svfjw" : ak.flatten(vars_noCut["svfjw"][cut]),
                            "pred1_evtw" : vars_noCut["pred1_evtw"][cut],
                            "pred2_evtw" : vars_noCut["pred2_evtw"][cut],
                            "pred3_evtw" : vars_noCut["pred3_evtw"][cut],
                            "pred4_evtw" : vars_noCut["pred4_evtw"][cut],
                    }

                    if len(events) > 0:
                        ## filling histograms
                        for histName, varDetail in variables(self.jNVar).items():
                            vX = vars_noCut[varDetail.varXName][cut]
                            vY = vars_noCut[varDetail.varYName][cut] if varDetail.dim == 2 else None
                            weight = weights["evtw"]
                            wKey = varDetail.weightName

                            # properly flatten certain inputs
                            if varDetail.flattenInfo >= 1:
                                vX = ak.flatten(vX)
                                vY = ak.flatten(vY) if varDetail.dim == 2 else None

                            # make sure the correct weights are applied
                            if wKey in weights.keys():
                                hW = weights[wKey]
                            elif wKey == "w1":
                                hW = np.ones(len(weight))
                            else:
                                hW = weight

                            if len(vX) > 0:
                                if   varDetail.dim == 1:                        
                                    output['h_{}{}'.format(histName,cutName)].fill(x=vX, weight=hW)
                                elif varDetail.dim == 2:
                                    output['h_{}{}'.format(histName,cutName)].fill(x=vX, y=vY, weight=hW)
                        
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
