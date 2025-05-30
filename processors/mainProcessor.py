from coffea import processor
import hist as h
import numpy as np
import awkward as ak
from utils import utility as utl
from utils import baseline as bl
from utils.variables import variables
from utils.inferenceParticleNet import runJetTagger, create_pn_related_variables
from utils.inferenceWNAE import create_wnae_related_variables
from utils.runEventTagger import runEventTagger
import uproot
import torch
from datetime import datetime

class MainProcessor(processor.ProcessorABC):
        def __init__(self,**kwargs):
                self._accumulator = processor.dict_accumulator({})
                self.setupHistos = None
                self.jNVar = kwargs["jNVar"]
                self.fakerateHisto = self.getHistoFromFile("fakerate.root", "jPt_Fakerate_SR;1") 
                self.hemPeriod = kwargs["hemPeriod"]
                self.hemStudy = kwargs["hemStudy"]
                self.trgEffStudy = kwargs["trgEffStudy"]
                self.evtTaggerDict = kwargs["evtTaggerDict"]
                self.sFactor = kwargs["sFactor"]
                self.skimSource = kwargs["skimSource"]
                self.skimCut = kwargs["skimCut"]
                self.runJetTag = kwargs["runJetTag"]
                self.runEvtClass = kwargs["runEvtClass"]
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
                    for histName, histDetail in variables(self.jNVar,runJetTag=self.runJetTag,runEvtClass=self.runEvtClass).items():
                        if   histDetail.dim == 1:
                            histograms['h_{}{}'.format(histName,cutName)] = h.Hist(histDetail.xbins, storage="weight")
                        elif histDetail.dim == 2:                        
                            histograms['h_{}{}'.format(histName,cutName)] = h.Hist(histDetail.xbins, histDetail.ybins, storage="weight")

                self._accumulator = histograms
                self.setupHistos = True

        def process(self, events):
                # cut loop
                ## objects used for cuts
                dataset = events.metadata['dataset']
                vars_noCut = utl.baselineVar(dataset,events,self.hemPeriod,self.sFactor,self.skimSource,self.runJetTag)
                if self.runJetTag:
                    if self.skimSource:
                        create_pn_related_variables(vars_noCut, self.fakerateHisto, vars_noCut["fjets"], vars_noCut["JetsAK8_pNetJetTaggerScore"][vars_noCut["JetsAK8_isGood"]])
                        create_wnae_related_variables(vars_noCut, vars_noCut["fjets"], 
                                                      svjWNAEPt0To200Loss = vars_noCut["JetsAK8_WNAEPt0To200Loss"][vars_noCut["JetsAK8_isGood"]],   wpt0To200 = 25.156,
                                                      svjWNAEPt200To300Loss = vars_noCut["JetsAK8_WNAEPt200To300Loss"][vars_noCut["JetsAK8_isGood"]], wpt200To300 = 18.284,
                                                      svjWNAEPt300To400Loss = vars_noCut["JetsAK8_WNAEPt300To400Loss"][vars_noCut["JetsAK8_isGood"]], wpt300To400 = 20.383,
                                                      svjWNAEPt400To500Loss = vars_noCut["JetsAK8_WNAEPt400To500Loss"][vars_noCut["JetsAK8_isGood"]], wpt400To500 = 21.941,
                                                      svjWNAEPt500ToInfLoss = vars_noCut["JetsAK8_WNAEPt500ToInfLoss"][vars_noCut["JetsAK8_isGood"]], wpt500ToInf = 16.370) 
                        # these WNAE working points are based on a 20% background jet rejection rate
                    else:
                        runJetTagger(events,vars_noCut,self.fakerateHisto)

                utl.varGetter(dataset,events,vars_noCut,np.ones(len(events),dtype=bool),self.jNVar)
                if self.runEvtClass:
                    runEventTagger(events, vars_noCut, self.skimSource, self.evtTaggerDict)                    
                cuts = bl.cutList(dataset,events,vars_noCut,self.hemStudy,self.trgEffStudy,self.hemPeriod,self.skimCut,self.skimSource,self.runJetTag)
                # setup histograms
                if self.setupHistos is None:
                    self.setupHistogram(cuts)
                output = self.accumulator
                # run cut loop
                for cutName,cut in cuts.items():
                    # print("cutName = {} \n cut = {}".format(cutName,cut))
                    # defining objects
                    weights = {
                            "evtw" : vars_noCut["evtw"][cut],
                            "jw"   : ak.flatten(vars_noCut["jw"][cut]),
                            "fjw"  : ak.flatten(vars_noCut["fjw"][cut]),
                            "ew"   : ak.flatten(vars_noCut["ew"][cut]),
                            "mw"   : ak.flatten(vars_noCut["mw"][cut]),
                            "crew" : ak.flatten(vars_noCut["crew"][cut]),
                            "crmw" : ak.flatten(vars_noCut["crmw"][cut]),
                            "nimw" : ak.flatten(vars_noCut["nimw"][cut]),
                            # "svfjw" : ak.flatten(vars_noCut["svfjw"][cut]),
                            # "pred1_evtw" : vars_noCut["pred1_evtw"][cut],
                            # "pred2_evtw" : vars_noCut["pred2_evtw"][cut],
                            # "pred3_evtw" : vars_noCut["pred3_evtw"][cut],
                            # "pred4_evtw" : vars_noCut["pred4_evtw"][cut],
                    }
                    if len(events) > 0:
                        ## filling histograms
                        for histName, varDetail in variables(self.jNVar,runJetTag=self.runJetTag,runEvtClass=self.runEvtClass).items():   
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
                                    if self.jNVar:
                                        finiteMask = np.isfinite(vX)
                                        vX = vX[finiteMask]
                                        hW = hW[finiteMask]
                                    output['h_{}{}'.format(histName,cutName)].fill(x=vX, weight=hW)
                                elif varDetail.dim == 2:
                                    if self.jNVar:
                                        finiteMask = np.isfinite(vX) & np.isfinite(vY)
                                        vX = vX[finiteMask]
                                        vY = vY[finiteMask]
                                        hW = hW[finiteMask]
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
