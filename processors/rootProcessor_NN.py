from coffea import hist, processor
import numpy as np
import awkward as ak
from utils import utility as utl
import utils.objects as ob
from utils import baseline as bl
from utils.variables import variables
from utils.inferenceParticleNet import runJetTagger
import uproot

def col_accumulator(a):
    return processor.column_accumulator(np.array(a))

class MainProcessor(processor.ProcessorABC):
        def __init__(self,**kwargs):
            self._accumulator = processor.dict_accumulator({})
            self.setupNPArr = None
            self.jNVar = kwargs["jNVar"]
            self.tcut = kwargs["tcut"]
            self.hemPeriod = ""
            self.fakerateHisto = self.getHistoFromFile("fakerate.root", "jPt_Fakerate_SR;1") 
            self.sFactor = kwargs["sFactor"]
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

        def setupNPArray(self,variables,maxNJets):
            varDict = {}
            for v,d in variables.items():
                if d.npzInfo == 1:
                    varDict[v] = processor.column_accumulator(np.empty((0,maxNJets)))
                elif d.npzInfo == 2:
                    varDict[v] = processor.column_accumulator(np.zeros(shape=(0)))
            self._accumulator = processor.dict_accumulator(varDict)
            self.setupNPArr = True

        def process(self, events):
                ## objects used for cuts
                dataset = events.metadata['dataset']
                maxNJets = 10
                if self.setupNPArr is None:
                    self.setupNPArray(variables(self.jNVar),maxNJets)
                output = self.accumulator.identity()
                vars_noCut = utl.baselineVar(dataset,events,self.hemPeriod,self.sFactor)
                preCut = bl.preCut(dataset,events,vars_noCut,self.hemPeriod) # preselection cut
                events = events[preCut]
                if len(events) > 0:
                    vars_noCut = utl.baselineVar(dataset,events,self.hemPeriod,self.sFactor)
                    runJetTagger(events,vars_noCut,self.fakerateHisto)
                    # Our preselection
                    cuts = bl.cutList(dataset,events,vars_noCut,self.hemPeriod,SVJCut=False)
                    # run cut loop
                    cut = cuts[self.tcut]
                    # saving the number of events before any cut for consistency check
                    if (len(events) > 0) and (np.any(cut)):
                        eventsCut = events[cut]
                        if (len(eventsCut) > 0):
                            utl.varGetter(dataset,eventsCut,vars_noCut,cut,self.jNVar)
                            for varName,varDetail in variables(self.jNVar).items():
                                # only store jetAK8 variables
                                if varDetail.npzInfo == 1:
                                    if varName == "nnOutput":
                                        hIn = vars_noCut[varName][cut]
                                    else:
                                        hIn = vars_noCut[varName]
                                    hIn = ak.to_numpy(ak.fill_none(
                                    ak.pad_none(hIn, maxNJets, axis=-1, clip=True),
                                    0))
                                    hIn = ak.nan_to_num(hIn,nan=0.0,posinf=0.0,neginf=0.0)
                                    output['{}'.format(varName)] += col_accumulator(hIn)
                                elif varDetail.npzInfo == 2:
                                    hIn = vars_noCut[varName]
                                    hIn = ak.nan_to_num(hIn,nan=0.0,posinf=0.0,neginf=0.0)
                                    output['{}'.format(varName)] += col_accumulator(hIn)
                return output

        def postprocess(self, accumulator):
                return accumulator
