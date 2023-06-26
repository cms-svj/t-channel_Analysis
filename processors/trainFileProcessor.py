from coffea import processor
import numpy as np
import awkward as ak
from utils import utility as utl
from utils import baseline as bl
from utils.python import jetutils as ju
from utils.variables import variables
from utils.python.svjgnntagger import SVJGNNTagger

def col_accumulator(a):
    return processor.column_accumulator(np.array(a))

class MainProcessor(processor.ProcessorABC):
        def __init__(self,dataset,sf,*args):
            self._accumulator = processor.dict_accumulator({})
            self.dataset = dataset
            self.setupNPArr = None
            self.scaleFactor = sf
            self.hemPeriod = ""
        @property
        def accumulator(self):
                return self._accumulator

        def setupNPArray(self,inputVars,inputShapes,numOfJetClasses):
            varDict = {}
            for v in inputVars:
                inputShape = inputShapes[v]
                varDict[v] = processor.column_accumulator(np.empty((0,inputShape[1],inputShape[2])))
            varDict["pT"] = processor.column_accumulator(np.empty((0)))
            varDict["signal"] = processor.column_accumulator(np.empty((0,numOfJetClasses)))
            varDict["weight"] = processor.column_accumulator(np.empty((0)))
            self._accumulator = processor.dict_accumulator(varDict)
            self.setupNPArr = True

        def process(self, events):
                ## objects used for cuts
                vars_noCut = utl.baselineVar(self.dataset,events,self.hemPeriod,self.scaleFactor)
                gnn = SVJGNNTagger( model_structure='utils.data.GNNTagger.SVJTagger',
                                    model_inputs='./utils/data/GNNTagger/svj.yaml')
                inputVars = gnn.data_config.input_names
                inputShapes = gnn.data_config.input_shapes
                numOfJetClasses = len(gnn.data_config.label_value)
                # Our preselection
                cuts = bl.cutList(self.dataset,events,vars_noCut,self.hemPeriod,SVJCut=False)
                if self.setupNPArr is None:                 
                    self.setupNPArray(inputVars,inputShapes,numOfJetClasses)
                output = self.accumulator.identity()

                # run cut loop
                cut = cuts["_qual_trg_st_1PJ"]
                if (len(events) > 0) and (np.any(cut)):
                    events = events[cut]
                    fjets = vars_noCut["fjets"][cut]
                    evtw = vars_noCut["evtw"][cut]
                    jets_in = ju.run_jet_constituent_matching(events, fjets)
                    jets_in = ak.flatten(jets_in)
                    feature_map = gnn.get_feature_map(jets_in)
                    X = gnn.structure_X(jets_in,feature_map)
                    for inputLabel in inputVars:
                        output[inputLabel] += col_accumulator(X[inputLabel])
                    weight = utl.awkwardReshape(fjets,evtw)
                    flattened_pT = ak.flatten(fjets.pt)
                    output["pT"] += col_accumulator(flattened_pT)
                    output["weight"] += col_accumulator(ak.flatten(weight))

                    # signal label
                    signalLabels = {
                        "QCD": 0,
                        "TTJets": 1,
                        "mMed": 2,
                    }
                    sigLabIndex = -1
                    for sigLab in signalLabels.keys():
                        if sigLab in self.dataset:
                            sigLabIndex = signalLabels[sigLab]
                    signal = np.zeros(numOfJetClasses)
                    signal[sigLabIndex] = 1
                    output["signal"] += col_accumulator(np.tile(signal,(len(flattened_pT),1)))

                return output

        def postprocess(self, accumulator):
                return accumulator
