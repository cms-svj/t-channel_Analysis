import numpy as np
import awkward as ak
import torch.utils.data as udata
import torch
import pandas as pd
from .variables import variables
from torch.nn import functional as f
from utils import utility as u
from utils.poibin import PoiBin
import random
from utils.data.DNNEventClassifier.models import DNN

#seed = 12345
#random.seed(12345)

def normalize(df,normMean,normStd):
    return (df-normMean)/normStd

def get_all_vars(varsIn,eventVar,jetVar,numOfJetsToUse):
    dSets = []
    dataSet = pd.DataFrame()
    ########### Normalization ###############
    for evar in eventVar:
        inputArr = np.array(varsIn[evar])
        if evar in ["njetsAK8","mT"]:
            dataSet[evar] = np.log(inputArr)
        else:
            dataSet[evar] = inputArr
    for jvar in jetVar:
        inputArr = varsIn[jvar]
        for i in range(numOfJetsToUse):
            jetiInput = np.array(u.jetVar_i(inputArr,i,padValue=0))
            if jvar in ["jPtAK8","jEAK8"]:
                dataSet["{}_{}".format(jvar,i)] = np.log(jetiInput)
            elif jvar == "jEtaAK8":
                dataSet["{}_{}".format(jvar,i)] = abs(jetiInput)
            else:
                dataSet["{}_{}".format(jvar,i)] = jetiInput
    ###########################################
    return dataSet

class RootDataset(udata.Dataset):
    def __init__(self, varsIn, eventVar, jetVar, numOfJetsToUse):
        dataSet = get_all_vars(varsIn, eventVar, jetVar, numOfJetsToUse)
        self.vars = dataSet.astype(float).values

    def __len__(self):
        return len(self.vars)

    def __getitem__(self, idx):
        data_np = self.vars[idx].copy()
        data  = torch.from_numpy(data_np)
        return data

def getNNOutput(dataset, model, num_classes):
    output_tag = np.array([])
    if dataset.__len__() > 0:
        loader = udata.DataLoader(dataset=dataset, batch_size=dataset.__len__(), num_workers=0)
        d = next(iter(loader))
        data = d.float()
        out_tag = model(data)
        signalIndex = num_classes - 1
        output_tag = f.softmax(out_tag,dim=1)[:,signalIndex].detach().numpy()
        output_tag = np.nan_to_num(output_tag,nan=-1)
    return output_tag,data

def runEventTagger(varsIn,evtTaggerDict):
    hyper = evtTaggerDict["hyper"]
    features = evtTaggerDict["features"]
    evtTaggerLocation = evtTaggerDict["evtTaggerLocation"]
    numOfJetsToUse = features.numOfJetsToKeep
    eventVariables = features.eventVariables
    jetVariables = features.jetVariables
    num_classes = hyper.num_classes
    device = torch.device('cpu')
    evtTagger = DNN(n_var=len(eventVariables)+len(jetVariables)*numOfJetsToUse, n_layers=hyper.num_of_layers, n_nodes=hyper.num_of_nodes, n_outputs=hyper.num_classes, drop_out_p=hyper.dropout).to(device=device)
    evtTagger.load_state_dict(torch.load("{}/model.pth".format(evtTaggerLocation),map_location=device))
    evtTagger.eval()
    evtTagger.to('cpu')
    dataset = RootDataset(varsIn=varsIn,eventVar=eventVariables, jetVar=jetVariables, numOfJetsToUse=numOfJetsToUse)
    nnOutput,data = getNNOutput(dataset, evtTagger, num_classes)    
    varsIn['nnEventOutput'] = nnOutput
    varsIn['nnEventOutputrMET'] = nnOutput/varsIn["met"]
    varsIn['nnEventOutputrST'] = nnOutput/varsIn["st"]
    return data
