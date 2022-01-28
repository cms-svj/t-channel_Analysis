import numpy as np
import awkward1 as ak
import torch.utils.data as udata
import torch
import pandas as pd
from .variables import variables
from torch.nn import functional as f
import matplotlib.pyplot as plt
from coffea.analysis_objects import JaggedCandidateArray

def normalize(df,normMean,normStd):
    return (df-normMean)/normStd

def get_all_vars(varsIn,varSet,normMean,normStd):
    dSets = []
    dataSet = pd.DataFrame()
    for var in varSet:
        inputArr = varsIn[var][0]
        if variables[var][4] == 2:
            inputArr = np.repeat(ak.to_numpy(inputArr),ak.to_numpy(varsIn["njetsAK8"][0]))
        if variables[var][5] == 1:
            inputArr = inputArr.flatten()
        elif variables[var][5] == 2:
            inputArr = ak.flatten(inputArr)
        dataSet[var] = inputArr
    dataSet = normalize(dataSet,normMean,normStd)
    return dataSet

class RootDataset(udata.Dataset):
    def __init__(self, varsIn, varSet, normMean, normStd):
        dataSet = get_all_vars(varsIn, varSet, normMean, normStd)
        self.vars = dataSet.astype(float).values

    def __len__(self):
        return len(self.vars)

    def __getitem__(self, idx):
        data_np = self.vars[idx].copy()
        data  = torch.from_numpy(data_np)
        return data

def getNNOutput(dataset, model):
    output_tag = np.array([])
    if dataset.__len__() > 0:
        loader = udata.DataLoader(dataset=dataset, batch_size=dataset.__len__(), num_workers=0)
        d = next(iter(loader))
        data = d.float()
        model.eval()
        out_tag, out_pTClass = model(data)
        output_tag = f.softmax(out_tag,dim=1)[:,1].detach().numpy()
        output_tag = np.nan_to_num(output_tag,nan=-1)
    return output_tag

def runNN(model,varsIn,varSet,normMean,normStd):
    dataset = RootDataset(varsIn=varsIn,varSet=varSet, normMean=normMean, normStd=normStd)
    nnOutput = getNNOutput(dataset, model)
    fjets = varsIn["fjets"]
    svjJetsAK8 = JaggedCandidateArray.candidatesfromcounts(
        fjets.counts,
        pt=fjets.pt.flatten(),
        eta=fjets.eta.flatten(),
        phi=fjets.phi.flatten(),
        mass=fjets.mass.flatten(),
        nnOutput = nnOutput
    )
    wpt = 0.5
    darksvjJetsAK8 = svjJetsAK8[svjJetsAK8.nnOutput >= wpt]
    varsIn['nsvjJetsAK8'] = [darksvjJetsAK8.counts,"evtw"]
    varsIn['nnOutput'] = [svjJetsAK8.nnOutput,"fjw"]
