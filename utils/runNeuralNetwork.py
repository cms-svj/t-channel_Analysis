import numpy as np
import awkward as ak
import torch.utils.data as udata
import torch
import pandas as pd
from .variables import variables
from torch.nn import functional as f
from utils import utility as u

def normalize(df,normMean,normStd):
    return (df-normMean)/normStd

def get_all_vars(varsIn,varSet,normMean,normStd):
    dSets = []
    dataSet = pd.DataFrame()
    for var in varSet:
        inputArr = varsIn[var]
        if variables[var][4] == 2:
            inputArr = np.repeat(ak.to_numpy(inputArr),ak.to_numpy(varsIn["njetsAK8"])
        if variables[var][5] == 1:
            inputArr = ak.flatten(inputArr)
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
    counts = ak.num(fjets.pt)
    svjJetsAK8 = ak.unflatten(nnOutput, counts)

    wpt = 0.5
    darksvjJetsAK8 = fjets[svjJetsAK8 >= wpt]
    varsIn['nsvjJetsAK8'] = ak.num(darksvjJetsAK8)
    varsIn['nnOutput'] = svjJetsAK8
