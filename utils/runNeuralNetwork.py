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
#seed = 12345
#random.seed(12345)

def normalize(df,normMean,normStd):
    return (df-normMean)/normStd

def get_all_vars(varsIn,varSet,normMean,normStd):
    dSets = []
    dataSet = pd.DataFrame()
    for var in varSet:
        inputArr = varsIn[var]
        if variables()[var].npzInfo == 2:
            inputArr = np.repeat(ak.to_numpy(inputArr),ak.to_numpy(varsIn["njetsAK8"]))
        if variables()[var].flattenInfo == 1:
            inputArr = ak.flatten(inputArr)
        elif variables()[var].flattenInfo == 2:
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

def extrapolateNTaggedJets(inputArray, nRank, counts):
    # Extrapolate N+nRank tagged jets from N tagged jets

    prediction = ak.to_list(ak.zeros_like(inputArray))
    for i, a in enumerate(inputArray):
        pNum = [0.0]
        pDen = [1.0]
        factor = 1.0
        if len(a) >= nRank: 
            n = nRank + counts[i]
            k = counts[i]
            choose = np.math.factorial(n)/(np.math.factorial(k)*np.math.factorial(n-k))
            factor = 1.0 / choose
            pb = PoiBin(ak.to_numpy(a))
            pNum = pb.pmf([nRank])
            pDen =  pb.pmf([0])

        prediction[i] = factor*(pNum[0]/pDen[0])
    
    return ak.Array(prediction)

def findFakeRate(fakerateHisto, bgroundJetsAK8):
    pt = bgroundJetsAK8.pt
    eta = bgroundJetsAK8.eta
    value = pt

    if fakerateHisto is None:
        return ak.ones_like(value)
    x, y = fakerateHisto.values()                

    fakerate = ak.to_list(ak.zeros_like(value))
    for i, a1 in enumerate(value):
        for j, v in enumerate(a1):
            idx = np.absolute(x-v).argmin()
            fakerate[i][j] = y[idx]

    return ak.Array(fakerate)

def getFlatScore(nnOutput):
    output = ak.to_list(ak.zeros_like(nnOutput))
    for i, a in enumerate(nnOutput):
        output[i] = random.uniform(0, 1)
    return output

def runNN(model,varsIn,varSet,normMean,normStd,fakerateHisto):
    dataset = RootDataset(varsIn=varsIn,varSet=varSet, normMean=normMean, normStd=normStd)
    nnOutput = getNNOutput(dataset, model)
    #nnOutput = getFlatScore(nnOutput)
    fjets = varsIn["fjets"]
    counts = ak.num(fjets.pt)
    svjJetsAK8 = ak.unflatten(nnOutput, counts)

    wpt = 0.7
    darksvjJetsAK8 = fjets[svjJetsAK8 >= wpt]
    bgroundJetsAK8 = fjets[svjJetsAK8 < wpt]
    nsvjJetsAK8 = ak.num(darksvjJetsAK8)
    varsIn['nsvjJetsAK8'] = nsvjJetsAK8
    varsIn['nnOutput'] = svjJetsAK8
    varsIn['svfjw'] = u.awkwardReshape(darksvjJetsAK8,varsIn['evtw'])
    varsIn['svjPtAK8'] = darksvjJetsAK8.pt
    varsIn['svjEtaAK8'] = darksvjJetsAK8.eta

    #######################################################
    # Extrapolate number of tag jets from low to high
    #######################################################
    fakerate = findFakeRate(fakerateHisto, bgroundJetsAK8)
    #fakerate = 0.3*ak.ones_like(svjJetsAK8[svjJetsAK8 < wpt])

    nsvjJetsAK8_pred1jets = extrapolateNTaggedJets(fakerate, 1, nsvjJetsAK8)
    nsvjJetsAK8_pred2jets = extrapolateNTaggedJets(fakerate, 2, nsvjJetsAK8)
    nsvjJetsAK8_pred3jets = extrapolateNTaggedJets(fakerate, 3, nsvjJetsAK8)
    nsvjJetsAK8_pred4jets = extrapolateNTaggedJets(fakerate, 4, nsvjJetsAK8)

    varsIn['nsvjJetsAK8_pred1jets'] = nsvjJetsAK8_pred1jets
    varsIn['nsvjJetsAK8_pred2jets'] = nsvjJetsAK8_pred2jets
    varsIn['nsvjJetsAK8_pred3jets'] = nsvjJetsAK8_pred3jets
    varsIn['nsvjJetsAK8_pred4jets'] = nsvjJetsAK8_pred4jets
    varsIn['pred1_evtw'] = varsIn['evtw']*nsvjJetsAK8_pred1jets
    varsIn['pred2_evtw'] = varsIn['evtw']*nsvjJetsAK8_pred2jets
    varsIn['pred3_evtw'] = varsIn['evtw']*nsvjJetsAK8_pred3jets
    varsIn['pred4_evtw'] = varsIn['evtw']*nsvjJetsAK8_pred4jets
    varsIn['nsvjJetsAK8Plus1'] = ak.num(darksvjJetsAK8) + 1.0
    varsIn['nsvjJetsAK8Plus2'] = ak.num(darksvjJetsAK8) + 2.0
    varsIn['nsvjJetsAK8Plus3'] = ak.num(darksvjJetsAK8) + 3.0
    varsIn['nsvjJetsAK8Plus4'] = ak.num(darksvjJetsAK8) + 4.0

