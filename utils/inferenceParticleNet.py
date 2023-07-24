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
from utils.python import jetutils as ju
from utils.python.svjgnntagger import SVJGNNTagger
from scipy.special import softmax

#seed = 12345
#random.seed(12345)

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

def runJetTagger(events,varsIn,fakerateHisto):
    gnn_triton = SVJGNNTagger(score_tag='score',
            triton_path='triton+grpc://triton.apps.okddev.fnal.gov:443/svj_tch_gnn/1',
            model_structure='utils.data.GNNTagger.SVJTagger',
            model_inputs='./utils/data/GNNTagger/svj.yaml',
            dec_thresh=0.999)
    # initialize model if not already done
    gnn_triton.use_triton = True
    if gnn_triton.model == None:
        gnn_triton.initialize_model()
    #nnOutput = getFlatScore(nnOutput)
    fjets = varsIn["fjets"]
    jets_in = ju.run_jet_constituent_matching(events, fjets)
    jets_in = ak.flatten(jets_in)
    batch_size = 1024
    nnOutput = np.array([])
    for ii in range(0,len(jets_in),batch_size):
        try:
            jets_eval = jets_in[ii:ii+batch_size]
        except:
            jets_eval = jets_in[ii:-1]
        # put data in proper format
        feature_map = gnn_triton.get_feature_map(jets_eval)
        X = gnn_triton.structure_X(jets_eval,feature_map)
        #inference with triton 
        outputs = gnn_triton.model(X)
        nnOutput = np.append(nnOutput,softmax(outputs, axis=-1)[:, 2]) # 2 is the label for SVJ_Dark, 0 = QCD, 1 = TTJets
    counts = ak.num(fjets.pt)
    svjJetsAK8 = ak.unflatten(nnOutput, counts)

    wpt = 0.4
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

