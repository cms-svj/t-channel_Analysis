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

def jetVar_i(var,i,padValue=np.Inf):
    paddedVar = ak.fill_none(ak.pad_none(var,i+1),padValue)
    return paddedVar[:,i]

#seed = 12345
#random.seed(12345)
def get_data(varsIn, fjets):
    j_wnae_pt_0_200 = varsIn["JetsAK8_WNAEPt0To200Loss"]
    j_wnae_pt_200_300 = varsIn["JetsAK8_WNAEPt200To300Loss"]
    j_wnae_pt_300_400 = varsIn["JetsAK8_WNAEPt300To400Loss"]
    j_wnae_pt_400_500 = varsIn["JetsAK8_WNAEPt400To500Loss"]
    j_wnae_pt_500_inf = varsIn["JetsAK8_WNAEPt500ToInfLoss"]

    j0_pt = jetVar_i(fjets.pt, 0)
    j1_pt = jetVar_i(fjets.pt, 1)
    j2_pt = jetVar_i(fjets.pt, 2)
    j3_pt = jetVar_i(fjets.pt, 3)

    j0_wnae_pt_0_200 = jetVar_i(j_wnae_pt_0_200, 0)
    j1_wnae_pt_0_200 = jetVar_i(j_wnae_pt_0_200, 1)
    j2_wnae_pt_0_200 = jetVar_i(j_wnae_pt_0_200, 2)
    j3_wnae_pt_0_200 = jetVar_i(j_wnae_pt_0_200, 3)

    j0_wnae_pt_200_300 = jetVar_i(j_wnae_pt_200_300, 0)
    j1_wnae_pt_200_300 = jetVar_i(j_wnae_pt_200_300, 1)
    j2_wnae_pt_200_300 = jetVar_i(j_wnae_pt_200_300, 2)
    j3_wnae_pt_200_300 = jetVar_i(j_wnae_pt_200_300, 3)

    j0_wnae_pt_300_400 = jetVar_i(j_wnae_pt_300_400, 0)
    j1_wnae_pt_300_400 = jetVar_i(j_wnae_pt_300_400, 1)
    j2_wnae_pt_300_400 = jetVar_i(j_wnae_pt_300_400, 2)
    j3_wnae_pt_300_400 = jetVar_i(j_wnae_pt_300_400, 3)

    j0_wnae_pt_400_500 = jetVar_i(j_wnae_pt_400_500, 0)
    j1_wnae_pt_400_500 = jetVar_i(j_wnae_pt_400_500, 1)
    j2_wnae_pt_400_500 = jetVar_i(j_wnae_pt_400_500, 2)
    j3_wnae_pt_400_500 = jetVar_i(j_wnae_pt_400_500, 3)

    j0_wnae_pt_500_inf = jetVar_i(j_wnae_pt_500_inf, 0)
    j1_wnae_pt_500_inf = jetVar_i(j_wnae_pt_500_inf, 1)
    j2_wnae_pt_500_inf = jetVar_i(j_wnae_pt_500_inf, 2)
    j3_wnae_pt_500_inf = jetVar_i(j_wnae_pt_500_inf, 3)

    j0_is_tagged = (j0_pt < 200)*(j0_wnae_pt_0_200 > 25.156) + \
                (j0_pt >= 200)*(j0_pt < 300)*(j0_wnae_pt_200_300 > 18.284) + \
                (j0_pt >= 300)*(j0_pt < 400)*(j0_wnae_pt_300_400 > 20.383) + \
                (j0_pt >= 400)*(j0_pt < 500)*(j0_wnae_pt_400_500 > 21.941) + \
                (j0_pt >= 500)*(j0_wnae_pt_500_inf > 16.370)

    j1_is_tagged = (j1_pt < 200)*(j1_wnae_pt_0_200 > 25.156) + \
                (j1_pt >= 200)*(j1_pt < 300)*(j1_wnae_pt_200_300 > 18.284) + \
                (j1_pt >= 300)*(j1_pt < 400)*(j1_wnae_pt_300_400 > 20.383) + \
                (j1_pt >= 400)*(j1_pt < 500)*(j1_wnae_pt_400_500 > 21.941) + \
                (j1_pt >= 500)*(j1_wnae_pt_500_inf > 16.370)

    j2_is_tagged = (j2_pt < 200)*(j2_wnae_pt_0_200 > 25.156) + \
                (j2_pt >= 200)*(j2_pt < 300)*(j2_wnae_pt_200_300 > 18.284) + \
                (j2_pt >= 300)*(j2_pt < 400)*(j2_wnae_pt_300_400 > 20.383) + \
                (j2_pt >= 400)*(j2_pt < 500)*(j2_wnae_pt_400_500 > 21.941) + \
                (j2_pt >= 500)*(j2_wnae_pt_500_inf > 16.370)

    j3_is_tagged = (j3_pt < 200)*(j3_wnae_pt_0_200 > 25.156) + \
                (j3_pt >= 200)*(j3_pt < 300)*(j3_wnae_pt_200_300 > 18.284) + \
                (j3_pt >= 300)*(j3_pt < 400)*(j3_wnae_pt_300_400 > 20.383) + \
                (j3_pt >= 400)*(j3_pt < 500)*(j3_wnae_pt_400_500 > 21.941) + \
                (j3_pt >= 500)*(j3_wnae_pt_500_inf > 16.370)

    n_tagged_jets = j0_is_tagged.astype(int) + j1_is_tagged.astype(int) + j2_is_tagged.astype(int) + j3_is_tagged.astype(int)

    return n_tagged_jets

def create_wnae_related_variables(varsIn, fjets):
    varsIn['nsvjJetsAK8WNAE'] = get_data(varsIn, fjets)

    # wpt = 0.7 # wpt 0.8, 0.045 fakerate
    # darksvjJetsAK8 = fjets[svjJetsAK8 >= wpt]
    # bgroundJetsAK8 = fjets[svjJetsAK8 < wpt]
    # nsvjJetsAK8 = ak.num(darksvjJetsAK8)
    # varsIn['nsvjJetsAK8WNAE'] = n_tagged_jets
    # varsIn['pNetJetTaggerScore'] = svjJetsAK8
    # varsIn['JetsAK8_pNetJetTaggerScore'] = svjJetsAK8 # needed for event classifier part to work properly; should probably unify the naming conventions
    # varsIn['svfjw'] = u.awkwardReshape(darksvjJetsAK8,varsIn['evtw'])
    # varsIn['svjPtAK8'] = darksvjJetsAK8.pt
    # varsIn['svjEtaAK8'] = darksvjJetsAK8.eta

    # varsIn['nsvjJetsAK8Plus1'] = ak.num(darksvjJetsAK8) + 1.0
    # varsIn['nsvjJetsAK8Plus2'] = ak.num(darksvjJetsAK8) + 2.0
    # varsIn['nsvjJetsAK8Plus3'] = ak.num(darksvjJetsAK8) + 3.0
    # varsIn['nsvjJetsAK8Plus4'] = ak.num(darksvjJetsAK8) + 4.0
