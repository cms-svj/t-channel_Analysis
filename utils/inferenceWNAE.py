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
def get_data(varsIn, fjets, year, category):
    scores = varsIn[""] #scores_file["Scores"]["score1"].array(library="np")
    training_file = uproot.open(f"test_dataset_{year}.root")

    j0_pt = training_file["Events"]["GoodJetsAK80_pt"].array(library="np")
    j1_pt = training_file["Events"]["GoodJetsAK81_pt"].array(library="np")
    j2_pt = training_file["Events"]["GoodJetsAK82_pt"].array(library="np")
    j3_pt = training_file["Events"]["GoodJetsAK83_pt"].array(library="np")

    j0_wnae_pt_0_200 = training_file["Events"]["GoodJetsAK80_WNAEPt0To200Loss"].array(library="np")
    j1_wnae_pt_0_200 = training_file["Events"]["GoodJetsAK81_WNAEPt0To200Loss"].array(library="np")
    j2_wnae_pt_0_200 = training_file["Events"]["GoodJetsAK82_WNAEPt0To200Loss"].array(library="np")
    j3_wnae_pt_0_200 = training_file["Events"]["GoodJetsAK83_WNAEPt0To200Loss"].array(library="np")

    j0_wnae_pt_200_300 = training_file["Events"]["GoodJetsAK80_WNAEPt200To300Loss"].array(library="np")
    j1_wnae_pt_200_300 = training_file["Events"]["GoodJetsAK81_WNAEPt200To300Loss"].array(library="np")
    j2_wnae_pt_200_300 = training_file["Events"]["GoodJetsAK82_WNAEPt200To300Loss"].array(library="np")
    j3_wnae_pt_200_300 = training_file["Events"]["GoodJetsAK83_WNAEPt200To300Loss"].array(library="np")

    j0_wnae_pt_300_400 = training_file["Events"]["GoodJetsAK80_WNAEPt300To400Loss"].array(library="np")
    j1_wnae_pt_300_400 = training_file["Events"]["GoodJetsAK81_WNAEPt300To400Loss"].array(library="np")
    j2_wnae_pt_300_400 = training_file["Events"]["GoodJetsAK82_WNAEPt300To400Loss"].array(library="np")
    j3_wnae_pt_300_400 = training_file["Events"]["GoodJetsAK83_WNAEPt300To400Loss"].array(library="np")

    j0_wnae_pt_400_500 = training_file["Events"]["GoodJetsAK80_WNAEPt400To500Loss"].array(library="np")
    j1_wnae_pt_400_500 = training_file["Events"]["GoodJetsAK81_WNAEPt400To500Loss"].array(library="np")
    j2_wnae_pt_400_500 = training_file["Events"]["GoodJetsAK82_WNAEPt400To500Loss"].array(library="np")
    j3_wnae_pt_400_500 = training_file["Events"]["GoodJetsAK83_WNAEPt400To500Loss"].array(library="np")

    j0_wnae_pt_500_inf = training_file["Events"]["GoodJetsAK80_WNAEPt500ToInfLoss"].array(library="np")
    j1_wnae_pt_500_inf = training_file["Events"]["GoodJetsAK81_WNAEPt500ToInfLoss"].array(library="np")
    j2_wnae_pt_500_inf = training_file["Events"]["GoodJetsAK82_WNAEPt500ToInfLoss"].array(library="np")
    j3_wnae_pt_500_inf = training_file["Events"]["GoodJetsAK83_WNAEPt500ToInfLoss"].array(library="np")

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
    year = "2016"
    category = "test"
    n_tagged_jets = get_data(varsIn, fjets, year, category)

    # wpt = 0.7 # wpt 0.8, 0.045 fakerate
    # darksvjJetsAK8 = fjets[svjJetsAK8 >= wpt]
    # bgroundJetsAK8 = fjets[svjJetsAK8 < wpt]
    # nsvjJetsAK8 = ak.num(darksvjJetsAK8)
    varsIn['nsvjJetsAK8WNAE'] = n_tagged_jets
    # varsIn['pNetJetTaggerScore'] = svjJetsAK8
    # varsIn['JetsAK8_pNetJetTaggerScore'] = svjJetsAK8 # needed for event classifier part to work properly; should probably unify the naming conventions
    # varsIn['svfjw'] = u.awkwardReshape(darksvjJetsAK8,varsIn['evtw'])
    # varsIn['svjPtAK8'] = darksvjJetsAK8.pt
    # varsIn['svjEtaAK8'] = darksvjJetsAK8.eta

    # varsIn['nsvjJetsAK8Plus1'] = ak.num(darksvjJetsAK8) + 1.0
    # varsIn['nsvjJetsAK8Plus2'] = ak.num(darksvjJetsAK8) + 2.0
    # varsIn['nsvjJetsAK8Plus3'] = ak.num(darksvjJetsAK8) + 3.0
    # varsIn['nsvjJetsAK8Plus4'] = ak.num(darksvjJetsAK8) + 4.0
