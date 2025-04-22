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

def create_wnae_related_variables(varsIn, fjets, 
                                  svjWNAEPt0To200Loss, wpt0To200, 
                                  svjWNAEPt200To300Loss, wpt200To300, 
                                  svjWNAEPt300To400Loss, wpt300To400, 
                                  svjWNAEPt400To500Loss, wpt400To500, 
                                  svjWNAEPt500ToInfLoss, wpt500ToInf):
    jetsAK8_pT = fjets.pt

    ptCon0to200 = (jetsAK8_pT >= 0) & (jetsAK8_pT < 200)
    ptCon200to300 = (jetsAK8_pT >= 200) & (jetsAK8_pT < 300)
    ptCon300to400 = (jetsAK8_pT >= 300) & (jetsAK8_pT < 400)
    ptCon400to500 = (jetsAK8_pT >= 400) & (jetsAK8_pT < 500)
    ptCon500toInf = (jetsAK8_pT >= 500)

    jetsAK8_wnaeScore = (
                            ak.where(ptCon0to200, svjWNAEPt0To200Loss, 0) 
                            + ak.where(ptCon200to300, svjWNAEPt200To300Loss, 0)
                            + ak.where(ptCon300to400, svjWNAEPt300To400Loss, 0)
                            + ak.where(ptCon400to500, svjWNAEPt400To500Loss, 0)
                            + ak.where(ptCon500toInf, svjWNAEPt500ToInfLoss, 0)
                        )
    varsIn['wnaeJetTaggerScore'] = jetsAK8_wnaeScore
    darksvjWNAEJetsAK8 = ak.concatenate(
                                        [
                                            fjets[ptCon0to200 & (svjWNAEPt0To200Loss > wpt0To200)],
                                            fjets[ptCon200to300 & (svjWNAEPt200To300Loss > wpt200To300)],
                                            fjets[ptCon300to400 & (svjWNAEPt300To400Loss > wpt300To400)],
                                            fjets[ptCon400to500 & (svjWNAEPt400To500Loss > wpt400To500)],
                                            fjets[ptCon500toInf & (svjWNAEPt500ToInfLoss > wpt500ToInf)],
                                        ], axis = 1
                                       )
    varsIn['nsvjWNAE'] = ak.num(darksvjWNAEJetsAK8)
    varsIn['svfWNAEjw'] = u.awkwardReshape(darksvjWNAEJetsAK8,varsIn['evtw'])
    varsIn['svjWNAEPtAK8'] = darksvjWNAEJetsAK8.pt
    varsIn['svjWNAEEtaAK8'] = darksvjWNAEJetsAK8.eta


