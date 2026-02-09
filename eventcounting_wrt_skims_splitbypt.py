#!/usr/bin/env python3
"""
Reproduce MET plots from skims and compare to t-channel output.
python3 eventcounting_wrt_skims.py   --year 2018   --sample 2016_QCD   --skimCut t_channel_pre_selection_WNAE   --tchan-root /us
cms/home/ashrivas/nobackup/Dark_Sector/t-channel_Analysis/output/Current_Model_wp90_w2016/2016_QCD.root   --tchan-hist h_MET_pre   --applyPU   --xmin 0 --xmax 4000 --nbins 120   --outdir ./ABCD_hists3/met_compare_2018QCD


You asked for a FULL script that makes 3 plots:
  (1) RAW skims MET (from Events.MET) overlaid with t-channel output hist 'h_MET_pre'
  (2) SCALED skims MET only (apply t-channel-like normalization via event weights)
  (3) SCALED skims MET overlaid with t-channel output hist 'h_MET_pre'

Also writes:
  - met_from_skims_raw.txt  (one MET value per event, concatenated across ALL skim files)
  - met_from_skims_w.txt    (one weight per event, aligned with MET list; for debugging)
  - per-file event counters + sumw summary in skim_summary.txt

Notes:
- This reads MET and Weight (and optionally puWeight) directly from the skim ROOT TTrees.
- "Scaling the same way tchannel does" in coffea path is per-event:
    evtw = luminosity * events.Weight * scaleFactor * kFactor * (puWeight for bkg)
  Here we implement the same structure:
    w = lumi * Weight * kFactor * scaleFactor * puWeight
  with defaults: scaleFactor=1 unless you pass --scaleFactor.
- For 2018, luminosity and kFactor depend on hemPeriod. Use --hemPeriod All/PreHEM/PostHEM.

Run:
  python met_skims_vs_tchannel.py \
    --tchan-root /path/to/2018_QCD.root \
    --tchan-hist h_MET_pre \
    --outdir ./met_compare_2018QCD \
    --hemPeriod All

If XRootD is slow, you can test with --maxFiles 5.
"""

import os
import sys
import math
import argparse
import numpy as np
import uproot
import matplotlib.pyplot as plt

# ---------------------------
# CONFIG (default candidates)
# ---------------------------
TREE_CANDIDATES = ["Events", "EventTree", "tree", "ntuple/tree"]
MET_BRANCH_CANDIDATES = ["MET", "MET_pt", "met", "met_pt", "PuppiMET_pt", "PuppiMET", "puppiMET_pt"]
WEIGHT_BRANCH_CANDIDATES = ["Weight", "weight", "genWeight", "Generator_weight"]
PUWEIGHT_BRANCH_CANDIDATES = ["puWeight", "PUWeight", "puweight", "PileupWeight"]

# ---------------------------
# ALL SKIM FILES (EXACTLY as you pasted)
# ---------------------------
SKIM_FILES_2018 = [
    # 300-470
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_300to470/part-0.root",

    # 470-600
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_470to600/part-0.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_470to600/part-1.root",

    # 600-800
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-0.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-1.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-2.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-3.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-4.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-5.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-6.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-7.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-8.root",

    # 800-1000
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_800to1000/part-0.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_800to1000/part-1.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_800to1000/part-10.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_800to1000/part-11.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_800to1000/part-2.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_800to1000/part-3.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_800to1000/part-4.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_800to1000/part-5.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_800to1000/part-6.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_800to1000/part-7.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_800to1000/part-8.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_800to1000/part-9.root",

    # 1000-1400
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1000to1400/part-0.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1000to1400/part-1.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1000to1400/part-10.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1000to1400/part-11.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1000to1400/part-12.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1000to1400/part-2.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1000to1400/part-3.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1000to1400/part-4.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1000to1400/part-5.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1000to1400/part-6.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1000to1400/part-7.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1000to1400/part-8.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1000to1400/part-9.root",

    # 1400-1800
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-0.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-1.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-10.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-11.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-12.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-13.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-14.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-15.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-2.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-3.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-4.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-5.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-6.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-7.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-8.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-9.root",

    # 1800-2400
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-0.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-1.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-10.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-11.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-12.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-13.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-2.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-3.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-4.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-5.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-6.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-7.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-8.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-9.root",

    # 2400-3200
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_2400to3200/part-0.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_2400to3200/part-1.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_2400to3200/part-10.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_2400to3200/part-11.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_2400to3200/part-12.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_2400to3200/part-2.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_2400to3200/part-3.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_2400to3200/part-4.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_2400to3200/part-5.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_2400to3200/part-6.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_2400to3200/part-7.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_2400to3200/part-8.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_2400to3200/part-9.root",

    # 3200-Inf
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_3200toInf/part-0.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_3200toInf/part-1.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_3200toInf/part-2.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_3200toInf/part-3.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_3200toInf/part-4.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_3200toInf/part-5.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_3200toInf/part-6.root",
]

SKIM_FILES_2016 = [
    # 300-470
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_300to470/part-0.root",

    # 470-600
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_470to600/part-0.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_470to600/part-1.root",

    # 600-800
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-0.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-1.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-2.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-3.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-4.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-5.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-6.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-7.root",

    # 800-1000
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_800to1000/part-0.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_800to1000/part-1.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_800to1000/part-10.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_800to1000/part-11.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_800to1000/part-2.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_800to1000/part-3.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_800to1000/part-4.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_800to1000/part-5.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_800to1000/part-6.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_800to1000/part-7.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_800to1000/part-8.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_800to1000/part-9.root",

    # 1000-1400
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1000to1400/part-0.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1000to1400/part-1.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1000to1400/part-10.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1000to1400/part-11.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1000to1400/part-12.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1000to1400/part-2.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1000to1400/part-3.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1000to1400/part-4.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1000to1400/part-5.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1000to1400/part-6.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1000to1400/part-7.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1000to1400/part-8.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1000to1400/part-9.root",

    # 1400-1800
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-0.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-1.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-10.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-11.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-12.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-13.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-14.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-2.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-3.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-4.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-5.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-6.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-7.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-8.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1400to1800/part-9.root",

    # 1800-2400
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-0.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-1.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-10.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-11.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-12.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-2.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-3.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-4.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-5.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-6.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-7.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-8.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_1800to2400/part-9.root",

    # 2400-3200
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_2400to3200/part-0.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_2400to3200/part-1.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_2400to3200/part-10.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_2400to3200/part-11.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_2400to3200/part-2.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_2400to3200/part-3.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_2400to3200/part-4.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_2400to3200/part-5.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_2400to3200/part-6.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_2400to3200/part-7.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_2400to3200/part-8.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_2400to3200/part-9.root",

    # 3200-Inf
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_3200toInf/part-0.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_3200toInf/part-1.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_3200toInf/part-2.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_3200toInf/part-3.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_3200toInf/part-4.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_3200toInf/part-5.root",
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_3200toInf/part-6.root",
]

def build_skim_file_list_from_tchannel(sample: str, skimCut: str, startFile: int = 0, nFiles: int = -1, verbose: bool = True):
    """
    Build the skim file list exactly the way t-channel does, from input/sampleJSONs/skims/<skimCut>/...
    Returns:
      files_flat: list[str] of root:// files (all subsamples concatenated)
      fileset: dict[sampleKey -> list[str]] (kept for optional debugging)
    """
    from utils import samples as s

    fileset = s.getFileset(
        sample=sample,
        skimCut=skimCut,
        startFile=startFile,
        nFiles=nFiles,
        skimSource=True,
        verbose=verbose,
    )

    # flatten in a deterministic order
    keys_sorted = sorted(fileset.keys())
    files_flat = []
    for k in keys_sorted:
        files_flat.extend(fileset[k])

    if len(files_flat) == 0:
        raise RuntimeError(f"[FILES] empty fileset for sample='{sample}' skimCut='{skimCut}'. Check your JSONs path.")

    if verbose:
        print(f"[FILES] sample='{sample}' skimCut='{skimCut}'")
        print(f"[FILES] subsamples={len(keys_sorted)}  total_files={len(files_flat)}")
        print(f"[FILES] first_subsample='{keys_sorted[0]}'  nfiles={len(fileset[keys_sorted[0]])}")
    return files_flat, fileset

# ---------------------------
# Scaling numbers copied from utils/utility.py logic you pasted
# ---------------------------
# New numbers (use these)
__lumi_per_year = {
    "2016": 36.31e3,
    "2017": 42.07e3,
    "2018": 59.56e3,
}
# new lumi
# def lumi_and_kfactor(year: str, dataset_is_qcd: bool, hemPeriod: str):
#     """
#     Returns (luminosity, kFactor).

#     Update: lumi comes from __lumi_per_year for All-years.
#     For 2018 PreHEM/PostHEM we keep the legacy split values ONLY if you explicitly request them,
#     otherwise hemPeriod=All uses __lumi_per_year["2018"].

#     kFactor logic in utility.py :
#       - 2016 QCD: 1.2428
#       - 2017 QCD: 1.4164
#       - 2018 PreHEM QCD: 1.4623
#       - 2018 PostHEM QCD: 1.5551
#       - 2018 All: (legacy code used kFactor=1.0)
#     """
#     year = str(year)
#     hemPeriod = hemPeriod.strip()

#     if year not in __lumi_per_year:
#         raise ValueError(f"Unsupported year '{year}' (expected 2016/2017/2018)")

#     # default lumi = your new numbers
#     luminosity = float(__lumi_per_year[year])

#     # kFactors (copied from your utility.py behavior)
#     kFactor = 1.0
#     if dataset_is_qcd:
#         if year == "2016":
#             kFactor = 1.2428
#         elif year == "2017":
#             kFactor = 1.4164
#         elif year == "2018":
#             # only apply the special QCD kFactors if you are explicitly splitting by HEM
#             if hemPeriod == "PreHEM":
#                 kFactor = 1.4623
#                 # optional: if you *also* want the lumi split, uncomment next line
#                 # luminosity = 21071.460
#             elif hemPeriod == "PostHEM":
#                 kFactor = 1.5551
#                 # optional: if you *also* want the lumi split, uncomment next line
#                 # luminosity = 38621.232
#             else:
#                 # "All" case: keep kFactor=1.0 to mirror the old logic for All
#                 kFactor = 1.0

#     return float(luminosity), float(kFactor)



def lumi_and_kfactor(year: str, dataset_is_qcd: bool, hemPeriod: str):
    """
    Mimic utility.baselineVar lumi + QCD kFactor blocks.
    Returns (luminosity, kFactor).
    """
    year = str(year)
    hemPeriod = hemPeriod.strip()

    kFactor = 1.0
    if year == "2016":
        luminosity = 35921.036
        if dataset_is_qcd:
            kFactor = 1.2428
    elif year == "2017":
        luminosity = 41521.331
        if dataset_is_qcd:
            kFactor = 1.4164
    elif year == "2018":
        if hemPeriod == "PreHEM":
            luminosity = 21071.460
            if dataset_is_qcd:
                kFactor = 1.4623
        elif hemPeriod == "PostHEM":
            luminosity = 38621.232
            if dataset_is_qcd:
                kFactor = 1.5551
        else:
            #"All" case in your pasted code: lumi=59692.692, kFactor stays default 1.0
            luminosity = 59692.692
    else:
        raise ValueError(f"Unsupported year '{year}' (expected 2016/2017/2018)")
    return float(luminosity), float(kFactor)



# ---------------------------
# Helpers for ROOT reading
# ---------------------------
def find_tree(fileobj):
    keys = fileobj.keys()
    for t in TREE_CANDIDATES:
        if t in keys:
            return t
    # fallback: first TTree-like object
    for k in keys:
        try:
            obj = fileobj[k]
            if isinstance(obj, uproot.behaviors.TTree.TTree):
                return k
        except Exception:
            pass
    raise RuntimeError(f"No TTree found. Top-level keys (sample): {list(keys)[:50]}")

def choose_branch(branches, candidates, label, required=True):
    branches = list(branches)
    # exact match
    for b in candidates:
        if b in branches:
            return b
    # fuzzy
    if label.lower().startswith("met"):
        met_like = [b for b in branches if "met" in b.lower()]
        if met_like:
            for b in met_like:
                if "pt" in b.lower():
                    return b
            return met_like[0]
    if not required:
        return None
    raise RuntimeError(
        f"Could not find {label} branch. "
        f"Candidates={candidates}. Example branches={branches[:50]}"
    )

def read_th1(rootfile, histname):
    with uproot.open(rootfile) as f:
        if histname not in f:
            raise KeyError(
                f"Histogram '{histname}' not found in {rootfile}. "
                f"Available keys (sample): {list(f.keys())[:50]}"
            )
        h = f[histname]
        vals = np.asarray(h.values(flow=False), dtype=float)
        edges = np.asarray(h.axis().edges(), dtype=float)
    return vals, edges

def rebin_hist_by_center(src_vals, src_edges, target_edges):
    """Rebin a TH1-like histogram to target_edges by treating each src bin content as weight at bin center."""
    centers = 0.5 * (src_edges[:-1] + src_edges[1:])
    rebinned, _ = np.histogram(centers, bins=target_edges, weights=src_vals)
    return rebinned

# ---------------------------
# Read MET (+ weights) from skims
# ---------------------------
def read_met_and_weights_from_skims(files, max_files=None):
    met_list = []
    w_list = []
    sumw_total = 0.0
    n_total = 0

    chosen_tree = None
    b_met = None
    b_w = None
    b_pu = None

    per_file = []

    use_files = files if max_files is None else files[:max_files]

    for i, fn in enumerate(use_files, 1):
        print(f"[{i:03d}/{len(use_files):03d}] opening {fn}")
        with uproot.open(fn) as f:
            tname = find_tree(f)
            tree = f[tname]
            branches = tree.keys()

            met_branch = choose_branch(branches, MET_BRANCH_CANDIDATES, "MET", required=True)
            w_branch = choose_branch(branches, WEIGHT_BRANCH_CANDIDATES, "Weight", required=True)
            pu_branch = choose_branch(branches, PUWEIGHT_BRANCH_CANDIDATES, "puWeight", required=False)

            if chosen_tree is None:
                chosen_tree = tname
                b_met, b_w, b_pu = met_branch, w_branch, pu_branch
                print(f"      -> using tree='{chosen_tree}', MET='{b_met}', Weight='{b_w}', puWeight='{b_pu}'")

            # Read arrays (numpy)
            met = np.asarray(tree[met_branch].array(library="np")).reshape(-1)
            w = np.asarray(tree[w_branch].array(library="np")).reshape(-1)

            if pu_branch is not None:
                pu = np.asarray(tree[pu_branch].array(library="np")).reshape(-1)
            else:
                pu = np.ones_like(w)

            # Basic sanity
            if len(met) != len(w) or len(w) != len(pu):
                raise RuntimeError(
                    f"Length mismatch in {fn}: len(MET)={len(met)} len(Weight)={len(w)} len(puWeight)={len(pu)}"
                )

            met_list.append(met)
            # store un-normalized per-event base weight and puWeight separately? we'll multiply later
            w_list.append(np.stack([w, pu], axis=1))

            n_total += len(met)
            sumw_file = float(np.sum(w * pu))
            sumw_total += sumw_file
            per_file.append((fn, len(met), sumw_file))

    met_all = np.concatenate(met_list) if met_list else np.array([], dtype=float)
    wpairs = np.concatenate(w_list) if w_list else np.zeros((0, 2), dtype=float)  # columns: [Weight, puWeight]
    return met_all, wpairs, chosen_tree, (b_met, b_w, b_pu), per_file, n_total, sumw_total

# ---------------------------
# Plotting
# ---------------------------
def plot_step(ax, bins, counts, label, total):
    centers = 0.5 * (bins[:-1] + bins[1:])
    ax.step(centers, counts, where="mid", linewidth=1.8, label=f"{label} (total={total:.6g})")

def make_plot(outpath_base, bins, series, x_range, title, logy=True):
    """
    series = list of tuples: (counts, label, total)
    """
    fig = plt.figure(figsize=(10, 7))
    ax = plt.gca()

    for counts, label, total in series:
        plot_step(ax, bins, counts, label, total)
    if logy:
        ax.set_yscale("log")
    ax.set_xlim(x_range)
    width = (bins[1] - bins[0])
    ax.set_xlabel("MET [GeV]")
    ax.set_ylabel(f"Events / {width:.1f} GeV")
    ax.grid(True, which="both", linestyle=":", linewidth=0.8)
    ax.set_title(title)
    ax.legend()

    fig.tight_layout()
    fig.savefig(outpath_base + ".png")
    fig.savefig(outpath_base + ".pdf")
    plt.close(fig)

def make_simple_hist(outpath_base, values, bins, xlabel, title, logy=True):
    import numpy as np   # <-- MUST be first

    fig = plt.figure(figsize=(9, 6))
    ax = plt.gca()

    counts, edges = np.histogram(values, bins=bins)
    centers = 0.5 * (edges[:-1] + edges[1:])

    ax.step(centers, counts, where="mid", linewidth=1.8)
    if logy:
        ax.set_yscale("log")

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Entries")
    ax.set_title(title)
    ax.grid(True, which="both", linestyle=":", linewidth=0.8)

    fig.tight_layout()
    fig.savefig(outpath_base + ".png")
    fig.savefig(outpath_base + ".pdf")
    plt.close(fig)


def rebin_hist_exact(src_edges, src_counts, dst_edges):
    """
    Rebin a histogram with piecewise-constant bin contents exactly by overlap.
    Conserves total counts (integral), assuming counts are per-bin yields (not densities).

    src_edges: (N+1,)
    src_counts: (N,)
    dst_edges: (M+1,)
    returns: dst_counts (M,)
    """
    src_edges = np.asarray(src_edges, dtype=float)
    src_counts = np.asarray(src_counts, dtype=float)
    dst_edges = np.asarray(dst_edges, dtype=float)

    src_w = np.diff(src_edges)
    if np.any(src_w <= 0):
        raise ValueError("src_edges not strictly increasing")
    if np.any(np.diff(dst_edges) <= 0):
        raise ValueError("dst_edges not strictly increasing")

    # Treat each source bin as uniform density over its width
    src_density = src_counts / src_w

    dst_counts = np.zeros(len(dst_edges) - 1, dtype=float)

    i = 0  # src bin index
    for j in range(len(dst_counts)):
        lo, hi = dst_edges[j], dst_edges[j+1]
        if hi <= src_edges[0] or lo >= src_edges[-1]:
            continue

        # advance src bin pointer until it might overlap
        while i < len(src_counts) and src_edges[i+1] <= lo:
            i += 1

        k = i
        while k < len(src_counts) and src_edges[k] < hi:
            overlap_lo = max(lo, src_edges[k])
            overlap_hi = min(hi, src_edges[k+1])
            if overlap_hi > overlap_lo:
                dst_counts[j] += src_density[k] * (overlap_hi - overlap_lo)
            if src_edges[k+1] >= hi:
                break
            k += 1

    return dst_counts

def read_met_and_weights_by_subsample(fileset_dict):
    """
    Reads files grouped by their keys (subsamples) to track contributions.
    """
    subsample_data = {}
    
    for label, files in fileset_dict.items():
        print(f"--- Processing Subsample: {label} ({len(files)} files) ---")
        met_list = []
        w_list = []
        
        # Reuse your existing logic but scoped to this label
        for fn in files:
            with uproot.open(fn) as f:
                tname = find_tree(f)
                tree = f[tname]
                branches = tree.keys()
                
                met_branch = choose_branch(branches, MET_BRANCH_CANDIDATES, "MET")
                w_branch = choose_branch(branches, WEIGHT_BRANCH_CANDIDATES, "Weight")
                pu_branch = choose_branch(branches, PUWEIGHT_BRANCH_CANDIDATES, "puWeight", required=False)
                
                met = np.asarray(tree[met_branch].array(library="np")).reshape(-1)
                w = np.asarray(tree[w_branch].array(library="np")).reshape(-1)
                pu = np.asarray(tree[pu_branch].array(library="np")).reshape(-1) if pu_branch else np.ones_like(w)
                
                met_list.append(met)
                w_list.append(np.stack([w, pu], axis=1))
        
        subsample_data[label] = {
            "met": np.concatenate(met_list) if met_list else np.array([]),
            "wpairs": np.concatenate(w_list) if w_list else np.zeros((0, 2))
        }
    return subsample_data
def build_fileset_dict(file_list):
    """
    Groups a flat list of files into a dictionary keyed by the pT bin name
    found in the directory path (e.g., 'QCD_Pt_600to800').
    """
    from collections import defaultdict
    d = defaultdict(list)
    for f in file_list:
        # Extract the directory name that contains the pT range
        # Path format: .../nominal//QCD_Pt_XXXtoYYY/part-N.root
        parts = f.split('/')
        pt_bin = "Unknown"
        for p in parts:
            if "QCD_Pt_" in p:
                pt_bin = p
                break
        d[pt_bin].append(f)
    return dict(d)

    
# ---------------------------
# Main
# ---------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tchan-root", required=True)
    ap.add_argument("--tchan-hist", default="h_MET_pre")

    ap.add_argument("--sample", required=True)
    ap.add_argument("--skimCut", required=True)
    ap.add_argument("--startFile", type=int, default=0)
    ap.add_argument("--nFiles", type=int, default=-1)

    ap.add_argument("--outdir", default="./ABCD_hists3/met_compare")
    ap.add_argument("--nbins", type=int, default=120)
    ap.add_argument("--xmax", type=float, default=2000.0)
    ap.add_argument("--xmin", type=float, default=0.0)

    ap.add_argument("--year", required=True)
    ap.add_argument("--hemPeriod", default="All", choices=["All", "PreHEM", "PostHEM"])

    ap.add_argument("--scaleFactor", type=float, default=1.0)
    ap.add_argument("--applyPU", action="store_true")

    ap.add_argument("--maxFiles", type=int, default=None)

    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Pick the correct list based on year
    full_list = SKIM_FILES_2018 if args.year == "2018" else SKIM_FILES_2016
    
    # NEW: Create the dictionary required by read_met_and_weights_by_subsample
    fileset_dict = build_fileset_dict(full_list)

    if args.maxFiles is not None and args.maxFiles > 0:
        print(f"Limiting to {args.maxFiles} files per pT bin...")
        for k in fileset_dict:
            fileset_dict[k] = fileset_dict[k][:args.maxFiles]

    # Histogram settings
    bins = np.linspace(args.xmin, args.xmax, args.nbins + 1)
    
    # --------------------------------------------------
    # 1. Read data grouped by pT bin
    # --------------------------------------------------
    # This calls the new function you added previously
    subsample_results = read_met_and_weights_by_subsample(fileset_dict)

    # Prepare for plotting
    all_raw_series = []
    all_scaled_series = []
    lumi, kfactor = lumi_and_kfactor(args.year, True, args.hemPeriod)

    # --------------------------------------------------
    # 2. Process each bin and prepare for breakdown plots
    # --------------------------------------------------
    for label, data in subsample_results.items():
        met = data["met"]
        w_base = data["wpairs"][:, 0]
        pu = data["wpairs"][:, 1]
        
        # Calculate event weights
        if args.applyPU:
            w_evt = lumi * w_base * args.scaleFactor * kfactor * pu
        else:
            w_evt = lumi * w_base * args.scaleFactor * kfactor
            
        # Create Histograms for this specific pT bin
        raw_counts, _ = np.histogram(met, bins=bins)
        scaled_counts, _ = np.histogram(met, bins=bins, weights=w_evt)
        
        all_raw_series.append((raw_counts, label, float(np.sum(raw_counts))))
        all_scaled_series.append((scaled_counts, label, float(np.sum(scaled_counts))))

    # --------------------------------------------------
    # 3. Generate Breakdown Plots
    # --------------------------------------------------
    # Plot: Raw Breakdown
    make_plot(
        os.path.join(args.outdir, "04_raw_breakdown_by_pt"),
        bins,
        all_raw_series,
        (args.xmin, args.xmax),
        f"Raw MET Breakdown by pT bin ({args.year})",
        logy=True
    )

    # Plot: Scaled Breakdown
    make_plot(
        os.path.join(args.outdir, "05_scaled_breakdown_by_pt"),
        bins,
        all_scaled_series,
        (args.xmin, args.xmax),
        f"Scaled MET Breakdown by pT bin ({args.year})",
        logy=True
    )

    print(f"Done! Check output in: {args.outdir}")

if __name__ == "__main__":
    main()