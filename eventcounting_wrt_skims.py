#!/usr/bin/env python3
"""
Reproduce MET plots from skims and compare to t-channel output.

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
SKIM_FILES = [
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

# ---------------------------
# Scaling numbers copied from utils/utility.py logic you pasted
# ---------------------------
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
            # "All" case in your pasted code: lumi=59692.692, kFactor stays default 1.0
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


    import numpy as np

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

# ---------------------------
# Main
# ---------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tchan-root", required=True, help="Path to t-channel output ROOT file (e.g. 2018_QCD.root)")
    ap.add_argument("--tchan-hist", default="h_MET_pre", help="Histogram name inside t-channel root")
    ap.add_argument("--outdir", default="./ABCD_hists3/met_compare_2018QCD", help="Output directory")
    ap.add_argument("--nbins", type=int, default=120, help="Overlay binning")
    ap.add_argument("--xmax", type=float, default=4000.0, help="Max MET for plots")
    ap.add_argument("--xmin", type=float, default=100.0, help="Min MET for plots")
    ap.add_argument("--year", default="2018", help="Year string: 2016/2017/2018")
    ap.add_argument("--hemPeriod", default="All", choices=["All", "PreHEM", "PostHEM"], help="2018 HEM period choice")
    ap.add_argument("--scaleFactor", type=float, default=1.0, help="Extra multiplicative scaleFactor (default 1)")
    ap.add_argument("--applyPU", action="store_true", help="Multiply by puWeight (recommended for bkg; matches utility.py)")
    ap.add_argument("--maxFiles", type=int, default=None, help="If set, only process first N skim files")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Build bins
    bins = np.linspace(args.xmin, args.xmax, args.nbins + 1)

    # ---- Read skims
    met, wpairs, tname, (b_met, b_w, b_pu), per_file, n_total, sumw_base = read_met_and_weights_from_skims(
        SKIM_FILES, max_files=args.maxFiles
    )
    w_base = wpairs[:, 0]  # events.Weight
    pu = wpairs[:, 1]      # events.puWeight or ones

    print(f"\n[SKIMS] total events read = {n_total}")
    print(f"[SKIMS] tree='{tname}', MET='{b_met}', Weight='{b_w}', puWeight='{b_pu}'")
    print(f"[SKIMS] sum(Weight*puWeight) across all files = {sumw_base:.6g}")

    # Save MET txt
    # met_txt = os.path.join(args.outdir, "met_from_skims_raw.txt")
    # np.savetxt(met_txt, met, fmt="%.6f")
    # print(f"[OK] wrote {met_txt}")

    # Save base weights txt (two columns)
    # w_txt = os.path.join(args.outdir, "met_from_skims_w.txt")
    # np.savetxt(w_txt, np.column_stack([w_base, pu]), fmt="%.8e", header="Weight puWeight")
    # print(f"[OK] wrote {w_txt}")

    # Summary txt
    summary_txt = os.path.join(args.outdir, "skim_summary.txt")
    with open(summary_txt, "w") as fp:
        fp.write(f"tree={tname}\n")
        fp.write(f"MET={b_met}\n")
        fp.write(f"Weight={b_w}\n")
        fp.write(f"puWeight={b_pu}\n")
        fp.write(f"total_events={n_total}\n")
        fp.write(f"sum_weight_times_pu={sumw_base:.9g}\n")
        fp.write("\nPer-file:\n")
        for fn, n, sw in per_file:
            fp.write(f"{n:10d}  {sw:15.9g}  {fn}\n")
    print(f"[OK] wrote {summary_txt}")

    # ---- Build RAW skim histogram (unweighted)
    raw_counts, _ = np.histogram(met, bins=bins)
    raw_total = float(np.sum(raw_counts))

    # ---- Build SCALED skim histogram (weighted like utility.py)
    lumi, kfactor = lumi_and_kfactor(args.year, dataset_is_qcd=True, hemPeriod=args.hemPeriod)
    if args.applyPU:
        w_evt = lumi * w_base * args.scaleFactor * kfactor * pu
    else:
        w_evt = lumi * w_base * args.scaleFactor * kfactor

    scaled_counts, _ = np.histogram(met, bins=bins, weights=w_evt)
    scaled_total = float(np.sum(scaled_counts))

    print(f"\n[SCALE] year={args.year} hemPeriod={args.hemPeriod}")
    print(f"[SCALE] lumi={lumi:.6f}  kFactor={kfactor:.6f}  scaleFactor={args.scaleFactor:.6g}  applyPU={args.applyPU}")
    print(f"[SCALE] total (sum of weights in range) = {scaled_total:.6g}")

    # ---- Read t-channel output histogram (TH1)
    tc_vals, tc_edges = read_th1(args.tchan_root, args.tchan_hist)

    # dump axis info (what you asked for)
    print(f"\n[TCHAN AXIS] hist='{args.tchan_hist}'")
    print(f"[TCHAN AXIS] n_edges={len(tc_edges)}  nbins={len(tc_edges)-1}")
    print(f"[TCHAN AXIS] xmin={float(tc_edges[0])}  xmax={float(tc_edges[-1])}")
    print(f"[TCHAN AXIS] first10={tc_edges[:10]}")
    print(f"[TCHAN AXIS] last10 ={tc_edges[-10:]}")

    # Rebin to our bins if needed (EXACT overlap-based rebin)
    same_edges = (len(tc_edges) == len(bins)) and np.allclose(tc_edges, bins)
    if same_edges:
        tc_reb = tc_vals.copy()
        print(f"[TCHAN] '{args.tchan_hist}' edges match overlay binning.")
    else:
        tc_reb = rebin_hist_exact(tc_edges, tc_vals, bins)
        print(f"[TCHAN] '{args.tchan_hist}' edges differ; rebinned EXACTLY by bin-overlap.")

        # sanity check: conserve integral in [xmin,xmax] overlap
        # (full-integral conservation holds if bins fully cover src range; here we usually do)
        print(f"[TCHAN] integral(original) = {float(np.sum(tc_vals)):.6g}")
        print(f"[TCHAN] integral(rebinned) = {float(np.sum(tc_reb)):.6g}")

    tc_total = float(np.sum(tc_reb))
    print(f"[TCHAN] integral(used for plots) = {tc_total:.6g}")

    # ---- Make 3 plots
    x_range = (args.xmin, args.xmax)

    # (1) raw skims + tchan overlay
    make_plot(
        outpath_base=os.path.join(args.outdir, "01_raw_skims_overlay_tchannel"),
        bins=bins,
        series=[
            (raw_counts, f"Skims MET raw (N={n_total})", raw_total),
            (tc_reb, f"t-channel {args.tchan_hist}", tc_total),
        ],
        x_range=x_range,
        title=f"MET: raw skims vs t-channel ({args.year} QCD)",
        logy=True,
    )

    # (2) scaled skims only
    make_plot(
        outpath_base=os.path.join(args.outdir, "02_scaled_skims_only"),
        bins=bins,
        series=[
            (scaled_counts, f"Skims MET scaled (total={scaled_total:.0f})", scaled_total),
        ],
        x_range=x_range,
        title=f"MET: scaled skims only ({args.year} QCD, hem={args.hemPeriod})",
        logy=True,
    )

    # (3) scaled skims + tchan overlay
    make_plot(
        outpath_base=os.path.join(args.outdir, "03_scaled_skims_overlay_tchannel"),
        bins=bins,
        series=[
            (scaled_counts, f"Skims MET scaled (hem={args.hemPeriod}, total={scaled_total:.0f})", scaled_total),
            (tc_reb, f"t-channel {args.tchan_hist} (total={tc_total:.0f})", tc_total),
        ],
        x_range=x_range,
        title=f"MET: scaled skims vs t-channel ({args.year} QCD, hem={args.hemPeriod})",
        logy=True,
    )

    print("\n[OK] Done. Outputs written to:")
    print(f"  {args.outdir}")
    print("  - 01_raw_skims_overlay_tchannel.(png/pdf)")
    print("  - 02_scaled_skims_only.(png/pdf)")
    print("  - 03_scaled_skims_overlay_tchannel.(png/pdf)")
    print("  - met_from_skims_raw.txt, met_from_skims_w.txt, skim_summary.txt")
    print("\nTip: If your t-channel ROOT is from PreHEM/PostHEM, re-run with --hemPeriod PreHEM or PostHEM for closer agreement.")

if __name__ == "__main__":
    main()
