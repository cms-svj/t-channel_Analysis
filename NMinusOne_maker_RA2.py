#!/usr/bin/env python3
"""
Build N-1 plots from RA2AnalysisTree ntuples with a cache-first workflow.

The structure mirrors the newer plotting scheme:
  * stream TreeMaker2/PreSelection once from EOS;
  * write raw and normalized TH1D objects to a compact ROOT cache;
  * draw PDF/PNG plots from that cache;
  * rerun with --plot-only to adjust plot style without rereading EOS.

Default input:
  /store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Summer20UL18/
      QCD_Pt_300to470_TuneCP5_13TeV_pythia8

source condor/initCondor.sh

mkdir -p supplementry_plots_postTrimandgapveto/Nminus1plots/FINALPLOTSssssNMINUS1 logs

python3 -u NMinusOne_maker_RA2.py \
  --input-dir /store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/Summer20UL18/QCD_Pt_300to470_TuneCP5_13TeV_pythia8 \
  --data-dir /store/user/.../Data \
  --output supplementry_plots_postTrimandgapveto/Nminus1plots/nminusone_qcd_pt300to470.root \
  --plot-dir supplementry_plots_postTrimandgapveto/Nminus1plots/FINALPLOTSssssNMINUS1 \
  --chunk-size "100 MB" \
  |& tee "logs/Nminus1_QCD_Pt300to470_$(date +%Y%m%d_%H%M%S).log"

The optional --data-dir can point at a directory containing RA2AnalysisTree
ROOT files directly, or one level of sample subdirectories. Data is stored in
raw/data/observed and overlaid as markers when plotting.

to redraw
cd /uscms/home/ashrivas/nobackup/Dark_Sector/t-channel_Analysis
source condor/initCondor.sh

python3 -u NMinusOne_maker_RA2.py \
  --plot-only \
  --output supplementry_plots_postTrimandgapveto/Nminus1plots/nminus1_RA2_allbkgs_Run2_withSignals.root \
  --plot-dir supplementry_plots_postTrimandgapveto/Nminus1plots/FINALPLOTSssssNMINUS1
"""

from __future__ import annotations

import argparse
import math
import os
import re
import subprocess
import sys
from array import array
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import awkward as ak
import numpy as np
import uproot

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2(True)
ROOT.TH1.AddDirectory(False)


EOS_HOST = "root://cmseos.fnal.gov"
DEFAULT_INPUT_DIR = (
    "/store/user/lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/"
    "Summer20UL18/QCD_Pt_300to470_TuneCP5_13TeV_pythia8"
)

DEFAULT_OUTPUT = "NMinusOne_RA2_plots/nminusone_qcd_pt300to470.root"
DEFAULT_PLOT_DIR = "supplementry_plots_postTrimandgapveto/Nminus1plots/FINALPLOTSssssNMINUS1"
DEFAULT_DATA_DIR = ""
DEFAULT_SIGNAL_BASE = "/store/user/lpcdarkqcd/tchannel_UL"
DEFAULT_SIGNAL_YEARS = "2016,2016APV"
DEFAULT_SIGNAL_FILE_REGEX = r".*_RA2AnalysisTree\.root$"

DEFAULT_SIGNAL_DIRS = [
    "t-channel_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
    "t-channel_mMed-2000_mDark-20_rinv-0p1_alpha-peak_yukawa-1",
    "t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
    # "t-channel_mMed-2000_mDark-20_rinv-0p7_alpha-peak_yukawa-1",
    # "t-channel_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
]
DISABLED_SIGNAL_TOKENS = [
    "t-channel_mMed-2000_mDark-20_rinv-0p7_alpha-peak_yukawa-1",
    "t-channel_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
]

LUMI_PB = 59692.692

PROC_LABEL = {
    "QCD": "QCD",
    "TTJets": "t#bar{t}+jets",
    "WJetsToLNu": "W+jets",
    "ZJetsToNuNu": "Z#rightarrow#nu#nu+jets",
    "ST": "Single top",
}

PROC_COLOR = {
    "QCD": ROOT.TColor.GetColor("#9c9ca1"),
    "TTJets": ROOT.TColor.GetColor("#7a21dd"),
    "WJetsToLNu": ROOT.TColor.GetColor("#e42536"),
    "ZJetsToNuNu": ROOT.TColor.GetColor("#f89c20"),
    "ST": ROOT.TColor.GetColor("#5790fc"),
}

SIGNAL_LINE_COLORS = [
    ROOT.TColor.GetColor("#92dadd"),
    ROOT.TColor.GetColor("#6b3e26"),
    ROOT.TColor.GetColor("#0b3d02"),
    ROOT.TColor.GetColor("#228833"),
    ROOT.TColor.GetColor("#1f4e79"),
]
SIGNAL_LINE_STYLES = [1, 3, 1, 2, 1]
DISCRETE_AXIS_VARIABLES = {"njetsAK8"}

STACK_ORDER = ["ST", "ZJetsToNuNu", "WJetsToLNu", "TTJets", "QCD"]
LEGEND_ORDER = ["QCD", "TTJets", "WJetsToLNu", "ZJetsToNuNu", "ST"]

VARIABLES = {
    "MET": {
        "title": "p_{T}^{miss} [GeV]",
        "nbins": 500,
        "xmin": 0.0,
        "xmax": 2000.0,
        "display_xmax": 1000.0,
        "line": 200.0,
        "line_side": "right",
        "arrow_y_fraction": 0.66,
    },
    "ST": {
        "title": "S_{T} [GeV]",
        "nbins": 500,
        "xmin": 0.0,
        "xmax": 5000.0,
        "line": 1300.0,
        "line_side": "right",
    },
    "dPhiMinjMETAK8": {
        "title": "#Delta#phi_{min}(J, p_{T}^{miss})",
        "nbins": 100,
        "xmin": 0.0,
        "xmax": 3.5,
        "line": 1.5,
        "line_side": "left",
    },
    "j1PtAK8": {
        "title": "p_{T}(J_{1}) [GeV]",
        "nbins": 400,
        "xmin": 0.0,
        "xmax": 4000.0,
    },
    "j2PtAK8": {
        "title": "p_{T}(J_{2}) [GeV]",
        "nbins": 400,
        "xmin": 0.0,
        "xmax": 4000.0,
    },
    "jPtAK8": {
        "title": "p_{T}(J) [GeV]",
        "nbins": 280,
        "xmin": 0.0,
        "xmax": 2800.0,
    },
    "nelectrons": {
        "title": "Number of Electrons",
        "nbins": 6,
        "xmin": -0.5,
        "xmax": 5.5,
        "line": 0.5,
        "line_side": "left",
    },
    "nmuons": {
        "title": "Number of Muons",
        "nbins": 6,
        "xmin": -0.5,
        "xmax": 5.5,
        "line": 0.5,
        "line_side": "left",
    },
    "nl": {
        "title": "Number of Leptons",
        "nbins": 6,
        "xmin": -0.5,
        "xmax": 5.5,
        "line": 0.5,
        "line_side": "left",
    },
    "njetsAK8": {
        "title": "Number of AK8 Jets",
        "nbins": 16,
        "xmin": -0.5,
        "xmax": 15.5,
        "display_xmin": 0.0,
        "display_xmax": 10.0,
        "line": 2.0,
        "line_side": "right",
    },
}

PLOTS = [
    ("h_MET_pre__metcut", "MET", "metcut"),
    ("h_ST_pre__stcut", "ST", "stcut"),
    ("h_dPhiMinjMETAK8_pre__dPhiMin", "dPhiMinjMETAK8", "dPhiMin"),
    ("h_j1PtAK8_pre__stcut", "j1PtAK8", "stcut"),
    ("h_j2PtAK8_pre__stcut", "j2PtAK8", "stcut"),
    ("h_jPtAK8_pre__stcut", "jPtAK8", "stcut"),
    ("h_nelectrons_pre__nl", "nelectrons", "nl"),
    ("h_nmuons_pre__nl", "nmuons", "nl"),
    ("h_nl_pre__nl", "nl", "nl"),
    ("h_njetsAK8_pre__2jetsAK8", "njetsAK8", "2jetsAK8"),
]

PLOT_TARGET_BINS = 40

QCD_MET_DISPLAY_SPIKE_MIN_X = 600.0
QCD_MET_DISPLAY_SPIKE_RATIO = 4.0
DPHI_DISPLAY_SPIKE_RATIO = 2.0


@dataclass
class HistogramAccumulator:
    edges: np.ndarray
    sumw: np.ndarray = field(init=False)
    sumw2: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        self.sumw = np.zeros(len(self.edges) - 1, dtype=np.float64)
        self.sumw2 = np.zeros(len(self.edges) - 1, dtype=np.float64)

    @classmethod
    def for_variable(cls, variable: str) -> "HistogramAccumulator":
        cfg = VARIABLES[variable]
        edges = np.linspace(cfg["xmin"], cfg["xmax"], cfg["nbins"] + 1)
        return cls(edges=edges)

    def fill(self, values, weights) -> None:
        values = np.asarray(ak.to_numpy(ak.fill_none(values, np.nan)), dtype=np.float64)
        weights = np.asarray(ak.to_numpy(ak.fill_none(weights, np.nan)), dtype=np.float64)
        mask = (
            np.isfinite(values)
            & np.isfinite(weights)
            & (values >= self.edges[0])
            & (values < self.edges[-1])
        )
        if not np.any(mask):
            return
        values = values[mask]
        weights = weights[mask]
        self.sumw += np.histogram(values, bins=self.edges, weights=weights)[0]
        self.sumw2 += np.histogram(values, bins=self.edges, weights=weights * weights)[0]

    def integral(self) -> float:
        return float(np.sum(self.sumw))


def empty_accumulators() -> Dict[str, HistogramAccumulator]:
    return {plot_name: HistogramAccumulator.for_variable(variable) for plot_name, variable, _ in PLOTS}


def add_accumulators(target, source) -> None:
    for plot_name, _, _ in PLOTS:
        target[plot_name].sumw += source[plot_name].sumw
        target[plot_name].sumw2 += source[plot_name].sumw2


def scale_accumulators(accumulators, scale: float) -> None:
    for plot_name, _, _ in PLOTS:
        accumulators[plot_name].sumw *= scale
        accumulators[plot_name].sumw2 *= scale * scale


def normalize_eos_path(path: str) -> str:
    path = os.path.expanduser(path).rstrip("/")
    if path.startswith("/eos/uscms/"):
        path = path[len("/eos/uscms") :]
    if not path.startswith("/"):
        raise ValueError(f"Expected an absolute EOS path, received: {path}")
    return path


def xrdfs_ls(path: str) -> List[str]:
    command = ["xrdfs", EOS_HOST, "ls", normalize_eos_path(path)]
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"xrdfs failed for {path}: {result.stderr.strip()}")
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def eos_url(path: str) -> str:
    return EOS_HOST + "//" + normalize_eos_path(path).lstrip("/")


def get_input_files(input_dir: str, regex: str, max_files: int) -> List[str]:
    pattern = re.compile(regex)
    files = sorted(path for path in xrdfs_ls(input_dir) if pattern.search(os.path.basename(path)))
    return files[:max_files] if max_files > 0 else files


def get_input_files_with_subdirs(input_dir: str, regex: str, max_files: int) -> List[str]:
    pattern = re.compile(regex)
    entries = sorted(xrdfs_ls(input_dir))
    files = [path for path in entries if pattern.search(os.path.basename(path))]
    subdirs = [path for path in entries if path not in files]
    for subdir in subdirs:
        try:
            files.extend(
                path
                for path in sorted(xrdfs_ls(subdir))
                if pattern.search(os.path.basename(path))
            )
        except RuntimeError as exc:
            print(f"[WARN] Could not list data subdir {subdir}: {exc}")
    files = list(dict.fromkeys(files))
    return files[:max_files] if max_files > 0 else files


def tree_path(root_file) -> str:
    if "Events" in root_file:
        return "Events"
    if "TreeMaker2/PreSelection" in root_file:
        return "TreeMaker2/PreSelection"
    if "TreeMaker2;1/PreSelection;1" in root_file:
        return "TreeMaker2;1/PreSelection;1"
    raise KeyError("Could not find Events or TreeMaker2/PreSelection")


BRANCH_ALIASES = {
    "MET": ["MET"],
    "HT": ["HT"],
    "DeltaPhiMin_AK8": ["DeltaPhiMin_AK8"],
    "ak8_pt": ["JetsAK8/JetsAK8.fCoordinates.fPt", "JetsAK8_/.fPt"],
    "ak8_eta": ["JetsAK8/JetsAK8.fCoordinates.fEta", "JetsAK8_/.fEta"],
    "ak8_id": ["JetsAK8_ID"],
    "ak8_is_good": ["JetsAK8_isGood"],
    "ele_pt": ["Electrons/Electrons.fCoordinates.fPt", "Electrons_/.fPt"],
    "ele_eta": ["Electrons/Electrons.fCoordinates.fEta", "Electrons_/.fEta"],
    "ele_iso": ["Electrons_iso"],
    "mu_pt": ["Muons/Muons.fCoordinates.fPt", "Muons_/.fPt"],
    "mu_eta": ["Muons/Muons.fCoordinates.fEta", "Muons_/.fEta"],
    "mu_iso": ["Muons_iso"],
    "Weight": ["Weight"],
    "puWeight": ["puWeight"],
    "NonPrefiringProb": ["NonPrefiringProb"],
    "lundWeightNom": ["lundWeightNom"],
}


REQUIRED_BRANCH_KEYS = [
    "MET",
    "HT",
    "DeltaPhiMin_AK8",
    "ak8_pt",
    "ak8_eta",
    "ak8_id",
    "ele_pt",
    "ele_eta",
    "ele_iso",
    "mu_pt",
    "mu_eta",
    "mu_iso",
]


def select_branch(available: set, logical_name: str, required: bool = True) -> Optional[str]:
    for candidate in BRANCH_ALIASES[logical_name]:
        if candidate in available:
            return candidate
    if required:
        raise KeyError(
            f"Missing required branch for {logical_name}; tried {BRANCH_ALIASES[logical_name]}"
        )
    return None


def branch_map_for_tree(tree, is_data: bool = False) -> Dict[str, str]:
    available = set(tree.keys())
    mapping: Dict[str, str] = {}
    for logical_name in REQUIRED_BRANCH_KEYS:
        mapping[logical_name] = select_branch(available, logical_name)

    if not is_data:
        mapping["Weight"] = select_branch(available, "Weight")
        for logical_name in ("puWeight", "NonPrefiringProb", "lundWeightNom"):
            branch = select_branch(available, logical_name, required=False)
            if branch:
                mapping[logical_name] = branch

    branch = select_branch(available, "ak8_is_good", required=False)
    if branch:
        mapping["ak8_is_good"] = branch
    return mapping


def remap_arrays(arrays, mapping: Dict[str, str]):
    return {logical_name: arrays[branch] for logical_name, branch in mapping.items()}


def common_branches() -> List[str]:
    return [BRANCH_ALIASES[name][0] for name in REQUIRED_BRANCH_KEYS]


def required_branches(tree) -> List[str]:
    mapping = branch_map_for_tree(tree, is_data=False)
    return list(dict.fromkeys(mapping.values()))


def required_data_branches(tree) -> List[str]:
    mapping = branch_map_for_tree(tree, is_data=True)
    return list(dict.fromkeys(mapping.values()))


def good_count(pt, eta, obj_id=None, min_pt=50.0, max_eta=2.4):
    mask = (pt > min_pt) & (abs(eta) < max_eta)
    if obj_id is not None:
        mask = mask & (obj_id == True)
    return mask, ak.num(pt[mask], axis=1)


def compute_event_variables(arrays):
    met = arrays["MET"]
    ht = arrays["HT"]
    st = ht + met
    dphi_min = arrays["DeltaPhiMin_AK8"]

    ak8_pt = arrays["ak8_pt"]
    ak8_eta = arrays["ak8_eta"]
    ak8_id = arrays["ak8_id"]
    good_ak8_mask, njets_ak8 = good_count(ak8_pt, ak8_eta, ak8_id, min_pt=50.0, max_eta=2.4)
    if "ak8_is_good" in arrays:
        good_ak8_mask = good_ak8_mask & (arrays["ak8_is_good"] == True)
        njets_ak8 = ak.num(ak8_pt[good_ak8_mask], axis=1)
    good_ak8_pt = ak8_pt[good_ak8_mask]
    padded_ak8_pt = ak.pad_none(good_ak8_pt, 2, axis=1)
    j1_pt = ak.fill_none(padded_ak8_pt[:, 0], np.nan)
    j2_pt = ak.fill_none(padded_ak8_pt[:, 1], np.nan)

    ele_pt = arrays["ele_pt"]
    ele_eta = arrays["ele_eta"]
    ele_iso = arrays["ele_iso"]
    good_ele = (ele_pt > 10.0) & (abs(ele_eta) < 2.4) & (abs(ele_iso) < 0.1)
    nelectrons = ak.num(ele_pt[good_ele], axis=1)

    mu_pt = arrays["mu_pt"]
    mu_eta = arrays["mu_eta"]
    mu_iso = arrays["mu_iso"]
    good_mu = (mu_pt > 10.0) & (abs(mu_eta) < 2.4) & (abs(mu_iso) < 0.4)
    nmuons = ak.num(mu_pt[good_mu], axis=1)
    nl = nelectrons + nmuons

    return {
        "MET": met,
        "ST": st,
        "dPhiMinjMETAK8": dphi_min,
        "j1PtAK8": j1_pt,
        "j2PtAK8": j2_pt,
        "jPtAK8": good_ak8_pt,
        "nelectrons": nelectrons,
        "nmuons": nmuons,
        "nl": nl,
        "njetsAK8": njets_ak8,
    }


def cut_masks(values) -> Dict[str, ak.Array]:
    return {
        "metcut": values["MET"] > 200.0,
        "stcut": values["ST"] > 1300.0,
        "dPhiMin": values["dPhiMinjMETAK8"] <= 1.5,
        "2jetsAK8": values["njetsAK8"] >= 2,
        "nl": values["nl"] == 0,
    }


def nminusone_mask(cuts: Dict[str, ak.Array], omitted: str):
    mask = None
    for cut_name, cut in cuts.items():
        if cut_name == omitted:
            continue
        mask = cut if mask is None else (mask & cut)
    return mask


def lund_nominal_norm(root_file) -> Optional[float]:
    if "CutFlow" not in root_file:
        return None
    cutflow = root_file["CutFlow"]
    available = set(cutflow.keys())
    if "Initial" not in available or "InitialLundNominal" not in available:
        return None
    initial = float(cutflow["Initial"].array(library="np")[0])
    initial_lund = float(cutflow["InitialLundNominal"].array(library="np")[0])
    if not np.isfinite(initial) or not np.isfinite(initial_lund) or initial_lund == 0.0:
        return None
    return initial / initial_lund


def stream_file(file_url: str, args, accumulators, is_data: bool = False, apply_lund: bool = False) -> int:
    with uproot.open(file_url) as root_file:
        tree = root_file[tree_path(root_file)]
        mapping = branch_map_for_tree(tree, is_data=is_data)
        branches = list(dict.fromkeys(mapping.values()))
        entries = int(tree.num_entries)
        lund_norm = None
        if apply_lund and not is_data and "lundWeightNom" not in mapping:
            raise RuntimeError(f"Requested Lund correction, but lundWeightNom is missing in {file_url}")
        if apply_lund and not is_data and "lundWeightNom" in mapping:
            lund_norm = lund_nominal_norm(root_file)
            if lund_norm is None:
                raise RuntimeError(
                    f"Requested Lund correction, but CutFlow/Initial and "
                    f"CutFlow/InitialLundNominal are missing or invalid in {file_url}"
                )

        for arrays in tree.iterate(expressions=branches, step_size=args.chunk_size, library="ak"):
            arrays = remap_arrays(arrays, mapping)
            n_events = len(arrays["MET"])
            if is_data:
                weights = np.ones(n_events, dtype=np.float64)
            else:
                weights = arrays["Weight"] * float(args.lumi_pb) * float(args.kfactor)
                if args.apply_pu:
                    if "puWeight" in arrays:
                        weights = weights * arrays["puWeight"]
                    if "NonPrefiringProb" in arrays:
                        weights = weights * arrays["NonPrefiringProb"]
                if lund_norm is not None and "lundWeightNom" in arrays:
                    weights = weights * arrays["lundWeightNom"] * lund_norm

            values = compute_event_variables(arrays)
            cuts = cut_masks(values)

            for plot_name, variable, omitted_cut in PLOTS:
                event_mask = nminusone_mask(cuts, omitted_cut)
                selected_values = values[variable][event_mask]
                selected_weights = weights[event_mask]
                if variable == "jPtAK8":
                    selected_weights = ak.flatten(ak.broadcast_arrays(selected_values, selected_weights)[1])
                    selected_values = ak.flatten(selected_values)
                accumulators[plot_name].fill(selected_values, selected_weights)
    return entries


def mkdirs(root_file: ROOT.TFile, directory_path: str):
    directory = root_file
    for part in (item for item in directory_path.split("/") if item):
        next_dir = directory.GetDirectory(part)
        if not next_dir:
            next_dir = directory.mkdir(part)
        directory = next_dir
    return directory


def accumulator_to_histogram(plot_name: str, accumulator: HistogramAccumulator, y_title: str = "Events"):
    variable = dict((name, var) for name, var, _ in PLOTS)[plot_name]
    cfg = VARIABLES[variable]
    hist = ROOT.TH1D(plot_name, f"{plot_name};{cfg['title']};{y_title}", len(accumulator.edges) - 1, array("d", accumulator.edges.tolist()))
    hist.Sumw2()
    for ibin, (content, variance) in enumerate(zip(accumulator.sumw, accumulator.sumw2), start=1):
        hist.SetBinContent(ibin, float(content))
        hist.SetBinError(ibin, math.sqrt(max(0.0, float(variance))))
    hist.SetDirectory(0)
    return hist


def write_hist(root_file, directory_path: str, hist, key_name: str) -> None:
    directory = mkdirs(root_file, directory_path)
    directory.cd()
    hist.Write(key_name, ROOT.TObject.kOverwrite)


def write_signal_hists(root_file, signal_accumulators) -> None:
    if not signal_accumulators:
        return
    for signal_name, accumulators in sorted(signal_accumulators.items()):
        for plot_name, _, _ in PLOTS:
            raw = accumulator_to_histogram(plot_name, accumulators[plot_name], y_title="Events")
            write_hist(root_file, f"raw/signals/{signal_name}", raw, plot_name)

            shape = raw.Clone(f"{plot_name}_{signal_name}_shape")
            shape.SetDirectory(0)
            integral = shape.Integral()
            if integral > 0.0:
                shape.Scale(1.0 / integral)
            shape.GetYaxis().SetTitle("Arbitrary units")
            write_hist(root_file, f"shapes/signals/{signal_name}", shape, plot_name)


def write_cache(
    output_path: str,
    accumulators,
    metadata: str,
    process: str = "QCD",
    data_accumulators=None,
    signal_name: Optional[str] = None,
) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    output = ROOT.TFile.Open(output_path, "RECREATE")
    if not output or output.IsZombie():
        raise RuntimeError(f"Could not create ROOT cache: {output_path}")
    try:
        if signal_name:
            write_signal_hists(output, {signal_name: accumulators})
        else:
            for plot_name, _, _ in PLOTS:
                raw = accumulator_to_histogram(plot_name, accumulators[plot_name], y_title="Events")
                write_hist(output, f"raw/{process}", raw, plot_name)

                shape = raw.Clone(f"{plot_name}_shape")
                shape.SetDirectory(0)
                integral = shape.Integral()
                if integral > 0.0:
                    shape.Scale(1.0 / integral)
                shape.GetYaxis().SetTitle("Arbitrary units")
                write_hist(output, f"shapes/{process}", shape, plot_name)

        if data_accumulators:
            for plot_name, _, _ in PLOTS:
                raw = accumulator_to_histogram(plot_name, data_accumulators[plot_name], y_title="Events")
                write_hist(output, "raw/data/observed", raw, plot_name)

                shape = raw.Clone(f"{plot_name}_data_shape")
                shape.SetDirectory(0)
                integral = shape.Integral()
                if integral > 0.0:
                    shape.Scale(1.0 / integral)
                shape.GetYaxis().SetTitle("Arbitrary units")
                write_hist(output, "shapes/data/observed", shape, plot_name)

        meta_dir = mkdirs(output, "metadata")
        meta_dir.cd()
        ROOT.TNamed("inputs", metadata).Write("inputs", ROOT.TObject.kOverwrite)
    finally:
        output.Close()


def read_accumulators_from_cache(cache_path: str, process: str) -> Optional[Dict[str, HistogramAccumulator]]:
    root_file = ROOT.TFile.Open(cache_path, "READ")
    if not root_file or root_file.IsZombie():
        print(f"[WARN] Could not open cache: {cache_path}")
        return None
    try:
        accumulators = empty_accumulators()
        any_found = False
        for plot_name, variable, _ in PLOTS:
            hist = root_file.Get(f"raw/groups/{process}/{plot_name}")
            if not hist:
                hist = root_file.Get(f"raw/{process}/{plot_name}")
            if not hist:
                continue
            any_found = True
            acc = HistogramAccumulator.for_variable(variable)
            for ibin in range(1, hist.GetNbinsX() + 1):
                content = hist.GetBinContent(ibin)
                error = hist.GetBinError(ibin)
                acc.sumw[ibin - 1] = content
                acc.sumw2[ibin - 1] = error * error
            accumulators[plot_name] = acc
        return accumulators if any_found else None
    finally:
        root_file.Close()


def hist_to_accumulator(hist, variable: str) -> HistogramAccumulator:
    acc = HistogramAccumulator.for_variable(variable)
    for ibin in range(1, hist.GetNbinsX() + 1):
        content = hist.GetBinContent(ibin)
        error = hist.GetBinError(ibin)
        acc.sumw[ibin - 1] = content
        acc.sumw2[ibin - 1] = error * error
    return acc


def read_data_accumulators_from_cache(cache_path: str) -> Optional[Dict[str, HistogramAccumulator]]:
    root_file = ROOT.TFile.Open(cache_path, "READ")
    if not root_file or root_file.IsZombie():
        print(f"[WARN] Could not open cache: {cache_path}")
        return None
    try:
        accumulators = empty_accumulators()
        any_found = False
        for plot_name, variable, _ in PLOTS:
            hist = root_file.Get(f"raw/data/observed/{plot_name}")
            if not hist:
                continue
            any_found = True
            acc = HistogramAccumulator.for_variable(variable)
            for ibin in range(1, hist.GetNbinsX() + 1):
                content = hist.GetBinContent(ibin)
                error = hist.GetBinError(ibin)
                acc.sumw[ibin - 1] = content
                acc.sumw2[ibin - 1] = error * error
            accumulators[plot_name] = acc
        return accumulators if any_found else None
    finally:
        root_file.Close()


def read_signal_accumulators_from_cache(cache_path: str) -> Dict[str, Dict[str, HistogramAccumulator]]:
    root_file = ROOT.TFile.Open(cache_path, "READ")
    if not root_file or root_file.IsZombie():
        print(f"[WARN] Could not open signal cache: {cache_path}")
        return {}
    try:
        directory = root_file.GetDirectory("raw/signals")
        if not directory:
            return {}
        output = {}
        for key in directory.GetListOfKeys():
            signal_name = key.GetName()
            accumulators = empty_accumulators()
            any_found = False
            for plot_name, variable, _ in PLOTS:
                hist = root_file.Get(f"raw/signals/{signal_name}/{plot_name}")
                if not hist:
                    continue
                accumulators[plot_name] = hist_to_accumulator(hist, variable)
                any_found = True
            if any_found:
                output[signal_name] = accumulators
        return output
    finally:
        root_file.Close()


def infer_process_from_cache(cache_path: str) -> Optional[str]:
    root_file = ROOT.TFile.Open(cache_path, "READ")
    if not root_file or root_file.IsZombie():
        return None
    try:
        directory = root_file.GetDirectory("raw/groups")
        if directory:
            keys = [key.GetName() for key in directory.GetListOfKeys()]
            if keys:
                return keys[0]
        directory = root_file.GetDirectory("raw")
        if directory:
            keys = [key.GetName() for key in directory.GetListOfKeys()]
            keys = [key for key in keys if key in PROC_LABEL]
            if keys:
                return keys[0]
        return None
    finally:
        root_file.Close()


def list_processes_from_cache(cache_path: str) -> List[str]:
    root_file = ROOT.TFile.Open(cache_path, "READ")
    if not root_file or root_file.IsZombie():
        return []
    try:
        directory = root_file.GetDirectory("raw/groups")
        if directory:
            return [
                key.GetName()
                for key in directory.GetListOfKeys()
                if key.GetName() in PROC_LABEL
            ]
        directory = root_file.GetDirectory("raw")
        if directory:
            return [
                key.GetName()
                for key in directory.GetListOfKeys()
                if key.GetName() in PROC_LABEL
            ]
        return []
    finally:
        root_file.Close()


def parse_merge_inputs(patterns: str) -> List[str]:
    import glob

    paths = []
    for token in [item.strip() for item in patterns.split(",") if item.strip()]:
        matches = sorted(glob.glob(token))
        paths.extend(matches if matches else [token])
    return list(dict.fromkeys(paths))


def parse_csv(value: str) -> List[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def ordered_signal_names(names: Sequence[str]) -> List[str]:
    default_index = {name: index for index, name in enumerate(DEFAULT_SIGNAL_DIRS)}
    return sorted(names, key=lambda name: (default_index.get(name, 10_000), name))


def signal_dir_name(signal_year: str, signal_name: str) -> str:
    return f"SVJ_UL{signal_year}_{signal_name}_13TeV-madgraphMLM-pythia8_n-1000"


def signal_name_from_training_dir(dirname: str, signal_year: str) -> Optional[str]:
    prefix = f"SVJ_UL{signal_year}_"
    suffix = "_13TeV-madgraphMLM-pythia8_n-1000"
    if dirname.startswith(prefix) and dirname.endswith(suffix):
        return dirname[len(prefix) : -len(suffix)]
    return None


def get_signal_directories(
    signal_base: str,
    signal_years: Sequence[str],
    requested: Sequence[str],
) -> List[Tuple[str, str, str]]:
    available: Dict[str, Dict[str, str]] = {}
    for signal_year in signal_years:
        private_samples = f"{normalize_eos_path(signal_base)}/{signal_year}/Full/PrivateSamples"
        entries = xrdfs_ls(private_samples)
        for entry in entries:
            dirname = os.path.basename(entry.rstrip("/"))
            signal_name = signal_name_from_training_dir(dirname, signal_year)
            if signal_name:
                available.setdefault(signal_name, {})[signal_year] = entry

    if len(requested) == 1 and requested[0].lower() in {"all", "*"}:
        selected_names = ordered_signal_names(available)
    else:
        selected_names = []
        seen = set()
        for token in requested:
            exact_matches = [name for name in available if token.lower() == name.lower()]
            matches = exact_matches or [
                name
                for name in available
                if token.lower() in name.lower()
            ]
            if not matches:
                print(f"[WARN] no training signal matches '{token}'")
                continue
            for name in ordered_signal_names(matches):
                if name in seen:
                    continue
                selected_names.append(name)
                seen.add(name)

    selected_names = [
        name
        for name in selected_names
        if not any(token.lower() in name.lower() for token in DISABLED_SIGNAL_TOKENS)
    ]

    selected = []
    for signal_name in selected_names:
        for signal_year in signal_years:
            signal_dir = available.get(signal_name, {}).get(signal_year)
            if not signal_dir:
                print(f"[WARN] Missing signal {signal_name} for {signal_year}")
                continue
            selected.append((signal_year, signal_name, signal_dir))
    return selected


def get_signal_files(sample_dir: str, regex: str, max_files: int = -1) -> List[str]:
    pattern = re.compile(regex)
    files = sorted(
        path
        for path in xrdfs_ls(sample_dir)
        if pattern.search(os.path.basename(path))
    )
    return files[:max_files] if max_files > 0 else files


def build_signal_accumulators(args) -> Dict[str, Dict[str, HistogramAccumulator]]:
    requested = parse_csv(args.signals)
    if args.no_signals or not requested:
        return {}
    signal_accumulators = {}
    signal_years = parse_csv(args.signal_years)
    for signal_year, signal_name, signal_dir in get_signal_directories(args.signal_base, signal_years, requested):
        files = get_signal_files(signal_dir, args.signal_file_regex, args.signal_max_files)
        if not files:
            print(f"[WARN] No signal ROOT files found for {signal_year}/{signal_name}")
            continue
        accumulators = signal_accumulators.setdefault(signal_name, empty_accumulators())
        entries = 0
        print(f"[SIGNAL] {signal_year}/{signal_name}: {len(files)} file(s)")
        for index, path in enumerate(files, start=1):
            entries += stream_file(eos_url(path), args, accumulators, is_data=False, apply_lund=args.apply_signal_lund)
            if index % 25 == 0 or index == len(files):
                print(f"[SIGNAL] {signal_year}/{signal_name}: processed {index}/{len(files)} files")
        print(f"[SIGNAL] {signal_year}/{signal_name}: entries={entries}")
    return signal_accumulators


def add_signal_accumulators(target, source) -> None:
    for signal_name, accumulators in source.items():
        if signal_name not in target:
            target[signal_name] = empty_accumulators()
        add_accumulators(target[signal_name], accumulators)


def write_merged_cache(output_path: str, per_process, metadata: str, data_accumulators=None, signal_accumulators=None) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    output = ROOT.TFile.Open(output_path, "RECREATE")
    if not output or output.IsZombie():
        raise RuntimeError(f"Could not create merged ROOT cache: {output_path}")
    try:
        for process, accumulators in sorted(per_process.items()):
            for plot_name, _, _ in PLOTS:
                raw = accumulator_to_histogram(plot_name, accumulators[plot_name], y_title="Events")
                write_hist(output, f"raw/groups/{process}", raw, plot_name)

        total_per_plot = {}
        for plot_name, variable, _ in PLOTS:
            total = HistogramAccumulator.for_variable(variable)
            for process in per_process:
                total.sumw += per_process[process][plot_name].sumw
                total.sumw2 += per_process[process][plot_name].sumw2
            total_per_plot[plot_name] = total

        for process, accumulators in sorted(per_process.items()):
            for plot_name, _, _ in PLOTS:
                shape = accumulator_to_histogram(plot_name, accumulators[plot_name], y_title="Arbitrary units")
                denom = total_per_plot[plot_name].integral()
                if denom > 0.0:
                    shape.Scale(1.0 / denom)
                write_hist(output, f"shapes/groups/{process}", shape, plot_name)

        if data_accumulators:
            for plot_name, _, _ in PLOTS:
                raw = accumulator_to_histogram(plot_name, data_accumulators[plot_name], y_title="Events")
                write_hist(output, "raw/data/observed", raw, plot_name)

                shape = raw.Clone(f"{plot_name}_data_shape")
                shape.SetDirectory(0)
                integral = shape.Integral()
                if integral > 0.0:
                    shape.Scale(1.0 / integral)
                shape.GetYaxis().SetTitle("Arbitrary units")
                write_hist(output, "shapes/data/observed", shape, plot_name)

        write_signal_hists(output, signal_accumulators)

        meta_dir = mkdirs(output, "metadata")
        meta_dir.cd()
        ROOT.TNamed("inputs", metadata).Write("inputs", ROOT.TObject.kOverwrite)
    finally:
        output.Close()


def merge_caches(args) -> None:
    paths = parse_merge_inputs(args.merge_inputs)
    if not paths:
        raise RuntimeError("No merge input caches found")

    per_process = {}
    data_accumulators = None
    signal_accumulators = {}
    lumi_scale = float(args.lumi_pb) / float(args.cache_lumi_pb)
    signal_lumi_scale = float(args.lumi_pb) / float(args.signal_cache_lumi_pb)
    for path in paths:
        processes = list_processes_from_cache(path)
        for process in processes:
            accumulators = read_accumulators_from_cache(path, process)
            if accumulators is not None:
                if abs(lumi_scale - 1.0) > 1.0e-12:
                    scale_accumulators(accumulators, lumi_scale)
                if process not in per_process:
                    per_process[process] = empty_accumulators()
                add_accumulators(per_process[process], accumulators)
        data_from_cache = read_data_accumulators_from_cache(path)
        if data_from_cache is not None:
            if data_accumulators is None:
                data_accumulators = empty_accumulators()
            add_accumulators(data_accumulators, data_from_cache)
        cached_signals = {} if args.ignore_input_signals else read_signal_accumulators_from_cache(path)
        if cached_signals:
            if abs(lumi_scale - 1.0) > 1.0e-12:
                for accumulators in cached_signals.values():
                    scale_accumulators(accumulators, lumi_scale)
            add_signal_accumulators(signal_accumulators, cached_signals)
        if not processes and data_from_cache is None and not cached_signals:
            print(f"[WARN] Could not infer process or data from {path}; skipping")

    signal_cache_paths = parse_merge_inputs(args.signal_cache_inputs)
    for path in signal_cache_paths:
        cached_signals = read_signal_accumulators_from_cache(path)
        if not cached_signals:
            print(f"[WARN] No signals found in cache {path}")
            continue
        if abs(signal_lumi_scale - 1.0) > 1.0e-12:
            for accumulators in cached_signals.values():
                scale_accumulators(accumulators, signal_lumi_scale)
        add_signal_accumulators(signal_accumulators, cached_signals)

    direct_signals = build_signal_accumulators(args)
    add_signal_accumulators(signal_accumulators, direct_signals)

    if not per_process:
        raise RuntimeError("No process histograms were found in merge inputs")

    metadata = (
        f"merged_inputs={len(paths)}; processes={','.join(sorted(per_process))}; "
        f"cache_lumi_pb={args.cache_lumi_pb}; output_lumi_pb={args.lumi_pb}; "
        f"lumi_scale={lumi_scale}; has_data={data_accumulators is not None}; "
        f"signal_cache_lumi_pb={args.signal_cache_lumi_pb}; "
        f"signal_lumi_scale={signal_lumi_scale}; "
        f"ignore_input_signals={args.ignore_input_signals}; "
        f"signals={','.join(ordered_signal_names(signal_accumulators))}"
    )
    write_merged_cache(
        args.output,
        per_process,
        metadata,
        data_accumulators=data_accumulators,
        signal_accumulators=signal_accumulators,
    )
    if not args.no_plot:
        plot_cache(args)


def divisor_generator(number: int) -> Iterable[int]:
    large = []
    for divisor in range(1, int(math.sqrt(number)) + 1):
        if number % divisor == 0:
            yield divisor
            if divisor * divisor != number:
                large.append(number // divisor)
    yield from reversed(large)


def rebin_factor(variable: str) -> int:
    nbins = int(VARIABLES[variable]["nbins"])
    divisors = list(divisor_generator(nbins))
    desired = nbins / float(PLOT_TARGET_BINS)
    return max(1, int(min(divisors, key=lambda value: abs(value - desired))))


def clone_hist(hist, name: str):
    out = hist.Clone(name)
    out.SetDirectory(0)
    return out


def get_hist(root_file, path: str, name: str):
    obj = root_file.Get(path)
    if not obj or not obj.InheritsFrom("TH1"):
        return None
    return clone_hist(obj, name)


def style_background(hist, process: str) -> None:
    color = PROC_COLOR.get(process, ROOT.kGray)
    hist.SetFillColor(color)
    hist.SetFillStyle(1001)
    hist.SetLineColor(color)
    hist.SetLineWidth(0)
    hist.SetMarkerSize(0)


def style_data(hist) -> None:
    hist.SetLineColor(ROOT.kBlack)
    hist.SetMarkerColor(ROOT.kBlack)
    hist.SetMarkerStyle(ROOT.kFullCircle)
    hist.SetMarkerSize(0.85)
    hist.SetLineWidth(2)
    hist.SetFillStyle(0)


def style_signal(hist, index: int) -> None:
    hist.SetLineColor(SIGNAL_LINE_COLORS[index % len(SIGNAL_LINE_COLORS)])
    hist.SetLineStyle(SIGNAL_LINE_STYLES[index % len(SIGNAL_LINE_STYLES)])
    hist.SetLineWidth(3)
    hist.SetFillStyle(0)
    hist.SetMarkerSize(0)


def signal_legend_label(signal_name: str) -> str:
    mmed = re.search(r"mMed-([0-9]+(?:p[0-9]+)?)", signal_name, re.IGNORECASE)
    rinv = re.search(r"rinv-([0-9]+(?:p[0-9]+)?)", signal_name, re.IGNORECASE)
    if mmed and rinv:
        return (
            f"m_{{#phi}} = {mmed.group(1).replace('p', '.')} GeV, "
            f"r_{{inv}} = {rinv.group(1).replace('p', '.')}"
        )
    return signal_name


def list_signal_names(root_file, cache_dir: str) -> List[str]:
    directory = root_file.GetDirectory(cache_dir)
    if not directory:
        return []
    return ordered_signal_names([key.GetName() for key in directory.GetListOfKeys()])


def signal_is_requested(signal_name: str, requested: Sequence[str]) -> bool:
    if any(token.lower() in signal_name.lower() for token in DISABLED_SIGNAL_TOKENS):
        return False
    if not requested:
        return False
    if any(token.lower() == "all" for token in requested):
        return True
    return any(token in signal_name for token in requested)


def suppress_qcd_met_display_spikes(hist) -> None:
    for ibin in range(2, hist.GetNbinsX()):
        xlow = hist.GetBinLowEdge(ibin)
        if xlow < QCD_MET_DISPLAY_SPIKE_MIN_X:
            continue

        content = hist.GetBinContent(ibin)
        left = hist.GetBinContent(ibin - 1)
        right = hist.GetBinContent(ibin + 1)
        if left <= 0.0 or right <= 0.0:
            continue
        if content < QCD_MET_DISPLAY_SPIKE_RATIO * left:
            continue
        if content < QCD_MET_DISPLAY_SPIKE_RATIO * right:
            continue

        replacement = 0.5 * (left + right)
        replacement_error = 0.5 * math.hypot(hist.GetBinError(ibin - 1), hist.GetBinError(ibin + 1))
        print(
            f"[INFO] Plot-side QCD MET spike cleanup: "
            f"{xlow:g}-{xlow + hist.GetBinWidth(ibin):g} GeV "
            f"{content:.6g} -> {replacement:.6g}"
        )
        hist.SetBinContent(ibin, replacement)
        hist.SetBinError(ibin, replacement_error)


def suppress_isolated_display_spikes(hist, variable: str, ratio_threshold: float = DPHI_DISPLAY_SPIKE_RATIO) -> None:
    updates = []
    contents = [hist.GetBinContent(ibin) for ibin in range(hist.GetNbinsX() + 2)]
    errors = [hist.GetBinError(ibin) for ibin in range(hist.GetNbinsX() + 2)]
    for ibin in range(2, hist.GetNbinsX()):
        content = hist.GetBinContent(ibin)
        left = contents[ibin - 1]
        right = contents[ibin + 1]
        if left <= 0.0 or right <= 0.0:
            continue

        window = [
            contents[index]
            for index in range(max(1, ibin - 2), min(hist.GetNbinsX(), ibin + 2) + 1)
            if index != ibin and contents[index] > 0.0
        ]
        if not window:
            continue
        local_level = float(np.median(window))
        if content < ratio_threshold * local_level:
            continue
        replacement = local_level
        replacement_error = float(np.median([errors[index] for index in range(max(1, ibin - 2), min(hist.GetNbinsX(), ibin + 2) + 1) if index != ibin]))
        updates.append((ibin, replacement, replacement_error))

    for ibin, replacement, replacement_error in updates:
        xlow = hist.GetBinLowEdge(ibin)
        print(
            f"[INFO] Plot-side {variable} spike cleanup: "
            f"{xlow:g}-{xlow + hist.GetBinWidth(ibin):g} "
            f"{hist.GetBinContent(ibin):.6g} -> {replacement:.6g}"
        )
        hist.SetBinContent(ibin, replacement)
        hist.SetBinError(ibin, replacement_error)


def pad_hist_xmax(hist, target_xmax: float, name: str):
    xaxis = hist.GetXaxis()
    old_xmax = xaxis.GetXmax()
    if old_xmax >= target_xmax - 1.0e-9:
        return hist

    old_nbins = hist.GetNbinsX()
    xmin = xaxis.GetXmin()
    bin_width = hist.GetBinWidth(old_nbins)
    extra_bins = int(math.ceil((target_xmax - old_xmax) / bin_width))
    new_xmax = old_xmax + extra_bins * bin_width
    padded = ROOT.TH1D(name, hist.GetTitle(), old_nbins + extra_bins, xmin, new_xmax)
    padded.Sumw2()
    padded.SetDirectory(0)
    for ibin in range(0, old_nbins + 2):
        padded.SetBinContent(ibin, hist.GetBinContent(ibin))
        padded.SetBinError(ibin, hist.GetBinError(ibin))
    return padded


def pad_hist_for_display(hist, variable: str, name: str):
    if variable in DISCRETE_AXIS_VARIABLES:
        xmin = VARIABLES[variable].get("display_xmin", 0.0)
        xmax = VARIABLES[variable].get("display_xmax", VARIABLES[variable]["xmax"])
        nbins = int(round(xmax - xmin))
        shifted = ROOT.TH1D(name, hist.GetTitle(), nbins, xmin, xmax)
        shifted.Sumw2()
        shifted.SetDirectory(0)
        for source_bin in range(1, hist.GetNbinsX() + 1):
            count_value = hist.GetBinCenter(source_bin) + 0.5
            if count_value < xmin or count_value >= xmax:
                continue
            target_bin = shifted.FindBin(count_value)
            shifted.SetBinContent(target_bin, shifted.GetBinContent(target_bin) + hist.GetBinContent(source_bin))
            shifted.SetBinError(
                target_bin,
                math.hypot(shifted.GetBinError(target_bin), hist.GetBinError(source_bin)),
            )
        return shifted
    target_xmax = VARIABLES[variable].get("display_xmax", VARIABLES[variable]["xmax"])
    return pad_hist_xmax(hist, target_xmax, name)


def draw_cms_label(
    lumi_fb: float,
    normalized: bool,
    has_data: bool = False,
    preliminary: bool = False,
) -> None:
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextColor(ROOT.kBlack)
    latex.SetTextFont(61)
    latex.SetTextSize(0.060)
    latex.SetTextAlign(11)
    latex.DrawLatex(0.16, 0.94, "CMS")
    if not has_data:
        latex.SetTextFont(52)
        latex.SetTextSize(0.040)
        latex.DrawLatex(0.290, 0.94, "Simulation")
        if preliminary:
            latex.DrawLatex(0.465, 0.94, "Preliminary")
    elif preliminary:
        latex.SetTextFont(52)
        latex.SetTextSize(0.040)
        latex.DrawLatex(0.255, 0.94, "Preliminary")
    latex.SetTextFont(42)
    latex.SetTextSize(0.045)
    latex.SetTextAlign(31)
    latex.DrawLatex(0.95, 0.94, "(13 TeV)")


def draw_cut_marker(axis, variable: str):
    if "line" not in VARIABLES[variable]:
        return None, None

    x = VARIABLES[variable]["line"]
    y1 = axis.GetMinimum()
    y2 = axis.GetMaximum()
    line = ROOT.TLine(x, y1, x, y2)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineStyle(ROOT.kDashed)
    line.SetLineWidth(3)
    line.Draw("same")

    xmin = VARIABLES[variable]["xmin"]
    xmax = VARIABLES[variable].get("display_xmax", VARIABLES[variable]["xmax"])
    width = max(1.0e-9, xmax - xmin)
    side = VARIABLES[variable].get("line_side", "right")
    if side == "left":
        arrow_symbol = "#leftarrow"
        arrow_align = 32
    else:
        arrow_symbol = "#rightarrow"
        arrow_align = 12

    arrow_y_fraction = VARIABLES[variable].get("arrow_y_fraction", 0.80)
    y_arrow = y1 + arrow_y_fraction * (y2 - y1)
    if ROOT.gPad.GetLogy() and y1 > 0.0 and y2 > y1:
        y_arrow = math.exp(math.log(y1) + arrow_y_fraction * (math.log(y2) - math.log(y1)))
    arrow = ROOT.TLatex()
    arrow.SetTextColor(ROOT.kBlack)
    arrow.SetTextFont(42)
    arrow.SetTextSize(0.055)
    arrow.SetTextAlign(arrow_align)
    arrow.DrawLatex(x, y_arrow, arrow_symbol)
    return line, arrow


def unique_legend_entries(entries):
    seen = set()
    unique = []
    for hist, label, option in entries:
        if label in seen:
            continue
        seen.add(label)
        unique.append((hist, label, option))
    return unique


def style_legend(legend, text_size: float = 0.031, margin: float = 0.25) -> None:
    legend.SetBorderSize(0)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetFillStyle(1001)
    legend.SetTextSize(text_size)
    legend.SetTextFont(42)
    legend.SetEntrySeparation(0.030)
    legend.SetMargin(margin)


def format_discrete_xaxis(hist, variable: str) -> None:
    if variable not in DISCRETE_AXIS_VARIABLES:
        return
    xaxis = hist.GetXaxis()
    xaxis.SetNdivisions(510 if variable == "njetsAK8" else 505, False)
    xaxis.SetLabelSize(0.040)


def draw_plot(root_file, plot_name: str, variable: str, plot_dir: str, args) -> None:
    grouped_cache_dir = "shapes/groups" if args.normalized else "raw/groups"
    legacy_cache_dir = "shapes/QCD" if args.normalized else "raw/QCD"
    data_cache_dir = "shapes/data/observed" if args.normalized else "raw/data/observed"
    signal_cache_dir = "shapes/signals" if args.normalized else "raw/signals"
    factor = rebin_factor(variable)

    process_hists = {}
    for process in STACK_ORDER:
        hist = get_hist(root_file, f"{grouped_cache_dir}/{process}/{plot_name}", f"display_{process}_{plot_name}")
        if not hist:
            continue
        if factor > 1:
            hist = hist.Rebin(factor, f"display_{process}_{plot_name}_rebin{factor}")
            hist.SetDirectory(0)
        hist = pad_hist_for_display(hist, variable, f"display_{process}_{plot_name}_padded")
        if process == "QCD" and plot_name == "h_MET_pre__metcut":
            suppress_qcd_met_display_spikes(hist)
        if variable == "dPhiMinjMETAK8":
            suppress_isolated_display_spikes(hist, variable)
        if hist.Integral() <= 0.0:
            continue
        style_background(hist, process)
        process_hists[process] = hist

    if not process_hists:
        hist = get_hist(root_file, f"{legacy_cache_dir}/{plot_name}", f"display_{plot_name}")
        if not hist:
            print(f"[WARN] Missing cached histogram: {grouped_cache_dir}/*/{plot_name}")
            return
        if factor > 1:
            hist = hist.Rebin(factor, f"display_{plot_name}_rebin{factor}")
            hist.SetDirectory(0)
        hist = pad_hist_for_display(hist, variable, f"display_{plot_name}_padded")
        if variable == "dPhiMinjMETAK8":
            suppress_isolated_display_spikes(hist, variable)
        style_background(hist, "QCD")
        process_hists["QCD"] = hist

    data_hist = get_hist(root_file, f"{data_cache_dir}/{plot_name}", f"display_data_{plot_name}")
    if data_hist:
        if factor > 1:
            data_hist = data_hist.Rebin(factor, f"display_data_{plot_name}_rebin{factor}")
            data_hist.SetDirectory(0)
        data_hist = pad_hist_for_display(data_hist, variable, f"display_data_{plot_name}_padded")
        if data_hist.Integral() <= 0.0:
            data_hist = None
        else:
            style_data(data_hist)

    requested_signals = parse_csv(args.signals)
    signal_hists = []
    for signal_name in list_signal_names(root_file, signal_cache_dir):
        if not signal_is_requested(signal_name, requested_signals):
            continue
        index = len(signal_hists)
        hist = get_hist(
            root_file,
            f"{signal_cache_dir}/{signal_name}/{plot_name}",
            f"display_signal_{index}_{plot_name}",
        )
        if not hist:
            continue
        if factor > 1:
            hist = hist.Rebin(factor, f"display_signal_{index}_{plot_name}_rebin{factor}")
            hist.SetDirectory(0)
        hist = pad_hist_for_display(hist, variable, f"display_signal_{index}_{plot_name}_padded")
        if not args.normalized and args.signal_scale != 1.0:
            hist.Scale(float(args.signal_scale))
        if variable == "dPhiMinjMETAK8":
            suppress_isolated_display_spikes(hist, variable)
        if hist.Integral() <= 0.0:
            continue
        style_signal(hist, index)
        signal_hists.append((signal_name, hist))

    total = None
    stack = ROOT.THStack(f"stack_{plot_name}", "")
    for process in STACK_ORDER:
        hist = process_hists.get(process)
        if not hist:
            continue
        stack.Add(hist)
        if total is None:
            total = clone_hist(hist, f"total_{plot_name}")
        else:
            total.Add(hist)

    canvas = ROOT.TCanvas(f"c_{plot_name}", "", 800, 800)
    canvas.cd()
    ROOT.gPad.SetLeftMargin(0.16)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.08)
    ROOT.gPad.SetBottomMargin(0.12)
    ROOT.gPad.SetTicks(1, 1)
    ROOT.gPad.SetLogy(True)

    if variable in DISCRETE_AXIS_VARIABLES:
        axis = ROOT.TH1D(
            f"axis_{plot_name}",
            "",
            100,
            VARIABLES[variable]["display_xmin"],
            VARIABLES[variable]["display_xmax"],
        )
        axis.SetDirectory(0)
    else:
        axis = clone_hist(total, f"axis_{plot_name}")
        axis.Reset("ICESM")
    axis.SetStats(0)
    axis.SetTitle("")
    axis.GetXaxis().SetTitle(VARIABLES[variable]["title"])
    axis.GetYaxis().SetTitle("Arbitrary units" if args.normalized else "Events")
    axis.GetXaxis().SetTitleSize(0.05)
    axis.GetXaxis().SetLabelSize(0.04)
    axis.GetYaxis().SetTitleSize(0.05)
    axis.GetYaxis().SetLabelSize(0.04)
    axis.GetYaxis().SetLabelOffset(0.01)
    axis.GetYaxis().SetTitleOffset(1.30)
    if variable == "ST":
        axis.GetXaxis().SetNdivisions(505)
    format_discrete_xaxis(axis, variable)
    if variable != "njetsAK8":
        axis.GetXaxis().SetRangeUser(
            VARIABLES[variable].get("display_xmin", VARIABLES[variable]["xmin"]),
            VARIABLES[variable].get("display_xmax", VARIABLES[variable]["xmax"]),
        )
    axis.SetMinimum(1.0e-5 if args.normalized else 1.0e-1)
    ymax = max(total.GetMaximum(), 1.0)
    if data_hist:
        ymax = max(ymax, data_hist.GetMaximum())
    for _, hist in signal_hists:
        ymax = max(ymax, hist.GetMaximum())
    axis.SetMaximum(ymax * 60.0)
    axis.Draw("hist")
    stack.Draw("hist same")
    for _, hist in signal_hists:
        hist.Draw("hist same")
    if data_hist:
        data_hist.Draw("E1 X0 same")

    cut_line, cut_arrow = draw_cut_marker(axis, variable)

    legend_y1 = 0.695 if signal_hists else 0.690
    if signal_hists and variable == "MET":
        legend_x1 = 0.40
    elif signal_hists and variable == "ST":
        legend_x1 = 0.42
    elif signal_hists and variable == "dPhiMinjMETAK8":
        legend_x1 = 0.45
    elif signal_hists and variable == "njetsAK8":
        legend_x1 = 0.42
    elif signal_hists and variable in ("nl", "nelectrons", "nmuons"):
        legend_x1 = 0.45
    else:
        legend_x1 = 0.48 if signal_hists else 0.72
    background_entries = [
        (process_hists[process], PROC_LABEL.get(process, process), "F")
        for process in LEGEND_ORDER
        if process in process_hists
    ]
    if data_hist:
        background_entries.insert(0, (data_hist, "Data", "EP"))
    signal_entries = unique_legend_entries(
        [
            (
                hist,
                signal_legend_label(signal_name)
                + (f" (#times {args.signal_scale:g})" if not args.normalized and args.signal_scale != 1.0 else ""),
                "L",
            )
            for signal_name, hist in signal_hists
        ]
    )
    if signal_hists:
        signal_entries.append(("", "m_{dark} = 20 GeV, #kern[-0.18]{#lambda} = 1", ""))
    if signal_hists:
        if variable == "dPhiMinjMETAK8":
            legend = ROOT.TLegend(0.58, 0.525, 0.98, 0.885)
            style_legend(legend, text_size=0.026, margin=0.14)
            for entry in background_entries:
                legend.AddEntry(*entry)
            for entry in signal_entries:
                legend.AddEntry(*entry)
            legend.Draw()
        else:
            background_legend = ROOT.TLegend(legend_x1, legend_y1 - 0.020, legend_x1 + 0.20, 0.885)
            style_legend(background_legend)
            for entry in background_entries:
                background_legend.AddEntry(*entry)
            background_legend.Draw()

            if variable == "MET":
                signal_x1 = 0.59
            elif variable == "ST":
                signal_x1 = 0.59
            elif variable == "njetsAK8":
                signal_x1 = 0.63
            else:
                signal_x1 = max(0.63, legend_x1 + 0.21)
            signal_legend = ROOT.TLegend(signal_x1, 0.690, 0.995, 0.885)
            style_legend(signal_legend, text_size=0.027, margin=0.14)
            for entry in signal_entries:
                signal_legend.AddEntry(*entry)
            signal_legend.Draw()
    else:
        legend = ROOT.TLegend(legend_x1, legend_y1, 0.98, 0.895)
        style_legend(legend, text_size=0.034)
        for entry in background_entries:
            legend.AddEntry(*entry)
        legend.Draw()

    draw_cms_label(
        args.lumi_pb / 1000.0,
        args.normalized,
        has_data=bool(data_hist),
        preliminary=args.preliminary_label,
    )
    ROOT.gPad.RedrawAxis()

    os.makedirs(plot_dir, exist_ok=True)
    tag = "_normalized" if args.normalized else ""
    outbase = os.path.join(plot_dir, f"{plot_name}{tag}")
    canvas.SaveAs(outbase + ".pdf")
    canvas.SaveAs(outbase + ".png")
    canvas.Close()
    print(f"[OK] Wrote {outbase}.pdf/.png")


def plot_cache(args) -> None:
    root_file = ROOT.TFile.Open(args.output, "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Could not open ROOT cache: {args.output}")
    try:
        for plot_name, variable, _ in PLOTS:
            draw_plot(root_file, plot_name, variable, args.plot_dir, args)
    finally:
        root_file.Close()


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input-dir", default=DEFAULT_INPUT_DIR)
    parser.add_argument("--file-regex", default=r".*_RA2AnalysisTree\.root$")
    parser.add_argument("--process", default="QCD", help="process label used when writing per-sample caches")
    parser.add_argument("--sample-name", default="", help="optional sample label stored in cache metadata")
    parser.add_argument("--as-signal", action="store_true", help="write this input cache under raw/signals/<sample-name>")
    parser.add_argument("--data-dir", default=DEFAULT_DATA_DIR, help="optional EOS data directory to overlay; may contain ROOT files or one level of sample subdirs")
    parser.add_argument("--data-file-regex", default=r".*_RA2AnalysisTree\.root$")
    parser.add_argument("--data-max-files", type=int, default=-1)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    parser.add_argument("--plot-dir", default=DEFAULT_PLOT_DIR)
    parser.add_argument("--chunk-size", default="100 MB")
    parser.add_argument("--max-files", type=int, default=-1)
    parser.add_argument("--lumi-pb", type=float, default=LUMI_PB)
    parser.add_argument(
        "--cache-lumi-pb",
        type=float,
        default=LUMI_PB,
        help="lumi used when existing per-sample caches were produced; merge mode rescales to --lumi-pb",
    )
    parser.add_argument(
        "--signal-cache-lumi-pb",
        type=float,
        default=LUMI_PB,
        help="lumi used when caches passed through --signal-cache-inputs were produced",
    )
    parser.add_argument("--kfactor", type=float, default=1.0)
    parser.add_argument("--apply-pu", action="store_true", default=True)
    parser.add_argument("--no-pu", action="store_false", dest="apply_pu")
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--no-plot", action="store_true", help="build/merge cache but do not draw plots")
    parser.add_argument("--merge-inputs", default="", help="comma-separated ROOT cache paths or glob patterns to merge")
    parser.add_argument(
        "--normalized",
        action="store_true",
        default=True,
        help="draw cached shapes normalized like Figure2: stacked backgrounds integrate to one, signals unit area",
    )
    parser.add_argument(
        "--raw-yields",
        action="store_false",
        dest="normalized",
        help="draw raw event-yield histograms instead of normalized shapes",
    )
    parser.add_argument("--signal-base", default=DEFAULT_SIGNAL_BASE, help="base containing YEAR/Full/PrivateSamples signal dirs")
    parser.add_argument("--signal-years", default=DEFAULT_SIGNAL_YEARS, help="comma-separated signal years to read from --signal-base")
    parser.add_argument("--signals", default=",".join(DEFAULT_SIGNAL_DIRS), help="comma-separated signal names/tokens, or all")
    parser.add_argument("--signal-cache-inputs", default="", help="comma-separated ROOT caches containing raw/signals/* histograms")
    parser.add_argument("--ignore-input-signals", action="store_true", help="ignore raw/signals found in --merge-inputs")
    parser.add_argument("--signal-max-files", type=int, default=-1, help="limit files per signal while testing")
    parser.add_argument("--signal-file-regex", default=DEFAULT_SIGNAL_FILE_REGEX, help="regex selecting signal ROOT files inside each signal directory")
    parser.add_argument("--apply-signal-lund", action="store_true", help="apply lundWeightNom normalized by Initial/InitialLundNominal to directly read signals when those branches exist")
    parser.add_argument("--no-signal-lund", action="store_false", dest="apply_signal_lund")
    parser.add_argument("--signal-scale", type=float, default=1.0, help="scale raw signal overlays at draw time")
    parser.add_argument("--no-signals", action="store_true", help="do not read or overlay signals")
    parser.add_argument("--preliminary-label", action="store_true", help="draw CMS Preliminary or CMS Simulation Preliminary")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.merge_inputs:
        merge_caches(args)
        return 0
    if args.plot_only:
        plot_cache(args)
        return 0

    files = get_input_files(args.input_dir, args.file_regex, args.max_files)
    if not files:
        raise RuntimeError(f"No input files found in {args.input_dir}")

    accumulators = empty_accumulators()
    failures = []
    entries = 0

    print(f"[INFO] Input files: {len(files)}")
    for index, remote_file in enumerate(files, start=1):
        try:
            entries += stream_file(
                eos_url(remote_file),
                args,
                accumulators,
                apply_lund=args.as_signal and args.apply_signal_lund,
            )
        except Exception as exc:
            message = f"{os.path.basename(remote_file)}: {exc}"
            failures.append(message)
            print(f"[WARN] {message}")
        if index % 25 == 0 or index == len(files):
            print(f"[INFO] processed {index:5d}/{len(files):5d} files")

    data_accumulators = None
    data_entries = 0
    data_files = []
    if args.data_dir:
        data_files = get_input_files_with_subdirs(args.data_dir, args.data_file_regex, args.data_max_files)
        if not data_files:
            print(f"[WARN] No data files found in {args.data_dir}")
        else:
            data_accumulators = empty_accumulators()
            print(f"[INFO] Data files: {len(data_files)}")
            for index, remote_file in enumerate(data_files, start=1):
                try:
                    data_entries += stream_file(eos_url(remote_file), args, data_accumulators, is_data=True)
                except Exception as exc:
                    message = f"data {os.path.basename(remote_file)}: {exc}"
                    failures.append(message)
                    print(f"[WARN] {message}")
                if index % 25 == 0 or index == len(data_files):
                    print(f"[INFO] processed data {index:5d}/{len(data_files):5d} files")

    metadata = (
        f"input_dir={args.input_dir}; file_regex={args.file_regex}; "
        f"process={args.process}; sample_name={args.sample_name}; "
        f"files={len(files)}; entries={entries}; lumi_pb={args.lumi_pb}; "
        f"data_dir={args.data_dir}; data_files={len(data_files)}; data_entries={data_entries}; "
        f"kfactor={args.kfactor}; apply_pu={args.apply_pu}; "
        "cuts: MET>200, ST>1300, nGoodAK8>=2, nGoodLeptons==0, DeltaPhiMin_AK8<=1.5"
    )
    signal_name = args.sample_name if args.as_signal else None
    if args.as_signal and not signal_name:
        raise RuntimeError("--as-signal requires --sample-name")
    write_cache(
        args.output,
        accumulators,
        metadata,
        process=args.process,
        data_accumulators=data_accumulators,
        signal_name=signal_name,
    )
    if not args.no_plot:
        plot_cache(args)

    if failures:
        print(f"[WARN] Completed with {len(failures)} failed file(s). First few:")
        for failure in failures[:10]:
            print(f"  - {failure}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[FATAL] {exc}", file=sys.stderr)
        raise
