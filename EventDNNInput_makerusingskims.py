#!/usr/bin/env python3
"""
Build normalized Event-DNN input plots from skim Events trees.

This follows the Figure2_makerusingskims.py pattern:
  1. stream EOS skim trees into a compact ROOT cache;
  2. write raw and normalized TH1D objects;
  3. draw legacy-style MC-only stack + signal overlays from the cache.

Example:
  source condor/initCondor.sh
  python3 -u EventDNNInput_makerusingskims.py \
    --background-base /store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/data_mc_noWNAE \
    --data-base /store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/data_mc_noWNAE \
    --signal-base /store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto \
    --years 2016,2017,2018 \
    --output EventDNNInput_plots/event_dnn_inputs.root \
    --plot-dir EventDNNInput_plots \
    --chunk-size "100 MB"



python3 -u EventDNNInput_makerusingskims.py \
  --background-base /store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/data_mc_WNAE \
  --data-base /store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/data_mc_WNAE \
  --signal-base /store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto \
  --years 2016,2017,2018 \
  --output supplementry_plots_postTrimandgapveto/DNN/EventDNNInput_plots_WNAE/event_dnn_inputs.root \
  --plot-dir supplementry_plots_postTrimandgapveto/DNN/EventDNNInput_plots_WNAE
  
to redraw

python3 -u EventDNNInput_makerusingskims.py \
  --plot-only \
  --output supplementry_plots_postTrimandgapveto/DNN/EventDNNInput_plots/event_dnn_inputs.root \
  --plot-dir supplementry_plots_postTrimandgapveto/DNN/EventDNNInput_plots

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

import Figure2_makerusingskims as f2

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetLineStyleString(4, "28 8 4 8")
ROOT.TH1.SetDefaultSumw2(True)
ROOT.TH1.AddDirectory(False)


EOS_HOST = "root://cmseos.fnal.gov"
DEFAULT_BACKGROUND_BASE = "/store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/data_mc_WNAE"
DEFAULT_DATA_BASE = DEFAULT_BACKGROUND_BASE
DEFAULT_SIGNAL_BASE = "/store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto"

LUMI_PB = {
    "2016": 36.31 * 1000.0,
    "2017": 42.07 * 1000.0,
    "2018": 59.56 * 1000.0,
}

VARIABLES = {
    "MET": {
        "source": "MET",
        "title": "p_{T}^{miss} [GeV]",
        "nbins": 500,
        "xmin": 200.0,
        "xmax": 2000.0,
    },
    "HT": {
        "source": "HT",
        "title": "H_{T} [GeV]",
        "nbins": 500,
        "xmin": 0.0,
        "xmax": 5000.0,
    },
    "ST": {
        "source": "ST",
        "title": "S_{T} [GeV]",
        "nbins": 500,
        "xmin": 1300.0,
        "xmax": 5000.0,
    },
    "mT": {
        "source": "MT_AK8",
        "title": "m_{T} [GeV]",
        "nbins": 500,
        "xmin": 0.0,
        "xmax": 6000.0,
    },
    "dPhiMinjMETAK8": {
        "source": "DeltaPhiMinGoodJetsAK8",
        "title": "#Delta#phi_{min}(J, p_{T}^{miss})",
        "nbins": 100,
        "xmin": 0.0,
        "xmax": 2.0,
    },
    "nsvjJetsAK8": {
        "source": "nsvjJetsAK8",
        "title": "Number of SVJ AK8 jets",
        "nbins": 20,
        "xmin": 0.0,
        "xmax": 5.0,
    },
    "dnnEventClassScore": {
        "source": "dnnEventClassScore",
        "title": "Event Classifier Score",
        "nbins": 100,
        "xmin": 0.0,
        "xmax": 1.0,
        "wp": 0.85,
        "wp_side": "right",
    },
    "LundJetPlaneZ01GoodJetsAK8": {
        "source": "LundJetPlaneZ01GoodJetsAK8",
        "title": "z(J_{1}, J_{2})",
        "nbins": 40,
        "xmin": 0.0,
        "xmax": 0.5,
        "log_ymin": 10.0,
    },
    "GoodJetsAK80_LundJetPlaneZ": {
        "source": "JetsAK8_LundJetPlaneZ",
        "jet_index": 0,
        "title": "z(J_{1}, p_{T}^{miss})",
        "nbins": 40,
        "xmin": 0.0,
        "xmax": 0.5,
    },
    "GoodJetsAK81_LundJetPlaneZ": {
        "source": "JetsAK8_LundJetPlaneZ",
        "jet_index": 1,
        "title": "z(J_{2}, p_{T}^{miss})",
        "nbins": 40,
        "xmin": 0.0,
        "xmax": 0.5,
    },
    "DeltaPhi01GoodJetsAK8": {
        "source": "DeltaPhi01GoodJetsAK8",
        "title": "#Delta#phi(J_{1}, J_{2})",
        "nbins": 40,
        "xmin": -3.2,
        "xmax": 3.2,
    },
    "GoodJetsAK80_deltaPhiMET": {
        "source": "JetsAK8_deltaPhiMET",
        "jet_index": 0,
        "title": "#Delta#phi(J_{1}, p_{T}^{miss})",
        "nbins": 40,
        "xmin": -3.2,
        "xmax": 3.2,
    },
    "GoodJetsAK81_deltaPhiMET": {
        "source": "JetsAK8_deltaPhiMET",
        "jet_index": 1,
        "title": "#Delta#phi(J_{2}, p_{T}^{miss})",
        "nbins": 64,
        "xmin": -3.2,
        "xmax": 3.2,
    },
    "DeltaEta01GoodJetsAK8": {
        "source": "DeltaEta01GoodJetsAK8",
        "title": "#Delta#eta(J_{1}, J_{2})",
        "nbins": 40,
        "xmin": -5.0,
        "xmax": 5.0,
        "log_ymin": 10.0,
    },
    "DeltaR01GoodJetsAK8": {
        "source": "DeltaR01GoodJetsAK8",
        "title": "#DeltaR(J_{1}, J_{2})",
        "nbins": 40,
        "xmin": 0.5,
        "xmax": 5.5,
        "log_ymin": 10.0,
    },
    "DijetMass01GoodJetsAK8Log": {
        "log_source": "DijetMass01GoodJetsAK8",
        "title": "log m(J_{1}, J_{2})",
        "nbins": 60,
        "xmin": 5.0,
        "xmax": 10.0,
    },
    "GoodJetsAK80_MTMETLog": {
        "source": "JetsAK8_MTMET",
        "jet_index": 0,
        "title": "m_{T}(J_{1}, p_{T}^{miss}) [GeV]",
        "nbins": 40,
        "xmin": 20.0,
        "xmax": 8000.0,
        "log_bins": True,
        "log_x": True,
        "display_rebin_factor_override": 1,
    },
    "GoodJetsAK81_MTMETLog": {
        "source": "JetsAK8_MTMET",
        "jet_index": 1,
        "title": "m_{T}(J_{2}, p_{T}^{miss}) [GeV]",
        "nbins": 40,
        "xmin": 20.0,
        "xmax": 8000.0,
        "log_bins": True,
        "log_x": True,
        "display_rebin_factor_override": 1,
    },
}

JET_VARIABLES = {
    "PtAK8": ("JetsAK8_/.fPt", "p_{T}(J_{idx}) [GeV]", 400, 0.0, 4000.0),
    "EtaAK8": ("JetsAK8_/.fEta", "#eta(J_{idx})", 80, -4.0, 4.0),
    "PhiAK8": ("JetsAK8_/.fPhi", "#phi(J_{idx})", 64, -3.2, 3.2),
    "MassAK8": ("JetsAK8_mass", "m(J_{idx}) [GeV]", 200, 0.0, 900.0),
    "MTMETAK8": ("JetsAK8_MTMET", "M_{T}(J_{idx}, p_{T}^{miss}) [GeV]", 200, 0.0, 6000.0),
    "DeltaPhiMETAK8": ("JetsAK8_deltaPhiMET", "#Delta#phi(J_{idx}, p_{T}^{miss})", 64, -3.2, 3.2),
    "LundJetPlaneZAK8": ("JetsAK8_LundJetPlaneZ", "z(J_{idx}, p_{T}^{miss})", 60, 0.0, 0.5),
    "Tau1AK8": ("JetsAK8_NsubjettinessTau1", "#tau_{1}(J_{idx})", 60, 0.0, 1.5),
    "Tau2AK8": ("JetsAK8_NsubjettinessTau2", "#tau_{2}(J_{idx})", 60, 0.0, 1.5),
    "Tau3AK8": ("JetsAK8_NsubjettinessTau3", "#tau_{3}(J_{idx})", 60, 0.0, 1.5),
    "Tau4AK8": ("JetsAK8_NsubjettinessTau4", "#tau_{4}(J_{idx})", 60, 0.0, 1.5),
    "Tau5AK8": ("JetsAK8_NsubjettinessTau5", "#tau_{5}(J_{idx})", 60, 0.0, 1.5),
    "NConstituentsSoftDropAK8": ("JetsAK8_nConstituentsSoftDrop", "N_{const}^{SD}(J_{idx})", 80, 0.0, 200.0),
    "PNetScoreAK8": ("JetsAK8_pNetJetTaggerScore", "ParticleNet score(J_{idx})", 100, 0.0, 1.0),
    "GirthAK8": ("JetsAK8_girth", "girth(J_{idx})", 60, 0.0, 1.0),
    "PtDAK8": ("JetsAK8_ptD", "p_{T}D(J_{idx})", 60, 0.0, 1.0),
    "AxisMajorAK8": ("JetsAK8_axismajor", "axis major(J_{idx})", 60, 0.0, 1.0),
    "AxisMinorAK8": ("JetsAK8_axisminor", "axis minor(J_{idx})", 60, 0.0, 1.0),
    "EcfC2b1AK8": ("JetsAK8_ecfC2b1", "ECF C_{2}^{#beta=1}(J_{idx})", 60, 0.0, 1.0),
    "EcfC2b2AK8": ("JetsAK8_ecfC2b2", "ECF C_{2}^{#beta=2}(J_{idx})", 60, 0.0, 1.0),
    "EcfD2b1AK8": ("JetsAK8_ecfD2b1", "ECF D_{2}^{#beta=1}(J_{idx})", 60, 0.0, 10.0),
    "EcfD2b2AK8": ("JetsAK8_ecfD2b2", "ECF D_{2}^{#beta=2}(J_{idx})", 60, 0.0, 10.0),
    "EcfN2b1AK8": ("JetsAK8_ecfN2b1", "ECF N_{2}^{#beta=1}(J_{idx})", 60, 0.0, 1.0),
    "EcfN2b2AK8": ("JetsAK8_ecfN2b2", "ECF N_{2}^{#beta=2}(J_{idx})", 60, 0.0, 1.0),
}

for jet_index in (1, 2, 3):
    for suffix, (source, title, nbins, xmin, xmax) in JET_VARIABLES.items():
        VARIABLES[f"j{jet_index}{suffix}"] = {
            "source": source,
            "jet_index": jet_index - 1,
            "title": title.replace("{idx}", str(jet_index)),
            "nbins": nbins,
            "xmin": xmin,
            "xmax": xmax,
        }

for jet_index in (1, 2, 3):
    for suffix, numerator, denominator, title in (
        ("Tau21AK8", "JetsAK8_NsubjettinessTau2", "JetsAK8_NsubjettinessTau1", "#tau_{21}(J_{idx})"),
        ("Tau32AK8", "JetsAK8_NsubjettinessTau3", "JetsAK8_NsubjettinessTau2", "#tau_{32}(J_{idx})"),
        ("Tau43AK8", "JetsAK8_NsubjettinessTau4", "JetsAK8_NsubjettinessTau3", "#tau_{43}(J_{idx})"),
    ):
        VARIABLES[f"j{jet_index}{suffix}"] = {
            "numerator": numerator,
            "denominator": denominator,
            "jet_index": jet_index - 1,
            "title": title.replace("{idx}", str(jet_index)),
            "nbins": 40,
            "xmin": 0.0,
            "xmax": 1.0,
        }

PLOT_VARIABLES = [
    "dnnEventClassScore",
    "DeltaEta01GoodJetsAK8",
    "DeltaR01GoodJetsAK8",
    "DijetMass01GoodJetsAK8Log",
    "GoodJetsAK80_deltaPhiMET",
    "GoodJetsAK80_LundJetPlaneZ",
    "GoodJetsAK80_MTMETLog",
    "GoodJetsAK81_deltaPhiMET",
    "GoodJetsAK81_LundJetPlaneZ",
    "GoodJetsAK81_MTMETLog",
    "j3DeltaPhiMETAK8",
    "j3LundJetPlaneZAK8",
    "j3MTMETAK8",
    "LundJetPlaneZ01GoodJetsAK8",
]

VARIABLES = {variable: VARIABLES[variable] for variable in PLOT_VARIABLES}


def figure2_variable_config() -> Dict[str, Dict[str, object]]:
    converted: Dict[str, Dict[str, object]] = {}
    for name, cfg in VARIABLES.items():
        out = dict(cfg)
        out["x_title"] = str(out.get("x_title", out.get("title", name)))
        out.setdefault("title", out["x_title"])
        out["branch"] = str(out.get("branch", out.get("source", name)))
        out["display_rebin_factor"] = DISPLAY_BIN_REDUCTION
        out["ratio_ymin"] = 0.0
        out["ratio_ymax"] = 2.0
        out["signal_scale_factors"] = SIGNAL_SCALE_FACTORS
        if name == "dnnEventClassScore":
            out["x_title"] = "Event classifier score"
            out["display_bin_width"] = 0.05
        converted[name] = out
    return converted


def configure_figure2_plot_style() -> None:
    f2.VARIABLES = figure2_variable_config()
    f2.SIGNAL_SCALE_FACTORS = SIGNAL_SCALE_FACTORS


PLOT_TARGET_BINS = 40
DISPLAY_BIN_REDUCTION = 2

DEFAULT_SIGNAL_DIRS = [
    "t-channel_mMed-2000_mDark-20_rinv-0p1_alpha-peak_yukawa-1",
    "t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
    "t-channel_mMed-2000_mDark-20_rinv-0p7_alpha-peak_yukawa-1",
    "t-channel_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
    "t-channel_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
]

SIGNAL_SCALE_FACTORS = {
    ("600", "0.3"): 10.0,
    ("2000", "0.1"): 500.0,
    ("2000", "0.3"): 500.0,
    ("2000", "0.7"): 1000.0,
    ("4000", "0.3"): 1000.0,
}

SAMPLE_MANIFEST = {
    "QCD": [
        "QCD_Pt_1000to1400",
        "QCD_Pt_1400to1800",
        "QCD_Pt_170to300",
        "QCD_Pt_1800to2400",
        "QCD_Pt_2400to3200",
        "QCD_Pt_300to470",
        "QCD_Pt_3200toInf",
        "QCD_Pt_470to600",
        "QCD_Pt_600to800",
        "QCD_Pt_800to1000",
    ],
    "ST": [
        "ST_s-channel_4f_hadronicDecays",
        "ST_s-channel_4f_leptonDecays",
        "ST_t-channel_antitop_5f_InclusiveDecays",
        "ST_t-channel_top_5f_InclusiveDecays",
        "ST_tW_antitop_5f_inclusiveDecays",
        "ST_tW_top_5f_inclusiveDecays",
    ],
    "TTJets": [
        "TTJets",
        "TTJets_DiLept",
        "TTJets_DiLept_genMET-150",
        "TTJets_HT-1200to2500",
        "TTJets_HT-2500toInf",
        "TTJets_HT-600to800",
        "TTJets_HT-800to1200",
        "TTJets_SingleLeptFromT",
        "TTJets_SingleLeptFromT_genMET-150",
        "TTJets_SingleLeptFromTbar",
        "TTJets_SingleLeptFromTbar_genMET-150",
    ],
    "WJetsToLNu": [
        "WJetsToLNu_HT-1200To2500",
        "WJetsToLNu_HT-2500ToInf",
        "WJetsToLNu_HT-400To600",
        "WJetsToLNu_HT-600To800",
        "WJetsToLNu_HT-800To1200",
    ],
    "ZJetsToNuNu": [
        "ZJetsToNuNu_HT-1200To2500",
        "ZJetsToNuNu_HT-2500ToInf",
        "ZJetsToNuNu_HT-400To600",
        "ZJetsToNuNu_HT-600To800",
        "ZJetsToNuNu_HT-800To1200",
    ],
}

PROCESS_ORDER = ["QCD", "ST", "TTJets", "WJetsToLNu", "ZJetsToNuNu"]
STACK_ORDER = ["ST", "ZJetsToNuNu", "WJetsToLNu", "TTJets", "QCD"]
LEGEND_ORDER = ["QCD", "TTJets", "WJetsToLNu", "ZJetsToNuNu", "ST"]
DEFAULT_DATA_SAMPLES = ["HTMHT", "JetHT", "MET"]

PROC_LABEL = {
    "QCD": "QCD",
    "TTJets": "t#bar{t}+jets",
    "ZJetsToNuNu": "Z#rightarrow#nu#nu+jets",
    "WJetsToLNu": "W#rightarrowl#nu+jets",
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
    ROOT.TColor.GetColor("#118AB2"),  # teal-blue (colorblind-safe; was low-contrast light cyan #92dadd)
    ROOT.TColor.GetColor("#6b3e26"),
    ROOT.TColor.GetColor("#0b3d02"),
    ROOT.TColor.GetColor("#228833"),
    ROOT.TColor.GetColor("#1f4e79"),
]
SIGNAL_LINE_STYLES = [4, 1, 2, 1, 2]


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
        if cfg.get("log_bins"):
            return cls(np.geomspace(cfg["xmin"], cfg["xmax"], cfg["nbins"] + 1))
        return cls(np.linspace(cfg["xmin"], cfg["xmax"], cfg["nbins"] + 1))

    def fill(self, values, weights) -> None:
        values = np.asarray(ak.to_numpy(ak.fill_none(values, np.nan)), dtype=np.float64)
        weights = np.asarray(weights, dtype=np.float64)
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

    def add(self, other: "HistogramAccumulator") -> None:
        self.sumw += other.sumw
        self.sumw2 += other.sumw2

    def integral(self) -> float:
        return float(np.sum(self.sumw))


def empty_accumulators() -> Dict[str, HistogramAccumulator]:
    return {variable: HistogramAccumulator.for_variable(variable) for variable in VARIABLES}


def add_accumulators(target, source) -> None:
    for variable in VARIABLES:
        target[variable].add(source[variable])


def parse_csv(value: str) -> List[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def normalize_eos_path(path: str) -> str:
    path = os.path.expanduser(path).rstrip("/")
    if path.startswith("/eos/uscms/"):
        path = path[len("/eos/uscms") :]
    if not path.startswith("/"):
        raise ValueError(f"Expected absolute EOS path, got {path}")
    return path


def xrdfs_ls(path: str, recursive: bool = False) -> List[str]:
    command = ["xrdfs", EOS_HOST, "ls"]
    if recursive:
        command.append("-R")
    command.append(normalize_eos_path(path))
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"xrdfs failed for {path}: {result.stderr.strip()}")
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def eos_url(path: str) -> str:
    return EOS_HOST + "//" + normalize_eos_path(path).lstrip("/")


def get_part_files(sample_dir: str, max_files: int) -> List[str]:
    files = sorted(path for path in xrdfs_ls(sample_dir, recursive=True) if re.search(r"/part-[0-9]+\.root$", path))
    return files[:max_files] if max_files > 0 else files


def nominal_dir(base: str, year: str) -> str:
    return f"{normalize_eos_path(base)}/{year}/t_channel_pre_selection/nominal"


def variable_sources(variable: str) -> List[str]:
    cfg = VARIABLES[variable]
    if "numerator" in cfg:
        return [cfg["numerator"], cfg["denominator"]]
    if "log_source" in cfg:
        return [cfg["log_source"]]
    return [cfg["source"]]


def required_branches(tree, is_data: bool = False) -> List[str]:
    available = set(tree.keys())
    required = set() if is_data else {"Weight"}
    for optional in ("puWeight", "NonPrefiringProb", "JetsAK8_isGood"):
        if optional == "JetsAK8_isGood" and optional in available:
            required.add(optional)
        elif not is_data and optional in available:
            required.add(optional)
    for variable in VARIABLES:
        for branch in variable_sources(variable):
            if branch in available:
                required.add(branch)
    missing = sorted(branch for branch in required if branch not in available)
    if missing:
        raise KeyError(f"Events tree lacks required branch(es): {missing}")
    return sorted(required)


def first_good_jets(jagged, arrays, index: int):
    if "JetsAK8_isGood" in arrays.fields:
        jagged = jagged[arrays["JetsAK8_isGood"]]
    return ak.fill_none(ak.pad_none(jagged, index + 1, axis=1)[:, index], np.nan)


def log_positive(values):
    values = ak.values_astype(values, np.float64)
    return np.log(ak.where(values > 0.0, values, np.nan))


def variable_values(variable: str, arrays):
    cfg = VARIABLES[variable]
    if "numerator" in cfg:
        numerator = first_good_jets(arrays[cfg["numerator"]], arrays, cfg["jet_index"])
        denominator = first_good_jets(arrays[cfg["denominator"]], arrays, cfg["jet_index"])
        safe_denominator = ak.where(denominator != 0.0, denominator, np.nan)
        with np.errstate(divide="ignore", invalid="ignore"):
            return numerator / safe_denominator
    if "log_source" in cfg:
        return log_positive(arrays[cfg["log_source"]])
    if "jet_index" in cfg:
        values = first_good_jets(arrays[cfg["source"]], arrays, cfg["jet_index"])
        return log_positive(values) if cfg.get("log") else values
    values = arrays[cfg["source"]]
    return log_positive(values) if cfg.get("log") else values


def can_compute_variable(variable: str, arrays) -> bool:
    return all(branch in arrays.fields for branch in variable_sources(variable))


def stream_file(file_url: str, year: str, step_size: str, accumulators, is_data: bool = False) -> int:
    with uproot.open(file_url) as root_file:
        if "Events" not in root_file:
            raise KeyError("No Events tree found")
        tree = root_file["Events"]
        branches = required_branches(tree, is_data=is_data)
        entries = int(tree.num_entries)
        for arrays in tree.iterate(expressions=branches, step_size=step_size, library="ak"):
            if is_data:
                weights = np.ones(len(arrays), dtype=np.float64)
            else:
                weights = ak.to_numpy(arrays["Weight"]) * LUMI_PB[year]
                if "puWeight" in arrays.fields:
                    weights *= ak.to_numpy(arrays["puWeight"])
                if "NonPrefiringProb" in arrays.fields:
                    weights *= ak.to_numpy(arrays["NonPrefiringProb"])
            for variable in VARIABLES:
                if can_compute_variable(variable, arrays):
                    accumulators[variable].fill(variable_values(variable, arrays), weights)
    return entries


def ordered_signal_names(names: Sequence[str]) -> List[str]:
    default_index = {name: index for index, name in enumerate(DEFAULT_SIGNAL_DIRS)}
    return sorted(names, key=lambda name: (default_index.get(name, 10000), name))


def get_signal_directories(signal_base: str, year: str, requested: Sequence[str], strict: bool):
    base = nominal_dir(signal_base, year)
    entries = xrdfs_ls(base)
    available = {
        os.path.basename(entry.rstrip("/")): entry
        for entry in entries
        if os.path.basename(entry.rstrip("/")).startswith("t-channel_")
    }
    if len(requested) == 1 and requested[0].lower() in {"all", "*"}:
        return [(name, available[name]) for name in ordered_signal_names(available)]

    selected = []
    seen = set()
    for token in requested:
        matches = [name for name in available if token.lower() == name.lower() or token.lower() in name.lower()]
        if not matches:
            message = f"{year}: no t-channel signal matches '{token}'"
            if strict:
                raise RuntimeError(message)
            print(f"[WARN] {message}")
            continue
        for name in ordered_signal_names(matches):
            if name not in seen:
                selected.append((name, available[name]))
                seen.add(name)
    return selected


def mkdirs(root_file: ROOT.TFile, directory_path: str):
    directory = root_file
    for part in (item for item in directory_path.split("/") if item):
        next_dir = directory.GetDirectory(part)
        if not next_dir:
            next_dir = directory.mkdir(part)
        directory = next_dir
    return directory


def accumulator_to_histogram(variable: str, accumulator, name: Optional[str] = None, y_title: str = "Events"):
    cfg = VARIABLES[variable]
    hist = ROOT.TH1D(
        name or variable,
        f"{name or variable};{cfg['title']};{y_title}",
        len(accumulator.edges) - 1,
        array("d", accumulator.edges.tolist()),
    )
    hist.Sumw2()
    for ibin, (content, variance) in enumerate(zip(accumulator.sumw, accumulator.sumw2), start=1):
        hist.SetBinContent(ibin, float(content))
        hist.SetBinError(ibin, math.sqrt(max(0.0, float(variance))))
    hist.SetDirectory(0)
    return hist


def write_histogram(root_file: ROOT.TFile, directory_path: str, histogram, key_name: str) -> None:
    directory = mkdirs(root_file, directory_path)
    directory.cd()
    histogram.Write(key_name, ROOT.TObject.kOverwrite)


def write_cache(
    output_path,
    sample_hists,
    group_hists,
    run2_hists,
    signal_hists,
    run2_signal_hists,
    data_hists,
    run2_data_hists,
    years,
):
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    output = ROOT.TFile.Open(output_path, "RECREATE")
    if not output or output.IsZombie():
        raise RuntimeError(f"Could not create ROOT output: {output_path}")
    try:
        for (year, process, sample), accumulators in sorted(sample_hists.items()):
            for variable in VARIABLES:
                write_histogram(output, f"raw/samples/{year}/{process}/{sample}", accumulator_to_histogram(variable, accumulators[variable]), variable)
        for (year, process), accumulators in sorted(group_hists.items()):
            for variable in VARIABLES:
                write_histogram(output, f"raw/groups/{year}/{process}", accumulator_to_histogram(variable, accumulators[variable]), variable)
        for process, accumulators in sorted(run2_hists.items()):
            for variable in VARIABLES:
                write_histogram(output, f"raw/groups/Run2/{process}", accumulator_to_histogram(variable, accumulators[variable]), variable)
        for (year, signal_name), accumulators in sorted(signal_hists.items()):
            for variable in VARIABLES:
                write_histogram(output, f"raw/signals/{year}/{signal_name}", accumulator_to_histogram(variable, accumulators[variable]), variable)
        for signal_name, accumulators in sorted(run2_signal_hists.items()):
            for variable in VARIABLES:
                write_histogram(output, f"raw/signals/Run2/{signal_name}", accumulator_to_histogram(variable, accumulators[variable]), variable)

        for (year, sample), accumulators in sorted(data_hists.items()):
            for variable in VARIABLES:
                write_histogram(output, f"raw/data/{year}/{sample}", accumulator_to_histogram(variable, accumulators[variable]), variable)

        for year in years:
            year_data = empty_accumulators()
            found = False
            for (data_year, _), accumulators in data_hists.items():
                if data_year == year:
                    add_accumulators(year_data, accumulators)
                    found = True
            if found:
                for variable in VARIABLES:
                    write_histogram(output, f"raw/data/{year}/observed", accumulator_to_histogram(variable, year_data[variable]), variable)

        if run2_data_hists:
            for variable in VARIABLES:
                write_histogram(output, "raw/data/Run2/observed", accumulator_to_histogram(variable, run2_data_hists[variable]), variable)

        eras = {year: {process: group_hists[(year, process)] for process in PROCESS_ORDER if (year, process) in group_hists} for year in years}
        eras["Run2"] = run2_hists
        for era, processes in eras.items():
            for variable in VARIABLES:
                total = sum(accs[variable].integral() for accs in processes.values())
                if total <= 0.0:
                    continue
                for process, accs in processes.items():
                    hist = accumulator_to_histogram(variable, accs[variable], y_title="Arbitrary units")
                    hist.Scale(1.0 / total)
                    write_histogram(output, f"shapes/{era}/groups/{process}", hist, variable)

        era_signals = {
            year: {name: accs for (sig_year, name), accs in signal_hists.items() if sig_year == year}
            for year in years
        }
        era_signals["Run2"] = run2_signal_hists
        for era, signals in era_signals.items():
            for signal_name, accs in signals.items():
                for variable in VARIABLES:
                    hist = accumulator_to_histogram(variable, accs[variable], y_title="Arbitrary units")
                    integral = hist.Integral()
                    if integral <= 0.0:
                        continue
                    hist.Scale(1.0 / integral)
                    write_histogram(output, f"shapes/{era}/signals/{signal_name}", hist, variable)

        for era in list(years) + ["Run2"]:
            source = run2_data_hists if era == "Run2" else None
            if era != "Run2":
                source = empty_accumulators()
                found = False
                for (data_year, _), accs in data_hists.items():
                    if data_year == era:
                        add_accumulators(source, accs)
                        found = True
                if not found:
                    source = None
            if not source:
                continue
            for variable in VARIABLES:
                hist = accumulator_to_histogram(variable, source[variable], y_title="Arbitrary units")
                integral = hist.Integral()
                if integral <= 0.0:
                    continue
                hist.Scale(1.0 / integral)
                write_histogram(output, f"shapes/{era}/data/observed", hist, variable)
    finally:
        output.Close()


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
    desired = nbins / float(PLOT_TARGET_BINS)
    divisors = list(divisor_generator(nbins))
    factor = max(1, int(min(divisors, key=lambda value: abs(value - desired))))
    visible_bins = nbins // factor if factor > 0 else nbins
    if visible_bins > 10 and nbins % (factor * DISPLAY_BIN_REDUCTION) == 0:
        factor *= DISPLAY_BIN_REDUCTION
    return factor


def clone_hist(histogram, name: str):
    out = histogram.Clone(name)
    out.SetDirectory(0)
    return out


def get_hist(root_file, path: str, clone_name: str):
    obj = root_file.Get(path)
    if not obj or not obj.InheritsFrom("TH1"):
        return None
    return clone_hist(obj, clone_name)


def rebin_for_display(histogram, variable: str, name: str):
    out = clone_hist(histogram, name)
    configured_bins = int(VARIABLES[variable]["nbins"])
    current_bins = int(out.GetNbinsX())
    factor = rebin_factor(variable)
    if current_bins != configured_bins and current_bins % configured_bins == 0:
        factor *= current_bins // configured_bins
    if current_bins % factor != 0:
        factor = 1
    if factor > 1:
        out = out.Rebin(factor, f"{name}_rebin{factor}")
        out.SetDirectory(0)
    return out


def signal_label(signal_name: str) -> str:
    mmed = re.search(r"mMed-([0-9]+(?:p[0-9]+)?)", signal_name, re.IGNORECASE)
    rinv = re.search(r"rinv-([0-9]+(?:p[0-9]+)?)", signal_name, re.IGNORECASE)
    if not (mmed and rinv):
        return signal_name
    return f"m_{{#phi}} = {mmed.group(1).replace('p', '.')} GeV, r_{{inv}} = {rinv.group(1).replace('p', '.')}"


def signal_scale_factor(signal_name: str) -> float:
    mmed = re.search(r"mMed-([0-9]+(?:p[0-9]+)?)", signal_name, re.IGNORECASE)
    rinv = re.search(r"rinv-([0-9]+(?:p[0-9]+)?)", signal_name, re.IGNORECASE)
    if not (mmed and rinv):
        return 1.0
    key = (mmed.group(1).replace("p", "."), rinv.group(1).replace("p", "."))
    return SIGNAL_SCALE_FACTORS.get(key, 1.0)


def scaled_signal_label(signal_name: str, normalized: bool) -> str:
    label = signal_label(signal_name)
    scale = signal_scale_factor(signal_name)
    if normalized or abs(scale - 1.0) < 1.0e-9:
        return label
    return f"{label} #times {scale:g}"


def style_background(histogram, process: str) -> None:
    color = PROC_COLOR[process]
    histogram.SetFillColor(color)
    histogram.SetFillStyle(1001)
    histogram.SetLineColor(color)
    histogram.SetLineWidth(0)
    histogram.SetMarkerSize(0)


def style_signal(histogram, index: int) -> None:
    histogram.SetLineColor(SIGNAL_LINE_COLORS[index % len(SIGNAL_LINE_COLORS)])
    histogram.SetLineStyle(SIGNAL_LINE_STYLES[index % len(SIGNAL_LINE_STYLES)])
    histogram.SetLineWidth(3)
    histogram.SetFillStyle(0)
    histogram.SetMarkerSize(0)


def style_data(histogram) -> None:
    histogram.SetLineColor(ROOT.kBlack)
    histogram.SetMarkerColor(ROOT.kBlack)
    histogram.SetMarkerStyle(ROOT.kFullCircle)
    histogram.SetMarkerSize(0.85)
    histogram.SetLineWidth(2)
    histogram.SetFillStyle(0)


def build_ratio_histograms(data_histogram, mc_histogram, name: str):
    ratio = clone_hist(data_histogram, f"{name}_ratio")
    ratio.Reset("ICESM")
    band = clone_hist(mc_histogram, f"{name}_ratio_band")
    band.Reset("ICESM")
    for ibin in range(1, data_histogram.GetNbinsX() + 1):
        data_value = data_histogram.GetBinContent(ibin)
        data_error = data_histogram.GetBinError(ibin)
        mc_value = mc_histogram.GetBinContent(ibin)
        mc_error = mc_histogram.GetBinError(ibin)
        if mc_value <= 0.0:
            ratio.SetBinContent(ibin, 0.0)
            ratio.SetBinError(ibin, 0.0)
            band.SetBinContent(ibin, 0.0)
            band.SetBinError(ibin, 0.0)
            continue
        ratio.SetBinContent(ibin, data_value / mc_value)
        ratio.SetBinError(ibin, data_error / mc_value)
        band.SetBinContent(ibin, 1.0)
        band.SetBinError(ibin, mc_error / mc_value)

    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetMarkerColor(ROOT.kBlack)
    ratio.SetMarkerStyle(ROOT.kFullCircle)
    ratio.SetMarkerSize(0.75)
    ratio.SetLineWidth(2)
    ratio.SetFillStyle(0)

    band.SetLineColor(ROOT.kGray + 2)
    band.SetFillColor(ROOT.kGray + 1)
    band.SetFillStyle(3004)
    band.SetMarkerSize(0)
    return ratio, band


def lumi_label(era: str) -> Optional[str]:
    if era == "Run2":
        return "138 fb^{-1} (13 TeV)"
    if era in LUMI_PB:
        return f"{LUMI_PB[era] / 1000.0:.1f} fb^{{-1}} (13 TeV)"
    return None


DRAW_PRELIMINARY_LABEL = False


def draw_cms_label(simulation: bool = True, era: Optional[str] = None, show_lumi: bool = False) -> None:
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextColor(ROOT.kBlack)
    latex.SetTextFont(61)
    latex.SetTextSize(0.060)
    latex.SetTextAlign(11)
    label_y = 0.918
    latex.DrawLatex(0.16, label_y, "CMS")
    if simulation:
        latex.SetTextFont(52)
        latex.SetTextSize(0.040)
        latex.DrawLatex(0.290, label_y, "Simulation")
        if DRAW_PRELIMINARY_LABEL:
            latex.DrawLatex(0.455, label_y, "Preliminary")
    elif DRAW_PRELIMINARY_LABEL:
        latex.SetTextFont(52)
        latex.SetTextSize(0.040)
        latex.DrawLatex(0.255, label_y, "Preliminary")
    latex.SetTextFont(42)
    latex.SetTextSize(0.045)
    latex.SetTextAlign(31)
    latex.DrawLatex(0.95, label_y, lumi_label(era) if show_lumi and era else "(13 TeV)")


def draw_common_signal_text(y_position: float = 0.650) -> None:
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.030)
    latex.SetTextColor(ROOT.kBlack)
    latex.SetTextAlign(12)
    latex.DrawLatex(0.565, y_position, "m_{dark} = 20 GeV, #kern[-0.18]{#lambda} = 1")


def interpolate_axis_y(ymin: float, ymax: float, fraction: float, log_y: bool) -> float:
    fraction = max(0.0, min(1.0, fraction))
    if log_y and ymin > 0.0 and ymax > ymin:
        return math.exp(math.log(ymin) + fraction * (math.log(ymax) - math.log(ymin)))
    return ymin + fraction * (ymax - ymin)


def draw_working_point_marker(variable: str, axis, log_y: bool):
    cfg = VARIABLES[variable]
    if "wp" not in cfg:
        return None
    x = float(cfg["wp"])
    xmin = float(cfg["xmin"])
    xmax = float(cfg["xmax"])
    if x <= xmin or x >= xmax:
        return None
    ymin = axis.GetMinimum()
    ymax = axis.GetMaximum()
    line_top = interpolate_axis_y(ymin, ymax, float(cfg.get("wp_line_y_fraction", 0.82)), log_y)
    arrow_y = interpolate_axis_y(ymin, ymax, float(cfg.get("wp_arrow_y_fraction", 0.74)), log_y)

    line = ROOT.TLine(x, ymin, x, line_top)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineStyle(ROOT.kDashed)
    line.SetLineWidth(2)
    line.Draw("same")

    arrow = ROOT.TLatex()
    arrow.SetTextColor(ROOT.kBlack)
    arrow.SetTextFont(42)
    arrow.SetTextSize(0.040)
    if cfg.get("wp_side", "right") == "left":
        arrow.SetTextAlign(32)
        arrow.DrawLatex(x, arrow_y, "#leftarrow")
    else:
        arrow.SetTextAlign(12)
        arrow.DrawLatex(x, arrow_y, "#rightarrow")
    return line, arrow


def cached_signal_names(root_file, era: str) -> List[str]:
    directory = root_file.GetDirectory(f"shapes/{era}/signals")
    if not directory:
        return []
    return ordered_signal_names([key.GetName() for key in directory.GetListOfKeys()])


def choose_signals(root_file, era: str, requested: Sequence[str]) -> List[str]:
    available = set(cached_signal_names(root_file, era))
    chosen = []
    for token in requested:
        if token in available:
            chosen.append(token)
            continue
        matches = sorted(name for name in available if token.lower() in name.lower())
        if len(matches) == 1:
            chosen.append(matches[0])
        elif not matches:
            print(f"[WARN] {era}: cached signal not found: {token}")
        else:
            print(f"[WARN] {era}: token '{token}' matches {len(matches)} signals; use full name")
    out = []
    seen = set()
    for name in chosen:
        if name not in seen:
            out.append(name)
            seen.add(name)
    return out


def draw_plot(
    root_file,
    era: str,
    variable: str,
    requested_signals,
    output_base: str,
    also_linear: bool,
    include_data: bool = False,
    normalized: bool = True,
) -> None:
    directory_kind = "shapes" if normalized else "raw"
    y_title = "Arbitrary units" if normalized else "Events"
    group_prefix = f"shapes/{era}/groups" if normalized else f"raw/groups/{era}"
    signal_prefix = f"shapes/{era}/signals" if normalized else f"raw/signals/{era}"
    data_prefix = f"shapes/{era}/data" if normalized else f"raw/data/{era}"
    process_hists = {}
    for process in STACK_ORDER:
        raw = get_hist(root_file, f"{group_prefix}/{process}/{variable}", f"{era}_{directory_kind}_{process}_{variable}")
        if not raw:
            continue
        hist = rebin_for_display(raw, variable, f"display_{era}_{process}_{variable}")
        if hist.Integral() <= 0.0:
            continue
        style_background(hist, process)
        process_hists[process] = hist
    if not process_hists:
        print(f"[WARN] {era}/{variable}: no {directory_kind} background histograms found")
        return

    signal_hists = []
    for index, signal_name in enumerate(choose_signals(root_file, era, requested_signals)):
        raw = get_hist(root_file, f"{signal_prefix}/{signal_name}/{variable}", f"{era}_{directory_kind}_{signal_name}_{variable}")
        if not raw:
            continue
        hist = rebin_for_display(raw, variable, f"display_{era}_signal_{index}_{variable}")
        if not normalized:
            hist.Scale(signal_scale_factor(signal_name))
        if hist.Integral() <= 0.0:
            continue
        style_signal(hist, index)
        signal_hists.append((signal_name, hist))

    data_hist = None
    if include_data:
        raw_data = get_hist(root_file, f"{data_prefix}/observed/{variable}", f"{era}_{directory_kind}_data_{variable}")
        if raw_data and raw_data.Integral() > 0.0:
            data_hist = rebin_for_display(raw_data, variable, f"display_{era}_data_{variable}")
            style_data(data_hist)
        else:
            print(f"[WARN] {era}/{variable}: no {directory_kind} data histogram found")

    for log_y in ([True, False] if also_linear else [True]):
        scale_tag = "normalized" if normalized else "raw"
        data_tag = "_wdata_ratio" if data_hist is not None else ""
        tag = f"{'log' if log_y else 'linear'}_{scale_tag}{data_tag}"
        stack = ROOT.THStack(f"stack_{era}_{variable}_{tag}", "")
        total = None
        for process in STACK_ORDER:
            hist = process_hists.get(process)
            if not hist:
                continue
            stack.Add(hist)
            if total is None:
                total = clone_hist(hist, f"total_{era}_{variable}_{tag}")
            else:
                total.Add(hist)

        has_ratio = data_hist is not None
        canvas = ROOT.TCanvas(f"c_{era}_{variable}_{tag}", "", 800, 900 if has_ratio else 800)
        if has_ratio:
            top_pad = ROOT.TPad(f"top_{era}_{variable}_{tag}", "", 0.0, 0.30, 1.0, 1.0)
            bottom_pad = ROOT.TPad(f"ratio_{era}_{variable}_{tag}", "", 0.0, 0.0, 1.0, 0.30)
            top_pad.SetLeftMargin(0.16)
            top_pad.SetRightMargin(0.05)
            top_pad.SetTopMargin(0.09)
            top_pad.SetBottomMargin(0.03)
            bottom_pad.SetLeftMargin(0.16)
            bottom_pad.SetRightMargin(0.05)
            bottom_pad.SetTopMargin(0.02)
            bottom_pad.SetBottomMargin(0.34)
            top_pad.SetTicks(1, 1)
            bottom_pad.SetTicks(1, 1)
            canvas.cd()
            top_pad.Draw()
            bottom_pad.Draw()
            top_pad.cd()
            ROOT.gPad.SetLogy(log_y)
        else:
            canvas.cd()
            ROOT.gPad.SetLeftMargin(0.16)
            ROOT.gPad.SetRightMargin(0.05)
            ROOT.gPad.SetTopMargin(0.08)
            ROOT.gPad.SetBottomMargin(0.12)
            ROOT.gPad.SetTicks(1, 1)
            ROOT.gPad.SetLogy(log_y)

        axis = clone_hist(total, f"axis_{era}_{variable}_{tag}")
        axis.Reset("ICESM")
        axis.SetStats(0)
        axis.SetTitle("")
        axis.GetXaxis().SetTitle(VARIABLES[variable]["title"])
        axis.GetYaxis().SetTitle(y_title)
        axis.GetXaxis().SetTitleSize(0.0 if has_ratio else 0.05)
        axis.GetXaxis().SetLabelSize(0.0 if has_ratio else 0.04)
        axis.GetYaxis().SetTitleSize(0.05 if has_ratio else 0.05)
        axis.GetYaxis().SetLabelSize(0.04 if has_ratio else 0.04)
        axis.GetYaxis().SetLabelOffset(0.01)
        axis.GetYaxis().SetTitleOffset(1.30)
        axis.GetXaxis().SetRangeUser(VARIABLES[variable]["xmin"], VARIABLES[variable]["xmax"])
        ymax = max([total.GetMaximum()] + [hist.GetMaximum() for _, hist in signal_hists] + ([data_hist.GetMaximum()] if data_hist else []) + [1.0])
        axis.SetMinimum(1.0 if log_y and not normalized else 1.0e-2)
        y_max = ymax * (100.0 if log_y else 2.6)
        if log_y and not normalized:
            y_max = max(y_max, 1.0e8)
        axis.SetMaximum(y_max)
        axis.Draw("hist")
        stack.Draw("hist same")
        for _, hist in signal_hists:
            hist.Draw("hist same")
        if data_hist:
            data_hist.Draw("E1 X0 same")
        wp_marker = draw_working_point_marker(variable, axis, log_y)

        if log_y and not normalized:
            legend_y1 = 0.615
            legend = ROOT.TLegend(0.23, legend_y1, 0.97, 0.855)
        else:
            legend_y1 = 0.660
            legend = ROOT.TLegend(0.23, legend_y1, 0.97, 0.875)
        legend.SetNColumns(2)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.034)
        legend.SetTextFont(42)
        legend.SetColumnSeparation(0.035 if log_y and not normalized else 0.050)
        legend.SetEntrySeparation(0.044 if log_y and not normalized else 0.046)
        bkg_entries = [(process_hists[p], PROC_LABEL[p], "F") for p in LEGEND_ORDER if p in process_hists]
        if data_hist:
            bkg_entries.insert(0, (data_hist, "Data", "EP"))
        sig_entries = [(hist, scaled_signal_label(name, normalized), "L") for name, hist in signal_hists]
        if signal_hists:
            sig_entries.append(("", "m_{dark} = 20 GeV, #kern[-0.18]{#lambda} = 1", ""))
        for row in range(max(len(bkg_entries), len(sig_entries))):
            legend.AddEntry(*(bkg_entries[row] if row < len(bkg_entries) else ("", "", "")))
            legend.AddEntry(*(sig_entries[row] if row < len(sig_entries) else ("", "", "")))
        legend.Draw()
        draw_cms_label(simulation=not include_data, era=era, show_lumi=True)
        ROOT.gPad.RedrawAxis()

        if has_ratio:
            ratio_hist, ratio_band = build_ratio_histograms(data_hist, total, f"{era}_{variable}_{tag}")
            bottom_pad.cd()
            ratio_hist.SetStats(0)
            ratio_hist.SetTitle("")
            ratio_hist.GetXaxis().SetTitle(VARIABLES[variable]["title"])
            ratio_hist.GetYaxis().SetTitle("Data/Sim")
            ratio_hist.GetXaxis().SetTitleSize(0.12)
            ratio_hist.GetXaxis().SetLabelSize(0.10)
            ratio_hist.GetYaxis().SetTitleSize(0.10)
            ratio_hist.GetYaxis().SetLabelSize(0.08)
            ratio_hist.GetYaxis().SetTitleOffset(0.58)
            ratio_hist.GetYaxis().SetNdivisions(505)
            ratio_hist.GetXaxis().SetRangeUser(VARIABLES[variable]["xmin"], VARIABLES[variable]["xmax"])
            ratio_ymin, ratio_ymax = (0.0, 2.0)
            ratio_hist.SetMinimum(ratio_ymin)
            ratio_hist.SetMaximum(ratio_ymax)
            ratio_hist.Draw("E1 X0")
            ratio_band.Draw("E2 same")
            ratio_hist.Draw("E1 X0 same")
            line = ROOT.TLine(VARIABLES[variable]["xmin"], 1.0, VARIABLES[variable]["xmax"], 1.0)
            line.SetLineColor(ROOT.kBlack)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            line.Draw("same")
            bottom_pad.RedrawAxis()

        os.makedirs(os.path.dirname(output_base), exist_ok=True)
        canvas.SaveAs(f"{output_base}_{tag}.pdf")
        canvas.SaveAs(f"{output_base}_{tag}.png")
        canvas.Close()
        print(f"[OK] Wrote {output_base}_{tag}.pdf/.png")


def plot_cache(
    output_path: str,
    plot_dir: str,
    eras: Sequence[str],
    signals: Sequence[str],
    also_linear: bool,
    draw_raw: bool,
    raw_only: bool = False,
) -> None:
    configure_figure2_plot_style()
    root_file = ROOT.TFile.Open(output_path, "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Could not open cache ROOT file: {output_path}")
    try:
        for era in eras:
            for variable in VARIABLES:
                if not raw_only:
                    f2.draw_stack(
                        root_file,
                        era,
                        variable,
                        os.path.join(plot_dir, era, "mc_signal_only", variable),
                        signals,
                        normalized=True,
                        also_linear=also_linear,
                        include_data=False,
                        include_ratio=False,
                        preliminary_label=DRAW_PRELIMINARY_LABEL,
                    )
                    f2.draw_stack(
                        root_file,
                        era,
                        variable,
                        os.path.join(plot_dir, era, "with_data", variable),
                        signals,
                        normalized=True,
                        also_linear=also_linear,
                        include_data=True,
                        include_ratio=False,
                        preliminary_label=DRAW_PRELIMINARY_LABEL,
                    )
                if draw_raw or raw_only:
                    f2.draw_stack(
                        root_file,
                        era,
                        variable,
                        os.path.join(plot_dir, era, "raw_wdata_ratio", variable),
                        signals,
                        include_data=True,
                        normalized=False,
                        also_linear=also_linear,
                        include_ratio=True,
                        preliminary_label=DRAW_PRELIMINARY_LABEL,
                    )
    finally:
        root_file.Close()


def read_sample(sample_dir, display_name, year, step_size, max_files, strict, dry_run, failures, is_data: bool = False):
    try:
        files = get_part_files(sample_dir, max_files)
    except Exception as exc:
        message = f"{display_name}: {exc}"
        if strict:
            raise RuntimeError(message) from exc
        print(f"[WARN] {message}")
        return None
    print(f"       part files: {len(files)}")
    if not files:
        return None
    if dry_run:
        return empty_accumulators()

    accumulators = empty_accumulators()
    entries = 0
    successful = 0
    for index, remote_file in enumerate(files, start=1):
        try:
            entries += stream_file(eos_url(remote_file), year, step_size, accumulators, is_data=is_data)
            successful += 1
        except Exception as exc:
            message = f"{display_name}/{os.path.basename(remote_file)}: {exc}"
            failures.append(message)
            print(f"[WARN] {message}")
        if index % 25 == 0 or index == len(files):
            print(f"       processed {index:4d}/{len(files):4d} files")
    if strict and successful != len(files):
        raise RuntimeError(f"{display_name}: {len(files) - successful} input file(s) failed")
    yield_label = "events" if is_data else "yield"
    print(f"       Events read: {entries:,} | z12 {yield_label}: {accumulators['LundJetPlaneZ01GoodJetsAK8'].integral():.8g}")
    return accumulators


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--background-base", default=DEFAULT_BACKGROUND_BASE)
    parser.add_argument("--data-base", default=DEFAULT_DATA_BASE)
    parser.add_argument("--data-samples", default=",".join(DEFAULT_DATA_SAMPLES))
    parser.add_argument("--no-data", action="store_true", help="skip data streaming and only draw MC/signal plots")
    parser.add_argument("--signal-base", default=DEFAULT_SIGNAL_BASE)
    parser.add_argument("--years", default="2016,2017,2018")
    parser.add_argument("--signals", default=",".join(DEFAULT_SIGNAL_DIRS))
    parser.add_argument("--output", default="EventDNNInput_plots/event_dnn_inputs.root")
    parser.add_argument("--plot-dir", default="EventDNNInput_plots")
    parser.add_argument("--chunk-size", default="100 MB")
    parser.add_argument("--max-files", type=int, default=-1)
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--strict", action="store_true")
    parser.add_argument("--also-linear", action="store_true")
    parser.add_argument("--draw-raw", action="store_true", help="also draw unnormalized weighted-yield stacks from the raw cache")
    parser.add_argument("--raw-only", action="store_true", help="only draw raw data/MC-ratio plots")
    parser.add_argument("--eras", default="", help="plot-only eras; default is years plus Run2")
    parser.add_argument("--preliminary-label", action="store_true", help="draw CMS Preliminary or CMS Simulation Preliminary")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    global DRAW_PRELIMINARY_LABEL
    DRAW_PRELIMINARY_LABEL = args.preliminary_label
    years = parse_csv(args.years)
    signals = parse_csv(args.signals)
    data_samples = parse_csv(args.data_samples)
    eras = parse_csv(args.eras) if args.eras else list(years) + ["Run2"]

    if args.plot_only:
        plot_cache(args.output, args.plot_dir, eras, signals, args.also_linear, args.draw_raw, args.raw_only)
        return 0

    sample_hists = {}
    group_hists = {}
    run2_hists = {process: empty_accumulators() for process in PROCESS_ORDER}
    signal_hists = {}
    run2_signal_hists = {}
    data_hists = {}
    run2_data_hists = empty_accumulators() if not args.no_data and data_samples else None
    failures = []

    for year in years:
        print(f"\n[INFO] Backgrounds {year}")
        for process in PROCESS_ORDER:
            group_hists[(year, process)] = empty_accumulators()
            for sample in SAMPLE_MANIFEST[process]:
                sample_dir = f"{nominal_dir(args.background_base, year)}/{sample}"
                print(f"  {process}/{sample}")
                accs = read_sample(sample_dir, f"{year}/{process}/{sample}", year, args.chunk_size, args.max_files, args.strict, args.dry_run, failures)
                if accs is None:
                    continue
                sample_hists[(year, process, sample)] = accs
                add_accumulators(group_hists[(year, process)], accs)
                add_accumulators(run2_hists[process], accs)

        print(f"\n[INFO] Signals {year}")
        for signal_name, signal_dir in get_signal_directories(args.signal_base, year, signals, args.strict):
            print(f"  {signal_name}")
            accs = read_sample(signal_dir, f"{year}/{signal_name}", year, args.chunk_size, args.max_files, args.strict, args.dry_run, failures)
            if accs is None:
                continue
            signal_hists[(year, signal_name)] = accs
            if signal_name not in run2_signal_hists:
                run2_signal_hists[signal_name] = empty_accumulators()
            add_accumulators(run2_signal_hists[signal_name], accs)

        if not args.no_data and data_samples:
            print(f"\n[INFO] Data {year}")
            for sample in data_samples:
                sample_dir = f"{nominal_dir(args.data_base, year)}/{sample}"
                print(f"  {sample}")
                accs = read_sample(
                    sample_dir,
                    f"{year}/data/{sample}",
                    year,
                    args.chunk_size,
                    args.max_files,
                    args.strict,
                    args.dry_run,
                    failures,
                    is_data=True,
                )
                if accs is None:
                    continue
                data_hists[(year, sample)] = accs
                if run2_data_hists is not None:
                    add_accumulators(run2_data_hists, accs)

    write_cache(
        args.output,
        sample_hists,
        group_hists,
        run2_hists,
        signal_hists,
        run2_signal_hists,
        data_hists,
        run2_data_hists,
        years,
    )
    plot_cache(args.output, args.plot_dir, eras, signals, args.also_linear, args.draw_raw, args.raw_only)

    if failures:
        print(f"\n[WARN] Completed with {len(failures)} file failure(s). First few:")
        for failure in failures[:10]:
            print(f"  - {failure}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[FATAL] {exc}", file=sys.stderr)
        raise
