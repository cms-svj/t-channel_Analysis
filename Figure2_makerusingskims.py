#!/usr/bin/env python3
"""
Create a compact ROOT cache and normalized log-y MC+signal stacks.

Background source:
  /store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/data_mc_WNAE/
      YEAR/t_channel_pre_selection/nominal/<background-sample>/part-*.root

Signal source:
  /store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/
      YEAR/t_channel_pre_selection/nominal/<t-channel signal>/part-*.root

No event trees are copied. The output ROOT file contains only TH1D
histograms, both raw luminosity-scaled histograms and the normalized shapes
used for the log plots.

Raw cache:
  raw/samples/<year>/<process>/<sample>/{MET,ST}
  raw/groups/<year>/<process>/{MET,ST}
  raw/groups/Run2/<process>/{MET,ST}
  raw/signals/<year>/<signal>/{MET,ST}
  raw/signals/Run2/<signal>/{MET,ST}

Normalized plotting cache:
  shapes/<era>/groups/<process>/{MET,ST}
  shapes/<era>/signals/<signal>/{MET,ST}

The normalized convention intentionally matches the legacy plotStack2 MC-only
signal mode:
  * Each background component is divided by the total background integral.
    Therefore the whole stacked MC distribution integrates to one.
  * Each signal is normalized to unit area independently.
  * Plots use log y, "Arbitrary units", and the same two-column legend style.

Once the ROOT cache exists, rerun with --plot-only to regenerate PDF/PNG plots
without reading EOS skims again.

 python3 -u Figure2_makerusingskims.py \
  --background-base /store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/data_mc_WNAE \
  --data-base /store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/data_mc_WNAE \
  --signal-base /store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto \
  --years 2016,2017,2018 \
  --output Figure2_plots/WNAE_selected_MC_vars_Run2.root \
  --plot-dir Figure2_plots \
  --chunk-size "100 MB" \
  "$@" \
  |& tee "logs/Figure2_MET_ST_$(date +%Y%m%d_%H%M%S).log"


PLOT ONLY 
cd /uscms/home/ashrivas/nobackup/Dark_Sector/t-channel_Analysis
source condor/initCondor.sh

python3 -u Figure2_makerusingskims.py \
  --plot-only \
  --years 2016,2017,2018 \
  --output Figure2_plots/noWNAE_selected_MC_vars_Run2.root \
  --plot-dir Figure2_plots
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


# =============================================================================
# Input locations and physics configuration
# =============================================================================
EOS_HOST = "root://cmseos.fnal.gov"

DEFAULT_BACKGROUND_BASE = (
    "/store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/data_mc_WNAE"
)
DEFAULT_DATA_BASE = DEFAULT_BACKGROUND_BASE

# Signals are NOT under the data/MC background base. They are directly under this base.
DEFAULT_SIGNAL_BASE = "/store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto"
DEFAULT_WNAE_SIGNAL_BASE = "/store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/signals_WNAE"

LUMI_PB = {
    "2016": 36.31 * 1000.0,
    "2017": 42.07 * 1000.0,
    "2018": 59.56 * 1000.0,
}

VARIABLES = {
    "MET": {
        "branch": "MET",
        "title": "p_{T}^{miss} [GeV]",
        "x_title": "p_{T}^{miss} [GeV]",
        "nbins": 500,
        "xmin": 200.0,
        "xmax": 2000.0,
    },
    "HT": {
        "branch": "HT",
        "title": "H_{T} [GeV]",
        "x_title": "H_{T} [GeV]",
        "nbins": 500,
        "xmin": 0.0,
        "xmax": 5000.0,
    },
    "ST": {
        "branch": "ST",
        "title": "S_{T} [GeV]",
        "x_title": "S_{T} [GeV]",
        "nbins": 500,
        "xmin": 1300.0,
        "xmax": 5000.0,
    },
    "mT": {
        "branch": "MT_AK8",
        "title": "m_{T} [GeV]",
        "x_title": "m_{T} [GeV]",
        "nbins": 500,
        "xmin": 0.0,
        "xmax": 6000.0,
    },
    "dPhiMinjMETAK8": {
        "branch": "DeltaPhiMinGoodJetsAK8",
        "title": "#Delta#phi_{min}(J, p_{T}^{miss})",
        "x_title": "#Delta#phi_{min}(J, p_{T}^{miss})",
        "nbins": 100,
        "xmin": 0.0,
        "xmax": 2.0,
        "display_xmax": 1.5,
    },
    "nsvjJetsAK8": {
        "branch": "nsvjJetsAK8",
        "title": "Number of SVJ AK8 jets",
        "x_title": "Number of SVJ AK8 jets",
        "nbins": 20,
        "xmin": 0.0,
        "xmax": 5.0,
    },
    "dnnEventClassScore": {
        "branch": "dnnEventClassScore",
        "title": "Event classifier score",
        "x_title": "Event classifier score",
        "nbins": 100,
        "xmin": 0.0,
        "xmax": 1.0,
        "wp": 0.85,
        "wp_side": "right",
    },
}

JET_VARIABLES = {
    "PtAK8": ("JetsAK8_/.fPt", "p_{T}(J_{{idx}}) [GeV]", 400, 0.0, 4000.0),
    "EtaAK8": ("JetsAK8_/.fEta", "#eta(J_{{idx}})", 80, -3.0, 3.0),
    "PhiAK8": ("JetsAK8_/.fPhi", "#phi(J_{{idx}})", 64, -3.2, 3.2),
    "MassAK8": ("JetsAK8_mass", "m(J_{{idx}}) [GeV]", 200, 0.0, 900.0),
    "SoftDropMassAK8": ("JetsAK8_softDropMass", "m_{SD}(J_{{idx}}) [GeV]", 200, 0.0, 900.0),
    "MTMETAK8": ("JetsAK8_MTMET", "M_{T}(J_{{idx}}, p_{T}^{miss}) [GeV]", 100, 0.0, 3000.0),
    "DeltaPhiMETAK8": ("JetsAK8_deltaPhiMET", "#Delta#phi(J_{{idx}}, p_{T}^{miss})", 64, -3.2, 3.2),
    "LundJetPlaneZAK8": ("JetsAK8_LundJetPlaneZ", "z(J_{{idx}}, p_{T}^{miss})", 60, 0.0, 0.5),
    "Tau1AK8": ("JetsAK8_NsubjettinessTau1", "#tau_{1}(J_{{idx}})", 40, 0.0, 0.7),
    "Tau2AK8": ("JetsAK8_NsubjettinessTau2", "#tau_{2}(J_{{idx}})", 40, 0.0, 0.6),
    "Tau3AK8": ("JetsAK8_NsubjettinessTau3", "#tau_{3}(J_{{idx}})", 40, 0.0, 0.4),
    "Tau4AK8": ("JetsAK8_NsubjettinessTau4", "#tau_{4}(J_{{idx}})", 40, 0.0, 0.4),
    "Tau5AK8": ("JetsAK8_NsubjettinessTau5", "#tau_{5}(J_{{idx}})", 40, 0.0, 0.3),
    "NConstituentsSoftDropAK8": ("JetsAK8_nConstituentsSoftDrop", "N_{const}^{SD}(J_{{idx}})", 80, 0.0, 200.0),
    "PNetScoreAK8": ("JetsAK8_pNetJetTaggerScore", "ParticleNet score (J_{{idx}})", 10, 0.0, 1.0),
    "GirthAK8": ("JetsAK8_girth", "girth(J_{{idx}})", 40, 0.0, 0.6),
    "PtDAK8": ("JetsAK8_ptD", "p_{T}^{D}(J_{{idx}})", 40, 0.05, 1.0),
    "AxisMajorAK8": ("JetsAK8_axismajor", "axis major(J_{{idx}})", 40, 0.0, 0.5),
    "AxisMinorAK8": ("JetsAK8_axisminor", "axis minor(J_{{idx}})", 40, 0.0, 0.4),
    "EcfC2b1AK8": ("JetsAK8_ecfC2b1", "ECF C_{2}^{#beta=1}(J_{{idx}})", 40, 0.0, 0.5),
    "EcfC2b2AK8": ("JetsAK8_ecfC2b2", "ECF C_{2}^{#beta=2}(J_{{idx}})", 40, 0.0, 0.5),
    "EcfD2b1AK8": ("JetsAK8_ecfD2b1", "ECF D_{2}^{#beta=1}(J_{{idx}})", 40, 0.0, 5.0),
    "EcfD2b2AK8": ("JetsAK8_ecfD2b2", "ECF D_{2}^{#beta=2}(J_{{idx}})", 40, 0.0, 7.0),
    "EcfN2b1AK8": ("JetsAK8_ecfN2b1", "ECF N_{2}^{#beta=1}(J_{{idx}})", 40, 0.0, 0.5),
    "EcfN2b2AK8": ("JetsAK8_ecfN2b2", "ECF N_{2}^{#beta=2}(J_{{idx}})", 40, 0.0, 0.4),
}

WNAE_PT_LOSS_VARIABLES = {
    "WNAEPt0To200LossAK8": (
        "JetsAK8_WNAEPt0To200Loss",
        "WNAE (0 < p_{T} < 200 GeV) score (J_{{idx}})",
        80,
        0.0,
        200.0,
    ),
    "WNAEPt200To300LossAK8": (
        "JetsAK8_WNAEPt200To300Loss",
        "WNAE (200 < p_{T} < 300 GeV) score (J_{{idx}})",
        80,
        0.0,
        200.0,
    ),
    "WNAEPt300To400LossAK8": (
        "JetsAK8_WNAEPt300To400Loss",
        "WNAE (300 < p_{T} < 400 GeV) score (J_{{idx}})",
        80,
        0.0,
        200.0,
    ),
    "WNAEPt400To500LossAK8": (
        "JetsAK8_WNAEPt400To500Loss",
        "WNAE (400 < p_{T} < 500 GeV) score (J_{{idx}})",
        80,
        0.0,
        200.0,
    ),
    "WNAEPt500ToInfLossAK8": (
        "JetsAK8_WNAEPt500ToInfLoss",
        "WNAE (p_{T} > 500 GeV) score (J_{{idx}})",
        80,
        0.0,
        200.0,
    ),
}

WNAE_LOSS_BRANCHES = [
    ("JetsAK8_WNAEPt0To200Loss", 0.0, 200.0),
    ("JetsAK8_WNAEPt200To300Loss", 200.0, 300.0),
    ("JetsAK8_WNAEPt300To400Loss", 300.0, 400.0),
    ("JetsAK8_WNAEPt400To500Loss", 400.0, 500.0),
    ("JetsAK8_WNAEPt500ToInfLoss", 500.0, None),
]

WNAE_MC_WORKING_POINTS = {
    "WNAEPt0To200LossAK8": 25.156,
    "WNAEPt200To300LossAK8": 18.284,
    "WNAEPt300To400LossAK8": 20.383,
    "WNAEPt400To500LossAK8": 21.941,
    "WNAEPt500ToInfLossAK8": 16.370,
}

def is_wnae_score_variable(variable: str) -> bool:
    return "WNAE" in variable


def has_wnae_overflow_display(variable: str) -> bool:
    return is_wnae_score_variable(variable) and not VARIABLES[variable].get("wnae_combined_loss")


def jet_variable_label(title: str, jet_label: str) -> str:
    if not jet_label:
        return title.replace(" (J_{{idx}})", "").replace("(J_{{idx}})", "").replace("  ", " ").strip()
    return title.replace("(J_{{idx}})", f"({jet_label})").replace(
        "score(J_{{idx}})", f"score({jet_label})"
    ).replace("{idx}", jet_label)


for jet_index in (1, 2):
    for suffix, (branch, title, nbins, xmin, xmax) in JET_VARIABLES.items():
        variable_config = {
            "branch": branch,
            "jet_index": jet_index - 1,
            "title": title.replace("{idx}", str(jet_index)),
            "x_title": title.replace("{idx}", str(jet_index)),
            "nbins": nbins,
            "xmin": xmin,
            "xmax": xmax,
        }
        if suffix == "PNetScoreAK8":
            variable_config["wp"] = 0.90
            variable_config["wp_side"] = "right"
            variable_config["wp_line_y_fraction"] = 0.900
            variable_config["wp_arrow_y_fraction"] = 0.560
        VARIABLES[f"j{jet_index}{suffix}"] = variable_config

for prefix, jet_label, selection_config in (
    ("j12", "J_{1}+J_{2}", {"leading_good_jets": 2}),
    ("all", "", {"all_good_jets": True}),
):
    for suffix, (branch, title, nbins, xmin, xmax) in JET_VARIABLES.items():
        variable_config = {
            "branch": branch,
            **selection_config,
            "title": jet_variable_label(title, jet_label),
            "x_title": jet_variable_label(title, jet_label),
            "nbins": nbins,
            "xmin": xmin,
            "xmax": xmax,
        }
        if suffix == "PNetScoreAK8":
            variable_config["wp"] = 0.90
            variable_config["wp_side"] = "right"
            variable_config["wp_line_y_fraction"] = 0.900
            variable_config["wp_arrow_y_fraction"] = 0.560
        VARIABLES[f"{prefix}{suffix}"] = variable_config

for jet_index in (1, 2):
    for suffix, (branch, title, nbins, xmin, xmax) in WNAE_PT_LOSS_VARIABLES.items():
        variable_config = {
            "branch": branch,
            "jet_index": jet_index - 1,
            "title": title.replace("{idx}", str(jet_index)),
            "x_title": title.replace("{idx}", str(jet_index)),
            "nbins": nbins,
            "xmin": xmin,
            "xmax": xmax,
            "display_xmin": 5.0,
            "display_xmax": 60.0,
            "display_bin_width": 5.0,
        }
        if suffix in WNAE_MC_WORKING_POINTS:
            variable_config["wp"] = WNAE_MC_WORKING_POINTS[suffix]
            variable_config["wp_side"] = "right"
            variable_config["wp_line_y_fraction"] = 0.620
            variable_config["wp_arrow_y_fraction"] = 0.580
        VARIABLES[f"j{jet_index}{suffix}"] = variable_config
    VARIABLES[f"j{jet_index}WNAECombinedLossAK8"] = {
        "wnae_combined_loss": True,
        "pt_branch": "JetsAK8_/.fPt",
        "jet_index": jet_index - 1,
        "title": f"p_{{T}}-selected WNAE score (J_{{{jet_index}}})",
        "x_title": f"p_{{T}}-selected WNAE score (J_{{{jet_index}}})",
        "nbins": 80,
        "xmin": 0.0,
        "xmax": 200.0,
        "display_xmin": 5.0,
        "display_xmax": 60.0,
        "display_bin_width": 5.0,
    }

for prefix, jet_label, selection_config in (
    ("j12", "J_{1}+J_{2}", {"leading_good_jets": 2}),
    ("all", "", {"all_good_jets": True}),
):
    for suffix, (branch, title, nbins, xmin, xmax) in WNAE_PT_LOSS_VARIABLES.items():
        variable_config = {
            "branch": branch,
            **selection_config,
            "title": jet_variable_label(title, jet_label),
            "x_title": jet_variable_label(title, jet_label),
            "nbins": nbins,
            "xmin": xmin,
            "xmax": xmax,
            "display_xmin": 5.0,
            "display_xmax": 60.0,
            "display_bin_width": 5.0,
        }
        if suffix in WNAE_MC_WORKING_POINTS:
            variable_config["wp"] = WNAE_MC_WORKING_POINTS[suffix]
            variable_config["wp_side"] = "right"
            variable_config["wp_line_y_fraction"] = 0.620
            variable_config["wp_arrow_y_fraction"] = 0.580
        VARIABLES[f"{prefix}{suffix}"] = variable_config
    VARIABLES[f"{prefix}WNAECombinedLossAK8"] = {
        "wnae_combined_loss": True,
        "pt_branch": "JetsAK8_/.fPt",
        **selection_config,
        "title": f"p_{{T}}-selected WNAE score{f' ({jet_label})' if jet_label else ''}",
        "x_title": f"p_{{T}}-selected WNAE score{f' ({jet_label})' if jet_label else ''}",
        "nbins": 80,
        "xmin": 0.0,
        "xmax": 200.0,
        "display_xmin": 5.0,
        "display_xmax": 60.0,
        "display_bin_width": 5.0,
    }

for jet_index in (1, 2):
    for suffix, numerator, denominator, title in (
        ("Tau21AK8", "JetsAK8_NsubjettinessTau2", "JetsAK8_NsubjettinessTau1", "#tau_{21}(J_{{idx}})"),
        ("Tau32AK8", "JetsAK8_NsubjettinessTau3", "JetsAK8_NsubjettinessTau2", "#tau_{32}(J_{{idx}})"),
        ("Tau43AK8", "JetsAK8_NsubjettinessTau4", "JetsAK8_NsubjettinessTau3", "#tau_{43}(J_{{idx}})"),
    ):
        VARIABLES[f"j{jet_index}{suffix}"] = {
            "numerator": numerator,
            "denominator": denominator,
            "jet_index": jet_index - 1,
            "title": title.replace("{idx}", str(jet_index)),
            "x_title": title.replace("{idx}", str(jet_index)),
            "nbins": 40,
            "xmin": 0.0,
            "xmax": 1.0,
        }

for prefix, jet_label, selection_config in (
    ("j12", "J_{1}+J_{2}", {"leading_good_jets": 2}),
    ("all", "", {"all_good_jets": True}),
):
    for suffix, numerator, denominator, title in (
        ("Tau21AK8", "JetsAK8_NsubjettinessTau2", "JetsAK8_NsubjettinessTau1", "#tau_{21}(J_{{idx}})"),
        ("Tau32AK8", "JetsAK8_NsubjettinessTau3", "JetsAK8_NsubjettinessTau2", "#tau_{32}(J_{{idx}})"),
        ("Tau43AK8", "JetsAK8_NsubjettinessTau4", "JetsAK8_NsubjettinessTau3", "#tau_{43}(J_{{idx}})"),
    ):
        VARIABLES[f"{prefix}{suffix}"] = {
            "numerator": numerator,
            "denominator": denominator,
            **selection_config,
            "title": jet_variable_label(title, jet_label),
            "x_title": jet_variable_label(title, jet_label),
            "nbins": 40,
            "xmin": 0.0,
            "xmax": 1.0,
        }


# The original plotter calls rebinCalc(totalBin, 40): 500 -> factor 10 -> 50 bins.
PLOT_TARGET_BINS = 40
DISPLAY_BIN_REDUCTION = 2
JET_AND_WNAE_DISPLAY_BIN_REDUCTION = 2

# These are the five signal benchmarks from the original plotStack2 usage
# translated to the real directory naming convention on EOS.
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

SIGNAL_SCALE_FACTORS = {
    ("600", "0.3"): 10.0,
    ("2000", "0.1"): 500.0,
    ("2000", "0.3"): 500.0,
    ("2000", "0.7"): 1000.0,
    ("4000", "0.3"): 1000.0,
}

# Explicitly requested background directory manifest.
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

# Retained from the supplied plotting script.
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


# =============================================================================
# Histogram accumulation
# =============================================================================
@dataclass
class HistogramAccumulator:
    """Numerically stable sumw/sumw2 storage for one fixed-binning TH1D."""

    edges: np.ndarray
    sumw: np.ndarray = field(init=False)
    sumw2: np.ndarray = field(init=False)
    entries_total: int = 0
    entries_in_range: int = 0

    def __post_init__(self) -> None:
        self.sumw = np.zeros(len(self.edges) - 1, dtype=np.float64)
        self.sumw2 = np.zeros(len(self.edges) - 1, dtype=np.float64)

    @classmethod
    def for_variable(cls, variable: str) -> "HistogramAccumulator":
        cfg = VARIABLES[variable]
        edges = np.linspace(
            cfg["xmin"], cfg["xmax"], cfg["nbins"] + 1, dtype=np.float64
        )
        return cls(edges=edges)

    def fill(self, values: np.ndarray, weights: np.ndarray) -> None:
        values_ak = ak.Array(values)
        weights = np.asarray(weights, dtype=np.float64)
        try:
            counts = np.asarray(ak.to_numpy(ak.num(values_ak, axis=1)), dtype=np.int64)
            weights = np.repeat(weights, counts)
            values_ak = ak.flatten(ak.fill_none(values_ak, np.nan), axis=None)
        except (ValueError, np.AxisError):
            values_ak = ak.fill_none(values_ak, np.nan)

        values = np.asarray(
            ak.to_numpy(values_ak),
            dtype=np.float64,
        )
        self.entries_total += int(values.size)

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
        self.entries_in_range += int(values.size)
        self.sumw += np.histogram(values, bins=self.edges, weights=weights)[0]
        self.sumw2 += np.histogram(
            values, bins=self.edges, weights=weights * weights
        )[0]

    def add(self, other: "HistogramAccumulator") -> None:
        if len(self.edges) != len(other.edges) or not np.allclose(
            self.edges, other.edges
        ):
            raise ValueError("Attempted to add incompatible histograms.")
        self.sumw += other.sumw
        self.sumw2 += other.sumw2
        self.entries_total += other.entries_total
        self.entries_in_range += other.entries_in_range

    def integral(self) -> float:
        return float(np.sum(self.sumw))


def empty_accumulators() -> Dict[str, HistogramAccumulator]:
    return {variable: HistogramAccumulator.for_variable(variable) for variable in VARIABLES}


def add_accumulator_dict(
    target: Dict[str, HistogramAccumulator],
    source: Dict[str, HistogramAccumulator],
) -> None:
    for variable in VARIABLES:
        target[variable].add(source[variable])


def replace_accumulator_variables(
    target: Dict[str, HistogramAccumulator],
    source: Dict[str, HistogramAccumulator],
    variables: Sequence[str],
) -> None:
    for variable in variables:
        target[variable] = source[variable]


def wnae_score_variables() -> List[str]:
    loss_branches = {branch for branch, _, _ in WNAE_LOSS_BRANCHES}
    selected = []
    for variable, cfg in VARIABLES.items():
        if cfg.get("wnae_combined_loss"):
            selected.append(variable)
            continue
        if any(branch in loss_branches for branch in variable_sources(variable)):
            selected.append(variable)
    return selected


# =============================================================================
# Legacy binning helpers
# =============================================================================
def divisor_generator(number: int) -> Iterable[int]:
    large_divisors: List[int] = []
    for divisor in range(1, int(math.sqrt(number)) + 1):
        if number % divisor == 0:
            yield divisor
            if divisor * divisor != number:
                large_divisors.append(number // divisor)
    yield from reversed(large_divisors)


def rebin_calc(nbins: int, target_bins: int) -> int:
    """Identical selection logic to the legacy plotStack2 rebinCalc function."""
    if nbins <= 0:
        return 1
    desired = nbins / float(target_bins)
    divisors = list(divisor_generator(nbins))
    if not divisors:
        return 1
    return max(1, int(min(divisors, key=lambda item: abs(item - desired))))


def plot_rebin_factor(variable: str) -> int:
    nbins = int(VARIABLES[variable]["nbins"])
    factor = rebin_calc(nbins, PLOT_TARGET_BINS)
    visible_bins = nbins // factor if factor > 0 else nbins
    if visible_bins > 10 and nbins % (factor * DISPLAY_BIN_REDUCTION) == 0:
        factor *= DISPLAY_BIN_REDUCTION
    visible_bins = nbins // factor if factor > 0 else nbins
    if (
        visible_bins > 10
        and (variable.startswith(("j1", "j2", "j12", "all")) or "WNAE" in variable)
        and nbins % (factor * JET_AND_WNAE_DISPLAY_BIN_REDUCTION) == 0
    ):
        factor *= JET_AND_WNAE_DISPLAY_BIN_REDUCTION
    return factor


# =============================================================================
# EOS helpers and tree streaming
# =============================================================================
def normalize_eos_path(path: str) -> str:
    """Allow /store/... and mounted /eos/uscms/store/... forms."""
    path = os.path.expanduser(path).rstrip("/")
    mounted_prefix = "/eos/uscms"
    if path.startswith(mounted_prefix + "/"):
        path = path[len(mounted_prefix) :]
    if not path.startswith("/"):
        raise ValueError(f"Expected an absolute EOS path, received: {path}")
    return path


def xrdfs_ls(path: str, recursive: bool = False) -> List[str]:
    command = ["xrdfs", EOS_HOST, "ls"]
    if recursive:
        command.append("-R")
    command.append(normalize_eos_path(path))

    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"xrdfs listing failed for {path}\n"
            f"command: {' '.join(command)}\n"
            f"stderr: {result.stderr.strip()}"
        )
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def eos_url(path: str) -> str:
    return EOS_HOST + "//" + normalize_eos_path(path).lstrip("/")


def get_part_files(sample_dir: str, max_files: int) -> List[str]:
    paths = xrdfs_ls(sample_dir, recursive=True)
    files = sorted(
        path for path in paths if re.search(r"/part-[0-9]+\.root$", path)
    )
    return files[:max_files] if max_files > 0 else files


def variable_sources(variable: str) -> List[str]:
    cfg = VARIABLES[variable]
    if "numerator" in cfg:
        return [cfg["numerator"], cfg["denominator"]]
    if cfg.get("wnae_combined_loss"):
        return [cfg["pt_branch"]] + [branch for branch, _, _ in WNAE_LOSS_BRANCHES]
    return [cfg["branch"]]


def available_weight_branches(tree: uproot.behaviors.TTree.TTree, is_data: bool = False) -> List[str]:
    available = set(tree.keys())
    required = [] if is_data else ["Weight"]
    missing = [branch for branch in required if branch not in available]
    if missing:
        raise KeyError(f"Events tree lacks required branch(es): {missing}")

    for branch in ("puWeight", "NonPrefiringProb", "lundWeightNom"):
        if not is_data and branch in available:
            required.append(branch)
    for variable in VARIABLES:
        for branch in variable_sources(variable):
            if branch in available and branch not in required:
                required.append(branch)
    if "JetsAK8_isGood" in available:
        required.append("JetsAK8_isGood")
    return required


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


def can_compute_variable(variable: str, arrays) -> bool:
    return all(branch in arrays.fields for branch in variable_sources(variable))


def good_jets(jagged, arrays):
    if "JetsAK8_isGood" in arrays.fields:
        jagged = jagged[arrays["JetsAK8_isGood"]]
    return jagged


def first_good_jets(jagged, arrays, index: int):
    jagged = good_jets(jagged, arrays)
    return ak.fill_none(ak.pad_none(jagged, index + 1, axis=1)[:, index], np.nan)


def selected_good_jets(jagged, arrays, cfg):
    jagged = good_jets(jagged, arrays)
    if cfg.get("all_good_jets"):
        return jagged
    if "leading_good_jets" in cfg:
        return ak.pad_none(
            jagged,
            int(cfg["leading_good_jets"]),
            axis=1,
            clip=True,
        )
    return ak.fill_none(
        ak.pad_none(jagged, int(cfg["jet_index"]) + 1, axis=1)[:, int(cfg["jet_index"])],
        np.nan,
    )


def variable_values(variable: str, arrays):
    cfg = VARIABLES[variable]
    if "numerator" in cfg:
        numerator = selected_good_jets(arrays[cfg["numerator"]], arrays, cfg)
        denominator = selected_good_jets(arrays[cfg["denominator"]], arrays, cfg)
        safe_denominator = ak.where(denominator != 0.0, denominator, np.nan)
        with np.errstate(divide="ignore", invalid="ignore"):
            return numerator / safe_denominator
    if cfg.get("wnae_combined_loss"):
        pt = selected_good_jets(arrays[cfg["pt_branch"]], arrays, cfg)
        combined = ak.zeros_like(pt)
        for branch, pt_min, pt_max in WNAE_LOSS_BRANCHES:
            loss = selected_good_jets(arrays[branch], arrays, cfg)
            if pt_max is None:
                in_bin = pt >= pt_min
            else:
                in_bin = (pt >= pt_min) & (pt < pt_max)
            combined = combined + ak.where(in_bin, loss, 0.0)
        return ak.where(np.isfinite(pt), combined, np.nan)
    if "jet_index" in cfg or cfg.get("all_good_jets") or "leading_good_jets" in cfg:
        return selected_good_jets(arrays[cfg["branch"]], arrays, cfg)
    return arrays[cfg["branch"]]


def stream_file(
    file_url: str,
    year: str,
    step_size: str,
    accumulators: Dict[str, HistogramAccumulator],
    is_data: bool = False,
) -> int:
    """Stream nominal MET/ST and apply the standard MC event normalization."""
    with uproot.open(file_url) as root_file:
        if "Events" not in root_file:
            raise KeyError("No Events tree found.")

        tree = root_file["Events"]
        requested = available_weight_branches(tree, is_data=is_data)
        lund_norm = lund_nominal_norm(root_file) if not is_data and "lundWeightNom" in requested else None
        entries = int(tree.num_entries)

        for arrays in tree.iterate(expressions=requested, step_size=step_size, library="ak"):
            if is_data:
                weights = np.ones(len(arrays), dtype=np.float64)
            else:
                weights = ak.to_numpy(arrays["Weight"]) * LUMI_PB[year]
                if "puWeight" in arrays.fields:
                    weights *= ak.to_numpy(arrays["puWeight"])
                if "NonPrefiringProb" in arrays.fields:
                    weights *= ak.to_numpy(arrays["NonPrefiringProb"])
                if lund_norm is not None and "lundWeightNom" in arrays.fields:
                    weights *= ak.to_numpy(arrays["lundWeightNom"]) * lund_norm

            for variable in VARIABLES:
                if can_compute_variable(variable, arrays):
                    accumulators[variable].fill(variable_values(variable, arrays), weights)
    return entries


# =============================================================================
# Signal discovery and labels
# =============================================================================
def parse_csv(value: str) -> List[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def ordered_signal_names(names: Sequence[str]) -> List[str]:
    """Keep the legacy benchmark order first, then alphabetical extras."""
    default_index = {name: index for index, name in enumerate(DEFAULT_SIGNAL_DIRS)}
    return sorted(names, key=lambda name: (default_index.get(name, 10_000), name))


def is_disabled_signal_name(signal_name: str) -> bool:
    lowered = signal_name.lower()
    return any(token.lower() in lowered for token in DISABLED_SIGNAL_TOKENS)


def signal_matches_request(signal_name: str, token: str) -> bool:
    """Full t-channel directory tokens are exact; short tokens remain substrings."""
    name_lower = signal_name.lower()
    token_lower = token.lower()
    if token_lower.startswith("t-channel_"):
        return token_lower == name_lower
    return token_lower in name_lower


def get_signal_directories(
    signal_base: str,
    year: str,
    requested: Sequence[str],
    strict: bool,
) -> List[Tuple[str, str]]:
    """
    Find requested direct signal subdirectories under either:
      signal_base/YEAR/t_channel_pre_selection/nominal/
      signal_base/YEAR/nominal/

    A requested token can be a full directory name or any unique substring.
    The literal special value 'all' selects every t-channel_* directory.
    """
    signal_base = normalize_eos_path(signal_base)
    nominal_candidates = [
        f"{signal_base}/{year}/t_channel_pre_selection/nominal",
        f"{signal_base}/{year}/nominal",
    ]
    entries = None
    nominal = None
    errors = []
    for candidate in nominal_candidates:
        try:
            entries = xrdfs_ls(candidate, recursive=False)
            nominal = candidate
            break
        except RuntimeError as exc:
            errors.append(str(exc))

    if entries is None or nominal is None:
        message = (
            f"{year}: could not list signal directory under any expected layout: "
            + ", ".join(nominal_candidates)
        )
        if strict:
            raise RuntimeError(message + "\n" + "\n".join(errors))
        print(f"[WARN] {message}")
        return []

    available = {
        os.path.basename(entry.rstrip("/")): entry
        for entry in entries
        if os.path.basename(entry.rstrip("/")).startswith("t-channel_")
    }

    available = {
        name: path
        for name, path in available.items()
        if not is_disabled_signal_name(name)
    }

    if len(requested) == 1 and requested[0].lower() in {"all", "*"}:
        return [(name, available[name]) for name in ordered_signal_names(available)]

    selected: List[Tuple[str, str]] = []
    seen = set()

    for token in requested:
        matches = [
            name for name in available if signal_matches_request(name, token)
        ]

        if not matches:
            message = f"{year}: no t-channel signal matches '{token}'"
            if strict:
                raise RuntimeError(message)
            print(f"[WARN] {message}")
            continue

        # A full directory name selects one; a loose substring can select more.
        for name in ordered_signal_names(matches):
            if name not in seen:
                selected.append((name, available[name]))
                seen.add(name)

    return selected


def signal_legend_label(signal_name: str) -> str:
    """Use the m_phi / r_inv naming convention from the legacy plotting script."""
    mmed = re.search(r"mMed-([0-9]+(?:p[0-9]+)?)", signal_name, re.IGNORECASE)
    rinv = re.search(r"rinv-([0-9]+(?:p[0-9]+)?)", signal_name, re.IGNORECASE)

    if mmed and rinv:
        return (
            f"m_{{#phi}} = {mmed.group(1).replace('p', '.')} GeV, "
            f"r_{{inv}} = {rinv.group(1).replace('p', '.')}"
        )
    return signal_name


def signal_parameter_key(signal_name: str) -> Optional[Tuple[str, str]]:
    mmed = re.search(r"mMed-([0-9]+(?:p[0-9]+)?)", signal_name, re.IGNORECASE)
    rinv = re.search(r"rinv-([0-9]+(?:p[0-9]+)?)", signal_name, re.IGNORECASE)
    if not (mmed and rinv):
        return None
    return (mmed.group(1).replace("p", "."), rinv.group(1).replace("p", "."))


def signal_scale_factor(
    signal_name: str,
    scale_factors: Optional[Dict[Tuple[str, str], float]] = None,
) -> float:
    key = signal_parameter_key(signal_name)
    if key is None:
        return 1.0
    factors = SIGNAL_SCALE_FACTORS if scale_factors is None else scale_factors
    return factors.get(key, 1.0)


def scaled_signal_legend_label(
    signal_name: str,
    normalized: bool,
    scale_factors: Optional[Dict[Tuple[str, str], float]] = None,
) -> str:
    label = signal_legend_label(signal_name)
    scale = signal_scale_factor(signal_name, scale_factors)
    if normalized or abs(scale - 1.0) < 1.0e-9:
        return label
    return f"{label} (#times {scale:g})"


def unique_legend_entries(entries):
    unique = []
    seen = set()
    for hist, label, option in entries:
        if label in seen:
            continue
        unique.append((hist, label, option))
        seen.add(label)
    return unique


# =============================================================================
# ROOT I/O
# =============================================================================
def mkdirs(root_file: ROOT.TFile, directory_path: str) -> ROOT.TDirectory:
    directory: ROOT.TDirectory = root_file
    for part in (item for item in directory_path.split("/") if item):
        next_dir = directory.GetDirectory(part)
        if not next_dir:
            next_dir = directory.mkdir(part)
        directory = next_dir
    return directory


def accumulator_to_histogram(
    variable: str,
    accumulator: HistogramAccumulator,
    name: Optional[str] = None,
    y_title: str = "Events",
) -> ROOT.TH1D:
    cfg = VARIABLES[variable]
    hist_name = name or variable
    edges = array("d", accumulator.edges.tolist())
    hist = ROOT.TH1D(
        hist_name,
        f"{hist_name};{cfg['x_title']};{y_title}",
        len(edges) - 1,
        edges,
    )
    hist.Sumw2()

    for ibin, (content, variance) in enumerate(
        zip(accumulator.sumw, accumulator.sumw2), start=1
    ):
        hist.SetBinContent(ibin, float(content))
        hist.SetBinError(ibin, math.sqrt(max(0.0, float(variance))))

    hist.SetDirectory(0)
    return hist


def write_histogram(
    root_file: ROOT.TFile,
    directory_path: str,
    histogram: ROOT.TH1,
    key_name: str,
) -> None:
    directory = mkdirs(root_file, directory_path)
    directory.cd()
    histogram.Write(key_name, ROOT.TObject.kOverwrite)


def write_cache(
    output_path: str,
    sample_hists: Dict[Tuple[str, str, str], Dict[str, HistogramAccumulator]],
    group_hists: Dict[Tuple[str, str], Dict[str, HistogramAccumulator]],
    run2_hists: Dict[str, Dict[str, HistogramAccumulator]],
    signal_hists: Dict[Tuple[str, str], Dict[str, HistogramAccumulator]],
    run2_signal_hists: Dict[str, Dict[str, HistogramAccumulator]],
    data_hists: Dict[Tuple[str, str], Dict[str, HistogramAccumulator]],
    run2_data_hists: Optional[Dict[str, HistogramAccumulator]],
    years: Sequence[str],
) -> None:
    """Write raw and normalized MET/ST TH1D objects, and no event-level trees."""
    parent = os.path.dirname(os.path.abspath(output_path))
    os.makedirs(parent, exist_ok=True)

    output = ROOT.TFile.Open(output_path, "RECREATE")
    if not output or output.IsZombie():
        raise RuntimeError(f"Could not create ROOT output: {output_path}")

    try:
        # Raw per-sample histograms.
        for (year, process, sample), accumulators in sorted(sample_hists.items()):
            for variable in VARIABLES:
                write_histogram(
                    output,
                    f"raw/samples/{year}/{process}/{sample}",
                    accumulator_to_histogram(variable, accumulators[variable]),
                    variable,
                )

        # Raw per-year process histograms.
        for (year, process), accumulators in sorted(group_hists.items()):
            for variable in VARIABLES:
                write_histogram(
                    output,
                    f"raw/groups/{year}/{process}",
                    accumulator_to_histogram(variable, accumulators[variable]),
                    variable,
                )

        # Raw Run-2 process histograms.
        for process, accumulators in sorted(run2_hists.items()):
            for variable in VARIABLES:
                write_histogram(
                    output,
                    f"raw/groups/Run2/{process}",
                    accumulator_to_histogram(variable, accumulators[variable]),
                    variable,
                )

        # Raw signal histograms.
        for (year, signal_name), accumulators in sorted(signal_hists.items()):
            for variable in VARIABLES:
                write_histogram(
                    output,
                    f"raw/signals/{year}/{signal_name}",
                    accumulator_to_histogram(variable, accumulators[variable]),
                    variable,
                )

        for signal_name, accumulators in sorted(run2_signal_hists.items()):
            for variable in VARIABLES:
                write_histogram(
                    output,
                    f"raw/signals/Run2/{signal_name}",
                    accumulator_to_histogram(variable, accumulators[variable]),
                    variable,
                )

        for (year, sample), accumulators in sorted(data_hists.items()):
            for variable in VARIABLES:
                write_histogram(
                    output,
                    f"raw/data/{year}/{sample}",
                    accumulator_to_histogram(variable, accumulators[variable]),
                    variable,
                )

        for year in years:
            year_data = empty_accumulators()
            found = False
            for (data_year, _), accumulators in data_hists.items():
                if data_year == year:
                    add_accumulator_dict(year_data, accumulators)
                    found = True
            if found:
                for variable in VARIABLES:
                    write_histogram(
                        output,
                        f"raw/data/{year}/observed",
                        accumulator_to_histogram(variable, year_data[variable]),
                        variable,
                    )

        if run2_data_hists:
            for variable in VARIABLES:
                write_histogram(
                    output,
                    "raw/data/Run2/observed",
                    accumulator_to_histogram(variable, run2_data_hists[variable]),
                    variable,
                )

        # Normalized shape cache for every year plus combined Run-2.
        era_backgrounds: Dict[str, Dict[str, Dict[str, HistogramAccumulator]]] = {}
        for year in years:
            era_backgrounds[year] = {
                process: group_hists[(year, process)]
                for process in PROCESS_ORDER
                if (year, process) in group_hists
            }
        era_backgrounds["Run2"] = run2_hists

        era_signals: Dict[str, Dict[str, Dict[str, HistogramAccumulator]]] = {}
        for year in years:
            era_signals[year] = {
                signal_name: accumulators
                for (signal_year, signal_name), accumulators in signal_hists.items()
                if signal_year == year
            }
        era_signals["Run2"] = run2_signal_hists

        for era, processes in era_backgrounds.items():
            for variable in VARIABLES:
                total_integral = sum(
                    accumulators[variable].integral()
                    for accumulators in processes.values()
                )
                if total_integral <= 0.0:
                    continue

                for process, accumulators in processes.items():
                    hist = accumulator_to_histogram(
                        variable,
                        accumulators[variable],
                        y_title="Arbitrary units",
                    )
                    hist.Scale(1.0 / total_integral)
                    write_histogram(
                        output,
                        f"shapes/{era}/groups/{process}",
                        hist,
                        variable,
                    )

        for era, signals in era_signals.items():
            for signal_name, accumulators in signals.items():
                for variable in VARIABLES:
                    hist = accumulator_to_histogram(
                        variable,
                        accumulators[variable],
                        y_title="Arbitrary units",
                    )
                    integral = hist.Integral()
                    if integral <= 0.0:
                        continue
                    hist.Scale(1.0 / integral)
                    write_histogram(
                        output,
                        f"shapes/{era}/signals/{signal_name}",
                        hist,
                        variable,
                    )

        for era in list(years) + ["Run2"]:
            source = run2_data_hists if era == "Run2" else None
            if era != "Run2":
                source = empty_accumulators()
                found = False
                for (data_year, _), accumulators in data_hists.items():
                    if data_year == era:
                        add_accumulator_dict(source, accumulators)
                        found = True
                if not found:
                    source = None
            if not source:
                continue
            for variable in VARIABLES:
                hist = accumulator_to_histogram(
                    variable,
                    source[variable],
                    y_title="Arbitrary units",
                )
                integral = hist.Integral()
                if integral <= 0.0:
                    continue
                hist.Scale(1.0 / integral)
                write_histogram(
                    output,
                    f"shapes/{era}/data/observed",
                    hist,
                    variable,
                )

        reference_variable = next(iter(VARIABLES))
        metadata = (
            "Raw weights: Weight*lumi_pb*puWeight*NonPrefiringProb when branches exist; "
            "files with lundWeightNom and CutFlow/InitialLundNominal additionally use "
            "lundWeightNom*Initial/InitialLundNominal. "
            "Data weights: one per event. "
            "Shapes: stack MC normalized by total MC integral; each signal and data normalized independently. "
            f"Stored variables: {', '.join(VARIABLES.keys())}. "
            f"Display rebin reference ({reference_variable}): rebinCalc(500, {PLOT_TARGET_BINS})="
            f"{plot_rebin_factor(reference_variable)}."
        )
        meta_dir = mkdirs(output, "metadata")
        meta_dir.cd()
        ROOT.TNamed("normalization", metadata).Write(
            "normalization", ROOT.TObject.kOverwrite
        )
    finally:
        output.Close()


# =============================================================================
# Plot style and ROOT-cache plotting
# =============================================================================
def clone_histogram(histogram: ROOT.TH1, name: str) -> ROOT.TH1:
    clone = histogram.Clone(name)
    clone.SetDirectory(0)
    return clone


def clone_wnae_overflow_display(histogram: ROOT.TH1, variable: str, name: str) -> ROOT.TH1:
    xmin = float(VARIABLES[variable].get("display_xmin", VARIABLES[variable]["xmin"]))
    xmax = float(VARIABLES[variable].get("display_xmax", 60.0))
    bin_width = float(VARIABLES[variable].get("display_bin_width", 2.5))
    wp = VARIABLES[variable].get("wp")
    if wp is not None and xmin < float(wp) < xmax:
        edge_start = float(wp) - math.ceil((float(wp) - xmin) / bin_width) * bin_width
    else:
        edge_start = xmin
    edge_stop = edge_start + math.ceil((xmax - edge_start) / bin_width) * bin_width
    edges = np.arange(edge_start, edge_stop + 0.5 * bin_width, bin_width)
    if wp is not None and xmin < float(wp) < xmax:
        edges[np.argmin(np.abs(edges - float(wp)))] = float(wp)
    edges = np.asarray([round(float(edge), 6) for edge in edges], dtype=np.float64)
    rebinned = ROOT.TH1D(f"{name}_overflow", histogram.GetTitle(), len(edges) - 1, array("d", edges))
    rebinned.SetDirectory(0)
    rebinned.Sumw2()

    last_bin = rebinned.GetNbinsX()
    for ibin in range(1, histogram.GetNbinsX() + 1):
        content = histogram.GetBinContent(ibin)
        error = histogram.GetBinError(ibin)
        if content == 0.0 and error == 0.0:
            continue
        center = histogram.GetBinCenter(ibin)
        if center < xmin:
            target_bin = 1
        elif center >= xmax:
            target_bin = last_bin
        else:
            target_bin = rebinned.FindBin(center)
        if target_bin > last_bin:
            target_bin = last_bin
        rebinned.SetBinContent(target_bin, rebinned.GetBinContent(target_bin) + content)
        rebinned.SetBinError(target_bin, math.hypot(rebinned.GetBinError(target_bin), error))
    return rebinned


def clone_rebin(histogram: ROOT.TH1, variable: str, name: str) -> ROOT.TH1:
    if has_wnae_overflow_display(variable):
        return clone_wnae_overflow_display(histogram, variable, name)
    clone = clone_histogram(histogram, name)
    configured_bins = int(VARIABLES[variable]["nbins"])
    current_bins = int(clone.GetNbinsX())
    factor = plot_rebin_factor(variable)
    if current_bins != configured_bins and current_bins % configured_bins == 0:
        factor *= current_bins // configured_bins
    display_rebin_factor = int(VARIABLES[variable].get("display_rebin_factor", 1))
    if display_rebin_factor > 1 and current_bins % (factor * display_rebin_factor) == 0:
        factor *= display_rebin_factor
    if current_bins % factor != 0:
        factor = 1
    if factor > 1:
        rebinned = clone.Rebin(factor, f"{name}_rebin{factor}")
        rebinned.SetDirectory(0)
        return rebinned
    return clone


def style_background(histogram: ROOT.TH1, process: str) -> None:
    color = PROC_COLOR[process]
    histogram.SetFillColor(color)
    histogram.SetFillStyle(1001)
    histogram.SetLineColor(color)
    histogram.SetLineWidth(0)
    histogram.SetMarkerSize(0)


def style_signal(histogram: ROOT.TH1, index: int) -> None:
    histogram.SetLineColor(SIGNAL_LINE_COLORS[index % len(SIGNAL_LINE_COLORS)])
    histogram.SetLineWidth(2)
    histogram.SetLineStyle(SIGNAL_LINE_STYLES[index % len(SIGNAL_LINE_STYLES)])
    histogram.SetFillStyle(0)
    histogram.SetMarkerSize(0)


def style_data(histogram: ROOT.TH1) -> None:
    histogram.SetLineColor(ROOT.kBlack)
    histogram.SetMarkerColor(ROOT.kBlack)
    histogram.SetMarkerStyle(ROOT.kFullCircle)
    histogram.SetMarkerSize(0.85)
    histogram.SetLineWidth(2)
    histogram.SetFillStyle(0)


def style_uncertainty_band(histogram: ROOT.TH1) -> None:
    histogram.SetFillColor(ROOT.kGray + 1)
    histogram.SetFillStyle(3354)
    histogram.SetLineColor(ROOT.kGray + 2)
    histogram.SetLineWidth(1)
    histogram.SetMarkerSize(0)


def clone_uncertainty_band(histogram: ROOT.TH1, name: str) -> ROOT.TH1:
    band = clone_histogram(histogram, name)
    style_uncertainty_band(band)
    return band


def lumi_label(era: str) -> Optional[str]:
    if era == "Run2":
        return "138 fb^{-1} (13 TeV)"
    elif era in LUMI_PB:
        lumi_fb = LUMI_PB[era] / 1000.0
    else:
        return None
    return f"{lumi_fb:.1f} fb^{{-1}} (13 TeV)"


def draw_cms_label(
    simulation: bool = True,
    era: Optional[str] = None,
    show_lumi: bool = False,
    preliminary: bool = False,
) -> None:
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
        latex.SetTextAlign(11)
        latex.DrawLatex(0.255, label_y, "Simulation")
        if preliminary:
            latex.DrawLatex(0.420, label_y, "Preliminary")
    elif preliminary:
        latex.SetTextFont(52)
        latex.SetTextSize(0.040)
        latex.SetTextAlign(11)
        latex.DrawLatex(0.255, label_y, "Preliminary")

    latex.SetTextFont(42)
    latex.SetTextSize(0.045)
    latex.SetTextAlign(31)
    latex.DrawLatex(0.95, label_y, lumi_label(era) if show_lumi and era else "(13 TeV)")


def draw_cms_simulation_label() -> None:
    draw_cms_label(simulation=True)


def draw_common_signal_text(y_position: float = 0.700) -> None:
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.030)
    latex.SetTextColor(ROOT.kBlack)
    latex.SetTextAlign(12)
    latex.DrawLatex(0.565, y_position, "m_{dark} = 20 GeV, #kern[-0.18]{#lambda} = 1")


def draw_ratio_y_title(title: str = "Data/Sim") -> None:
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.12)
    latex.SetTextAngle(90)
    latex.SetTextAlign(31)
    latex.DrawLatex(0.085, 0.985, title)


def interpolate_axis_y(ymin: float, ymax: float, fraction: float, log_y: bool) -> float:
    fraction = max(0.0, min(1.0, fraction))
    if log_y and ymin > 0.0 and ymax > ymin:
        return math.exp(math.log(ymin) + fraction * (math.log(ymax) - math.log(ymin)))
    return ymin + fraction * (ymax - ymin)


def draw_working_point_marker(variable: str, axis: ROOT.TH1, log_y: bool):
    cfg = VARIABLES[variable]
    if "wp" not in cfg:
        return None

    x = float(cfg["wp"])
    xmin = float(cfg.get("display_xmin", cfg["xmin"]))
    xmax = float(cfg.get("display_xmax", cfg["xmax"]))
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


def draw_working_point_arrow(variable: str, axis: ROOT.TH1, log_y: bool):
    cfg = VARIABLES[variable]
    if "wp" not in cfg:
        return None

    x = float(cfg["wp"])
    xmin = float(cfg.get("display_xmin", cfg["xmin"]))
    xmax = float(cfg.get("display_xmax", cfg["xmax"]))
    if x <= xmin or x >= xmax:
        return None

    arrow_y = interpolate_axis_y(
        axis.GetMinimum(),
        axis.GetMaximum(),
        float(cfg.get("wp_arrow_y_fraction", 0.74)),
        log_y,
    )
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
    return arrow


def configure_axis_exponent(variable: str) -> None:
    cfg = VARIABLES[variable]
    ROOT.TGaxis.SetMaxDigits(int(cfg.get("y_max_digits", 3)))
    if cfg.get("y_exponent_left", False):
        x_offset, y_offset = cfg.get("y_exponent_offset", (-0.070, -0.015))
        ROOT.TGaxis.SetExponentOffset(float(x_offset), float(y_offset), "y")
    else:
        ROOT.TGaxis.SetExponentOffset(0.0, 0.0, "y")


def get_root_histogram(
    root_file: ROOT.TFile,
    object_path: str,
    clone_name: str,
) -> Optional[ROOT.TH1]:
    obj = root_file.Get(object_path)
    if not obj:
        return None
    if not obj.InheritsFrom("TH1"):
        return None
    return clone_histogram(obj, clone_name)


def list_signal_shape_names(root_file: ROOT.TFile, era: str) -> List[str]:
    directory = root_file.GetDirectory(f"shapes/{era}/signals")
    if not directory:
        return []
    return ordered_signal_names(
        [key.GetName() for key in directory.GetListOfKeys()]
    )


def selected_cached_signal_names(
    root_file: ROOT.TFile,
    era: str,
    requested_signals: Sequence[str],
) -> List[str]:
    available = [
        name
        for name in list_signal_shape_names(root_file, era)
        if not is_disabled_signal_name(name)
    ]
    if not requested_signals:
        return []
    if len(requested_signals) == 1 and requested_signals[0].lower() in {"all", "*"}:
        return available

    selected: List[str] = []
    seen = set()
    for token in requested_signals:
        matches = [
            name
            for name in available
            if signal_matches_request(name, token)
        ]
        if not matches:
            print(f"[WARN] {era}: cached signal not found for token '{token}'")
            continue
        for name in ordered_signal_names(matches):
            if name not in seen:
                selected.append(name)
                seen.add(name)
    return selected


def draw_stack(
    root_file: ROOT.TFile,
    era: str,
    variable: str,
    output_base: str,
    requested_signals: Sequence[str],
    normalized: bool = True,
    also_linear: bool = False,
    include_data: bool = False,
    include_ratio: bool = False,
    force_log_y: Optional[bool] = None,
    draw_uncertainty_band: bool = True,
    preliminary_label: bool = False,
) -> None:
    """Draw stack(s) directly from cached histograms."""
    directory_kind = "shapes" if normalized else "raw"
    y_title = "Arbitrary units" if normalized else "Events"
    group_prefix = f"shapes/{era}/groups" if normalized else f"raw/groups/{era}"
    signal_prefix = f"shapes/{era}/signals" if normalized else f"raw/signals/{era}"
    data_prefix = f"shapes/{era}/data" if normalized else f"raw/data/{era}"
    process_hists: Dict[str, ROOT.TH1] = {}
    for process in STACK_ORDER:
        hist = get_root_histogram(
            root_file,
            f"{group_prefix}/{process}/{variable}",
            f"{era}_{process}_{variable}",
        )
        if not hist:
            continue
        display = clone_rebin(hist, variable, f"display_{era}_{process}_{variable}")
        if display.Integral() <= 0.0:
            continue
        style_background(display, process)
        process_hists[process] = display

    if not process_hists:
        print(f"[WARN] No cached {directory_kind} {variable} backgrounds for {era}.")
        return

    scale_factors = VARIABLES[variable].get("signal_scale_factors")
    signal_hists: List[Tuple[str, ROOT.TH1]] = []
    for index, signal_name in enumerate(selected_cached_signal_names(root_file, era, requested_signals)):
        if is_disabled_signal_name(signal_name):
            continue
        hist = get_root_histogram(
            root_file,
            f"{signal_prefix}/{signal_name}/{variable}",
            f"{era}_{signal_name}_{variable}",
        )
        if not hist:
            continue
        display = clone_rebin(
            hist,
            variable,
            f"display_{era}_signal_{index}_{variable}",
        )
        if not normalized:
            display.Scale(signal_scale_factor(signal_name, scale_factors))
        if display.Integral() <= 0.0:
            continue
        style_signal(display, index)
        signal_hists.append((signal_name, display))

    data_hist = None
    if include_data:
        raw_data = get_root_histogram(
            root_file,
            f"{data_prefix}/observed/{variable}",
            f"{era}_data_{variable}",
        )
        if raw_data and raw_data.Integral() > 0.0:
            data_hist = clone_rebin(raw_data, variable, f"display_{era}_data_{variable}")
            style_data(data_hist)
        else:
            print(f"[WARN] No cached {directory_kind} {variable} data for {era}.")

    if include_ratio and data_hist is None:
        print(f"[WARN] Cannot draw Data/Sim ratio for {era}/{variable} without cached data.")
        return

    draw_modes = [force_log_y] if force_log_y is not None else ([True, False] if also_linear else [True])
    for log_y in draw_modes:
        scale_tag = "normalized" if normalized else "raw"
        ratio_tag = "_ratio" if include_ratio else ""
        tag = f"{'log' if log_y else 'linear'}_{scale_tag}{ratio_tag}"
        canvas = ROOT.TCanvas(
            f"canvas_{era}_{variable}_{tag}",
            f"canvas_{era}_{variable}_{tag}",
            800,
            900 if include_ratio else 800,
        )
        canvas.cd()

        if include_ratio:
            pad_top = ROOT.TPad(f"pad_top_{era}_{variable}_{tag}", "", 0.0, 0.30, 1.0, 1.0)
            pad_bottom = ROOT.TPad(f"pad_bottom_{era}_{variable}_{tag}", "", 0.0, 0.0, 1.0, 0.30)
            pad_top.SetLeftMargin(0.16)
            pad_top.SetRightMargin(0.05)
            pad_top.SetTopMargin(0.10)
            pad_top.SetBottomMargin(0.02)
            pad_top.SetTicks(1, 1)
            pad_top.SetLogy(log_y)
            pad_bottom.SetLeftMargin(0.16)
            pad_bottom.SetRightMargin(0.05)
            pad_bottom.SetTopMargin(0.04)
            pad_bottom.SetBottomMargin(0.38)
            pad_bottom.SetTicks(1, 1)
            pad_bottom.SetGridy(True)
            pad_top.Draw()
            pad_bottom.Draw()
            pad_top.cd()
        else:
            ROOT.gPad.SetLeftMargin(0.16)
            ROOT.gPad.SetRightMargin(0.05)
            ROOT.gPad.SetTopMargin(0.08)
            ROOT.gPad.SetBottomMargin(0.12)
            ROOT.gPad.SetTicks(1, 1)
            ROOT.gPad.SetLogy(log_y)

        stack = ROOT.THStack(f"stack_{era}_{variable}_{tag}", "")
        total = None
        for process in STACK_ORDER:
            hist = process_hists.get(process)
            if not hist:
                continue
            stack.Add(hist)
            if total is None:
                total = clone_histogram(hist, f"total_{era}_{variable}_{tag}")
            else:
                total.Add(hist)

        if total is None:
            canvas.Close()
            continue

        axis = clone_histogram(total, f"axis_{era}_{variable}_{tag}")
        axis.Reset("ICESM")
        axis.SetStats(0)
        axis.SetTitle("")
        axis.GetXaxis().SetTitle(VARIABLES[variable]["x_title"])
        axis.GetYaxis().SetTitle(y_title)
        axis.GetXaxis().SetTitleSize(0.0 if include_ratio else 0.05)
        axis.GetXaxis().SetLabelSize(0.0 if include_ratio else 0.04)
        axis.GetYaxis().SetTitleSize(0.055 if include_ratio else 0.05)
        axis.GetYaxis().SetLabelSize(0.045 if include_ratio else 0.04)
        axis.GetYaxis().SetLabelOffset(0.01)
        axis.GetYaxis().SetTitleOffset(1.18 if include_ratio else 1.30)
        axis.GetXaxis().SetRangeUser(
            VARIABLES[variable].get("display_xmin", VARIABLES[variable]["xmin"]),
            VARIABLES[variable].get("display_xmax", VARIABLES[variable]["xmax"]),
        )

        ymax = total.GetMaximum()
        for _, hist in signal_hists:
            ymax = max(ymax, hist.GetMaximum())
        if data_hist:
            ymax = max(ymax, data_hist.GetMaximum())
        ymax = max(ymax, 1.0)

        if log_y and not normalized and ("PNetScore" in variable or is_wnae_score_variable(variable)):
            axis.SetMinimum(1.0e2)
        else:
            axis.SetMinimum(1.0 if log_y and not normalized else 1.0e-2)
        if not log_y and "linear_ymax" in VARIABLES[variable]:
            y_max = float(VARIABLES[variable]["linear_ymax"])
        else:
            y_max_factor = 100.0 if log_y else float(VARIABLES[variable].get("linear_ymax_factor", 2.6))
            y_max = ymax * y_max_factor
        if log_y and not normalized:
            y_max = max(y_max, 1.0e9)
        axis.SetMaximum(y_max)

        configure_axis_exponent(variable)
        axis.Draw("hist")
        stack.Draw("hist same")
        total_unc = None
        if draw_uncertainty_band:
            total_unc = clone_uncertainty_band(total, f"total_unc_{era}_{variable}_{tag}")
            total_unc.Draw("E2 same")
        for _, hist in signal_hists:
            hist.Draw("hist same")
        if data_hist:
            data_hist.Draw("E1 same")
        wp_marker = draw_working_point_marker(variable, axis, log_y)

        # Exact MC-only two-column arrangement from the supplied stack plotter:
        # background entries on the left, signal overlays on the right.
        if include_data:
            if is_wnae_score_variable(variable):
                bkg_legend_coords = (0.28, 0.43)
                sig_legend_coords = (0.48, 0.99)
            else:
                bkg_legend_coords = (0.21, 0.41)
                sig_legend_coords = (0.45, 0.98)
            if include_ratio:
                if log_y and not normalized:
                    legend_y1 = 0.605 if is_wnae_score_variable(variable) else 0.600
                else:
                    legend_y1 = 0.615
                legend_text_size = 0.039
                legend_entry_separation = 0.047 if log_y and not normalized and is_wnae_score_variable(variable) else (0.050 if log_y and not normalized else 0.052)
                legend_y2 = 0.865 if log_y and not normalized and is_wnae_score_variable(variable) else (0.855 if log_y and not normalized else 0.875)
            else:
                if log_y and not normalized:
                    legend_y1 = 0.625 if is_wnae_score_variable(variable) else 0.620
                else:
                    legend_y1 = 0.660
                legend_text_size = 0.038
                legend_entry_separation = 0.047 if log_y and not normalized and is_wnae_score_variable(variable) else (0.050 if log_y and not normalized else 0.050)
                legend_y2 = 0.865 if log_y and not normalized and is_wnae_score_variable(variable) else (0.855 if log_y and not normalized else 0.910)
        else:
            bkg_legend_coords = (0.21, 0.41)
            sig_legend_coords = (0.45, 0.98)
            legend_y1 = 0.740
            legend_text_size = 0.041
            legend_entry_separation = 0.034
            legend_y2 = 0.900

        background_entries = [
            (process_hists[process], PROC_LABEL[process], "F")
            for process in LEGEND_ORDER
            if process in process_hists
        ]
        if data_hist:
            background_entries.insert(0, (data_hist, "Data", "EP"))
        if total_unc:
            background_entries.append((total_unc, "Sim. unc.", "F"))
        signal_entries = unique_legend_entries(
            [
                (hist, scaled_signal_legend_label(signal_name, normalized, scale_factors), "L")
                for signal_name, hist in signal_hists
            ]
        )

        if signal_hists:
            signal_entries.append(("", "m_{dark} = 20 GeV, #kern[-0.18]{#lambda} = 1", ""))

        def configure_legend(legend, margin: float = 0.25):
            legend.SetNColumns(1)
            legend.SetBorderSize(0)
            legend.SetFillColor(ROOT.kWhite)
            legend.SetFillStyle(1001)
            legend.SetTextSize(legend_text_size)
            legend.SetTextFont(42)
            legend.SetEntrySeparation(legend_entry_separation)
            legend.SetMargin(margin)

        if signal_hists:
            bkg_legend = ROOT.TLegend(bkg_legend_coords[0], legend_y1, bkg_legend_coords[1], legend_y2)
            sig_legend = ROOT.TLegend(sig_legend_coords[0], legend_y1, sig_legend_coords[1], legend_y2)
            configure_legend(bkg_legend)
            bkg_width = bkg_legend_coords[1] - bkg_legend_coords[0]
            sig_width = sig_legend_coords[1] - sig_legend_coords[0]
            configure_legend(sig_legend, margin=0.25 * bkg_width / sig_width if sig_width > 0.0 else 0.25)
            for entry in background_entries:
                bkg_legend.AddEntry(*entry)
            for entry in signal_entries:
                sig_legend.AddEntry(*entry)
            bkg_legend.Draw()
            sig_legend.Draw()
            if wp_marker and is_wnae_score_variable(variable):
                draw_working_point_arrow(variable, axis, log_y)
        else:
            legend = ROOT.TLegend(0.21 if include_data else 0.23, legend_y1, 0.98 if include_data else 0.97, legend_y2)
            configure_legend(legend)
            for entry in background_entries:
                legend.AddEntry(*entry)
            legend.Draw()
            if wp_marker and is_wnae_score_variable(variable):
                draw_working_point_arrow(variable, axis, log_y)

        draw_cms_label(
            simulation=not include_data,
            era=era,
            show_lumi=True,
            preliminary=preliminary_label,
        )
        ROOT.gPad.RedrawAxis()

        if include_ratio:
            pad_bottom.cd()
            ratio = clone_histogram(data_hist, f"ratio_{era}_{variable}_{tag}")
            ratio.Divide(total)
            ratio.SetStats(0)
            ratio.SetTitle("")
            ratio.SetMinimum(0.0)
            ratio.SetMaximum(2.0)
            ratio.GetXaxis().SetTitle(VARIABLES[variable]["x_title"])
            ratio.GetYaxis().SetTitle("")
            ratio.GetYaxis().CenterTitle(False)
            ratio.GetYaxis().SetNdivisions(505)
            ratio.GetXaxis().SetTitleSize(0.13)
            ratio.GetXaxis().SetLabelSize(0.11)
            ratio.GetXaxis().SetTitleOffset(1.05)
            ratio.GetYaxis().SetLabelSize(0.10)
            ratio.GetYaxis().SetTitleOffset(0.40)
            ratio.GetYaxis().SetTitleSize(0.12)
            ratio.GetYaxis().SetLabelOffset(0.01)
            ratio.GetXaxis().SetRangeUser(
                VARIABLES[variable].get("display_xmin", VARIABLES[variable]["xmin"]),
                VARIABLES[variable].get("display_xmax", VARIABLES[variable]["xmax"]),
            )
            ratio.Draw("PE")

            if draw_uncertainty_band:
                ratio_unc = clone_histogram(total, f"ratio_unc_{era}_{variable}_{tag}")
                for ibin in range(1, ratio_unc.GetNbinsX() + 1):
                    mc = total.GetBinContent(ibin)
                    err = total.GetBinError(ibin)
                    ratio_unc.SetBinContent(ibin, 1.0 if mc > 0.0 else 0.0)
                    ratio_unc.SetBinError(ibin, err / mc if mc > 0.0 else 0.0)
                style_uncertainty_band(ratio_unc)
                ratio_unc.Draw("E2 same")
                ratio.Draw("PE same")

            line_xmin = VARIABLES[variable].get("display_xmin", VARIABLES[variable]["xmin"])
            line_xmax = VARIABLES[variable].get("display_xmax", VARIABLES[variable]["xmax"])
            unity = ROOT.TLine(line_xmin, 1.0, line_xmax, 1.0)
            unity.SetLineStyle(ROOT.kDashed)
            unity.Draw("same")
            draw_ratio_y_title()
            ROOT.gPad.RedrawAxis()

        os.makedirs(os.path.dirname(output_base), exist_ok=True)
        canvas.SaveAs(f"{output_base}_{tag}.pdf")
        canvas.SaveAs(f"{output_base}_{tag}.png")
        canvas.Close()
        print(f"[OK] Wrote {output_base}_{tag}.pdf/.png")


def site_subdir_for_variable(variable: str) -> str:
    has_jet_prefix = variable.startswith(("j1", "j2", "j12", "all"))
    has_wnae = "WNAE" in variable
    has_pnet = "PNetScore" in variable
    if has_pnet and has_jet_prefix:
        return os.path.join("ParticleNET_Supervised_tagger", "ScoreDistribution")
    if variable == "HT":
        return os.path.join("WNAE", "Preselection", "PostPreselecton_Eventlevel")
    if has_wnae and has_jet_prefix:
        return os.path.join("WNAE", "WNAEScores")
    if has_jet_prefix:
        return os.path.join("WNAE", "InputVariables")
    return os.path.join("WNAE", "Preselection", "PostPreselecton_Eventlevel")


def variable_selected_for_site_class(variable: str, site_variable_class: str) -> bool:
    has_wnae = "WNAE" in variable
    if site_variable_class == "wnae-only":
        return has_wnae
    if site_variable_class == "non-wnae":
        return not has_wnae
    if site_variable_class == "softdrop-only":
        return variable in {
            "j1SoftDropMassAK8",
            "j2SoftDropMassAK8",
            "j12SoftDropMassAK8",
            "allSoftDropMassAK8",
        }
    return True


def plot_cache(
    output_path: str,
    plot_dir: str,
    years: Sequence[str],
    requested_signals: Sequence[str],
    also_linear: bool,
    draw_raw: bool = False,
    site_layout: bool = False,
    site_variable_class: str = "all",
    raw_only: bool = False,
    draw_uncertainty_band: bool = True,
    preliminary_label: bool = False,
) -> None:
    root_file = ROOT.TFile.Open(output_path, "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Could not open cache ROOT file: {output_path}")

    try:
        for era in list(years) + ["Run2"]:
            for variable in VARIABLES:
                if not variable_selected_for_site_class(variable, site_variable_class):
                    continue
                variable_plot_dir = (
                    os.path.join(plot_dir, site_subdir_for_variable(variable))
                    if site_layout
                    else plot_dir
                )
                if not raw_only:
                    output_base = os.path.join(
                        variable_plot_dir,
                        era,
                        "mc_signal_only",
                        variable,
                    )
                    draw_stack(
                        root_file,
                        era,
                        variable,
                        output_base,
                        requested_signals,
                        normalized=True,
                        also_linear=also_linear,
                        draw_uncertainty_band=draw_uncertainty_band,
                        preliminary_label=preliminary_label,
                    )
                    output_base = os.path.join(
                        variable_plot_dir,
                        era,
                        "with_data",
                        variable,
                    )
                    draw_stack(
                        root_file,
                        era,
                        variable,
                        output_base,
                        requested_signals,
                        normalized=True,
                        also_linear=also_linear,
                        include_data=True,
                        draw_uncertainty_band=draw_uncertainty_band,
                        preliminary_label=preliminary_label,
                    )
                    output_base = os.path.join(
                        variable_plot_dir,
                        era,
                        "wdata_ratio",
                        variable,
                    )
                    draw_stack(
                        root_file,
                        era,
                        variable,
                        output_base,
                        requested_signals,
                        normalized=True,
                        also_linear=False,
                        include_data=True,
                        include_ratio=True,
                        draw_uncertainty_band=draw_uncertainty_band,
                        preliminary_label=preliminary_label,
                    )
                if draw_raw or raw_only:
                    output_base = os.path.join(
                        variable_plot_dir,
                        era,
                        "raw_wdata_ratio",
                        variable,
                    )
                    draw_stack(
                        root_file,
                        era,
                        variable,
                        output_base,
                        requested_signals,
                        normalized=False,
                        also_linear=also_linear,
                        include_data=True,
                        include_ratio=True,
                        draw_uncertainty_band=draw_uncertainty_band,
                        preliminary_label=preliminary_label,
                    )
    finally:
        root_file.Close()


# =============================================================================
# Build driver
# =============================================================================
def read_sample(
    sample_dir: str,
    display_name: str,
    year: str,
    step_size: str,
    max_files: int,
    strict: bool,
    dry_run: bool,
    failures: List[str],
    is_data: bool = False,
) -> Optional[Dict[str, HistogramAccumulator]]:
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
        message = f"{display_name}: no part-*.root files found"
        if strict:
            raise RuntimeError(message)
        print(f"[WARN] {message}")
        return None

    if dry_run:
        return empty_accumulators()

    accumulators = empty_accumulators()
    entries = 0
    successful = 0

    for index, remote_file in enumerate(files, start=1):
        try:
            entries += stream_file(
                eos_url(remote_file),
                year,
                step_size,
                accumulators,
                is_data=is_data,
            )
            successful += 1
        except Exception as exc:
            message = (
                f"{display_name}/{os.path.basename(remote_file)}: {exc}"
            )
            failures.append(message)
            print(f"[WARN] {message}")

        if index % 25 == 0 or index == len(files):
            print(f"       processed {index:4d}/{len(files):4d} files")

    if strict and successful != len(files):
        raise RuntimeError(
            f"{display_name}: {len(files) - successful} input file(s) failed"
        )

    print(
        f"       Events read: {entries:,} | "
        f"MET yield: {accumulators['MET'].integral():.8g} | "
        f"ST yield: {accumulators['ST'].integral():.8g}"
    )
    return accumulators


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build ROOT histogram cache from selected WNAE backgrounds "
            "and direct t-channel signal skims; make legacy-style normalized stacks."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--background-base",
        default=DEFAULT_BACKGROUND_BASE,
        help="EOS directory immediately before YEAR/t_channel_pre_selection/nominal for backgrounds",
    )
    parser.add_argument(
        "--data-base",
        default=DEFAULT_DATA_BASE,
        help="EOS directory immediately before YEAR/t_channel_pre_selection/nominal for data",
    )
    parser.add_argument(
        "--signal-base",
        default=DEFAULT_SIGNAL_BASE,
        help=(
            "EOS directory immediately before YEAR/t_channel_pre_selection/nominal "
            "or YEAR/nominal for signals"
        ),
    )
    parser.add_argument(
        "--wnae-signal-base",
        default=DEFAULT_WNAE_SIGNAL_BASE,
        help=(
            "Signal EOS directory containing WNAE loss branches. WNAE-score "
            "signal histograms are read from here and replace those variables "
            "from --signal-base."
        ),
    )
    parser.add_argument(
        "--signals",
        default=",".join(DEFAULT_SIGNAL_DIRS),
        help=(
            "Comma-separated full signal directory names or unique substrings. "
            "Use --signals all to process every t-channel signal."
        ),
    )
    parser.add_argument(
        "--no-signals",
        action="store_true",
        help="Skip all signal EOS reads and omit signal overlays",
    )
    parser.add_argument(
        "--data-samples",
        default=",".join(DEFAULT_DATA_SAMPLES),
        help="Comma-separated data sample directories under each nominal skim directory",
    )
    parser.add_argument(
        "--no-data",
        action="store_true",
        help="Skip data EOS reads and do not write with_data plots",
    )
    parser.add_argument(
        "--years",
        default="2016,2017,2018",
        help="Comma-separated eras to process",
    )
    parser.add_argument(
        "--output",
        default="Figure2_plots/WNAE_selected_MC_vars_Run2.root",
        help="Compact ROOT output containing selected TH1 histograms",
    )
    parser.add_argument(
        "--plot-dir",
        default="Figure2_plots",
        help="Directory for normalized PDF/PNG stacks",
    )
    parser.add_argument(
        "--plot-only",
        action="store_true",
        help="Use the cached ROOT histograms to recreate plots without touching EOS",
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Build/update ROOT cache only; do not write PDF/PNG plots",
    )
    parser.add_argument(
        "--also-linear",
        action="store_true",
        help="Additionally make linear-y normalized versions",
    )
    parser.add_argument(
        "--draw-raw",
        action="store_true",
        help="Also draw unnormalized weighted-yield stacks from the raw cache",
    )
    parser.add_argument(
        "--raw-only",
        action="store_true",
        help="Only draw raw data/MC-ratio plots into raw_wdata_ratio subdirectories",
    )
    parser.add_argument(
        "--no-uncertainty-band",
        action="store_true",
        help="Do not draw the simulated-background uncertainty band or its legend entry",
    )
    parser.add_argument(
        "--preliminary-label",
        action="store_true",
        help="Draw CMS Preliminary or CMS Simulation Preliminary",
    )
    parser.add_argument(
        "--site-layout",
        action="store_true",
        help="Route Figure2 variables into site subdirectories by variable type",
    )
    parser.add_argument(
        "--site-variable-class",
        choices=("all", "wnae-only", "non-wnae", "softdrop-only"),
        default="all",
        help="Restrict which variables are drawn when using cached site plots",
    )
    parser.add_argument(
        "--chunk-size",
        default="100 MB",
        help="uproot streaming step size",
    )
    parser.add_argument(
        "--max-files-per-sample",
        type=int,
        default=0,
        help="Testing only: use at most N ROOT files in each sample directory; 0 uses all",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="List matched inputs but do not read trees, write ROOT, or draw plots",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail on a missing requested directory or a failed input ROOT file",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    years = [year.strip() for year in args.years.split(",") if year.strip()]
    invalid = [year for year in years if year not in LUMI_PB]
    if invalid:
        raise SystemExit(
            f"Unsupported year(s): {invalid}. Supported values: {sorted(LUMI_PB)}"
        )

    if args.plot_only:
        if args.no_plots:
            raise SystemExit("--plot-only and --no-plots cannot be used together.")
        requested_signals = [] if args.no_signals else parse_csv(args.signals)
        print(f"[INFO] Plotting only from cached ROOT file: {args.output}")
        plot_cache(
            args.output,
            args.plot_dir,
            years,
            requested_signals,
            args.also_linear,
            args.draw_raw,
            args.site_layout,
            args.site_variable_class,
            args.raw_only,
            not args.no_uncertainty_band,
            args.preliminary_label,
        )
        return

    requested_signals = [] if args.no_signals else parse_csv(args.signals)
    requested_data = [] if args.no_data else parse_csv(args.data_samples)
    failures: List[str] = []
    sample_hists: Dict[
        Tuple[str, str, str], Dict[str, HistogramAccumulator]
    ] = {}
    group_hists: Dict[
        Tuple[str, str], Dict[str, HistogramAccumulator]
    ] = {}
    signal_hists: Dict[
        Tuple[str, str], Dict[str, HistogramAccumulator]
    ] = {}
    data_hists: Dict[
        Tuple[str, str], Dict[str, HistogramAccumulator]
    ] = {}

    print("=" * 100)
    print("WNAE variable cache and normalized log-y Figure 2 stacks")
    print(f"Background base : {normalize_eos_path(args.background_base)}")
    print(f"Data base       : {normalize_eos_path(args.data_base)}")
    print(f"Signal base     : {normalize_eos_path(args.signal_base)}")
    print(f"WNAE signal base: {normalize_eos_path(args.wnae_signal_base)}")
    print(f"Years           : {', '.join(years)}")
    print(
        "Signals         : "
        + (", ".join(requested_signals) if requested_signals else "disabled")
    )
    print(
        "Data            : "
        + (", ".join(requested_data) if requested_data else "disabled")
    )
    print(f"ROOT output     : {args.output}")
    print(f"Plot directory  : {args.plot_dir}")
    print(
        "Weight          : Weight * lumi_pb * puWeight * NonPrefiringProb; "
        "if lundWeightNom and CutFlow are present, multiply "
        "lundWeightNom * Initial / InitialLundNominal"
    )
    print(
        f"Stored variables: {len(VARIABLES)} one-dimensional histograms; "
        f"display rebin factor = {plot_rebin_factor('MET')} (50 visible bins)"
    )
    print(
        "[WARNING] TTJets group deliberately sums every requested TTJets "
        "directory, including potentially overlapping inclusive/enriched samples."
    )
    print("=" * 100)

    # ----------------------------- Backgrounds -----------------------------
    for year in years:
        nominal = (
            f"{normalize_eos_path(args.background_base)}/{year}/"
            "t_channel_pre_selection/nominal"
        )
        print(f"\n{'#' * 100}\nBACKGROUND {year}\n{nominal}\n{'#' * 100}")

        for process in PROCESS_ORDER:
            for sample in SAMPLE_MANIFEST[process]:
                sample_dir = f"{nominal}/{sample}"
                display_name = f"{year}/{process}/{sample}"
                print(f"\n[READ] {display_name}")

                accumulators = read_sample(
                    sample_dir,
                    display_name,
                    year,
                    args.chunk_size,
                    args.max_files_per_sample,
                    args.strict,
                    args.dry_run,
                    failures,
                )
                if accumulators is None:
                    continue

                sample_hists[(year, process, sample)] = accumulators
                group_hists.setdefault((year, process), empty_accumulators())
                add_accumulator_dict(group_hists[(year, process)], accumulators)

    # ------------------------------- Signals -------------------------------
    if requested_signals:
        for year in years:
            print(
                f"\n{'#' * 100}\nSIGNALS {year}\n"
                f"{normalize_eos_path(args.signal_base)}/{year}/"
                f"[t_channel_pre_selection/nominal or nominal]\n{'#' * 100}"
            )

            matches = get_signal_directories(
                args.signal_base,
                year,
                requested_signals,
                args.strict,
            )
            if not matches:
                print("[INFO] No requested signal directories found for this year.")
                continue

            for signal_name, sample_dir in matches:
                display_name = f"{year}/signal/{signal_name}"
                print(f"\n[READ] {display_name}")

                accumulators = read_sample(
                    sample_dir,
                    display_name,
                    year,
                    args.chunk_size,
                    args.max_files_per_sample,
                    args.strict,
                    args.dry_run,
                    failures,
                )
                if accumulators is not None:
                    signal_hists[(year, signal_name)] = accumulators

            wnae_variables = wnae_score_variables()
            if args.wnae_signal_base and wnae_variables:
                wnae_matches = get_signal_directories(
                    args.wnae_signal_base,
                    year,
                    requested_signals,
                    args.strict,
                )
                for signal_name, sample_dir in wnae_matches:
                    if (year, signal_name) not in signal_hists:
                        continue
                    display_name = f"{year}/signal-wnae/{signal_name}"
                    print(f"\n[READ] {display_name}")

                    accumulators = read_sample(
                        sample_dir,
                        display_name,
                        year,
                        args.chunk_size,
                        args.max_files_per_sample,
                        args.strict,
                        args.dry_run,
                        failures,
                    )
                    if accumulators is not None:
                        replace_accumulator_variables(
                            signal_hists[(year, signal_name)],
                            accumulators,
                            wnae_variables,
                        )

    # -------------------------------- Data ---------------------------------
    if requested_data:
        for year in years:
            nominal = (
                f"{normalize_eos_path(args.data_base)}/{year}/"
                "t_channel_pre_selection/nominal"
            )
            print(f"\n{'#' * 100}\nDATA {year}\n{nominal}\n{'#' * 100}")

            for sample in requested_data:
                sample_dir = f"{nominal}/{sample}"
                display_name = f"{year}/data/{sample}"
                print(f"\n[READ] {display_name}")

                accumulators = read_sample(
                    sample_dir,
                    display_name,
                    year,
                    args.chunk_size,
                    args.max_files_per_sample,
                    args.strict,
                    args.dry_run,
                    failures,
                    is_data=True,
                )
                if accumulators is not None:
                    data_hists[(year, sample)] = accumulators

    if args.dry_run:
        print("\n[DRY RUN] Nothing was read from ROOT trees and no output was written.")
        return

    if not group_hists:
        raise RuntimeError("No background histograms were made.")

    # --------------------------- Run-2 combinations ------------------------
    run2_hists: Dict[str, Dict[str, HistogramAccumulator]] = {}
    for (_, process), accumulators in group_hists.items():
        run2_hists.setdefault(process, empty_accumulators())
        add_accumulator_dict(run2_hists[process], accumulators)

    run2_signal_hists: Dict[str, Dict[str, HistogramAccumulator]] = {}
    for (_, signal_name), accumulators in signal_hists.items():
        run2_signal_hists.setdefault(signal_name, empty_accumulators())
        add_accumulator_dict(run2_signal_hists[signal_name], accumulators)

    run2_data_hists = empty_accumulators() if data_hists else None
    if run2_data_hists:
        for accumulators in data_hists.values():
            add_accumulator_dict(run2_data_hists, accumulators)

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
    print(f"\n[OK] Wrote compact ROOT histogram cache:\n  {args.output}")

    if not args.no_plots:
        # Reopen the actual stored cache. Plotting therefore validates the
        # no-reread workflow and later --plot-only runs are identical.
        plot_cache(
            args.output,
            args.plot_dir,
            years,
            requested_signals,
            args.also_linear,
            args.draw_raw,
            args.site_layout,
            args.site_variable_class,
            args.raw_only,
            not args.no_uncertainty_band,
            args.preliminary_label,
        )

    print("\n[RUN-2 RAW INTEGRALS]")
    for process in PROCESS_ORDER:
        if process in run2_hists:
            print(
                f"  {process:12s}"
                f" MET={run2_hists[process]['MET'].integral():.10g}"
                f" ST={run2_hists[process]['ST'].integral():.10g}"
            )

    if run2_signal_hists:
        print("\n[RUN-2 RAW SIGNAL INTEGRALS]")
        for signal_name in ordered_signal_names(run2_signal_hists.keys()):
            print(
                f"  {signal_name}\n"
                f"    MET={run2_signal_hists[signal_name]['MET'].integral():.10g}"
                f" ST={run2_signal_hists[signal_name]['ST'].integral():.10g}"
            )

    if failures:
        print("\n[WARNING] Files skipped after I/O errors:")
        for failure in failures:
            print(f"  - {failure}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        raise SystemExit("\nStopped by user.")
    except Exception as exc:
        print(f"\n[FATAL] {exc}", file=sys.stderr)
        raise
