#!/usr/bin/env python3
"""
Build WNAE-vs-PNET jet-score heatmaps in WNAE pT bins.

The script streams the WNAE skims, chooses the WNAE score branch that matches
each AK8 jet pT bin, and fills 2D histograms of

    x = JetsAK8_pNetJetTaggerScore
    y = WNAE Score for the corresponding pT bin

Outputs:
  * WNAE_PNET_comparison_histograms.root
  * WNAE_PNET_comparison_<pt-bin>.png/pdf
  * WNAE_PNET_comparison_pass_fail.png/pdf

Example: Signal
python3 -u WNAE_PNET_comparison.py \
  --base /store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/signals_WNAE \
  --pnet-friend-base /store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto \
  --years 2016,2017,2018 \
  --processes Signal \
  --sample-dirs all \
  --signal-mdark 20 \
  --signal-yukawa 1 \
  --out-dir WNAE_PNET_comparison_Run2_Signal_mDark20_yukawa1 \
  --log-z \
  --log-matrix-z \
  --chunk-size "100 MB"

python3 -u WNAE_PNET_comparison.py \
  --plot-only \
  --cache WNAE_PNET_comparison_Run2_AllSignals/WNAE_PNET_comparison_histograms.root \
  --out-dir WNAE_PNET_comparison_Run2_AllSignals \
  --log-z \
  --log-matrix-z


example data
python3 -u WNAE_PNET_comparison.py \
  --years 2016,2017,2018 \
  --processes Data \
  --out-dir WNAE_PNET_comparison_Run2_Data \
  --chunk-size "100 MB"

python3 -u WNAE_PNET_comparison.py \
  --plot-only \
  --plot-ymax 100 \
  --log-z \
  --cache WNAE_PNET_comparison_Run2_Data/WNAE_PNET_comparison_histograms.root \
  --out-dir WNAE_PNET_comparison_Run2_Data
Quick test:
  python3 -u WNAE_PNET_comparison.py --years 2018 --processes QCD \
    --max-files-per-sample 1 --max-events 20000

"""

from __future__ import annotations

import argparse
from contextlib import ExitStack
import json
import os
import re
import subprocess
import sys
from dataclasses import asdict, dataclass
from typing import Dict, List, Optional, Sequence, Tuple


def require_runtime_packages():
    missing = []
    try:
        import awkward as ak  # noqa: F401
    except Exception:
        missing.append("awkward")
    try:
        import numpy as np  # noqa: F401
    except Exception:
        missing.append("numpy")
    try:
        import uproot  # noqa: F401
    except Exception:
        missing.append("uproot")
    try:
        import matplotlib  # noqa: F401
    except Exception:
        missing.append("matplotlib")

    if missing:
        raise RuntimeError(
            "Missing required Python package(s): "
            + ", ".join(missing)
            + ". Source the analysis/condor environment before running this script."
        )


EOS_HOST = "root://cmseos.fnal.gov"
DEFAULT_BASE = "/store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/data_mc_WNAE"
TREE_NAME = "Events"

PT_BRANCH = "JetsAK8_/.fPt"
PNET_BRANCH = "JetsAK8_pNetJetTaggerScore"
GOOD_JET_BRANCH = "JetsAK8_isGood"
SVJ_CATEGORY_KEYS = ["0SVJ", "1SVJ", "2SVJ", "3PSVJ"]
SVJ_CATEGORY_LABELS = ["0", "1", "2", "3+"]

LUMI_PB = {
    "2016": 36.31 * 1000.0,
    "2017": 42.07 * 1000.0,
    "2018": 59.56 * 1000.0,
}


@dataclass(frozen=True)
class WNAEPtBin:
    key: str
    label: str
    pt_min: float
    pt_max: Optional[float]
    loss_branch: str
    mc_wp: float
    data_wp: float


WNAE_PT_BINS = [
    WNAEPtBin(
        key="pt0to200",
        label="0 <= pT < 200 GeV",
        pt_min=0.0,
        pt_max=200.0,
        loss_branch="JetsAK8_WNAEPt0To200Loss",
        mc_wp=25.156,
        data_wp=24.974,
    ),
    WNAEPtBin(
        key="pt200to300",
        label="200 <= pT < 300 GeV",
        pt_min=200.0,
        pt_max=300.0,
        loss_branch="JetsAK8_WNAEPt200To300Loss",
        mc_wp=18.284,
        data_wp=18.445,
    ),
    WNAEPtBin(
        key="pt300to400",
        label="300 <= pT < 400 GeV",
        pt_min=300.0,
        pt_max=400.0,
        loss_branch="JetsAK8_WNAEPt300To400Loss",
        mc_wp=20.383,
        data_wp=21.031,
    ),
    WNAEPtBin(
        key="pt400to500",
        label="400 <= pT < 500 GeV",
        pt_min=400.0,
        pt_max=500.0,
        loss_branch="JetsAK8_WNAEPt400To500Loss",
        mc_wp=21.941,
        data_wp=22.131,
    ),
    WNAEPtBin(
        key="pt500toInf",
        label="pT >= 500 GeV",
        pt_min=500.0,
        pt_max=None,
        loss_branch="JetsAK8_WNAEPt500ToInfLoss",
        mc_wp=16.370,
        data_wp=16.535,
    ),
]

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
    "Data": ["HTMHT", "JetHT", "MET"],
    "Signal": ["t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1"],
}


class Hist2D:
    def __init__(self, x_edges, y_edges):
        import numpy as np

        self.x_edges = np.asarray(x_edges, dtype="float64")
        self.y_edges = np.asarray(y_edges, dtype="float64")
        self.sumw = np.zeros((len(self.x_edges) - 1, len(self.y_edges) - 1), dtype="float64")
        self.sumw2 = np.zeros_like(self.sumw)
        self.entries_total = 0
        self.entries_in_range = 0

    def fill(self, x, y, weights):
        import numpy as np

        x = np.asarray(x, dtype="float64")
        y = np.asarray(y, dtype="float64")
        weights = np.asarray(weights, dtype="float64")
        self.entries_total += int(x.size)

        mask = (
            np.isfinite(x)
            & np.isfinite(y)
            & np.isfinite(weights)
            & (x >= self.x_edges[0])
            & (x < self.x_edges[-1])
            & (y >= self.y_edges[0])
            & (y < self.y_edges[-1])
        )
        if not np.any(mask):
            return

        x = x[mask]
        y = y[mask]
        weights = weights[mask]
        self.entries_in_range += int(x.size)
        self.sumw += np.histogram2d(x, y, bins=(self.x_edges, self.y_edges), weights=weights)[0]
        self.sumw2 += np.histogram2d(
            x, y, bins=(self.x_edges, self.y_edges), weights=weights * weights
        )[0]

    def add(self, other):
        import numpy as np

        if not np.allclose(self.x_edges, other.x_edges) or not np.allclose(
            self.y_edges, other.y_edges
        ):
            raise ValueError("Cannot add incompatible 2D histograms.")
        self.sumw += other.sumw
        self.sumw2 += other.sumw2
        self.entries_total += other.entries_total
        self.entries_in_range += other.entries_in_range

    def integral(self):
        import numpy as np

        return float(np.sum(self.sumw))


class PassFailAccumulator:
    def __init__(self):
        import numpy as np

        self.sumw = np.zeros((2, 2), dtype="float64")
        self.sumw2 = np.zeros((2, 2), dtype="float64")
        self.entries = np.zeros((2, 2), dtype="int64")

    def fill(self, pnet_pass, wnae_pass, weights):
        import numpy as np

        pnet_pass = np.asarray(pnet_pass, dtype=bool)
        wnae_pass = np.asarray(wnae_pass, dtype=bool)
        weights = np.asarray(weights, dtype="float64")
        mask = np.isfinite(weights)
        if not np.any(mask):
            return

        x_index = pnet_pass[mask].astype("int64")
        y_index = wnae_pass[mask].astype("int64")
        weights = weights[mask]

        for y in (0, 1):
            for x in (0, 1):
                cell = (y_index == y) & (x_index == x)
                if np.any(cell):
                    self.sumw[y, x] += float(np.sum(weights[cell]))
                    self.sumw2[y, x] += float(np.sum(weights[cell] * weights[cell]))
                    self.entries[y, x] += int(np.count_nonzero(cell))


class EventSVJAccumulator:
    def __init__(self):
        import numpy as np

        self.sumw = np.zeros((4, 4), dtype="float64")
        self.sumw2 = np.zeros((4, 4), dtype="float64")
        self.entries = np.zeros((4, 4), dtype="int64")

    def fill(self, pnet_counts, wnae_counts, weights):
        import numpy as np

        pnet_bins = np.minimum(np.asarray(pnet_counts, dtype="int64"), 3)
        wnae_bins = np.minimum(np.asarray(wnae_counts, dtype="int64"), 3)
        weights = np.asarray(weights, dtype="float64")
        mask = np.isfinite(weights)
        if not np.any(mask):
            return

        pnet_bins = pnet_bins[mask]
        wnae_bins = wnae_bins[mask]
        weights = weights[mask]

        for y in range(4):
            for x in range(4):
                cell = (wnae_bins == y) & (pnet_bins == x)
                if np.any(cell):
                    self.sumw[y, x] += float(np.sum(weights[cell]))
                    self.sumw2[y, x] += float(np.sum(weights[cell] * weights[cell]))
                    self.entries[y, x] += int(np.count_nonzero(cell))


def signal_group_from_sample(sample: str) -> str:
    match = re.search(r"(rinv-[^_]+)", sample)
    if match:
        return match.group(1)
    return "signal_other"


def sample_matches_signal_filters(sample: str, args) -> bool:
    if args.signal_mdark and not re.search(rf"(?:^|_)mDark-{re.escape(args.signal_mdark)}(?:_|$)", sample):
        return False
    if args.signal_yukawa and not re.search(rf"(?:^|_)yukawa-{re.escape(args.signal_yukawa)}(?:_|$)", sample):
        return False
    return True


def filename_safe(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def signal_mmed_label(metadata: dict) -> str:
    label = str(metadata.get("signal_mmed_label") or "").strip()
    if not label:
        return ""
    if label.lower() in {"all", "all_mphi", "all_mmed"}:
        return r"all $m_{\phi}$"
    return r"$m_{\phi} = " + label + r"$ GeV"


def signal_page_label(page_label: str, metadata: dict) -> str:
    pieces = []
    is_signal = any(
        str(process).lower() == "signal" for process in metadata.get("processes", [])
    ) or bool(metadata.get("event_svj_group_keys"))
    mmed_label = signal_mmed_label(metadata)
    if mmed_label:
        pieces.append(mmed_label)
    signal_mdark = metadata.get("signal_mdark") or ("20" if is_signal else "")
    signal_yukawa = metadata.get("signal_yukawa") or ("1" if is_signal else "")
    if signal_mdark:
        pieces.append(r"$m_{\mathrm{dark}} = " + str(signal_mdark) + r"$ GeV")
    if page_label and page_label != "inclusive":
        rinv_value = page_label.replace("rinv-", "").replace("p", ".")
        pieces.append(r"$r_{\mathrm{inv}} = " + rinv_value + r"$")
    if signal_yukawa:
        pieces.append(r"$\lambda = " + str(signal_yukawa) + r"$")
    return ", ".join(pieces)


def is_signal_comparison(metadata: dict) -> bool:
    return any(
        str(process).lower() == "signal" for process in metadata.get("processes", [])
    ) or bool(metadata.get("event_svj_group_keys"))


def run2_lumi_label(metadata: dict) -> str:
    years = set(str(year) for year in metadata.get("years", []))
    if years == {"2016", "2017", "2018"}:
        return r"138 fb$^{-1}$ (13 TeV)"
    lumi_pb = sum(LUMI_PB.get(year, 0.0) for year in years)
    if lumi_pb <= 0.0:
        return "(13 TeV)"
    lumi_fb = lumi_pb / 1000.0
    if abs(lumi_fb - round(lumi_fb)) < 0.05:
        return rf"{int(round(lumi_fb))} fb$^{{-1}}$ (13 TeV)"
    return rf"{lumi_fb:.1f} fb$^{{-1}}$ (13 TeV)"


def add_cms_label(
    ax,
    metadata: dict,
    has_mplhep: bool,
    page_label: str = "",
    matrix_style: bool = False,
) -> None:
    is_signal = is_signal_comparison(metadata)
    if matrix_style:
        header_y = 1.005
        ax.text(
            0.0,
            header_y,
            "CMS",
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=34,
            fontweight="bold",
        )
        if is_signal:
            ax.text(
                0.185,
                header_y,
                "Simulation",
                transform=ax.transAxes,
                ha="left",
                va="bottom",
                fontsize=22,
                style="italic",
            )
        right_label = run2_lumi_label(metadata)
        ax.text(
            1.0,
            header_y,
            right_label,
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=24,
        )
        return

    if is_signal:
        header_y = 1.018
        cms_size = 20
        sim_size = 16
        tev_size = 16
        sim_x = 0.13
        ax.text(
            0.0,
            header_y,
            "CMS",
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=cms_size,
            fontweight="bold",
        )
        ax.text(
            sim_x,
            header_y,
            "Simulation Preliminary",
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=max(12, sim_size - 2),
            style="italic",
        )
        ax.text(
            1.0,
            header_y,
            "(13 TeV)",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=tev_size,
        )
    elif has_mplhep and not metadata.get("all_data", False):
        import mplhep as hep

        hep.cms.label("Preliminary", data=False, ax=ax, loc=0, fontsize=12)
    else:
        ax.text(
            0.0,
            1.005,
            "CMS",
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=14,
            fontweight="bold",
        )
    if metadata.get("all_data", False):
        ax.text(
            1.0,
            1.005,
            run2_lumi_label(metadata),
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=12,
        )
    if not is_signal:
        extra_label = signal_page_label(page_label, metadata)
        if extra_label:
            ax.text(
                0.52,
                1.005,
                extra_label,
                transform=ax.transAxes,
                ha="left",
                va="bottom",
                fontsize=11,
            )


def apply_figure_layout(fig, metadata: dict) -> None:
    if is_signal_comparison(metadata):
        fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.965), pad=0.8)
    else:
        fig.tight_layout(pad=0.8)


def apply_matrix_layout(fig) -> None:
    fig.subplots_adjust(left=0.105, right=0.875, bottom=0.165, top=0.865)


def binned_weighted_kendall_tau(hist: Hist2D) -> Optional[float]:
    import numpy as np

    weights = np.asarray(hist.sumw, dtype="float64")
    if weights.size == 0:
        return None
    weights = np.where(np.isfinite(weights) & (weights > 0.0), weights, 0.0)
    total = float(np.sum(weights))
    if total <= 0.0:
        return None

    nx, ny = weights.shape
    prefix = np.zeros((nx + 1, ny + 1), dtype="float64")
    prefix[1:, 1:] = np.cumsum(np.cumsum(weights, axis=0), axis=1)
    concordant = 0.0
    discordant = 0.0

    for ix in range(nx):
        for iy in range(ny):
            weight = float(weights[ix, iy])
            if weight <= 0.0:
                continue
            lower_left = prefix[ix, iy]
            upper_right = total - prefix[ix + 1, ny] - prefix[nx, iy + 1] + prefix[ix + 1, iy + 1]
            upper_left = prefix[ix, ny] - prefix[ix, iy + 1]
            lower_right = prefix[nx, iy] - prefix[ix + 1, iy]
            concordant += weight * (lower_left + upper_right)
            discordant += weight * (upper_left + lower_right)

    denom = concordant + discordant
    if denom <= 0.0:
        return None
    return float((concordant - discordant) / denom)


def add_score_info_box(ax, hist: Hist2D, metadata: dict, page_label: str = "") -> None:
    if not is_signal_comparison(metadata):
        return

    lines = []
    signal_label = signal_page_label(page_label, metadata)
    if signal_label:
        lines.append(signal_label)
    tau = binned_weighted_kendall_tau(hist)
    if tau is not None:
        lines.append(r"Kendall $\tau = " + f"{tau:.3f}" + r"$")
    if not lines:
        return

    ax.text(
        0.025,
        0.965,
        "\n".join(lines),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10.5,
        color="black",
        bbox={
            "boxstyle": "square,pad=0.25",
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.82,
        },
        zorder=10,
    )


def add_wp_label_box(
    ax,
    *,
    pnet_wp: Optional[float],
    wnae_wp: Optional[float],
    is_quantile_axis: bool,
    x_edges: np.ndarray,
    plot_ymin: float,
    plot_ymax: float,
) -> None:
    lines = []
    if wnae_wp is not None:
        if is_quantile_axis:
            lines.append(f"WNAE WP q = {wnae_wp:.2f}")
        else:
            lines.append(f"WNAE WP {wnae_wp:.2f}")
    if pnet_wp is not None:
        lines.append(f"PNET WP {pnet_wp:.2f}")
    if not lines:
        return

    x_span = x_edges[-1] - x_edges[0]
    y_span = plot_ymax - plot_ymin
    x = x_edges[0] + 0.72 * x_span
    y = plot_ymin + 0.78 * y_span
    if pnet_wp is not None and wnae_wp is not None:
        x = min(max(x_edges[0] + 0.05 * x_span, pnet_wp - 0.03 * x_span), x_edges[-1] - 0.05 * x_span)
        y = min(max(plot_ymin + 0.06 * y_span, wnae_wp + 0.06 * y_span), plot_ymax - 0.06 * y_span)

    ax.text(
        x,
        y,
        "\n".join(lines),
        color="black",
        fontsize=11,
        ha="right",
        va="bottom",
        bbox={
            "boxstyle": "square,pad=0.25",
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.88,
        },
        zorder=12,
    )


def add_signal_matrix_info_box(ax, metadata: dict, page_label: str = "") -> None:
    if not is_signal_comparison(metadata):
        return

    mmed_label = signal_mmed_label(metadata) or r"$m_{\phi} = 2000$ GeV"
    signal_mdark = metadata.get("signal_mdark") or "20"
    signal_yukawa = metadata.get("signal_yukawa") or "1"
    top_line = mmed_label + r", $m_{\mathrm{dark}} = " + str(signal_mdark) + r"$ GeV"
    bottom_pieces = []
    if page_label and page_label != "inclusive":
        rinv_value = page_label.replace("rinv-", "").replace("p", ".")
        bottom_pieces.append(r"$r_{\mathrm{inv}} = " + rinv_value + r"$")
    bottom_pieces.append(r"$\lambda = " + str(signal_yukawa) + r"$")
    label = top_line + "\n" + ", ".join(bottom_pieces)
    ax.text(
        0.0,
        -0.085,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=16,
        color="black",
        bbox={
            "boxstyle": "square,pad=0.25",
            "facecolor": "white",
            "edgecolor": "none",
            "alpha": 0.82,
        },
        zorder=10,
    )


def formatted_cell_value(value: float, metadata: dict) -> str:
    if metadata.get("all_data", False):
        return f"{value:.0f}"
    return f"{value:.1f}"


def parse_csv(value: str) -> List[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def normalize_eos_path(path: str) -> str:
    path = os.path.expanduser(path).rstrip("/")
    mounted_prefix = "/eos/uscms"
    if path.startswith(mounted_prefix + "/"):
        path = path[len(mounted_prefix) :]
    if path.startswith("root://"):
        match = re.match(r"^root://[^/]+/+(.*)$", path)
        if not match:
            raise ValueError(f"Could not normalize xrootd path: {path}")
        path = "/" + match.group(1).lstrip("/")
    if not path.startswith("/"):
        raise ValueError(f"Expected an absolute EOS path, received: {path}")
    return path


def eos_url(path: str) -> str:
    return EOS_HOST + "//" + normalize_eos_path(path).lstrip("/")


def xrdfs_ls(path: str, recursive: bool = False) -> List[str]:
    command = ["xrdfs", EOS_HOST, "ls"]
    if recursive:
        command.append("-R")
    command.append(normalize_eos_path(path))
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"xrdfs listing failed for {path}\n"
            f"command: {' '.join(command)}\n"
            f"stderr: {result.stderr.strip()}"
        )
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def year_nominal_dirs(base: str, year: str) -> List[str]:
    normalized = normalize_eos_path(base)
    return [
        f"{normalized}/{year}/t_channel_pre_selection/nominal",
        f"{normalized}/{year}/nominal",
    ]


def resolve_sample_dirs(base: str, year: str, samples: Sequence[str]) -> List[Tuple[str, str]]:
    resolved = []
    for nominal_dir in year_nominal_dirs(base, year):
        try:
            entries = xrdfs_ls(nominal_dir, recursive=False)
        except Exception:
            continue

        available = {
            os.path.basename(entry.rstrip("/")): entry.rstrip("/")
            for entry in entries
        }

        if len(samples) == 1 and samples[0].lower() in {"all", "*"}:
            for sample_name in sorted(available):
                if sample_name.startswith("t-channel_"):
                    resolved.append((sample_name, available[sample_name]))
            if resolved:
                return resolved

        for sample in samples:
            if sample in available:
                resolved.append((sample, available[sample]))
            else:
                matches = [name for name in available if sample.lower() in name.lower()]
                for name in sorted(matches):
                    resolved.append((name, available[name]))

        if resolved:
            return resolved

    return resolved


def discover_files_from_base(
    base: str,
    years: Sequence[str],
    processes: Sequence[str],
    sample_dirs: Sequence[str],
    max_files_per_sample: int,
) -> List[Tuple[str, str, str, str]]:
    discovered = []
    for year in years:
        for process in processes:
            samples = list(sample_dirs) if sample_dirs else SAMPLE_MANIFEST.get(process)
            if not samples:
                raise ValueError(
                    f"No samples known for process '{process}'. "
                    "Use --sample-dirs for custom sample directory names."
                )

            for sample, sample_dir in resolve_sample_dirs(base, year, samples):
                try:
                    paths = xrdfs_ls(sample_dir, recursive=True)
                except Exception as exc:
                    print(f"[WARN] {year} {process}/{sample}: {exc}")
                    continue

                root_files = sorted(
                    path for path in paths if re.search(r"/part-[0-9]+\.root$", path)
                )
                if max_files_per_sample > 0:
                    root_files = root_files[:max_files_per_sample]
                for path in root_files:
                    discovered.append((year, process, sample, eos_url(path)))

    return discovered


def files_from_text(path: str) -> List[Tuple[str, str, str, str]]:
    out = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            raw = raw.strip()
            if not raw or raw.startswith("#"):
                continue
            parts = [part.strip() for part in raw.split(",")]
            if len(parts) == 1:
                out.append(("unknown", "input", "input", parts[0]))
            elif len(parts) >= 4:
                out.append((parts[0], parts[1], parts[2], parts[3]))
            else:
                raise ValueError(
                    f"Bad file-list line '{raw}'. Use either FILE or YEAR,PROCESS,SAMPLE,FILE."
                )
    return out


def friend_file_url(friend_base: str, year: str, sample: str, file_url: str) -> str:
    part_name = os.path.basename(normalize_eos_path(file_url))
    nominal_dir = year_nominal_dirs(friend_base, year)[0]
    return eos_url(f"{nominal_dir}/{sample}/{part_name}")


def selected_weight_branches(tree, is_data: bool, no_weights: bool, no_lund_correction: bool) -> List[str]:
    if is_data or no_weights:
        return []
    available = set(tree.keys())
    branches = []
    if "Weight" in available:
        branches.append("Weight")
    for name in ("puWeight", "NonPrefiringProb"):
        if name in available:
            branches.append(name)
    if not no_lund_correction and "lundWeightNom" in available:
        branches.append("lundWeightNom")
    return branches


def lund_nominal_norm(root_file):
    import numpy as np

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


def chunk_event_weights(arrays, year: str, is_data: bool, no_weights: bool, lund_norm=None):
    import numpy as np

    if is_data or no_weights or "Weight" not in arrays.fields:
        return np.ones(len(arrays), dtype="float64")

    import awkward as ak

    weights = np.asarray(ak.to_numpy(arrays["Weight"]), dtype="float64")
    if year in LUMI_PB:
        weights = weights * LUMI_PB[year]
    if "puWeight" in arrays.fields:
        weights = weights * np.asarray(ak.to_numpy(arrays["puWeight"]), dtype="float64")
    if "NonPrefiringProb" in arrays.fields:
        weights = weights * np.asarray(ak.to_numpy(arrays["NonPrefiringProb"]), dtype="float64")
    if lund_norm is not None and "lundWeightNom" in arrays.fields:
        weights = weights * np.asarray(ak.to_numpy(arrays["lundWeightNom"]), dtype="float64") * lund_norm
    return weights


def finite_flatten(values):
    import awkward as ak
    import numpy as np

    return np.asarray(ak.to_numpy(ak.flatten(values, axis=None)), dtype="float64")


def stream_file(
    file_url: str,
    year: str,
    process: str,
    sample: str,
    args,
    histograms: Dict[str, Hist2D],
    histogram_groups: Optional[Dict[str, Dict[str, Hist2D]]] = None,
    pass_fail: Optional[Dict[str, PassFailAccumulator]] = None,
    event_svj: Optional[EventSVJAccumulator] = None,
    event_svj_groups: Optional[Dict[str, EventSVJAccumulator]] = None,
    event_svj_group_key: Optional[str] = None,
) -> int:
    import awkward as ak
    import numpy as np
    import uproot

    is_data = process.lower() == "data"
    entries_seen = 0

    with ExitStack() as stack:
        root_file = stack.enter_context(uproot.open(file_url))
        if args.tree not in root_file:
            raise KeyError(f"No {args.tree} tree found.")
        tree = root_file[args.tree]
        available = set(tree.keys())

        required = [PT_BRANCH]
        missing = [branch for branch in required if branch not in available]
        if missing:
            raise KeyError(f"Missing required branch(es): {missing}")

        friend_tree = None
        if PNET_BRANCH not in available:
            if not args.pnet_friend_base:
                raise KeyError(f"Missing required branch(es): ['{PNET_BRANCH}']")
            friend_url = friend_file_url(args.pnet_friend_base, year, sample=sample, file_url=file_url)
            friend_root = stack.enter_context(uproot.open(friend_url))
            if args.tree not in friend_root:
                raise KeyError(f"No {args.tree} tree found in PNET friend file: {friend_url}")
            friend_tree = friend_root[args.tree]
            if PNET_BRANCH not in friend_tree.keys():
                raise KeyError(f"PNET friend file lacks {PNET_BRANCH}: {friend_url}")
            if friend_tree.num_entries != tree.num_entries:
                print(
                    f"[WARN] friend entry mismatch for {file_url}: "
                    f"WNAE={tree.num_entries}, PNET={friend_tree.num_entries}"
                )

        loss_branches = [pt_bin.loss_branch for pt_bin in WNAE_PT_BINS]
        present_loss_branches = [branch for branch in loss_branches if branch in available]
        if not present_loss_branches:
            raise KeyError("No WNAE loss branches found in Events tree.")

        expressions = [PT_BRANCH] + present_loss_branches
        if PNET_BRANCH in available:
            expressions.append(PNET_BRANCH)
        if GOOD_JET_BRANCH in available:
            expressions.append(GOOD_JET_BRANCH)
        expressions.extend(
            selected_weight_branches(
                tree,
                is_data=is_data,
                no_weights=args.no_weights,
                no_lund_correction=args.no_lund_correction,
            )
        )
        expressions = list(dict.fromkeys(expressions))
        lund_norm = lund_nominal_norm(root_file) if not is_data and "lundWeightNom" in expressions else None

        for arrays in tree.iterate(
            expressions=expressions,
            step_size=args.chunk_size,
            library="ak",
            entry_stop=args.max_events if args.max_events > 0 else None,
            report=True,
        ):
            arrays, report = arrays
            entries_seen += len(arrays)
            event_weights = chunk_event_weights(
                arrays, year=year, is_data=is_data, no_weights=args.no_weights, lund_norm=lund_norm
            )

            pt = arrays[PT_BRANCH]
            if friend_tree is not None:
                friend_arrays = friend_tree.arrays(
                    [PNET_BRANCH],
                    entry_start=report.tree_entry_start,
                    entry_stop=report.tree_entry_stop,
                    library="ak",
                )
                pnet = friend_arrays[PNET_BRANCH]
            else:
                pnet = arrays[PNET_BRANCH]
            if GOOD_JET_BRANCH in arrays.fields:
                good = arrays[GOOD_JET_BRANCH]
                pt = pt[good]
                pnet = pnet[good]
            else:
                good = None

            if event_svj is not None or (
                event_svj_groups is not None and event_svj_group_key is not None
            ):
                pnet_counts = ak.to_numpy(ak.sum(pnet >= args.pnet_wp, axis=1))
                wnae_counts = np.zeros(len(arrays), dtype="int64")
                for count_pt_bin in WNAE_PT_BINS:
                    if count_pt_bin.loss_branch not in arrays.fields:
                        continue
                    count_loss = arrays[count_pt_bin.loss_branch]
                    if good is not None:
                        count_loss = count_loss[good]

                    if count_pt_bin.pt_max is None:
                        count_in_pt = pt >= count_pt_bin.pt_min
                    else:
                        count_in_pt = (pt >= count_pt_bin.pt_min) & (pt < count_pt_bin.pt_max)

                    count_wp = count_pt_bin.data_wp if is_data else count_pt_bin.mc_wp
                    count_tags = count_in_pt & (count_loss >= count_wp)
                    wnae_counts += ak.to_numpy(ak.sum(count_tags, axis=1))

                if event_svj is not None:
                    event_svj.fill(pnet_counts, wnae_counts, event_weights)
                if event_svj_groups is not None and event_svj_group_key is not None:
                    event_svj_groups.setdefault(
                        event_svj_group_key, EventSVJAccumulator()
                    ).fill(pnet_counts, wnae_counts, event_weights)

            for pt_bin in WNAE_PT_BINS:
                if pt_bin.loss_branch not in arrays.fields:
                    continue

                loss = arrays[pt_bin.loss_branch]
                if good is not None:
                    loss = loss[good]

                if pt_bin.pt_max is None:
                    in_pt = pt >= pt_bin.pt_min
                else:
                    in_pt = (pt >= pt_bin.pt_min) & (pt < pt_bin.pt_max)

                selected_pnet = pnet[in_pt]
                selected_loss = loss[in_pt]

                if pass_fail is not None:
                    _, pf_jet_weights = ak.broadcast_arrays(selected_pnet, ak.Array(event_weights))
                    pf_x = finite_flatten(selected_pnet)
                    pf_y = finite_flatten(selected_loss)
                    pf_w = finite_flatten(pf_jet_weights)
                    wnae_wp = pt_bin.data_wp if is_data else pt_bin.mc_wp
                    pass_fail[pt_bin.key].fill(pf_x >= args.pnet_wp, pf_y >= wnae_wp, pf_w)
                    pass_fail["inclusive"].fill(pf_x >= args.pnet_wp, pf_y >= wnae_wp, pf_w)

                if args.tag_region != "all":
                    pnet_tag = selected_pnet >= args.pnet_wp
                    wnae_wp = pt_bin.data_wp if is_data else pt_bin.mc_wp
                    wnae_tag = selected_loss >= wnae_wp
                    if args.tag_region == "pnet":
                        tag_mask = pnet_tag
                    elif args.tag_region == "wnae":
                        tag_mask = wnae_tag
                    elif args.tag_region == "both":
                        tag_mask = pnet_tag & wnae_tag
                    elif args.tag_region == "either":
                        tag_mask = pnet_tag | wnae_tag
                    else:
                        raise ValueError(f"Unknown tag region: {args.tag_region}")
                    selected_pnet = selected_pnet[tag_mask]
                    selected_loss = selected_loss[tag_mask]

                _, jet_weights = ak.broadcast_arrays(selected_pnet, ak.Array(event_weights))

                x = finite_flatten(selected_pnet)
                y = finite_flatten(selected_loss)
                w = finite_flatten(jet_weights)

                histograms[pt_bin.key].fill(x, y, w)
                histograms["inclusive"].fill(x, y, w)
                if histogram_groups is not None and event_svj_group_key is not None:
                    group_histograms = histogram_groups.setdefault(
                        event_svj_group_key,
                        make_empty_histograms(x_edges=histograms["inclusive"].x_edges, y_edges=histograms["inclusive"].y_edges),
                    )
                    group_histograms[pt_bin.key].fill(x, y, w)
                    group_histograms["inclusive"].fill(x, y, w)

            if args.max_events > 0 and entries_seen >= args.max_events:
                break

    return entries_seen


def make_empty_histograms(x_edges, y_edges) -> Dict[str, Hist2D]:
    keys = ["inclusive"] + [pt_bin.key for pt_bin in WNAE_PT_BINS]
    return {key: Hist2D(x_edges, y_edges) for key in keys}


def save_cache(
    path: str,
    histograms: Dict[str, Hist2D],
    histogram_groups: Dict[str, Dict[str, Hist2D]],
    pass_fail: Dict[str, PassFailAccumulator],
    event_svj: EventSVJAccumulator,
    event_svj_groups: Dict[str, EventSVJAccumulator],
    metadata: dict,
) -> None:
    import numpy as np
    import uproot

    if path.endswith(".npz"):
        payload = {"metadata": np.array(json.dumps(metadata, indent=2))}
        for key, hist in histograms.items():
            payload[f"{key}_sumw"] = hist.sumw
            payload[f"{key}_sumw2"] = hist.sumw2
            payload[f"{key}_x_edges"] = hist.x_edges
            payload[f"{key}_y_edges"] = hist.y_edges
            payload[f"{key}_entries_total"] = np.array(hist.entries_total)
            payload[f"{key}_entries_in_range"] = np.array(hist.entries_in_range)
        for group_key, group_histograms in histogram_groups.items():
            safe_group = group_key.replace("-", "_").replace(".", "p")
            for key, hist in group_histograms.items():
                payload[f"hist_group_{safe_group}_{key}_sumw"] = hist.sumw
                payload[f"hist_group_{safe_group}_{key}_sumw2"] = hist.sumw2
                payload[f"hist_group_{safe_group}_{key}_x_edges"] = hist.x_edges
                payload[f"hist_group_{safe_group}_{key}_y_edges"] = hist.y_edges
                payload[f"hist_group_{safe_group}_{key}_entries_total"] = np.array(hist.entries_total)
                payload[f"hist_group_{safe_group}_{key}_entries_in_range"] = np.array(hist.entries_in_range)
        for key, matrix in pass_fail.items():
            payload[f"pass_fail_{key}_sumw"] = matrix.sumw
            payload[f"pass_fail_{key}_sumw2"] = matrix.sumw2
            payload[f"pass_fail_{key}_entries"] = matrix.entries
        payload["event_svj_sumw"] = event_svj.sumw
        payload["event_svj_sumw2"] = event_svj.sumw2
        payload["event_svj_entries"] = event_svj.entries
        for group_key, group_matrix in event_svj_groups.items():
            safe_key = group_key.replace("-", "_").replace(".", "p")
            payload[f"event_svj_group_{safe_key}_sumw"] = group_matrix.sumw
            payload[f"event_svj_group_{safe_key}_sumw2"] = group_matrix.sumw2
            payload[f"event_svj_group_{safe_key}_entries"] = group_matrix.entries
        np.savez_compressed(path, **payload)
        return

    with uproot.recreate(path) as root_file:
        root_file["metadata"] = {"json": np.array([json.dumps(metadata, indent=2)])}
        for key, hist in histograms.items():
            root_file[f"hist_{key}"] = {
                "sumw": hist.sumw.reshape(1, -1),
                "sumw2": hist.sumw2.reshape(1, -1),
                "x_edges": hist.x_edges.reshape(1, -1),
                "y_edges": hist.y_edges.reshape(1, -1),
                "shape": np.array([[hist.sumw.shape[0], hist.sumw.shape[1]]], dtype=np.int32),
                "entries_total": np.array([hist.entries_total], dtype=np.int64),
                "entries_in_range": np.array([hist.entries_in_range], dtype=np.int64),
            }
        for group_key, group_histograms in histogram_groups.items():
            for key, hist in group_histograms.items():
                root_file[f"hist_group_{group_key}_{key}"] = {
                    "sumw": hist.sumw.reshape(1, -1),
                    "sumw2": hist.sumw2.reshape(1, -1),
                    "x_edges": hist.x_edges.reshape(1, -1),
                    "y_edges": hist.y_edges.reshape(1, -1),
                    "shape": np.array([[hist.sumw.shape[0], hist.sumw.shape[1]]], dtype=np.int32),
                    "entries_total": np.array([hist.entries_total], dtype=np.int64),
                    "entries_in_range": np.array([hist.entries_in_range], dtype=np.int64),
                }
        for key, matrix in pass_fail.items():
            root_file[f"pass_fail_{key}"] = {
                "sumw": matrix.sumw.reshape(1, -1),
                "sumw2": matrix.sumw2.reshape(1, -1),
                "entries": matrix.entries.reshape(1, -1),
            }
        root_file["event_svj"] = {
            "sumw": event_svj.sumw.reshape(1, -1),
            "sumw2": event_svj.sumw2.reshape(1, -1),
            "entries": event_svj.entries.reshape(1, -1),
        }
        for group_key, group_matrix in event_svj_groups.items():
            root_file[f"event_svj_group_{group_key}"] = {
                "sumw": group_matrix.sumw.reshape(1, -1),
                "sumw2": group_matrix.sumw2.reshape(1, -1),
                "entries": group_matrix.entries.reshape(1, -1),
            }


def validate_raw_cache(metadata: dict, histograms: Dict[str, Hist2D], path: str) -> None:
    if metadata.get("wnae_transform") == "quantile":
        raise RuntimeError(
            f"{path} was built with the old quantile WNAE axis. "
            "Rebuild the cache with the current raw-only script before using --plot-only."
        )

    inclusive = histograms.get("inclusive")
    if inclusive is not None and len(inclusive.y_edges) > 1 and inclusive.y_edges[-1] <= 1.01:
        raise RuntimeError(
            f"{path} has a 0-1 WNAE y-axis, so it is not a raw WNAE Score cache. "
            "Rebuild the cache without the old quantile option."
        )


def load_cache(
    path: str,
) -> Tuple[
    Dict[str, Hist2D],
    Dict[str, Dict[str, Hist2D]],
    Dict[str, PassFailAccumulator],
    EventSVJAccumulator,
    Dict[str, EventSVJAccumulator],
    dict,
]:
    import numpy as np

    if not path.endswith(".npz"):
        import uproot

        with uproot.open(path) as root_file:
            metadata = json.loads(root_file["metadata"].arrays(library="np")["json"][0])
            histograms = {}
            for key in metadata["histogram_keys"]:
                arrays = root_file[f"hist_{key}"].arrays(library="np")
                shape = tuple(int(x) for x in arrays["shape"][0])
                hist = Hist2D(arrays["x_edges"][0], arrays["y_edges"][0])
                hist.sumw = arrays["sumw"][0].reshape(shape)
                hist.sumw2 = arrays["sumw2"][0].reshape(shape)
                hist.entries_total = int(arrays["entries_total"][0])
                hist.entries_in_range = int(arrays["entries_in_range"][0])
                histograms[key] = hist

            histogram_groups = {}
            for group_key in metadata.get("score_histogram_group_keys", []):
                group_histograms = {}
                for key in metadata["histogram_keys"]:
                    tree_name = f"hist_group_{group_key}_{key}"
                    if tree_name not in root_file:
                        continue
                    arrays = root_file[tree_name].arrays(library="np")
                    shape = tuple(int(x) for x in arrays["shape"][0])
                    hist = Hist2D(arrays["x_edges"][0], arrays["y_edges"][0])
                    hist.sumw = arrays["sumw"][0].reshape(shape)
                    hist.sumw2 = arrays["sumw2"][0].reshape(shape)
                    hist.entries_total = int(arrays["entries_total"][0])
                    hist.entries_in_range = int(arrays["entries_in_range"][0])
                    group_histograms[key] = hist
                if group_histograms:
                    histogram_groups[group_key] = group_histograms

            pass_fail = {}
            for key in metadata.get("pass_fail_keys", ["inclusive"]):
                matrix = PassFailAccumulator()
                tree_name = f"pass_fail_{key}"
                if tree_name in root_file:
                    arrays = root_file[tree_name].arrays(library="np")
                    matrix.sumw = arrays["sumw"][0].reshape(2, 2)
                    matrix.sumw2 = arrays["sumw2"][0].reshape(2, 2)
                    matrix.entries = arrays["entries"][0].reshape(2, 2)
                pass_fail[key] = matrix

            event_svj = EventSVJAccumulator()
            if "event_svj" in root_file:
                arrays = root_file["event_svj"].arrays(library="np")
                event_svj.sumw = arrays["sumw"][0].reshape(4, 4)
                event_svj.sumw2 = arrays["sumw2"][0].reshape(4, 4)
                event_svj.entries = arrays["entries"][0].reshape(4, 4)

            event_svj_groups = {}
            for group_key in metadata.get("event_svj_group_keys", []):
                tree_name = f"event_svj_group_{group_key}"
                if tree_name not in root_file:
                    continue
                arrays = root_file[tree_name].arrays(library="np")
                matrix = EventSVJAccumulator()
                matrix.sumw = arrays["sumw"][0].reshape(4, 4)
                matrix.sumw2 = arrays["sumw2"][0].reshape(4, 4)
                matrix.entries = arrays["entries"][0].reshape(4, 4)
                event_svj_groups[group_key] = matrix

        validate_raw_cache(metadata, histograms, path)
        return histograms, histogram_groups, pass_fail, event_svj, event_svj_groups, metadata

    loaded = np.load(path, allow_pickle=False)
    metadata = json.loads(loaded["metadata"].item())
    histograms = {}
    for key in metadata["histogram_keys"]:
        hist = Hist2D(loaded[f"{key}_x_edges"], loaded[f"{key}_y_edges"])
        hist.sumw = loaded[f"{key}_sumw"]
        hist.sumw2 = loaded[f"{key}_sumw2"]
        hist.entries_total = int(loaded[f"{key}_entries_total"])
        hist.entries_in_range = int(loaded[f"{key}_entries_in_range"])
        histograms[key] = hist

    histogram_groups = {}
    for group_key in metadata.get("score_histogram_group_keys", []):
        safe_group = group_key.replace("-", "_").replace(".", "p")
        group_histograms = {}
        for key in metadata["histogram_keys"]:
            prefix = f"hist_group_{safe_group}_{key}"
            if f"{prefix}_sumw" not in loaded:
                continue
            hist = Hist2D(loaded[f"{prefix}_x_edges"], loaded[f"{prefix}_y_edges"])
            hist.sumw = loaded[f"{prefix}_sumw"]
            hist.sumw2 = loaded[f"{prefix}_sumw2"]
            hist.entries_total = int(loaded[f"{prefix}_entries_total"])
            hist.entries_in_range = int(loaded[f"{prefix}_entries_in_range"])
            group_histograms[key] = hist
        if group_histograms:
            histogram_groups[group_key] = group_histograms

    pass_fail = {}
    for key in metadata.get("pass_fail_keys", ["inclusive"]):
        matrix = PassFailAccumulator()
        if f"pass_fail_{key}_sumw" in loaded:
            matrix.sumw = loaded[f"pass_fail_{key}_sumw"]
            matrix.sumw2 = loaded[f"pass_fail_{key}_sumw2"]
            matrix.entries = loaded[f"pass_fail_{key}_entries"]
        pass_fail[key] = matrix

    event_svj = EventSVJAccumulator()
    if "event_svj_sumw" in loaded:
        event_svj.sumw = loaded["event_svj_sumw"]
        event_svj.sumw2 = loaded["event_svj_sumw2"]
        event_svj.entries = loaded["event_svj_entries"]
    event_svj_groups = {}
    for group_key in metadata.get("event_svj_group_keys", []):
        safe_key = group_key.replace("-", "_").replace(".", "p")
        sumw_key = f"event_svj_group_{safe_key}_sumw"
        if sumw_key not in loaded:
            continue
        matrix = EventSVJAccumulator()
        matrix.sumw = loaded[sumw_key]
        matrix.sumw2 = loaded[f"event_svj_group_{safe_key}_sumw2"]
        matrix.entries = loaded[f"event_svj_group_{safe_key}_entries"]
        event_svj_groups[group_key] = matrix
    validate_raw_cache(metadata, histograms, path)
    return histograms, histogram_groups, pass_fail, event_svj, event_svj_groups, metadata


def normalized_z(hist: Hist2D, mode: str):
    import numpy as np

    z = hist.sumw.astype("float64").copy()
    if mode == "none":
        return z
    if mode == "area":
        total = np.sum(z)
        return z / total if total > 0 else z
    if mode == "x":
        denom = np.sum(z, axis=1, keepdims=True)
        return np.divide(z, denom, out=np.zeros_like(z), where=denom > 0)
    if mode == "y":
        denom = np.sum(z, axis=0, keepdims=True)
        return np.divide(z, denom, out=np.zeros_like(z), where=denom > 0)
    raise ValueError(f"Unknown normalization mode: {mode}")


def quantile_value_from_projection(y_edges, y_projection, raw_value: float) -> float:
    import numpy as np

    y_edges = np.asarray(y_edges, dtype="float64")
    y_projection = np.asarray(y_projection, dtype="float64")
    total = float(np.sum(y_projection))
    if total <= 0.0:
        return 0.0
    if raw_value <= y_edges[0]:
        return 0.0
    if raw_value >= y_edges[-1]:
        return 1.0

    bin_index = int(np.searchsorted(y_edges, raw_value, side="right") - 1)
    bin_index = max(0, min(bin_index, len(y_projection) - 1))
    bin_width = y_edges[bin_index + 1] - y_edges[bin_index]
    fraction_in_bin = 0.0 if bin_width <= 0.0 else (raw_value - y_edges[bin_index]) / bin_width
    cumulative = float(np.sum(y_projection[:bin_index])) + fraction_in_bin * float(y_projection[bin_index])
    return max(0.0, min(1.0, cumulative / total))


def quantile_histogram(hist: Hist2D, y_bins: int) -> Hist2D:
    import numpy as np

    y_projection = np.sum(hist.sumw, axis=0)
    y_edges_quantile = np.linspace(0.0, 1.0, y_bins + 1)
    transformed = Hist2D(hist.x_edges, y_edges_quantile)
    transformed.entries_total = hist.entries_total
    transformed.entries_in_range = hist.entries_in_range

    if float(np.sum(y_projection)) <= 0.0:
        return transformed

    cumulative = np.concatenate(([0.0], np.cumsum(y_projection))) / float(np.sum(y_projection))
    for y_index in range(hist.sumw.shape[1]):
        if y_projection[y_index] <= 0.0:
            continue
        q_low = cumulative[y_index]
        q_high = cumulative[y_index + 1]
        q_width = q_high - q_low
        if q_width <= 0.0:
            continue
        first_target = max(0, int(np.searchsorted(y_edges_quantile, q_low, side="right") - 1))
        last_target = min(
            transformed.sumw.shape[1] - 1,
            int(np.searchsorted(y_edges_quantile, q_high, side="left")),
        )
        for target_y in range(first_target, last_target + 1):
            overlap = max(
                0.0,
                min(q_high, y_edges_quantile[target_y + 1]) - max(q_low, y_edges_quantile[target_y]),
            )
            if overlap <= 0.0:
                continue
            fraction = overlap / q_width
            transformed.sumw[:, target_y] += hist.sumw[:, y_index] * fraction
            transformed.sumw2[:, target_y] += hist.sumw2[:, y_index] * fraction * fraction

    return transformed


def quantile_heatmap_inputs(histograms: Dict[str, Hist2D], metadata: dict, y_bins: int):
    transformed = {}
    quantile_wp_values = {}
    wp_kind_default = "data" if metadata.get("all_data", False) else "mc"

    for key, hist in histograms.items():
        transformed[key] = quantile_histogram(hist, y_bins)
        pt_bin = pt_bin_for_key(key)
        if pt_bin is None:
            continue

        y_projection = hist.sumw.sum(axis=0)
        raw_wp_values = metadata.get("wnae_wp_plot_values", {}).get(key, {})
        raw_wp = raw_wp_values.get(wp_kind_default)
        if raw_wp is None:
            raw_wp = pt_bin.data_wp if metadata.get("all_data", False) else pt_bin.mc_wp
        quantile_wp_values.setdefault(key, {})[wp_kind_default] = quantile_value_from_projection(
            hist.y_edges,
            y_projection,
            raw_wp,
        )

    quantile_metadata = json.loads(json.dumps(metadata))
    quantile_metadata["wnae_transform"] = "quantile_plot"
    quantile_metadata["wnae_axis_label"] = "WNAE Score quantile"
    quantile_metadata["wnae_wp_plot_values"] = quantile_wp_values
    return transformed, quantile_metadata


def pt_bin_for_key(key: str) -> Optional[WNAEPtBin]:
    for pt_bin in WNAE_PT_BINS:
        if pt_bin.key == key:
            return pt_bin
    return None


def plot_heatmap(key: str, hist: Hist2D, metadata: dict, args, page_label: str = "") -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.colors import LogNorm

    try:
        import mplhep as hep

        plt.style.use(hep.style.CMS)
        has_mplhep = True
    except Exception:
        has_mplhep = False

    z = normalized_z(hist, args.normalize)
    x_edges = hist.x_edges
    y_edges = hist.y_edges
    is_quantile_axis = metadata.get("wnae_transform") == "quantile_plot"
    plot_ymin = y_edges[0] if args.plot_ymin is None else args.plot_ymin
    if args.plot_ymax is None:
        plot_ymax = y_edges[-1] if is_quantile_axis else 100.0
    else:
        plot_ymax = args.plot_ymax

    fig, ax = plt.subplots(figsize=(8.5, 7.2))

    plot_z = z.T
    positive = plot_z[plot_z > 0]
    norm = None
    if args.log_z and positive.size:
        norm = LogNorm(vmin=max(float(np.min(positive)), 1e-12), vmax=float(np.max(positive)))

    mesh = ax.pcolormesh(
        x_edges,
        y_edges,
        plot_z,
        cmap=args.cmap,
        shading="auto",
        norm=norm,
        edgecolors="none",
        linewidth=0.0,
        antialiased=False,
        rasterized=True,
    )
    ax.grid(False)
    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label(colorbar_label(args.normalize))
    cbar.ax.tick_params(labelsize=12)
    cbar.ax.yaxis.label.set_size(14)
    wp_color = "black" if args.log_z else "white"
    wnae_wp_for_label = None

    pt_bin = pt_bin_for_key(key)
    if pt_bin is not None:
        wp_kind = "data" if metadata.get("all_data", False) else "mc"
        wnae_wp = metadata.get("wnae_wp_plot_values", {}).get(key, {}).get(wp_kind)
        if wnae_wp is None:
            wnae_wp = pt_bin.data_wp if metadata.get("all_data", False) else pt_bin.mc_wp
        ax.axhline(wnae_wp, color=wp_color, linewidth=1.8, linestyle="-", alpha=0.95)
        wnae_wp_for_label = wnae_wp

    if args.pnet_wp is not None:
        ax.axvline(args.pnet_wp, color=wp_color, linewidth=1.8, linestyle="-", alpha=0.95)
    add_wp_label_box(
        ax,
        pnet_wp=args.pnet_wp,
        wnae_wp=wnae_wp_for_label,
        is_quantile_axis=is_quantile_axis,
        x_edges=x_edges,
        plot_ymin=plot_ymin,
        plot_ymax=plot_ymax,
    )

    ax.set_xlabel("ParticleNet SVJ tagger score", fontsize=16)
    ax.set_ylabel(metadata.get("wnae_axis_label", "WNAE Score"), fontsize=16)
    ax.tick_params(axis="both", which="major", labelsize=13)
    ax.set_xlim(x_edges[0], x_edges[-1])
    ax.set_ylim(plot_ymin, plot_ymax)

    add_score_info_box(ax, hist, metadata, page_label=page_label)
    add_cms_label(ax, metadata, has_mplhep, page_label=page_label)

    apply_figure_layout(fig, metadata)
    for ext in ("png", "pdf"):
        out_path = os.path.join(args.out_dir, f"WNAE_PNET_comparison_{key}.{ext}")
        fig.savefig(out_path, dpi=180 if ext == "png" else None, bbox_inches="tight")
        print(f"  -> {out_path}")
    plt.close(fig)


def colorbar_label(normalize: str) -> str:
    if normalize == "none":
        return "Events"
    if normalize == "area":
        return "Fraction of jets / bin"
    if normalize == "x":
        return "Fraction in PNET-score slice"
    if normalize == "y":
        return "Fraction in WNAE Score slice"
    return normalize


def write_summary(path: str, histograms: Dict[str, Hist2D]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("pt_bin,weighted_integral,entries_total,entries_in_range\n")
        for key in ["inclusive"] + [pt_bin.key for pt_bin in WNAE_PT_BINS]:
            hist = histograms[key]
            handle.write(
                f"{key},{hist.integral():.10g},{hist.entries_total},{hist.entries_in_range}\n"
            )


def plot_heatmap_set(
    histograms: Dict[str, Hist2D],
    metadata: dict,
    args,
    page_label: str = "",
) -> None:
    for key in ["inclusive"] + [pt_bin.key for pt_bin in WNAE_PT_BINS]:
        plot_heatmap(key, histograms[key], metadata, args, page_label=page_label)


def rinv_sort_value(group_key: str) -> float:
    match = re.search(r"rinv-([0-9p.]+)", group_key)
    if not match:
        return 999.0
    try:
        return float(match.group(1).replace("p", "."))
    except ValueError:
        return 999.0


def rinv_curve_label(group_key: str) -> str:
    match = re.search(r"rinv-([0-9p.]+)", group_key)
    if not match:
        return group_key
    value = match.group(1).replace("p", ".")
    return r"$r_{\mathrm{inv}} = " + value + r"$"


def mean_wnae_profile(hist: Hist2D):
    import numpy as np

    weights = np.asarray(hist.sumw, dtype="float64")
    weights = np.where(np.isfinite(weights) & (weights > 0.0), weights, 0.0)
    x_centers = 0.5 * (hist.x_edges[:-1] + hist.x_edges[1:])
    y_centers = 0.5 * (hist.y_edges[:-1] + hist.y_edges[1:])
    denom = np.sum(weights, axis=1)
    numer = np.sum(weights * y_centers[np.newaxis, :], axis=1)
    mean = np.divide(numer, denom, out=np.full_like(denom, np.nan), where=denom > 0.0)
    return x_centers, mean, denom


def plot_mean_score_overlay_by_rinv(
    histogram_groups: Dict[str, Dict[str, Hist2D]],
    metadata: dict,
    args,
    *,
    out_dir: str,
    quantile: bool = False,
) -> None:
    if len(histogram_groups) < 2:
        return

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    try:
        import mplhep as hep

        plt.style.use(hep.style.CMS)
        has_mplhep = True
    except Exception:
        has_mplhep = False

    os.makedirs(out_dir, exist_ok=True)
    plot_groups = histogram_groups
    plot_metadata = metadata
    if quantile:
        plot_groups = {}
        for group_key, group_histograms in histogram_groups.items():
            q_histograms, q_metadata = quantile_heatmap_inputs(
                group_histograms,
                metadata,
                args.quantile_wnae_bins,
            )
            plot_groups[group_key] = q_histograms
            plot_metadata = q_metadata

    keys = ["inclusive"] + [pt_bin.key for pt_bin in WNAE_PT_BINS]
    sorted_groups = sorted(plot_groups, key=rinv_sort_value)
    colors = plt.cm.coolwarm(np.linspace(0.08, 0.92, len(sorted_groups)))

    for key in keys:
        fig, ax = plt.subplots(figsize=(8.5, 7.2))
        profile_means = []
        for color, group_key in zip(colors, sorted_groups):
            hist = plot_groups[group_key].get(key)
            if hist is None:
                continue
            x_centers, mean, denom = mean_wnae_profile(hist)
            mask = np.isfinite(mean) & (denom > 0.0)
            if not np.any(mask):
                continue
            profile_means.append(mean[mask])
            ax.plot(
                x_centers[mask],
                mean[mask],
                marker="o",
                markersize=4.2,
                linewidth=2.0,
                color=color,
                label=rinv_curve_label(group_key),
            )

        if args.pnet_wp is not None:
            ax.axvline(args.pnet_wp, color="black", linewidth=1.6, linestyle="--", alpha=0.85)

        pt_bin = pt_bin_for_key(key)
        if pt_bin is not None:
            wp_kind = "data" if plot_metadata.get("all_data", False) else "mc"
            wnae_wp = plot_metadata.get("wnae_wp_plot_values", {}).get(key, {}).get(wp_kind)
            if wnae_wp is None:
                wnae_wp = pt_bin.data_wp if plot_metadata.get("all_data", False) else pt_bin.mc_wp
            ax.axhline(wnae_wp, color="black", linewidth=1.6, linestyle="--", alpha=0.85)

        ax.set_xlim(args.pnet_min, args.pnet_max)
        if quantile:
            ax.set_ylim(0.0, 1.0)
            ax.set_ylabel("Mean WNAE score quantile", fontsize=16)
        else:
            ymax = 1.18 * max([float(np.nanmax(mean)) for mean in profile_means] or [1.0])
            ax.set_ylim(0.0, max(ymax, 1.0))
            ax.set_ylabel("Mean WNAE score", fontsize=16)
        ax.set_xlabel("ParticleNet SVJ tagger score", fontsize=16)
        ax.tick_params(axis="both", which="major", labelsize=13)
        ax.grid(False)
        legend = ax.legend(
            loc="lower right",
            frameon=True,
            facecolor="white",
            edgecolor="none",
            framealpha=0.9,
            fontsize=11,
            title="Signal",
            title_fontsize=12,
            ncol=2,
            columnspacing=1.0,
            handlelength=1.6,
        )
        if legend:
            legend.set_zorder(20)

        signal_label = signal_page_label("", plot_metadata)
        if signal_label:
            ax.text(
                0.025,
                0.965,
                signal_label,
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=12,
                color="black",
                bbox={
                    "boxstyle": "square,pad=0.25",
                    "facecolor": "white",
                    "edgecolor": "none",
                    "alpha": 0.82,
                },
                zorder=10,
            )
        ax.text(
            0.0,
            1.018,
            "CMS",
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=20,
            fontweight="bold",
        )
        ax.text(
            0.13,
            1.018,
            "Simulation",
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=14,
            style="italic",
        )
        ax.text(
            1.0,
            1.018,
            "(13 TeV)",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=16,
        )
        apply_figure_layout(fig, plot_metadata)

        suffix = "_quantile" if quantile else ""
        for ext in ("png", "pdf"):
            out_path = os.path.join(out_dir, f"WNAE_PNET_mean_wnae_score_by_rinv_{key}{suffix}.{ext}")
            fig.savefig(out_path, dpi=180 if ext == "png" else None, bbox_inches="tight")
            print(f"  -> {out_path}")
        plt.close(fig)


def plot_score_histograms_by_rinv(
    histogram_groups: Dict[str, Dict[str, Hist2D]],
    metadata: dict,
    args,
) -> None:
    if not histogram_groups:
        return

    original_out_dir = args.out_dir
    original_plot_ymin = args.plot_ymin
    original_plot_ymax = args.plot_ymax
    by_rinv_dir = os.path.join(original_out_dir, "by_rinv")
    print(f"Making score heatmaps by rinv in {by_rinv_dir}")
    for group_key in sorted(histogram_groups):
        args.out_dir = os.path.join(by_rinv_dir, filename_safe(group_key))
        os.makedirs(args.out_dir, exist_ok=True)
        plot_heatmap_set(histogram_groups[group_key], metadata, args, page_label=group_key)

    if not args.no_quantile_plots and args.quantile_subdir:
        quantile_by_rinv_dir = os.path.join(original_out_dir, args.quantile_subdir, "by_rinv")
        args.plot_ymin = 0.0
        args.plot_ymax = 1.0
        print(f"Making quantile score heatmaps by rinv in {quantile_by_rinv_dir}")
        for group_key in sorted(histogram_groups):
            args.out_dir = os.path.join(quantile_by_rinv_dir, filename_safe(group_key))
            os.makedirs(args.out_dir, exist_ok=True)
            quantile_histograms, quantile_metadata = quantile_heatmap_inputs(
                histogram_groups[group_key],
                metadata,
                args.quantile_wnae_bins,
            )
            plot_heatmap_set(quantile_histograms, quantile_metadata, args, page_label=group_key)

    mean_score_dir = os.path.join(original_out_dir, "mean_score_by_rinv")
    print(f"Making mean-score overlays by rinv in {mean_score_dir}")
    plot_mean_score_overlay_by_rinv(
        histogram_groups,
        metadata,
        args,
        out_dir=mean_score_dir,
        quantile=False,
    )

    if not args.no_quantile_plots and args.quantile_subdir:
        quantile_mean_score_dir = os.path.join(
            original_out_dir, args.quantile_subdir, "mean_score_by_rinv"
        )
        print(f"Making quantile mean-score overlays by rinv in {quantile_mean_score_dir}")
        plot_mean_score_overlay_by_rinv(
            histogram_groups,
            metadata,
            args,
            out_dir=quantile_mean_score_dir,
            quantile=True,
        )

    args.out_dir = original_out_dir
    args.plot_ymin = original_plot_ymin
    args.plot_ymax = original_plot_ymax


def write_pass_fail_summary(path: str, pass_fail: Dict[str, PassFailAccumulator]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("scope,wnae,pnet,weighted_count,percent,raw_entries\n")
        for scope, matrix in pass_fail.items():
            total = float(matrix.sumw.sum())
            for y_index, wnae_label in enumerate(("Fail", "Pass")):
                for x_index, pnet_label in enumerate(("Fail", "Pass")):
                    value = float(matrix.sumw[y_index, x_index])
                    percent = 100.0 * value / total if total else 0.0
                    entries = int(matrix.entries[y_index, x_index])
                    handle.write(
                        f"{scope},{wnae_label},{pnet_label},"
                        f"{value:.10g},{percent:.6g},{entries}\n"
                    )


def plot_pass_fail(
    scope: str,
    matrix: PassFailAccumulator,
    metadata: dict,
    args,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    import numpy as np

    try:
        import mplhep as hep

        plt.style.use(hep.style.CMS)
        has_mplhep = True
    except Exception:
        has_mplhep = False

    values = matrix.sumw.astype("float64")
    total = float(values.sum())
    percentages = 100.0 * values / total if total else np.zeros_like(values)
    plot_values = values.copy()

    fig, ax = plt.subplots(figsize=(7.2, 6.4))
    positive = plot_values[plot_values > 0.0]
    norm = None
    vmin = 0.0
    vmax = float(np.max(np.abs(plot_values))) if plot_values.size else 1.0
    if vmax <= 0.0:
        vmax = 1.0
    if args.log_matrix_z and positive.size:
        vmin = float(np.min(positive))
        norm = LogNorm(vmin=vmin, vmax=vmax)
        plot_values = np.ma.masked_less_equal(plot_values, 0.0)
    imshow_kwargs = {"origin": "lower", "cmap": args.pass_fail_cmap}
    if norm is None:
        imshow_kwargs.update({"vmin": vmin, "vmax": vmax})
    else:
        imshow_kwargs["norm"] = norm
    image = ax.imshow(plot_values, **imshow_kwargs)
    cbar = fig.colorbar(image, ax=ax, pad=0.03, fraction=0.046)
    cbar.set_label("")

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Fail", "Pass"])
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["Fail", "Pass"])
    ax.set_xlabel("PNET", fontsize=16)
    ax.set_ylabel("WNAE", fontsize=16)
    ax.tick_params(axis="both", which="major", labelsize=14)
    ax.set_aspect("equal")

    for y in (0, 1):
        for x in (0, 1):
            value = values[y, x]
            percent = percentages[y, x]
            text = f"{formatted_cell_value(value, metadata)}\n({percent:.2f}%)"
            ax.text(
                x,
                y,
                text,
                ha="center",
                va="center",
                fontsize=15,
                color="white",
                fontweight="bold",
            )

    ax.grid(False)
    ax.axvline(0.5, color="black", linewidth=1.4)
    ax.axhline(0.5, color="black", linewidth=1.4)
    for spine in ax.spines.values():
        spine.set_linewidth(1.6)

    add_cms_label(ax, metadata, has_mplhep)

    apply_figure_layout(fig, metadata)
    suffix = "pass_fail" if scope == "inclusive" else f"pass_fail_{scope}"
    for ext in ("png", "pdf"):
        out_path = os.path.join(args.out_dir, f"WNAE_PNET_comparison_{suffix}.{ext}")
        fig.savefig(out_path, dpi=180 if ext == "png" else None, bbox_inches="tight")
        print(f"  -> {out_path}")
    plt.close(fig)


def write_event_svj_summary(path: str, event_svj: EventSVJAccumulator) -> None:
    total = float(event_svj.sumw.sum())
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("wnae_nsvj,pnet_nsvj,weighted_count,percent,raw_entries\n")
        for y, wnae_label in enumerate(SVJ_CATEGORY_LABELS):
            for x, pnet_label in enumerate(SVJ_CATEGORY_LABELS):
                value = float(event_svj.sumw[y, x])
                percent = 100.0 * value / total if total else 0.0
                entries = int(event_svj.entries[y, x])
                handle.write(
                    f"{wnae_label},{pnet_label},{value:.10g},{percent:.6g},{entries}\n"
                )


def write_event_svj_group_summary(
    path: str,
    event_svj: EventSVJAccumulator,
    event_svj_groups: Dict[str, EventSVJAccumulator],
) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("group,wnae_nsvj,pnet_nsvj,weighted_count,percent,raw_entries\n")
        groups = [("inclusive", event_svj)] + [
            (group_key, event_svj_groups[group_key]) for group_key in sorted(event_svj_groups)
        ]
        for group_key, matrix in groups:
            total = float(matrix.sumw.sum())
            for y, wnae_label in enumerate(SVJ_CATEGORY_LABELS):
                for x, pnet_label in enumerate(SVJ_CATEGORY_LABELS):
                    value = float(matrix.sumw[y, x])
                    percent = 100.0 * value / total if total else 0.0
                    entries = int(matrix.entries[y, x])
                    handle.write(
                        f"{group_key},{wnae_label},{pnet_label},"
                        f"{value:.10g},{percent:.6g},{entries}\n"
                    )


def draw_event_svj_matrix(ax, fig, event_svj: EventSVJAccumulator, metadata: dict, args, page_label: str = "") -> None:
    from matplotlib.colors import LogNorm
    import matplotlib.patheffects as pe
    import numpy as np

    values = event_svj.sumw.astype("float64")
    total = float(values.sum())
    percentages = 100.0 * values / total if total else np.zeros_like(values)
    plot_values = values.copy()
    positive = plot_values[plot_values > 0.0]
    norm = None
    vmin = 0.0
    vmax = float(np.max(plot_values)) if plot_values.size else 1.0
    if vmax <= 0.0:
        vmax = 1.0
    if args.log_matrix_z and positive.size:
        vmin = float(np.min(positive))
        norm = LogNorm(vmin=vmin, vmax=vmax)
        plot_values = np.ma.masked_less_equal(plot_values, 0.0)

    imshow_kwargs = {"origin": "lower", "cmap": args.cmap}
    if norm is None:
        imshow_kwargs.update({"vmin": vmin, "vmax": vmax})
    else:
        imshow_kwargs["norm"] = norm
    image = ax.imshow(plot_values, **imshow_kwargs)
    cbar = fig.colorbar(image, ax=ax, pad=0.03, fraction=0.046)
    cbar.set_label("Events", fontsize=22)
    cbar.ax.tick_params(labelsize=16)

    ax.set_xticks(range(4))
    ax.set_xticklabels(SVJ_CATEGORY_LABELS)
    ax.set_yticks(range(4))
    ax.set_yticklabels(SVJ_CATEGORY_LABELS)
    ax.set_xlabel(r"$n_{\mathrm{SVJ}}^{\mathrm{PN}}$", fontsize=24, labelpad=24)
    ax.set_ylabel(r"$n_{\mathrm{SVJ}}^{\mathrm{WNAE}}$", fontsize=24)
    ax.tick_params(axis="both", which="major", labelsize=18)
    ax.set_aspect("equal")

    for y in range(4):
        for x in range(4):
            value = values[y, x]
            percent = percentages[y, x]
            text = f"{formatted_cell_value(value, metadata)}\n({percent:.2f}%)"
            ax.text(
                x,
                y,
                text,
                ha="center",
                va="center",
                fontsize=13.5,
                color="white",
                fontweight="bold",
                linespacing=0.90,
                path_effects=[pe.withStroke(linewidth=2.4, foreground="black")],
            )

    for edge in np.arange(-0.5, 4.0, 1.0):
        ax.axvline(edge, color="black", linewidth=1.0)
        ax.axhline(edge, color="black", linewidth=1.0)
    for spine in ax.spines.values():
        spine.set_linewidth(1.6)


def plot_event_svj_matrix(event_svj: EventSVJAccumulator, metadata: dict, args) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    try:
        import mplhep as hep

        plt.style.use(hep.style.CMS)
        has_mplhep = True
    except Exception:
        has_mplhep = False

    fig, ax = plt.subplots(figsize=(10.2, 9.1))
    draw_event_svj_matrix(ax, fig, event_svj, metadata, args)

    add_signal_matrix_info_box(ax, metadata)
    add_cms_label(ax, metadata, has_mplhep, matrix_style=True)

    apply_matrix_layout(fig)
    for ext in ("png", "pdf"):
        out_path = os.path.join(args.out_dir, f"WNAE_PNET_comparison_event_nsvj_matrix.{ext}")
        fig.savefig(out_path, dpi=180 if ext == "png" else None)
        print(f"  -> {out_path}")
    plt.close(fig)


def plot_event_svj_matrix_multipage(
    event_svj: EventSVJAccumulator,
    event_svj_groups: Dict[str, EventSVJAccumulator],
    metadata: dict,
    args,
) -> None:
    if not event_svj_groups:
        return

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    try:
        import mplhep as hep

        plt.style.use(hep.style.CMS)
        has_mplhep = True
    except Exception:
        has_mplhep = False

    out_path = os.path.join(args.out_dir, "WNAE_PNET_comparison_event_nsvj_matrix_by_rinv.pdf")
    with PdfPages(out_path) as pdf:
        pages = [("inclusive", event_svj)] + [
            (group_key, event_svj_groups[group_key]) for group_key in sorted(event_svj_groups)
        ]
        for page_label, matrix in pages:
            fig, ax = plt.subplots(figsize=(10.2, 9.1))
            draw_event_svj_matrix(ax, fig, matrix, metadata, args, page_label=page_label)
            add_signal_matrix_info_box(ax, metadata, page_label=page_label)
            add_cms_label(ax, metadata, has_mplhep, page_label=page_label, matrix_style=True)
            apply_matrix_layout(fig)
            pdf.savefig(fig)
            if page_label != "inclusive":
                suffix = filename_safe(page_label)
                for ext in ("png", "pdf"):
                    page_out = os.path.join(
                        args.out_dir, f"WNAE_PNET_comparison_event_nsvj_matrix_{suffix}.{ext}"
                    )
                    fig.savefig(page_out, dpi=180 if ext == "png" else None)
                    print(f"  -> {page_out}")
            plt.close(fig)
    print(f"  -> {out_path}")


def wnae_wp_plot_values() -> Dict[str, Dict[str, float]]:
    values = {}
    for pt_bin in WNAE_PT_BINS:
        values[pt_bin.key] = {"mc": pt_bin.mc_wp, "data": pt_bin.data_wp}
    return values


def build_histograms(args):
    import numpy as np

    x_edges = np.linspace(args.pnet_min, args.pnet_max, args.pnet_bins + 1)
    y_edges = np.linspace(args.wnae_min, args.wnae_max, args.wnae_bins + 1)
    keys = ["inclusive"] + [pt_bin.key for pt_bin in WNAE_PT_BINS]
    histograms = {key: Hist2D(x_edges, y_edges) for key in keys}
    histogram_groups: Dict[str, Dict[str, Hist2D]] = {}
    pass_fail = {key: PassFailAccumulator() for key in keys}
    event_svj = EventSVJAccumulator()
    event_svj_groups: Dict[str, EventSVJAccumulator] = {}

    if args.file_list:
        files = files_from_text(args.file_list)
    else:
        files = discover_files_from_base(
            base=args.base,
            years=parse_csv(args.years),
            processes=parse_csv(args.processes),
            sample_dirs=parse_csv(args.sample_dirs) if args.sample_dirs else [],
            max_files_per_sample=args.max_files_per_sample,
        )

    if not files:
        raise RuntimeError("No input ROOT files were found.")
    files = [
        file_info
        for file_info in files
        if file_info[1].lower() != "signal" or sample_matches_signal_filters(file_info[2], args)
    ]
    if not files:
        raise RuntimeError("No input ROOT files survived the requested signal filters.")

    print(f"Streaming {len(files)} ROOT files")
    total_entries = 0
    failures = 0
    for index, (year, process, sample, file_url) in enumerate(files, start=1):
        print(f"[{index}/{len(files)}] {year} {process}/{sample}: {file_url}")
        event_group_key = signal_group_from_sample(sample) if process.lower() == "signal" else None
        try:
            entries = stream_file(
                file_url,
                year,
                process,
                sample,
                args,
                histograms,
                histogram_groups=histogram_groups,
                pass_fail=pass_fail,
                event_svj=event_svj,
                event_svj_groups=event_svj_groups,
                event_svj_group_key=event_group_key,
            )
            total_entries += entries
        except Exception as exc:
            failures += 1
            print(f"[WARN] skipped file: {exc}")
            if args.fail_fast:
                raise

    metadata = {
        "script": os.path.basename(__file__),
        "base": args.base,
        "pnet_friend_base": args.pnet_friend_base,
        "years": parse_csv(args.years),
        "processes": parse_csv(args.processes),
        "sample_dirs": parse_csv(args.sample_dirs) if args.sample_dirs else [],
        "tree": args.tree,
        "pt_branch": PT_BRANCH,
        "pnet_branch": PNET_BRANCH,
        "good_jet_branch": GOOD_JET_BRANCH,
        "wnae_pt_bins": [asdict(pt_bin) for pt_bin in WNAE_PT_BINS],
        "histogram_keys": keys,
        "score_histogram_group_keys": sorted(histogram_groups),
        "n_files": len(files),
        "failed_files": failures,
        "total_entries_seen": total_entries,
        "all_data": bool(files) and all(process.lower() == "data" for _, process, _, _ in files),
        "tag_region": args.tag_region,
        "pnet_wp": args.pnet_wp,
        "no_weights": args.no_weights,
        "no_lund_correction": args.no_lund_correction,
        "weight_formula": (
            "data or --no-weights: 1; MC: Weight*lumi_pb*puWeight*NonPrefiringProb "
            "when branches exist"
            + (
                "; Lund correction disabled by --no-lund-correction"
                if args.no_lund_correction
                else (
                    "; signal-like files with lundWeightNom and "
                    "CutFlow/InitialLundNominal additionally use "
                    "lundWeightNom*Initial/InitialLundNominal"
                )
            )
        ),
        "wnae_wp_plot_values": wnae_wp_plot_values(),
        "pass_fail_keys": keys,
        "svj_category_keys": SVJ_CATEGORY_KEYS,
        "svj_category_labels": SVJ_CATEGORY_LABELS,
        "event_svj_group_keys": sorted(event_svj_groups),
        "signal_mdark": args.signal_mdark,
        "signal_yukawa": args.signal_yukawa,
        "signal_mmed_label": args.signal_mmed_label,
    }
    return histograms, histogram_groups, pass_fail, event_svj, event_svj_groups, metadata


def parse_args(argv: Optional[Sequence[str]] = None):
    parser = argparse.ArgumentParser(
        description="Make WNAE Score vs PNET-score 2D heatmaps in WNAE pT bins."
    )
    parser.add_argument("--base", default=DEFAULT_BASE, help="EOS base before YEAR/t_channel_pre_selection/nominal")
    parser.add_argument("--years", default="2016,2017,2018", help="Comma-separated years")
    parser.add_argument(
        "--processes",
        default="QCD,TTJets,WJetsToLNu,ZJetsToNuNu,ST",
        help="Comma-separated process keys from the built-in manifest. Include Data for data.",
    )
    parser.add_argument(
        "--sample-dirs",
        default="",
        help="Optional comma-separated sample directories to use for every requested process.",
    )
    parser.add_argument(
        "--signal-mdark",
        default="",
        help="Optional signal sample filter, e.g. 20 keeps sample names containing mDark-20.",
    )
    parser.add_argument(
        "--signal-yukawa",
        default="",
        help="Optional signal sample filter, e.g. 1 keeps sample names containing yukawa-1.",
    )
    parser.add_argument(
        "--signal-mmed-label",
        default="",
        help="Optional label for signal mMed/mphi, e.g. 2000 or all.",
    )
    parser.add_argument(
        "--file-list",
        default="",
        help="Optional text file with FILE or YEAR,PROCESS,SAMPLE,FILE per line.",
    )
    parser.add_argument(
        "--pnet-friend-base",
        default="",
        help=(
            "Optional EOS base containing matching files with "
            "JetsAK8_pNetJetTaggerScore, used when the main input has WNAE "
            "branches but no PNET branch."
        ),
    )
    parser.add_argument("--tree", default=TREE_NAME, help="Input TTree name")
    parser.add_argument("--out-dir", default="WNAE_PNET_comparison", help="Output directory")
    parser.add_argument(
        "--cache",
        default="",
        help="ROOT cache path. Defaults to OUT_DIR/WNAE_PNET_comparison_histograms.root",
    )
    parser.add_argument("--plot-only", action="store_true", help="Read --cache and only remake plots")
    parser.add_argument("--no-plots", action="store_true", help="Build cache only")
    parser.add_argument("--chunk-size", default="100 MB", help="uproot iterate step size")
    parser.add_argument("--max-files-per-sample", type=int, default=0, help="Testing: cap files per sample; 0 means all")
    parser.add_argument("--max-events", type=int, default=0, help="Testing: cap events per input file; 0 means all")
    parser.add_argument("--fail-fast", action="store_true", help="Stop on the first bad file")
    parser.add_argument("--no-weights", action="store_true", help="Use unit event weights")
    parser.add_argument(
        "--no-lund-correction",
        action="store_true",
        help="Do not apply lundWeightNom*Initial/InitialLundNominal even if available.",
    )

    parser.add_argument("--pnet-bins", type=int, default=25)
    parser.add_argument("--pnet-min", type=float, default=0.0)
    parser.add_argument("--pnet-max", type=float, default=1.0)
    parser.add_argument("--wnae-bins", type=int, default=200)
    parser.add_argument("--wnae-min", type=float, default=0.0)
    parser.add_argument("--wnae-max", type=float, default=200.0)
    parser.add_argument(
        "--plot-ymin",
        type=float,
        default=None,
        help="Plot-only y-axis minimum override for heatmaps.",
    )
    parser.add_argument(
        "--plot-ymax",
        type=float,
        default=None,
        help="Plot-only y-axis maximum override for heatmaps.",
    )
    parser.add_argument(
        "--quantile-subdir",
        default="quantile",
        help="Subdirectory under --out-dir for additional WNAE 0-1 quantile heatmaps.",
    )
    parser.add_argument(
        "--quantile-wnae-bins",
        type=int,
        default=200,
        help="Number of WNAE quantile-axis bins for the additional 0-1 quantile heatmaps.",
    )
    parser.add_argument(
        "--no-quantile-plots",
        action="store_true",
        help="Do not write the additional WNAE 0-1 quantile heatmaps.",
    )
    parser.add_argument("--pnet-wp", type=float, default=0.90, help="PNET WP line and tag-region threshold")
    parser.add_argument(
        "--tag-region",
        choices=["all", "pnet", "wnae", "both", "either"],
        default="all",
        help="Fill all jets or only jets passing one/both tag definitions.",
    )
    parser.add_argument(
        "--normalize",
        choices=["none", "area", "x", "y"],
        default="none",
        help="Plot normalization only; cache always stores raw sumw.",
    )
    parser.add_argument("--log-z", action="store_true", default=False, help="Use log color scale")
    parser.add_argument("--linear-z", dest="log_z", action="store_false", help="Use linear color scale")
    parser.add_argument(
        "--log-matrix-z",
        action="store_true",
        default=False,
        help="Use log color scale for pass/fail and event nSVJ matrix plots.",
    )
    parser.add_argument("--cmap", default="viridis")
    parser.add_argument(
        "--pass-fail-cmap",
        default=None,
        help="Colormap for the pass/fail map. Defaults to --cmap.",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    require_runtime_packages()
    if args.pass_fail_cmap is None:
        args.pass_fail_cmap = args.cmap

    os.makedirs(args.out_dir, exist_ok=True)
    if not args.cache:
        args.cache = os.path.join(args.out_dir, "WNAE_PNET_comparison_histograms.root")

    if args.plot_only:
        histograms, histogram_groups, pass_fail, event_svj, event_svj_groups, metadata = load_cache(args.cache)
        if args.signal_mdark:
            metadata["signal_mdark"] = args.signal_mdark
        if args.signal_yukawa:
            metadata["signal_yukawa"] = args.signal_yukawa
        if args.signal_mmed_label:
            metadata["signal_mmed_label"] = args.signal_mmed_label
    else:
        histograms, histogram_groups, pass_fail, event_svj, event_svj_groups, metadata = build_histograms(args)
        save_cache(args.cache, histograms, histogram_groups, pass_fail, event_svj, event_svj_groups, metadata)
        print(f"Saved histogram cache: {args.cache}")

    summary_path = os.path.join(args.out_dir, "WNAE_PNET_comparison_summary.csv")
    write_summary(summary_path, histograms)
    print(f"Saved summary: {summary_path}")

    pass_fail_summary_path = os.path.join(args.out_dir, "WNAE_PNET_comparison_pass_fail.csv")
    write_pass_fail_summary(pass_fail_summary_path, pass_fail)
    print(f"Saved pass/fail summary: {pass_fail_summary_path}")

    event_svj_summary_path = os.path.join(args.out_dir, "WNAE_PNET_comparison_event_nsvj_matrix.csv")
    write_event_svj_summary(event_svj_summary_path, event_svj)
    print(f"Saved event nSVJ matrix summary: {event_svj_summary_path}")
    if event_svj_groups:
        event_svj_group_summary_path = os.path.join(
            args.out_dir, "WNAE_PNET_comparison_event_nsvj_matrix_by_rinv.csv"
        )
        write_event_svj_group_summary(event_svj_group_summary_path, event_svj, event_svj_groups)
        print(f"Saved event nSVJ by-rinv summary: {event_svj_group_summary_path}")

    if not args.no_plots:
        print("Making heatmaps")
        plot_heatmap_set(histograms, metadata, args)
        if not args.no_quantile_plots and args.quantile_subdir:
            quantile_dir = os.path.join(args.out_dir, args.quantile_subdir)
            os.makedirs(quantile_dir, exist_ok=True)
            quantile_histograms, quantile_metadata = quantile_heatmap_inputs(
                histograms,
                metadata,
                args.quantile_wnae_bins,
            )
            original_out_dir = args.out_dir
            original_plot_ymin = args.plot_ymin
            original_plot_ymax = args.plot_ymax
            args.out_dir = quantile_dir
            args.plot_ymin = 0.0
            args.plot_ymax = 1.0
            print(f"Making quantile heatmaps in {quantile_dir}")
            plot_heatmap_set(quantile_histograms, quantile_metadata, args)
            args.out_dir = original_out_dir
            args.plot_ymin = original_plot_ymin
            args.plot_ymax = original_plot_ymax
        plot_score_histograms_by_rinv(histogram_groups, metadata, args)
        print("Making pass/fail map")
        plot_pass_fail("inclusive", pass_fail["inclusive"], metadata, args)
        print("Making event nSVJ matrix")
        plot_event_svj_matrix(event_svj, metadata, args)
        plot_event_svj_matrix_multipage(event_svj, event_svj_groups, metadata, args)

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(1)
