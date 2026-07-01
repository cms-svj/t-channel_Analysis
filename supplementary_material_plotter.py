#!/usr/bin/env python3
"""
supplementary_material_plotter_FIXED_v2.py

LPC-ready supplementary material plotter for SVJ DNN and ParticleNet/PNet plots.

Default behavior is intentionally FAST:
  - DNN reference score distribution
  - DNN reference ROC
  - PNet reference score distribution
  - PNet reference ROC

Full signal scans are OFF by default. Turn them on only with:
  --make-auc-maps
  --make-eff-maps
  --make-auc-tables
  --make-selected-roc
  --make-eff-vs-scorecut

Important fixes relative to the old script:
  - PNet reads Events/JetsAK8_pNetJetTaggerScore from skim files.
  - PNet event score is max(JetsAK8_pNetJetTaggerScore) over AK8 jets per event.
  - PNet signal points are directories under the EOS nominal directory, not ABCD_scores files.
  - ROC axes are standard: False Positive Rate (FPR), True Positive Rate (TPR).
  - Score distributions use black for QCD/background and red for signal.
  - Score distributions have no title like "DNN score distribution".
  - DNN exact nominal file selection avoids pdfUp/pdfDown/etc.
  - AUC tables are only written with --make-auc-tables.

Quick run:
  python3 supplementary_material_plotter_FIXED_v2.py --all --out-dir supplementary_material_plots

Full maps later:
  python3 supplementary_material_plotter_FIXED_v2.py --all --make-auc-maps --make-eff-maps --out-dir supplementary_material_plots
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import uproot

try:
    import awkward as ak
except Exception:
    ak = None

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe

try:
    import mplhep as hep
    plt.style.use(hep.style.CMS)
    HAS_MPLHEP = True
except Exception:
    hep = None
    HAS_MPLHEP = False

from sklearn.metrics import roc_auc_score, roc_curve


# =============================================================================
# DNN defaults from the supplied AUC-maker
# =============================================================================

DEFAULT_DNN_SIG_DIR = (
    "/uscms/home/nparmar/nobackup/SVJ/eventLevelTagger/output/dataset11/"
    "training_set2/application_results/"
    "sdt_allBkgs_disco_0p001_closure_0p06_damp_1_net_64_32_16_8_"
    "trainingAllYears_met250_dnn85_eval2018/v_5/"
    "ML_fit_results_grid_optimization_binned/"
    "scores_signals_GapJetVeto/2018/nominal"
)

DEFAULT_DNN_BKG_FILE = (
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/"
    "tchannel_UL/skims_gapJetVeto/"
    "DNN_scores_for_background_with_PNET/"
    "scores_backgrounds_GapJetVeto/2018/nominal/"
    "ABCD_scores_QCD.root"
)


# =============================================================================
# PNet defaults from Akshat's skim path and old PNet script
# =============================================================================

DEFAULT_PNET_SIG_BASE = (
    "/store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/"
    "2018/t_channel_pre_selection/nominal"
)

PNET_TREE = "Events"
PNET_SCORE_BRANCH = "JetsAK8_pNetJetTaggerScore"

DEFAULT_PNET_QCD_BASE = (
    "/store/user/ashrivas/tchannel_UL/skims_gapJetVeto/"
    "data_mc_withWNAE_METtrim/trimmed_root"
)
DEFAULT_PNET_YEAR = "2018"


DNN_TRAINING_SAMPLES = [
    {"mMed": 600,  "rinv": 0.3},
    {"mMed": 800,  "rinv": 0.3},
    {"mMed": 1000, "rinv": 0.3},
    {"mMed": 1500, "rinv": 0.3},
    {"mMed": 2000, "rinv": 0.1},
    {"mMed": 2000, "rinv": 0.3},
    {"mMed": 2000, "rinv": 0.5},
    {"mMed": 2000, "rinv": 0.7},
    {"mMed": 3000, "rinv": 0.3},
    {"mMed": 4000, "rinv": 0.3},
]

PNET_TRAINING_SAMPLES = [
    {"mMed": 2000, "mDark": 20,  "rinv": 0.3},
    {"mMed": 600,  "mDark": 20,  "rinv": 0.3},
    {"mMed": 800,  "mDark": 20,  "rinv": 0.3},
    {"mMed": 1000, "mDark": 20,  "rinv": 0.3},
    {"mMed": 1500, "mDark": 20,  "rinv": 0.3},
    {"mMed": 3000, "mDark": 20,  "rinv": 0.3},
    {"mMed": 4000, "mDark": 20,  "rinv": 0.3},
    {"mMed": 2000, "mDark": 1,   "rinv": 0.3},
    {"mMed": 2000, "mDark": 100, "rinv": 0.3},
    {"mMed": 2000, "mDark": 20,  "rinv": 0.1},
    {"mMed": 2000, "mDark": 20,  "rinv": 0.5},
    {"mMed": 2000, "mDark": 20,  "rinv": 0.7},
]


# =============================================================================
# General utilities
# =============================================================================

def mkdir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def clean_name(s: str) -> str:
    s = str(s).strip().replace("/", "_").replace(".", "p")
    return re.sub(r"[^A-Za-z0-9_+\-]+", "_", s).strip("_")


def save_fig(fig: plt.Figure, prefix: str, formats: Sequence[str] = ("pdf", "png")) -> None:
    for ext in formats:
        path = f"{prefix}.{ext}"
        fig.savefig(path, bbox_inches="tight", dpi=180 if ext == "png" else None)
        print(f"  -> {path}")
    plt.close(fig)


def cms_label(ax: plt.Axes, fontsize: int = 24, com: Optional[int] = 13) -> None:
    if HAS_MPLHEP:
        try:
            hep.cms.label(data=False, ax=ax, loc=0, fontsize=fontsize, com=com)
        except TypeError:
            hep.cms.label(data=False, ax=ax, loc=0, fontsize=fontsize)
    else:
        ax.text(0.02, 0.98, "CMS Simulation", transform=ax.transAxes,
                ha="left", va="top", fontsize=fontsize, fontweight="bold")


def file_exists(path: str) -> bool:
    return bool(path) and os.path.isfile(path)


def dir_exists(path: str) -> bool:
    return bool(path) and os.path.isdir(path)


def eos_to_xrootd(path: str) -> str:
    """Convert /eos/uscms or /store paths to root://cmseos.fnal.gov paths."""
    if path.startswith("root://"):
        return path
    if path.startswith("/eos/uscms/"):
        return "root://cmseos.fnal.gov//" + path.replace("/eos/uscms/", "", 1)
    if path.startswith("/store/"):
        return "root://cmseos.fnal.gov/" + path
    return path


def xrootd_to_xrdfs_path(path: str) -> str:
    if path.startswith("root://cmseos.fnal.gov//"):
        return "/" + path.split("root://cmseos.fnal.gov//", 1)[1]
    if path.startswith("root://cmseos.fnal.gov/"):
        return path.split("root://cmseos.fnal.gov", 1)[1]
    if path.startswith("/eos/uscms/"):
        return "/" + path.replace("/eos/uscms/", "", 1)
    return path


def list_dir(path: str) -> List[str]:
    """List a local/EOS/xrootd directory and return full paths."""
    if os.path.isdir(path):
        return [os.path.join(path, x) for x in sorted(os.listdir(path))]

    # If user passed /store or root://, use xrdfs.
    xrdfs_path = xrootd_to_xrdfs_path(path)
    cmd = ["xrdfs", "root://cmseos.fnal.gov", "ls", xrdfs_path]
    try:
        out = subprocess.check_output(cmd, text=True, stderr=subprocess.STDOUT)
        lines = [x.strip() for x in out.splitlines() if x.strip()]
        # Use xrootd paths for uproot if local path does not exist.
        return ["root://cmseos.fnal.gov/" + x if x.startswith("/store/") else x for x in lines]
    except Exception as e:
        print(f"  [WARN] Could not list directory: {path}\n         {e}")
        return []


def root_files_in_dir(path: str) -> List[str]:
    if os.path.isdir(path):
        return sorted(glob.glob(os.path.join(path, "*.root")))
    return [x for x in list_dir(path) if x.endswith(".root")]


def qcd_sort_key(name_or_path: str) -> Tuple[int, int]:
    """Sort QCD_Pt_170to300 through QCD_Pt_3200toInf numerically."""
    name = os.path.basename(os.path.normpath(name_or_path))
    name = re.sub(r"^\d{4}_", "", name)
    match = re.match(r"QCD_Pt_(\d+)to(\d+|Inf)$", name)
    if match is None:
        return (10**9, 10**9)
    low = int(match.group(1))
    high = 10**9 if match.group(2) == "Inf" else int(match.group(2))
    return (low, high)


def discover_trimmed_qcd_files(qcd_base: str, year: str) -> List[str]:
    """Discover all ROOT files in <year>_QCD_Pt_* children of a trimmed-QCD base."""
    qcd_dirs = []
    prefix = f"{year}_QCD_Pt_"
    for entry in list_dir(qcd_base):
        name = os.path.basename(os.path.normpath(entry))
        if name.startswith(prefix):
            qcd_dirs.append(entry)

    files: List[str] = []
    for qcd_dir in sorted(qcd_dirs, key=qcd_sort_key):
        files.extend(root_files_in_dir(qcd_dir))
    return files


# =============================================================================
# ROOT loading: DNN score files and generic input branches
# =============================================================================

def get_tree(fin: uproot.ReadOnlyDirectory,
             preferred: Sequence[str] = ("ABCD_scores", "Events", "tree", "PreSelection", "TreeMaker2/PreSelection")):
    for key in preferred:
        try:
            obj = fin[key]
            if isinstance(obj, uproot.behaviors.TTree.TTree):
                return obj
        except Exception:
            pass
    for key in fin.keys():
        kclean = key.split(";")[0]
        try:
            obj = fin[kclean]
            if isinstance(obj, uproot.behaviors.TTree.TTree):
                return obj
        except Exception:
            pass
    for key in fin.keys():
        kclean = key.split(";")[0]
        try:
            obj = fin[kclean]
            if hasattr(obj, "keys"):
                for sub in obj.keys():
                    subclean = sub.split(";")[0]
                    maybe = obj[subclean]
                    if isinstance(maybe, uproot.behaviors.TTree.TTree):
                        return maybe
        except Exception:
            pass
    return None


def tree_keys(path: str) -> List[str]:
    with uproot.open(eos_to_xrootd(path) if path.startswith("/store/") else path) as fin:
        tree = get_tree(fin)
        if tree is None:
            return []
        return list(tree.keys())


def load_branch(path: str, branch: str, max_events: int = -1):
    open_path = eos_to_xrootd(path) if path.startswith("/store/") else path
    with uproot.open(open_path) as fin:
        tree = get_tree(fin)
        if tree is None:
            raise RuntimeError(f"No TTree found in {path}")
        if branch not in tree.keys():
            raise KeyError(f"Branch '{branch}' not found in {path}")
        entry_stop = None if max_events is None or max_events < 0 else max_events
        return tree[branch].array(library="ak" if ak is not None else "np", entry_stop=entry_stop)


def to_flat_numpy(arr) -> np.ndarray:
    if ak is not None and isinstance(arr, ak.Array):
        try:
            # Flatten both scalar and jagged safely.
            if "var" in str(arr.type):
                arr = ak.flatten(arr, axis=None)
        except Exception:
            pass
        arr = ak.to_numpy(arr)
    arr = np.asarray(arr)
    if arr.dtype.kind not in "biufc":
        return np.array([], dtype=float)
    arr = arr.astype(float).reshape(-1)
    return arr[np.isfinite(arr)]


def load_dnn_scores_weights(path: str, score_branch: str = "scores", weight_branch: str = "weights",
                            max_events: int = -1) -> Tuple[np.ndarray, np.ndarray]:
    if not path.startswith("root://") and not path.startswith("/store/") and not file_exists(path):
        raise FileNotFoundError(path)
    open_path = eos_to_xrootd(path) if path.startswith("/store/") else path
    with uproot.open(open_path) as fin:
        tree = get_tree(fin, preferred=("ABCD_scores", "Events", "tree"))
        if tree is None:
            raise RuntimeError(f"No TTree found in {path}")
        keys = list(tree.keys())
        if score_branch not in keys:
            raise KeyError(f"Score branch '{score_branch}' not found in {path}. Available keys include: {keys[:40]}")
        entry_stop = None if max_events is None or max_events < 0 else max_events
        scores = tree[score_branch].array(library="np", entry_stop=entry_stop).astype(float).reshape(-1)
        if weight_branch and weight_branch in keys:
            weights = tree[weight_branch].array(library="np", entry_stop=entry_stop).astype(float).reshape(-1)
        else:
            weights = np.ones_like(scores, dtype=float)
    mask = np.isfinite(scores) & np.isfinite(weights) & (weights >= 0)
    return scores[mask], weights[mask]


def discover_numeric_scalar_branches(path_a: str, path_b: str, exclude: Sequence[str], max_vars: int) -> List[str]:
    try:
        keys_a = set(tree_keys(path_a))
        keys_b = set(tree_keys(path_b))
    except Exception as e:
        print(f"  [WARN] Could not discover branches: {e}")
        return []
    common = sorted(keys_a & keys_b)
    out: List[str] = []
    exclude_set = set(exclude)
    for br in common:
        if br in exclude_set:
            continue
        if any(x.lower() in br.lower() for x in ("weight", "score", "event", "run", "lumi")):
            continue
        try:
            flat = to_flat_numpy(load_branch(path_a, br, max_events=1000))
            if len(flat) == 0:
                continue
            out.append(br)
            if len(out) >= max_vars:
                break
        except Exception:
            continue
    return out


# =============================================================================
# DNN signal file parsing
# =============================================================================

def parse_signal_params_from_name(name_or_path: str) -> Optional[Dict[str, object]]:
    name = os.path.basename(os.path.normpath(name_or_path))
    m = re.search(
        r"mMed-(\d+)_mDark-(\d+)_rinv-([0-9p.]+)_alpha-([^_]+)_yukawa-([0-9p.]+)",
        name,
    )
    if not m:
        return None
    return {
        "mMed": int(m.group(1)),
        "mDark": int(m.group(2)),
        "rinv": float(m.group(3).replace("p", ".")),
        "alpha": m.group(4),
        "yukawa": float(str(m.group(5)).replace("p", ".")),
        "path": name_or_path,
    }


def is_exact_dnn_nominal_file(path: str) -> bool:
    """Accept only the true nominal DNN score file, not pdfUp/pdfDown/etc."""
    base = os.path.basename(path)
    if not base.endswith(".root"):
        return False
    if not base.startswith("ABCD_scores_nominal_t-channel_"):
        return False
    return base.endswith("_nominal.root") or base.endswith("_lundWeightNom.root")


def find_dnn_signal_files(sig_dir: str, mdark: int, alpha: str, yukawa: float) -> List[str]:
    if not dir_exists(sig_dir):
        return []
    files = sorted(glob.glob(os.path.join(sig_dir, "*.root")))
    out: List[str] = []
    for f in files:
        if not is_exact_dnn_nominal_file(f):
            continue
        p = parse_signal_params_from_name(f)
        if p is None:
            continue
        if int(p["mDark"]) != int(mdark):
            continue
        if str(p["alpha"]) != str(alpha):
            continue
        if abs(float(p["yukawa"]) - float(yukawa)) > 1e-9:
            continue
        out.append(f)
    return out


def find_dnn_reference_signal(sig_dir: str, mmed: int, mdark: int, rinv: float, alpha: str, yukawa: float) -> Optional[str]:
    cands: List[str] = []
    for f in find_dnn_signal_files(sig_dir, mdark, alpha, yukawa):
        p = parse_signal_params_from_name(f)
        if p and int(p["mMed"]) == int(mmed) and abs(float(p["rinv"]) - float(rinv)) < 1e-6:
            cands.append(f)
    if not cands:
        return None
    # Prefer pure _nominal.root over _lundWeightNom.root.
    cands = sorted(cands, key=lambda x: (not os.path.basename(x).endswith("_nominal.root"), len(os.path.basename(x))))
    return cands[0]


# =============================================================================
# PNet signal directory parsing and score loading
# =============================================================================

def is_pnet_signal_dir(path: str) -> bool:
    base = os.path.basename(os.path.normpath(path))
    return base.startswith("t-channel_mMed-") and parse_signal_params_from_name(base) is not None


def find_pnet_signal_dirs(base_dir: str, mdark: int, alpha: str, yukawa: float) -> List[str]:
    entries = list_dir(base_dir)
    dirs = []
    for e in entries:
        if not is_pnet_signal_dir(e):
            continue
        p = parse_signal_params_from_name(e)
        if p is None:
            continue
        if int(p["mDark"]) != int(mdark):
            continue
        if str(p["alpha"]) != str(alpha):
            continue
        if abs(float(p["yukawa"]) - float(yukawa)) > 1e-9:
            continue
        dirs.append(e)
    return sorted(dirs, key=lambda x: (parse_signal_params_from_name(x)["mMed"], parse_signal_params_from_name(x)["rinv"]))


def find_pnet_reference_dir(base_dir: str, mmed: int, mdark: int, rinv: float, alpha: str, yukawa: float) -> Optional[str]:
    for d in find_pnet_signal_dirs(base_dir, mdark, alpha, yukawa):
        p = parse_signal_params_from_name(d)
        if p and int(p["mMed"]) == int(mmed) and abs(float(p["rinv"]) - float(rinv)) < 1e-6:
            return d
    return None


def load_pnet_event_scores(files: Sequence[str], branch: str = PNET_SCORE_BRANCH, tree_name: str = PNET_TREE,
                           max_events: int = -1, step_size: int = 50_000) -> Tuple[np.ndarray, np.ndarray]:
    """Return event-level max PNet score over AK8 jets and unit weights."""
    if ak is None:
        raise RuntimeError("awkward is required for PNet jagged jet-score branches")
    chunks: List[np.ndarray] = []
    total = 0
    for fpath in files:
        if max_events > 0 and total >= max_events:
            break
        open_path = eos_to_xrootd(fpath) if fpath.startswith("/store/") else fpath
        try:
            with uproot.open(open_path) as fin:
                if tree_name in fin:
                    tree = fin[tree_name]
                else:
                    tree = get_tree(fin, preferred=(tree_name, "Events"))
                if tree is None:
                    print(f"    [SKIP PNet] no tree in {fpath}")
                    continue
                if branch not in tree.keys():
                    print(f"    [SKIP PNet] missing branch {branch} in {fpath}")
                    continue
                for arrays in tree.iterate([branch], step_size=step_size, library="ak"):
                    raw = arrays[branch]
                    # Branch is jagged: events -> AK8 jets. Keep max per event.
                    if "var" in str(raw.type):
                        raw = ak.where(raw < -9.0, np.nan, raw)
                        score = ak.to_numpy(ak.fill_none(ak.max(raw, axis=1), np.nan))
                    else:
                        score = ak.to_numpy(raw)
                    score = np.asarray(score, dtype=float).reshape(-1)
                    score = score[np.isfinite(score)]
                    if len(score) == 0:
                        continue
                    if max_events > 0 and total + len(score) > max_events:
                        score = score[: max_events - total]
                    chunks.append(score)
                    total += len(score)
                    if max_events > 0 and total >= max_events:
                        break
        except Exception as e:
            print(f"    [SKIP PNet] {fpath}: {e}")
            continue
    if not chunks:
        return np.array([], dtype=float), np.array([], dtype=float)
    scores = np.concatenate(chunks).astype(float)
    weights = np.ones_like(scores, dtype=float)
    return scores, weights


# =============================================================================
# Metrics
# =============================================================================

def weighted_quantile(values: np.ndarray, quantile: float, weights: Optional[np.ndarray] = None) -> float:
    values = np.asarray(values, dtype=float)
    mask = np.isfinite(values)
    values = values[mask]
    if len(values) == 0:
        return np.nan
    if weights is None:
        return float(np.quantile(values, quantile))
    weights = np.asarray(weights, dtype=float)[mask]
    maskw = np.isfinite(weights) & (weights >= 0)
    values, weights = values[maskw], weights[maskw]
    if len(values) == 0 or np.sum(weights) <= 0:
        return np.nan
    order = np.argsort(values)
    values, weights = values[order], weights[order]
    cdf = np.cumsum(weights)
    cdf = cdf / cdf[-1]
    return float(np.interp(quantile, cdf, values))


def weighted_efficiency(scores: np.ndarray, weights: np.ndarray, cut: float, greater: bool = True) -> float:
    scores = np.asarray(scores, dtype=float)
    weights = np.asarray(weights, dtype=float)
    mask = np.isfinite(scores) & np.isfinite(weights) & (weights >= 0)
    scores, weights = scores[mask], weights[mask]
    den = np.sum(weights)
    if den <= 0:
        return np.nan
    passing = scores >= cut if greater else scores <= cut
    return float(np.sum(weights[passing]) / den)


def compute_auc_roc(sig_scores: np.ndarray, sig_w: np.ndarray, bkg_scores: np.ndarray, bkg_w: np.ndarray) -> Dict[str, object]:
    labels = np.concatenate([np.ones(len(sig_scores), dtype=int), np.zeros(len(bkg_scores), dtype=int)])
    scores = np.concatenate([sig_scores, bkg_scores]).astype(float)
    weights = np.concatenate([sig_w, bkg_w]).astype(float)
    mask = np.isfinite(scores) & np.isfinite(weights) & (weights >= 0)
    labels, scores, weights = labels[mask], scores[mask], weights[mask]
    if len(np.unique(labels)) < 2:
        raise RuntimeError("Need both signal and background labels for ROC/AUC")
    raw_auc = float(roc_auc_score(labels, scores, sample_weight=weights))
    if raw_auc < 0.5:
        roc_scores = -scores
        auc = 1.0 - raw_auc
        flipped = True
    else:
        roc_scores = scores
        auc = raw_auc
        flipped = False
    fpr, tpr, thr = roc_curve(labels, roc_scores, sample_weight=weights)
    return {"auc": auc, "raw_auc": raw_auc, "fpr": fpr, "tpr": tpr, "thresholds": thr, "flipped": flipped}


def compute_auc_only(sig_scores: np.ndarray, sig_w: np.ndarray, bkg_scores: np.ndarray, bkg_w: np.ndarray) -> Dict[str, object]:
    labels = np.concatenate([np.ones(len(sig_scores), dtype=int), np.zeros(len(bkg_scores), dtype=int)])
    scores = np.concatenate([sig_scores, bkg_scores]).astype(float)
    weights = np.concatenate([sig_w, bkg_w]).astype(float)
    mask = np.isfinite(scores) & np.isfinite(weights) & (weights >= 0)
    labels, scores, weights = labels[mask], scores[mask], weights[mask]
    if len(np.unique(labels)) < 2:
        raise RuntimeError("Need both signal and background labels for AUC")
    raw_auc = float(roc_auc_score(labels, scores, sample_weight=weights))
    if raw_auc < 0.5:
        return {"auc": 1.0 - raw_auc, "raw_auc": raw_auc, "flipped": True}
    return {"auc": raw_auc, "raw_auc": raw_auc, "flipped": False}


# =============================================================================
# Plotting helpers
# =============================================================================

def common_hist_bins(arrays: Sequence[np.ndarray], nbins: int = 50,
                     force_range: Optional[Tuple[float, float]] = None) -> Optional[np.ndarray]:
    vals: List[np.ndarray] = []
    for a in arrays:
        a = np.asarray(a, dtype=float).reshape(-1)
        a = a[np.isfinite(a)]
        if len(a):
            vals.append(a)
    if not vals:
        return None
    allv = np.concatenate(vals)
    if force_range is not None:
        lo, hi = force_range
    else:
        lo, hi = np.percentile(allv, [0.5, 99.5])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = float(np.nanmin(allv)), float(np.nanmax(allv))
    if hi <= lo:
        hi = lo + 1.0
    return np.linspace(lo, hi, nbins + 1)


def normalized_step(ax: plt.Axes, values: np.ndarray, weights: Optional[np.ndarray], bins: np.ndarray,
                    label: str, color: str, linewidth: float = 3.0) -> None:
    values = np.asarray(values, dtype=float).reshape(-1)
    mask = np.isfinite(values)
    values = values[mask]
    if weights is not None:
        weights = np.asarray(weights, dtype=float).reshape(-1)[mask]
        maskw = np.isfinite(weights) & (weights >= 0)
        values, weights = values[maskw], weights[maskw]
        if len(weights) and np.sum(weights) > 0:
            weights = weights / np.sum(weights)
    if len(values) == 0:
        return
    ax.hist(values, bins=bins, weights=weights, histtype="step", density=(weights is None),
            color=color, linewidth=linewidth, label=label)


def plot_score_distribution(tagger: str, sig_scores: np.ndarray, sig_w: np.ndarray,
                            bkg_scores: np.ndarray, bkg_w: np.ndarray,
                            out_dir: str, ref_label: str, nbins: int = 50) -> None:
    mkdir(out_dir)
    bins = common_hist_bins([sig_scores, bkg_scores], nbins=nbins, force_range=(0.0, 1.0))
    if bins is None:
        return
    for logy in (False, True):
        fig, ax = plt.subplots(figsize=(10, 8))
        normalized_step(ax, bkg_scores, bkg_w, bins, "QCD background", color="black", linewidth=3.0)
        normalized_step(ax, sig_scores, sig_w, bins, ref_label, color="red", linewidth=3.0)
        ax.set_xlabel(f"{tagger} score", fontsize=24)
        ax.set_ylabel("Normalized events / bin", fontsize=24)
        ax.set_xlim(0.0, 1.0)
        ax.tick_params(labelsize=18)
        ax.legend(fontsize=15, framealpha=0.95)
        ax.grid(True, axis="y", alpha=0.25)
        if logy:
            ax.set_yscale("log")
            ax.set_ylim(bottom=1e-5)
        cms_label(ax, fontsize=26)
        # Intentionally no title.
        plt.tight_layout()
        suffix = "logy" if logy else "linear"
        save_fig(fig, os.path.join(out_dir, f"{tagger.lower()}_score_distribution_{suffix}"))


def plot_reference_roc(tagger: str, roc_info: Dict[str, object], out_dir: str, ref_label: str) -> None:
    mkdir(out_dir)
    auc = float(roc_info["auc"])
    fpr = np.asarray(roc_info["fpr"], dtype=float)
    tpr = np.asarray(roc_info["tpr"], dtype=float)
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.plot(fpr, tpr, color="red", linewidth=3.0, label=f"{ref_label}  (AUC = {auc:.3f})")
    ax.plot([0, 1], [0, 1], "--", color="black", linewidth=1.5, alpha=0.7, label="Random")
    ax.set_xlabel("False Positive Rate (FPR)", fontsize=24)
    ax.set_ylabel("True Positive Rate (TPR)", fontsize=24)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.tick_params(labelsize=18)
    ax.legend(fontsize=14, loc="lower right", framealpha=0.95)
    ax.grid(True, alpha=0.25)
    cms_label(ax, fontsize=26)
    # Intentionally no title.
    plt.tight_layout()
    save_fig(fig, os.path.join(out_dir, f"{tagger.lower()}_roc_reference"))


def paper_grid_from_dict(values: Dict[Tuple[int, float], float], z_label: str, out_prefix: str,
                         training_samples: Optional[List[Dict[str, object]]] = None,
                         vmin: float = 0.0, vmax: float = 1.0,
                         fmt: str = ".2f", cmap: str = "viridis",
                         title: Optional[str] = None) -> None:
    if not values:
        print(f"  [SKIP] No values for {z_label}")
        return
    mmed_vals = sorted({int(k[0]) for k in values.keys()})
    rinv_vals = sorted({float(k[1]) for k in values.keys()})
    Z = np.full((len(rinv_vals), len(mmed_vals)), np.nan)
    for (mmed, rinv), val in values.items():
        Z[rinv_vals.index(float(rinv)), mmed_vals.index(int(mmed))] = val

    fig, ax = plt.subplots(figsize=(16, 12))
    im = ax.imshow(Z, origin="lower", aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax,
                   extent=[-0.5, len(mmed_vals) - 0.5, -0.5, len(rinv_vals) - 0.5])
    outline = [pe.Stroke(linewidth=4, foreground="black"), pe.Normal()]
    for i in range(len(rinv_vals)):
        for j in range(len(mmed_vals)):
            val = Z[i, j]
            if np.isfinite(val):
                ax.text(j, i + 0.12, f"{val:{fmt}}", ha="center", va="center",
                        fontsize=23, color="white", path_effects=outline)
            else:
                ax.text(j, i, "—", ha="center", va="center", fontsize=24, color="white", path_effects=outline)
    if training_samples is not None:
        for s in training_samples:
            m, r = int(s["mMed"]), float(s["rinv"])
            if m in mmed_vals and r in rinv_vals:
                ax.scatter(mmed_vals.index(m), rinv_vals.index(r) - 0.24,
                           marker="*", s=500, facecolor="white", edgecolor="black",
                           linewidth=1.5, zorder=10)
    ax.set_xticks(range(len(mmed_vals)))
    ax.set_xticklabels([str(m) for m in mmed_vals], fontsize=28)
    ax.set_yticks(range(len(rinv_vals)))
    ax.set_yticklabels([f"{r:.2g}" if abs(r - round(r, 1)) > 1e-8 else f"{r:.1f}" for r in rinv_vals], fontsize=28)
    ax.set_xlabel(r"$m_{\Phi}\ [\mathrm{GeV}]$", fontsize=44, labelpad=22)
    ax.set_ylabel(r"$r_{\mathrm{inv}}$", fontsize=44, labelpad=22)
    cb = fig.colorbar(im, ax=ax, pad=0.03)
    cb.set_label(z_label, fontsize=36, labelpad=25)
    cb.ax.tick_params(labelsize=28)
    cms_label(ax, fontsize=38, com=13)
    if title:
        ax.set_title(title, fontsize=28, pad=24)
    plt.tight_layout()
    save_fig(fig, out_prefix)


def write_auc_tables(rows: List[Dict[str, object]], out_prefix: str) -> None:
    mkdir(os.path.dirname(out_prefix))
    csv_path = f"{out_prefix}.csv"
    txt_path = f"{out_prefix}.txt"
    rows = sorted(rows, key=lambda r: (int(r["mMed"]), float(r["rinv"]), int(r.get("mDark", 0))))
    fields = ["mMed", "mDark", "rinv", "alpha", "yukawa", "auc"]
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})
    with open(txt_path, "w") as f:
        f.write(f"{'mMed':>8} {'mDark':>8} {'rinv':>8} {'alpha':>8} {'yukawa':>8} {'AUC':>10}\n")
        f.write("-" * 62 + "\n")
        for r in rows:
            f.write(f"{int(r['mMed']):8d} {int(r.get('mDark', 0)):8d} {float(r['rinv']):8.2f} "
                    f"{str(r.get('alpha', '')):>8} {float(r.get('yukawa', 1.0)):8.3g} {float(r['auc']):10.5f}\n")
    print(f"  -> {csv_path}")
    print(f"  -> {txt_path}")


def plot_selected_roc_overlay(tagger: str, results: Dict[Tuple[int, float], Dict[str, object]],
                              out_dir: str, max_curves: int = 10) -> None:
    if not results:
        return
    mkdir(out_dir)
    keys = sorted(results.keys(), key=lambda x: (x[0], x[1]))
    preferred = [k for k in keys if abs(k[1] - 0.3) < 1e-6] + [k for k in keys if k[0] == 2000]
    selected: List[Tuple[int, float]] = []
    seen = set()
    for k in preferred + keys:
        if k not in seen:
            selected.append(k)
            seen.add(k)
        if len(selected) >= max_curves:
            break
    if not selected:
        return
    fig, ax = plt.subplots(figsize=(10, 9))
    cmap = plt.cm.viridis
    norm = mcolors.Normalize(vmin=min(k[0] for k in selected), vmax=max(k[0] for k in selected))
    for k in selected:
        info = results[k]
        ax.plot(info["fpr"], info["tpr"], linewidth=2.4, color=cmap(norm(k[0])),
                label=rf"$m_\Phi={k[0]}$, $r_{{\rm inv}}={k[1]:.1f}$, AUC={info['auc']:.3f}")
    ax.plot([0, 1], [0, 1], "--", color="black", linewidth=1.4, alpha=0.65)
    ax.set_xlabel("False Positive Rate (FPR)", fontsize=24)
    ax.set_ylabel("True Positive Rate (TPR)", fontsize=24)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.tick_params(labelsize=18)
    ax.legend(fontsize=11, loc="lower right", framealpha=0.94)
    ax.grid(True, alpha=0.25)
    cms_label(ax, fontsize=26)
    plt.tight_layout()
    save_fig(fig, os.path.join(out_dir, f"{tagger.lower()}_roc_selected_points"))


def plot_efficiency_vs_scorecut(tagger: str, selected_scores: Dict[Tuple[int, float], Tuple[np.ndarray, np.ndarray]],
                                bkg_scores: np.ndarray, bkg_w: np.ndarray, out_dir: str) -> None:
    if not selected_scores:
        return
    mkdir(out_dir)
    cuts = np.linspace(0.0, 1.0, 101)
    fig, ax = plt.subplots(figsize=(10, 8))
    bkg_eff = [weighted_efficiency(bkg_scores, bkg_w, c, greater=True) for c in cuts]
    ax.plot(cuts, bkg_eff, "--", color="black", linewidth=2.5, label="QCD background")
    for (mmed, rinv), (scores, weights) in sorted(selected_scores.items(), key=lambda x: (x[0][0], x[0][1])):
        eff = [weighted_efficiency(scores, weights, c, greater=True) for c in cuts]
        ax.plot(cuts, eff, linewidth=2.4, label=rf"$m_\Phi={mmed}$, $r_{{\rm inv}}={rinv:.1f}$")
    ax.set_xlabel(f"{tagger} score cut", fontsize=24)
    ax.set_ylabel("Efficiency", fontsize=24)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.05)
    ax.tick_params(labelsize=18)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=12, framealpha=0.95, ncol=2)
    cms_label(ax, fontsize=26)
    plt.tight_layout()
    save_fig(fig, os.path.join(out_dir, f"{tagger.lower()}_efficiency_vs_scorecut_selected_points"))


def reference_label(mmed: int, mdark: int, rinv: float) -> str:
    return (
        rf"Signal $m_\Phi={mmed}$ GeV, $m_{{\rm dark}}={mdark}$ GeV"
        "\n"
        rf"$\Lambda_{{\rm dark}}=1$ , $r_{{\rm inv}}={rinv:.1f}$"
    )
# =============================================================================
# DNN input-variable plots
# =============================================================================

def plot_input_variables(sig_file: str, bkg_file: str, branches: List[str], out_dir: str, max_events: int) -> None:
    mkdir(out_dir)
    made = 0
    for br in branches:
        try:
            s = to_flat_numpy(load_branch(sig_file, br, max_events=max_events))
            b = to_flat_numpy(load_branch(bkg_file, br, max_events=max_events))
        except Exception as e:
            print(f"    [SKIP input] {br}: {e}")
            continue
        if len(s) == 0 or len(b) == 0:
            print(f"    [SKIP input] {br}: empty/non-numeric")
            continue
        bins = common_hist_bins([s, b], nbins=45)
        if bins is None:
            continue
        fig, ax = plt.subplots(figsize=(10, 8))
        normalized_step(ax, b, None, bins, f"QCD background (N={len(b):,})", color="black", linewidth=3.0)
        normalized_step(ax, s, None, bins, f"Signal (N={len(s):,})", color="red", linewidth=3.0)
        ax.set_xlabel(br, fontsize=24)
        ax.set_ylabel("Normalized events / bin", fontsize=24)
        ax.tick_params(labelsize=18)
        ax.legend(fontsize=15, framealpha=0.95)
        ax.grid(True, axis="y", alpha=0.25)
        cms_label(ax, fontsize=26)
        plt.tight_layout()
        save_fig(fig, os.path.join(out_dir, f"dnn_input_{clean_name(br)}"))
        made += 1
    print(f"  Made {made} DNN input-variable plots")


# =============================================================================
# Processing functions
# =============================================================================

@dataclass
class ScanResult:
    auc_map: Dict[Tuple[int, float], float]
    eff_maps: Dict[float, Dict[Tuple[int, float], float]]
    roc_results: Dict[Tuple[int, float], Dict[str, object]]
    selected_scores: Dict[Tuple[int, float], Tuple[np.ndarray, np.ndarray]]
    table_rows: List[Dict[str, object]]


def scan_requested(args: argparse.Namespace) -> bool:
    return any([
        args.make_auc_maps,
        args.make_eff_maps,
        args.make_auc_tables,
        args.make_selected_roc,
        args.make_eff_vs_scorecut,
    ])


def selected_keys(args: argparse.Namespace) -> set:
    out = set()
    for token in args.selected_points:
        try:
            mm, rr = token.split(":")
            out.add((int(mm), float(rr)))
        except Exception:
            pass
    return out


def scan_dnn(args: argparse.Namespace, bkg_scores: np.ndarray, bkg_w: np.ndarray) -> ScanResult:
    sig_files = find_dnn_signal_files(args.dnn_sig_dir, args.ref_mdark, args.ref_alpha, args.ref_yukawa)
    print(f"  DNN scan files after exact nominal filtering: {len(sig_files):,}")
    return scan_file_like_points(
        tagger="DNN",
        point_paths=sig_files,
        path_to_scores=lambda p: load_dnn_scores_weights(p, args.dnn_score_branch, args.dnn_weight_branch, args.max_events_scores),
        path_to_meta=lambda p: parse_signal_params_from_name(p),
        bkg_scores=bkg_scores,
        bkg_w=bkg_w,
        args=args,
    )


def scan_pnet(args: argparse.Namespace, bkg_scores: np.ndarray, bkg_w: np.ndarray) -> ScanResult:
    sig_dirs = find_pnet_signal_dirs(args.pnet_sig_base, args.ref_mdark, args.ref_alpha, args.ref_yukawa)
    print(f"  PNet scan directories: {len(sig_dirs):,}")

    def load_dir_scores(d: str) -> Tuple[np.ndarray, np.ndarray]:
        files = root_files_in_dir(d)
        return load_pnet_event_scores(files, branch=args.pnet_score_branch, tree_name=args.pnet_tree,
                                      max_events=args.max_events_scores, step_size=args.pnet_step_size)

    return scan_file_like_points(
        tagger="PNet",
        point_paths=sig_dirs,
        path_to_scores=load_dir_scores,
        path_to_meta=lambda p: parse_signal_params_from_name(p),
        bkg_scores=bkg_scores,
        bkg_w=bkg_w,
        args=args,
    )


def scan_file_like_points(tagger: str,
                          point_paths: Sequence[str],
                          path_to_scores,
                          path_to_meta,
                          bkg_scores: np.ndarray,
                          bkg_w: np.ndarray,
                          args: argparse.Namespace) -> ScanResult:
    bkg_cuts: Dict[float, float] = {}
    if args.make_eff_maps:
        for beff in args.bkg_eff_targets:
            q = max(0.0, min(1.0, 1.0 - beff))
            cut = weighted_quantile(bkg_scores, q, bkg_w)
            actual = weighted_efficiency(bkg_scores, bkg_w, cut, greater=True)
            bkg_cuts[beff] = cut
            print(f"  {tagger} score cut for FPR/QCD eff {beff:g}: cut={cut:.6g}, actual={actual:.6g}")

    want_auc = args.make_auc_maps or args.make_auc_tables or args.make_selected_roc
    want_roc = args.make_selected_roc
    want_selected_scores = args.make_eff_vs_scorecut
    wanted_keys = selected_keys(args)

    auc_map: Dict[Tuple[int, float], float] = {}
    eff_maps: Dict[float, Dict[Tuple[int, float], float]] = {beff: {} for beff in args.bkg_eff_targets}
    roc_results: Dict[Tuple[int, float], Dict[str, object]] = {}
    selected_scores: Dict[Tuple[int, float], Tuple[np.ndarray, np.ndarray]] = {}
    table_rows: List[Dict[str, object]] = []

    n = len(point_paths)
    for i, path in enumerate(point_paths, start=1):
        meta = path_to_meta(path)
        if meta is None:
            continue
        key = (int(meta["mMed"]), float(meta["rinv"]))
        try:
            sig_scores, sig_w = path_to_scores(path)
            if len(sig_scores) == 0:
                continue
            info = None
            if want_auc:
                info = (
                    compute_auc_roc(sig_scores, sig_w, bkg_scores, bkg_w)
                    if want_roc
                    else compute_auc_only(sig_scores, sig_w, bkg_scores, bkg_w)
                )
                auc_map[key] = float(info["auc"])
                if want_roc:
                    roc_results[key] = info
                table_rows.append({
                    "mMed": int(meta["mMed"]),
                    "mDark": int(meta["mDark"]),
                    "rinv": float(meta["rinv"]),
                    "alpha": str(meta["alpha"]),
                    "yukawa": float(meta["yukawa"]),
                    "auc": float(info["auc"]),
                })
            if args.make_eff_maps:
                for beff, cut in bkg_cuts.items():
                    eff_maps[beff][key] = weighted_efficiency(sig_scores, sig_w, cut, greater=True)
            if want_selected_scores and key in wanted_keys:
                selected_scores[key] = (sig_scores, sig_w)
            if args.verbose_scan:
                msg = f"    [{i:4d}/{n}] mMed={key[0]:5d} rinv={key[1]:.2f}"
                if info is not None:
                    msg += f" AUC={info['auc']:.4f}"
                print(msg)
            elif i % max(1, n // 10) == 0 or i == n:
                print(f"    scanned {i:4d}/{n}")
        except Exception as e:
            if args.verbose_scan:
                print(f"    [SKIP] {path}: {e}")
            continue

    return ScanResult(auc_map, eff_maps, roc_results, selected_scores, table_rows)


def process_dnn(args: argparse.Namespace) -> None:
    if not args.run_dnn:
        print("\n[DNN] skipped")
        return
    print("\n" + "=" * 100)
    print("Processing DNN")
    print("=" * 100)
    tag_dir = os.path.join(args.out_dir, "DNN")
    score_dir = os.path.join(tag_dir, "01_score_distributions")
    roc_dir = os.path.join(tag_dir, "02_roc_curves")
    map_dir = os.path.join(tag_dir, "03_maps")
    cut_dir = os.path.join(tag_dir, "04_efficiency_vs_scorecut")
    table_dir = os.path.join(tag_dir, "05_tables")
    for d in (score_dir, roc_dir, map_dir, cut_dir, table_dir):
        mkdir(d)

    ref = args.dnn_ref_sig_file or find_dnn_reference_signal(
        args.dnn_sig_dir, args.ref_mmed, args.ref_mdark, args.ref_rinv, args.ref_alpha, args.ref_yukawa
    )
    if not ref:
        print("  [ERROR] DNN reference signal not found")
        return
    print(f"  Background: {args.dnn_bkg_file}")
    print(f"  Reference signal: {ref}")
    print(f"  Score branch: {args.dnn_score_branch}")

    bkg_scores, bkg_w = load_dnn_scores_weights(args.dnn_bkg_file, args.dnn_score_branch, args.dnn_weight_branch, args.max_events_scores)
    sig_scores, sig_w = load_dnn_scores_weights(ref, args.dnn_score_branch, args.dnn_weight_branch, args.max_events_scores)
    print(f"  Loaded QCD scores: {len(bkg_scores):,}")
    print(f"  Loaded reference signal scores: {len(sig_scores):,}")

    ref_lab = reference_label(args.ref_mmed, args.ref_mdark, args.ref_rinv)
    plot_score_distribution("DNN", sig_scores, sig_w, bkg_scores, bkg_w, score_dir, ref_lab, args.score_bins)
    roc_info = compute_auc_roc(sig_scores, sig_w, bkg_scores, bkg_w)
    print(f"  Reference DNN AUC: {roc_info['auc']:.5f} (raw={roc_info['raw_auc']:.5f}, flipped={roc_info['flipped']})")
    plot_reference_roc("DNN", roc_info, roc_dir, ref_lab)

    if not scan_requested(args):
        print("  DNN full scan skipped. Add --make-auc-maps/--make-eff-maps/--make-auc-tables if needed.")
        return

    res = scan_dnn(args, bkg_scores, bkg_w)
    if args.make_auc_maps and res.auc_map:
        paper_grid_from_dict(res.auc_map, "DNN AUC", os.path.join(map_dir, "dnn_auc_grid"),
                             training_samples=DNN_TRAINING_SAMPLES, vmin=0.5, vmax=1.0, fmt=".2f")
    if args.make_eff_maps:
        for beff, emap in res.eff_maps.items():
            if emap:
                paper_grid_from_dict(emap, f"DNN signal efficiency @ FPR {beff:g}",
                                     os.path.join(map_dir, f"dnn_signal_eff_map_fpr_{str(beff).replace('.', 'p')}"),
                                     training_samples=DNN_TRAINING_SAMPLES, vmin=0.0, vmax=1.0, fmt=".2f")
    if args.make_auc_tables and res.table_rows:
        write_auc_tables(res.table_rows, os.path.join(table_dir, "dnn_auc_table"))
    if args.make_selected_roc:
        plot_selected_roc_overlay("DNN", res.roc_results, roc_dir, args.max_roc_curves)
    if args.make_eff_vs_scorecut:
        plot_efficiency_vs_scorecut("DNN", res.selected_scores, bkg_scores, bkg_w, cut_dir)


def process_pnet(args: argparse.Namespace) -> None:
    if not args.run_pnet:
        print("\n[PNet] skipped")
        return
    print("\n" + "=" * 100)
    print("Processing PNet")
    print("=" * 100)
    tag_dir = os.path.join(args.out_dir, "PNet")
    score_dir = os.path.join(tag_dir, "01_score_distributions")
    roc_dir = os.path.join(tag_dir, "02_roc_curves")
    map_dir = os.path.join(tag_dir, "03_maps")
    cut_dir = os.path.join(tag_dir, "04_efficiency_vs_scorecut")
    table_dir = os.path.join(tag_dir, "05_tables")
    for d in (score_dir, roc_dir, map_dir, cut_dir, table_dir):
        mkdir(d)

    ref_dir = args.pnet_ref_sig_dir or find_pnet_reference_dir(
        args.pnet_sig_base, args.ref_mmed, args.ref_mdark, args.ref_rinv, args.ref_alpha, args.ref_yukawa
    )
    if not ref_dir:
        print("  [ERROR] PNet reference signal directory not found")
        print(f"          base={args.pnet_sig_base}")
        return
    ref_files = root_files_in_dir(ref_dir)
    if not ref_files:
        print(f"  [ERROR] No ROOT files in PNet reference dir: {ref_dir}")
        return
    print(f"  PNet signal base: {args.pnet_sig_base}")
    print(f"  Reference signal dir: {ref_dir}")
    print(f"  Reference signal files: {len(ref_files)}")
    print(f"  Tree: {args.pnet_tree}")
    print(f"  Branch: {args.pnet_score_branch}")
    print("  Event score: max branch value over AK8 jets")

    qcd_files = (
        args.pnet_qcd_files
        if args.pnet_qcd_files
        else discover_trimmed_qcd_files(args.pnet_qcd_base, args.pnet_year)
    )
    if not args.pnet_qcd_files:
        print(f"  PNet QCD base: {args.pnet_qcd_base}")
        print(f"  PNet QCD year: {args.pnet_year}")
    print(f"  PNet QCD files: {len(qcd_files)}")
    if args.pnet_qcd_cache and os.path.exists(args.pnet_qcd_cache) and not args.pnet_no_cache:
        print(f"  Loading cached PNet QCD scores: {args.pnet_qcd_cache}")
        bkg_scores = np.load(args.pnet_qcd_cache)
        bkg_w = np.ones_like(bkg_scores, dtype=float)
    else:
        bkg_scores, bkg_w = load_pnet_event_scores(qcd_files, branch=args.pnet_score_branch, tree_name=args.pnet_tree,
                                                   max_events=args.max_events_scores, step_size=args.pnet_step_size)
        if args.pnet_qcd_cache and len(bkg_scores) and not args.pnet_no_cache:
            np.save(args.pnet_qcd_cache, bkg_scores)
            print(f"  Saved PNet QCD cache: {args.pnet_qcd_cache}")

    sig_scores, sig_w = load_pnet_event_scores(ref_files, branch=args.pnet_score_branch, tree_name=args.pnet_tree,
                                               max_events=args.max_events_scores, step_size=args.pnet_step_size)
    print(f"  Loaded QCD PNet scores: {len(bkg_scores):,}")
    print(f"  Loaded reference signal PNet scores: {len(sig_scores):,}")
    if len(bkg_scores) == 0 or len(sig_scores) == 0:
        print("  [ERROR] Empty PNet scores")
        return

    ref_lab = reference_label(args.ref_mmed, args.ref_mdark, args.ref_rinv)
    plot_score_distribution("PNet", sig_scores, sig_w, bkg_scores, bkg_w, score_dir, ref_lab, args.score_bins)
    roc_info = compute_auc_roc(sig_scores, sig_w, bkg_scores, bkg_w)
    print(f"  Reference PNet AUC: {roc_info['auc']:.5f} (raw={roc_info['raw_auc']:.5f}, flipped={roc_info['flipped']})")
    plot_reference_roc("PNet", roc_info, roc_dir, ref_lab)

    if not scan_requested(args):
        print("  PNet full scan skipped. Add --make-auc-maps/--make-eff-maps/--make-auc-tables if needed.")
        return

    res = scan_pnet(args, bkg_scores, bkg_w)
    if args.make_auc_maps and res.auc_map:
        paper_grid_from_dict(res.auc_map, "PNet AUC", os.path.join(map_dir, "pnet_auc_grid"),
                             training_samples=PNET_TRAINING_SAMPLES, vmin=0.5, vmax=1.0, fmt=".2f")
    if args.make_eff_maps:
        for beff, emap in res.eff_maps.items():
            if emap:
                paper_grid_from_dict(emap, f"PNet signal efficiency @ FPR {beff:g}",
                                     os.path.join(map_dir, f"pnet_signal_eff_map_fpr_{str(beff).replace('.', 'p')}"),
                                     training_samples=PNET_TRAINING_SAMPLES, vmin=0.0, vmax=1.0, fmt=".2f")
    if args.make_auc_tables and res.table_rows:
        write_auc_tables(res.table_rows, os.path.join(table_dir, "pnet_auc_table"))
    if args.make_selected_roc:
        plot_selected_roc_overlay("PNet", res.roc_results, roc_dir, args.max_roc_curves)
    if args.make_eff_vs_scorecut:
        plot_efficiency_vs_scorecut("PNet", res.selected_scores, bkg_scores, bkg_w, cut_dir)


def main() -> None:
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Make supplementary SVJ DNN and PNet plots. Fast by default; scans require flags."
    )
    ap.add_argument("--out-dir", default="./supplementry_plots_postTrimandgapveto")

    ap.add_argument("--run-dnn", action="store_true")
    ap.add_argument("--run-pnet", action="store_true")
    ap.add_argument("--all", action="store_true")

    ap.add_argument("--ref-mmed", type=int, default=2000)
    ap.add_argument("--ref-mdark", type=int, default=20)
    ap.add_argument("--ref-rinv", type=float, default=0.3)
    ap.add_argument("--ref-alpha", default="peak")
    ap.add_argument("--ref-yukawa", type=float, default=1.0)

    ap.add_argument("--dnn-sig-dir", default=DEFAULT_DNN_SIG_DIR)
    ap.add_argument("--dnn-bkg-file", default=DEFAULT_DNN_BKG_FILE)
    ap.add_argument("--dnn-ref-sig-file", default="")
    ap.add_argument("--dnn-score-branch", default="scores")
    ap.add_argument("--dnn-weight-branch", default="weights")

    ap.add_argument("--pnet-sig-base", default=DEFAULT_PNET_SIG_BASE,
                    help="Parent dir containing t-channel_mMed-* PNet signal directories")
    ap.add_argument("--pnet-ref-sig-dir", default="")
    ap.add_argument("--pnet-qcd-base", default=DEFAULT_PNET_QCD_BASE,
                    help="Parent dir containing YEAR_QCD_Pt_* trimmed-QCD directories")
    ap.add_argument("--pnet-year", default=DEFAULT_PNET_YEAR)
    ap.add_argument("--pnet-qcd-files", nargs="*", default=[])
    ap.add_argument("--pnet-qcd-cache", default="./pnet_qcd_scores_2018_postTrimandgapveto.npy")
    ap.add_argument("--pnet-no-cache", action="store_true")
    ap.add_argument("--pnet-tree", default=PNET_TREE)
    ap.add_argument("--pnet-score-branch", default=PNET_SCORE_BRANCH)
    ap.add_argument("--pnet-step-size", type=int, default=50_000)

    ap.add_argument("--make-dnn-inputs", action="store_true")
    ap.add_argument("--dnn-input-sig-file", default="")
    ap.add_argument("--dnn-input-bkg-file", default="")
    ap.add_argument("--dnn-input-vars", nargs="*", default=[])
    ap.add_argument("--auto-input-vars", action="store_true")
    ap.add_argument("--max-input-plots", type=int, default=30)

    # Full-scan outputs: all off by default.
    ap.add_argument("--make-auc-maps", action="store_true")
    ap.add_argument("--make-eff-maps", action="store_true")
    ap.add_argument("--make-auc-tables", action="store_true")
    ap.add_argument("--make-selected-roc", action="store_true")
    ap.add_argument("--make-eff-vs-scorecut", action="store_true")

    ap.add_argument("--max-events-scores", type=int, default=-1)
    ap.add_argument("--max-events-inputs", type=int, default=200_000)
    ap.add_argument("--score-bins", type=int, default=50)
    ap.add_argument("--bkg-eff-targets", type=float, nargs="+", default=[0.10, 0.01, 0.001])
    ap.add_argument("--selected-points", nargs="+", default=[
        "600:0.3", "1000:0.3", "2000:0.3", "3000:0.3", "4000:0.3",
        "2000:0.1", "2000:0.5", "2000:0.7"
    ])
    ap.add_argument("--max-roc-curves", type=int, default=10)
    ap.add_argument("--verbose-scan", action="store_true")
    ap.add_argument("--print-branches", action="store_true")

    args = ap.parse_args()

    if args.all:
        args.run_dnn = True
        args.run_pnet = True
    if not args.run_dnn and not args.run_pnet:
        args.run_dnn = True

    mkdir(args.out_dir)
    print("\n" + "=" * 100)
    print("SVJ supplementary material plotter FIXED_v2")
    print("=" * 100)
    print(f"Output directory: {args.out_dir}")
    print(f"Reference point: mMed={args.ref_mmed}, mDark={args.ref_mdark}, rinv={args.ref_rinv}, "
          f"alpha={args.ref_alpha}, yukawa={args.ref_yukawa:g}")
    print(f"Full scan requested: {scan_requested(args)}")

    if args.print_branches:
        dnn_ref = args.dnn_ref_sig_file or find_dnn_reference_signal(
            args.dnn_sig_dir, args.ref_mmed, args.ref_mdark, args.ref_rinv, args.ref_alpha, args.ref_yukawa
        ) or ""
        for label, path in [("DNN background", args.dnn_bkg_file), ("DNN reference signal", dnn_ref)]:
            if path and (file_exists(path) or path.startswith("root://") or path.startswith("/store/")):
                print(f"\n[{label}] {path}")
                for br in tree_keys(path):
                    print(f"  {br}")
            else:
                print(f"\n[{label}] no file found/provided")
        # PNet branch is known, but also print one ref file's branches if available.
        pnet_ref = args.pnet_ref_sig_dir or find_pnet_reference_dir(
            args.pnet_sig_base, args.ref_mmed, args.ref_mdark, args.ref_rinv, args.ref_alpha, args.ref_yukawa
        ) or ""
        pnet_files = root_files_in_dir(pnet_ref)[:1] if pnet_ref else []
        if pnet_files:
            print(f"\n[PNet reference example] {pnet_files[0]}")
            for br in tree_keys(pnet_files[0]):
                print(f"  {br}")
        else:
            print("\n[PNet reference example] no file found")
        return

    if args.make_dnn_inputs:
        print("\n" + "=" * 100)
        print("DNN input-variable plots")
        print("=" * 100)
        dnn_ref = args.dnn_ref_sig_file or find_dnn_reference_signal(
            args.dnn_sig_dir, args.ref_mmed, args.ref_mdark, args.ref_rinv, args.ref_alpha, args.ref_yukawa
        ) or ""
        input_sig = args.dnn_input_sig_file or dnn_ref
        input_bkg = args.dnn_input_bkg_file or args.dnn_bkg_file
        branches = list(args.dnn_input_vars)
        if not branches and args.auto_input_vars and input_sig and input_bkg:
            branches = discover_numeric_scalar_branches(
                input_sig, input_bkg, exclude=[args.dnn_score_branch, args.dnn_weight_branch], max_vars=args.max_input_plots
            )
            print(f"  Auto-discovered DNN input branches: {branches}")
        if not branches:
            print("  [WARN] No DNN input variables requested/discovered.")
            print("         ABCD_scores files usually only have scores/weights; pass skim files and --dnn-input-vars for real inputs.")
        else:
            plot_input_variables(input_sig, input_bkg, branches,
                                 os.path.join(args.out_dir, "DNN", "00_input_variables"), args.max_events_inputs)

    process_dnn(args)
    process_pnet(args)

    print("\n" + "=" * 100)
    print(f"Done. Plots saved under: {args.out_dir}")
    print("=" * 100 + "\n")


if __name__ == "__main__":
    main()
