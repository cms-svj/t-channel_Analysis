#!/usr/bin/env python3
"""
DNN AUC grid, mDark vs mMed (rinv=0.3, yukawa=1 fixed), against both QCD and
ttbar -- the DNN counterpart to makeAUCtables_pnet_extra_axes.py's PNET
mDark-vs-mMed grids.

Unlike ParticleNet, the event-level DNN score isn't read from the raw EOS
skims -- it's read from pre-computed per-event score+weight ROOT files
(the same source makeAUCtables_DNN.py's mMed-vs-rinv grid uses), generalized
here from its hardcoded mDark=20 to scan the full mDark axis:
  Signal: .../scores_signals_GapJetVeto/2018/nominal/
          ABCD_scores_nominal_t-channel_mMed-*_mDark-*_rinv-0p3_alpha-peak_yukawa-1_{nominal,lundWeightNom}.root
          (both suffixes are identical duplicates for every point; "nominal" is used)
  Background: .../DNN_scores_and_ABCD_histograms/DNN_scores/scores_backgrounds_GapJetVeto/2018/nominal/
          ABCD_scores_{QCD,TTJets}.root
Both signal and background files carry pre-computed "scores"/"weights"
branches (weights already include the full per-event physics weight), so no
additional Lund-plane reweighting is needed here -- it's baked into the
signal weights by the score-production pipeline.

Two-step, mirroring the other AUC/signal-efficiency scripts in this repo:
  --compute   Reads the score files, computes weighted AUC per grid point vs
              QCD and vs ttbar. Writes dnn_mdark_scan_auc_qcd.csv and
              dnn_mdark_scan_auc_ttbar.csv.
  (default)   Reads those CSVs and draws the two heat maps.

Usage:
  source condor/initCondor.sh
  python3 makeAUCtables_DNN_extra_axes.py --compute
  python3 makeAUCtables_DNN_extra_axes.py
"""
import argparse
import csv
import glob
import hashlib
import os
import re
import sys

import numpy as np
import uproot
from sklearn.metrics import roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import makeAUCtables_pnet as auc_plot

OUT_DIR = os.path.join(HERE, "dnn_auc_extra_axes")
QCD_CSV = os.path.join(OUT_DIR, "dnn_mdark_scan_auc_qcd.csv")
TTBAR_CSV = os.path.join(OUT_DIR, "dnn_mdark_scan_auc_ttbar.csv")
PLOT_DIR = os.path.join(
    HERE, "SVJ-tchannel-Run2_site", "supplemetry", "assets_v2", "auc"
)

SIG_DIR = (
    "/uscms/home/nparmar/nobackup/SVJ/eventLevelTagger/output/dataset11/training_set2/"
    "application_results/sdt_allBkgs_disco_0p001_closure_0p06_damp_1_net_64_32_16_8_"
    "trainingAllYears_met250_dnn85_eval2018/v_5/ML_fit_results_grid_optimization_binned/"
    "scores_signals_GapJetVeto/2018/nominal"
)
BKG_DIR = (
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/"
    "DNN_scores_and_ABCD_histograms/DNN_scores/scores_backgrounds_GapJetVeto/2018/nominal"
)
BACKGROUND_FILES = {
    "qcd": f"{BKG_DIR}/ABCD_scores_QCD.root",
    "ttbar": f"{BKG_DIR}/ABCD_scores_TTJets.root",
}

FIXED_RINV = 0.3
FIXED_YUKAWA = 1.0
AUC_CAP = 500_000
SEED = 42


def get_arrays(filepath):
    with uproot.open(filepath) as f:
        tree = f["ABCD_scores"]
        return tree["scores"].array(library="np"), tree["weights"].array(library="np")


def rinv_token(rinv: float) -> str:
    return "0" if rinv == 0 else str(rinv).replace(".", "p")


def discover_mdark_scan_files():
    """One (mMed, mDark) -> file path. Both suffixes are identical duplicates
    for every point (verified directly); "nominal" is preferred deterministically."""
    pattern = os.path.join(
        SIG_DIR,
        f"ABCD_scores_nominal_t-channel_mMed-*_mDark-*_rinv-{rinv_token(FIXED_RINV)}_alpha-peak_yukawa-1_*.root",
    )
    all_files = glob.glob(pattern)
    name_re = re.compile(
        r"mMed-(\d+)_mDark-(\d+)_rinv-([\dp]+)_alpha-(\w+)_yukawa-([\dp]+)_(nominal|lundWeightNom)\.root$"
    )
    by_point = {}
    for f in all_files:
        m = name_re.search(os.path.basename(f))
        if not m:
            continue
        point = (int(m.group(1)), int(m.group(2)))
        suffix = m.group(6)
        if point not in by_point or suffix == "nominal":
            by_point[point] = f
    return by_point


def fingerprint(*parts) -> str:
    joined = "|".join(map(str, parts))
    return hashlib.sha1(joined.encode("utf-8")).hexdigest()[:16]


def weighted_auc(sig_scores, sig_weights, bkg_scores, bkg_weights, rng):
    sig_idx = rng.choice(len(sig_scores), min(len(sig_scores), AUC_CAP), replace=False) if len(sig_scores) else np.array([], dtype=int)
    bkg_idx = rng.choice(len(bkg_scores), min(len(bkg_scores), AUC_CAP), replace=False) if len(bkg_scores) else np.array([], dtype=int)

    scores = np.concatenate([sig_scores[sig_idx], bkg_scores[bkg_idx]])
    weights = np.concatenate([sig_weights[sig_idx], bkg_weights[bkg_idx]])
    labels = np.concatenate([np.ones(len(sig_idx)), np.zeros(len(bkg_idx))])
    weights = np.clip(weights, 0.0, None)

    auc = float(roc_auc_score(labels, scores, sample_weight=weights))
    return auc if auc >= 0.5 else 1.0 - auc


def compute() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)

    print("Discovering DNN mDark-scan signal score files …")
    points = discover_mdark_scan_files()
    print(f"  {len(points)} (mMed, mDark) points found")

    print("Loading QCD and ttbar background scores …")
    backgrounds = {}
    for name, path in BACKGROUND_FILES.items():
        sc, wt = get_arrays(path)
        backgrounds[name] = (sc, wt)
        print(f"  {name}: {len(sc):,} events")

    rows = {"qcd": [], "ttbar": []}
    for i, (point, path) in enumerate(sorted(points.items()), start=1):
        mmed, mdark = point
        sig_scores, sig_weights = get_arrays(path)
        rng_seed = SEED + int(hashlib.sha1(f"{mmed}_{mdark}".encode()).hexdigest()[:8], 16)
        for bkg_name, (bkg_scores, bkg_weights) in backgrounds.items():
            rng = np.random.default_rng(rng_seed)
            auc = weighted_auc(sig_scores, sig_weights, bkg_scores, bkg_weights, rng)
            rows[bkg_name].append({"mMed": float(mmed), "mDark": float(mdark), "auc": auc})
        print(f"  [{i:03d}/{len(points):03d}] mMed={mmed} mDark={mdark}  "
              f"AUC(qcd)={rows['qcd'][-1]['auc']:.4f}  AUC(ttbar)={rows['ttbar'][-1]['auc']:.4f}")

    for bkg_name, csv_path in (("qcd", QCD_CSV), ("ttbar", TTBAR_CSV)):
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["mMed", "mDark", "auc"])
            writer.writeheader()
            writer.writerows(rows[bkg_name])
        print(f"[OK] wrote {csv_path}")


def load_scan_csv(csv_path):
    rows = []
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            rows.append({"mMed": float(row["mMed"]), "mDark": float(row["mDark"]), "auc": float(row["auc"])})
    mMed_vals = sorted({r["mMed"] for r in rows})
    mdark_vals = sorted({r["mDark"] for r in rows})
    cells = {(r["mMed"], r["mDark"]): r["auc"] for r in rows}
    raw = [[cells.get((m, d), float("nan")) for d in mdark_vals] for m in mMed_vals]
    mMed_vals = [int(m) if float(m).is_integer() else m for m in mMed_vals]
    return mMed_vals, mdark_vals, raw


def plot() -> None:
    os.makedirs(PLOT_DIR, exist_ok=True)

    for bkg_name, csv_path, z_label, out_name in (
        ("qcd", QCD_CSV, "Event DNN AUC vs. QCD", "dnn_auc_grid_mdark_vs_mmed"),
        ("ttbar", TTBAR_CSV, "Event DNN AUC vs. ttbar", "dnn_auc_grid_mdark_vs_mmed_ttbar"),
    ):
        if not os.path.isfile(csv_path):
            raise SystemExit(f"{csv_path} not found -- run with --compute first.")
        mMed_vals, mdark_vals, raw = load_scan_csv(csv_path)
        auc_plot.plot_auc_grid(
            mMed_vals=mMed_vals, rinv_vals=mdark_vals, raw=raw,
            z_label=z_label,
            out_prefix=os.path.join(PLOT_DIR, out_name),
            y_label=r"$m_{\mathrm{dark}}\ [\mathrm{GeV}]$",
            y_tick_fmt=lambda v: f"{v:g}",
        )
    print(f"[OK] wrote heat maps under {PLOT_DIR}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--compute", action="store_true", help="(re)compute the AUC CSVs from the score files")
    args = ap.parse_args()
    if args.compute:
        compute()
    plot()


if __name__ == "__main__":
    main()
