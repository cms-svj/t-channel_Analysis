#!/usr/bin/env python3
"""
Extra ParticleNet AUC grids, in addition to the existing mMed-vs-rinv grid
in makeAUCtables_pnet.py:

  1. mDark vs mMed   (rinv=0.3, yukawa=1 fixed) vs QCD
  2. yukawa vs mMed  (rinv=0.3, mDark=20 fixed) vs QCD
  3. mDark vs mMed   (rinv=0.3, yukawa=1 fixed) vs ttbar

Background: 2018 QCD/TTJets from the same base Figure2_makerusingskims.py
uses,
  /store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/data_mc_WNAE
  (year)/t_channel_pre_selection/nominal/{QCD_Pt_*,TTJets*}
Signal: 2018 t-channel skims from
  /store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto
  (year)/t_channel_pre_selection/nominal/<signal>

Weighted AUC (sklearn roc_auc_score sample_weight), since the pNet score is
correlated with jet substructure and the signal MC needs the Lund-plane
correction to model that substructure correctly:
  QCD event weight    = Weight * puWeight * NonPrefiringProb
  Signal event weight = Weight * puWeight * NonPrefiringProb * lundWeightNom
                         * (CutFlow/Initial / CutFlow/InitialLundNominal)
matching the weight formula and Initial/InitialLundNominal renormalization
convention used throughout Figure2_makerusingskims.py.

Two-step, mirroring the other AUC/signal-efficiency scripts in this repo:
  --compute   Slow. Streams ParticleNet scores + weights from EOS, computes
              weighted AUC per grid point. Writes mdark_scan_auc.csv and
              yukawa_scan_auc.csv.
  (default)   Fast. Reads those CSVs and draws the two heat maps.

Usage:
  source condor/initCondor.sh
  python3 makeAUCtables_pnet_extra_axes.py --compute
  python3 makeAUCtables_pnet_extra_axes.py
"""
import argparse
import csv
import hashlib
import os
import re
import sys

import awkward as ak
import numpy as np
import uproot
from sklearn.metrics import roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import calculateAUCtables_pnetPOSTtrimming as auc_calc
import makeAUCtables_pnet as auc_plot

OUT_DIR = os.path.join(HERE, "pnet_auc_extra_axes")
MDARK_CSV = os.path.join(OUT_DIR, "mdark_scan_auc.csv")
YUKAWA_CSV = os.path.join(OUT_DIR, "yukawa_scan_auc.csv")
MDARK_TTBAR_CSV = os.path.join(OUT_DIR, "mdark_scan_auc_ttbar.csv")
PLOT_DIR = os.path.join(
    HERE, "SVJ-tchannel-Run2_site", "supplemetry", "assets_v2", "auc"
)

SIGNAL_EOS_BASE = "/store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/2018/t_channel_pre_selection/nominal"
BACKGROUND_EOS_BASE = "/store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/data_mc_WNAE/2018/t_channel_pre_selection/nominal"
QCD_EOS_BASE = BACKGROUND_EOS_BASE
BACKGROUND_PREFIXES = {
    "qcd": ("QCD_Pt_",),
    "ttbar": ("TTJets",),
}

PNET_BRANCH = "JetsAK8_pNetJetTaggerScore"
TREE = "Events"
INVALID_SCORE_CUT = -9.0
STEP_SIZE = 50_000
AUC_CAP = 500_000
SEED = 42

FIXED_RINV = 0.3
FIXED_YUKAWA_FOR_MDARK_SCAN = 1.0
FIXED_MDARK_FOR_YUKAWA_SCAN = 20.0
ALPHA = "peak"


def discover_axis_scan(signal_base, *, scan="mdark"):
    """Discover signal points varying mDark (yukawa=1 fixed) or yukawa (mDark=20 fixed), rinv=0.3 fixed."""
    samples = []
    for entry in auc_calc.xrdfs_ls(signal_base):
        meta = auc_calc.parse_signal_metadata(entry)
        if meta is None:
            continue
        if not np.isclose(float(meta["rinv"]), FIXED_RINV, atol=1e-10):
            continue
        if str(meta["alpha"]) != ALPHA:
            continue
        if scan == "mdark":
            if not np.isclose(float(meta["yukawa"]), FIXED_YUKAWA_FOR_MDARK_SCAN, atol=1e-10):
                continue
        else:
            if not np.isclose(float(meta["mDark"]), FIXED_MDARK_FOR_YUKAWA_SCAN, atol=1e-10):
                continue

        files = tuple(auc_calc.discover_root_files(entry))
        if not files:
            print(f"[WARN] no ROOT files for {entry}")
            continue
        samples.append(
            auc_calc.SignalSample(name=str(meta["signal_name"]), eos_dir=entry, files=files, metadata=meta)
        )
    return sorted(samples, key=lambda s: (float(s.metadata["mMed"]), float(s.metadata["mDark"]), float(s.metadata["yukawa"])))


def discover_background_samples(background_base, prefixes):
    samples = []
    for entry in auc_calc.xrdfs_ls(background_base):
        name = os.path.basename(entry.rstrip("/"))
        if not name.startswith(prefixes):
            continue
        files = tuple(auc_calc.discover_root_files(entry))
        if not files:
            print(f"[WARN] no ROOT files for {entry}")
            continue
        samples.append((name, entry, files))
    return samples


def lund_norm_factor(fin) -> float:
    """Initial / InitialLundNominal, the same renormalization used in Figure2_makerusingskims.py
    so that CutFlow/Initial remains valid regardless of the Lund-plane reweighting."""
    cutflow = fin["CutFlow"]
    initial = float(cutflow["Initial"].array(library="np")[0])
    initial_lund = float(cutflow["InitialLundNominal"].array(library="np")[0])
    if not np.isfinite(initial) or not np.isfinite(initial_lund) or initial_lund == 0.0:
        return 1.0
    return initial / initial_lund


def fingerprint(*parts) -> str:
    joined = "|".join(map(str, parts))
    return hashlib.sha1(joined.encode("utf-8")).hexdigest()[:16]


def safe_filename(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text)


def load_scores_and_weights(files, *, label, is_signal, cache_path, cache_key):
    """Per-event max valid pNet score, and the matching physics weight
    (Weight*puWeight*NonPrefiringProb, times lundWeightNom*Initial/InitialLundNominal for signal)."""
    if os.path.exists(cache_path):
        with np.load(cache_path, allow_pickle=False) as cache:
            if str(cache["cache_key"].item()) == cache_key:
                return np.asarray(cache["scores"]), np.asarray(cache["weights"])

    base_branches = [PNET_BRANCH, "Weight", "puWeight", "NonPrefiringProb"]
    branches = base_branches + (["lundWeightNom"] if is_signal else [])

    score_chunks, weight_chunks = [], []
    for file_index, file_path in enumerate(files, start=1):
        with uproot.open(file_path) as fin:
            norm = lund_norm_factor(fin) if is_signal else 1.0
            tree = fin[TREE]
            for arrays in tree.iterate(branches, step_size=STEP_SIZE, library="ak"):
                raw_scores = arrays[PNET_BRANCH]
                valid_scores = ak.mask(raw_scores, raw_scores >= INVALID_SCORE_CUT)
                event_scores = ak.to_numpy(
                    ak.fill_none(ak.max(valid_scores, axis=1, mask_identity=True), np.nan)
                ).astype(np.float64, copy=False)

                weight = (
                    ak.to_numpy(arrays["Weight"]).astype(np.float64)
                    * ak.to_numpy(arrays["puWeight"]).astype(np.float64)
                    * ak.to_numpy(arrays["NonPrefiringProb"]).astype(np.float64)
                )
                if is_signal:
                    weight = weight * ak.to_numpy(arrays["lundWeightNom"]).astype(np.float64) * norm

                score_chunks.append(event_scores)
                weight_chunks.append(weight)
        print(f"  {label}: {file_index}/{len(files)} files read")

    scores = np.concatenate(score_chunks) if score_chunks else np.empty(0)
    weights = np.concatenate(weight_chunks) if weight_chunks else np.empty(0)

    finite = np.isfinite(scores)
    scores, weights = scores[finite], weights[finite]

    np.savez_compressed(cache_path, scores=scores, weights=weights, cache_key=np.asarray(cache_key))
    return scores, weights


def weighted_auc(sig_scores, sig_weights, bkg_scores, bkg_weights, rng):
    sig_idx = rng.choice(len(sig_scores), min(len(sig_scores), AUC_CAP), replace=False) if len(sig_scores) else np.array([], dtype=int)
    bkg_idx = rng.choice(len(bkg_scores), min(len(bkg_scores), AUC_CAP), replace=False) if len(bkg_scores) else np.array([], dtype=int)

    scores = np.concatenate([sig_scores[sig_idx], bkg_scores[bkg_idx]])
    weights = np.concatenate([sig_weights[sig_idx], bkg_weights[bkg_idx]])
    labels = np.concatenate([np.ones(len(sig_idx)), np.zeros(len(bkg_idx))])

    # Negative weights (rare MC negative-weight events) break roc_auc_score's
    # sample_weight handling; clip like the rest of the analysis treats them.
    weights = np.clip(weights, 0.0, None)
    return float(roc_auc_score(labels, scores, sample_weight=weights))


def compute_for_background(background_name, scan_configs) -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    prefixes = BACKGROUND_PREFIXES[background_name]

    print(f"Discovering 2018 {background_name} background from {BACKGROUND_EOS_BASE} (prefixes={prefixes}) …")
    bkg_samples = discover_background_samples(BACKGROUND_EOS_BASE, prefixes)
    if not bkg_samples:
        raise RuntimeError(f"No {background_name} ROOT files discovered.")
    print(f"  {len(bkg_samples)} {background_name} samples, {sum(len(files) for _, _, files in bkg_samples)} files")

    bkg_score_chunks, bkg_weight_chunks = [], []
    for name, eos_dir, files in bkg_samples:
        key = fingerprint(auc_calc.normalise_eos_path(eos_dir), tuple(files), PNET_BRANCH, INVALID_SCORE_CUT)
        cache = os.path.join(OUT_DIR, f"{background_name}_pnet_scores_weights_{safe_filename(name)}_{key}.npz")
        scores, weights = load_scores_and_weights(files, label=name, is_signal=False, cache_path=cache, cache_key=key)
        bkg_score_chunks.append(scores)
        bkg_weight_chunks.append(weights)
        print(f"  {background_name} ready: {name:<20} events={len(scores):,}")
    combined_bkg_scores = np.concatenate(bkg_score_chunks)
    combined_bkg_weights = np.concatenate(bkg_weight_chunks)

    for scan, fixed_label, csv_path, extra_field in scan_configs:
        print(f"\n=== {background_name} / {scan} scan ({fixed_label}) ===")
        samples = discover_axis_scan(SIGNAL_EOS_BASE, scan=scan)
        print(f"Discovered {len(samples)} signal points for the {scan} scan")

        rows = []
        for i, sample in enumerate(samples, start=1):
            key = fingerprint(auc_calc.normalise_eos_path(sample.eos_dir), tuple(sample.files), PNET_BRANCH, INVALID_SCORE_CUT)
            cache = os.path.join(OUT_DIR, f"signal_pnet_scores_weights_{safe_filename(sample.name)}_{key}.npz")
            sig_scores, sig_weights = load_scores_and_weights(
                sample.files, label=sample.name, is_signal=True, cache_path=cache, cache_key=key
            )
            rng = np.random.default_rng(SEED + int(hashlib.sha1(sample.name.encode()).hexdigest()[:8], 16))
            auc = weighted_auc(sig_scores, sig_weights, combined_bkg_scores, combined_bkg_weights, rng)
            rows.append({
                "mMed": float(sample.metadata["mMed"]),
                extra_field: float(sample.metadata[extra_field]),
                "auc": auc,
            })
            print(
                f"  [{i:03d}/{len(samples):03d}] mMed={sample.metadata['mMed']:.0f} "
                f"{extra_field}={sample.metadata[extra_field]:g}  AUC={auc:.4f}"
            )

        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["mMed", extra_field, "auc"])
            writer.writeheader()
            writer.writerows(rows)
        print(f"[OK] wrote {csv_path}")


def compute() -> None:
    compute_for_background("qcd", [
        ("mdark", f"rinv={FIXED_RINV}, yukawa={FIXED_YUKAWA_FOR_MDARK_SCAN}", MDARK_CSV, "mDark"),
        ("yukawa", f"rinv={FIXED_RINV}, mDark={FIXED_MDARK_FOR_YUKAWA_SCAN}", YUKAWA_CSV, "yukawa"),
    ])
    compute_for_background("ttbar", [
        ("mdark", f"rinv={FIXED_RINV}, yukawa={FIXED_YUKAWA_FOR_MDARK_SCAN}", MDARK_TTBAR_CSV, "mDark"),
    ])


def load_scan_csv(csv_path, axis_field):
    rows = []
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            rows.append({"mMed": float(row["mMed"]), axis_field: float(row[axis_field]), "auc": float(row["auc"])})
    mMed_vals = sorted({r["mMed"] for r in rows})
    axis_vals = sorted({r[axis_field] for r in rows})
    cells = {(r["mMed"], r[axis_field]): r["auc"] for r in rows}
    raw = [[cells.get((m, a), float("nan")) for a in axis_vals] for m in mMed_vals]
    mMed_vals = [int(m) if float(m).is_integer() else m for m in mMed_vals]
    return mMed_vals, axis_vals, raw


def plot() -> None:
    os.makedirs(PLOT_DIR, exist_ok=True)

    if not os.path.isfile(MDARK_CSV):
        raise SystemExit(f"{MDARK_CSV} not found -- run with --compute first.")
    mMed_vals, mdark_vals, raw = load_scan_csv(MDARK_CSV, "mDark")
    auc_plot.plot_auc_grid(
        mMed_vals=mMed_vals, rinv_vals=mdark_vals, raw=raw,
        z_label="ParticleNet AUC vs. QCD",
        out_prefix=os.path.join(PLOT_DIR, "pnet_auc_grid_mdark_vs_mmed"),
        y_label=r"$m_{\mathrm{dark}}\ [\mathrm{GeV}]$",
        y_tick_fmt=lambda v: f"{v:g}",
    )

    if not os.path.isfile(YUKAWA_CSV):
        raise SystemExit(f"{YUKAWA_CSV} not found -- run with --compute first.")
    mMed_vals, yukawa_vals, raw = load_scan_csv(YUKAWA_CSV, "yukawa")
    auc_plot.plot_auc_grid(
        mMed_vals=mMed_vals, rinv_vals=yukawa_vals, raw=raw,
        z_label="ParticleNet AUC vs. QCD",
        out_prefix=os.path.join(PLOT_DIR, "pnet_auc_grid_yukawa_vs_mmed"),
        y_label=r"$\lambda$",
        y_tick_fmt=lambda v: f"{v:g}",
    )

    if not os.path.isfile(MDARK_TTBAR_CSV):
        raise SystemExit(f"{MDARK_TTBAR_CSV} not found -- run with --compute first.")
    mMed_vals, mdark_vals, raw = load_scan_csv(MDARK_TTBAR_CSV, "mDark")
    auc_plot.plot_auc_grid(
        mMed_vals=mMed_vals, rinv_vals=mdark_vals, raw=raw,
        z_label="ParticleNet AUC vs. ttbar",
        out_prefix=os.path.join(PLOT_DIR, "pnet_auc_grid_mdark_vs_mmed_ttbar"),
        y_label=r"$m_{\mathrm{dark}}\ [\mathrm{GeV}]$",
        y_tick_fmt=lambda v: f"{v:g}",
    )
    print(f"[OK] wrote heat maps under {PLOT_DIR}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--compute", action="store_true", help="(re)compute the AUC CSVs from EOS")
    args = ap.parse_args()
    if args.compute:
        compute()
    plot()


if __name__ == "__main__":
    main()
