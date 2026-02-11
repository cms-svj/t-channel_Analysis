#!/usr/bin/env python3
"""
trimQCD_skims.py

Trim QCD skim ROOT files by a MET upper cut, write trimmed ROOT files LOCALLY
(so you can move to EOS later), and make PRE/POST MET plots.

Key design choices (to avoid ROOT segfault/OOM you hit):
  - Read + write with uproot in CHUNKS (streaming) so memory stays flat.
  - Drop only branches matching --exclude regex(es) (default: none).
    For your broken Nano-style leaflists, use:  --exclude '^.*/\.f'
    (this drops GenJets_/.f*, Electrons_/.f*, etc).

Example:
  python3 trimQCD_skims.py --met-branch MET --plot-dir ./QCD_plots/trimming_skims \
    --ymin 1e-4 --overwrite --exclude '^.*/\.f' --year 2016

Optional:
  --keepPtBins QCD_Pt_470to600 QCD_Pt_600to800
"""

import os
import re
import sys
import argparse
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional

import numpy as np
import uproot
import awkward as ak

import matplotlib.pyplot as plt
import mplhep as hep
from matplotlib.backends.backend_pdf import PdfPages

plt.style.use(hep.style.CMS)

# =============================================================================
# USER CONFIG: dataset file lists + per-ptbin MET cuts
# =============================================================================

# Add/extend years/bins as needed. These are the ones from your logs.
from typing import Dict, List

SKIM_FILES: Dict[str, Dict[str, List[str]]] = {
    "2016": {
        "QCD_Pt_470to600": [
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_470to600/part-0.root",
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_470to600/part-1.root",
        ],
        "QCD_Pt_600to800": [
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-0.root",
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-1.root",
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-2.root",
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-3.root",
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-4.root",
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-5.root",
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-6.root",
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2016/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-7.root",
        ],
    },

    "2017": {
        "QCD_Pt_470to600": [
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2017/t_channel_pre_selection/nominal//QCD_Pt_470to600/part-0.root",
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2017/t_channel_pre_selection/nominal//QCD_Pt_470to600/part-1.root",
        ],
        "QCD_Pt_600to800": [
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2017/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-0.root",
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2017/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-1.root",
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2017/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-2.root",
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2017/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-3.root",
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2017/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-4.root",
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2017/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-5.root",
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2017/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-6.root",
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2017/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-7.root",
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2017/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-8.root",
        ],
    },

    "2018": {
        "QCD_Pt_470to600": [
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_470to600/part-0.root",
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_470to600/part-1.root",
        ],
        "QCD_Pt_600to800": [
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-0.root",
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-1.root",
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-2.root",
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-3.root",
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-4.root",
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-5.root",
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-6.root",
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-7.root",
            "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-8.root",
        ],
    },
}


# Per pt-bin MET max cut (from your printouts)
MET_MAX_CUT: Dict[str, float] = {
    "QCD_Pt_470to600": 750.0,
    "QCD_Pt_600to800": 900.0,
    # Add others as needed...
}

# =============================================================================
# Core helpers
# =============================================================================

@dataclass
class HistConfig:
    xmin: float
    xmax: float
    nbins: int

def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)

def compile_excludes(patterns: List[str]) -> List[re.Pattern]:
    return [re.compile(p) for p in patterns]

def keep_branch(name: str, exclude_res: List[re.Pattern]) -> bool:
    return not any(r.search(name) for r in exclude_res)

def build_hist_from_files(
    files: List[str],
    tree_name: str,
    met_branch: str,
    hist_cfg: HistConfig,
    step_size: int,
    met_max_for_count: Optional[float] = None,
) -> Tuple[np.ndarray, np.ndarray, int, int]:
    """
    Stream MET from many files and return (counts, bin_edges, total_events, events_above_cut).
    Uses uproot so it doesn't trigger the ROOT leaflist crash.
    """
    bins = np.linspace(hist_cfg.xmin, hist_cfg.xmax, hist_cfg.nbins + 1)
    counts = np.zeros(hist_cfg.nbins, dtype=np.float64)
    total = 0
    above = 0

    for f in files:
        with uproot.open(f) as fin:
            tree = fin[tree_name]
            if met_branch not in tree.keys():
                raise RuntimeError(f"MET branch '{met_branch}' not found in {f}")

            for arrays in tree.iterate(
                expressions=[met_branch],
                step_size=step_size,
                library="np",
            ):
                met = arrays[met_branch]
                total += met.shape[0]
                if met_max_for_count is not None:
                    above += int(np.sum(met > met_max_for_count))
                h, _ = np.histogram(met, bins=bins)
                counts += h

    return counts, bins, total, above
def trim_root_file_uproot(
    infile: str,
    outfile: str,
    tree_name: str,
    met_branch: str,
    met_max: float,
    exclude_patterns: List[str],
    step_size: int,
    overwrite: bool,
) -> Tuple[int, int]:

    if (not overwrite) and os.path.exists(outfile):
        with uproot.open(outfile) as f:
            t = f[tree_name]
            kept_existing = int(t.num_entries)
        return kept_existing, 0

    exclude_res = compile_excludes(exclude_patterns)
    ensure_dir(os.path.dirname(outfile))

    # If an old empty/bad file exists, remove it first when overwriting
    if overwrite and os.path.exists(outfile):
        os.remove(outfile)

    with uproot.open(infile) as fin:
        tree = fin[tree_name]
        all_branches = list(tree.keys())

        branches = [b for b in all_branches if keep_branch(b, exclude_res)]

        # Ensure MET exists for filtering
        if met_branch not in branches:
            if met_branch in all_branches:
                branches.append(met_branch)
            else:
                raise RuntimeError(f"MET branch '{met_branch}' not found in {infile}")

        kept_total = 0
        dropped_total = 0
        wrote_anything = False

        with uproot.recreate(outfile) as fout:
            for arrays in tree.iterate(
                expressions=branches,
                step_size=step_size,
                library="ak",
            ):
                met = arrays[met_branch]
                mask = met <= met_max

                n_chunk = int(len(met))
                n_keep = int(ak.sum(mask))
                n_drop = n_chunk - n_keep

                kept_total += n_keep
                dropped_total += n_drop

                if n_keep == 0:
                    continue

                out_chunk = {k: arrays[k][mask] for k in arrays.fields}

                if not wrote_anything:
                    # Create the tree with inferred types on first write
                    fout[tree_name] = out_chunk
                    wrote_anything = True
                else:
                    fout[tree_name].extend(out_chunk)

        # If nothing passed the cut, remove empty output file (prevents POST crash)
        if not wrote_anything:
            os.remove(outfile)

    return kept_total, dropped_total

def plot_hist_pair(
    pdf: PdfPages,
    title: str,
    bins: np.ndarray,
    pre_counts: np.ndarray,
    post_counts: Optional[np.ndarray],
    xlabel: str,
    ymin: Optional[float],
    logy: bool = True,
) -> None:
    fig, ax = plt.subplots(figsize=(12, 9))

    centers = 0.5 * (bins[:-1] + bins[1:])
    width = (bins[1] - bins[0])

    # totals
    n_pre  = int(np.sum(pre_counts))
    n_post = int(np.sum(post_counts)) if post_counts is not None else None

    # labels with counts
    pre_label = f"PRE (input), N = {n_pre:,}"
    post_label = f"POST (trimmed), N = {n_post:,}" if n_post is not None else None

    ax.step(
        centers,
        pre_counts,
        where="mid",
        linewidth=2.0,
        label=pre_label,
        color="C0",
    )

    if post_counts is not None:
        ax.step(
            centers,
            post_counts,
            where="mid",
            linewidth=2.0,
            label=post_label,
            color="C1",
        )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(f"Events / {width:.1f} GeV")
    ax.set_title(title)

    if logy:
        ax.set_yscale("log")
    if ymin is not None:
        ax.set_ylim(bottom=ymin)

    ax.set_xlim(bins[0], bins[-1])
    ax.legend(frameon=False)
    #hep.cms.label("Work in progress", ax=ax)

    pdf.savefig(fig)
    plt.close(fig)


# =============================================================================
# Main
# =============================================================================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--year", default="2018", help="Year key in SKIM_FILES (e.g. 2016)")
    ap.add_argument("--tree", default="Events", help="TTree name")
    ap.add_argument("--met-branch", default="MET", help="MET branch name (e.g. MET or MET_pt)")
    ap.add_argument("--plot-dir", default="./QCD_plots/trimming_skims", help="Directory for plots + outputs")
    ap.add_argument("--trim-subdir", default="trimmed_root", help="Subdir under plot-dir for trimmed ROOT files")
    ap.add_argument("--step-size", type=int, default=50_000, help="Chunk size for uproot.iterate")
    ap.add_argument("--xmin", type=float, default=0.0)
    ap.add_argument("--xmax", type=float, default=4000.0)
    ap.add_argument("--nbins", type=int, default=120)
    ap.add_argument("--ymin", type=float, default=None, help="Force y-axis min (e.g. 1e-4)")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite local outputs")
    ap.add_argument("--exclude", action="append", default=[],
                    help="Regex for branches to exclude. Repeatable. "
                         "For broken leaflists use: --exclude '^.*/\\.f'")
    ap.add_argument("--keepPtBins", nargs="*", default=None,
                    help="Only process these pt bins (e.g. QCD_Pt_470to600 QCD_Pt_600to800)")
    ap.add_argument("--plots-only",action="store_true",help="Only build plots from existing trimmed ROOT files (no trimming)")

    args = ap.parse_args()

    if args.year not in SKIM_FILES:
        print(f"[ERROR] year '{args.year}' not found in SKIM_FILES. Available: {list(SKIM_FILES.keys())}")
        sys.exit(1)

    ensure_dir(args.plot_dir)
    trimmed_root_base = os.path.join(args.plot_dir, args.trim_subdir)
    ensure_dir(trimmed_root_base)

    hist_cfg = HistConfig(args.xmin, args.xmax, args.nbins)

    # Select pt bins
    ptbins = list(SKIM_FILES[args.year].keys())
    if args.keepPtBins is not None and len(args.keepPtBins) > 0:
        ptbins = [p for p in ptbins if p in set(args.keepPtBins)]

    # PDF with all plots for that year
    pdf_path = os.path.join(args.plot_dir, f"trimQCD_{args.year}.pdf")

    print("\n=== LOCAL trim + plots (uproot streaming) ===")
    print(f"Year         : {args.year}")
    print(f"Plot dir     : {args.plot_dir}")
    print(f"Trimmed root : {trimmed_root_base}")
    print(f"Tree         : {args.tree}")
    print(f"MET branch   : {args.met_branch}")
    print(f"Step size    : {args.step_size}")
    print(f"Hist range   : [{args.xmin}, {args.xmax}] with {args.nbins} bins")
    print(f"Exclude      : {args.exclude}\n")

    grand_kept = 0
    grand_dropped = 0

    with PdfPages(pdf_path) as pdf:
        for pt in ptbins:
            if pt not in MET_MAX_CUT:
                print(f"[WARN] No MET cut configured for {pt} in MET_MAX_CUT. Skipping.")
                continue

            files = SKIM_FILES[args.year][pt]
            metmax = MET_MAX_CUT[pt]

            tag = f"{args.year}_{pt}"
            out_dir = os.path.join(trimmed_root_base, tag)
            ensure_dir(out_dir)

            print(f"--- {tag}  (drop {args.met_branch} > {metmax}) ---")
            print("  Building PRE histogram (remote inputs)...")

            pre_counts, bins, total_pre, above_pre = build_hist_from_files(
                files=files,
                tree_name=args.tree,
                met_branch=args.met_branch,
                hist_cfg=hist_cfg,
                step_size=args.step_size,
                met_max_for_count=metmax,
            )
            print(f"    Total events (pre): {total_pre}")
            print(f"    Events above cut (pre): {above_pre}\n")

            # Trim each file
            kept_dataset = 0
            dropped_dataset = 0
            out_files_local = []

            for f in files:
                base = os.path.basename(f)
                out_f = os.path.join(out_dir, base)
                out_files_local.append(out_f)

                print(f"  Input : {f}")
                print(f"  Out   : {out_f}")

                try:
                    kept, dropped = trim_root_file_uproot(
                        infile=f,
                        outfile=out_f,
                        tree_name=args.tree,
                        met_branch=args.met_branch,
                        met_max=metmax,
                        exclude_patterns=args.exclude,
                        step_size=args.step_size,
                        overwrite=args.overwrite,
                    )
                    kept_dataset += kept
                    dropped_dataset += dropped
                    print(f"    kept={kept}  dropped={dropped}\n")
                except Exception as e:
                    print(f"[ERROR] Failed on {f}\n  -> {e}\n")

            grand_kept += kept_dataset
            grand_dropped += dropped_dataset

            # POST histogram from local outputs that exist
            good_outs = []
            for x in out_files_local:
                if not os.path.exists(x):
                    continue
                try:
                    with uproot.open(x) as f:
                        if args.tree in f:
                            good_outs.append(x)
                except Exception:
                    pass

            post_counts = None
            if len(good_outs) == 0:
                print("  [WARN] No outputs written for this dataset. Skipping POST + plot.\n")
                plot_hist_pair(
                    pdf=pdf,
                    title=f"{tag}: PRE only (no POST written)",
                    bins=bins,
                    pre_counts=pre_counts,
                    post_counts=None,
                    xlabel=args.met_branch,
                    ymin=args.ymin,
                    logy=True,
                )
                continue

            print("  Building POST histogram (local trimmed outputs)...")
            post_counts, _, total_post, above_post = build_hist_from_files(
                files=good_outs,
                tree_name=args.tree,
                met_branch=args.met_branch,
                hist_cfg=hist_cfg,
                step_size=args.step_size,
                met_max_for_count=metmax,  # should be 0 above cut ideally
            )
            print(f"    Total events (post): {total_post}")
            print(f"    Events above cut (post): {above_post} (should be ~0)\n")

            plot_hist_pair(
                pdf=pdf,
                title=f"{tag}: MET (cut {metmax})",
                bins=bins,
                pre_counts=pre_counts,
                post_counts=post_counts,
                xlabel=args.met_branch,
                ymin=args.ymin,
                logy=True,
            )

    print("=== Done ===")
    print(f"PDF saved: {pdf_path}")
    print(f"Grand kept   : {grand_kept}")
    print(f"Grand dropped: {grand_dropped}")

if __name__ == "__main__":
    main()
