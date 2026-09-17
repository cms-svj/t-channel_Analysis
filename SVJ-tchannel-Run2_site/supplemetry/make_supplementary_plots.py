#!/usr/bin/env python3
"""
Single entry point for the supplementary-material plots: N-1, ParticleNet
inputs, WNAE inputs, DNN inputs, background-composition pies, and signal
efficiency maps. Histogram-style plots are drawn from one merged cache
(supplemetry_cache.root, built by this same script).

This does not reimplement any histogram-filling or drawing logic: it reuses
the existing, already-validated functions from Figure2_makerusingskims.py,
NMinusOne_maker_RA2.py, ParticleNetInput_makerusingskims.py,
EventDNNInput_makerusingskims.py and combinehistplotter.py exactly as they
run today, just pointed at namespaced subdirectories of one merged file
instead of four separate cache files scattered across the repo.

Usage:
  source condor/initCondor.sh
  python3 make_supplementary_plots.py                 # merge + draw
  python3 make_supplementary_plots.py --skip-merge     # reuse existing cache, just redraw
  python3 make_supplementary_plots.py --merge-only     # only rebuild the merged cache
"""
import argparse
import os
import subprocess
import sys
import types

import ROOT

ROOT.gROOT.SetBatch(True)

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(os.path.dirname(HERE))  # .../t-channel_Analysis
sys.path.insert(0, BASE)

import Figure2_makerusingskims as f2  # noqa: E402
import NMinusOne_maker_RA2 as n1_mod  # noqa: E402
import ParticleNetInput_makerusingskims as pnet_mod  # noqa: E402
import EventDNNInput_makerusingskims as dnn_mod  # noqa: E402
import combinehistplotter as pie_mod  # noqa: E402

CACHE = os.path.join(HERE, "supplemetry_cache.root")
ASSETS = os.path.join(HERE, "assets_v2")

# Each source cache keeps its own internal layout (raw/groups/<era>/<process>/
# <var>, raw/data/..., raw/signals/..., shapes/... etc.) untouched; it is
# copied wholesale under one top-level namespace per section (wnae/, dnn/,
# pnet/, nminus1/) so the four don't collide. Downstream drawing code
# (Figure2's draw_stack, etc.) takes a ROOT directory-like object and does
# `directory.Get("raw/groups/...")` -- a TDirectory supports that exactly
# like a TFile does, so pointing it at merged_cache.Get("wnae") works with
# zero changes to the existing per-section maker scripts.
CACHE_SOURCES = {
    "wnae": os.path.join(BASE, "Figure2_plots", "WNAE_alljets_dataMC_vars_Run2_wlundcorrection_v4.root"),
    "dnn": os.path.join(BASE, "supplementry_plots_postTrimandgapveto", "DNN", "EventDNNInput_plots_noWNAE", "event_dnn_inputs_v3.root"),
    "pnet": os.path.join(BASE, "ParticleNetInput_plots", "pnet_input_vars_wlundcorrection.root"),
    "nminus1": os.path.join(BASE, "supplementry_plots_postTrimandgapveto", "Nminus1plots", "nminus1_RA2_allbkgs_Run2_withSignals.root"),
}


def copy_directory(src_dir: "ROOT.TDirectory", dst_dir: "ROOT.TDirectory") -> int:
    count = 0
    for key in src_dir.GetListOfKeys():
        name = key.GetName()
        obj = key.ReadObj()
        if obj.InheritsFrom("TDirectory"):
            next_dst = dst_dir.GetDirectory(name)
            if not next_dst:
                next_dst = dst_dir.mkdir(name)
            count += copy_directory(obj, next_dst)
        elif obj.InheritsFrom("TH1"):
            dst_dir.cd()
            obj.Write(name, ROOT.TObject.kOverwrite)
            count += 1
    return count


def merge_caches() -> None:
    out = ROOT.TFile.Open(CACHE, "RECREATE")
    if not out or out.IsZombie():
        raise RuntimeError(f"Could not create {CACHE}")
    try:
        for label, path in CACHE_SOURCES.items():
            if not os.path.isfile(path):
                print(f"[WARN] Missing source for '{label}': {path}")
                continue
            src = ROOT.TFile.Open(path, "READ")
            if not src or src.IsZombie():
                print(f"[WARN] Could not open source for '{label}': {path}")
                continue
            top = out.mkdir(label)
            n = copy_directory(src, top)
            src.Close()
            print(f"[OK] {label}: copied {n} histograms from {path}")
    finally:
        out.Write()
        out.Close()
    print(f"[OK] wrote merged cache: {CACHE}")

# Variables actually referenced by supplementryv2.tex, per section.
WNAE_VARIABLES = [
    "allPtDAK8",
    "allSoftDropMassAK8",
    "allTau43AK8",
    "allTau32AK8",
    "allTau21AK8",
    "allEcfD2b1AK8",
    "allEcfC2b1AK8",
    "allEcfC2b2AK8",
    "allEcfN2b1AK8",
]

PNET_VARIABLES = ["del_r", "log_e_jete", "log_e", "log_pt_jetpt", "log_pt"]

DNN_VARIABLES = [
    "DeltaEta01GoodJetsAK8",
    "DeltaR01GoodJetsAK8",
    "LundJetPlaneZ01GoodJetsAK8",
    "GoodJetsAK80_deltaPhiMET",
    "GoodJetsAK81_deltaPhiMET",
    "GoodJetsAK80_LundJetPlaneZ",
    "GoodJetsAK81_LundJetPlaneZ",
    "GoodJetsAK80_MTMETLog",
    "GoodJetsAK81_MTMETLog",
]

N1_PLOTS = [
    ("h_j1PtAK8_pre__stcut", "j1PtAK8"),
    ("h_j2PtAK8_pre__stcut", "j2PtAK8"),
    ("h_jPtAK8_pre__stcut", "jPtAK8"),
    ("h_njetsAK8_pre__2jetsAK8", "njetsAK8"),
]


def draw_wnae_inputs(cache: ROOT.TFile, out_dir: str) -> None:
    section = cache.Get("wnae")
    os.makedirs(out_dir, exist_ok=True)
    for variable in WNAE_VARIABLES:
        f2.draw_stack(
            section,
            "Run2",
            variable,
            os.path.join(out_dir, variable),
            f2.DEFAULT_SIGNAL_DIRS,
            normalized=False,
            include_data=True,
            include_ratio=True,
            force_log_y=True,
        )


def draw_pnet_inputs(cache: ROOT.TFile, out_dir: str) -> None:
    pnet_mod.configure_figure2_globals()
    section = cache.Get("pnet")
    os.makedirs(out_dir, exist_ok=True)
    for variable in PNET_VARIABLES:
        f2.draw_stack(
            section,
            "Run2",
            variable,
            os.path.join(out_dir, variable),
            pnet_mod.PARTICLENET_SIGNAL_DIRS,
            normalized=False,
            include_data=True,
            include_ratio=True,
            force_log_y=True,
        )


def draw_dnn_inputs(cache: ROOT.TFile, out_dir: str) -> None:
    f2.VARIABLES = dnn_mod.figure2_variable_config()
    f2.SIGNAL_SCALE_FACTORS = dnn_mod.SIGNAL_SCALE_FACTORS
    section = cache.Get("dnn")
    os.makedirs(out_dir, exist_ok=True)
    for variable in DNN_VARIABLES:
        f2.draw_stack(
            section,
            "Run2",
            variable,
            os.path.join(out_dir, variable),
            dnn_mod.DEFAULT_SIGNAL_DIRS,
            normalized=False,
            include_data=True,
            include_ratio=True,
            force_log_y=True,
        )
        # figure2_variable_config()/draw_stack write "<variable>_log_raw_ratio.pdf";
        # the tex expects the "_wdata_ratio" naming already used elsewhere in
        # this document, so mirror that here.
        src = os.path.join(out_dir, f"{variable}_log_raw_ratio.pdf")
        dst = os.path.join(out_dir, f"{variable}_log_raw_wdata_ratio.pdf")
        if os.path.exists(src):
            os.replace(src, dst)
        src_png = os.path.join(out_dir, f"{variable}_log_raw_ratio.png")
        dst_png = os.path.join(out_dir, f"{variable}_log_raw_wdata_ratio.png")
        if os.path.exists(src_png):
            os.replace(src_png, dst_png)


def draw_nminus1(cache: ROOT.TFile, out_dir: str) -> None:
    section = cache.Get("nminus1")
    os.makedirs(out_dir, exist_ok=True)
    args = types.SimpleNamespace(
        normalized=True,
        signal_scale=1.0,
        lumi_pb=n1_mod.LUMI_PB,
        preliminary_label=False,
        signals=",".join(n1_mod.DEFAULT_SIGNAL_DIRS),
    )
    for plot_name, variable in N1_PLOTS:
        n1_mod.draw_plot(section, plot_name, variable, out_dir, args)


# Owned by the stat-inference step, not duplicated into our cache; both the
# pie chart and the signal region-A table read it directly.
COMBINE_FILE = (
    "/uscms/home/nparmar/nobackup/SVJ/stat_inference_unblind_gapJetveto/"
    "stat_histograms/MET_mMed-fullScan_run2_pNet_V5_DNN85_WP90_MET250_BD_lWN_"
    "allsyst_min2p_distortion_alternative_trigger_unblind_data_gapJetveto_0to3PSVJ.root"
)


def draw_background_composition(out_dir: str) -> None:
    os.makedirs(out_dir, exist_ok=True)
    import uproot

    with uproot.open(COMBINE_FILE) as f:
        years = pie_mod.detect_years(f)
    total_yields = {
        proc: {region: {svj: 0.0 for svj in pie_mod.SVJ_ORDER} for region in pie_mod.PLOT_REGIONS}
        for proc in pie_mod.BKG_ORDER
    }
    for year in years:
        with uproot.open(COMBINE_FILE) as f:
            for proc in pie_mod.BKG_ORDER:
                for region_plot in pie_mod.PLOT_REGIONS:
                    region_root = pie_mod.ROOT_REGION_FOR[region_plot]
                    for svj in pie_mod.SVJ_ORDER:
                        hist_path = f"{svj}Y{year}_Run2/{region_root}/{proc}"
                        if hist_path not in f:
                            continue
                        values, _ = f[hist_path].to_numpy()
                        total_yields[proc][region_plot][svj] += float(values.sum())
    pie_mod.write_fraction_pie_charts(
        "Run2_Combined",
        total_yields,
        out_dir,
        prelim=False,
        tagger_label="PN",
    )


def write_signal_acceptance_table(out_dir: str) -> None:
    os.makedirs(out_dir, exist_ok=True)
    import uproot

    with uproot.open(COMBINE_FILE) as f:
        years = pie_mod.detect_years(f)
    acceptance = pie_mod.compute_signal_region_a_acceptance(COMBINE_FILE, years)
    pie_mod.write_signal_region_a_table(
        acceptance, os.path.join(out_dir, "signal_region_a_acceptance.tex")
    )


def combine_name_to_skim_dir(name: str) -> str:
    """'mMed600_mDark20_rinv0p3_yukawa1' -> the skim-directory name
    NMinusOne_maker_RA2.py/Figure2_makerusingskims.py use to find the
    corresponding signal sample on EOS."""
    import re

    m = re.match(r"mMed(\d+)_mDark(\d+)_rinv([\dp]+)_yukawa([\dp]+)", name)
    if not m:
        raise ValueError(f"Unrecognized signal name format: {name}")
    mmed, mdark, rinv, yukawa = m.groups()
    return f"t-channel_mMed-{mmed}_mDark-{mdark}_rinv-{rinv}_alpha-peak_yukawa-{yukawa}"


def compute_signal_generated_totals(signals, years=("2016", "2017", "2018")) -> dict:
    """
    Run2-combined, luminosity-scaled generated (pre-selection) yield for each
    signal, read from CutFlow/Initial in each sample's skim files -- the same
    quantity NMinusOne_maker_RA2.lund_nominal_norm reads for its Lund-plane
    correction. This only opens the small CutFlow tree per file (one entry),
    not the full Events tree, so it's cheap even though it touches EOS.
    """
    import uproot

    totals = {}
    for signal_name in signals:
        skim_dir = combine_name_to_skim_dir(signal_name)
        total = 0.0
        for year in years:
            try:
                entries = f2.get_signal_directories(f2.DEFAULT_SIGNAL_BASE, year, [skim_dir], strict=False)
            except RuntimeError:
                entries = []
            if not entries:
                print(f"[WARN] {signal_name}/{year}: signal directory not found, skipping")
                continue
            _, sample_dir = entries[0]
            try:
                files = [e for e in f2.xrdfs_ls(sample_dir, recursive=True) if e.endswith(".root")]
            except RuntimeError as exc:
                print(f"[WARN] {signal_name}/{year}: could not list files ({exc})")
                continue
            year_initial = 0.0
            for path in files:
                url = f2.EOS_HOST + "/" + f2.normalize_eos_path(path)
                try:
                    with uproot.open(url) as fin:
                        year_initial += float(fin["CutFlow"]["Initial"].array(library="np")[0])
                except Exception as exc:
                    print(f"[WARN] {signal_name}/{year}: could not read CutFlow from {path} ({exc})")
            total += year_initial * f2.LUMI_PB[year]
        totals[signal_name] = total
    return totals


def write_signal_acceptance_vs_generated_table(out_dir: str) -> None:
    os.makedirs(out_dir, exist_ok=True)
    import uproot

    signals = [name for name, _ in pie_mod.REFERENCE_SIGNALS]
    with uproot.open(COMBINE_FILE) as f:
        years = pie_mod.detect_years(f)
    region_a_yields = pie_mod.compute_signal_region_yields(COMBINE_FILE, years, signals)
    generated = compute_signal_generated_totals(signals, years)
    pie_mod.write_signal_region_a_vs_generated_table(
        region_a_yields, generated, os.path.join(out_dir, "signal_region_a_vs_generated.tex")
    )


# --- postfits, limits, AUC, tagger comparisons -----------------------------
#
# These four sections aren't drawn from supplemetry_cache.root: postfits come
# from combine's FitDiagnostics output (a full re-fit is out of scope for a
# lightweight redraw), limits come from a CSV of already-computed combine
# scan points, the ttbar AUC grids are hardcoded arrays checked directly into
# makeAUCtables_pnet.py, and the tagger-comparison plots come from their own
# EOS-built ROOT cache. Each already has its own stable, version-controlled
# source, so rather than force them into supplemetry_cache.root's histogram
# shape, this just points at those sources directly -- one script, but not
# one physical file, for the parts that aren't natively histogram-shaped.

POSTFIT_SRC = os.path.join(BASE, "t-channel_plotting_scripts", "postfits")
LIMITS_SRC = os.path.join(BASE, "t-channel_plotting_scripts", "limits")
TAGGER_COMPARISON_DIRS = {
    "data": os.path.join(BASE, "WNAE_PNET_comparison_Run2_Data"),
    "signal": os.path.join(BASE, "WNAE_PNET_comparison_Run2_Signal_mDark20_yukawa1"),
}


def copy_postfits(out_dir_pnet: str, out_dir_wnae: str) -> None:
    import shutil

    os.makedirs(out_dir_pnet, exist_ok=True)
    os.makedirs(out_dir_wnae, exist_ok=True)
    for region in ("A", "B", "C", "D"):
        for tagger, src_sub, dst in (
            ("pnet", "plots_pnet", out_dir_pnet),
            ("wnae", "plots_wnae", out_dir_wnae),
        ):
            src = os.path.join(POSTFIT_SRC, src_sub, "postfit_bonly_stitched", f"plot_postfit_bonly_stitched_{region}_Run2.pdf")
            if not os.path.isfile(src):
                print(f"[WARN] missing postfit source: {src}")
                continue
            shutil.copy2(src, os.path.join(dst, os.path.basename(src)))
    print(f"[OK] copied postfit plots (regions A-D, both taggers) to {out_dir_pnet} / {out_dir_wnae}")


def generate_signal_contamination(out_dir: str) -> None:
    subprocess.run(
        [sys.executable, os.path.join(HERE, "make_signal_contamination.py"), "--out-dir", out_dir],
        check=True,
        cwd=BASE,
    )


def generate_limits(out_dir: str) -> None:
    import subprocess

    os.makedirs(out_dir, exist_ok=True)
    for tagger, csv_name in (
        ("pnet", "unblinded_gapVeto_perCateYearNonclosure_PNET-2.csv"),
        ("wnae", "unblinded_gapVeto_perCateYearNonclosure_WNAE-2.csv"),
    ):
        for suffix, observed, color_mode in (
            ("", "True", "observed_xsec"),
            ("_expected", "False", "expected_xsec"),
            ("_theory_over_expected", "False", "theory_over_expected"),
        ):
            out_pdf = os.path.join(out_dir, f"limits2d_{tagger}_mMed_rinv{suffix}.pdf")
            cmd = [
                sys.executable, os.path.join(LIMITS_SRC, "limitPlot2D-griddata_interpolation.py"),
                "--poi-x", "mMed", "--poi-y", "rinv",
                "--fixed", "mDark=20", "yukawa=1",
                "--output-csv", csv_name,
                "--output-pdf", out_pdf,
                "--observed", observed,
                "--color-mode", color_mode,
                "--tagger", tagger,
                "--lumi", r"138 fb$^{-1}$ (13 TeV)",
            ]
            subprocess.run(cmd, cwd=LIMITS_SRC, check=True)
    print(f"[OK] wrote limit plots to {out_dir}")


def generate_auc_ttbar(out_dir: str) -> None:
    sys.path.insert(0, BASE)
    import makeAUCtables_pnet as auc_mod

    os.makedirs(out_dir, exist_ok=True)
    auc_mod.plot_auc_grid(
        mMed_vals=auc_mod.TTBAR_MMED_PNET, rinv_vals=auc_mod.TTBAR_RINV_PNET, raw=auc_mod.TTBAR_PNET_RAW,
        z_label="ParticleNet AUC vs. ttbar", out_prefix=os.path.join(out_dir, "pnet_auc_grid_ttbar"),
        training_samples=auc_mod.PNET_TRAINING_SAMPLES,
    )
    auc_mod.plot_auc_grid(
        mMed_vals=auc_mod.MMED_DNN, rinv_vals=auc_mod.RINV_DNN, raw=auc_mod.DNN_RAW_TTBAR,
        z_label="Event DNN AUC vs. ttbar", out_prefix=os.path.join(out_dir, "dnn_auc_grid_ttbar"),
        training_samples=auc_mod.DNN_TRAINING_SAMPLES,
    )
    print(f"[OK] wrote ttbar AUC grids to {out_dir}")


def generate_signal_efficiency_maps() -> None:
    import subprocess

    scripts = [
        [sys.executable, os.path.join(HERE, "make_signal_efficiency_maps.py")],
        [sys.executable, os.path.join(HERE, "make_signal_efficiency_maps_extra_axes.py")],
        [sys.executable, os.path.join(HERE, "make_signal_efficiency_maps_wnae.py"), "--compute"],
    ]
    for cmd in scripts:
        subprocess.run(cmd, cwd=HERE, check=True)
    print("[OK] generated PNET and WNAE signal-efficiency maps")


def copy_tagger_comparisons(out_dir: str) -> None:
    import shutil

    os.makedirs(out_dir, exist_ok=True)
    signal_lund_dir = os.path.join(TAGGER_COMPARISON_DIRS["signal"], "wlundcorrection")
    # (source dir, source relpath, destination filename) -- the tex references
    # renamed copies (e.g. data vs. signal event_nsvj_matrix plots, and the
    # quantile-subdirectory variant of pt200to300), not the raw source names.
    files = [
        (TAGGER_COMPARISON_DIRS["data"], "WNAE_PNET_comparison_event_nsvj_matrix.pdf", "WNAE_PNET_comparison_event_nsvj_matrix_data.pdf"),
        (signal_lund_dir, "WNAE_PNET_comparison_event_nsvj_matrix.pdf", "WNAE_PNET_comparison_event_nsvj_matrix_baseline_signal.pdf"),
        (signal_lund_dir, "WNAE_PNET_comparison_pt200to300.pdf", "WNAE_PNET_comparison_pt200to300.pdf"),
        (signal_lund_dir, os.path.join("quantile", "WNAE_PNET_comparison_pt200to300.pdf"), "quantile_WNAE_PNET_comparison_pt200to300.pdf"),
    ]
    for src_dir, relpath, name in files:
        src = os.path.join(src_dir, relpath)
        if not os.path.isfile(src):
            print(f"[WARN] missing tagger-comparison source: {src}")
            continue
        shutil.copy2(src, os.path.join(out_dir, name))
    print(f"[OK] copied tagger-comparison plots to {out_dir}")


def generate_lund_systematics(out_dir: str) -> None:
    script = os.path.join(HERE, "make_lund_plots.py")
    lund_root = os.path.join(HERE, "lund-histograms.root")
    distortion_root = os.path.join(HERE, "distortion-LJP.root")
    for path in (script, lund_root, distortion_root):
        if not os.path.isfile(path):
            raise FileNotFoundError(path)
    subprocess.run(
        [
            sys.executable,
            script,
            "--from-root",
            "--lund-root",
            lund_root,
            "--distortion-root",
            distortion_root,
            "--output-dir",
            out_dir,
        ],
        check=True,
        cwd=HERE,
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--skip-merge", action="store_true", help="reuse the existing merged cache instead of rebuilding it")
    ap.add_argument("--merge-only", action="store_true", help="only (re)build the merged cache; do not draw plots")
    ap.add_argument(
        "--skip-generated-table",
        action="store_true",
        help="skip the region-A-vs-generated table (it touches EOS to read CutFlow/Initial per signal file)",
    )
    args = ap.parse_args()

    if not args.skip_merge:
        merge_caches()
    elif not os.path.isfile(CACHE):
        raise SystemExit(f"Merged cache not found: {CACHE}\nRun without --skip-merge first.")

    if args.merge_only:
        return

    cache = ROOT.TFile.Open(CACHE, "READ")
    if not cache or cache.IsZombie():
        raise SystemExit(f"Could not open merged cache: {CACHE}")

    try:
        draw_wnae_inputs(cache, os.path.join(ASSETS, "wnae_inputs"))
        draw_pnet_inputs(cache, os.path.join(ASSETS, "pnet_inputs"))
        draw_dnn_inputs(cache, os.path.join(ASSETS, "dnn_inputs"))
        draw_nminus1(cache, os.path.join(ASSETS, "nminus1"))
    finally:
        cache.Close()

    draw_background_composition(os.path.join(ASSETS, "background_composition"))
    write_signal_acceptance_table(os.path.join(ASSETS, "signal_acceptance"))
    if not args.skip_generated_table:
        write_signal_acceptance_vs_generated_table(os.path.join(ASSETS, "signal_acceptance"))

    copy_postfits(os.path.join(ASSETS, "postfit_pnet"), os.path.join(ASSETS, "postfit_wnae"))
    generate_signal_contamination(os.path.join(ASSETS, "signal_contamination"))
    generate_signal_efficiency_maps()
    generate_limits(os.path.join(ASSETS, "limits"))
    generate_auc_ttbar(os.path.join(ASSETS, "auc"))
    copy_tagger_comparisons(os.path.join(ASSETS, "tagger_comparisons"))
    generate_lund_systematics(os.path.join(ASSETS, "lund_systematics"))

    print(f"\n[DONE] wrote plots under {ASSETS}")


if __name__ == "__main__":
    main()
