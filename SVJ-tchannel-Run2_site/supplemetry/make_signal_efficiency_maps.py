#!/usr/bin/env python3
"""
Signal efficiency (N_A / N_generated) heat maps over the mMed x rinv grid,
one per nSVJ category plus inclusive -- same grid and rendering style as the
existing AUC scans (makeAUCtables_pnet.py's PNET_GRID_MMED/PNET_GRID_RINV and
plot_auc_grid), axes laid out like
t-channel_plotting_scripts/limits/plots/limits2d_pnet_mMed_rinv.pdf (mPhi vs
rinv).

Two-step, mirroring how the AUC scripts already separate scanning from
plotting (load_auc_grid_from_csv / plot_auc_grid_from_csv), so a replot never
needs to touch EOS again:

  --compute   Slow. For every (mMed, rinv) grid point: reads N_A in each
              nSVJ bin from the ABCD combine file (same source as the
              region-A acceptance tables), and N_generated from
              CutFlow/Initial in the EOS skims -- the same normalization
              NMinusOne_maker_RA2.py uses for its Lund-plane correction, and
              the same mechanism make_supplementary_plots.py already uses
              for the 10-point table, just batched over the full grid.
              Writes signal_efficiency_grid.csv.
  (default)   Fast. Reads that CSV and draws the 5 heat maps.

Usage:
  source condor/initCondor.sh
  python3 make_signal_efficiency_maps.py --compute
  python3 make_signal_efficiency_maps.py
"""
import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(os.path.dirname(HERE))  # .../t-channel_Analysis
sys.path.insert(0, BASE)

CSV_PATH = os.path.join(HERE, "signal_efficiency_grid.csv")
PLOT_DIR = os.path.join(HERE, "assets_v2", "signal_efficiency_maps")

COMBINE_FILE = (
    "/uscms/home/nparmar/nobackup/SVJ/stat_inference_unblind_gapJetveto/"
    "stat_histograms/MET_mMed-fullScan_run2_pNet_V5_DNN85_WP90_MET250_BD_lWN_"
    "allsyst_min2p_distortion_alternative_trigger_unblind_data_gapJetveto_0to3PSVJ.root"
)

SVJ_ORDER = ["0SVJ", "1SVJ", "2SVJ", "3PSVJ"]
SVJ_TITLES = {
    "0SVJ": r"$n_{\mathrm{SVJ}}^{\mathrm{PN}} = 0$ signal efficiency [%]",
    "1SVJ": r"$n_{\mathrm{SVJ}}^{\mathrm{PN}} = 1$ signal efficiency [%]",
    "2SVJ": r"$n_{\mathrm{SVJ}}^{\mathrm{PN}} = 2$ signal efficiency [%]",
    "3PSVJ": r"$n_{\mathrm{SVJ}}^{\mathrm{PN}} \geq 3$ signal efficiency [%]",
    "Inclusive": "Inclusive signal efficiency [%]",
}

# Same grid as makeAUCtables_pnet.py's PNET_GRID_MMED / PNET_GRID_RINV.
GRID_MMED = [500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000, 3500, 4000]
GRID_RINV = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]


def rinv_token(rinv: float) -> str:
    return "0" if rinv == 0 else str(rinv).replace(".", "p")


def combine_name(mmed: int, rinv: float) -> str:
    return f"mMed{mmed}_mDark20_rinv{rinv_token(rinv)}_yukawa1"


def skim_name(mmed: int, rinv: float) -> str:
    return f"t-channel_mMed-{mmed}_mDark-20_rinv-{rinv_token(rinv)}_alpha-peak_yukawa-1"


def get_all_signal_directories(f2, signal_base: str, year: str, requested):
    """Same batched per-year xrdfs_ls as f2.get_signal_directories, but without
    f2.DISABLED_SIGNAL_TOKENS -- that filter only trims Figure2's default
    5-signal overlay and should not drop otherwise-valid grid points from a
    full mMed x rinv efficiency scan."""
    import os

    signal_base = f2.normalize_eos_path(signal_base)
    nominal = f"{signal_base}/{year}/t_channel_pre_selection/nominal"
    entries = f2.xrdfs_ls(nominal, recursive=False)
    available = {
        os.path.basename(entry.rstrip("/")): entry
        for entry in entries
        if os.path.basename(entry.rstrip("/")).startswith("t-channel_")
    }
    wanted = set(requested)
    return [(name, path) for name, path in available.items() if name in wanted]


def compute() -> None:
    import uproot
    import Figure2_makerusingskims as f2
    import combinehistplotter as pie_mod

    grid = [(m, r) for m in GRID_MMED for r in GRID_RINV]

    # --- N_A per nSVJ bin, from the ABCD combine file ---
    region_yields = {}  # (mmed, rinv) -> {svj: yield, "Inclusive": yield}
    with uproot.open(COMBINE_FILE) as f:
        years = pie_mod.detect_years(f)
        for mmed, rinv in grid:
            cname = combine_name(mmed, rinv)
            per_svj = {}
            total = 0.0
            for svj in SVJ_ORDER:
                y = 0.0
                for year in years:
                    y += pie_mod.read_region_signal(f, svj, year, "A", cname)
                per_svj[svj] = y
                total += y
            per_svj["Inclusive"] = total
            region_yields[(mmed, rinv)] = per_svj
    print(f"[OK] read region-A yields for {len(grid)} grid points from the combine file")

    # --- N_generated, from CutFlow/Initial in the EOS skims. One directory
    # listing per year for the whole grid (not one per point). ---
    sample_dirs = {}  # (mmed, rinv, year) -> sample_dir or None
    for year in ("2016", "2017", "2018"):
        requested = [skim_name(m, r) for m, r in grid]
        entries = get_all_signal_directories(f2, f2.DEFAULT_SIGNAL_BASE, year, requested)
        found = dict(entries)
        for m, r in grid:
            sample_dirs[(m, r, year)] = found.get(skim_name(m, r))
    n_missing = sum(1 for v in sample_dirs.values() if v is None)
    print(f"[OK] resolved signal directories for all three years ({n_missing}/{len(sample_dirs)} sample-years missing)")

    generated = {}
    for i, (mmed, rinv) in enumerate(grid, start=1):
        total = 0.0
        any_found = False
        for year in ("2016", "2017", "2018"):
            sample_dir = sample_dirs[(mmed, rinv, year)]
            if sample_dir is None:
                continue
            any_found = True
            try:
                files = [e for e in f2.xrdfs_ls(sample_dir, recursive=True) if e.endswith(".root")]
            except RuntimeError as exc:
                print(f"[WARN] mMed={mmed} rinv={rinv} {year}: could not list files ({exc})")
                continue
            year_initial = 0.0
            for path in files:
                url = f2.EOS_HOST + "/" + f2.normalize_eos_path(path)
                try:
                    with uproot.open(url) as fin:
                        year_initial += float(fin["CutFlow"]["Initial"].array(library="np")[0])
                except Exception as exc:
                    print(f"[WARN] mMed={mmed} rinv={rinv} {year}: could not read CutFlow from {path} ({exc})")
            total += year_initial * f2.LUMI_PB[year]
        generated[(mmed, rinv)] = total if any_found else None
        if i % 20 == 0 or i == len(grid):
            print(f"  ... {i}/{len(grid)} generated yields computed")
    print("[OK] computed generated yields")

    os.makedirs(os.path.dirname(CSV_PATH), exist_ok=True)
    fieldnames = (
        ["mMed", "rinv", "N_generated"]
        + [f"NA_{svj}" for svj in SVJ_ORDER]
        + ["NA_Inclusive"]
        + [f"eff_{svj}" for svj in SVJ_ORDER]
        + ["eff_Inclusive"]
    )
    with open(CSV_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for mmed, rinv in grid:
            gen = generated[(mmed, rinv)]
            yields = region_yields[(mmed, rinv)]
            row = {"mMed": mmed, "rinv": rinv, "N_generated": gen if gen is not None else ""}
            for svj in SVJ_ORDER:
                row[f"NA_{svj}"] = yields[svj]
            row["NA_Inclusive"] = yields["Inclusive"]
            for svj in SVJ_ORDER:
                row[f"eff_{svj}"] = (yields[svj] / gen) if gen else ""
            row["eff_Inclusive"] = (yields["Inclusive"] / gen) if gen else ""
            writer.writerow(row)
    print(f"[OK] wrote {CSV_PATH}")


def load_grid(value_col: str):
    rows = []
    with open(CSV_PATH, newline="") as f:
        for row in csv.DictReader(f):
            if row[value_col] == "":
                continue
            rows.append({"mMed": float(row["mMed"]), "rinv": float(row["rinv"]), "value": float(row[value_col])})
    mmed_vals = sorted({row["mMed"] for row in rows})
    rinv_vals = sorted({row["rinv"] for row in rows})
    cells = {(row["mMed"], row["rinv"]): row["value"] for row in rows}
    raw = [[cells.get((m, r), float("nan")) for r in rinv_vals] for m in mmed_vals]
    mmed_vals = [int(m) if float(m).is_integer() else m for m in mmed_vals]
    return mmed_vals, rinv_vals, raw


def plot() -> None:
    import makeAUCtables_pnet as auc_mod

    if not os.path.isfile(CSV_PATH):
        raise SystemExit(f"{CSV_PATH} not found -- run with --compute first.")

    os.makedirs(PLOT_DIR, exist_ok=True)
    grids = {}
    global_vmax = 1.0
    for svj in SVJ_ORDER + ["Inclusive"]:
        mmed_vals, rinv_vals, raw = load_grid(f"eff_{svj}")
        # plot_auc_grid's cell labels are hardcoded to "{v:.2f}", tuned for
        # AUC values around 0.5-1.0. Our efficiencies are all a few percent
        # at most, so pass percentages instead of fractions -- otherwise
        # every cell would render as an uninformative "0.01" or "0.02".
        raw_pct = [[v * 100.0 for v in row] for row in raw]
        grids[svj] = (mmed_vals, rinv_vals, raw_pct)
        global_vmax = max(global_vmax, max((v for row in raw_pct for v in row if v == v), default=0.0))
    for svj in SVJ_ORDER + ["Inclusive"]:
        mmed_vals, rinv_vals, raw_pct = grids[svj]
        auc_mod.plot_auc_grid(
            mMed_vals=mmed_vals,
            rinv_vals=rinv_vals,
            raw=raw_pct,
            z_label=SVJ_TITLES[svj],
            out_prefix=os.path.join(PLOT_DIR, f"signal_efficiency_region_a_{svj}"),
            vmin=0.0,
            vmax=global_vmax,
            cmap="viridis",
            training_samples=auc_mod.PNET_TRAINING_SAMPLES,
            signal_efficiency_template=True,
        )
    print(f"[OK] wrote heat maps under {PLOT_DIR}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--compute", action="store_true", help="(re)compute signal_efficiency_grid.csv from EOS + the ABCD combine file")
    args = ap.parse_args()

    if args.compute:
        compute()
    plot()


if __name__ == "__main__":
    main()
