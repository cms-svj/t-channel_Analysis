#!/usr/bin/env python3
"""
Two extra signal-efficiency-into-region-A grids (N_A / N_generated), in
addition to the existing mMed-vs-rinv grid in make_signal_efficiency_maps.py:

  1. mDark vs mMed   (rinv=0.3, yukawa=1 fixed)
  2. yukawa vs mMed  (rinv=0.3, mDark=20 fixed)

Same method as make_signal_efficiency_maps.py: N_A per nSVJ category (and
inclusive) from the ABCD combine file, N_generated from CutFlow/Initial in
the EOS skims (un-reweighted by the Lund correction, matching the existing
mMed-vs-rinv grid's convention -- CutFlow/Initial is the generator-level
total and stays valid regardless of downstream Lund reweighting). Signal
directory discovery bypasses Figure2's DISABLED_SIGNAL_TOKENS filter, same
fix as get_all_signal_directories() in make_signal_efficiency_maps.py.

Two-step:
  --compute   Slow. Writes mdark_scan_efficiency.csv and
              yukawa_scan_efficiency.csv.
  (default)   Fast. Reads those CSVs and draws the heat maps: one combined
              2x2 (0SVJ/1SVJ/2SVJ/3PSVJ) figure per scan, plus a separate
              Inclusive figure per scan.

Usage:
  source condor/initCondor.sh
  python3 make_signal_efficiency_maps_extra_axes.py --compute
  python3 make_signal_efficiency_maps_extra_axes.py
"""
import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(os.path.dirname(HERE))  # .../t-channel_Analysis
sys.path.insert(0, BASE)
sys.path.insert(0, HERE)

import make_signal_efficiency_maps as sigeff  # reuses COMBINE_FILE, SVJ_ORDER, SVJ_TITLES, get_all_signal_directories

OUT_DIR = HERE
MDARK_CSV = os.path.join(OUT_DIR, "signal_efficiency_mdark_scan.csv")
YUKAWA_CSV = os.path.join(OUT_DIR, "signal_efficiency_yukawa_scan.csv")
PLOT_DIR = os.path.join(HERE, "assets_v2", "signal_efficiency_maps")

FIXED_RINV = 0.3
FIXED_YUKAWA_FOR_MDARK_SCAN = 1.0
FIXED_MDARK_FOR_YUKAWA_SCAN = 20.0
ALPHA = "peak"

GRID_MMED = [500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000, 3500, 4000]
GRID_MDARK = [1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
GRID_YUKAWA = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5]


def yukawa_token(y: float) -> str:
    return str(int(y)) if float(y).is_integer() else str(y).replace(".", "p")


def combine_name_mdark_scan(mmed: int, mdark: int) -> str:
    return f"mMed{mmed}_mDark{mdark}_rinv{sigeff.rinv_token(FIXED_RINV)}_yukawa1"


def combine_name_yukawa_scan(mmed: int, yukawa: float) -> str:
    return f"mMed{mmed}_mDark20_rinv{sigeff.rinv_token(FIXED_RINV)}_yukawa{yukawa_token(yukawa)}"


def skim_name_mdark_scan(mmed: int, mdark: int) -> str:
    return f"t-channel_mMed-{mmed}_mDark-{mdark}_rinv-{sigeff.rinv_token(FIXED_RINV)}_alpha-{ALPHA}_yukawa-1"


def skim_name_yukawa_scan(mmed: int, yukawa: float) -> str:
    return f"t-channel_mMed-{mmed}_mDark-20_rinv-{sigeff.rinv_token(FIXED_RINV)}_alpha-{ALPHA}_yukawa-{yukawa_token(yukawa)}"


def compute() -> None:
    import uproot
    import Figure2_makerusingskims as f2
    import combinehistplotter as pie_mod

    for scan, axis_field, axis_vals, combine_name_fn, skim_name_fn, csv_path in (
        ("mdark", "mDark", GRID_MDARK, combine_name_mdark_scan, skim_name_mdark_scan, MDARK_CSV),
        ("yukawa", "yukawa", GRID_YUKAWA, combine_name_yukawa_scan, skim_name_yukawa_scan, YUKAWA_CSV),
    ):
        print(f"\n=== {scan} scan ===")
        grid = [(m, a) for m in GRID_MMED for a in axis_vals]

        # --- N_A per nSVJ bin, from the ABCD combine file ---
        # read_region_signal() returns 0.0 both when a signal genuinely has no
        # events in a region AND when the signal process is entirely absent
        # from the combine file (e.g. yukawa=0.1 isn't simulated there at
        # all). Those two cases must not be conflated -- an absent signal
        # should leave the grid cell blank, not silently render as 0.00.
        import re

        with uproot.open(sigeff.COMBINE_FILE) as f:
            available_signal_tokens = set()
            for key in f.keys():
                m = re.search(r"mMed\d+_mDark\d+_rinv[\dp]+_yukawa[\dp]+", key)
                if m:
                    available_signal_tokens.add(m.group(0))

            region_yields = {}
            n_absent = 0
            years = pie_mod.detect_years(f)
            for mmed, axis_val in grid:
                cname = combine_name_fn(mmed, axis_val)
                if cname not in available_signal_tokens:
                    region_yields[(mmed, axis_val)] = None
                    n_absent += 1
                    continue
                per_svj = {}
                total = 0.0
                for svj in sigeff.SVJ_ORDER:
                    y = 0.0
                    for year in years:
                        y += pie_mod.read_region_signal(f, svj, year, "A", cname)
                    per_svj[svj] = y
                    total += y
                per_svj["Inclusive"] = total
                region_yields[(mmed, axis_val)] = per_svj
        print(f"[OK] read region-A yields for {len(grid)} grid points from the combine file ({n_absent} absent from the file entirely)")

        # --- N_generated, from CutFlow/Initial in the EOS skims ---
        sample_dirs = {}
        for year in ("2016", "2017", "2018"):
            requested = [skim_name_fn(m, a) for m, a in grid]
            entries = sigeff.get_all_signal_directories(f2, f2.DEFAULT_SIGNAL_BASE, year, requested)
            found = dict(entries)
            for m, a in grid:
                sample_dirs[(m, a, year)] = found.get(skim_name_fn(m, a))
        n_missing = sum(1 for v in sample_dirs.values() if v is None)
        print(f"[OK] resolved signal directories for all three years ({n_missing}/{len(sample_dirs)} sample-years missing)")

        generated = {}
        for i, (mmed, axis_val) in enumerate(grid, start=1):
            total = 0.0
            any_found = False
            for year in ("2016", "2017", "2018"):
                sample_dir = sample_dirs[(mmed, axis_val, year)]
                if sample_dir is None:
                    continue
                any_found = True
                try:
                    files = [e for e in f2.xrdfs_ls(sample_dir, recursive=True) if e.endswith(".root")]
                except RuntimeError as exc:
                    print(f"[WARN] mMed={mmed} {axis_field}={axis_val} {year}: could not list files ({exc})")
                    continue
                year_initial = 0.0
                for path in files:
                    url = f2.EOS_HOST + "/" + f2.normalize_eos_path(path)
                    try:
                        with uproot.open(url) as fin:
                            year_initial += float(fin["CutFlow"]["Initial"].array(library="np")[0])
                    except Exception as exc:
                        print(f"[WARN] mMed={mmed} {axis_field}={axis_val} {year}: could not read CutFlow from {path} ({exc})")
                total += year_initial * f2.LUMI_PB[year]
            generated[(mmed, axis_val)] = total if any_found else None
            if i % 20 == 0 or i == len(grid):
                print(f"  ... {i}/{len(grid)} generated yields computed")
        print("[OK] computed generated yields")

        fieldnames = (
            ["mMed", axis_field, "N_generated"]
            + [f"NA_{svj}" for svj in sigeff.SVJ_ORDER]
            + ["NA_Inclusive"]
            + [f"eff_{svj}" for svj in sigeff.SVJ_ORDER]
            + ["eff_Inclusive"]
        )
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for mmed, axis_val in grid:
                gen = generated[(mmed, axis_val)]
                yields = region_yields[(mmed, axis_val)]
                row = {"mMed": mmed, axis_field: axis_val, "N_generated": gen if gen is not None else ""}
                if yields is None:
                    for svj in sigeff.SVJ_ORDER:
                        row[f"NA_{svj}"] = ""
                        row[f"eff_{svj}"] = ""
                    row["NA_Inclusive"] = ""
                    row["eff_Inclusive"] = ""
                else:
                    for svj in sigeff.SVJ_ORDER:
                        row[f"NA_{svj}"] = yields[svj]
                    row["NA_Inclusive"] = yields["Inclusive"]
                    for svj in sigeff.SVJ_ORDER:
                        row[f"eff_{svj}"] = (yields[svj] / gen) if gen else ""
                    row["eff_Inclusive"] = (yields["Inclusive"] / gen) if gen else ""
                writer.writerow(row)
        print(f"[OK] wrote {csv_path}")


def load_grid(csv_path, axis_field, value_col):
    rows = []
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            if row[value_col] == "":
                continue
            rows.append({"mMed": float(row["mMed"]), axis_field: float(row[axis_field]), "value": float(row[value_col])})
    mmed_vals = sorted({row["mMed"] for row in rows})
    axis_vals = sorted({row[axis_field] for row in rows})
    cells = {(row["mMed"], row[axis_field]): row["value"] for row in rows}
    raw = [[cells.get((m, a), float("nan")) for a in axis_vals] for m in mmed_vals]
    mmed_vals = [int(m) if float(m).is_integer() else m for m in mmed_vals]
    return mmed_vals, axis_vals, raw


def plot_scan(csv_path, axis_field, y_label, y_tick_fmt, out_tag) -> None:
    import makeAUCtables_pnet as auc_mod

    if not os.path.isfile(csv_path):
        raise SystemExit(f"{csv_path} not found -- run with --compute first.")

    os.makedirs(PLOT_DIR, exist_ok=True)

    grids = {}
    global_vmax = 1.0
    for svj in sigeff.SVJ_ORDER + ["Inclusive"]:
        mmed_vals, axis_vals, raw = load_grid(csv_path, axis_field, f"eff_{svj}")
        raw_pct = [[v * 100.0 for v in row] for row in raw]
        grids[svj] = (mmed_vals, axis_vals, raw_pct)
        global_vmax = max(global_vmax, max((v for row in raw_pct for v in row if v == v), default=0.0))
    for svj in sigeff.SVJ_ORDER + ["Inclusive"]:
        mmed_vals, axis_vals, raw_pct = grids[svj]
        z_label = sigeff.SVJ_TITLES[svj]
        auc_mod.plot_auc_grid(
            mMed_vals=mmed_vals, rinv_vals=axis_vals, raw=raw_pct,
            z_label=z_label,
            out_prefix=os.path.join(PLOT_DIR, f"signal_efficiency_region_a_{out_tag}_{svj}"),
            vmin=0.0, vmax=global_vmax, cmap="viridis",
            y_label=y_label, y_tick_fmt=y_tick_fmt,
            signal_efficiency_template=True,
        )
    print(f"[OK] wrote {out_tag} heat maps under {PLOT_DIR}")


def plot() -> None:
    plot_scan(MDARK_CSV, "mDark", r"$m_{\mathrm{dark}}\ [\mathrm{GeV}]$", lambda v: f"{v:g}", "mdark")
    plot_scan(YUKAWA_CSV, "yukawa", r"$\lambda$", lambda v: f"{v:g}", "yukawa")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--compute", action="store_true", help="(re)compute the CSVs from EOS + the ABCD combine file")
    args = ap.parse_args()
    if args.compute:
        compute()
    plot()


if __name__ == "__main__":
    main()
