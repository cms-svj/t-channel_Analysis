#!/usr/bin/env python3
"""
Build WNAE signal-efficiency-into-region-A plots, defined as

    efficiency = N_A / N_generated

where N_A is read from signal_background_ratios.csv, summed over Run 2 years,
and N_generated is taken from the already-computed PNET signal-efficiency CSV
for the same physical signal sample. The CSV region labels are swapped relative
to the ABCD convention used in the plots: A <-> B and C <-> D.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, BASE)
sys.path.insert(0, HERE)

import make_signal_efficiency_maps as pnet_eff
import make_signal_efficiency_maps_extra_axes as pnet_eff_extra

PNET_GRID_CSV = os.path.join(HERE, "signal_efficiency_grid.csv")
PNET_MDARK_CSV = os.path.join(HERE, "signal_efficiency_mdark_scan.csv")
PNET_YUKAWA_CSV = os.path.join(HERE, "signal_efficiency_yukawa_scan.csv")
WNAE_RATIO_CSV = os.path.join(HERE, "signal_background_ratios.csv")

GRID_CSV = os.path.join(HERE, "wnae_signal_efficiency_grid.csv")
MDARK_CSV = os.path.join(HERE, "wnae_signal_efficiency_mdark_scan.csv")
YUKAWA_CSV = os.path.join(HERE, "wnae_signal_efficiency_yukawa_scan.csv")
SUMMARY_CSV = os.path.join(HERE, "wnae_signal_efficiency_baseline.csv")
PLOT_DIR = os.path.join(HERE, "assets_v2", "wnae_signal_efficiency_maps")

BASELINE = {"mMed": 2000, "rinv": 0.3, "mDark": 20, "yukawa": 1.0}
SVJ_ORDER = ["0SVJ", "1SVJ", "2SVJ", "3PSVJ"]
SVJ_TITLES = {
    "0SVJ": r"$n_{\mathrm{SVJ}}^{\mathrm{WNAE}} = 0$ signal efficiency [%]",
    "1SVJ": r"$n_{\mathrm{SVJ}}^{\mathrm{WNAE}} = 1$ signal efficiency [%]",
    "2SVJ": r"$n_{\mathrm{SVJ}}^{\mathrm{WNAE}} = 2$ signal efficiency [%]",
    "3PSVJ": r"$n_{\mathrm{SVJ}}^{\mathrm{WNAE}} \geq 3$ signal efficiency [%]",
    "Inclusive": "Inclusive signal efficiency [%]",
}
REGION_LABEL_SWAP = {"A": "B", "B": "A", "C": "D", "D": "C"}


def _as_float(value: str):
    if value == "" or value is None:
        return None
    return float(value)


def _read_csv_value(path: str, match: dict[str, float], column: str):
    with open(path, newline="") as fin:
        for row in csv.DictReader(fin):
            ok = True
            for key, want in match.items():
                if abs(float(row[key]) - float(want)) > 1e-9:
                    ok = False
                    break
            if ok:
                if row[column] == "":
                    return None
                return float(row[column])
    raise RuntimeError(f"Could not find {match} in {path}")


def _generated_lookup(path: str, keys: list[str]) -> dict[tuple, float]:
    lookup = {}
    with open(path, newline="") as fin:
        for row in csv.DictReader(fin):
            generated = _as_float(row.get("N_generated", ""))
            if generated is None:
                continue
            lookup[tuple(float(row[key]) for key in keys)] = generated
    return lookup


def _scan_yields(match: dict[str, float], x_field: str, y_field: str) -> dict[tuple[float, float], dict[str, float]]:
    """Run-2 WNAE N_A yields from signal_background_ratios.csv.

    The CSV producer labels the ABCD regions opposite to the plotted convention
    for the two axes, so true region A is stored under region B.
    """
    csv_signal_region = REGION_LABEL_SWAP["A"]
    yields = defaultdict(lambda: defaultdict(float))
    with open(WNAE_RATIO_CSV, newline="") as fin:
        for row in csv.DictReader(fin):
            if row["region"] != csv_signal_region:
                continue
            if row["category"] not in SVJ_ORDER:
                continue
            keep = True
            for key, want in match.items():
                if abs(float(row[key]) - float(want)) > 1e-9:
                    keep = False
                    break
            if not keep:
                continue
            key = (float(row[x_field]), float(row[y_field]))
            # Efficiencies must use the unit-cross-section numerator. The
            # scaled yield includes the signal cross section, which can exceed
            # the generated-yield normalization by orders of magnitude for the
            # high-yukawa scan points.
            yields[key][row["category"]] += float(row["signal_yield_unit_xsec"])
    output = {}
    for key, per_svj in yields.items():
        output[key] = {svj: per_svj.get(svj, 0.0) for svj in SVJ_ORDER}
        output[key]["Inclusive"] = sum(output[key][svj] for svj in SVJ_ORDER)
    return output


def _write_scan_csv(path: str, y_field: str, x_vals, y_vals, generated_lookup, yields_lookup) -> None:
    fieldnames = (
        ["mMed", y_field, "N_generated"]
        + [f"NA_{svj}" for svj in SVJ_ORDER]
        + ["NA_Inclusive"]
        + [f"eff_{svj}" for svj in SVJ_ORDER]
        + ["eff_Inclusive"]
    )
    with open(path, "w", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=fieldnames)
        writer.writeheader()
        for mmed in x_vals:
            for y_val in y_vals:
                key = (float(mmed), float(y_val))
                row = {"mMed": mmed, y_field: y_val, "N_generated": ""}
                if key in yields_lookup and key in generated_lookup:
                    fill_row(row, yields_lookup[key], generated_lookup[key])
                writer.writerow(row)


def fill_row(row: dict, region_a: dict[str, float], generated: float) -> None:
    row["N_generated"] = generated
    for svj in SVJ_ORDER:
        row[f"NA_{svj}"] = region_a[svj]
        row[f"eff_{svj}"] = region_a[svj] / generated
    row["NA_Inclusive"] = region_a["Inclusive"]
    row["eff_Inclusive"] = region_a["Inclusive"] / generated


def _write_summary_csv(region_a: dict[str, float], generated: float) -> None:
    with open(SUMMARY_CSV, "w", newline="") as fout:
        fieldnames = ["category", "N_A", "N_generated", "efficiency", "efficiency_percent"]
        writer = csv.DictWriter(fout, fieldnames=fieldnames)
        writer.writeheader()
        for svj in SVJ_ORDER + ["Inclusive"]:
            writer.writerow(
                {
                    "category": svj,
                    "N_A": region_a[svj],
                    "N_generated": generated,
                    "efficiency": region_a[svj] / generated,
                    "efficiency_percent": 100.0 * region_a[svj] / generated,
                }
            )


def compute() -> None:
    grid_generated = _generated_lookup(PNET_GRID_CSV, ["mMed", "rinv"])
    mdark_generated = _generated_lookup(PNET_MDARK_CSV, ["mMed", "mDark"])
    yukawa_generated = _generated_lookup(PNET_YUKAWA_CSV, ["mMed", "yukawa"])

    grid_yields = _scan_yields(
        {"mDark": BASELINE["mDark"], "yukawa": BASELINE["yukawa"]},
        "mMed",
        "rinv",
    )
    mdark_yields = _scan_yields(
        {"rinv": BASELINE["rinv"], "yukawa": BASELINE["yukawa"]},
        "mMed",
        "mDark",
    )
    yukawa_yields = _scan_yields(
        {"rinv": BASELINE["rinv"], "mDark": BASELINE["mDark"]},
        "mMed",
        "yukawa",
    )

    _write_scan_csv(GRID_CSV, "rinv", pnet_eff.GRID_MMED, pnet_eff.GRID_RINV, grid_generated, grid_yields)
    _write_scan_csv(MDARK_CSV, "mDark", pnet_eff_extra.GRID_MMED, pnet_eff_extra.GRID_MDARK, mdark_generated, mdark_yields)
    _write_scan_csv(YUKAWA_CSV, "yukawa", pnet_eff_extra.GRID_MMED, pnet_eff_extra.GRID_YUKAWA, yukawa_generated, yukawa_yields)

    baseline_key = (float(BASELINE["mMed"]), float(BASELINE["rinv"]))
    if baseline_key in grid_yields and baseline_key in grid_generated:
        _write_summary_csv(grid_yields[baseline_key], grid_generated[baseline_key])
    else:
        raise RuntimeError("Could not find baseline WNAE point in the signal/background ratio CSV and generated-yield CSV")

    print(f"[OK] wrote {GRID_CSV}")
    print(f"[OK] wrote {MDARK_CSV}")
    print(f"[OK] wrote {YUKAWA_CSV}")
    print(f"[OK] wrote {SUMMARY_CSV}")


def load_grid(csv_path: str, x_field: str, y_field: str, value_col: str):
    axis_rows = []
    valued_rows = []
    with open(csv_path, newline="") as fin:
        for row in csv.DictReader(fin):
            axis_rows.append({x_field: float(row[x_field]), y_field: float(row[y_field])})
            if row[value_col] == "":
                continue
            valued_rows.append(
                {
                    x_field: float(row[x_field]),
                    y_field: float(row[y_field]),
                    "value": float(row[value_col]),
                }
            )
    x_vals = sorted({row[x_field] for row in axis_rows})
    y_vals = sorted({row[y_field] for row in axis_rows})
    cells = {(row[x_field], row[y_field]): row["value"] for row in valued_rows}
    raw = [[cells.get((x, y), float("nan")) for y in y_vals] for x in x_vals]
    x_vals = [int(x) if float(x).is_integer() else x for x in x_vals]
    return x_vals, y_vals, raw


def reference_vmax(csv_path: str, x_field: str, y_field: str) -> float:
    if not os.path.isfile(csv_path):
        return 1.0
    vmax = 1.0
    for svj in SVJ_ORDER + ["Inclusive"]:
        _, _, raw = load_grid(csv_path, x_field, y_field, f"eff_{svj}")
        raw_pct = [[v * 100.0 for v in row] for row in raw]
        vmax = max(vmax, max((v for row in raw_pct for v in row if v == v), default=0.0))
    return vmax


def plot_scan(csv_path: str, y_field: str, y_label: str, y_tick_fmt, out_tag: str, reference_csv=None) -> None:
    import makeAUCtables_pnet as auc_mod

    os.makedirs(PLOT_DIR, exist_ok=True)
    grids = {}
    global_vmax = reference_vmax(reference_csv, "mMed", y_field) if reference_csv else 1.0
    for svj in SVJ_ORDER + ["Inclusive"]:
        mmed_vals, y_vals, raw = load_grid(csv_path, "mMed", y_field, f"eff_{svj}")
        raw_pct = [[v * 100.0 for v in row] for row in raw]
        finite = [v for row in raw_pct for v in row if v == v]
        global_vmax = max(global_vmax, max(finite) if finite else 0.0)
        grids[svj] = (mmed_vals, y_vals, raw_pct)
    for svj in SVJ_ORDER + ["Inclusive"]:
        mmed_vals, y_vals, raw_pct = grids[svj]
        auc_mod.plot_auc_grid(
            mMed_vals=mmed_vals,
            rinv_vals=y_vals,
            raw=raw_pct,
            z_label=SVJ_TITLES[svj],
            out_prefix=os.path.join(PLOT_DIR, f"wnae_signal_efficiency_region_a_{out_tag}_{svj}"),
            vmin=0.0,
            vmax=global_vmax,
            cmap="viridis",
            training_samples=None,
            y_label=y_label,
            y_tick_fmt=y_tick_fmt,
            signal_efficiency_template=True,
        )


def plot_yukawa_summary() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import mplhep as hep
    from matplotlib.ticker import AutoMinorLocator, MultipleLocator

    selected_yukawas = (0.5, 1.0, 2.0)
    line_styles = {
        0.5: (0, (1.2, 1.6)),
        1.0: "-",
        2.0: (0, (6.0, 2.2)),
    }
    points = {yukawa: [] for yukawa in selected_yukawas}
    with open(YUKAWA_CSV, newline="") as fin:
        for row in csv.DictReader(fin):
            yukawa = float(row["yukawa"])
            if yukawa not in points or not row["eff_Inclusive"]:
                continue
            points[yukawa].append(
                (float(row["mMed"]), 100.0 * float(row["eff_Inclusive"]))
            )

    missing = [yukawa for yukawa, values in points.items() if not values]
    if missing:
        raise RuntimeError(f"Missing WNAE inclusive efficiencies for Yukawa values: {missing}")

    hep.style.use("CMS")
    fig, ax = plt.subplots(figsize=(11, 9))
    for yukawa in selected_yukawas:
        values = sorted(points[yukawa])
        masses = [value[0] for value in values]
        efficiencies = [value[1] for value in values]
        ax.plot(
            masses,
            efficiencies,
            color="#118AB2",
            linestyle=line_styles[yukawa],
            linewidth=3.2,
            marker="o",
            markersize=6.5,
            markerfacecolor="white",
            markeredgewidth=1.8,
            label=rf"$\lambda = {yukawa:g}$",
        )

    ax.set_xlim(500, 4000)
    ax.set_ylim(bottom=0.0)
    ax.set_xticks(list(range(500, 4001, 500)))
    ax.xaxis.set_minor_locator(MultipleLocator(100))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel(r"$m_{\Phi}$ [GeV]")
    ax.set_ylabel("WNAE inclusive signal efficiency [%]")
    ax.tick_params(which="major", length=10)
    ax.tick_params(which="minor", length=5)

    hep.cms.label(data=False, ax=ax, loc=0, com=13)
    ax.text(
        0.52,
        1.012,
        r"$m_{\mathrm{dark}} = 20$ GeV, $r_{\mathrm{inv}} = 0.3$",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=18,
    )
    ax.legend(loc="upper right", frameon=False, fontsize=21, handlelength=3.2)

    fig.tight_layout()
    os.makedirs(PLOT_DIR, exist_ok=True)
    out_prefix = os.path.join(
        PLOT_DIR,
        "wnae_signal_efficiency_region_a_yukawa_inclusive_summary",
    )
    fig.savefig(out_prefix + ".pdf", bbox_inches="tight", pad_inches=0.12)
    fig.savefig(out_prefix + ".png", bbox_inches="tight", pad_inches=0.12, dpi=180)
    plt.close(fig)
    print(f"[OK] wrote {out_prefix}.pdf/.png")


def plot() -> None:
    if not os.path.isfile(GRID_CSV):
        raise SystemExit(f"{GRID_CSV} not found -- run with --compute first.")
    plot_scan(GRID_CSV, "rinv", r"$r_{\mathrm{inv}}$", lambda v: f"{v:.1f}", "mmed_rinv", PNET_GRID_CSV)
    plot_scan(MDARK_CSV, "mDark", r"$m_{\mathrm{dark}}\ [\mathrm{GeV}]$", lambda v: f"{v:g}", "mdark", PNET_MDARK_CSV)
    plot_scan(YUKAWA_CSV, "yukawa", r"$\lambda$", lambda v: f"{v:g}", "yukawa", PNET_YUKAWA_CSV)
    plot_yukawa_summary()
    print(f"[OK] wrote WNAE signal-efficiency maps under {PLOT_DIR}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--compute", action="store_true", help="read WNAE fitDiagnostics and write WNAE efficiency CSVs")
    args = parser.parse_args()
    if args.compute:
        compute()
    plot()


if __name__ == "__main__":
    main()
