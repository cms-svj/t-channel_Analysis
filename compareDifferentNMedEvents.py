#!/usr/bin/env python3
"""Plot direct, associated, and pair production fractions versus mediator mass.

The input files are the unselected 2018 TreeMaker t-channel samples generated
at lambda=1. The mediator multiplicity is reconstructed exactly as in
utils/utility.py by counting generator particles with |PDG ID| in
4900001--4900006. Fractions at other Yukawa couplings are obtained from the
leading-order topology dependence: direct scales as lambda^4, associated as
lambda^2, and QCD pair production as lambda^0.

The first pass writes compact CSV and ROOT caches. Later style-only reruns can
use --plot-only and never touch EOS.

Examples:
  source condor/initCondor.sh
  python3 compareDifferentNMedEvents.py
  python3 compareDifferentNMedEvents.py --plot-only
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List

import awkward as ak
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import uproot


HERE = Path(__file__).resolve().parent
SUPPLEMENT_DIR = HERE / "SVJ-tchannel-Run2_site" / "supplemetry"
DEFAULT_OUTPUT_DIR = SUPPLEMENT_DIR / "assets_v2" / "yukawa_nmed"
DEFAULT_CSV_CACHE = SUPPLEMENT_DIR / "yukawa_nmed_fractions.csv"
DEFAULT_ROOT_CACHE = SUPPLEMENT_DIR / "yukawa_nmed_fractions.root"

sys.path.insert(0, str(HERE))
import Figure2_makerusingskims as figure2


DEFAULT_TREE_MAKER_BASE = "/store/user/lpcdarkqcd/tchannel_UL/2018/Full/PrivateSamples"
DEFAULT_MASSES = [600, 800, 1000, 1500, 2000, 3000, 4000]
DEFAULT_YUKAWAS = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5]
MEDIATOR_PDG_IDS = np.asarray([4900001, 4900002, 4900003, 4900004, 4900005, 4900006])
TOPOLOGIES = (
    ("direct", "Direct", 0, 4, "#e42536"),
    ("associated", "Associated", 1, 2, "#5790fc"),
    ("pair", "Pair", 2, 0, "#009e73"),
)


def parse_csv_numbers(value: str, cast=float) -> List:
    return [cast(item.strip()) for item in value.split(",") if item.strip()]


def yukawa_token(value: float) -> str:
    return str(int(value)) if float(value).is_integer() else str(value).replace(".", "p")


def sample_name(mass: int, yukawa: float = 1.0) -> str:
    return (
        f"t-channel_mMed-{mass}_mDark-20_rinv-0p3_"
        f"alpha-peak_yukawa-{yukawa_token(yukawa)}"
    )


def tree_maker_directory_name(mass: int) -> str:
    return (
        f"SVJ_UL2018_{sample_name(mass)}_"
        "13TeV-madgraphMLM-pythia8_n-1000"
    )


def discover_samples(tree_maker_base: str, masses: Iterable[int]) -> Dict[int, str]:
    base = figure2.normalize_eos_path(tree_maker_base)
    available = {
        os.path.basename(entry.rstrip("/")): entry
        for entry in figure2.xrdfs_ls(base, recursive=False)
    }
    return {
        mass: available[tree_maker_directory_name(mass)]
        for mass in masses
        if tree_maker_directory_name(mass) in available
    }


def tree_path(root_file) -> str:
    if "TreeMaker2/PreSelection" in root_file:
        return "TreeMaker2/PreSelection"
    if "TreeMaker2;1/PreSelection;1" in root_file:
        return "TreeMaker2;1/PreSelection;1"
    raise KeyError("Could not find TreeMaker2/PreSelection")


def tree_maker_files(sample_dir: str, max_files: int) -> List[str]:
    pattern = re.compile(r".*_RA2AnalysisTree\.root$")
    files = sorted(
        path
        for path in figure2.xrdfs_ls(sample_dir, recursive=False)
        if pattern.search(os.path.basename(path))
    )
    return files[:max_files] if max_files > 0 else files


def accumulate_file(file_url: str) -> dict:
    result = {
        "n_events": 0,
        "n_overflow": 0,
        "direct_sumw": 0.0,
        "direct_sumw2": 0.0,
        "associated_sumw": 0.0,
        "associated_sumw2": 0.0,
        "pair_sumw": 0.0,
        "pair_sumw2": 0.0,
    }
    with uproot.open(file_url) as root_file:
        tree = root_file[tree_path(root_file)]
        if "GenParticles_PdgId" not in tree:
            raise KeyError(f"{file_url} lacks GenParticles_PdgId")
        abs_pdg_id = abs(tree["GenParticles_PdgId"].array(library="ak"))
        is_mediator = ak.zeros_like(abs_pdg_id, dtype=bool)
        for pdg_id in MEDIATOR_PDG_IDS:
            is_mediator = is_mediator | (abs_pdg_id == int(pdg_id))
        n_mediator = ak.to_numpy(ak.sum(is_mediator, axis=1))

        result["n_events"] = int(len(n_mediator))
        result["n_overflow"] = int(np.count_nonzero(n_mediator > 2))
        for key, _, multiplicity, _, _ in TOPOLOGIES:
            count = float(np.count_nonzero(n_mediator == multiplicity))
            result[f"{key}_sumw"] = count
            result[f"{key}_sumw2"] = count
    return result


def merge_result(target: dict, source: dict) -> None:
    for key, value in source.items():
        target[key] = target.get(key, 0) + value


def fraction_uncertainty(category_sumw: float, category_sumw2: float, total_sumw: float, total_sumw2: float) -> float:
    if total_sumw == 0.0:
        return float("nan")
    fraction = category_sumw / total_sumw
    other_sumw2 = max(0.0, total_sumw2 - category_sumw2)
    variance = ((1.0 - fraction) ** 2 * category_sumw2 + fraction**2 * other_sumw2) / total_sumw**2
    return math.sqrt(max(0.0, variance))


def finalize_row(mass: int, yukawa: float, accum: dict) -> dict:
    scaled = {}
    for key, _, _, coupling_power, _ in TOPOLOGIES:
        factor = yukawa**coupling_power
        scaled[f"{key}_sumw"] = accum[f"{key}_sumw"] * factor
        scaled[f"{key}_sumw2"] = accum[f"{key}_sumw2"] * factor * factor

    total_sumw = sum(scaled[f"{key}_sumw"] for key, _, _, _, _ in TOPOLOGIES)
    total_sumw2 = sum(scaled[f"{key}_sumw2"] for key, _, _, _, _ in TOPOLOGIES)
    row = {
        "mPhi": mass,
        "yukawa": yukawa,
        "n_events": accum["n_events"],
        "n_overflow": accum["n_overflow"],
        "total_sumw": total_sumw,
        "total_sumw2": total_sumw2,
    }
    for key, _, _, _, _ in TOPOLOGIES:
        sumw = scaled[f"{key}_sumw"]
        sumw2 = scaled[f"{key}_sumw2"]
        row[f"{key}_source_count"] = int(accum[f"{key}_sumw"])
        row[f"{key}_sumw"] = sumw
        row[f"{key}_sumw2"] = sumw2
        row[f"{key}_fraction"] = sumw / total_sumw if total_sumw else float("nan")
        row[f"{key}_uncertainty"] = fraction_uncertainty(sumw, sumw2, total_sumw, total_sumw2)
    return row


def compute(args) -> List[dict]:
    directories = discover_samples(args.tree_maker_base, args.masses)
    missing = sorted(set(args.masses) - set(directories))
    if missing:
        raise RuntimeError(f"Missing unselected TreeMaker samples for masses: {missing}")

    baseline = {}
    for mass in args.masses:
        files = tree_maker_files(directories[mass], args.max_files)
        if not files:
            raise RuntimeError(f"No TreeMaker ROOT files found under {directories[mass]}")
        print(f"[READ] lambda=1, mPhi={mass}: {len(files)} unselected TreeMaker files")
        accum = {}
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = {
                executor.submit(accumulate_file, figure2.eos_url(path)): path
                for path in files
            }
            for index, future in enumerate(as_completed(futures), start=1):
                merge_result(accum, future.result())
                if index % 100 == 0 or index == len(futures):
                    print(f"[READ] mPhi={mass}: processed {index}/{len(futures)} files")
        baseline[mass] = accum

    rows = []
    for yukawa in args.yukawas:
        for mass in args.masses:
            row = finalize_row(mass, yukawa, baseline[mass])
            rows.append(row)
            print(
                f"[OK] lambda={yukawa:g}, mPhi={mass}: "
                f"Direct={row['direct_fraction']:.4f}, "
                f"Associated={row['associated_fraction']:.4f}, "
                f"Pair={row['pair_fraction']:.4f}"
            )
    return rows


def write_caches(rows: List[dict], csv_path: Path, root_path: Path) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with csv_path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    root_path.parent.mkdir(parents=True, exist_ok=True)
    integer_fields = {
        "mPhi",
        "n_events",
        "n_overflow",
        *(f"{key}_source_count" for key, _, _, _, _ in TOPOLOGIES),
    }
    columns = {
        key: np.asarray(
            [row[key] for row in rows],
            dtype=np.int64 if key in integer_fields else np.float64,
        )
        for key in fieldnames
    }
    with uproot.recreate(root_path) as output:
        output["yukawa_nmed_fractions"] = columns
    print(f"[OK] wrote {csv_path}")
    print(f"[OK] wrote {root_path}")


def read_cache(csv_path: Path) -> List[dict]:
    integer_fields = {
        "mPhi",
        "n_events",
        "n_overflow",
        *(f"{key}_source_count" for key, _, _, _, _ in TOPOLOGIES),
    }
    rows = []
    with csv_path.open(newline="") as source:
        for raw in csv.DictReader(source):
            rows.append(
                {
                    key: int(float(value)) if key in integer_fields else float(value)
                    for key, value in raw.items()
                }
            )
    return rows


def draw_plot(rows: List[dict], yukawa: float, output_dir: Path, preliminary: bool) -> None:
    selected = sorted(
        (row for row in rows if np.isclose(row["yukawa"], yukawa)),
        key=lambda row: row["mPhi"],
    )
    if not selected:
        raise RuntimeError(f"No cached rows found for lambda={yukawa:g}")

    hep.style.use("CMS")
    fig, ax = plt.subplots(figsize=(9.2, 7.4))
    masses = np.asarray([row["mPhi"] for row in selected], dtype=float)
    for key, label, _, _, color in TOPOLOGIES:
        fractions = np.asarray([row[f"{key}_fraction"] for row in selected])
        uncertainties = np.asarray([row[f"{key}_uncertainty"] for row in selected])
        ax.errorbar(
            masses,
            fractions,
            yerr=uncertainties,
            color=color,
            marker="o",
            markersize=5.5,
            linewidth=2.8,
            capsize=2.5,
            label=label,
        )

    ax.set_xscale("log")
    ax.set_xlim(450, 4500)
    ax.set_ylim(0.0, 1.04)
    ax.set_xticks(masses)
    ax.set_xticklabels([str(int(value)) for value in masses], rotation=55, ha="right")
    ax.set_xlabel(r"$m_{\Phi}$ [GeV]")
    ax.set_ylabel("Fraction of events")
    ax.minorticks_on()
    ax.tick_params(which="major", length=9)
    ax.tick_params(which="minor", length=4)

    hep.cms.label("Preliminary" if preliminary else None, data=False, com=13, ax=ax)
    ax.text(
        0.04,
        0.91,
        rf"$m_{{\mathrm{{dark}}}}=20$ GeV, $r_{{\mathrm{{inv}}}}=0.3$, $\lambda={yukawa:g}$",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=18,
    )
    ax.legend(loc="center right", frameon=False, fontsize=19, handlelength=2.0)
    fig.tight_layout()

    output_dir.mkdir(parents=True, exist_ok=True)
    token = yukawa_token(yukawa)
    for extension in ("pdf", "png"):
        path = output_dir / f"nmed_fractions_vs_mphi_yukawa-{token}.{extension}"
        fig.savefig(path, bbox_inches="tight", pad_inches=0.12, dpi=180 if extension == "png" else None)
        print(f"[OK] wrote {path}")
    plt.close(fig)


def draw_combined_plot(rows: List[dict], output_dir: Path, preliminary: bool) -> None:
    from matplotlib.lines import Line2D
    from matplotlib.ticker import AutoMinorLocator, MultipleLocator

    selected_yukawas = (0.5, 1.0, 2.0)
    line_styles = {
        0.5: (0, (1.2, 1.6)),
        1.0: "-",
        2.0: (0, (6.0, 2.2)),
    }
    cms_colors = {
        "direct": "#5790fc",
        "associated": "#f89c20",
        "pair": "#e42536",
    }
    markers = {"direct": "o", "associated": "s", "pair": "^"}
    by_yukawa = {
        yukawa: sorted(
            (row for row in rows if np.isclose(row["yukawa"], yukawa)),
            key=lambda row: row["mPhi"],
        )
        for yukawa in selected_yukawas
    }
    missing = [yukawa for yukawa, values in by_yukawa.items() if not values]
    if missing:
        raise RuntimeError(f"Missing cached mediator fractions for Yukawa values: {missing}")

    hep.style.use("CMS")
    fig, ax = plt.subplots(figsize=(11.5, 9.0))
    ax.axhline(1.0, color="#4d4d4d", linestyle="--", linewidth=1.4, alpha=0.75, zorder=0)
    for yukawa in selected_yukawas:
        selected = by_yukawa[yukawa]
        masses = np.asarray([row["mPhi"] for row in selected], dtype=float)
        for key, _, _, _, _ in TOPOLOGIES:
            fractions = np.asarray([row[f"{key}_fraction"] for row in selected])
            ax.plot(
                masses,
                fractions,
                color=cms_colors[key],
                linestyle=line_styles[yukawa],
                linewidth=2.8,
                marker=markers[key],
                markersize=5.8,
                markerfacecolor="white",
                markeredgewidth=1.5,
            )

    ax.set_xlim(500, 4000)
    ax.set_ylim(0.0, 1.24)
    ax.set_xticks(list(range(500, 4001, 500)))
    ax.xaxis.set_minor_locator(MultipleLocator(100))
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel(r"$m_{\Phi}$ [GeV]")
    ax.set_ylabel("Fraction of events")
    ax.tick_params(which="major", length=10)
    ax.tick_params(which="minor", length=5)

    hep.cms.label("Preliminary" if preliminary else None, data=False, com=13, ax=ax)
    ax.text(
        0.96,
        0.70,
        r"$m_{\mathrm{dark}} = 20$ GeV, $r_{\mathrm{inv}} = 0.3$",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=18,
    )

    topology_handles = [
        Line2D(
            [],
            [],
            color=cms_colors[key],
            marker=markers[key],
            markerfacecolor="white",
            markeredgewidth=1.5,
            linewidth=2.8,
            label=label,
        )
        for key, label, _, _, _ in TOPOLOGIES
    ]
    topology_legend = ax.legend(
        handles=topology_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.985),
        ncol=3,
        frameon=False,
        fontsize=18,
        handlelength=2.2,
        columnspacing=1.6,
    )
    ax.add_artist(topology_legend)

    yukawa_handles = [
        Line2D(
            [],
            [],
            color="black",
            linestyle=line_styles[yukawa],
            linewidth=2.8,
            label=rf"$\lambda = {yukawa:g}$",
        )
        for yukawa in selected_yukawas
    ]
    ax.legend(
        handles=yukawa_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.91),
        ncol=3,
        frameon=False,
        fontsize=18,
        handlelength=2.8,
        columnspacing=2.0,
    )

    fig.tight_layout()
    output_dir.mkdir(parents=True, exist_ok=True)
    for extension in ("pdf", "png"):
        path = output_dir / f"nmed_fractions_vs_mphi_yukawa-combined.{extension}"
        fig.savefig(path, bbox_inches="tight", pad_inches=0.12, dpi=180 if extension == "png" else None)
        print(f"[OK] wrote {path}")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot-only", action="store_true", help="read the CSV cache and skip EOS")
    parser.add_argument("--masses", default=",".join(str(value) for value in DEFAULT_MASSES))
    parser.add_argument("--yukawas", default=",".join(str(value) for value in DEFAULT_YUKAWAS))
    parser.add_argument("--tree-maker-base", default=DEFAULT_TREE_MAKER_BASE)
    parser.add_argument("--max-files", type=int, default=0, help="files per sample; 0 means all")
    parser.add_argument("--workers", type=int, default=12, help="parallel TreeMaker file readers")
    parser.add_argument("--csv-cache", type=Path, default=DEFAULT_CSV_CACHE)
    parser.add_argument("--root-cache", type=Path, default=DEFAULT_ROOT_CACHE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--preliminary",
        action="store_true",
        help="add the Preliminary qualifier to the CMS label",
    )
    args = parser.parse_args()
    args.masses = parse_csv_numbers(args.masses, int)
    args.yukawas = parse_csv_numbers(args.yukawas, float)

    if args.plot_only:
        rows = read_cache(args.csv_cache)
    else:
        rows = compute(args)
        write_caches(rows, args.csv_cache, args.root_cache)
    for yukawa in args.yukawas:
        draw_plot(rows, yukawa, args.output_dir, preliminary=args.preliminary)
    draw_combined_plot(rows, args.output_dir, preliminary=args.preliminary)


if __name__ == "__main__":
    main()
