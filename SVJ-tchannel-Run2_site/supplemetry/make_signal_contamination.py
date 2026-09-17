#!/usr/bin/env python3
"""Make baseline-signal contamination tables for ABCD control regions."""
import argparse
import contextlib
import csv
import io
import os
import sys
from types import SimpleNamespace

import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np


HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(os.path.dirname(HERE))
POSTFITS = os.path.join(BASE, "t-channel_plotting_scripts", "postfits")
DEFAULT_OUT = os.path.join(HERE, "assets_v2", "signal_contamination")


def _load_postfit_helpers():
    sys.path.insert(0, POSTFITS)
    from postfit_plotter import (  # noqa: WPS433
        _cli_mapping,
        __combine_histograms_years,
        __get_bkg_histograms,
        __get_signal_histograms,
    )

    return _cli_mapping, __combine_histograms_years, __get_bkg_histograms, __get_signal_histograms


def _collect_rows():
    _cli_mapping, combine_years, get_bkg, get_signal = _load_postfit_helpers()
    abcd = SimpleNamespace(remap_back_postfit=False, map_regions={region: region for region in ("A", "B", "C", "D")})
    rows = []

    for tagger_label, subdir in (("ParticleNet", "PNET"), ("WNAE", "WNAE")):
        fit_file = os.path.join(
            POSTFITS,
            subdir,
            "fitDiagnosticscombinedFitDiagnostics_mMed2000_mDark20_rinv0p3_yukawa1_Run2.root",
        )
        card = os.path.join(
            POSTFITS,
            subdir,
            "combinedFitDiagnostics_mMed2000_mDark20_rinv0p3_yukawa1_Run2.txt",
        )
        mapping = _cli_mapping(card)
        by_category = {}

        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            for category in mapping:
                base, year = (category.split("Y", 1) + ["Run2"])[:2] if "Y" in category else (category, "Run2")
                by_category.setdefault(base, {})[year] = {
                    "prefit_signal": get_signal(fit_file, "prefit", mapping, category, abcd),
                    "prefit_bkg": get_bkg(fit_file, "prefit", mapping, category, abcd),
                }

            combined = {
                category: {
                    "signal": combine_years(years, "prefit", "signal"),
                    "background": combine_years(years, "prefit", "bkg"),
                }
                for category, years in by_category.items()
            }

        categories = ["0SVJ", "1SVJ", "2SVJ", "3PSVJ"]
        for region in ("B", "C", "D"):
            inclusive_signal = 0.0
            inclusive_background = 0.0
            for category in categories:
                signal = combined[category]["signal"][region].Integral()
                background = combined[category]["background"][region].Integral()
                inclusive_signal += signal
                inclusive_background += background
                rows.append(_row(tagger_label, region, category, signal, background))
            rows.append(_row(tagger_label, region, "Inclusive", inclusive_signal, inclusive_background))

    return rows


def _row(tagger, region, category, signal, background):
    contamination = 100.0 * signal / background if background else 0.0
    return {
        "tagger": tagger,
        "region": region,
        "category": category,
        "signal": signal,
        "background": background,
        "contamination_percent": contamination,
    }


def _format_yield(value):
    if abs(value) >= 1.0e4:
        return f"{value:.2e}"
    if abs(value) >= 100.0:
        return f"{value:.1f}"
    return f"{value:.2f}"


def _category_tex(category):
    return {
        "0SVJ": r"$0$",
        "1SVJ": r"$1$",
        "2SVJ": r"$2$",
        "3PSVJ": r"$\geq 3$",
        "Inclusive": "Inclusive",
    }[category]


def _write_csv(rows, out_dir):
    path = os.path.join(out_dir, "baseline_signal_contamination_bcd_prefit.csv")
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return path


def _write_tex(rows, out_dir):
    path = os.path.join(out_dir, "baseline_signal_contamination_bcd_prefit.tex")
    with open(path, "w") as handle:
        handle.write(r"""\begin{table}[H]
\centering
\small
\resizebox{\textwidth}{!}{%
\begin{tabular}{ll l r r r}
\hline
Tagger & Region & $n_{\mathrm{SVJ}}$ category & Signal yield & Prefit background yield & $S/B$ [\%] \\
\hline
""")
        for tagger in ("ParticleNet", "WNAE"):
            for region in ("B", "C", "D"):
                selected = [row for row in rows if row["tagger"] == tagger and row["region"] == region]
                for index, row in enumerate(selected):
                    tagger_text = tagger if region == "B" and index == 0 else ""
                    region_text = region if index == 0 else ""
                    category = _category_tex(row["category"])
                    signal = _format_yield(row["signal"])
                    background = _format_yield(row["background"])
                    contamination = f"{row['contamination_percent']:.3f}"
                    handle.write(
                        f"{tagger_text} & {region_text} & {category} & {signal} & {background} & {contamination} \\\\\n"
                    )
                handle.write(r"\hline" + "\n")
        handle.write(r"""\end{tabular}
}
\caption{Expected signal contamination in the ABCD control regions for the baseline signal point $m_{\Phi}=2000$ GeV, $m_{\mathrm{dark}}=20$ GeV, $r_{\mathrm{inv}}=0.3$, and $\lambda=1$. Both the signal and background yields are the nominal prefit expectations from simulation. The contamination is shown as $S/B$.}
\label{tab:baseline-signal-contamination-bcd}
\end{table}
""")
    return path


def _write_plot(rows, out_dir):
    path = os.path.join(out_dir, "baseline_signal_contamination_bcd_prefit_by_nsvj.pdf")
    hep.style.use(hep.style.CMS)
    fig, ax = plt.subplots(figsize=(12.0, 7.2))

    regions = ["B", "C", "D"]
    categories = ["0SVJ", "1SVJ", "2SVJ", "3PSVJ"]
    category_labels = [r"$0$", r"$1$", r"$2$", r"$\geq 3$"]
    group_gap = 0.65
    x = []
    xticklabels = []
    for region_index, region in enumerate(regions):
        base = region_index * (len(categories) + group_gap)
        for category_index, label in enumerate(category_labels):
            x.append(base + category_index)
            xticklabels.append(label)
    x = np.asarray(x, dtype=float)
    width = 0.36
    colors = {"ParticleNet": "#118AB2", "WNAE": "#C1121F"}

    for offset, tagger in ((-width / 2, "ParticleNet"), (width / 2, "WNAE")):
        values = [
            next(
                row["contamination_percent"]
                for row in rows
                if row["tagger"] == tagger and row["region"] == region and row["category"] == "Inclusive"
            )
            if category == "Inclusive"
            else next(
                row["contamination_percent"]
                for row in rows
                if row["tagger"] == tagger and row["region"] == region and row["category"] == category
            )
            for region in regions
            for category in categories
        ]
        ax.bar(x + offset, values, width=width, label=tagger, color=colors[tagger], edgecolor="black", linewidth=1.0)

    ax.text(0.00, 1.005, "CMS", transform=ax.transAxes, ha="left", va="bottom", fontsize=28, fontweight="bold")
    ax.text(0.13, 1.005, "Simulation Preliminary", transform=ax.transAxes, ha="left", va="bottom", fontsize=22, style="italic")
    ax.text(1.00, 1.005, r"138 fb$^{-1}$ (13 TeV)", transform=ax.transAxes, ha="right", va="bottom", fontsize=22)
    ax.set_xticks(x)
    ax.set_xticklabels(xticklabels)
    ax.set_ylabel("Signal contamination $S/B$ [%]")
    ax.set_xlabel(r"$n_{\mathrm{SVJ}}$ category")
    ymax = max(row["contamination_percent"] for row in rows if row["category"] != "Inclusive")
    ax.set_ylim(0, ymax * 1.45 if ymax > 0 else 1.0)
    for region_index, region in enumerate(regions):
        center = region_index * (len(categories) + group_gap) + (len(categories) - 1) / 2.0
        ax.text(center, -0.16, f"Region {region}", transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=18)
        if region_index > 0:
            ax.axvline(region_index * (len(categories) + group_gap) - group_gap / 2.0, color="black", linewidth=1.0, alpha=0.35)
    ax.legend(loc="upper right", frameon=True, edgecolor="white", framealpha=1)
    ax.text(
        0.03,
        0.92,
        r"$m_{\Phi}=2000$ GeV, $m_{\mathrm{dark}}=20$ GeV" + "\n" + r"$r_{\mathrm{inv}}=0.3$, $\lambda=1$",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=16,
    )
    fig.subplots_adjust(left=0.11, right=0.98, bottom=0.20, top=0.90)
    fig.savefig(path)
    plt.close(fig)
    return path


def main():
    parser = argparse.ArgumentParser(description="Make baseline-signal contamination assets.")
    parser.add_argument("-o", "--out-dir", default=DEFAULT_OUT)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    rows = _collect_rows()
    print(f"[OK] wrote {_write_csv(rows, args.out_dir)}")
    print(f"[OK] wrote {_write_tex(rows, args.out_dir)}")
    print(f"[OK] wrote {_write_plot(rows, args.out_dir)}")


if __name__ == "__main__":
    main()
