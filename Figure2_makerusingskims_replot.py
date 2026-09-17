#!/usr/bin/env python3
"""Replot cached Figure-2 MET/ST normalized stacks without reading EOS skims.

The script reads only the existing ROOT cache produced by the skim builder:
  shapes/<era>/groups/<process>/{MET,ST}
  shapes/<era>/signals/<signal>/{MET,ST}

It intentionally draws exactly the five legacy signal benchmarks by default
and uses the exact single two-column MC-only legend layout from the legacy
plotStack2 script. No event trees and no EOS files are opened.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import sys
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetLineStyleString(4, "28 8 4 8")
ROOT.TH1.AddDirectory(False)


# -----------------------------------------------------------------------------
# Cache and display configuration
# -----------------------------------------------------------------------------
DEFAULT_CACHE = "Figure2_plots/noWNAE_selected_MC_MET_ST_Run2.root"
DEFAULT_PLOT_DIR = "Figure2_plots_legacy_legend"
DEFAULT_ERAS = ["2016", "2017", "2018", "Run2"]

# Exactly the five signals used by the legacy plotting command.  Extra cached
# signals are deliberately ignored unless the user explicitly passes --signals.
# Reference order, colors, and line styles match the approved layout:
# cyan 2000/0.1, brown dotted 2000/0.3, dark-green 2000/0.7,
# green dashed 4000/0.3, dark-blue 600/0.3.
DEFAULT_SIGNAL_DIRS = [
    "t-channel_mMed-2000_mDark-20_rinv-0p1_alpha-peak_yukawa-1",
    "t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
    "t-channel_mMed-2000_mDark-20_rinv-0p7_alpha-peak_yukawa-1",
    "t-channel_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
    "t-channel_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
]

VARIABLES = {
    "MET": {
        "title": "p_{T}^{miss} [GeV]",
        "nbins": 500,
        "xmin": 200.0,
        "xmax": 2000.0,
    },
    "ST": {
        "title": "S_{T} [GeV]",
        "nbins": 500,
        "xmin": 1300.0,
        "xmax": 5000.0,
    },
}

PROC_LABEL = {
    "QCD": "QCD",
    "TTJets": "t#bar{t}+jets",
    "WJetsToLNu": "W#rightarrowl#nu+jets",
    "ZJetsToNuNu": "Z#rightarrow#nu#nu+jets",
    "ST": "Single top",
}

PROC_COLOR = {
    "QCD": ROOT.TColor.GetColor("#9c9ca1"),
    "TTJets": ROOT.TColor.GetColor("#7a21dd"),
    "WJetsToLNu": ROOT.TColor.GetColor("#e42536"),
    "ZJetsToNuNu": ROOT.TColor.GetColor("#f89c20"),
    "ST": ROOT.TColor.GetColor("#5790fc"),
}

STACK_ORDER = ["ST", "ZJetsToNuNu", "WJetsToLNu", "TTJets", "QCD"]
LEGEND_ORDER = ["QCD", "TTJets", "WJetsToLNu", "ZJetsToNuNu", "ST"]

SIGNAL_LINE_COLORS = [
    ROOT.TColor.GetColor("#118AB2"),  # m_phi=2000, r_inv=0.1 -- teal-blue, colorblind-safe (was low-contrast light cyan #92dadd)
    ROOT.TColor.GetColor("#6b3e26"),  # m_phi=2000, r_inv=0.3
    ROOT.TColor.GetColor("#0b3d02"),  # m_phi=2000, r_inv=0.7
    ROOT.TColor.GetColor("#228833"),  # m_phi=4000, r_inv=0.3
    ROOT.TColor.GetColor("#1f4e79"),  # m_phi=600,  r_inv=0.3
]
SIGNAL_LINE_STYLES = [4, 1, 2, 1, 2]


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def parse_csv(text: str) -> List[str]:
    return [item.strip() for item in text.split(",") if item.strip()]


def divisor_generator(number: int) -> Iterable[int]:
    large: List[int] = []
    for divisor in range(1, int(math.sqrt(number)) + 1):
        if number % divisor == 0:
            yield divisor
            if divisor * divisor != number:
                large.append(number // divisor)
    yield from reversed(large)


def rebin_calc(nbins: int, target_bins: int = 40) -> int:
    """Same legacy rebinCalc logic: 500 bins -> factor 10 -> 50 display bins."""
    desired = nbins / float(target_bins)
    divisors = list(divisor_generator(nbins))
    return max(1, int(min(divisors, key=lambda value: abs(value - desired))))


def clone_hist(histogram: ROOT.TH1, name: str) -> ROOT.TH1:
    clone = histogram.Clone(name)
    clone.SetDirectory(0)
    return clone


def rebin_for_display(histogram: ROOT.TH1, variable: str, name: str) -> ROOT.TH1:
    clone = clone_hist(histogram, name)
    factor = rebin_calc(VARIABLES[variable]["nbins"], 40)
    if factor > 1:
        clone = clone.Rebin(factor, f"{name}_rebin{factor}")
        clone.SetDirectory(0)
    return clone


def get_hist(root_file: ROOT.TFile, path: str, clone_name: str) -> Optional[ROOT.TH1]:
    obj = root_file.Get(path)
    if not obj or not obj.InheritsFrom("TH1"):
        return None
    return clone_hist(obj, clone_name)


def cached_signal_names(root_file: ROOT.TFile, era: str) -> List[str]:
    directory = root_file.GetDirectory(f"shapes/{era}/signals")
    if not directory:
        return []
    return sorted(key.GetName() for key in directory.GetListOfKeys())


def choose_signals(
    root_file: ROOT.TFile,
    era: str,
    requested: Sequence[str],
) -> List[str]:
    """Select only exact legacy signal names unless the user overrides them."""
    available = set(cached_signal_names(root_file, era))
    chosen: List[str] = []

    for token in requested:
        if token in available:
            chosen.append(token)
            continue

        matches = sorted(name for name in available if token.lower() in name.lower())
        if len(matches) == 1:
            chosen.append(matches[0])
        elif not matches:
            print(f"[WARN] {era}: cached signal not found: {token}")
        else:
            print(
                f"[WARN] {era}: token '{token}' matches {len(matches)} cached signals; "
                "use the full directory name to select exactly one."
            )

    # Preserve requested order and remove any accidental duplicate token.
    out: List[str] = []
    seen = set()
    for name in chosen:
        if name not in seen:
            out.append(name)
            seen.add(name)
    return out


def signal_label(signal_name: str) -> str:
    mmed = re.search(r"mMed-([0-9]+(?:p[0-9]+)?)", signal_name, re.IGNORECASE)
    rinv = re.search(r"rinv-([0-9]+(?:p[0-9]+)?)", signal_name, re.IGNORECASE)
    if not (mmed and rinv):
        return signal_name
    return (
        f"m_{{#phi}} = {mmed.group(1).replace('p', '.')} GeV, "
        f"r_{{inv}} = {rinv.group(1).replace('p', '.')}"
    )


def style_background(histogram: ROOT.TH1, process: str) -> None:
    color = PROC_COLOR[process]
    histogram.SetFillColor(color)
    histogram.SetFillStyle(1001)
    histogram.SetLineColor(color)
    histogram.SetLineWidth(0)
    histogram.SetMarkerSize(0)


def style_signal(histogram: ROOT.TH1, index: int) -> None:
    histogram.SetLineColor(SIGNAL_LINE_COLORS[index])
    histogram.SetLineStyle(SIGNAL_LINE_STYLES[index])
    histogram.SetLineWidth(3)
    histogram.SetFillStyle(0)
    histogram.SetMarkerSize(0)


def draw_cms_label() -> None:
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextColor(ROOT.kBlack)

    latex.SetTextFont(61)
    latex.SetTextSize(0.060)
    latex.SetTextAlign(11)
    latex.DrawLatex(0.16, 0.94, "CMS")

    latex.SetTextFont(52)
    latex.SetTextSize(0.040)
    latex.DrawLatex(0.290, 0.94, "Simulation")

    latex.SetTextFont(42)
    latex.SetTextSize(0.045)
    latex.SetTextAlign(31)
    latex.DrawLatex(0.95, 0.94, "(13 TeV)")


def draw_common_signal_text() -> None:
    """Exact annotation placement from the legacy MC-only stack plotter."""
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.030)
    latex.SetTextColor(ROOT.kBlack)
    latex.SetTextAlign(12)
    latex.DrawLatex(0.605, 0.700, "m_{dark} = 20 GeV, #lambda = 1")


def make_legacy_mc_only_legend(
    process_hists: Dict[str, ROOT.TH1],
    signal_hists: List[Tuple[str, ROOT.TH1]],
) -> ROOT.TLegend:
    """
    Match the legacy `plot_paper_stack` MC-only legend exactly:
    one TLegend, two columns, with background entries in the left column and
    signal entries in the right column on the same rows.
    """
    legend = ROOT.TLegend(0.19, 0.730, 0.93, 0.895)
    legend.SetNColumns(2)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.030)
    legend.SetTextFont(42)
    legend.SetColumnSeparation(0.04)
    legend.SetEntrySeparation(0.028)

    background_entries = [
        (process_hists[process], PROC_LABEL[process], "F")
        for process in LEGEND_ORDER
        if process in process_hists
    ]
    signal_entries = [
        (hist, signal_label(signal_name), "L")
        for signal_name, hist in signal_hists
    ]

    rows = max(len(background_entries), len(signal_entries))
    for row in range(rows):
        legend.AddEntry(
            *(
                background_entries[row]
                if row < len(background_entries)
                else ("", "", "")
            )
        )
        legend.AddEntry(
            *(
                signal_entries[row]
                if row < len(signal_entries)
                else ("", "", "")
            )
        )

    return legend


def draw_plot(
    root_file: ROOT.TFile,
    era: str,
    variable: str,
    requested_signals: Sequence[str],
    output_base: str,
    make_linear: bool,
) -> None:
    process_hists: Dict[str, ROOT.TH1] = {}
    for process in STACK_ORDER:
        raw = get_hist(
            root_file,
            f"shapes/{era}/groups/{process}/{variable}",
            f"{era}_{process}_{variable}",
        )
        if not raw:
            continue
        hist = rebin_for_display(raw, variable, f"display_{era}_{process}_{variable}")
        if hist.Integral() <= 0.0:
            continue
        style_background(hist, process)
        process_hists[process] = hist

    if not process_hists:
        print(f"[WARN] {era}/{variable}: no normalized background histograms found.")
        return

    signal_hists: List[Tuple[str, ROOT.TH1]] = []
    selected = choose_signals(root_file, era, requested_signals)
    for index, signal_name in enumerate(selected):
        raw = get_hist(
            root_file,
            f"shapes/{era}/signals/{signal_name}/{variable}",
            f"{era}_{signal_name}_{variable}",
        )
        if not raw:
            continue
        hist = rebin_for_display(raw, variable, f"display_{era}_signal_{index}_{variable}")
        if hist.Integral() <= 0.0:
            continue
        style_signal(hist, index)
        signal_hists.append((signal_name, hist))

    modes = [True, False] if make_linear else [True]
    for log_y in modes:
        tag = "log_normalized" if log_y else "linear_normalized"
        canvas = ROOT.TCanvas(f"c_{era}_{variable}_{tag}", "", 800, 800)
        canvas.cd()
        ROOT.gPad.SetLeftMargin(0.16)
        ROOT.gPad.SetRightMargin(0.05)
        ROOT.gPad.SetTopMargin(0.08)
        ROOT.gPad.SetBottomMargin(0.12)
        ROOT.gPad.SetTicks(1, 1)
        ROOT.gPad.SetLogy(log_y)

        stack = ROOT.THStack(f"stack_{era}_{variable}_{tag}", "")
        total: Optional[ROOT.TH1] = None
        for process in STACK_ORDER:
            hist = process_hists.get(process)
            if not hist:
                continue
            stack.Add(hist)
            if total is None:
                total = clone_hist(hist, f"total_{era}_{variable}_{tag}")
            else:
                total.Add(hist)

        if total is None:
            canvas.Close()
            continue

        axis = clone_hist(total, f"axis_{era}_{variable}_{tag}")
        axis.Reset("ICESM")
        axis.SetStats(0)
        axis.SetTitle("")
        axis.GetXaxis().SetTitle(VARIABLES[variable]["title"])
        axis.GetYaxis().SetTitle("Arbitrary units")
        axis.GetXaxis().SetTitleSize(0.05)
        axis.GetXaxis().SetLabelSize(0.04)
        axis.GetYaxis().SetTitleSize(0.05)
        axis.GetYaxis().SetLabelSize(0.04)
        axis.GetYaxis().SetLabelOffset(0.01)
        axis.GetYaxis().SetTitleOffset(1.30)
        axis.GetXaxis().SetRangeUser(VARIABLES[variable]["xmin"], VARIABLES[variable]["xmax"])

        ymax = total.GetMaximum()
        for _, hist in signal_hists:
            ymax = max(ymax, hist.GetMaximum())
        ymax = max(ymax, 1.0)

        if log_y:
            axis.SetMinimum(1.0e-5)
            axis.SetMaximum(ymax * 60.0)
        else:
            axis.SetMinimum(0.0)
            axis.SetMaximum(ymax * 1.5)

        axis.Draw("hist")
        stack.Draw("hist same")
        for _, hist in signal_hists:
            hist.Draw("hist same")

        legend = make_legacy_mc_only_legend(process_hists, signal_hists)
        legend.Draw()
        if signal_hists:
            draw_common_signal_text()
        draw_cms_label()
        ROOT.gPad.RedrawAxis()

        os.makedirs(os.path.dirname(output_base), exist_ok=True)
        canvas.SaveAs(f"{output_base}_{tag}.pdf")
        canvas.SaveAs(f"{output_base}_{tag}.png")
        canvas.Close()
        print(f"[OK] Wrote {output_base}_{tag}.pdf/.png")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Replot cached Figure-2 MET/ST shapes with the legacy two-column MC-only legend; no EOS input is read.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", default=DEFAULT_CACHE, help="Existing compact Figure-2 ROOT cache")
    parser.add_argument("--plot-dir", default=DEFAULT_PLOT_DIR, help="Directory for PDF/PNG output")
    parser.add_argument("--eras", default=",".join(DEFAULT_ERAS), help="Comma-separated cached eras")
    parser.add_argument(
        "--signals",
        default=",".join(DEFAULT_SIGNAL_DIRS),
        help="Exactly selected cached signal directories, in desired drawing order",
    )
    parser.add_argument("--also-linear", action="store_true", help="Also write linear-y versions")
    parser.add_argument("--list-cache-signals", action="store_true", help="List cached signal directories and exit")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    eras = parse_csv(args.eras)
    signals = parse_csv(args.signals)

    root_file = ROOT.TFile.Open(args.input, "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Could not open ROOT cache: {args.input}")

    try:
        if args.list_cache_signals:
            for era in eras:
                print(f"\n[{era}] cached signals")
                for name in cached_signal_names(root_file, era):
                    print(f"  {name}")
            return 0

        for era in eras:
            for variable in ("MET", "ST"):
                output_base = os.path.join(args.plot_dir, era, "mc_signal_only", variable)
                draw_plot(
                    root_file=root_file,
                    era=era,
                    variable=variable,
                    requested_signals=signals,
                    output_base=output_base,
                    make_linear=args.also_linear,
                )
    finally:
        root_file.Close()
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[FATAL] {exc}", file=sys.stderr)
        raise
