#!/usr/bin/env python3
"""
Build DNN score plots from ABCD_scores ROOT files.

This is intentionally separate from EventDNNInput_makerusingskims.py so the
existing input-variable maker and any current stash/worktree state stay intact.

Example:
  source condor/initCondor.sh
  python3 -u EventDNNScore_makerusingscores.py \
    --background-base /eos/uscms/store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/DNN_scores_and_ABCD_histograms/DNN_scores/scores_backgrounds_GapJetVeto \
    --data-base /eos/uscms/store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/DNN_scores_and_ABCD_histograms/DNN_scores/scores_data_GapJetVeto \
    --signal-base /eos/uscms/store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/DNN_scores_and_ABCD_histograms/DNN_scores/scores_signals_GapJetVeto \
    --years 2016 \
    --output EventDNNScore_plots/dnn_scores.root \
    --plot-dir EventDNNScore_plots

To redraw from the ROOT cache:
  python3 -u EventDNNScore_makerusingscores.py \
    --plot-only \
    --output EventDNNScore_plots/dnn_scores.root \
    --plot-dir EventDNNScore_plots
"""

from __future__ import annotations

import argparse
import math
import os
import re
import sys
from array import array
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import uproot

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.SetDefaultSumw2(True)
ROOT.TH1.AddDirectory(False)


DEFAULT_SCORE_BASE = (
    "/eos/uscms/store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/"
    "DNN_scores_and_ABCD_histograms/DNN_scores"
)
DEFAULT_BACKGROUND_BASE = f"{DEFAULT_SCORE_BASE}/scores_backgrounds_GapJetVeto"
DEFAULT_DATA_BASE = f"{DEFAULT_SCORE_BASE}/scores_data_GapJetVeto"
DEFAULT_SIGNAL_BASE = f"{DEFAULT_SCORE_BASE}/scores_signals_GapJetVeto"
DEFAULT_SIGNAL_DIRS = [
    "t-channel_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
    "t-channel_mMed-2000_mDark-20_rinv-0p1_alpha-peak_yukawa-1",
    "t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
    # "t-channel_mMed-2000_mDark-20_rinv-0p7_alpha-peak_yukawa-1",
    # "t-channel_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
]
DISABLED_SIGNAL_TOKENS = [
    "t-channel_mMed-2000_mDark-20_rinv-0p7_alpha-peak_yukawa-1",
    "t-channel_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
]

SIGNAL_SCALE_FACTORS = {
    ("600", "0.3"): 10.0,
    ("2000", "0.1"): 500.0,
    ("2000", "0.3"): 500.0,
    ("2000", "0.7"): 1000.0,
    ("4000", "0.3"): 1000.0,
}

LUMI_PB = {
    "2016": 36.31 * 1000.0,
    "2017": 42.07 * 1000.0,
    "2018": 59.56 * 1000.0,
}

SCORE_VARIABLE = "DNNScore"
SCORE_FILL_BINNING = "equal20"
VARIABLES = {
    SCORE_VARIABLE: {
        "source": "scores",
        "weight": "weights",
        "title": "Event classifier score",
        "x_title": "Event classifier score",
        "nbins": 20,
        "xmin": 0.0,
        "xmax": 1.0,
        "wp": 0.85,
        "wp_side": "right",
        "wp_line_y_fraction": 0.900,
        "wp_arrow_y_fraction": 0.600,
    }
}
PLOT_TARGET_BINS = 40

PROCESS_ORDER = ["QCD", "ST", "TTJets", "WJetsToLNu", "ZJetsToNuNu"]
STACK_ORDER = ["ST", "ZJetsToNuNu", "WJetsToLNu", "TTJets", "QCD"]
LEGEND_ORDER = ["QCD", "TTJets", "WJetsToLNu", "ZJetsToNuNu", "ST"]
DEFAULT_DATA_SAMPLES = ["HTMHT", "JetHT", "MET"]

PROC_LABEL = {
    "QCD": "QCD",
    "TTJets": "t#bar{t}+jets",
    "ZJetsToNuNu": "Z#rightarrow#nu#nu+jets",
    "WJetsToLNu": "W#rightarrowl#nu+jets",
    "ST": "Single top",
}

PROC_COLOR = {
    "QCD": ROOT.TColor.GetColor("#9c9ca1"),
    "TTJets": ROOT.TColor.GetColor("#7a21dd"),
    "WJetsToLNu": ROOT.TColor.GetColor("#e42536"),
    "ZJetsToNuNu": ROOT.TColor.GetColor("#f89c20"),
    "ST": ROOT.TColor.GetColor("#5790fc"),
}

SIGNAL_LINE_COLORS = [
    ROOT.TColor.GetColor("#92dadd"),
    ROOT.TColor.GetColor("#6b3e26"),
    ROOT.TColor.GetColor("#0b3d02"),
    ROOT.TColor.GetColor("#228833"),
    ROOT.TColor.GetColor("#1f4e79"),
]
SIGNAL_LINE_STYLES = [1, 3, 1, 2, 1]


@dataclass
class HistogramAccumulator:
    edges: np.ndarray
    sumw: np.ndarray = field(init=False)
    sumw2: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        self.sumw = np.zeros(len(self.edges) - 1, dtype=np.float64)
        self.sumw2 = np.zeros(len(self.edges) - 1, dtype=np.float64)

    @classmethod
    def for_variable(cls, variable: str) -> "HistogramAccumulator":
        cfg = VARIABLES[variable]
        if variable == SCORE_VARIABLE and SCORE_FILL_BINNING == "wpedge20":
            return cls(np.asarray(list(np.linspace(0.0, 0.85, 20)) + [1.0], dtype=np.float64))
        return cls(np.linspace(cfg["xmin"], cfg["xmax"], cfg["nbins"] + 1))

    def fill(self, values, weights) -> None:
        values = np.asarray(values, dtype=np.float64)
        weights = np.asarray(weights, dtype=np.float64)
        mask = (
            np.isfinite(values)
            & np.isfinite(weights)
            & (values >= self.edges[0])
            & (values < self.edges[-1])
        )
        if not np.any(mask):
            return
        values = values[mask]
        weights = weights[mask]
        self.sumw += np.histogram(values, bins=self.edges, weights=weights)[0]
        self.sumw2 += np.histogram(values, bins=self.edges, weights=weights * weights)[0]

    def add(self, other: "HistogramAccumulator") -> None:
        if not np.array_equal(self.edges, other.edges):
            raise ValueError("Attempted to add incompatible histograms")
        self.sumw += other.sumw
        self.sumw2 += other.sumw2

    def integral(self) -> float:
        return float(np.sum(self.sumw))


def empty_accumulators() -> Dict[str, HistogramAccumulator]:
    return {variable: HistogramAccumulator.for_variable(variable) for variable in VARIABLES}


def add_accumulators(target, source) -> None:
    for variable in VARIABLES:
        target[variable].add(source[variable])


def parse_csv(value: str) -> List[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def local_eos_path(path: str) -> Path:
    path = os.path.expanduser(path).rstrip("/")
    if path.startswith("/store/"):
        path = "/eos/uscms" + path
    return Path(path)


def score_file_for_process(background_base: str, year: str, process: str) -> Path:
    return local_eos_path(background_base) / year / "nominal" / f"ABCD_scores_{process}.root"


def score_file_for_sample(base: str, year: str, sample: str) -> Path:
    return local_eos_path(base) / year / "nominal" / f"ABCD_scores_{sample}.root"


def stream_score_file(
    path: Path,
    accumulators,
    step_size: str,
    strict_positive_weights: bool = False,
    is_data: bool = False,
) -> int:
    with uproot.open(path) as root_file:
        if "ABCD_scores" not in root_file:
            raise KeyError(f"{path} has no ABCD_scores tree")
        tree = root_file["ABCD_scores"]
        available = set(tree.keys())
        missing = ["scores"] if "scores" not in available else []
        if not is_data and "weights" not in available:
            missing.append("weights")
        if missing:
            raise KeyError(f"{path} lacks required branch(es): {missing}")

        entries = int(tree.num_entries)
        expressions = ["scores"] + (["weights"] if "weights" in available else [])
        for arrays in tree.iterate(expressions=expressions, step_size=step_size, library="np"):
            scores = arrays["scores"]
            weights = arrays["weights"] if "weights" in arrays else np.ones(len(scores), dtype=np.float64)
            if strict_positive_weights:
                keep = weights > 0.0
                scores = scores[keep]
                weights = weights[keep]
            accumulators[SCORE_VARIABLE].fill(scores, weights)
    return entries


def read_backgrounds(background_base: str, years: Sequence[str], step_size: str, strict: bool, strict_positive_weights: bool):
    group_hists: Dict[Tuple[str, str], Dict[str, HistogramAccumulator]] = {}
    run2_hists = {process: empty_accumulators() for process in PROCESS_ORDER}
    failures = []

    for year in years:
        print(f"\n[INFO] Background score files {year}")
        for process in PROCESS_ORDER:
            path = score_file_for_process(background_base, year, process)
            group_hists[(year, process)] = empty_accumulators()
            if not path.exists():
                message = f"{year}/{process}: missing {path}"
                if strict:
                    raise RuntimeError(message)
                print(f"[WARN] {message}")
                continue
            try:
                entries = stream_score_file(path, group_hists[(year, process)], step_size, strict_positive_weights)
                add_accumulators(run2_hists[process], group_hists[(year, process)])
                yield_sum = group_hists[(year, process)][SCORE_VARIABLE].integral()
                print(f"  {process:12s} entries={entries:10,d} yield={yield_sum:.8g}")
            except Exception as exc:
                message = f"{year}/{process}: {exc}"
                failures.append(message)
                if strict:
                    raise RuntimeError(message) from exc
                print(f"[WARN] {message}")

    return group_hists, run2_hists, failures


def read_data(data_base: str, years: Sequence[str], samples: Sequence[str], step_size: str, strict: bool, strict_positive_weights: bool):
    data_hists: Dict[Tuple[str, str], Dict[str, HistogramAccumulator]] = {}
    run2_data_hists = empty_accumulators()
    failures = []
    found_any = False

    if not data_base or not samples:
        return data_hists, None, failures

    for year in years:
        print(f"\n[INFO] Data score files {year}")
        for sample in samples:
            path = score_file_for_sample(data_base, year, sample)
            if not path.exists():
                message = f"{year}/{sample}: missing {path}"
                if strict:
                    raise RuntimeError(message)
                print(f"[WARN] {message}")
                continue
            accs = empty_accumulators()
            try:
                entries = stream_score_file(
                    path,
                    accs,
                    step_size,
                    strict_positive_weights=strict_positive_weights,
                    is_data=True,
                )
                data_hists[(year, sample)] = accs
                add_accumulators(run2_data_hists, accs)
                found_any = True
                print(f"  {sample:12s} entries={entries:10,d} yield={accs[SCORE_VARIABLE].integral():.8g}")
            except Exception as exc:
                message = f"{year}/{sample}: {exc}"
                failures.append(message)
                if strict:
                    raise RuntimeError(message) from exc
                print(f"[WARN] {message}")

    return data_hists, run2_data_hists if found_any else None, failures


def signal_label(signal_name: str) -> str:
    mmed = re.search(r"mMed-([0-9]+(?:p[0-9]+)?)", signal_name, re.IGNORECASE)
    rinv = re.search(r"rinv-([0-9]+(?:p[0-9]+)?)", signal_name, re.IGNORECASE)
    if not (mmed and rinv):
        return signal_name
    return f"m_{{#phi}} = {mmed.group(1).replace('p', '.')} GeV, r_{{inv}} = {rinv.group(1).replace('p', '.')}"


def signal_scale_factor(signal_name: str) -> float:
    mmed = re.search(r"mMed-([0-9]+(?:p[0-9]+)?)", signal_name, re.IGNORECASE)
    rinv = re.search(r"rinv-([0-9]+(?:p[0-9]+)?)", signal_name, re.IGNORECASE)
    if not (mmed and rinv):
        return 1.0
    key = (mmed.group(1).replace("p", "."), rinv.group(1).replace("p", "."))
    return SIGNAL_SCALE_FACTORS.get(key, 1.0)


def scaled_signal_label(signal_name: str, normalized: bool) -> str:
    label = signal_label(signal_name)
    scale = signal_scale_factor(signal_name)
    if normalized or abs(scale - 1.0) < 1.0e-9:
        return label
    return f"{label} (#times {scale:g})"


def signal_name_from_file(path: Path) -> str:
    name = path.name
    name = re.sub(r"^ABCD_scores_(?:nominal_)?", "", name)
    name = re.sub(r"_(?:nominal|lundWeightNom)\.root$", "", name)
    name = re.sub(r"\.root$", "", name)
    return name


def signal_file_matches_weight_variation(path: Path, weight_variation: str) -> bool:
    if weight_variation == "auto":
        return bool(re.search(r"_(?:lundWeightNom|nominal)\.root$", path.name))
    return path.name.endswith(f"_{weight_variation}.root")


def central_signal_priority(path: Path) -> int:
    return 0 if path.name.endswith("_lundWeightNom.root") else 1


def find_signal_files(
    signal_base: str,
    year: str,
    requested: Sequence[str],
    weight_variation: str,
) -> List[Tuple[str, Path]]:
    if not signal_base:
        return []
    base = local_eos_path(signal_base) / year / "nominal"
    if not base.exists():
        print(f"[WARN] {year}: signal score directory not found: {base}")
        return []
    files = sorted(
        (
            path
            for path in base.glob("ABCD_scores*.root")
            if signal_file_matches_weight_variation(path, weight_variation)
        ),
        key=lambda path: (signal_name_from_file(path), central_signal_priority(path), path.name),
    )
    available = []
    seen_available = set()
    for path in files:
        name = signal_name_from_file(path)
        if name in seen_available:
            continue
        available.append((name, path))
        seen_available.add(name)
    if not requested:
        return []
    if len(requested) == 1 and requested[0].lower() in {"all", "*"}:
        return available

    selected = []
    seen = set()
    for token in requested:
        token_lower = token.lower()
        matches = [(name, path) for name, path in available if token_lower == name.lower()]
        if not matches:
            matches = [(name, path) for name, path in available if token_lower in name.lower()]
        if not matches:
            print(f"[WARN] {year}: no signal score file matches '{token}'")
            continue
        for name, path in matches:
            if name not in seen:
                selected.append((name, path))
                seen.add(name)
    return selected


def read_signals(
    signal_base: str,
    years: Sequence[str],
    requested: Sequence[str],
    step_size: str,
    strict_positive_weights: bool,
    signal_weight_variation: str,
):
    signal_hists: Dict[Tuple[str, str], Dict[str, HistogramAccumulator]] = {}
    run2_signal_hists: Dict[str, Dict[str, HistogramAccumulator]] = {}
    if not signal_base or not requested:
        return signal_hists, run2_signal_hists

    for year in years:
        print(f"\n[INFO] Signal score files {year}")
        for signal_name, path in find_signal_files(signal_base, year, requested, signal_weight_variation):
            accs = empty_accumulators()
            entries = stream_score_file(path, accs, step_size, strict_positive_weights)
            signal_hists[(year, signal_name)] = accs
            if signal_name not in run2_signal_hists:
                run2_signal_hists[signal_name] = empty_accumulators()
            add_accumulators(run2_signal_hists[signal_name], accs)
            print(f"  {signal_name:70s} entries={entries:10,d} yield={accs[SCORE_VARIABLE].integral():.8g}")
    return signal_hists, run2_signal_hists


def mkdirs(root_file: ROOT.TFile, directory_path: str):
    directory = root_file
    for part in (item for item in directory_path.split("/") if item):
        next_dir = directory.GetDirectory(part)
        if not next_dir:
            next_dir = directory.mkdir(part)
        directory = next_dir
    return directory


def accumulator_to_histogram(variable: str, accumulator, name: Optional[str] = None, y_title: str = "Events"):
    cfg = VARIABLES[variable]
    hist_name = name or variable
    hist = ROOT.TH1D(
        hist_name,
        f"{hist_name};{cfg['title']};{y_title}",
        len(accumulator.edges) - 1,
        array("d", accumulator.edges.tolist()),
    )
    hist.Sumw2()
    for ibin, (content, variance) in enumerate(zip(accumulator.sumw, accumulator.sumw2), start=1):
        hist.SetBinContent(ibin, float(content))
        hist.SetBinError(ibin, math.sqrt(max(0.0, float(variance))))
    hist.SetDirectory(0)
    return hist


def write_histogram(root_file: ROOT.TFile, directory_path: str, histogram, key_name: str) -> None:
    directory = mkdirs(root_file, directory_path)
    directory.cd()
    histogram.Write(key_name, ROOT.TObject.kOverwrite)


def write_cache(
    output_path: str,
    group_hists,
    run2_hists,
    signal_hists,
    run2_signal_hists,
    data_hists,
    run2_data_hists,
    years: Sequence[str],
) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    output = ROOT.TFile.Open(output_path, "RECREATE")
    if not output or output.IsZombie():
        raise RuntimeError(f"Could not create ROOT output: {output_path}")
    try:
        for (year, process), accumulators in sorted(group_hists.items()):
            for variable in VARIABLES:
                write_histogram(
                    output,
                    f"raw/groups/{year}/{process}",
                    accumulator_to_histogram(variable, accumulators[variable]),
                    variable,
                )
        for process, accumulators in sorted(run2_hists.items()):
            for variable in VARIABLES:
                write_histogram(
                    output,
                    f"raw/groups/Run2/{process}",
                    accumulator_to_histogram(variable, accumulators[variable]),
                    variable,
                )
        for (year, signal_name), accumulators in sorted(signal_hists.items()):
            for variable in VARIABLES:
                write_histogram(
                    output,
                    f"raw/signals/{year}/{signal_name}",
                    accumulator_to_histogram(variable, accumulators[variable]),
                    variable,
                )
        for signal_name, accumulators in sorted(run2_signal_hists.items()):
            for variable in VARIABLES:
                write_histogram(
                    output,
                    f"raw/signals/Run2/{signal_name}",
                    accumulator_to_histogram(variable, accumulators[variable]),
                    variable,
                )
        for (year, sample), accumulators in sorted(data_hists.items()):
            for variable in VARIABLES:
                write_histogram(
                    output,
                    f"raw/data/{year}/{sample}",
                    accumulator_to_histogram(variable, accumulators[variable]),
                    variable,
                )
        for year in years:
            year_data = empty_accumulators()
            found = False
            for (data_year, _), accumulators in data_hists.items():
                if data_year == year:
                    add_accumulators(year_data, accumulators)
                    found = True
            if found:
                for variable in VARIABLES:
                    write_histogram(
                        output,
                        f"raw/data/{year}/observed",
                        accumulator_to_histogram(variable, year_data[variable]),
                        variable,
                    )
        if run2_data_hists:
            for variable in VARIABLES:
                write_histogram(
                    output,
                    "raw/data/Run2/observed",
                    accumulator_to_histogram(variable, run2_data_hists[variable]),
                    variable,
                )

        eras = {
            year: {
                process: group_hists[(year, process)]
                for process in PROCESS_ORDER
                if (year, process) in group_hists
            }
            for year in years
        }
        eras["Run2"] = run2_hists
        for era, processes in eras.items():
            for variable in VARIABLES:
                total = sum(accs[variable].integral() for accs in processes.values())
                if total <= 0.0:
                    continue
                for process, accs in processes.items():
                    hist = accumulator_to_histogram(variable, accs[variable], y_title="Arbitrary units")
                    hist.Scale(1.0 / total)
                    write_histogram(output, f"shapes/{era}/groups/{process}", hist, variable)

        era_signals = {
            year: {name: accs for (signal_year, name), accs in signal_hists.items() if signal_year == year}
            for year in years
        }
        era_signals["Run2"] = run2_signal_hists
        for era, signals in era_signals.items():
            for signal_name, accs in signals.items():
                for variable in VARIABLES:
                    hist = accumulator_to_histogram(variable, accs[variable], y_title="Arbitrary units")
                    integral = hist.Integral()
                    if integral <= 0.0:
                        continue
                    hist.Scale(1.0 / integral)
                    write_histogram(output, f"shapes/{era}/signals/{signal_name}", hist, variable)

        for era in list(years) + ["Run2"]:
            source = run2_data_hists if era == "Run2" else None
            if era != "Run2":
                source = empty_accumulators()
                found = False
                for (data_year, _), accumulators in data_hists.items():
                    if data_year == era:
                        add_accumulators(source, accumulators)
                        found = True
                if not found:
                    source = None
            if not source:
                continue
            for variable in VARIABLES:
                hist = accumulator_to_histogram(variable, source[variable], y_title="Arbitrary units")
                integral = hist.Integral()
                if integral <= 0.0:
                    continue
                hist.Scale(1.0 / integral)
                write_histogram(output, f"shapes/{era}/data/observed", hist, variable)
    finally:
        output.Close()


def divisor_generator(number: int) -> Iterable[int]:
    large = []
    for divisor in range(1, int(math.sqrt(number)) + 1):
        if number % divisor == 0:
            yield divisor
            if divisor * divisor != number:
                large.append(number // divisor)
    yield from reversed(large)


def rebin_factor(variable: str) -> int:
    nbins = int(VARIABLES[variable]["nbins"])
    desired = nbins / float(PLOT_TARGET_BINS)
    divisors = list(divisor_generator(nbins))
    return max(1, int(min(divisors, key=lambda value: abs(value - desired))))


def clone_hist(histogram, name: str):
    out = histogram.Clone(name)
    out.SetDirectory(0)
    return out


def get_hist(root_file, path: str, clone_name: str):
    obj = root_file.Get(path)
    if not obj or not obj.InheritsFrom("TH1"):
        return None
    return clone_hist(obj, clone_name)


def variable_rebin(histogram, name: str, edges: Sequence[float]):
    rebinned = ROOT.TH1D(
        name,
        histogram.GetTitle(),
        len(edges) - 1,
        array("d", list(edges)),
    )
    rebinned.Sumw2()
    for ibin in range(1, histogram.GetNbinsX() + 1):
        center = histogram.GetXaxis().GetBinCenter(ibin)
        target = rebinned.FindBin(center)
        if target < 1 or target > rebinned.GetNbinsX():
            continue
        content = histogram.GetBinContent(ibin)
        error = histogram.GetBinError(ibin)
        rebinned.SetBinContent(target, rebinned.GetBinContent(target) + content)
        rebinned.SetBinError(target, math.hypot(rebinned.GetBinError(target), error))
    rebinned.SetDirectory(0)
    return rebinned


def rebin_for_display(histogram, variable: str, name: str, score_display_binning: str = "equal20"):
    out = clone_hist(histogram, name)
    if variable == SCORE_VARIABLE:
        if score_display_binning == "equal20":
            return out
        if score_display_binning == "wpedge10":
            edges = [0.1 * index for index in range(9)] + [0.85, 1.0]
            return variable_rebin(histogram, f"{name}_wpedge10", edges)
        if score_display_binning == "wpedge20":
            edges = list(np.linspace(0.0, 0.85, 20)) + [1.0]
            return variable_rebin(histogram, f"{name}_wpedge20", edges)
        raise ValueError(f"Unknown score display binning: {score_display_binning}")
    factor = rebin_factor(variable)
    if factor > 1:
        out = out.Rebin(factor, f"{name}_rebin{factor}")
        out.SetDirectory(0)
    return out


def style_background(histogram, process: str) -> None:
    color = PROC_COLOR[process]
    histogram.SetFillColor(color)
    histogram.SetFillStyle(1001)
    histogram.SetLineColor(color)
    histogram.SetLineWidth(0)
    histogram.SetMarkerSize(0)


def style_signal(histogram, index: int) -> None:
    histogram.SetLineColor(SIGNAL_LINE_COLORS[index % len(SIGNAL_LINE_COLORS)])
    histogram.SetLineStyle(SIGNAL_LINE_STYLES[index % len(SIGNAL_LINE_STYLES)])
    histogram.SetLineWidth(3)
    histogram.SetFillStyle(0)
    histogram.SetMarkerSize(0)


def style_data(histogram) -> None:
    histogram.SetLineColor(ROOT.kBlack)
    histogram.SetMarkerColor(ROOT.kBlack)
    histogram.SetMarkerStyle(ROOT.kFullCircle)
    histogram.SetMarkerSize(0.85)
    histogram.SetLineWidth(2)
    histogram.SetFillStyle(0)


def style_uncertainty_band(histogram) -> None:
    histogram.SetFillColor(ROOT.kGray + 1)
    histogram.SetFillStyle(3354)
    histogram.SetLineColor(ROOT.kGray + 2)
    histogram.SetLineWidth(1)
    histogram.SetMarkerSize(0)


def clone_uncertainty_band(histogram, name: str):
    band = clone_hist(histogram, name)
    style_uncertainty_band(band)
    return band


def lumi_label(era: str) -> Optional[str]:
    if era == "Run2":
        return "138 fb^{-1} (13 TeV)"
    elif era in LUMI_PB:
        lumi_fb = LUMI_PB[era] / 1000.0
    else:
        return None
    return f"{lumi_fb:.1f} fb^{{-1}} (13 TeV)"


def draw_cms_label(
    simulation: bool = True,
    era: Optional[str] = None,
    show_lumi: bool = False,
    preliminary: bool = False,
) -> None:
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextColor(ROOT.kBlack)
    latex.SetTextFont(61)
    latex.SetTextSize(0.060)
    latex.SetTextAlign(11)
    label_y = 0.918
    latex.DrawLatex(0.16, label_y, "CMS")
    if simulation:
        latex.SetTextFont(52)
        latex.SetTextSize(0.040)
        latex.DrawLatex(0.255, label_y, "Simulation")
        if preliminary:
            latex.DrawLatex(0.420, label_y, "Preliminary")
    elif preliminary:
        latex.SetTextFont(52)
        latex.SetTextSize(0.040)
        latex.DrawLatex(0.255, label_y, "Preliminary")
    latex.SetTextFont(42)
    latex.SetTextSize(0.045)
    latex.SetTextAlign(31)
    latex.DrawLatex(0.95, label_y, lumi_label(era) if show_lumi and era else "(13 TeV)")


def draw_common_signal_text(y_position: float = 0.700) -> None:
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.030)
    latex.SetTextColor(ROOT.kBlack)
    latex.SetTextAlign(12)
    latex.DrawLatex(0.565, y_position, "m_{dark} = 20 GeV, #kern[-0.18]{#lambda} = 1")


def draw_ratio_y_title(title: str = "Data/Sim") -> None:
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.12)
    latex.SetTextAngle(90)
    latex.SetTextAlign(31)
    latex.DrawLatex(0.085, 0.985, title)


def interpolate_axis_y(ymin: float, ymax: float, fraction: float, log_y: bool) -> float:
    fraction = max(0.0, min(1.0, fraction))
    if log_y and ymin > 0.0 and ymax > ymin:
        return math.exp(math.log(ymin) + fraction * (math.log(ymax) - math.log(ymin)))
    return ymin + fraction * (ymax - ymin)


def draw_working_point_marker(variable: str, axis, log_y: bool):
    cfg = VARIABLES[variable]
    if "wp" not in cfg:
        return None
    x = float(cfg["wp"])
    xmin = float(cfg["xmin"])
    xmax = float(cfg["xmax"])
    if x <= xmin or x >= xmax:
        return None
    ymin = axis.GetMinimum()
    ymax = axis.GetMaximum()
    line_top = interpolate_axis_y(ymin, ymax, float(cfg.get("wp_line_y_fraction", 0.82)), log_y)
    arrow_y = interpolate_axis_y(ymin, ymax, float(cfg.get("wp_arrow_y_fraction", 0.74)), log_y)

    line = ROOT.TLine(x, ymin, x, line_top)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineStyle(ROOT.kDashed)
    line.SetLineWidth(2)
    line.Draw("same")

    arrow = ROOT.TLatex()
    arrow.SetTextColor(ROOT.kBlack)
    arrow.SetTextFont(42)
    arrow.SetTextSize(0.040)
    if cfg.get("wp_side", "right") == "left":
        arrow.SetTextAlign(32)
        arrow.DrawLatex(x, arrow_y, "#leftarrow")
    else:
        arrow.SetTextAlign(12)
        arrow.DrawLatex(x, arrow_y, "#rightarrow")
    return line, arrow


def cached_signal_names(root_file, era: str) -> List[str]:
    directory = root_file.GetDirectory(f"shapes/{era}/signals")
    if not directory:
        return []
    return sorted(key.GetName() for key in directory.GetListOfKeys())


def signal_is_disabled(signal_name: str) -> bool:
    lower_name = signal_name.lower()
    return any(token.lower() in lower_name for token in DISABLED_SIGNAL_TOKENS)


def choose_signals(root_file, era: str, requested: Sequence[str]) -> List[str]:
    available = set(name for name in cached_signal_names(root_file, era) if not signal_is_disabled(name))
    if not requested:
        return []
    if len(requested) == 1 and requested[0].lower() in {"all", "*"}:
        return sorted(available)
    chosen = []
    for token in requested:
        token_lower = token.lower()
        matches = sorted(name for name in available if token_lower == name.lower())
        if not matches:
            matches = sorted(name for name in available if token_lower in name.lower())
        if len(matches) == 1:
            chosen.append(matches[0])
        elif not matches:
            print(f"[WARN] {era}: cached signal not found: {token}")
        else:
            print(f"[WARN] {era}: token '{token}' matches {len(matches)} signals; use full name")
    out = []
    seen = set()
    for name in chosen:
        if name not in seen:
            out.append(name)
            seen.add(name)
    return out


def draw_plot(
    root_file,
    era: str,
    variable: str,
    requested_signals: Sequence[str],
    plot_dir: str,
    normalized: bool,
    also_linear: bool,
    include_data: bool = False,
    include_ratio: bool = False,
    score_display_binning: str = "equal20",
    draw_uncertainty_band: bool = True,
    preliminary_label: bool = False,
) -> None:
    directory_kind = "shapes" if normalized else "raw"
    y_title = "Arbitrary units" if normalized else "Events"
    group_prefix = f"shapes/{era}/groups" if normalized else f"raw/groups/{era}"
    signal_prefix = f"shapes/{era}/signals" if normalized else f"raw/signals/{era}"
    data_prefix = f"shapes/{era}/data" if normalized else f"raw/data/{era}"
    process_hists = {}
    for process in STACK_ORDER:
        raw = get_hist(
            root_file,
            f"{group_prefix}/{process}/{variable}",
            f"{era}_{directory_kind}_{process}_{variable}",
        )
        if not raw:
            continue
        hist = rebin_for_display(
            raw,
            variable,
            f"display_{era}_{directory_kind}_{process}_{variable}",
            score_display_binning,
        )
        if hist.Integral() <= 0.0:
            continue
        style_background(hist, process)
        process_hists[process] = hist
    if not process_hists:
        print(f"[WARN] {era}/{variable}: no {directory_kind} background histograms found")
        return

    signal_hists = []
    for index, signal_name in enumerate(choose_signals(root_file, era, requested_signals)):
        raw = get_hist(
            root_file,
            f"{signal_prefix}/{signal_name}/{variable}",
            f"{era}_{directory_kind}_{signal_name}_{variable}",
        )
        if not raw:
            continue
        hist = rebin_for_display(
            raw,
            variable,
            f"display_{era}_{directory_kind}_signal_{index}_{variable}",
            score_display_binning,
        )
        if not normalized:
            hist.Scale(signal_scale_factor(signal_name))
        if hist.Integral() <= 0.0:
            continue
        style_signal(hist, index)
        signal_hists.append((signal_name, hist))

    data_hist = None
    if include_data:
        raw_data = get_hist(
            root_file,
            f"{data_prefix}/observed/{variable}",
            f"{era}_{directory_kind}_data_{variable}",
        )
        if raw_data and raw_data.Integral() > 0.0:
            data_hist = rebin_for_display(
                raw_data,
                variable,
                f"display_{era}_{directory_kind}_data_{variable}",
                score_display_binning,
            )
            style_data(data_hist)
        else:
            print(f"[WARN] {era}/{variable}: no {directory_kind} observed data histogram found")

    if include_ratio and data_hist is None:
        print(f"[WARN] {era}/{variable}: cannot draw Data/Sim ratio without data")
        return

    draw_modes = ([True, False] if also_linear else [True])
    for log_y in draw_modes:
        scale_tag = "normalized" if normalized else "raw"
        data_tag = "_wdata_ratio" if include_ratio else ("_with_data" if include_data else "")
        tag = f"{'log' if log_y else 'linear'}_{scale_tag}{data_tag}"
        canvas = ROOT.TCanvas(f"c_{era}_{variable}_{tag}", "", 800, 900 if include_ratio else 800)
        canvas.cd()
        if include_ratio:
            pad_top = ROOT.TPad(f"pad_top_{era}_{variable}_{tag}", "", 0.0, 0.30, 1.0, 1.0)
            pad_bottom = ROOT.TPad(f"pad_bottom_{era}_{variable}_{tag}", "", 0.0, 0.0, 1.0, 0.30)
            pad_top.SetLeftMargin(0.16)
            pad_top.SetRightMargin(0.05)
            pad_top.SetTopMargin(0.10)
            pad_top.SetBottomMargin(0.02)
            pad_top.SetTicks(1, 1)
            pad_top.SetLogy(log_y)
            pad_bottom.SetLeftMargin(0.16)
            pad_bottom.SetRightMargin(0.05)
            pad_bottom.SetTopMargin(0.04)
            pad_bottom.SetBottomMargin(0.38)
            pad_bottom.SetTicks(1, 1)
            pad_bottom.SetGridy(True)
            pad_top.Draw()
            pad_bottom.Draw()
            pad_top.cd()
        else:
            ROOT.gPad.SetLeftMargin(0.16)
            ROOT.gPad.SetRightMargin(0.05)
            ROOT.gPad.SetTopMargin(0.08)
            ROOT.gPad.SetBottomMargin(0.12)
            ROOT.gPad.SetTicks(1, 1)
            ROOT.gPad.SetLogy(log_y)

        stack = ROOT.THStack(f"stack_{era}_{variable}_{tag}", "")
        total = None
        for process in STACK_ORDER:
            hist = process_hists.get(process)
            if not hist:
                continue
            stack.Add(hist)
            if total is None:
                total = clone_hist(hist, f"total_{era}_{variable}_{tag}")
            else:
                total.Add(hist)

        axis = clone_hist(total, f"axis_{era}_{variable}_{tag}")
        axis.Reset("ICESM")
        axis.SetStats(0)
        axis.SetTitle("")
        axis.GetXaxis().SetTitle(VARIABLES[variable]["x_title"])
        axis.GetYaxis().SetTitle(y_title)
        axis.GetXaxis().SetTitleSize(0.0 if include_ratio else 0.05)
        axis.GetXaxis().SetLabelSize(0.0 if include_ratio else 0.04)
        axis.GetYaxis().SetTitleSize(0.055 if include_ratio else 0.05)
        axis.GetYaxis().SetLabelSize(0.045 if include_ratio else 0.04)
        axis.GetYaxis().SetLabelOffset(0.01)
        axis.GetYaxis().SetTitleOffset(1.18 if include_ratio else 1.30)
        axis.GetXaxis().SetRangeUser(VARIABLES[variable]["xmin"], VARIABLES[variable]["xmax"])

        maxima = [total.GetMaximum()] + [hist.GetMaximum() for _, hist in signal_hists] + [1.0]
        if data_hist:
            maxima.append(data_hist.GetMaximum())
        axis.SetMinimum(1.0 if log_y and not normalized else 1.0e-2)
        y_max = max(maxima) * (100.0 if log_y else 2.6)
        if log_y and not normalized:
            y_max = max(y_max, 1.0e9)
        axis.SetMaximum(y_max)
        axis.Draw("hist")
        stack.Draw("hist same")
        total_unc = None
        if draw_uncertainty_band:
            total_unc = clone_uncertainty_band(total, f"total_unc_{era}_{variable}_{tag}")
            total_unc.Draw("E2 same")
        for _, hist in signal_hists:
            hist.Draw("hist same")
        if data_hist:
            data_hist.Draw("E1 same")
        wp_marker = draw_working_point_marker(variable, axis, log_y)

        if include_data:
            if log_y and not normalized:
                legend_y1 = 0.600 if include_ratio else 0.620
            else:
                legend_y1 = 0.600 if include_ratio else 0.645
            legend_text_size = 0.039 if include_ratio else 0.038
            legend_entry_separation = 0.050 if log_y and not normalized else (0.052 if include_ratio else 0.050)
            legend_y2 = 0.855 if log_y and not normalized else (0.875 if include_ratio else 0.910)
        else:
            legend_y1 = 0.730
            legend_text_size = 0.041
            legend_entry_separation = 0.034
            legend_y2 = 0.895
        bkg_entries = [(process_hists[p], PROC_LABEL[p], "F") for p in LEGEND_ORDER if p in process_hists]
        if data_hist:
            bkg_entries.insert(0, (data_hist, "Data", "EP"))
        if total_unc:
            bkg_entries.append((total_unc, "Sim. unc.", "F"))
        sig_entries = [(hist, scaled_signal_label(name, normalized), "L") for name, hist in signal_hists]
        if signal_hists:
            sig_entries.append(("", "m_{dark} = 20 GeV, #kern[-0.18]{#lambda} = 1", ""))

        def configure_legend(legend, margin: float = 0.25):
            legend.SetNColumns(1)
            legend.SetBorderSize(0)
            legend.SetFillColor(ROOT.kWhite)
            legend.SetFillStyle(1001)
            legend.SetTextSize(legend_text_size)
            legend.SetTextFont(42)
            legend.SetEntrySeparation(legend_entry_separation)
            legend.SetMargin(margin)

        if signal_hists:
            bkg_legend = ROOT.TLegend(0.21, legend_y1, 0.41, legend_y2)
            sig_legend = ROOT.TLegend(0.45, legend_y1, 0.98, legend_y2)
            configure_legend(bkg_legend)
            configure_legend(sig_legend, margin=0.09)
            for entry in bkg_entries:
                bkg_legend.AddEntry(*entry)
            for entry in sig_entries:
                sig_legend.AddEntry(*entry)
            bkg_legend.Draw()
            sig_legend.Draw()
        else:
            legend = ROOT.TLegend(0.21 if include_data else 0.23, legend_y1, 0.98 if include_data else 0.97, legend_y2)
            configure_legend(legend)
            for entry in bkg_entries:
                legend.AddEntry(*entry)
            legend.Draw()
        draw_cms_label(
            simulation=not include_data,
            era=era,
            show_lumi=True,
            preliminary=preliminary_label,
        )
        ROOT.gPad.RedrawAxis()

        if include_ratio:
            pad_bottom.cd()
            ratio = clone_hist(data_hist, f"ratio_{era}_{variable}_{tag}")
            ratio.Divide(total)
            ratio.SetStats(0)
            ratio.SetTitle("")
            ratio.SetMinimum(0.75)
            ratio.SetMaximum(1.25)
            ratio.GetXaxis().SetTitle(VARIABLES[variable]["x_title"])
            ratio.GetYaxis().SetTitle("")
            ratio.GetYaxis().CenterTitle(False)
            ratio.GetYaxis().SetNdivisions(505)
            ratio.GetXaxis().SetTitleSize(0.13)
            ratio.GetXaxis().SetLabelSize(0.11)
            ratio.GetXaxis().SetTitleOffset(1.05)
            ratio.GetYaxis().SetLabelSize(0.10)
            ratio.GetYaxis().SetTitleOffset(0.40)
            ratio.GetYaxis().SetTitleSize(0.12)
            ratio.GetYaxis().SetLabelOffset(0.01)
            ratio.GetXaxis().SetRangeUser(VARIABLES[variable]["xmin"], VARIABLES[variable]["xmax"])
            ratio.Draw("PE")

            if draw_uncertainty_band:
                ratio_unc = clone_hist(total, f"ratio_unc_{era}_{variable}_{tag}")
                for ibin in range(1, ratio_unc.GetNbinsX() + 1):
                    mc = total.GetBinContent(ibin)
                    err = total.GetBinError(ibin)
                    ratio_unc.SetBinContent(ibin, 1.0 if mc > 0.0 else 0.0)
                    ratio_unc.SetBinError(ibin, err / mc if mc > 0.0 else 0.0)
                style_uncertainty_band(ratio_unc)
                ratio_unc.Draw("E2 same")
                ratio.Draw("PE same")

            unity = ROOT.TLine(VARIABLES[variable]["xmin"], 1.0, VARIABLES[variable]["xmax"], 1.0)
            unity.SetLineStyle(ROOT.kDashed)
            unity.Draw("same")
            draw_ratio_y_title()
            ROOT.gPad.RedrawAxis()

        output_kind = "raw_wdata_ratio" if (not normalized and include_data and include_ratio) else directory_kind
        out_prefix = os.path.join(plot_dir, era, output_kind, variable)
        os.makedirs(os.path.dirname(out_prefix), exist_ok=True)
        canvas.SaveAs(f"{out_prefix}_{tag}.pdf")
        canvas.SaveAs(f"{out_prefix}_{tag}.png")
        canvas.Close()
        print(f"[OK] Wrote {out_prefix}_{tag}.pdf/.png")


def plot_cache(
    output_path: str,
    plot_dir: str,
    eras: Sequence[str],
    signals: Sequence[str],
    also_linear: bool,
    draw_raw: bool,
    include_data: bool,
    raw_only: bool = False,
    score_display_binning: str = "equal20",
    draw_uncertainty_band: bool = True,
    preliminary_label: bool = False,
) -> None:
    root_file = ROOT.TFile.Open(output_path, "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Could not open cache ROOT file: {output_path}")
    try:
        for era in eras:
            for variable in VARIABLES:
                if not raw_only:
                    draw_plot(root_file, era, variable, signals, plot_dir, normalized=True, also_linear=also_linear, score_display_binning=score_display_binning, draw_uncertainty_band=draw_uncertainty_band, preliminary_label=preliminary_label)
                    if include_data:
                        draw_plot(root_file, era, variable, signals, plot_dir, normalized=True, also_linear=also_linear, include_data=True, score_display_binning=score_display_binning, draw_uncertainty_band=draw_uncertainty_band, preliminary_label=preliminary_label)
                        draw_plot(root_file, era, variable, signals, plot_dir, normalized=True, also_linear=also_linear, include_data=True, include_ratio=True, score_display_binning=score_display_binning, draw_uncertainty_band=draw_uncertainty_band, preliminary_label=preliminary_label)
                if draw_raw or raw_only:
                    if not raw_only:
                        draw_plot(root_file, era, variable, signals, plot_dir, normalized=False, also_linear=also_linear, score_display_binning=score_display_binning, draw_uncertainty_band=draw_uncertainty_band, preliminary_label=preliminary_label)
                        if include_data:
                            draw_plot(root_file, era, variable, signals, plot_dir, normalized=False, also_linear=also_linear, include_data=True, score_display_binning=score_display_binning, draw_uncertainty_band=draw_uncertainty_band, preliminary_label=preliminary_label)
                    if include_data:
                        draw_plot(root_file, era, variable, signals, plot_dir, normalized=False, also_linear=also_linear, include_data=True, include_ratio=True, score_display_binning=score_display_binning, draw_uncertainty_band=draw_uncertainty_band, preliminary_label=preliminary_label)
    finally:
        root_file.Close()


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--background-base", default=DEFAULT_BACKGROUND_BASE)
    parser.add_argument("--data-base", default=DEFAULT_DATA_BASE)
    parser.add_argument("--signal-base", default=DEFAULT_SIGNAL_BASE, help="optional base containing YEAR/nominal/ABCD_scores*.root")
    parser.add_argument(
        "--signal-weight-variation",
        choices=["lundWeightNom", "nominal", "auto"],
        default="lundWeightNom",
        help=(
            "Signal score-file suffix to use. The lundWeightNom files have the "
            "nominal Lund-corrected weights already written into their weights branch."
        ),
    )
    parser.add_argument("--signals", default=",".join(DEFAULT_SIGNAL_DIRS), help="comma-separated signal name tokens, or all")
    parser.add_argument("--no-signals", action="store_true", help="skip signal score files and overlays")
    parser.add_argument("--data-samples", default=",".join(DEFAULT_DATA_SAMPLES), help="comma-separated data score samples")
    parser.add_argument("--years", default="2016,2017,2018")
    parser.add_argument("--output", default="EventDNNScore_plots/dnn_scores.root")
    parser.add_argument("--plot-dir", default="EventDNNScore_plots")
    parser.add_argument("--chunk-size", default="100 MB")
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--strict", action="store_true")
    parser.add_argument("--strict-positive-weights", action="store_true", help="drop zero/negative weights while filling")
    parser.add_argument("--no-data", action="store_true", help="skip data score files and data overlay/ratio plots")
    parser.add_argument("--also-linear", action="store_true")
    parser.add_argument("--draw-raw", action="store_true", help="also draw unnormalized weighted-yield stacks")
    parser.add_argument("--raw-only", action="store_true", help="only draw raw data/MC-ratio plots")
    parser.add_argument(
        "--no-uncertainty-band",
        action="store_true",
        help="Do not draw the simulated-background uncertainty band or its legend entry",
    )
    parser.add_argument("--preliminary-label", action="store_true", help="draw CMS Preliminary or CMS Simulation Preliminary")
    parser.add_argument("--eras", default="", help="plot-only eras; default is years plus Run2")
    parser.add_argument(
        "--score-display-binning",
        choices=["equal20", "wpedge10", "wpedge20"],
        default="equal20",
        help=(
            "DNN score display binning. equal20 keeps 20 equal bins; "
            "wpedge10 uses 10 bins with 0.85 as an edge; wpedge20 uses "
            "20 bins with 0.85 as the final bin edge before the high-score tail."
        ),
    )
    parser.add_argument(
        "--score-fill-binning",
        choices=["equal20", "wpedge20"],
        default="equal20",
        help="DNN score cache binning used while reading score trees.",
    )
    return parser.parse_args()


def main() -> int:
    global SCORE_FILL_BINNING
    args = parse_args()
    SCORE_FILL_BINNING = args.score_fill_binning
    years = parse_csv(args.years)
    signals = [] if args.no_signals else parse_csv(args.signals)
    data_samples = [] if args.no_data else parse_csv(args.data_samples)
    eras = parse_csv(args.eras) if args.eras else list(years) + ["Run2"]

    if args.plot_only:
        plot_cache(
            args.output,
            args.plot_dir,
            eras,
            signals,
            args.also_linear,
            args.draw_raw,
            include_data=not args.no_data,
            raw_only=args.raw_only,
            score_display_binning=args.score_display_binning,
            draw_uncertainty_band=not args.no_uncertainty_band,
            preliminary_label=args.preliminary_label,
        )
        return 0

    group_hists, run2_hists, background_failures = read_backgrounds(
        args.background_base,
        years,
        args.chunk_size,
        args.strict,
        args.strict_positive_weights,
    )
    data_hists, run2_data_hists, data_failures = read_data(
        args.data_base,
        years,
        data_samples,
        args.chunk_size,
        args.strict,
        args.strict_positive_weights,
    )
    signal_hists, run2_signal_hists = read_signals(
        args.signal_base,
        years,
        signals,
        args.chunk_size,
        args.strict_positive_weights,
        args.signal_weight_variation,
    )
    write_cache(
        args.output,
        group_hists,
        run2_hists,
        signal_hists,
        run2_signal_hists,
        data_hists,
        run2_data_hists,
        years,
    )
    plot_cache(
        args.output,
        args.plot_dir,
        eras,
        signals,
        args.also_linear,
        args.draw_raw,
        include_data=bool(data_hists) and not args.no_data,
        raw_only=args.raw_only,
        score_display_binning="equal20",
        draw_uncertainty_band=not args.no_uncertainty_band,
        preliminary_label=args.preliminary_label,
    )

    failures = background_failures + data_failures
    if failures:
        print(f"\n[WARN] Completed with {len(failures)} file failure(s). First few:")
        for failure in failures[:10]:
            print(f"  - {failure}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[FATAL] {exc}", file=sys.stderr)
        raise
