#!/usr/bin/env python
"""Make cut-flow tables, plots, and a ROOT cache from skim CutFlow trees.

Run from the analysis area with the usual environment:

  source condor/initCondor.sh
  python makeCutFlowTables_fromSkims.py --years 2016,2017,2018

The skim files contain a one-entry TTree named ``CutFlow``. The values are
cumulative weighted cross sections in pb. This script sums those counters over
all files in each sample group and multiplies by the per-year luminosity to get
expected event yields. For multiple years, the table values are computed from
the summed expected-event counters.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.TH1.SetDefaultSumw2(True)


EOS_HOST = "root://cmseos.fnal.gov"
LUMI_PB = {
    "2016": 36.31 * 1000.0,
    "2017": 42.07 * 1000.0,
    "2018": 59.56 * 1000.0,
}

DEFAULT_BACKGROUND_BASE = (
    "/store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/data_mc_noWNAE"
)
DEFAULT_SIGNAL_BASE = "/store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto"

BACKGROUND_MANIFEST = {
    "QCD": [
        "QCD_Pt_1000to1400",
        "QCD_Pt_1400to1800",
        "QCD_Pt_170to300",
        "QCD_Pt_1800to2400",
        "QCD_Pt_2400to3200",
        "QCD_Pt_300to470",
        "QCD_Pt_3200toInf",
        "QCD_Pt_470to600",
        "QCD_Pt_600to800",
        "QCD_Pt_800to1000",
    ],
    "TTJets": [
        "TTJets",
        "TTJets_DiLept",
        "TTJets_DiLept_genMET-150",
        "TTJets_HT-1200to2500",
        "TTJets_HT-2500toInf",
        "TTJets_HT-600to800",
        "TTJets_HT-800to1200",
        "TTJets_SingleLeptFromT",
        "TTJets_SingleLeptFromT_genMET-150",
        "TTJets_SingleLeptFromTbar",
        "TTJets_SingleLeptFromTbar_genMET-150",
    ],
    "WJetsToLNu": [
        "WJetsToLNu_HT-1200To2500",
        "WJetsToLNu_HT-2500ToInf",
        "WJetsToLNu_HT-400To600",
        "WJetsToLNu_HT-600To800",
        "WJetsToLNu_HT-800To1200",
    ],
    "ZJetsToNuNu": [
        "ZJetsToNuNu_HT-1200To2500",
        "ZJetsToNuNu_HT-2500ToInf",
        "ZJetsToNuNu_HT-400To600",
        "ZJetsToNuNu_HT-600To800",
        "ZJetsToNuNu_HT-800To1200",
    ],
}

BACKGROUND_ORDER = ["QCD", "TTJets", "WJetsToLNu", "ZJetsToNuNu"]
BACKGROUND_TABLE_LABELS = {
    "QCD": "QCD",
    "TTJets": "ttbar",
    "WJetsToLNu": "W(->lv)+jets",
    "ZJetsToNuNu": "Z(->vv)+jets",
}

DEFAULT_SIGNAL_NAMES = [
    "t-channel_mMed-1000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
    "t-channel_mMed-2000_mDark-20_rinv-0p1_alpha-peak_yukawa-1",
    "t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
    "t-channel_mMed-2000_mDark-20_rinv-0p7_alpha-peak_yukawa-1",
    "t-channel_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
]

# The branch order is the real cumulative order in the skim CutFlow tree.
# By default the displayed table stops at PhiSpikeFilter to match the paper
# preselection table; the ROOT cache also stores GapJetVeto and Final.
ALL_CUTFLOW_STEPS = [
    ("Initial", "Initial"),
    ("Trigger", "Trigger"),
    ("STGt1300GeV", "ST > 1300 GeV"),
    ("GoodJetsAK8", "All good AK8 jets"),
    ("nJetsAK8Gt2", "n AK8 jets >= 2"),
    ("LeptonVeto", "Lepton Veto"),
    ("DeltaPhiMinLt1p5", "DeltaPhi_min(J,MET) < 1.5"),
    ("METGt200GeV", "MET > 200 GeV"),
    ("METFilters", "MET filters"),
    ("PhiSpikeFilter", "phi spike filter"),
    ("GapJetVeto", "Gap jet veto"),
    ("Final", "Final"),
]

DEFAULT_LAST_TABLE_BRANCH = "PhiSpikeFilter"


@dataclass
class CutFlowResult:
    name: str
    label: str
    files: List[str]
    sample_names: List[str]
    years: List[str]
    weighted_pb: Dict[str, float]
    expected_events_by_branch: Dict[str, float]
    raw_final_entries: int

    def value(self, branch: str) -> float:
        return float(self.expected_events_by_branch.get(branch, 0.0))

    def expected_events(self, branch: str) -> float:
        return self.value(branch)

    def pb_value(self, branch: str) -> float:
        return float(self.weighted_pb.get(branch, 0.0))


def normalize_eos_path(path: str) -> str:
    path = os.path.expanduser(path).rstrip("/")
    for prefix in ("/eos/uscms", "/eos/cms"):
        if path.startswith(prefix + "/"):
            return path[len(prefix) :]
    return path


def eos_url(path: str) -> str:
    return EOS_HOST + "//" + normalize_eos_path(path).lstrip("/")


def xrdfs_ls(path: str, recursive: bool = False) -> List[str]:
    command = ["xrdfs", EOS_HOST, "ls"]
    if recursive:
        command.append("-R")
    command.append(normalize_eos_path(path))
    proc = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"xrdfs failed for {path}\ncommand: {' '.join(command)}\n"
            f"stderr: {proc.stderr.strip()}"
        )
    return [line.strip() for line in proc.stdout.splitlines() if line.strip()]


def get_part_files(sample_dir: str, max_files: int = 0) -> List[str]:
    files = sorted(
        eos_url(path)
        for path in xrdfs_ls(sample_dir, recursive=True)
        if re.search(r"/part-[0-9]+\.root$", path)
    )
    return files[:max_files] if max_files > 0 else files


def selected_steps(include_gap_jet_veto: bool) -> List[Tuple[str, str]]:
    if include_gap_jet_veto:
        return [(branch, label) for branch, label in ALL_CUTFLOW_STEPS if branch != "Final"]

    out: List[Tuple[str, str]] = []
    for branch, label in ALL_CUTFLOW_STEPS:
        if branch == "Initial":
            continue
        out.append((branch, label))
        if branch == DEFAULT_LAST_TABLE_BRANCH:
            break
    return out


def root_steps() -> List[Tuple[str, str]]:
    return [(branch, label) for branch, label in ALL_CUTFLOW_STEPS if branch != "Final"]


def read_one_cutflow(file_url: str) -> Tuple[Dict[str, float], int]:
    with uproot.open(file_url) as root_file:
        if "CutFlow" not in root_file:
            raise KeyError(f"No CutFlow tree in {file_url}")

        cutflow = root_file["CutFlow"]
        arrays = cutflow.arrays(cutflow.keys(), library="np")
        values = {
            branch: float(np.asarray(arrays[branch])[0])
            for branch in cutflow.keys()
            if len(np.asarray(arrays[branch])) > 0
        }

        final_entries = 0
        if "Events" in root_file:
            final_entries = int(root_file["Events"].num_entries)
        return values, final_entries


def sum_cutflow(
    files: Sequence[str],
    name: str,
    label: str,
    sample_names: Sequence[str],
    year: str,
) -> CutFlowResult:
    totals_pb: Dict[str, float] = {}
    totals_events: Dict[str, float] = {}
    raw_final_entries = 0
    lumi_pb = LUMI_PB[year]

    for index, file_url in enumerate(files, start=1):
        try:
            values, entries = read_one_cutflow(file_url)
        except Exception as exc:
            print(f"[WARN] {name}: failed to read {file_url}: {exc}", file=sys.stderr)
            continue

        raw_final_entries += entries
        for branch, value in values.items():
            totals_pb[branch] = totals_pb.get(branch, 0.0) + value
            totals_events[branch] = totals_events.get(branch, 0.0) + value * lumi_pb

        if index == 1 or index == len(files) or index % 25 == 0:
            print(f"  {name}: {index}/{len(files)} files")

    return CutFlowResult(
        name=name,
        label=label,
        files=list(files),
        sample_names=list(sample_names),
        years=[year],
        weighted_pb=totals_pb,
        expected_events_by_branch=totals_events,
        raw_final_entries=raw_final_entries,
    )


def add_results(target: CutFlowResult, source: CutFlowResult) -> None:
    target.files.extend(source.files)
    target.sample_names.extend(source.sample_names)
    target.years.extend(source.years)
    target.raw_final_entries += source.raw_final_entries
    for branch, value in source.weighted_pb.items():
        target.weighted_pb[branch] = target.weighted_pb.get(branch, 0.0) + value
    for branch, value in source.expected_events_by_branch.items():
        target.expected_events_by_branch[branch] = (
            target.expected_events_by_branch.get(branch, 0.0) + value
        )


def aggregate_year_results(results_by_year: Sequence[CutFlowResult]) -> List[CutFlowResult]:
    aggregated: Dict[str, CutFlowResult] = {}
    order: List[str] = []
    for result in results_by_year:
        if result.name not in aggregated:
            aggregated[result.name] = CutFlowResult(
                name=result.name,
                label=result.label,
                files=[],
                sample_names=[],
                years=[],
                weighted_pb={},
                expected_events_by_branch={},
                raw_final_entries=0,
            )
            order.append(result.name)
        add_results(aggregated[result.name], result)
    return [aggregated[name] for name in order]


def discover_backgrounds(
    background_base: str,
    years: Sequence[str],
    max_files_per_sample: int,
) -> List[CutFlowResult]:
    all_results: List[CutFlowResult] = []
    for year in years:
        print(f"[INFO] Reading {year} background CutFlow trees")
        all_results.extend(
            discover_backgrounds_one_year(background_base, year, max_files_per_sample)
        )
    return aggregate_year_results(all_results)


def discover_backgrounds_one_year(
    background_base: str,
    year: str,
    max_files_per_sample: int,
) -> List[CutFlowResult]:
    nominal = (
        f"{normalize_eos_path(background_base)}/{year}/"
        "t_channel_pre_selection/nominal"
    )
    available = {
        os.path.basename(path.rstrip("/")): path
        for path in xrdfs_ls(nominal, recursive=False)
    }

    results: List[CutFlowResult] = []
    for process in BACKGROUND_ORDER:
        files: List[str] = []
        used_samples: List[str] = []
        for sample in BACKGROUND_MANIFEST[process]:
            if sample not in available:
                print(f"[WARN] {year}: missing {process} sample: {sample}", file=sys.stderr)
                continue
            sample_files = get_part_files(available[sample], max_files_per_sample)
            if not sample_files:
                print(f"[WARN] no ROOT files for {sample}", file=sys.stderr)
                continue
            files.extend(sample_files)
            used_samples.append(sample)

        print(f"[INFO] {year} {process}: {len(used_samples)} samples, {len(files)} files")
        results.append(
            sum_cutflow(
                files,
                name=process,
                label=BACKGROUND_TABLE_LABELS[process],
                sample_names=used_samples,
                year=year,
            )
        )
    return results


def discover_signals(
    signal_base: str,
    years: Sequence[str],
    signal_names: Sequence[str],
    max_files_per_sample: int,
) -> List[CutFlowResult]:
    all_results: List[CutFlowResult] = []
    for year in years:
        print(f"[INFO] Reading {year} signal CutFlow trees")
        all_results.extend(
            discover_signals_one_year(signal_base, year, signal_names, max_files_per_sample)
        )
    return aggregate_year_results(all_results)


def discover_signals_one_year(
    signal_base: str,
    year: str,
    signal_names: Sequence[str],
    max_files_per_sample: int,
) -> List[CutFlowResult]:
    nominal = (
        f"{normalize_eos_path(signal_base)}/{year}/"
        "t_channel_pre_selection/nominal"
    )
    available = {
        os.path.basename(path.rstrip("/")): path
        for path in xrdfs_ls(nominal, recursive=False)
    }

    results: List[CutFlowResult] = []
    for signal_name in signal_names:
        if signal_name not in available:
            print(f"[WARN] missing signal: {signal_name}", file=sys.stderr)
            continue

        files = get_part_files(available[signal_name], max_files_per_sample)
        print(f"[INFO] {signal_name}: {len(files)} files")
        results.append(
            sum_cutflow(
                files,
                name=signal_name,
                label=signal_label(signal_name),
                sample_names=[signal_name],
                year=year,
            )
        )
    return results


def token_to_text(token: str) -> str:
    return token.replace("p", ".")


def signal_label(signal_name: str) -> str:
    mmed = re.search(r"mMed-([0-9]+(?:p[0-9]+)?)", signal_name)
    rinv = re.search(r"rinv-([0-9]+(?:p[0-9]+)?)", signal_name)
    if not mmed or not rinv:
        return signal_name
    return f"mMed {token_to_text(mmed.group(1))}, rinv {token_to_text(rinv.group(1))}"


def signal_header_rows(signals: Sequence[CutFlowResult]) -> Tuple[List[str], List[str]]:
    mmeds: List[str] = []
    rinvs: List[str] = []
    for result in signals:
        mmed = re.search(r"mMed-([0-9]+(?:p[0-9]+)?)", result.name)
        rinv = re.search(r"rinv-([0-9]+(?:p[0-9]+)?)", result.name)
        mmeds.append(token_to_text(mmed.group(1)) if mmed else result.label)
        rinvs.append(token_to_text(rinv.group(1)) if rinv else "")
    return mmeds, rinvs


def relative_efficiency(result: CutFlowResult, branch: str) -> float:
    branches = [item[0] for item in ALL_CUTFLOW_STEPS]
    if branch not in branches or branch == "Initial":
        return float("nan")
    previous = branches[branches.index(branch) - 1]
    denominator = result.value(previous)
    if denominator == 0.0:
        return float("nan")
    return 100.0 * result.value(branch) / denominator


def percent_of_initial(result: CutFlowResult, branch: str) -> float:
    initial = result.value("Initial")
    if initial == 0.0:
        return float("nan")
    return 100.0 * result.value(branch) / initial


def fmt_percent(value: float) -> str:
    if not np.isfinite(value):
        return "nan"
    if abs(value - round(value)) < 0.05 and abs(value) >= 99.95:
        return f"{value:.0f}"
    if value < 0.01:
        return f"{value:.4g}"
    if value < 10.0:
        return f"{value:.2g}"
    return f"{value:.1f}"


def fmt_yield(value: float) -> str:
    if not np.isfinite(value):
        return "nan"
    rounded = int(round(value))
    if abs(rounded) >= 100000:
        return f"{rounded / 1000.0:.0f}k"
    return f"{rounded:d}"


def table_matrix(
    results: Sequence[CutFlowResult],
    steps: Sequence[Tuple[str, str]],
    era_label: str,
) -> List[List[str]]:
    rows: List[List[str]] = []
    for branch, label in steps:
        rows.append([label] + [fmt_percent(relative_efficiency(result, branch)) for result in results])

    last_branch = steps[-1][0]
    rows.append(
        ["Total absolute efficiency [%]"]
        + [fmt_percent(percent_of_initial(result, last_branch)) for result in results]
    )
    rows.append(
        [f"Expected yield ({era_label})"]
        + [fmt_yield(result.expected_events(last_branch)) for result in results]
    )
    return rows


def write_csv_table(path: str, headers: Sequence[str], rows: Sequence[Sequence[str]]) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        writer.writerows(rows)


def write_txt_table(path: str, headers: Sequence[str], rows: Sequence[Sequence[str]]) -> None:
    widths = [len(item) for item in headers]
    for row in rows:
        widths = [max(width, len(item)) for width, item in zip(widths, row)]

    def line(items: Sequence[str]) -> str:
        return "  ".join(item.rjust(width) if idx else item.ljust(width) for idx, (item, width) in enumerate(zip(items, widths)))

    with open(path, "w") as handle:
        handle.write(line(headers) + "\n")
        handle.write(line(["-" * width for width in widths]) + "\n")
        for row in rows:
            handle.write(line(row) + "\n")


def write_tex_table(path: str, headers: Sequence[str], rows: Sequence[Sequence[str]]) -> None:
    with open(path, "w") as handle:
        handle.write("\\begin{tabular}{l" + "c" * (len(headers) - 1) + "}\n")
        handle.write("\\hline\n")
        handle.write(" & ".join(headers) + " \\\\\n")
        handle.write("\\hline\n")
        for row in rows:
            handle.write(" & ".join(row) + " \\\\\n")
        handle.write("\\hline\n")
        handle.write("\\end{tabular}\n")


def write_detailed_csv(
    path: str,
    category: str,
    results: Sequence[CutFlowResult],
    steps: Sequence[Tuple[str, str]],
) -> None:
    fields = [
        "category",
        "name",
        "label",
        "branch",
        "cut_label",
        "weighted_pb",
        "expected_events",
        "relative_efficiency_percent",
        "percent_of_initial",
        "raw_final_skim_entries",
        "n_files",
        "years",
        "samples",
    ]
    with open(path, "a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        if handle.tell() == 0:
            writer.writeheader()
        for result in results:
            for branch, label in steps:
                writer.writerow(
                    {
                        "category": category,
                        "name": result.name,
                        "label": result.label,
                        "branch": branch,
                        "cut_label": label,
                        "weighted_pb": f"{result.pb_value(branch):.12g}",
                        "expected_events": f"{result.expected_events(branch):.12g}",
                        "relative_efficiency_percent": f"{relative_efficiency(result, branch):.12g}",
                        "percent_of_initial": f"{percent_of_initial(result, branch):.12g}",
                        "raw_final_skim_entries": result.raw_final_entries,
                        "n_files": len(result.files),
                        "years": ",".join(sorted(set(result.years))),
                        "samples": ",".join(result.sample_names),
                    }
                )


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text)


def draw_cutflow_plot(
    path_base: str,
    title: str,
    results: Sequence[CutFlowResult],
    steps: Sequence[Tuple[str, str]],
    era_label: str,
) -> None:
    labels = [label for _, label in steps]
    x = np.arange(len(steps))

    fig, (ax_yield, ax_pct) = plt.subplots(
        2,
        1,
        figsize=(max(9.5, 0.75 * len(steps) + 3.0), 8.0),
        sharex=True,
        gridspec_kw={"height_ratios": [1.2, 1.0]},
    )

    for result in results:
        yields = np.array([result.expected_events(branch) for branch, _ in steps], dtype=float)
        pct = np.array([percent_of_initial(result, branch) for branch, _ in steps], dtype=float)
        ax_yield.plot(x, yields, marker="o", linewidth=1.8, label=result.label)
        ax_pct.plot(x, pct, marker="o", linewidth=1.8, label=result.label)

    ax_yield.set_title(title)
    ax_yield.set_ylabel(f"Expected events, {era_label}")
    ax_yield.set_yscale("log")
    ax_yield.grid(True, which="both", alpha=0.25)
    ax_yield.legend(fontsize=8, ncol=2)

    ax_pct.set_ylabel("Percent of initial")
    ax_pct.set_ylim(bottom=0.0)
    ax_pct.grid(True, alpha=0.25)
    ax_pct.set_xticks(x)
    ax_pct.set_xticklabels(labels, rotation=35, ha="right")

    fig.tight_layout()
    fig.savefig(path_base + ".png", dpi=200)
    fig.savefig(path_base + ".pdf")
    plt.close(fig)


def make_hist(name: str, title: str, values: Sequence[float], steps: Sequence[Tuple[str, str]]) -> ROOT.TH1D:
    hist = ROOT.TH1D(name, title, len(steps), 0.5, len(steps) + 0.5)
    hist.SetDirectory(0)
    for idx, ((_, label), value) in enumerate(zip(steps, values), start=1):
        hist.GetXaxis().SetBinLabel(idx, label)
        hist.SetBinContent(idx, float(value))
        hist.SetBinError(idx, 0.0)
    return hist


def mkdirs(root_file: ROOT.TFile, path: str) -> ROOT.TDirectory:
    directory: ROOT.TDirectory = root_file
    for part in [item for item in path.split("/") if item]:
        next_dir = directory.GetDirectory(part)
        if not next_dir:
            next_dir = directory.mkdir(part)
        directory = next_dir
    return directory


def write_root_file(
    path: str,
    backgrounds: Sequence[CutFlowResult],
    signals: Sequence[CutFlowResult],
    steps: Sequence[Tuple[str, str]],
    years: Sequence[str],
    era_label: str,
    args: argparse.Namespace,
) -> None:
    root_file = ROOT.TFile.Open(path, "RECREATE")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Could not create ROOT file: {path}")

    try:
        meta = {
            "years": list(years),
            "era_label": era_label,
            "lumi_pb": {year: LUMI_PB[year] for year in years},
            "background_base": args.background_base,
            "signal_base": args.signal_base,
            "steps": [branch for branch, _ in steps],
            "note": "weighted_pb is summed over input files; expected_events is sum_year(lumi_pb_year * weighted_pb_year).",
        }
        ROOT.TNamed("metadata_json", json.dumps(meta, indent=2, sort_keys=True)).Write()

        for category, results in [("backgrounds", backgrounds), ("signals", signals)]:
            for result in results:
                directory = mkdirs(root_file, f"cutflow/{category}/{safe_name(result.name)}")
                directory.cd()
                pb = [result.pb_value(branch) for branch, _ in steps]
                yields = [result.expected_events(branch) for branch, _ in steps]
                rel = [relative_efficiency(result, branch) for branch, _ in steps]
                remain = [percent_of_initial(result, branch) for branch, _ in steps]

                make_hist("weighted_pb", "CutFlow weighted cross section;Cut;pb", pb, steps).Write()
                make_hist("expected_events", "CutFlow expected events;Cut;Events", yields, steps).Write()
                make_hist("relative_efficiency_percent", "Relative cut efficiency;Cut;%", rel, steps).Write()
                make_hist("percent_of_initial", "Percent of initial;Cut;%", remain, steps).Write()
                ROOT.TNamed("input_files", "\n".join(result.files)).Write()
                ROOT.TNamed("sample_names", "\n".join(result.sample_names)).Write()
                ROOT.TNamed("years", "\n".join(sorted(set(result.years)))).Write()

        root_file.Write()
    finally:
        root_file.Close()


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--background-base", default=DEFAULT_BACKGROUND_BASE)
    parser.add_argument("--signal-base", default=DEFAULT_SIGNAL_BASE)
    parser.add_argument("--signals", default=",".join(DEFAULT_SIGNAL_NAMES))
    parser.add_argument("--years", default="2018")
    parser.add_argument("--out-dir", default=None)
    parser.add_argument("--max-files-per-sample", type=int, default=0)
    parser.add_argument(
        "--include-gap-jet-veto",
        action="store_true",
        help="Include GapJetVeto in the printed tables/plots.",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    years = [item.strip() for item in args.years.split(",") if item.strip()]
    invalid_years = [year for year in years if year not in LUMI_PB]
    if invalid_years:
        raise ValueError(
            f"Unsupported year(s): {invalid_years}; choose from {sorted(LUMI_PB)}"
        )
    era_label = "Run2" if years == ["2016", "2017", "2018"] else "+".join(years)
    if args.out_dir is None:
        args.out_dir = (
            "cutflow_Run2_from_skims"
            if era_label == "Run2"
            else f"cutflow_{era_label}_from_skims"
        )
    os.makedirs(args.out_dir, exist_ok=True)

    signal_names = [item.strip() for item in args.signals.split(",") if item.strip()]
    table_steps = selected_steps(args.include_gap_jet_veto)
    cache_steps = root_steps()

    print(f"[INFO] Aggregating years: {', '.join(years)} ({era_label})")
    backgrounds = discover_backgrounds(args.background_base, years, args.max_files_per_sample)
    signals = discover_signals(args.signal_base, years, signal_names, args.max_files_per_sample)

    detailed_path = os.path.join(args.out_dir, "cutflow_detailed.csv")
    if os.path.exists(detailed_path):
        os.remove(detailed_path)
    write_detailed_csv(detailed_path, "background", backgrounds, cache_steps)
    write_detailed_csv(detailed_path, "signal", signals, cache_steps)

    background_headers = ["Background"] + [result.label for result in backgrounds]
    background_rows = table_matrix(backgrounds, table_steps, era_label)
    write_csv_table(os.path.join(args.out_dir, "cutflow_background_table.csv"), background_headers, background_rows)
    write_txt_table(os.path.join(args.out_dir, "cutflow_background_table.txt"), background_headers, background_rows)
    write_tex_table(os.path.join(args.out_dir, "cutflow_background_table.tex"), background_headers, background_rows)

    signal_mmeds, signal_rinvs = signal_header_rows(signals)
    signal_headers = ["Signal"] + [result.label for result in signals]
    signal_rows = (
        [["mMed [GeV]"] + signal_mmeds, ["rinv"] + signal_rinvs]
        + table_matrix(signals, table_steps, era_label)
    )
    write_csv_table(os.path.join(args.out_dir, "cutflow_signal_table.csv"), signal_headers, signal_rows)
    write_txt_table(os.path.join(args.out_dir, "cutflow_signal_table.txt"), signal_headers, signal_rows)
    write_tex_table(os.path.join(args.out_dir, "cutflow_signal_table.tex"), signal_headers, signal_rows)

    if backgrounds:
        draw_cutflow_plot(
            os.path.join(args.out_dir, f"cutflow_background_{era_label}"),
            f"{era_label} background cut flow",
            backgrounds,
            table_steps,
            era_label,
        )
    if signals:
        draw_cutflow_plot(
            os.path.join(args.out_dir, f"cutflow_signal_{era_label}"),
            f"{era_label} signal cut flow",
            signals,
            table_steps,
            era_label,
        )

    write_root_file(
        os.path.join(args.out_dir, f"cutflow_{era_label}.root"),
        backgrounds,
        signals,
        cache_steps,
        years,
        era_label,
        args,
    )

    print(f"[DONE] Wrote outputs under {args.out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
