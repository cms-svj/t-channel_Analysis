#!/usr/bin/env python3
"""Build ParticleNet PF-candidate input-variable stacks with Figure2 styling."""

from __future__ import annotations

import argparse
import os
from typing import Dict, List, Optional, Sequence, Tuple

import awkward as ak
import numpy as np
import uproot

import Figure2_makerusingskims as f2


VARIABLES = {
    "del_eta": {
        "title": "#Delta#eta",
        "x_title": "#Delta#eta",
        "nbins": 60,
        "xmin": -1.2,
        "xmax": 1.2,
        "display_xmin": -1.0,
        "display_xmax": 1.0,
    },
    "del_phi": {
        "title": "#Delta#phi",
        "x_title": "#Delta#phi",
        "nbins": 60,
        "xmin": -1.2,
        "xmax": 1.2,
        "display_xmin": -1.0,
        "display_xmax": 1.0,
    },
    "del_r": {
        "title": "Constituent#kern[0.30]{#Delta}R",
        "x_title": "Constituent#kern[0.30]{#Delta}R",
        "nbins": 60,
        "xmin": 0.0,
        "xmax": 0.8,
    },
    "log_pt": {
        "title": "Constituent log(p_{T})",
        "x_title": "Constituent log(p_{T})",
        "nbins": 60,
        "xmin": -2.0,
        "xmax": 8.0,
        "display_xmin": -1.0,
        "display_xmax": 5.0,
        "display_rebin_factor_override": 1,
    },
    "log_e": {
        "title": "Constituent log(E)",
        "x_title": "Constituent log(E)",
        "nbins": 60,
        "xmin": -2.0,
        "xmax": 8.0,
        "display_xmin": -1.0,
        "display_xmax": 6.0,
        "display_rebin_factor_override": 1,
    },
    "log_pt_jetpt": {
        "title": "Constituent log(p_{T}/p_{T}^{J})",
        "x_title": "Constituent log(p_{T}/p_{T}^{J})",
        "nbins": 60,
        "xmin": -10.0,
        "xmax": 0.0,
        "display_xmin": -8.0,
        "display_xmax": -1.0,
        "display_rebin_factor_override": 1,
    },
    "log_e_jete": {
        "title": "Constituent log(E/E_{J})",
        "x_title": "Constituent log(E/E_{J})",
        "nbins": 60,
        "xmin": -10.0,
        "xmax": 0.0,
        "display_xmin": -8.0,
        "display_xmax": -1.0,
        "display_rebin_factor_override": 1,
    },
}

PNET_INPUT_STYLE = {
    "linear_ymax_factor": 1.35,
    "linear_ymax": 20.0e6,
    "display_rebin_factor": 2,
    "y_exponent_left": True,
    "y_exponent_offset": (-0.085, 0.020),
    "y_max_digits": 3,
    "ratio_ymin": 0.75,
    "ratio_ymax": 1.25,
    "log_ymin": 1.0e3,
    "signal_scale_factors": {
        ("600", "0.3"): 8.0,
        ("2000", "0.1"): 550.0,
        ("2000", "0.3"): 250.0,
    },
}

PARTICLENET_SIGNAL_DIRS = [
    "t-channel_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
    "t-channel_mMed-2000_mDark-20_rinv-0p1_alpha-peak_yukawa-1",
    "t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
    # "t-channel_mMed-2000_mDark-20_rinv-0p7_alpha-peak_yukawa-1",
    # "t-channel_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
]

for _cfg in VARIABLES.values():
    _cfg.update(PNET_INPUT_STYLE)

for _variable in ("del_eta", "del_phi"):
    VARIABLES[_variable]["linear_ymax"] = 30.0e6
VARIABLES["del_r"]["linear_ymax"] = 18.0e6
for _variable in ("log_pt", "log_e", "log_pt_jetpt", "log_e_jete"):
    VARIABLES[_variable]["linear_ymax"] = 18.0e6

# PNET_INPUT_STYLE's blanket "display_rebin_factor": 2 would still coarsen
# these four beyond the "display_rebin_factor_override": 1 set above (the
# two keys stack in clone_rebin), which would round the requested display
# crop out to the next wide bin edge. Neutralize it just for these four so
# the crop lands exactly where requested.
for _variable in ("log_pt", "log_e", "log_pt_jetpt", "log_e_jete"):
    VARIABLES[_variable]["display_rebin_factor"] = 1


def halve_bins_above_threshold(variables: Dict[str, Dict[str, object]], threshold: int = 10) -> None:
    for cfg in variables.values():
        nbins = int(cfg.get("nbins", 0))
        if nbins > threshold:
            cfg["nbins"] = max(1, nbins // 2)


halve_bins_above_threshold(VARIABLES)

REQUIRED_BRANCHES = [
    "JetsAK8_constituentsIndex",
    "JetsAK8_constituentsIndexCounts",
    "JetsAK8_isGood",
    "JetsAK8_/.fPt",
    "JetsAK8_/.fEta",
    "JetsAK8_/.fPhi",
    "JetsAK8_/.fE",
    "JetsConstituents_/.fPt",
    "JetsConstituents_/.fEta",
    "JetsConstituents_/.fPhi",
    "JetsConstituents_/.fE",
]


def configure_figure2_globals() -> None:
    f2.VARIABLES = VARIABLES


def available_branches(tree: uproot.behaviors.TTree.TTree, is_data: bool = False) -> List[str]:
    available = set(tree.keys())
    required = list(REQUIRED_BRANCHES)
    if not is_data:
        required.append("Weight")
        for branch in ("puWeight", "NonPrefiringProb", "lundWeightNom"):
            if branch in available:
                required.append(branch)

    missing = [branch for branch in required if branch not in available]
    if missing:
        raise KeyError(f"Events tree lacks required branch(es): {missing}")
    return required


def delta_phi(phi_a, phi_b):
    dphi = phi_a - phi_b
    dphi = ak.where(dphi < -np.pi, dphi + 2.0 * np.pi, dphi)
    dphi = ak.where(dphi > np.pi, dphi - 2.0 * np.pi, dphi)
    return dphi


def pnet_candidate_values(arrays, event_weights: np.ndarray) -> Tuple[Dict[str, np.ndarray], np.ndarray]:
    flat_index = ak.flatten(arrays["JetsAK8_constituentsIndex"], axis=1)
    flat_counts = ak.flatten(arrays["JetsAK8_constituentsIndexCounts"], axis=None)
    index_by_jet_flat = ak.unflatten(flat_index, flat_counts, axis=0)
    index_by_event = ak.unflatten(
        index_by_jet_flat,
        ak.num(arrays["JetsAK8_constituentsIndexCounts"], axis=1),
        axis=0,
    )

    good = arrays["JetsAK8_isGood"]
    index_by_good_jet = index_by_event[good]
    candidate_counts = ak.to_numpy(ak.flatten(ak.num(index_by_good_jet, axis=2), axis=None))
    total_candidates_by_event = ak.to_numpy(ak.sum(ak.num(index_by_good_jet, axis=2), axis=1))
    if int(np.sum(candidate_counts)) == 0:
        return {name: np.array([], dtype=np.float64) for name in VARIABLES}, np.array([], dtype=np.float64)

    constituent_counts_by_event = ak.to_numpy(ak.num(arrays["JetsConstituents_/.fPt"], axis=1))
    event_offsets = np.concatenate([[0], np.cumsum(constituent_counts_by_event)[:-1]])
    flat_global_index = ak.to_numpy(
        ak.flatten(index_by_good_jet + ak.broadcast_arrays(index_by_good_jet, ak.Array(event_offsets))[1], axis=None)
    )

    cand_pt = ak.to_numpy(ak.flatten(arrays["JetsConstituents_/.fPt"], axis=None))[flat_global_index]
    cand_eta = ak.to_numpy(ak.flatten(arrays["JetsConstituents_/.fEta"], axis=None))[flat_global_index]
    cand_phi = ak.to_numpy(ak.flatten(arrays["JetsConstituents_/.fPhi"], axis=None))[flat_global_index]
    cand_e = ak.to_numpy(ak.flatten(arrays["JetsConstituents_/.fE"], axis=None))[flat_global_index]

    jet_pt = np.repeat(ak.to_numpy(ak.flatten(arrays["JetsAK8_/.fPt"][good], axis=None)), candidate_counts)
    jet_eta = np.repeat(ak.to_numpy(ak.flatten(arrays["JetsAK8_/.fEta"][good], axis=None)), candidate_counts)
    jet_phi = np.repeat(ak.to_numpy(ak.flatten(arrays["JetsAK8_/.fPhi"][good], axis=None)), candidate_counts)
    jet_e = np.repeat(ak.to_numpy(ak.flatten(arrays["JetsAK8_/.fE"][good], axis=None)), candidate_counts)
    weights = np.repeat(event_weights, total_candidates_by_event)

    d_eta = cand_eta - jet_eta
    d_phi = np.asarray(delta_phi(ak.Array(cand_phi), ak.Array(jet_phi)), dtype=np.float64)
    values = {
        "del_eta": d_eta,
        "del_phi": d_phi,
        "del_r": np.sqrt(d_eta * d_eta + d_phi * d_phi),
        "log_pt": np.log(np.where(cand_pt > 0.0, cand_pt, np.nan)),
        "log_e": np.log(np.where(cand_e > 0.0, cand_e, np.nan)),
        "log_pt_jetpt": np.log(np.where((cand_pt > 0.0) & (jet_pt > 0.0), cand_pt / jet_pt, np.nan)),
        "log_e_jete": np.log(np.where((cand_e > 0.0) & (jet_e > 0.0), cand_e / jet_e, np.nan)),
    }
    return values, weights


def stream_file(
    file_url: str,
    year: str,
    step_size: str,
    accumulators: Dict[str, f2.HistogramAccumulator],
    is_data: bool = False,
) -> int:
    with uproot.open(file_url) as root_file:
        if "Events" not in root_file:
            raise KeyError("No Events tree found.")

        tree = root_file["Events"]
        requested = available_branches(tree, is_data=is_data)
        lund_norm = f2.lund_nominal_norm(root_file) if not is_data and "lundWeightNom" in requested else None
        entries = int(tree.num_entries)

        for arrays in tree.iterate(expressions=requested, step_size=step_size, library="ak"):
            if is_data:
                event_weights = np.ones(len(arrays), dtype=np.float64)
            else:
                event_weights = ak.to_numpy(arrays["Weight"]) * f2.LUMI_PB[year]
                if "puWeight" in arrays.fields:
                    event_weights *= ak.to_numpy(arrays["puWeight"])
                if "NonPrefiringProb" in arrays.fields:
                    event_weights *= ak.to_numpy(arrays["NonPrefiringProb"])
                if lund_norm is not None and "lundWeightNom" in arrays.fields:
                    event_weights *= ak.to_numpy(arrays["lundWeightNom"]) * lund_norm

            values, candidate_weights = pnet_candidate_values(arrays, event_weights)
            for variable, variable_values in values.items():
                accumulators[variable].fill(variable_values, candidate_weights)
    return entries


def read_sample(
    sample_dir: str,
    display_name: str,
    year: str,
    step_size: str,
    max_files: int,
    strict: bool,
    dry_run: bool,
    failures: List[str],
    is_data: bool = False,
) -> Optional[Dict[str, f2.HistogramAccumulator]]:
    try:
        files = f2.get_part_files(sample_dir, max_files)
    except Exception as exc:
        message = f"{display_name}: {exc}"
        if strict:
            raise RuntimeError(message) from exc
        print(f"[WARN] {message}")
        return None

    print(f"       part files: {len(files)}")
    if not files:
        message = f"{display_name}: no part-*.root files found"
        if strict:
            raise RuntimeError(message)
        print(f"[WARN] {message}")
        return None
    if dry_run:
        return f2.empty_accumulators()

    accumulators = f2.empty_accumulators()
    entries = 0
    successful = 0
    for index, remote_file in enumerate(files, start=1):
        try:
            entries += stream_file(
                f2.eos_url(remote_file),
                year,
                step_size,
                accumulators,
                is_data=is_data,
            )
            successful += 1
        except Exception as exc:
            message = f"{display_name}/{os.path.basename(remote_file)}: {exc}"
            failures.append(message)
            print(f"[WARN] {message}")

        if index % 25 == 0 or index == len(files):
            print(f"       processed {index:4d}/{len(files):4d} files")

    if strict and successful != len(files):
        raise RuntimeError(f"{display_name}: {len(files) - successful} input file(s) failed")

    print(
        f"       Events read: {entries:,} | "
        f"PF candidates: {accumulators['del_r'].entries_in_range:,} | "
        f"weighted yield: {accumulators['del_r'].integral():.8g}"
    )
    return accumulators


def draw_mode(output_path: str, output_dir: str, years: Sequence[str], signals: Sequence[str], log_y: bool, preliminary_label: bool) -> None:
    root_file = f2.ROOT.TFile.Open(output_path, "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Could not open cache ROOT file: {output_path}")
    try:
        for era in list(years) + ["Run2"]:
            for variable in VARIABLES:
                base = os.path.join(output_dir, era, "raw_wdata_ratio", variable)
                f2.draw_stack(
                    root_file,
                    era,
                    variable,
                    base,
                    signals,
                    normalized=False,
                    also_linear=False,
                    include_data=True,
                    include_ratio=True,
                    force_log_y=log_y,
                    preliminary_label=preliminary_label,
                )
    finally:
        root_file.Close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--background-base", default=f2.DEFAULT_BACKGROUND_BASE)
    parser.add_argument("--data-base", default=f2.DEFAULT_DATA_BASE)
    parser.add_argument("--signal-base", default=f2.DEFAULT_SIGNAL_BASE)
    parser.add_argument("--signals", default=",".join(PARTICLENET_SIGNAL_DIRS))
    parser.add_argument("--no-signals", action="store_true")
    parser.add_argument("--data-samples", default=",".join(f2.DEFAULT_DATA_SAMPLES))
    parser.add_argument("--no-data", action="store_true")
    parser.add_argument("--years", default="2016,2017,2018")
    parser.add_argument("--output", default="ParticleNetInput_plots/pnet_input_vars.root")
    parser.add_argument("--plot-dir", default="SVJ-tchannel-Run2_site/ParticleNET_Supervised_tagger/Inputvars")
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--chunk-size", default="100 MB")
    parser.add_argument("--max-files-per-sample", type=int, default=0)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--strict", action="store_true")
    parser.add_argument("--preliminary-label", action="store_true", help="draw CMS Preliminary or CMS Simulation Preliminary")
    return parser.parse_args()


def main() -> None:
    configure_figure2_globals()
    args = parse_args()
    years = f2.parse_csv(args.years)
    signals = [] if args.no_signals else f2.parse_csv(args.signals)
    data_samples = [] if args.no_data else f2.parse_csv(args.data_samples)

    if args.plot_only:
        draw_mode(args.output, os.path.join(args.plot_dir, "log"), years, signals, log_y=True, preliminary_label=args.preliminary_label)
        draw_mode(args.output, os.path.join(args.plot_dir, "linear"), years, signals, log_y=False, preliminary_label=args.preliminary_label)
        return

    failures: List[str] = []
    sample_hists: Dict[Tuple[str, str, str], Dict[str, f2.HistogramAccumulator]] = {}
    group_hists: Dict[Tuple[str, str], Dict[str, f2.HistogramAccumulator]] = {}
    signal_hists: Dict[Tuple[str, str], Dict[str, f2.HistogramAccumulator]] = {}
    data_hists: Dict[Tuple[str, str], Dict[str, f2.HistogramAccumulator]] = {}

    print("=" * 100)
    print("ParticleNet PF-candidate input variable cache and Figure2-style stacks")
    print(f"Background base : {f2.normalize_eos_path(args.background_base)}")
    print(f"Data base       : {f2.normalize_eos_path(args.data_base)}")
    print(f"Signal base     : {f2.normalize_eos_path(args.signal_base)}")
    print(f"Years           : {', '.join(years)}")
    print(f"ROOT output     : {args.output}")
    print(f"Plot directory  : {args.plot_dir}")
    print("=" * 100)

    for year in years:
        nominal = f"{f2.normalize_eos_path(args.background_base)}/{year}/t_channel_pre_selection/nominal"
        print(f"\n{'#' * 100}\nBACKGROUND {year}\n{nominal}\n{'#' * 100}")
        for process in f2.PROCESS_ORDER:
            for sample in f2.SAMPLE_MANIFEST[process]:
                sample_dir = f"{nominal}/{sample}"
                display_name = f"{year}/{process}/{sample}"
                print(f"\n[READ] {display_name}")
                accumulators = read_sample(
                    sample_dir,
                    display_name,
                    year,
                    args.chunk_size,
                    args.max_files_per_sample,
                    args.strict,
                    args.dry_run,
                    failures,
                )
                if accumulators is None:
                    continue
                sample_hists[(year, process, sample)] = accumulators
                group_hists.setdefault((year, process), f2.empty_accumulators())
                f2.add_accumulator_dict(group_hists[(year, process)], accumulators)

    if signals:
        for year in years:
            print(f"\n{'#' * 100}\nSIGNALS {year}\n{'#' * 100}")
            for signal_name, sample_dir in f2.get_signal_directories(args.signal_base, year, signals, args.strict):
                display_name = f"{year}/signal/{signal_name}"
                print(f"\n[READ] {display_name}")
                accumulators = read_sample(
                    sample_dir,
                    display_name,
                    year,
                    args.chunk_size,
                    args.max_files_per_sample,
                    args.strict,
                    args.dry_run,
                    failures,
                )
                if accumulators is not None:
                    signal_hists[(year, signal_name)] = accumulators

    if data_samples:
        for year in years:
            nominal = f"{f2.normalize_eos_path(args.data_base)}/{year}/t_channel_pre_selection/nominal"
            print(f"\n{'#' * 100}\nDATA {year}\n{nominal}\n{'#' * 100}")
            for sample in data_samples:
                display_name = f"{year}/data/{sample}"
                print(f"\n[READ] {display_name}")
                accumulators = read_sample(
                    f"{nominal}/{sample}",
                    display_name,
                    year,
                    args.chunk_size,
                    args.max_files_per_sample,
                    args.strict,
                    args.dry_run,
                    failures,
                    is_data=True,
                )
                if accumulators is not None:
                    data_hists[(year, sample)] = accumulators

    if args.dry_run:
        print("\n[DRY RUN] Nothing was written.")
        return
    if not group_hists:
        raise RuntimeError("No background histograms were made.")

    run2_hists: Dict[str, Dict[str, f2.HistogramAccumulator]] = {}
    for (_, process), accumulators in group_hists.items():
        run2_hists.setdefault(process, f2.empty_accumulators())
        f2.add_accumulator_dict(run2_hists[process], accumulators)

    run2_signal_hists: Dict[str, Dict[str, f2.HistogramAccumulator]] = {}
    for (_, signal_name), accumulators in signal_hists.items():
        run2_signal_hists.setdefault(signal_name, f2.empty_accumulators())
        f2.add_accumulator_dict(run2_signal_hists[signal_name], accumulators)

    run2_data_hists = f2.empty_accumulators() if data_hists else None
    if run2_data_hists:
        for accumulators in data_hists.values():
            f2.add_accumulator_dict(run2_data_hists, accumulators)

    f2.write_cache(
        args.output,
        sample_hists,
        group_hists,
        run2_hists,
        signal_hists,
        run2_signal_hists,
        data_hists,
        run2_data_hists,
        years,
    )
    print(f"\n[OK] Wrote ParticleNet input cache: {args.output}")

    if not args.no_plots:
        draw_mode(args.output, os.path.join(args.plot_dir, "log"), years, signals, log_y=True)
        draw_mode(args.output, os.path.join(args.plot_dir, "linear"), years, signals, log_y=False)

    if failures:
        print(f"\n[WARN] Completed with {len(failures)} file failure(s). First few:")
        for failure in failures[:10]:
            print(f"  - {failure}")


if __name__ == "__main__":
    main()
