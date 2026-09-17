#!/usr/bin/env python3
"""Build WNAE signal-efficiency caches from signal WNAE skims.

For each signal point this computes region-A yields with
  MET > 250 GeV, event-DNN score > 0.85,
and WNAE event nSVJ categories 0, 1, 2, 3+.

The output ROOT file is the expensive cache. CSVs are also written in the same
format used by make_signal_efficiency_maps_wnae.py.
"""

from __future__ import annotations

import argparse
import csv
import importlib
import json
import math
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import awkward as ak
import joblib
import numpy as np
import pandas as pd
import torch
import uproot

HERE = Path(__file__).resolve().parent
ANALYSIS = HERE.parents[1]
sys.path.insert(0, str(ANALYSIS))
sys.path.insert(0, str(HERE))

import WNAE_PNET_comparison as wnae
import make_signal_efficiency_maps as pnet_eff
import make_signal_efficiency_maps_extra_axes as pnet_eff_extra

SVJ_ORDER = ["0SVJ", "1SVJ", "2SVJ", "3PSVJ"]
YEARS = ["2016", "2017", "2018"]
DEFAULT_SIGNAL_BASE = "/store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/signals_WNAE"
DEFAULT_MODEL = "utils/data/DNNEventClassifier/sdt_allBkgs_disco_0p001_closure_0p06_damp_1_net_64_32_16_8_trainingAllYears_eval2016"
DEFAULT_OUTPUT = HERE / "wnae_signal_efficiency_cache.root"


def rinv_token(rinv: float) -> str:
    return "0" if abs(rinv) < 1e-12 else f"{rinv:g}".replace(".", "p")


def val_token(value: float) -> str:
    if abs(value - round(value)) < 1e-12:
        return str(int(round(value)))
    return f"{value:g}".replace(".", "p")


def sample_name(mmed: int, mdark: int, rinv: float, yukawa: float) -> str:
    return (
        f"t-channel_mMed-{mmed}_mDark-{mdark}_rinv-{rinv_token(rinv)}"
        f"_alpha-peak_yukawa-{val_token(yukawa)}"
    )


@dataclass(frozen=True)
class Point:
    scan: str
    mmed: int
    mdark: int
    rinv: float
    yukawa: float

    @property
    def sample(self) -> str:
        return sample_name(self.mmed, self.mdark, self.rinv, self.yukawa)


def build_points() -> List[Point]:
    points: Dict[Tuple[str, int, int, float, float], Point] = {}
    for mmed in pnet_eff.GRID_MMED:
        for rinv in pnet_eff.GRID_RINV:
            p = Point("mmed_rinv", int(mmed), 20, float(rinv), 1.0)
            points[(p.scan, p.mmed, p.mdark, p.rinv, p.yukawa)] = p
    for mmed in pnet_eff.GRID_MMED:
        for mdark in pnet_eff_extra.GRID_MDARK:
            p = Point("mdark", int(mmed), int(mdark), 0.3, 1.0)
            points[(p.scan, p.mmed, p.mdark, p.rinv, p.yukawa)] = p
    for mmed in pnet_eff.GRID_MMED:
        for yukawa in pnet_eff_extra.GRID_YUKAWA:
            p = Point("yukawa", int(mmed), 20, 0.3, float(yukawa))
            points[(p.scan, p.mmed, p.mdark, p.rinv, p.yukawa)] = p
    return list(points.values())


def load_event_tagger(model_dir: str):
    model_path = Path(model_dir)
    module_name = str(model_path / "training_model").replace("/", ".")
    training_model = importlib.import_module(module_name)
    model = training_model.model
    checkpoint = torch.load(model_path / "model.pt", map_location=torch.device("cpu"))
    state = checkpoint["model_state_dict"] if isinstance(checkpoint, dict) and "model_state_dict" in checkpoint else checkpoint
    model.load_state_dict(state)
    model.eval()
    scaler = joblib.load(model_path / "scaler.joblib")
    return model, scaler


def pad_array(array, idx: int, pad_value: float = 0.0):
    counts = ak.count(array, axis=1)
    masked = ak.mask(array, counts > idx)
    return ak.fill_none(masked[:, idx], pad_value)


def build_dnn_inputs(arrays, scaler) -> pd.DataFrame:
    var: Dict[str, object] = {field: arrays[field] for field in arrays.fields}
    good = var["JetsAK8_isGood"]
    var["JetsAK8_MTMETLog"] = np.log(var["JetsAK8_MTMET"])

    variable_names = ["deltaPhiMET", "LundJetPlaneZ", "MTMETLog", "/.fPt", "/.fEta", "/.fPhi", "mass"]
    for variable_name in variable_names:
        if variable_name == "/.fPt":
            target = "pt"
        elif variable_name == "/.fEta":
            target = "eta"
        elif variable_name == "/.fPhi":
            target = "phi"
        else:
            target = variable_name
        source = f"JetsAK8_{variable_name}"
        good_values = var[source][good]
        for idx in range(8):
            name = f"GoodJetsAK8{idx}_{target}"
            var[name] = pad_array(good_values, idx)
            if variable_name == "deltaPhiMET":
                var[name] = ak.where(var[name] == 0, var["METPhi"], var[name])

    for idx_1 in range(4):
        for idx_2 in range(idx_1 + 1, 4):
            base = f"DijetMass{idx_1}{idx_2}GoodJetsAK8"
            dijet_mass = var[base]
            jet1_mass = var[f"GoodJetsAK8{idx_1}_mass"]
            jet2_mass = var[f"GoodJetsAK8{idx_2}_mass"]
            dijet_mass = ak.where(dijet_mass == 0, jet1_mass, dijet_mass)
            dijet_mass = ak.where(dijet_mass == 0, jet2_mass, dijet_mass)
            log_mass = np.log(dijet_mass)
            var[f"{base}Log"] = ak.nan_to_num(log_mass, nan=0.0, posinf=0.0, neginf=0.0)

    features = list(scaler.feature_names_in_)
    missing = [name for name in features if name not in var]
    if missing:
        raise KeyError(f"Missing DNN feature(s): {missing}")
    frame = pd.DataFrame({name: np.asarray(ak.to_numpy(var[name]), dtype="float64") for name in features})
    scaled = scaler.transform(frame)
    return pd.DataFrame(scaled, columns=features)


def evaluate_dnn(arrays, model, scaler) -> np.ndarray:
    inputs = build_dnn_inputs(arrays, scaler)
    with torch.no_grad():
        scores = model(torch.tensor(inputs.values, dtype=torch.float32)).detach().numpy().flatten()
    return np.asarray(scores, dtype="float64")


def required_branches(tree) -> List[str]:
    available = set(tree.keys())
    branches = [
        "MET",
        "METPhi",
        "JetsAK8_isGood",
        "JetsAK8_deltaPhiMET",
        "JetsAK8_LundJetPlaneZ",
        "JetsAK8_MTMET",
        "JetsAK8_/.fPt",
        "JetsAK8_/.fEta",
        "JetsAK8_/.fPhi",
        "JetsAK8_mass",
    ]
    for i in range(4):
        for j in range(i + 1, 4):
            branches.extend(
                [
                    f"DeltaEta{i}{j}GoodJetsAK8",
                    f"DeltaPhi{i}{j}GoodJetsAK8",
                    f"DeltaR{i}{j}GoodJetsAK8",
                    f"LundJetPlaneZ{i}{j}GoodJetsAK8",
                    f"DijetMass{i}{j}GoodJetsAK8",
                ]
            )
    branches.extend(pt_bin.loss_branch for pt_bin in wnae.WNAE_PT_BINS)
    branches.extend(
        wnae.selected_weight_branches(
            tree,
            is_data=False,
            no_weights=False,
            no_lund_correction=False,
        )
    )
    missing = [branch for branch in branches if branch not in available]
    if missing:
        raise KeyError(f"Missing required branch(es): {missing}")
    return list(dict.fromkeys(branches))


def read_initial(root_file) -> float:
    if "CutFlow" not in root_file:
        return 0.0
    cutflow = root_file["CutFlow"]
    if "Initial" not in cutflow.keys():
        return 0.0
    return float(cutflow["Initial"].array(library="np")[0])


def stream_file(file_url: str, year: str, model, scaler, chunk_size: str, max_events: int):
    sums = np.zeros(4, dtype="float64")
    sums2 = np.zeros(4, dtype="float64")
    raw = np.zeros(4, dtype="int64")
    entries_seen = 0
    entries_region = 0
    initial = 0.0

    with uproot.open(file_url) as root_file:
        initial = read_initial(root_file)
        tree = root_file["Events"]
        branches = required_branches(tree)
        lund_norm = wnae.lund_nominal_norm(root_file)
        stop = max_events if max_events and max_events > 0 else None
        for arrays in tree.iterate(branches, step_size=chunk_size, library="ak", entry_stop=stop):
            n = len(arrays)
            entries_seen += n
            scores = evaluate_dnn(arrays, model, scaler)
            event_weights = wnae.chunk_event_weights(
                arrays, year=year, is_data=False, no_weights=False, lund_norm=lund_norm
            )
            region_mask = (np.asarray(ak.to_numpy(arrays["MET"]), dtype="float64") > 250.0) & (scores > 0.85)
            entries_region += int(np.count_nonzero(region_mask))

            pt = arrays[wnae.PT_BRANCH]
            good = arrays[wnae.GOOD_JET_BRANCH]
            pt = pt[good]
            wnae_counts = np.zeros(n, dtype="int64")
            for pt_bin in wnae.WNAE_PT_BINS:
                loss = arrays[pt_bin.loss_branch][good]
                if pt_bin.pt_max is None:
                    in_pt = pt >= pt_bin.pt_min
                else:
                    in_pt = (pt >= pt_bin.pt_min) & (pt < pt_bin.pt_max)
                tags = in_pt & (loss >= pt_bin.mc_wp)
                wnae_counts += ak.to_numpy(ak.sum(tags, axis=1))

            cats = np.minimum(wnae_counts, 3)
            for idx in range(4):
                mask = region_mask & (cats == idx)
                if not np.any(mask):
                    continue
                selected_weights = event_weights[mask]
                sums[idx] += float(np.sum(selected_weights))
                sums2[idx] += float(np.sum(selected_weights * selected_weights))
                raw[idx] += int(np.count_nonzero(mask))

    return initial, sums, sums2, raw, entries_seen, entries_region


def write_root_cache(path: Path, rows: List[dict], metadata: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    arrays = {}
    fields = [
        "scan",
        "sample",
        "mMed",
        "mDark",
        "rinv",
        "yukawa",
        "N_generated",
        "entries_seen",
        "entries_region",
        "n_files",
        "n_failed",
    ]
    fields += [f"NA_{svj}" for svj in SVJ_ORDER] + ["NA_Inclusive"]
    fields += [f"sumw2_{svj}" for svj in SVJ_ORDER] + ["sumw2_Inclusive"]
    fields += [f"raw_{svj}" for svj in SVJ_ORDER] + ["raw_Inclusive"]
    fields += [f"eff_{svj}" for svj in SVJ_ORDER] + ["eff_Inclusive"]
    for field in fields:
        values = [row[field] for row in rows]
        if field in ("scan", "sample"):
            arrays[field] = np.asarray(values, dtype=object)
        elif field in ("mMed", "mDark", "entries_seen", "entries_region", "n_files", "n_failed") or field.startswith("raw_"):
            arrays[field] = np.asarray(values, dtype="int64")
        else:
            arrays[field] = np.asarray(values, dtype="float64")
    with uproot.recreate(path) as fout:
        fout["results"] = arrays
        fout["metadata"] = {"json": np.asarray([json.dumps(metadata, sort_keys=True)], dtype=object)}


def write_csvs(rows: List[dict]) -> None:
    row_by_scan_point = {
        (r["scan"], r["mMed"], r["mDark"], round(r["rinv"], 6), round(r["yukawa"], 6)): r
        for r in rows
    }

    def csv_row(row):
        out = {k: row[k] for k in ["N_generated"] + [f"NA_{s}" for s in SVJ_ORDER] + ["NA_Inclusive"] + [f"eff_{s}" for s in SVJ_ORDER] + ["eff_Inclusive"]}
        return out

    grid_path = HERE / "wnae_signal_efficiency_grid.csv"
    with open(grid_path, "w", newline="") as fout:
        fieldnames = ["mMed", "rinv", "N_generated"] + [f"NA_{s}" for s in SVJ_ORDER] + ["NA_Inclusive"] + [f"eff_{s}" for s in SVJ_ORDER] + ["eff_Inclusive"]
        writer = csv.DictWriter(fout, fieldnames=fieldnames)
        writer.writeheader()
        for mmed in pnet_eff.GRID_MMED:
            for rinv in pnet_eff.GRID_RINV:
                row = {"mMed": mmed, "rinv": rinv}
                cached = row_by_scan_point.get(("mmed_rinv", int(mmed), 20, round(float(rinv), 6), 1.0))
                if cached:
                    row.update(csv_row(cached))
                writer.writerow(row)

    mdark_path = HERE / "wnae_signal_efficiency_mdark_scan.csv"
    with open(mdark_path, "w", newline="") as fout:
        fieldnames = ["mMed", "mDark", "N_generated"] + [f"NA_{s}" for s in SVJ_ORDER] + ["NA_Inclusive"] + [f"eff_{s}" for s in SVJ_ORDER] + ["eff_Inclusive"]
        writer = csv.DictWriter(fout, fieldnames=fieldnames)
        writer.writeheader()
        for mmed in pnet_eff.GRID_MMED:
            for mdark in pnet_eff_extra.GRID_MDARK:
                row = {"mMed": mmed, "mDark": mdark}
                cached = row_by_scan_point.get(("mdark", int(mmed), int(mdark), 0.3, 1.0))
                if cached:
                    row.update(csv_row(cached))
                writer.writerow(row)

    yukawa_path = HERE / "wnae_signal_efficiency_yukawa_scan.csv"
    with open(yukawa_path, "w", newline="") as fout:
        fieldnames = ["mMed", "yukawa", "N_generated"] + [f"NA_{s}" for s in SVJ_ORDER] + ["NA_Inclusive"] + [f"eff_{s}" for s in SVJ_ORDER] + ["eff_Inclusive"]
        writer = csv.DictWriter(fout, fieldnames=fieldnames)
        writer.writeheader()
        for mmed in pnet_eff.GRID_MMED:
            for yukawa in pnet_eff_extra.GRID_YUKAWA:
                row = {"mMed": mmed, "yukawa": yukawa}
                cached = row_by_scan_point.get(("yukawa", int(mmed), 20, 0.3, round(float(yukawa), 6)))
                if cached:
                    row.update(csv_row(cached))
                writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", default=DEFAULT_SIGNAL_BASE)
    parser.add_argument("--years", default=",".join(YEARS))
    parser.add_argument("--model-dir", default=DEFAULT_MODEL)
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--chunk-size", default="100 MB")
    parser.add_argument("--max-events", type=int, default=-1)
    parser.add_argument("--max-files-per-sample", type=int, default=-1)
    parser.add_argument("--only-sample", action="append", default=[])
    args = parser.parse_args()

    years = [y.strip() for y in args.years.split(",") if y.strip()]
    model, scaler = load_event_tagger(args.model_dir)
    points = build_points()
    if args.only_sample:
        wanted = set(args.only_sample)
        points = [p for p in points if p.sample in wanted]

    rows = []
    for ip, point in enumerate(points, 1):
        discovered = wnae.discover_files_from_base(
            base=args.base,
            years=years,
            processes=["Signal"],
            sample_dirs=[point.sample],
            max_files_per_sample=args.max_files_per_sample,
        )
        sums = np.zeros(4, dtype="float64")
        sums2 = np.zeros(4, dtype="float64")
        raw = np.zeros(4, dtype="int64")
        entries_seen = 0
        entries_region = 0
        n_generated = 0.0
        n_failed = 0
        print(f"[{ip}/{len(points)}] {point.sample}: {len(discovered)} files")
        for year, _process, _sample, file_url in discovered:
            try:
                initial, s, s2, r, seen, region = stream_file(
                    file_url, year, model, scaler, args.chunk_size, args.max_events
                )
            except Exception as exc:
                n_failed += 1
                print(f"[WARN] failed {file_url}: {exc}", flush=True)
                continue
            n_generated += initial * wnae.LUMI_PB[year]
            sums += s
            sums2 += s2
            raw += r
            entries_seen += seen
            entries_region += region

        row = {
            "scan": point.scan,
            "sample": point.sample,
            "mMed": point.mmed,
            "mDark": point.mdark,
            "rinv": point.rinv,
            "yukawa": point.yukawa,
            "N_generated": n_generated,
            "entries_seen": entries_seen,
            "entries_region": entries_region,
            "n_files": len(discovered),
            "n_failed": n_failed,
        }
        for idx, svj in enumerate(SVJ_ORDER):
            row[f"NA_{svj}"] = sums[idx]
            row[f"sumw2_{svj}"] = sums2[idx]
            row[f"raw_{svj}"] = int(raw[idx])
            row[f"eff_{svj}"] = sums[idx] / n_generated if n_generated > 0 else math.nan
        row["NA_Inclusive"] = float(np.sum(sums))
        row["sumw2_Inclusive"] = float(np.sum(sums2))
        row["raw_Inclusive"] = int(np.sum(raw))
        row["eff_Inclusive"] = row["NA_Inclusive"] / n_generated if n_generated > 0 else math.nan
        rows.append(row)
        print(
            f"    NA incl={row['NA_Inclusive']:.6g}, eff incl={100.0 * row['eff_Inclusive']:.5g}%, "
            f"region entries={entries_region}",
            flush=True,
        )

    metadata = {
        "base": args.base,
        "years": years,
        "model_dir": args.model_dir,
        "selection": "MET > 250 and dnnEventClassScore > 0.85",
        "weights": "Weight*lumi*puWeight*NonPrefiringProb*lundWeightNom*(Initial/InitialLundNominal)",
        "wnae_pt_bins": [pt_bin.__dict__ for pt_bin in wnae.WNAE_PT_BINS],
    }
    write_root_cache(Path(args.output), rows, metadata)
    write_csvs(rows)
    print(f"[OK] wrote ROOT cache: {args.output}")
    print(f"[OK] wrote WNAE efficiency CSVs in: {HERE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
