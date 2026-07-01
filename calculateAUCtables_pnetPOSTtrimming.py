#!/usr/bin/env python3
"""
python3 calculateAUCtables_pnetPOSTtrimming.py  --no-cache

Compute event-level ParticleNet AUC for EVERY discovered 2018 t-channel signal
point satisfying the requested signal filters, against every discovered 2018
WNAE + MET-trimmed QCD pT bin and the combined QCD sample.

Default signal selection
------------------------
  signal base : /store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/
                signals/2018/nominal
  mDark       : 20
  yukawa      : 1
  alpha       : peak

With the supplied directory listing, this selects the 77 available
(mMed, rinv) signal points.

The event-level ParticleNet score is:

    max(JetsAK8_pNetJetTaggerScore)

across valid AK8 jets in each event. Scores below --invalid-score-cut are
interpreted as invalid sentinels. Events without any valid AK8 score are
excluded from the AUC calculation.

Outputs
-------
  <out-dir>/pnet_auc_table.csv
      Detailed table: each selected signal point vs each QCD pT bin plus
      ALL_QCD_COMBINED.

  <out-dir>/pnet_auc_combined_qcd.csv
      One row per selected signal point vs ALL_QCD_COMBINED.

  <out-dir>/pnet_auc_table.txt
      Human-readable summary, raw/oriented AUC grids in (mMed, rinv), and
      the detailed signal x QCD-bin table.

Examples
--------
Quick validation (read at most 10k events from each signal/QCD bin):

  python3 pnet_auc_all_signals_mDark20_yukawa1.py \\
    --out-dir pnet_auc_test_all_signals \\
    --max-events-per-sample 10000 \\
    --auc-cap 10000 \\
    --no-cache

Full run (read all events; use <=500k finite events/class per AUC):

  python3 pnet_auc_all_signals_mDark20_yukawa1.py \\
    --out-dir pnet_auc_2018_allSignals_mDark20_yukawa1

Use every finite event in every AUC:

python3 calculateAUCtables_pnetPOSTtrimming.py  --no-cache
Restrict to a small set of signal points for debugging:

  python3 pnet_auc_all_signals_mDark20_yukawa1.py \\
    --mmeds 1000 2000 \\
    --rinvs 0.1 0.3 \\
    --out-dir pnet_auc_debug
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import awkward as ak
import numpy as np
import uproot
from sklearn.metrics import roc_auc_score


# -----------------------------------------------------------------------------
# Defaults
# -----------------------------------------------------------------------------
EOS_HOST = "root://cmseos.fnal.gov"
DEFAULT_TREE = "Events"
DEFAULT_PNET_BRANCH = "JetsAK8_pNetJetTaggerScore"

# DEFAULT_SIGNAL_EOS_BASE = (
#     "/store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/"
#     "signals/2018/nominal"
# )

DEFAULT_SIGNAL_EOS_BASE=(
    "/store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/2018/t_channel_pre_selection/nominal"
)
DEFAULT_QCD_EOS_BASE = (
    "/store/user/ashrivas/tchannel_UL/skims_gapJetVeto/"
    "data_mc_withWNAE_METtrim/trimmed_root"
)
DEFAULT_TTBAR_EOS_BASE = (
    "/store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/"
    "data_mc_noWNAE"
)

NO_WNAE_BACKGROUND_CONFIG = {
    "ttjets": {
        "label": "TTJets",
        "combined_name": "ALL_TTJETS_COMBINED",
        "prefixes": ("TTJets",),
    },
    "st": {
        "label": "ST",
        "combined_name": "ALL_ST_COMBINED",
        "prefixes": ("ST_",),
    },
    "wjets": {
        "label": "WJets",
        "combined_name": "ALL_WJETS_COMBINED",
        "prefixes": ("WJetsToLNu",),
    },
    "zjets": {
        "label": "ZJets",
        "combined_name": "ALL_ZJETS_COMBINED",
        "prefixes": ("ZJetsToNuNu",),
    },
}

SIGNAL_DIR_RE = re.compile(
    r"^t-channel_"
    r"mMed-(?P<mmed>[\dp]+)_"
    r"mDark-(?P<mdark>[\dp]+)_"
    r"rinv-(?P<rinv>[\dp]+)_"
    r"alpha-(?P<alpha>[^_]+)_"
    r"yukawa-(?P<yukawa>[^_]+)$"
)


# -----------------------------------------------------------------------------
# EOS helpers
# -----------------------------------------------------------------------------
def normalise_eos_path(path: str) -> str:
    """Accept either /store/... or LPC-mounted /eos/uscms/store/... paths."""
    path = os.path.expanduser(path).rstrip("/")
    for mount in ("/eos/uscms", "/eos/cms"):
        if path == mount:
            return "/"
        if path.startswith(mount + "/"):
            return path[len(mount):]
    return path


def xrootd_url(eos_path: str) -> str:
    """Convert an EOS path to a canonical XRootD URL."""
    return f"{EOS_HOST}//{normalise_eos_path(eos_path).lstrip('/')}"


def xrdfs_ls(eos_path: str) -> List[str]:
    """List one EOS directory and return its canonical /store/... entries."""
    eos_path = normalise_eos_path(eos_path)
    command = ["xrdfs", EOS_HOST, "ls", eos_path]
    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"xrdfs failed for {eos_path}\n"
            f"Command: {' '.join(command)}\n"
            f"{result.stderr.strip()}"
        )
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def discover_root_files(eos_dir: str) -> List[str]:
    """Return direct ROOT-file children of one EOS directory."""
    return sorted(
        xrootd_url(entry)
        for entry in xrdfs_ls(eos_dir)
        if entry.endswith(".root")
    )


# -----------------------------------------------------------------------------
# Dataset metadata and discovery
# -----------------------------------------------------------------------------
def token_to_float(value: str) -> float:
    """Convert EOS naming tokens such as 0p05 into Python floats."""
    return float(value.replace("p", "."))


def normalise_float_token(value: float) -> str:
    """Stable human-readable floating-point text used in tables."""
    return f"{value:.12g}"


def parse_signal_metadata(signal_dir: str) -> Optional[Dict[str, object]]:
    """Parse parameters from one t-channel signal-directory basename."""
    signal_name = os.path.basename(normalise_eos_path(signal_dir))
    match = SIGNAL_DIR_RE.match(signal_name)
    if match is None:
        return None

    groups = match.groupdict()
    try:
        return {
            "signal_name": signal_name,
            "mMed": token_to_float(groups["mmed"]),
            "mDark": token_to_float(groups["mdark"]),
            "rinv": token_to_float(groups["rinv"]),
            "alpha": groups["alpha"],
            "yukawa": token_to_float(groups["yukawa"]),
        }
    except ValueError:
        return None


@dataclass(frozen=True)
class SignalSample:
    name: str
    eos_dir: str
    files: Tuple[str, ...]
    metadata: Dict[str, object]


@dataclass(frozen=True)
class QCDSample:
    name: str
    eos_dir: str
    files: Tuple[str, ...]


def signal_sort_key(sample: SignalSample) -> Tuple[float, float, str]:
    return (
        float(sample.metadata["mMed"]),
        float(sample.metadata["rinv"]),
        sample.name,
    )


def qcd_sort_key(name: str) -> Tuple[int, int]:
    """Sort QCD_Pt_170to300 through QCD_Pt_3200toInf numerically."""
    match = re.match(r"QCD_Pt_(\d+)to(\d+|Inf)$", name)
    if match is None:
        return (10**9, 10**9)
    low = int(match.group(1))
    high = 10**9 if match.group(2) == "Inf" else int(match.group(2))
    return (low, high)


def _matches_requested_float(value: float, requested: Optional[Sequence[float]]) -> bool:
    if not requested:
        return True
    return any(np.isclose(value, target, rtol=0.0, atol=1.0e-10) for target in requested)


def discover_signal_samples(
    signal_base: str,
    *,
    m_dark: float,
    yukawa: float,
    alpha: str,
    requested_mmeds: Optional[Sequence[float]],
    requested_rinvs: Optional[Sequence[float]],
) -> List[SignalSample]:
    """Discover signals and retain points that pass all requested filters."""
    samples: List[SignalSample] = []

    for entry in xrdfs_ls(signal_base):
        metadata = parse_signal_metadata(entry)
        if metadata is None:
            continue

        if not np.isclose(float(metadata["mDark"]), m_dark, rtol=0.0, atol=1.0e-10):
            continue
        if not np.isclose(float(metadata["yukawa"]), yukawa, rtol=0.0, atol=1.0e-10):
            continue
        if alpha.lower() != "any" and str(metadata["alpha"]) != alpha:
            continue
        if not _matches_requested_float(float(metadata["mMed"]), requested_mmeds):
            continue
        if not _matches_requested_float(float(metadata["rinv"]), requested_rinvs):
            continue

        files = tuple(discover_root_files(entry))
        if not files:
            print(f"[WARNING] No ROOT files found for signal: {entry}", file=sys.stderr)
            continue
        samples.append(
            SignalSample(
                name=str(metadata["signal_name"]),
                eos_dir=entry,
                files=files,
                metadata=metadata,
            )
        )

    return sorted(samples, key=signal_sort_key)


def discover_trimmed_qcd_samples(
    qcd_base: str,
    year: str,
    requested_bins: Optional[Sequence[str]],
) -> List[QCDSample]:
    """Discover <year>_QCD_Pt_* bins in the trimmed-QCD directory layout."""
    wanted = set(requested_bins or [])
    candidates: List[Tuple[str, str]] = []

    for entry in xrdfs_ls(qcd_base):
        basename = os.path.basename(entry.rstrip("/"))
        prefix = f"{year}_"
        if basename.startswith(prefix + "QCD_Pt_"):
            candidates.append((basename[len(prefix):], entry))

    samples: List[QCDSample] = []
    for qcd_name, eos_dir in sorted(candidates, key=lambda item: qcd_sort_key(item[0])):
        if wanted and qcd_name not in wanted:
            continue
        files = tuple(discover_root_files(eos_dir))
        if not files:
            print(f"[WARNING] No ROOT files found for QCD: {eos_dir}", file=sys.stderr)
            continue
        samples.append(QCDSample(name=qcd_name, eos_dir=eos_dir, files=files))

    if wanted:
        found = {sample.name for sample in samples}
        missing = sorted(wanted - found, key=qcd_sort_key)
        if missing:
            print(
                "[WARNING] Requested QCD bins not found: " + ", ".join(missing),
                file=sys.stderr,
            )
    return samples


def discover_no_wnae_background_samples(
    background_base: str,
    year: str,
    background_kind: str,
    requested_samples: Optional[Sequence[str]],
) -> List[QCDSample]:
    """Discover a no-WNAE MC background family in YEAR/t_channel_pre_selection/nominal."""
    config = NO_WNAE_BACKGROUND_CONFIG[background_kind]
    wanted = set(requested_samples or [])
    nominal_dir = os.path.join(
        normalise_eos_path(background_base),
        year,
        "t_channel_pre_selection",
        "nominal",
    )

    samples: List[QCDSample] = []
    for entry in xrdfs_ls(nominal_dir):
        sample_name = os.path.basename(entry.rstrip("/"))
        if not sample_name.startswith(config["prefixes"]):
            continue
        if wanted and sample_name not in wanted:
            continue
        files = tuple(discover_root_files(entry))
        if not files:
            print(f"[WARNING] No ROOT files found for {config['label']}: {entry}", file=sys.stderr)
            continue
        samples.append(QCDSample(name=sample_name, eos_dir=entry, files=files))

    if wanted:
        found = {sample.name for sample in samples}
        missing = sorted(wanted - found)
        if missing:
            print(
                f"[WARNING] Requested {config['label']} samples not found: "
                + ", ".join(missing),
                file=sys.stderr,
            )
    return sorted(samples, key=lambda sample: sample.name)


# -----------------------------------------------------------------------------
# ROOT score extraction
# -----------------------------------------------------------------------------
def find_tree(fin: uproot.ReadOnlyDirectory) -> Optional[str]:
    """Use common TTree names first, then find any available TTree."""
    for candidate in ("Events", "TreeMaker2/PreSelection", "PreSelection", "tree"):
        if candidate in fin:
            return candidate

    for key, obj in fin.items(recursive=True):
        if isinstance(obj, uproot.behaviors.TTree.TTree):
            return key.split(";")[0]
    return None


def load_event_max_scores(
    files: Sequence[str],
    *,
    branch: str,
    requested_tree: str,
    step_size: int,
    max_events: int,
    invalid_score_cut: float,
    label: str,
) -> np.ndarray:
    """Read one max-valid pNet score per event from direct ROOT file paths."""
    score_chunks: List[np.ndarray] = []
    events_read = 0
    files_total = len(files)

    for file_index, file_path in enumerate(files, start=1):
        if max_events > 0 and events_read >= max_events:
            break

        try:
            with uproot.open(file_path) as fin:
                tree_name = requested_tree if requested_tree in fin else find_tree(fin)
                if tree_name is None:
                    print(f"[WARNING] {label}: no TTree in {file_path}", file=sys.stderr)
                    continue

                tree = fin[tree_name]
                if branch not in tree:
                    print(
                        f"[WARNING] {label}: missing branch '{branch}' in {file_path}",
                        file=sys.stderr,
                    )
                    continue

                for arrays in tree.iterate([branch], step_size=step_size, library="ak"):
                    raw_scores = arrays[branch]

                    # A jagged event can contain a valid jet even when another
                    # jet has the invalid sentinel. Mask before the per-event max.
                    valid_scores = ak.mask(raw_scores, raw_scores >= invalid_score_cut)
                    event_scores = ak.to_numpy(
                        ak.fill_none(
                            ak.max(valid_scores, axis=1, mask_identity=True),
                            np.nan,
                        )
                    ).astype(np.float64, copy=False)

                    if max_events > 0:
                        remaining = max_events - events_read
                        event_scores = event_scores[:remaining]

                    score_chunks.append(event_scores)
                    events_read += len(event_scores)
                    if max_events > 0 and events_read >= max_events:
                        break

        except Exception as exc:
            # A transient EOS/XRootD failure should not discard a multi-hour job.
            print(f"[WARNING] {label}: failed to read {file_path}: {exc}", file=sys.stderr)

        if file_index == 1 or file_index == files_total or file_index % 10 == 0:
            print(f"  {label}: {file_index}/{files_total} files, {events_read:,} events read")

    if not score_chunks:
        return np.empty(0, dtype=np.float64)
    return np.concatenate(score_chunks)


def finite_scores(scores: np.ndarray) -> np.ndarray:
    scores = np.asarray(scores, dtype=np.float64)
    return scores[np.isfinite(scores)]


# -----------------------------------------------------------------------------
# Score cache helpers
# -----------------------------------------------------------------------------
def fingerprint(*parts: object) -> str:
    """Small stable cache key derived from the exact input/configuration."""
    joined = "|".join(map(str, parts))
    return hashlib.sha1(joined.encode("utf-8")).hexdigest()[:16]


def safe_filename(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text)


def load_or_build_scores(
    *,
    files: Sequence[str],
    label: str,
    cache_path: str,
    cache_key: str,
    args: argparse.Namespace,
) -> np.ndarray:
    """Load an exact cache match, otherwise stream scores and create one."""
    if not args.no_cache and os.path.exists(cache_path):
        try:
            with np.load(cache_path, allow_pickle=False) as cache:
                saved_key = str(cache["cache_key"].item())
                if saved_key == cache_key:
                    scores = np.asarray(cache["scores"], dtype=np.float64)
                    print(f"Loading cached {label} scores: {cache_path}")
                    return scores
                print(f"Cache key changed; rebuilding {label}: {cache_path}")
        except Exception as exc:
            print(f"[WARNING] Cannot reuse cache {cache_path}: {exc}", file=sys.stderr)

    scores = load_event_max_scores(
        files,
        branch=args.pnet_branch,
        requested_tree=args.tree,
        step_size=args.step_size,
        max_events=args.max_events_per_sample,
        invalid_score_cut=args.invalid_score_cut,
        label=label,
    )

    if len(scores) and not args.no_cache:
        np.savez_compressed(cache_path, scores=scores, cache_key=np.asarray(cache_key))
        print(f"Saved {label} cache: {cache_path}")
    return scores


# -----------------------------------------------------------------------------
# AUC calculation
# -----------------------------------------------------------------------------
def deterministic_rng(seed: int, label: str) -> np.random.Generator:
    """Independent deterministic random state for each output row."""
    digest = hashlib.sha1(label.encode("utf-8")).hexdigest()[:8]
    return np.random.default_rng(seed + int(digest, 16))


def maybe_downsample(values: np.ndarray, cap: int, rng: np.random.Generator) -> np.ndarray:
    if cap < 0 or len(values) <= cap:
        return values
    indices = rng.choice(len(values), size=cap, replace=False)
    return values[indices]


def compute_auc(
    signal_scores: np.ndarray,
    qcd_scores: np.ndarray,
    auc_cap: int,
    rng: np.random.Generator,
) -> Dict[str, object]:
    """Calculate raw and direction-independent AUC with event-count bookkeeping."""
    signal = finite_scores(signal_scores)
    qcd = finite_scores(qcd_scores)

    result: Dict[str, object] = {
        "n_signal_finite": int(len(signal)),
        "n_qcd_finite": int(len(qcd)),
        "n_signal_used": 0,
        "n_qcd_used": 0,
        "auc_raw_signal_high": float("nan"),
        "auc_oriented": float("nan"),
        "orientation": "",
    }
    if len(signal) == 0 or len(qcd) == 0:
        return result

    signal_used = maybe_downsample(signal, auc_cap, rng)
    qcd_used = maybe_downsample(qcd, auc_cap, rng)

    labels = np.concatenate(
        [np.ones(len(signal_used), dtype=np.int8), np.zeros(len(qcd_used), dtype=np.int8)]
    )
    values = np.concatenate([signal_used, qcd_used])
    raw_auc = float(roc_auc_score(labels, values))

    # The raw value retains the intended pNet convention: larger score is more
    # signal-like. The oriented value is just a separation diagnostic.
    orientation = "+score" if raw_auc >= 0.5 else "-score"
    oriented_auc = raw_auc if raw_auc >= 0.5 else 1.0 - raw_auc

    result.update(
        n_signal_used=int(len(signal_used)),
        n_qcd_used=int(len(qcd_used)),
        auc_raw_signal_high=raw_auc,
        auc_oriented=oriented_auc,
        orientation=orientation,
    )
    return result


# -----------------------------------------------------------------------------
# Output helpers
# -----------------------------------------------------------------------------
DETAIL_FIELDS = [
    "signal_name",
    "signal_eos_dir",
    "signal_files",
    "mMed",
    "mDark",
    "rinv",
    "alpha",
    "yukawa",
    "qcd_sample",
    "qcd_files",
    "n_signal_finite",
    "n_qcd_finite",
    "n_signal_used",
    "n_qcd_used",
    "auc_raw_signal_high",
    "auc_oriented",
    "orientation",
]

COMBINED_FIELDS = DETAIL_FIELDS.copy()


def format_auc(value: object, width: int = 0) -> str:
    number = float(value)
    text = f"{number:.6f}" if np.isfinite(number) else "nan"
    return f"{text:>{width}}" if width else text


def write_csv(path: str, fields: Sequence[str], rows: Sequence[Dict[str, object]]) -> None:
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields))
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def combined_grid(
    combined_rows: Sequence[Dict[str, object]],
    value_field: str,
) -> Tuple[List[float], List[float], Dict[Tuple[float, float], float]]:
    mmeds = sorted({float(row["mMed"]) for row in combined_rows})
    rinvs = sorted({float(row["rinv"]) for row in combined_rows})
    cells: Dict[Tuple[float, float], float] = {}
    for row in combined_rows:
        cells[(float(row["mMed"]), float(row["rinv"]))] = float(row[value_field])
    return mmeds, rinvs, cells


def write_auc_grid(
    handle: object,
    *,
    title: str,
    combined_rows: Sequence[Dict[str, object]],
    value_field: str,
) -> None:
    mmeds, rinvs, cells = combined_grid(combined_rows, value_field)
    handle.write(title + "\n")
    header = f"{'mMed [GeV]':>11} |" + "".join(f" {normalise_float_token(r):>9}" for r in rinvs)
    handle.write(header + "\n")
    handle.write("-" * len(header) + "\n")
    for mmed in mmeds:
        values: List[str] = []
        for rinv in rinvs:
            value = cells.get((mmed, rinv))
            values.append("       --" if value is None else f" {format_auc(value, 8)}")
        handle.write(f"{mmed:11.0f} |" + "".join(values) + "\n")
    handle.write("\n")


def write_text_table(
    *,
    path: str,
    args: argparse.Namespace,
    background_kind: str,
    background_label: str,
    combined_background_name: str,
    background_base: str,
    signal_samples: Sequence[SignalSample],
    qcd_samples: Sequence[QCDSample],
    combined_rows: Sequence[Dict[str, object]],
    detail_rows: Sequence[Dict[str, object]],
) -> None:
    total_qcd_files = sum(len(sample.files) for sample in qcd_samples)

    with open(path, "w") as handle:
        handle.write(f"ParticleNet event-level AUC table: all selected signals vs {background_label}\n")
        handle.write("=" * 150 + "\n")
        handle.write(f"Signal EOS base      : {normalise_eos_path(args.signal_eos_base)}\n")
        handle.write(f"Signal filters       : mDark={args.mdark:g}, yukawa={args.yukawa:g}, alpha={args.alpha}\n")
        handle.write(
            "Optional signal subsets: "
            f"mMed={args.mmeds if args.mmeds else 'all'}, "
            f"rinv={args.rinvs if args.rinvs else 'all'}\n"
        )
        handle.write(f"Selected signals     : {len(signal_samples)} points, {sum(len(s.files) for s in signal_samples)} ROOT files\n")
        handle.write(f"{background_label} EOS base      : {normalise_eos_path(background_base)}\n")
        handle.write(f"{background_label} samples       : {len(qcd_samples)} samples, {total_qcd_files} ROOT files\n")
        handle.write(f"Tree                  : {args.tree}\n")
        handle.write(f"ParticleNet branch    : {args.pnet_branch}\n")
        handle.write("Event score           : max(valid JetsAK8_pNetJetTaggerScore) per event\n")
        handle.write(f"Invalid-score cut     : values below {args.invalid_score_cut:g} are ignored\n")
        handle.write(
            "AUC cap               : "
            + ("all finite scores" if args.auc_cap < 0 else f"{args.auc_cap:,} finite events per class")
            + "\n\n"
        )

        handle.write(f"COMBINED-{background_label.upper()} AUC BY SIGNAL POINT\n")
        handle.write("-" * 150 + "\n")
        handle.write(
            f"{'mMed':>7} {'mDark':>7} {'rinv':>7} {'files':>7} {'Nsig(fin)':>12} {'Nqcd(fin)':>12} "
            f"{'Nsig(use)':>12} {'Nqcd(use)':>12} {'AUC raw':>10} {'AUC oriented':>13} {'direction':>10}  signal point\n"
        )
        handle.write("-" * 150 + "\n")
        for row in combined_rows:
            handle.write(
                f"{float(row['mMed']):7.0f} "
                f"{float(row['mDark']):7.3g} "
                f"{float(row['rinv']):7.3g} "
                f"{int(row['signal_files']):7d} "
                f"{int(row['n_signal_finite']):12,d} "
                f"{int(row['n_qcd_finite']):12,d} "
                f"{int(row['n_signal_used']):12,d} "
                f"{int(row['n_qcd_used']):12,d} "
                f"{format_auc(row['auc_raw_signal_high'], 10)} "
                f"{format_auc(row['auc_oriented'], 13)} "
                f"{str(row['orientation']):>10}  "
                f"{row['signal_name']}\n"
            )
        handle.write("\n")

        write_auc_grid(
            handle,
            title=(
                f"RAW AUC GRID: signal label=1, {background_label} label=0; "
                "larger pNet score assumed signal-like"
            ),
            combined_rows=combined_rows,
            value_field="auc_raw_signal_high",
        )
        write_auc_grid(
            handle,
            title="ORIENTED AUC GRID: max(raw AUC, 1 - raw AUC), for direction-independent separation only",
            combined_rows=combined_rows,
            value_field="auc_oriented",
        )

        handle.write(
            "DETAILED AUC TABLE: EVERY SIGNAL POINT VS EACH "
            f"{background_label} SAMPLE AND {combined_background_name}\n"
        )
        handle.write("-" * 150 + "\n")
        handle.write(
            f"{'mMed':>7} {'rinv':>7} {background_label + ' sample':<35} {'bkg files':>9} {'Nsig(fin)':>12} {'Nbkg(fin)':>12} "
            f"{'Nsig(use)':>12} {'Nqcd(use)':>12} {'AUC raw':>10} {'AUC oriented':>13} {'direction':>10}\n"
        )
        handle.write("-" * 150 + "\n")
        for row in detail_rows:
            handle.write(
                f"{float(row['mMed']):7.0f} "
                f"{float(row['rinv']):7.3g} "
                f"{str(row['qcd_sample']):<35} "
                f"{int(row['qcd_files']):9d} "
                f"{int(row['n_signal_finite']):12,d} "
                f"{int(row['n_qcd_finite']):12,d} "
                f"{int(row['n_signal_used']):12,d} "
                f"{int(row['n_qcd_used']):12,d} "
                f"{format_auc(row['auc_raw_signal_high'], 10)} "
                f"{format_auc(row['auc_oriented'], 13)} "
                f"{str(row['orientation']):>10}\n"
            )

        handle.write("\nNotes:\n")
        handle.write("  * Raw AUC is the physically conventional value for the selected pNet score direction.\n")
        handle.write("  * Oriented AUC is included only as a separation summary if the score convention is inverted.\n")
        handle.write(
            f"  * Per-sample {background_label} rows use the same score definition; "
            f"{combined_background_name} concatenates all available {background_label} samples.\n"
        )


# -----------------------------------------------------------------------------
# Main execution
# -----------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Compute ParticleNet AUC for all selected t-channel signals versus "
            "the selected 2018 background samples."
        ),
    )

    parser.add_argument("--signal-eos-base", default=DEFAULT_SIGNAL_EOS_BASE)
    parser.add_argument(
        "--background-kind",
        choices=("qcd", "ttjets", "st", "wjets", "zjets", "all"),
        default="qcd",
        help=(
            "Background sample family to compare against. "
            "'all' runs QCD, TTJets, ST, WJets, and ZJets into separate output subdirectories."
        ),
    )
    parser.add_argument("--qcd-eos-base", default=DEFAULT_QCD_EOS_BASE)
    parser.add_argument("--qcd-year", default="2018")
    parser.add_argument("--ttbar-eos-base", default=DEFAULT_TTBAR_EOS_BASE)
    parser.add_argument("--ttbar-year", default="2018")
    parser.add_argument(
        "--background-samples",
        nargs="+",
        default=None,
        metavar="SAMPLE",
        help="Optional subset for no-WNAE MC backgrounds such as ST, WJets, ZJets, or TTJets.",
    )

    parser.add_argument("--mdark", type=float, default=20.0)
    parser.add_argument("--yukawa", type=float, default=1.0)
    parser.add_argument(
        "--alpha",
        default="peak",
        help="Signal alpha token to retain. Pass 'any' to accept every alpha token.",
    )
    parser.add_argument(
        "--mmeds",
        nargs="+",
        type=float,
        default=None,
        metavar="MMED",
        help="Optional mediator-mass subset, e.g. --mmeds 1000 2000",
    )
    parser.add_argument(
        "--rinvs",
        nargs="+",
        type=float,
        default=None,
        metavar="RINV",
        help="Optional rinv subset, e.g. --rinvs 0.1 0.3",
    )
    parser.add_argument(
        "--qcd-bins",
        nargs="+",
        default=None,
        metavar="QCD_Pt_BIN",
        help="Optional QCD subset, e.g. QCD_Pt_600to800 QCD_Pt_800to1000",
    )
    parser.add_argument(
        "--ttbar-samples",
        nargs="+",
        default=None,
        metavar="TTJETS_SAMPLE",
        help="Optional TTJets subset, e.g. TTJets TTJets_HT-600to800",
    )

    parser.add_argument("--tree", default=DEFAULT_TREE)
    parser.add_argument("--pnet-branch", default=DEFAULT_PNET_BRANCH)
    parser.add_argument(
        "--invalid-score-cut",
        type=float,
        default=-9.0,
        help="Scores below this cut are treated as invalid sentinels.",
    )
    parser.add_argument("--step-size", type=int, default=50_000)
    parser.add_argument(
        "--max-events-per-sample",
        type=int,
        default=-1,
        help="Testing cap per signal and per background sample; -1 reads all available events.",
    )
    parser.add_argument(
        "--auc-cap",
        type=int,
        default=500_000,
        help="Max finite events used from each class in each AUC; -1 uses every finite event.",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--out-dir",
        default="/uscms/home/ashrivas/nobackup/Dark_Sector/t-channel_Analysis/pnet_auc_2018_GapVeto_PostQCDtrim",
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Do not read or write compressed per-sample score caches.",
    )
    return parser.parse_args()


def make_row(
    *,
    signal: SignalSample,
    qcd_name: str,
    qcd_files: int,
    auc_info: Dict[str, object],
) -> Dict[str, object]:
    metadata = signal.metadata
    return {
        "signal_name": signal.name,
        "signal_eos_dir": normalise_eos_path(signal.eos_dir),
        "signal_files": len(signal.files),
        "mMed": float(metadata["mMed"]),
        "mDark": float(metadata["mDark"]),
        "rinv": float(metadata["rinv"]),
        "alpha": str(metadata["alpha"]),
        "yukawa": float(metadata["yukawa"]),
        "qcd_sample": qcd_name,
        "qcd_files": qcd_files,
        **auc_info,
    }


def background_metadata(background_kind: str) -> Tuple[str, str]:
    if background_kind == "qcd":
        return "QCD", "ALL_QCD_COMBINED"
    config = NO_WNAE_BACKGROUND_CONFIG[background_kind]
    return str(config["label"]), str(config["combined_name"])


def output_path_for_background(args: argparse.Namespace, background_kind: str) -> str:
    if args.background_kind == "all":
        return os.path.join(args.out_dir, background_kind)
    return args.out_dir


def run_background(
    *,
    args: argparse.Namespace,
    signal_samples: Sequence[SignalSample],
    background_kind: str,
) -> None:
    out_dir = output_path_for_background(args, background_kind)
    os.makedirs(out_dir, exist_ok=True)
    background_label, combined_background_name = background_metadata(background_kind)

    print("=" * 100)
    print(f"ParticleNet AUC: all selected signals vs {background_label}")
    print("=" * 100)
    print(f"Signal EOS base     : {normalise_eos_path(args.signal_eos_base)}")
    print(f"Signal filters      : mDark={args.mdark:g}, yukawa={args.yukawa:g}, alpha={args.alpha}")
    if background_kind == "qcd":
        background_base = args.qcd_eos_base
        background_year = args.qcd_year
        print(f"QCD EOS base        : {normalise_eos_path(background_base)}")
    else:
        background_base = args.ttbar_eos_base
        background_year = args.ttbar_year
        print(f"{background_label} EOS base     : {normalise_eos_path(background_base)}")
    print(f"ParticleNet branch  : {args.pnet_branch}")
    print(f"AUC cap             : {'all finite events' if args.auc_cap < 0 else args.auc_cap}")

    if background_kind == "qcd":
        qcd_samples = discover_trimmed_qcd_samples(
            background_base,
            background_year,
            args.qcd_bins,
        )
    else:
        requested_samples = args.background_samples
        if background_kind == "ttjets" and requested_samples is None:
            requested_samples = args.ttbar_samples
        qcd_samples = discover_no_wnae_background_samples(
            background_base,
            background_year,
            background_kind,
            requested_samples,
        )
    if not qcd_samples:
        raise RuntimeError(f"No {background_label} ROOT files were discovered.")

    total_qcd_files = sum(len(sample.files) for sample in qcd_samples)
    print(
        f"\nStreaming {background_label} from EOS: "
        f"{len(qcd_samples)} samples, {total_qcd_files} ROOT files"
    )
    for sample in qcd_samples:
        print(f"  {background_year} {sample.name:<35} {len(sample.files):3d} files  {sample.eos_dir}")

    # Read/cache the background score arrays exactly once. They are re-used
    # against all signal points, avoiding repeated scans of the same ROOT files.
    qcd_scores_by_name: Dict[str, np.ndarray] = {}
    all_qcd_score_chunks: List[np.ndarray] = []
    for sample in qcd_samples:
        qcd_key = fingerprint(
            normalise_eos_path(sample.eos_dir),
            tuple(sample.files),
            args.tree,
            args.pnet_branch,
            args.invalid_score_cut,
            args.max_events_per_sample,
        )
        qcd_cache = os.path.join(
            out_dir,
            f"{background_kind}_pnet_scores_{safe_filename(sample.name)}_{qcd_key}.npz",
        )
        scores = load_or_build_scores(
            files=sample.files,
            label=sample.name,
            cache_path=qcd_cache,
            cache_key=qcd_key,
            args=args,
        )
        qcd_scores_by_name[sample.name] = scores
        all_qcd_score_chunks.append(scores)
        print(
            f"  {background_label} ready: {sample.name:<35} "
            f"events={len(scores):,}, finite={len(finite_scores(scores)):,}"
        )

    combined_qcd_scores = (
        np.concatenate(all_qcd_score_chunks)
        if all_qcd_score_chunks
        else np.empty(0, dtype=np.float64)
    )
    if len(finite_scores(combined_qcd_scores)) == 0:
        raise RuntimeError(f"No finite {background_label} ParticleNet scores were loaded.")

    detail_rows: List[Dict[str, object]] = []
    combined_rows: List[Dict[str, object]] = []

    print("\nComputing AUCs for every selected signal point")
    for signal_index, signal in enumerate(signal_samples, start=1):
        print("-" * 100)
        print(
            f"[{signal_index:02d}/{len(signal_samples):02d}] "
            f"mMed={float(signal.metadata['mMed']):.0f}, "
            f"rinv={float(signal.metadata['rinv']):g}: {signal.name}"
        )

        signal_key = fingerprint(
            normalise_eos_path(signal.eos_dir),
            tuple(signal.files),
            args.tree,
            args.pnet_branch,
            args.invalid_score_cut,
            args.max_events_per_sample,
        )
        signal_cache = os.path.join(
            out_dir,
            f"signal_pnet_scores_{safe_filename(signal.name)}_{signal_key}.npz",
        )
        signal_scores = load_or_build_scores(
            files=signal.files,
            label=signal.name,
            cache_path=signal_cache,
            cache_key=signal_key,
            args=args,
        )
        n_signal_finite = len(finite_scores(signal_scores))
        print(f"  Signal ready: events={len(signal_scores):,}, finite={n_signal_finite:,}")

        # Keep this point in the output even if a particular signal has no valid
        # scores. Its AUC fields will be nan, making the problem visible.
        combined_info = compute_auc(
            signal_scores,
            combined_qcd_scores,
            args.auc_cap,
            deterministic_rng(args.seed, signal.name + "|" + combined_background_name),
        )
        combined_row = make_row(
            signal=signal,
            qcd_name=combined_background_name,
            qcd_files=total_qcd_files,
            auc_info=combined_info,
        )
        combined_rows.append(combined_row)
        detail_rows.append(combined_row)
        print(
            f"  {combined_background_name:<35} "
            f"Nsig={int(combined_info['n_signal_finite']):,} "
            f"Nqcd={int(combined_info['n_qcd_finite']):,} "
            f"AUC(raw/oriented)={format_auc(combined_info['auc_raw_signal_high'])}/"
            f"{format_auc(combined_info['auc_oriented'])}"
        )

        for qcd_sample in qcd_samples:
            auc_info = compute_auc(
                signal_scores,
                qcd_scores_by_name[qcd_sample.name],
                args.auc_cap,
                deterministic_rng(args.seed, signal.name + "|" + qcd_sample.name),
            )
            row = make_row(
                signal=signal,
                qcd_name=qcd_sample.name,
                qcd_files=len(qcd_sample.files),
                auc_info=auc_info,
            )
            detail_rows.append(row)
            print(
                f"  {qcd_sample.name:<22} "
                f"AUC(raw/oriented)={format_auc(auc_info['auc_raw_signal_high'])}/"
                f"{format_auc(auc_info['auc_oriented'])}"
            )

    # Preserve a compact, useful order in every output file.
    qcd_rank = {combined_background_name: -1}
    qcd_rank.update({sample.name: index for index, sample in enumerate(qcd_samples)})
    detail_rows.sort(
        key=lambda row: (
            float(row["mMed"]),
            float(row["rinv"]),
            qcd_rank.get(str(row["qcd_sample"]), 10**6),
        )
    )
    combined_rows.sort(key=lambda row: (float(row["mMed"]), float(row["rinv"])))

    detailed_csv_path = os.path.join(out_dir, "pnet_auc_table.csv")
    combined_csv_path = os.path.join(out_dir, f"pnet_auc_combined_{background_kind}.csv")
    text_path = os.path.join(out_dir, "pnet_auc_table.txt")
    write_csv(detailed_csv_path, DETAIL_FIELDS, detail_rows)
    write_csv(combined_csv_path, COMBINED_FIELDS, combined_rows)
    write_text_table(
        path=text_path,
        args=args,
        background_kind=background_kind,
        background_label=background_label,
        combined_background_name=combined_background_name,
        background_base=background_base,
        signal_samples=signal_samples,
        qcd_samples=qcd_samples,
        combined_rows=combined_rows,
        detail_rows=detail_rows,
    )

    print("\nDone.")
    print(f"  Detailed CSV : {detailed_csv_path}")
    print(f"  Combined CSV : {combined_csv_path}")
    print(f"  TXT table    : {text_path}")


def main() -> None:
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    signal_samples = discover_signal_samples(
        args.signal_eos_base,
        m_dark=args.mdark,
        yukawa=args.yukawa,
        alpha=args.alpha,
        requested_mmeds=args.mmeds,
        requested_rinvs=args.rinvs,
    )
    if not signal_samples:
        raise RuntimeError("No signal points survived the requested filters.")

    print(f"\nDiscovered signals: {len(signal_samples)} selected points")
    for signal in signal_samples:
        print(
            f"  mMed={float(signal.metadata['mMed']):6.0f}  "
            f"rinv={float(signal.metadata['rinv']):5.2f}  "
            f"files={len(signal.files):3d}  {signal.name}"
        )

    background_kinds = (
        ("qcd", "ttjets", "st", "wjets", "zjets")
        if args.background_kind == "all"
        else (args.background_kind,)
    )
    for background_kind in background_kinds:
        run_background(
            args=args,
            signal_samples=signal_samples,
            background_kind=background_kind,
        )


if __name__ == "__main__":
    main()
