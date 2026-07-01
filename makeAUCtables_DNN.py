#!/usr/bin/env python3
import uproot
import numpy as np
from sklearn.metrics import roc_auc_score
import os
import re
import glob
import subprocess
from collections import defaultdict

# ---------------------------------------------------------
# Paths
# ---------------------------------------------------------



def eos_recursive_ls(eos_dir):
    """List EOS files recursively using xrdfs."""
    result = subprocess.run(
        ["xrdfs", EOS_HOST, "ls", "-R", eos_dir],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"Could not list EOS directory:\n{eos_dir}\n\n{result.stderr}"
        )

    return [
        line.strip()
        for line in result.stdout.splitlines()
        if line.strip()
    ]


def eos_to_xrootd(eos_path):
    """Convert /store/... into an uproot-readable root:// URL."""
    return f"{EOS_HOST}//{eos_path.lstrip('/')}"


def get_arrays(filepath):
    """Extract scores and weights from the ABCD_scores TTree."""
    with uproot.open(filepath) as f:
        if "ABCD_scores" not in f:
            available = [key.split(";")[0] for key in f.keys()]
            raise KeyError(
                f"No ABCD_scores tree in {filepath}\n"
                f"Available objects: {available}"
            )

        tree = f["ABCD_scores"]

        if "scores" not in tree.keys() or "weights" not in tree.keys():
            raise KeyError(
                f"Missing scores/weights branches in {filepath}\n"
                f"Branches: {list(tree.keys())}"
            )

        return (
            tree["scores"].array(library="np"),
            tree["weights"].array(library="np"),
        )


#SIG_DIR = "/uscms/home/nparmar/nobackup/SVJ/eventLevelTagger/output/dataset11/training_set2/application_results/sdt_allBkgs_disco_0p001_closure_0p06_damp_1_net_64_32_16_8_trainingAllYears_met250_dnn85_eval2018/v_5/ML_fit_results_grid_optimization_binned/scores_signals_skims_extended_distortion_checks_alternative_v2/2018/nominal"
SIG_DIR = "/uscms/home/nparmar/nobackup/SVJ/eventLevelTagger/output/dataset11/training_set2/application_results/sdt_allBkgs_disco_0p001_closure_0p06_damp_1_net_64_32_16_8_trainingAllYears_met250_dnn85_eval2018/v_5/ML_fit_results_grid_optimization_binned/scores_signals_GapJetVeto/2018/nominal"

#SIG_DIR = "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/signals/2018/nominal/"
#BKG_FILE = "/uscms/home/nparmar/nobackup/SVJ/eventLevelTagger/output/dataset11/training_set2/application_results/sdt_allBkgs_disco_0p001_closure_0p06_damp_1_net_64_32_16_8_trainingAllYears_met250_dnn85_eval2018/v_5/ML_fit_results_grid_optimization_binned/scores_backgrounds_scaler_keane_lumiscaled_2018metbinning_filteredQCD/2018/nominal/ABCD_scores_QCD.root"
BKG_DIR = (
    "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_gapJetVeto/"
    "DNN_scores_for_background_with_PNET/scores_backgrounds_GapJetVeto/2018/nominal"
)
BACKGROUND_FILES = {
    "QCD": f"{BKG_DIR}/ABCD_scores_QCD.root",
    "TTJets": f"{BKG_DIR}/ABCD_scores_TTJets.root",
    "ST": f"{BKG_DIR}/ABCD_scores_ST.root",
    "WJetsToLNu": f"{BKG_DIR}/ABCD_scores_WJetsToLNu.root",
    "ZJetsToNuNu": f"{BKG_DIR}/ABCD_scores_ZJetsToNuNu.root",
}


def print_auc_grid(background_name, results_dict):
    """Print a 2D AUC grid for one background family."""
    mMed_vals = sorted(list(set(k[0] for k in results_dict.keys())))
    rinv_vals = sorted(list(set(k[1] for k in results_dict.keys())))

    print("\n" + "=" * 95)
    print(f"DNN AUC GRID vs {background_name} (mDark=20, alpha=peak, yukawa=1)")
    print("=" * 95)

    header = f"{'mMed / rinv':>12} | " + " | ".join([f"{r:>5.2f}" for r in rinv_vals])
    print(header)
    print("-" * len(header))

    for m in mMed_vals:
        row_str = f"{m:>12} | "
        for r in rinv_vals:
            auc = results_dict.get((m, r), float("nan"))
            if np.isnan(auc):
                row_str += "  ---  | "
            else:
                row_str += f"{auc:>5.3f} | "
        print(row_str)
    print("=" * 95 + "\n")

def main():
    background_arrays = {}
    for background_name, background_file in BACKGROUND_FILES.items():
        print(f"Loading {background_name} background from:\n  {background_file}")
        try:
            b_sc, b_wt = get_arrays(background_file)
        except Exception as e:
            print(f"  -> Error loading {background_name}: {e}")
            continue
        background_arrays[background_name] = (b_sc, b_wt)
        print(f"  -> Loaded {len(b_sc):,} background events.\n")

    if not background_arrays:
        print("No background score files could be loaded.")
        return

    # Broaden search to catch all rinv values for yukawa-1
    search_pattern = os.path.join(SIG_DIR, "ABCD_scores_nominal_t-channel_mMed-*_mDark-20_rinv-*_alpha-peak_yukawa-1_*.root")
    all_files = glob.glob(search_pattern)

    # Filter to keep ONLY the nominal files (accounting for the inconsistent naming convention)
    sig_files = [f for f in all_files if f.endswith("_nominal.root") or f.endswith("_lundWeightNom.root")]

    if not sig_files:
        print("No nominal files found matching the pattern.")
        return

    results_by_background = {name: {} for name in background_arrays}

    print(f"Scanning {len(sig_files)} signal points and computing AUCs...")
    for sf in sig_files:
        # Extract parameters using regex
        fname = os.path.basename(sf)
        m = re.search(r"mMed-(\d+)_mDark-(\d+)_rinv-([\dp]+)_alpha-([\w]+)", fname)
        if not m: continue

        mMed = int(m.group(1))
        rinv = float(m.group(3).replace("p", "."))

        # Load signal arrays
        try:
            s_sc, s_wt = get_arrays(sf)
        except Exception as e:
            print(f"Error reading {fname}: {e}")
            continue

        if len(s_sc) == 0:
            continue

        for background_name, (b_sc, b_wt) in background_arrays.items():
            scores = np.concatenate([s_sc, b_sc])
            wts = np.concatenate([s_wt, b_wt])
            labels = np.concatenate([np.ones(len(s_sc)), np.zeros(len(b_sc))])

            auc = roc_auc_score(labels, scores, sample_weight=wts)
            auc = auc if auc > 0.5 else 1.0 - auc

            results_by_background[background_name][(mMed, rinv)] = auc

    for background_name, results_dict in results_by_background.items():
        print_auc_grid(background_name, results_dict)

if __name__ == "__main__":
    main()
