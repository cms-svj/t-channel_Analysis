# #!/usr/bin/env python3
# """
# pnet_auc_grid_eos.py

# Compute ParticleNet AUC for each signal point vs 2018 QCD background.
# Scans directly from an EOS directory using xrdfs.
# """

# import os, re, sys, subprocess, argparse
# from typing import Dict, List, Tuple, Optional

# import numpy as np
# import uproot
# import awkward as ak
# import matplotlib
# matplotlib.use("Agg")
# import matplotlib.pyplot as plt
# import matplotlib.colors as mcolors
# import mplhep as hep
# from sklearn.metrics import roc_auc_score, roc_curve

# plt.style.use(hep.style.CMS)

# # ─────────────────────────────────────────────────────────────────────────────
# # Constants
# # ─────────────────────────────────────────────────────────────────────────────

# PNET_BRANCH = "JetsAK8_pNetJetTaggerScore"
# TREE        = "Events"

# QCD_FILES_2018: Dict[str, List[str]] = {

#     "QCD_Pt_470to600": [
#             "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_470to600/part-0.root",
#             "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_470to600/part-1.root",
#         ],

#     "QCD_Pt_600to800": [
#         "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-0.root",
#         "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-1.root",
#         "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-2.root",
#         "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-3.root",
#         "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-4.root",
#         "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-5.root",
#         "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-6.root",
#         "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-7.root",
#         "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims_wnae/2018/t_channel_pre_selection/nominal//QCD_Pt_600to800/part-8.root",
#     ],
# }

# REF = dict(mMed=2000, rinv=0.3, mDark=20)


# # ─────────────────────────────────────────────────────────────────────────────
# # EOS Discovery & Parsing
# # ─────────────────────────────────────────────────────────────────────────────

# def ensure_dir(p: str) -> None:
#     os.makedirs(p, exist_ok=True)

# def parse_signal_dir(dir_path: str) -> Optional[Dict]:
#     """Parse parameters from the EOS directory name."""
#     name = os.path.basename(dir_path)
#     m = re.search(
#         r"mMed-(\d+)_mDark-(\d+)_rinv-([\dp]+)_alpha-([\w]+)_yukawa-([\dp]+)",
#         name,
#     )
#     if not m:
#         return None
#     return dict(
#         mMed   = int(m.group(1)),
#         mDark  = int(m.group(2)),
#         rinv   = float(m.group(3).replace("p", ".")),
#         alpha  = m.group(4),
#         yukawa = m.group(5),
#         eos_path = dir_path,
#         label  = (f"mMed={m.group(1)}, mDark={m.group(2)}, "
#                   f"rinv={m.group(3).replace('p','.')}, α={m.group(4)}"),
#     )

# def files_from_eos_dir(dir_path: str) -> List[str]:
#     """Uses xrdfs to list all .root files inside a specific EOS directory."""
#     cmd = f"xrdfs root://cmseos.fnal.gov ls {dir_path}"
#     try:
#         out = subprocess.check_output(cmd, shell=True, text=True)
#         lines = out.strip().split('\n')
#         # Prepend xrootd prefix for uproot
#         return [f"root://cmseos.fnal.gov/{l}" for l in lines if l.endswith(".root")]
#     except Exception as e:
#         print(f"  [ERROR] Failed to list files in {dir_path}: {e}")
#         return []


# # ─────────────────────────────────────────────────────────────────────────────
# # Score loading
# # ─────────────────────────────────────────────────────────────────────────────

# def find_tree(fin):
#     candidates = ["Events", "TreeMaker2/PreSelection", "PreSelection", "tree"]
#     for key in candidates:
#         if key in fin: return key
#     for k, v in fin.items():
#         if isinstance(v, uproot.behaviors.TTree.TTree): return k
#     return None

# def load_max_jet_scores(
#     files: List[str],
#     branch: str    = PNET_BRANCH,
#     tree: str      = TREE,
#     step_size: int = 50_000,
#     max_events: int = -1,
# ) -> np.ndarray:
#     chunks = []
#     total  = 0
#     for fpath in files:
#         if max_events > 0 and total >= max_events: break
#         try:
#             with uproot.open(fpath) as fin:
#                 use_tree = tree if tree in fin else find_tree(fin)
#                 if not use_tree: continue
                
#                 t = fin[use_tree]
#                 if branch not in t.keys(): continue
#                 for arrays in t.iterate([branch], step_size=step_size, library="ak"):
#                     raw = arrays[branch]
#                     raw = ak.where(raw < -9.0, np.nan, raw)
#                     mx  = ak.to_numpy(ak.fill_none(ak.max(raw, axis=1), np.nan))
#                     mx  = mx[np.isfinite(mx)]
#                     chunks.append(mx)
#                     total += len(mx)
#                     if max_events > 0 and total >= max_events: break
#         except Exception as e:
#             pass
#     return np.concatenate(chunks) if chunks else np.array([], dtype=np.float64)


# # ─────────────────────────────────────────────────────────────────────────────
# # AUC + ROC
# # ─────────────────────────────────────────────────────────────────────────────

# def compute_auc_roc(sig: np.ndarray, bkg: np.ndarray, cap: int = 500_000) -> Tuple[float, np.ndarray, np.ndarray]:
#     rng = np.random.default_rng(42)
#     s = rng.choice(sig, min(len(sig), cap), replace=False)
#     b = rng.choice(bkg, min(len(bkg), cap), replace=False)
#     scores = np.concatenate([s, b])
#     labels = np.concatenate([np.ones(len(s)), np.zeros(len(b))])
#     try:
#         auc = roc_auc_score(labels, scores)
#         fpr, tpr, _ = roc_curve(labels, scores)
#         if auc < 0.5:
#             auc = 1.0 - auc
#             fpr, tpr = tpr, fpr
#         return auc, fpr, tpr
#     except Exception:
#         return float("nan"), np.array([0., 1.]), np.array([0., 1.])


# # ─────────────────────────────────────────────────────────────────────────────
# # Plot helpers
# # ─────────────────────────────────────────────────────────────────────────────

# def _cms_label(ax):
#     hep.cms.label("Simulation Preliminary", data=False, ax=ax, loc=0, fontsize=12)
#     #ax.text(1.0, 1.01, "(13 TeV)", transform=ax.transAxes, ha="right", va="bottom", fontsize=11)

# def save(fig, prefix: str) -> None:
#     for ext in ("pdf", "png"):
#         path = f"{prefix}.{ext}"
#         fig.savefig(path, bbox_inches="tight", dpi=150 if ext == "png" else None)
#         print(f"  -> {path}")
#     plt.close(fig)


# # ─────────────────────────────────────────────────────────────────────────────
# # Plotting Functions (Grid, ROC, Dist)
# # ────--─────────────────────────────────────────────────────────────────────────

# def plot_auc_grid(results: List[Dict], out_prefix: str) -> None:
#     pts = [r for r in results if r["mDark"] == 20 and np.isfinite(r["auc"])]
#     if not pts: return

#     mMed_vals = sorted(set(r["mMed"] for r in pts))
#     rinv_vals = sorted(set(r["rinv"] for r in pts))

#     Z = np.full((len(mMed_vals), len(rinv_vals)), np.nan)
#     for r in pts:
#         i, j = mMed_vals.index(r["mMed"]), rinv_vals.index(r["rinv"])
#         Z[i, j] = r["auc"]

#     nrows, ncols = len(mMed_vals), len(rinv_vals)
#     fig, ax = plt.subplots(figsize=(max(6, ncols * 1.5 + 2), max(5, nrows * 0.9 + 2)))
#     im = ax.imshow(Z, origin="lower", aspect="auto", cmap="viridis", vmin=0.5, vmax=1.0, extent=[-0.5, ncols - 0.5, -0.5, nrows - 0.5])

#     for i in range(nrows):
#         for j in range(ncols):
#             v = Z[i, j]
#             if np.isfinite(v):
#                 txt_c = "white" #if (v - 0.5) / 0.5 < 0.55 else "black"
#                 ax.text(j, i, f"{v:.3f}", ha="center", va="center", fontsize=10, color=txt_c, fontweight="bold")
#             else:
#                 ax.text(j, i, "—", ha="center", va="center", fontsize=10, color="grey")

#     ax.set_xticks(range(ncols)); ax.set_xticklabels([f"{r:.1f}" for r in rinv_vals], fontsize=11)
#     ax.set_yticks(range(nrows)); ax.set_yticklabels([str(m) for m in mMed_vals], fontsize=11)
#     ax.set_xlabel(r"$r_{\mathrm{inv}}$", fontsize=13)
#     ax.set_ylabel(r"$m_{\Phi}\ [\mathrm{GeV}]$", fontsize=13)

#     cb = fig.colorbar(im, ax=ax, pad=0.02)
#     cb.set_label("AUC", fontsize=12)
#     _cms_label(ax)
#     #ax.set_title(r"ParticleNet SVJ Tagger AUC  |  $m_{\rm dark}=20$ GeV, $\alpha={\rm peak}$  |  2018 QCD bkg", fontsize=10, pad=36)
#     plt.tight_layout()
#     save(fig, out_prefix)

# def plot_roc_curves(results: List[Dict], out_prefix: str) -> None:
#     finite = [r for r in results if np.isfinite(r["auc"])]
#     if not finite: return

#     fig, ax = plt.subplots(figsize=(9, 8))
#     mMed_u = sorted(set(r["mMed"] for r in finite))
#     cmap = plt.cm.plasma
#     m_norm = mcolors.Normalize(vmin=min(mMed_u), vmax=max(mMed_u))
#     lstyles = ["-", "--", "-.", ":", (0,(3,1,1,1)), (0,(5,2))]
#     style_keys = sorted(set((r["mDark"], r["rinv"]) for r in finite))

#     for r in sorted(finite, key=lambda x: (x["mMed"], x["rinv"], x["mDark"])):
#         colour = cmap(m_norm(r["mMed"]))
#         ls_idx = style_keys.index((r["mDark"], r["rinv"])) % len(lstyles)
#         lbl = f"$m_{{\\Phi}}={r['mMed']}$, $r_{{\\rm inv}}={r['rinv']:.1f}$, $m_{{\\rm dark}}={r['mDark']}$ (AUC={r['auc']:.3f})"
#         ax.plot(r["fpr"], r["tpr"], color=colour, linestyle=lstyles[ls_idx], linewidth=1.4, label=lbl)

#     ax.plot([0,1],[0,1], "k--", linewidth=0.8, label="Random")
#     ax.set_xlabel("False Positive Rate  (QCD efficiency)", fontsize=13)
#     ax.set_ylabel("True Positive Rate  (Signal efficiency)", fontsize=13)
#     ax.set_xlim(0, 1); ax.set_ylim(0, 1)

#     sm = plt.cm.ScalarMappable(cmap=cmap, norm=m_norm)
#     sm.set_array([])
#     cb = fig.colorbar(sm, ax=ax, pad=0.02)
#     cb.set_label(r"$m_{\Phi}$ [GeV]", fontsize=12)

#     ax.legend(fontsize=6.5, loc="lower right", ncol=2, framealpha=0.7)
#     _cms_label(ax)
#     ax.set_title("ROC Curves — ParticleNet SVJ Tagger  |  2018 QCD bkg", fontsize=11, pad=36)
#     plt.tight_layout()
#     save(fig, out_prefix)

# def plot_score_dist(bkg_scores: np.ndarray, sig_scores: np.ndarray, ref: Dict, out_prefix: str, nbins: int = 50) -> None:
#     bins = np.linspace(0, 1, nbins + 1)
#     centres = 0.5 * (bins[:-1] + bins[1:])
#     width = bins[1] - bins[0]

#     def norm_hist(arr):
#         h, _ = np.histogram(arr, bins=bins)
#         s = h.sum() * width
#         return h / s if s > 0 else h.astype(float)

#     fig, ax = plt.subplots(figsize=(9, 7))
#     bkg_h = norm_hist(bkg_scores)
#     ax.bar(centres, bkg_h, width=width, alpha=0.45, color="C0", label=f"QCD background  (N={len(bkg_scores):,})", linewidth=0)
#     ax.step(centres, bkg_h, where="mid", color="C0", linewidth=2.0)

#     sig_h = norm_hist(sig_scores)
#     sig_label = f"Signal  $m_{{\\Phi}}={ref['mMed']}$ GeV, $r_{{\\rm inv}}={ref['rinv']}$, $m_{{\\rm dark}}={ref['mDark']}$ GeV\n(N={len(sig_scores):,})"
#     ax.bar(centres, sig_h, width=width, alpha=0.45, color="C3", label=sig_label, linewidth=0)
#     ax.step(centres, sig_h, where="mid", color="C3", linewidth=2.0)

#     ax.set_xlabel("ParticleNet SVJ Tagger Score  (max over jets per event)", fontsize=13)
#     ax.set_ylabel("Normalised events / bin", fontsize=13)
#     #ax.set_yscale("log")
#     ax.set_xlim(0, 1); ax.set_ylim(bottom=1e-4)
#     ax.legend(fontsize=11, framealpha=0.85)
#     _cms_label(ax)
#     ax.set_title("Score Distribution  |  2018 QCD vs Signal", fontsize=11, pad=36)
#     plt.tight_layout()
#     save(fig, out_prefix)


# # ─────────────────────────────────────────────────────────────────────────────
# # Main
# # ─────────────────────────────────────────────────────────────────────────────

# def main():
#     ap = argparse.ArgumentParser()
#     ap.add_argument("--signal-eos-dir", required=True, help="Base EOS path containing the signal subdirectories")
#     ap.add_argument("--pnet-branch", default=PNET_BRANCH)
#     ap.add_argument("--tree",        default=TREE)
#     ap.add_argument("--out-dir",     default="./pnet_auc_plots")
#     ap.add_argument("--step-size",   type=int, default=50_000)
#     ap.add_argument("--max-events",  type=int, default=-1)
#     ap.add_argument("--qcd-cache",   default="./qcd_pnet_scores.npy")
#     ap.add_argument("--no-cache",    action="store_true")
#     args = ap.parse_args()

#     ensure_dir(args.out_dir)

#     # 1. SCAN EOS DIRECTORY
#     print(f"\nScanning EOS directory: {args.signal_eos_dir}")
#     cmd = f"xrdfs root://cmseos.fnal.gov ls {args.signal_eos_dir}"
#     try:
#         out = subprocess.check_output(cmd, shell=True, text=True)
#         subdirs = [line for line in out.strip().split('\n') if "t-channel_mMed" in line]
#     except Exception as e:
#         print(f"[ERROR] Failed to run xrdfs on {args.signal_eos_dir}: {e}")
#         sys.exit(1)

#     signal_points = []
#     for d in subdirs:
#         meta = parse_signal_dir(d)
#         if meta: signal_points.append(meta)

#     if not signal_points:
#         print("[ERROR] No signal directories could be parsed.")
#         sys.exit(1)

#     print(f"Found {len(signal_points)} signal points. Examples:")
#     for sp in sorted(signal_points, key=lambda x: (x["mMed"], x["rinv"], x["mDark"]))[:5]:
#         print(f"  mMed={sp['mMed']:5d}  rinv={sp['rinv']:.1f}  mDark={sp['mDark']:4d}")

#     # 2. LOAD QCD BACKGROUND
#     cache = args.qcd_cache
#     if not args.no_cache and os.path.exists(cache):
#         print(f"\nLoading QCD scores from cache: {cache}")
#         bkg_scores = np.load(cache)
#     else:
#         print("\nStreaming QCD background from EOS …")
#         all_qcd = [f for fs in QCD_FILES_2018.values() for f in fs]
#         bkg_scores = load_max_jet_scores(all_qcd, branch=args.pnet_branch, tree=args.tree, step_size=args.step_size, max_events=args.max_events)
#         if len(bkg_scores) > 0 and not args.no_cache:
#             np.save(cache, bkg_scores)

#     if len(bkg_scores) == 0:
#         print("[ERROR] No QCD scores loaded.")
#         sys.exit(1)

#     # 3. COMPUTE AUC
#     print(f"\nComputing AUC for {len(signal_points)} points …")
#     results = []
#     ref_sig_scores = None

#     for sp in sorted(signal_points, key=lambda x: (x["mMed"], x["rinv"], x["mDark"])):
#         root_files = files_from_eos_dir(sp["eos_path"])
#         if not root_files:
#             continue
        
#         sig_scores = load_max_jet_scores(root_files, branch=args.pnet_branch, tree=args.tree, step_size=args.step_size, max_events=args.max_events)
        
#         if len(sig_scores) == 0:
#             results.append({**sp, "auc": float("nan"), "fpr": None, "tpr": None})
#             continue

#         auc, fpr, tpr = compute_auc_roc(sig_scores, bkg_scores)
#         results.append({**sp, "auc": auc, "fpr": fpr, "tpr": tpr})
#         print(f"  [{sp['label']}]  N_sig={len(sig_scores):>7,}  AUC={auc:.4f}")

#         if (sp["mMed"] == REF["mMed"] and abs(sp["rinv"] - REF["rinv"]) < 1e-4 and sp["mDark"] == REF["mDark"]):
#             ref_sig_scores = sig_scores

#     # 4. EXPORT AND PLOT
#     table_path = os.path.join(args.out_dir, "pnet_auc_table.txt")
#     with open(table_path, "w") as f:
#         f.write(f"{'mMed':>8}  {'mDark':>6}  {'rinv':>6}  {'alpha':>6}  {'AUC':>8}\n")
#         f.write("─" * 46 + "\n")
#         for r in sorted(results, key=lambda x: (x["mMed"], x["rinv"], x["mDark"])):
#             f.write(f"{r['mMed']:>8}  {r['mDark']:>6}  {r['rinv']:>6.1f}  {r['alpha']:>6}  {r['auc']:>8.4f}\n")

#     print("\nMaking plots …")
#     plot_auc_grid(results, out_prefix=os.path.join(args.out_dir, "pnet_auc_grid"))
#     plot_roc_curves(results, out_prefix=os.path.join(args.out_dir, "pnet_roc_curves"))

#     if ref_sig_scores is not None and len(ref_sig_scores) > 0:
#         plot_score_dist(bkg_scores, ref_sig_scores, ref=REF, out_prefix=os.path.join(args.out_dir, "pnet_score_dist"))

#     print(f"\n=== Done.  All outputs in: {args.out_dir}/ ===")



# if __name__ == "__main__":
#     main()




#!/usr/bin/env python3
"""
plot_dual_auc.py

Produce two separate AUC grid plots:
  1. ParticleNet SVJ Tagger  (pnet_auc_grid.png/pdf)
  2. WNAE Tagger             (wnae_auc_grid.png/pdf)
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import matplotlib.patheffects as pe
plt.style.use(hep.style.CMS)

# ─────────────────────────────────────────────────────────────────────────────
# PNet AUC values
# rows = mMed 500→4000, cols = rinv 0.0→0.9
# ─────────────────────────────────────────────────────────────────────────────

MMED_PNET = [500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000, 3500, 4000]
#MMED_PNET = [500, 600, 700, 800, 1000, 1500, 2000, 3000, 4000]
RINV_PNET = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

PNET_RAW = [
    [0.892, 0.914, 0.886, 0.819, 0.833, 0.802, 0.753, 0.698, 0.626, 0.520],  # 500
    [0.914, 0.907, 0.894, 0.832, 0.840, 0.811, 0.760, 0.706, 0.636, 0.547],  # 600
    [0.953, 0.926, 0.888, 0.840, 0.854, 0.813, 0.781, 0.724, 0.656, 0.546],  # 700
    [0.937, 0.920, 0.909, 0.850, 0.859, 0.823, 0.785, 0.733, 0.672, 0.563],  # 800
    [0.942, 0.936, 0.909, 0.856, 0.872, 0.832, 0.799, 0.740, 0.687, 0.580],  # 900
    [0.945, 0.931, 0.912, 0.860, 0.872, 0.835, 0.804, 0.741, 0.667, 0.574],  # 1000
    [0.944, 0.945, 0.928, 0.889, 0.878, 0.845, 0.804, 0.757, 0.681, 0.601],  # 1500
    [0.944, 0.953, 0.935, 0.915, 0.901, 0.866, 0.827, 0.775, 0.709, 0.594],  # 2000
    [0.950, 0.957, 0.943, 0.928, 0.911, 0.877, 0.828, 0.779, 0.707, 0.625],  # 2500
    [0.967, 0.963, 0.944, 0.931, 0.917, 0.888, 0.864, 0.804, 0.723, 0.638],  # 3000
    [0.954, 0.968, 0.952, 0.931, 0.917, 0.883, 0.850, 0.818, 0.717, 0.632],  # 3500
    [0.970, 0.958, 0.944, 0.932, 0.915, 0.896, 0.857, 0.814, 0.737, 0.637],  # 4000
]

# ─────────────────────────────────────────────────────────────────────────────
# WNAE AUC values
# rows = mMed 500→4000, cols = rinv 0.0→0.9
# ─────────────────────────────────────────────────────────────────────────────

#MMED_WNAE = [500, 600, 700, 800, 1000, 1500, 2000, 3000, 4000]
MMED_WNAE = [500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000, 3500, 4000]
RINV_WNAE = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

WNAE_RAW = [
    [0.886, 0.786, 0.764, 0.794, 0.744, 0.708, 0.708, 0.654, 0.580, 0.522],  # 500
    [0.780, 0.774, 0.772, 0.777, 0.734, 0.722, 0.710, 0.667, 0.634, 0.531],  # 600
    [0.643, 0.816, 0.791, 0.754, 0.762, 0.738, 0.701, 0.652, 0.610, 0.557],  # 700
    [0.804, 0.780, 0.802, 0.782, 0.756, 0.723, 0.707, 0.665, 0.627, 0.593],  # 800
    [0.821, 0.784, 0.793, 0.759, 0.769, 0.733, 0.730, 0.666, 0.633, 0.584],  # 900
    [0.828, 0.816, 0.783, 0.772, 0.757, 0.741, 0.727, 0.662, 0.649, 0.581],  # 1000
    [0.756, 0.793, 0.789, 0.756, 0.748, 0.741, 0.723, 0.681, 0.612, 0.544],  # 1500
    [0.759, 0.772, 0.765, 0.763, 0.741, 0.735, 0.733, 0.711, 0.674, 0.537],  # 2000
    [0.743, 0.768, 0.760, 0.767, 0.775, 0.749, 0.735, 0.689, 0.662, 0.561],  # 2500
    [0.734, 0.774, 0.767, 0.754, 0.766, 0.752, 0.727, 0.678, 0.630, 0.576],  # 3000
    [0.692, 0.783, 0.784, 0.754, 0.762, 0.751, 0.733, 0.673, 0.661, 0.580],  # 3500
    [0.672, 0.782, 0.772, 0.778, 0.774, 0.736, 0.730, 0.683, 0.692, 0.622],  # 4000
]




# ─────────────────────────────────────────────────────────────────────────────
# Third tagger AUC values  (rename tagger_name in main() as needed)
# rows = mMed 500→4000, cols = rinv 0.0→1.0  (11 columns)
# ─────────────────────────────────────────────────────────────────────────────
 
MMED_DNN   = [500, 600, 700, 800, 1000, 1500, 2000, 3000, 4000]
RINV_DNN   = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

 
DNN_RAW = [
    [0.792, 0.821, 0.818, 0.799, 0.790, 0.754, 0.743, 0.684, 0.659, 0.586, 0.521],  # 500
    [0.763, 0.839, 0.817, 0.809, 0.791, 0.788, 0.760, 0.718, 0.666, 0.609, 0.539],  # 600
    [0.820, 0.839, 0.830, 0.816, 0.811, 0.782, 0.762, 0.738, 0.697, 0.655, 0.548],  # 700
    [0.871, 0.838, 0.834, 0.830, 0.817, 0.787, 0.776, 0.755, 0.702, 0.644, 0.577],  # 800
    [0.814, 0.845, 0.833, 0.826, 0.816, 0.798, 0.774, 0.740, 0.697, 0.661, 0.555],  # 1000
    [0.832, 0.821, 0.812, 0.803, 0.783, 0.757, 0.724, 0.684, 0.640, 0.608, 0.554],  # 1500
    [0.822, 0.802, 0.791, 0.775, 0.763, 0.736, 0.715, 0.676, 0.620, 0.584, 0.535],  # 2000
    [0.771, 0.782, 0.760, 0.738, 0.732, 0.720, 0.683, 0.662, 0.607, 0.571, 0.550],  # 3000
    [0.736, 0.754, 0.747, 0.729, 0.718, 0.698, 0.662, 0.634, 0.604, 0.562, 0.530],  # 4000
]
# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def make_Z(mMed_vals, rinv_vals, raw):
    if len(raw) != len(mMed_vals):
        raise ValueError(
            f"Row mismatch: len(raw)={len(raw)} but len(mMed_vals)={len(mMed_vals)}"
        )
    for i, row in enumerate(raw):
        if len(row) != len(rinv_vals):
            raise ValueError(
                f"Column mismatch in row {i}: len(row)={len(row)} but len(rinv_vals)={len(rinv_vals)}"
            )

    Z = np.full((len(mMed_vals), len(rinv_vals)), np.nan)
    for i, row in enumerate(raw):
        for j, v in enumerate(row):
            Z[i, j] = v
    return Z


def save(fig, prefix):
    for ext in ("pdf", "png"):
        path = f"{prefix}.{ext}"
        fig.savefig(path, bbox_inches="tight", dpi=150 if ext == "png" else None)
        print(f"  -> {path}")
    plt.close(fig)


def plot_auc_grid(mMed_vals, rinv_vals, raw, tagger_name, out_prefix,
                  vmin=0.5, vmax=1.0, cmap="viridis"):
    """AUC heatmap with mPhi on x-axis and r_inv on y-axis."""
    Z = make_Z(mMed_vals, rinv_vals, raw)

    # Flip axes: rows become r_inv, columns become mPhi
    Zplot = Z.T
    nrows, ncols = len(rinv_vals), len(mMed_vals)

    fig_w = max(8, ncols * 1.25 + 3)
    fig_h = max(5, nrows * 0.75 + 2.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    im = ax.imshow(
        Zplot,
        origin="lower",
        aspect="auto",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        extent=[-0.5, ncols - 0.5, -0.5, nrows - 0.5],
    )

    outline = [
        pe.Stroke(linewidth=2.2, foreground="black"),
        pe.Normal()
    ]

    for i in range(nrows):
        for j in range(ncols):
            v = Zplot[i, j]
            if np.isfinite(v):
                ax.text(
                    j, i, f"{v:.3f}",
                    ha="center", va="center",
                    fontsize=13,
                    color="white",
                    fontweight="bold",
                    path_effects=outline,
                )
            else:
                ax.text(
                    j, i, "—",
                    ha="center", va="center",
                    fontsize=10,
                    color="white",
                    path_effects=outline,
                )

    ax.set_xticks(range(ncols))
    ax.set_xticklabels([str(m) for m in mMed_vals], fontsize=16)

    ax.set_yticks(range(nrows))
    ax.set_yticklabels([f"{r:.1f}" for r in rinv_vals], fontsize=16)

    ax.set_xlabel(r"$m_{\Phi}\ [\mathrm{GeV}]$", fontsize=28, labelpad=12)
    ax.set_ylabel(r"$r_{\mathrm{inv}}$", fontsize=28, labelpad=12)

    ax.tick_params(axis="both", labelsize=18)

    cb = fig.colorbar(im, ax=ax, pad=0.03)
    cb.set_label("AUC", fontsize=20, labelpad=14)
    cb.ax.tick_params(labelsize=14)

    hep.cms.label(
        data=False,
        ax=ax,
        loc=0,
        fontsize=26,
        com=13
    )

    plt.tight_layout()
    save(fig, out_prefix)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import os, argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="./paper_auc_plots")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    print("\nPlotting ParticleNet AUC grid …")
    plot_auc_grid(
        mMed_vals   = MMED_PNET,
        rinv_vals   = RINV_PNET,
        raw         = PNET_RAW,
        tagger_name = (r"ParticleNet SVJ Tagger  |  "
                       r"$m_{\rm dark}=20$ GeV, $\alpha={\rm peak}$  |  2018 QCD bkg"),
        out_prefix  = os.path.join(args.out_dir, "pnet_auc_grid"),
    )

    print("\nPlotting WNAE AUC grid …")
    plot_auc_grid(
        mMed_vals   = MMED_WNAE,
        rinv_vals   = RINV_WNAE,
        raw         = WNAE_RAW,
        tagger_name = (r"WNAE SVJ Tagger  |  "
                       r"$m_{\rm dark}=20$ GeV, $\alpha={\rm peak}$  |  2018 QCD bkg"),
        out_prefix  = os.path.join(args.out_dir, "wnae_auc_grid"),
    )

    print("\nPlotting DNN AUC grid …")
    plot_auc_grid(
        mMed_vals   = MMED_DNN,
        rinv_vals   = RINV_DNN,
        raw         = DNN_RAW,
        tagger_name = (r"DNN SVJ Tagger  |  "
                       r"$m_{\rm dark}=20$ GeV, $\alpha={\rm peak}$  |  2018 QCD bkg"),
        out_prefix  = os.path.join(args.out_dir, "dnn_auc_grid"),
    )
    print("\n=== Done ===")