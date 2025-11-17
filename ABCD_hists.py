#!/usr/bin/env python3
"""
ABCD stacked-region plot generator (PNET and WNAE).
Saves per-year / per-mode:
  {outdir}/{year}/{mode}_ABCD_counts.txt
  {outdir}/{year}/{mode}_ABCD_MC.pdf
  {outdir}/{year}/{mode}_ABCD_Data.pdf

Requirements: your project's utils.DataSetInfo (get2DHistoIntegral) available on PYTHONPATH.
"""
import os
import sys
import optparse
import numpy as np
import matplotlib.pyplot as plt
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

import utils.DataSetInfo as info

# ---------------------------
# CONFIG (global)
# ---------------------------
MET_CUT = 250.0
DNN_CUT = 0.85

ABCD_HIST_NAME = "h_METvsDNN"   # prefix inside ROOT files
# MAINCUT set per mode: "_pre_" or "_pre_WNAE_"

SVJ_ORDER = ["0SVJ", "1SVJ", "2SVJ", "3PSVJ"]
SVJ_XLABELS = ["0", "1", "2", "3+"]

# color mapping requested
COLOR_MAP = {
    "QCD": "royalblue",
    "TT": "gray",
    "WJets": "palevioletred",
    "ZJets": "darkorange",
    "ST": "forestgreen",
    "other": "gray",
}

LABEL_KEYWORDS = {
    "QCD": ["QCD"],
    "TT": ["TT", "TTJets", "t#bar{t}", "tt"],
    "WJets": ["WJets", "W+jets", "W+JETS", "Wjets", "W+jets"],
    "ZJets": ["ZJets", "Z#rightarrow#nu#nu", "Z->nu", "Z"],
    "ST": ["Single top", "ST"],
}

# ---------------------------
# Load datasets
# ---------------------------
def getData(path, scale=1.0, year="2018"):
    """
    Loads Data, Signal, Background DataSetInfo objects.
    Files are expected to be named like: {year}_WJets.root etc in 'path'.
    """
    Data = [
        info.DataSetInfo(basedir=path, fileName=year + "_Data.root", sys=-1.0, label="Data", scale=scale),
    ]
    bgData = [
        info.DataSetInfo(basedir=path, fileName=year + "_ST.root",     label="Single top", scale=scale, color=ROOT.TColor.GetColor("#5790fc")),
        info.DataSetInfo(basedir=path, fileName=year + "_TTJets.root", label="t#bar{t}",    scale=scale, color=ROOT.TColor.GetColor("#f89c20")),
        info.DataSetInfo(basedir=path, fileName=year + "_ZJets.root",  label="Z#rightarrow#nu#nu+jets", scale=scale, color=ROOT.TColor.GetColor("#e42536")),
        info.DataSetInfo(basedir=path, fileName=year + "_WJets.root",  label="W+jets",      scale=scale, color=ROOT.TColor.GetColor("#964a8b")),
        info.DataSetInfo(basedir=path, fileName=year + "_QCD.root",    label="QCD",         scale=scale, color=ROOT.TColor.GetColor("#9c9ca1")),
    ]
    sgData = [
        info.DataSetInfo(basedir=path, fileName=year + "_m2000_d20_r0p3_y1_N-1_M0_.root", label="mMed_2000", scale=scale, color=ROOT.kBlue + 1),
    ]
    return Data, sgData, bgData

def getMCOnly(dataset_dir, scale, year):
    """Return (Data=[], sgData, bgData) to handle the 2016-no-data case easily."""
    _, sgData, bgData = getData(dataset_dir, scale, year)
    return [], sgData, bgData

# ---------------------------
# Label/color utilities
# ---------------------------
def get_comp_label(comp):
    """Return a readable label for a DataSetInfo object."""
    for attr in ("label", "label_", "name", "fileName", "file_name"):
        if hasattr(comp, attr):
            val = getattr(comp, attr)
            if attr in ("fileName", "file_name") and isinstance(val, str):
                return os.path.splitext(val)[0]
            return str(val)
    return repr(comp)

def simplify_label_key(label):
    L = label.lower()
    for key, patterns in LABEL_KEYWORDS.items():
        for p in patterns:
            if p.lower() in L:
                return key
    return "other"

def color_for_label(label):
    key = simplify_label_key(label)
    return COLOR_MAP.get(key, COLOR_MAP["other"])

def clean_legend_label(raw):
    """Map a few raw labels to clean human-readable legend labels."""
    mapping = {
        "Single top": "Single top",
        "t#bar{t}": "ttbar",
        "Z#rightarrow#nu#nu+jets": "Z→νν",
        "W+jets": "W+jets",
        "QCD": "QCD"
    }
    return mapping.get(raw, raw)

# ---------------------------
# ABCD counting (robust)
# ---------------------------
def _safe_get2DHistoIntegral(data_obj, histName, xmin, xmax, ymin, ymax, showEvents=False):
    """
    Wrap data_obj.get2DHistoIntegral; be tolerant of differing return signatures and errors.
    Returns tuple (histo_or_None, integral_float, error_float)
    On failure returns (None, 0.0, 0.0) and prints a warning.
    """
    try:
        # Try named args (preferred)
        result = data_obj.get2DHistoIntegral(histName, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, showEvents=showEvents)
    except TypeError:
        # maybe older signature expects positional args
        try:
            result = data_obj.get2DHistoIntegral(histName, xmin, xmax, ymin, ymax, showEvents)
        except Exception as e:
            print(f"[WARN] get2DHistoIntegral failed for '{histName}' on {get_comp_label(data_obj)}: {e}")
            return None, 0.0, 0.0
    except Exception as e:
        print(f"[WARN] get2DHistoIntegral exception for '{histName}' on {get_comp_label(data_obj)}: {e}")
        return None, 0.0, 0.0

    # normalize different return shapes
    if result is None:
        return None, 0.0, 0.0

    # if result is a ROOT histogram (rare), treat as no integral
    try:
        # if result is iterable / tuple
        length = len(result)
    except Exception:
        # single value
        try:
            val = float(result)
            return None, val, 0.0
        except Exception:
            return None, 0.0, 0.0

    # result may be (histo, integral, error) or (integral, error)
    if length == 3:
        histo, integral, error = result
        try:
            integral_f = float(integral)
        except Exception:
            integral_f = 0.0
        try:
            error_f = float(error)
        except Exception:
            error_f = 0.0
        return histo, integral_f, error_f
    elif length == 2:
        # (integral, error)
        integral, error = result
        try:
            integral_f = float(integral)
        except Exception:
            integral_f = 0.0
        try:
            error_f = float(error)
        except Exception:
            error_f = 0.0
        return None, integral_f, error_f
    else:
        # unknown shape
        print(f"[WARN] Unexpected return from get2DHistoIntegral for '{histName}' on {get_comp_label(data_obj)}: got {result}")
        return None, 0.0, 0.0

def get_region_count_from_dataset(data_obj, svj_label, maincut, met_cut=MET_CUT, dnn_cut=DNN_CUT):
    """
    Returns dict {'A': (n,err), 'B':..., 'C':..., 'D':...}
    Builds histName as: ABCD_HIST_NAME + maincut + svj_label
    maincut should include leading/trailing underscores as needed, e.g. "_pre_" or "_pre_WNAE_"
    """
    histName = f"{ABCD_HIST_NAME}{maincut}{svj_label}"
    # Regions:
    # A: MET > met_cut, DNN > dnn_cut
    _, A_int, A_err = _safe_get2DHistoIntegral(data_obj, histName, xmin=met_cut, xmax=1e9, ymin=dnn_cut, ymax=1.0, showEvents=True)
    _, B_int, B_err = _safe_get2DHistoIntegral(data_obj, histName, xmin=0.0,   xmax=met_cut, ymin=dnn_cut, ymax=1.0, showEvents=True)
    _, C_int, C_err = _safe_get2DHistoIntegral(data_obj, histName, xmin=met_cut, xmax=1e9, ymin=0.0,     ymax=dnn_cut, showEvents=True)
    _, D_int, D_err = _safe_get2DHistoIntegral(data_obj, histName, xmin=0.0,   xmax=met_cut, ymin=0.0,     ymax=dnn_cut, showEvents=True)

    return {
        "A": (float(A_int), float(A_err)),
        "B": (float(B_int), float(B_err)),
        "C": (float(C_int), float(C_err)),
        "D": (float(D_int), float(D_err)),
    }

def build_bg_counts_by_SVJ(bgData, maincut, svj_list=SVJ_ORDER):
    """
    Return counts[(svj,label)] = (A,B,C,D), aggregated dict per svj, components list
    """
    counts = {}
    components = []
    for comp in bgData:
        label = get_comp_label(comp)
        components.append(label)
        for svj in svj_list:
            reg = get_region_count_from_dataset(comp, svj, maincut)
            counts[(svj, label)] = (reg["A"][0], reg["B"][0], reg["C"][0], reg["D"][0])

    aggregated = {}
    for svj in svj_list:
        A_tot = sum(counts[(svj, c)][0] for c in components)
        B_tot = sum(counts[(svj, c)][1] for c in components)
        C_tot = sum(counts[(svj, c)][2] for c in components)
        D_tot = sum(counts[(svj, c)][3] for c in components)
        aggregated[svj] = {"A": A_tot, "B": B_tot, "C": C_tot, "D": D_tot}
    return counts, aggregated, components

def build_data_counts_by_SVJ(dataList, maincut, svj_list=SVJ_ORDER):
    """
    dataList typically has a single Data DataSetInfo; return dict keyed by (svj,'Data')
    """
    data_counts = {}
    for comp in dataList:
        label = "Data"
        for svj in svj_list:
            reg = get_region_count_from_dataset(comp, svj, maincut)
            data_counts[(svj, label)] = (reg["A"][0], reg["B"][0], reg["C"][0], reg["D"][0])
    return data_counts

# ---------------------------
# TXT writer (nicely formatted)
# ---------------------------
def save_all_to_txt(counts, aggregated, components, outpath, mode_label, year):
    with open(outpath, "w") as f:
        f.write(f"# ABCD summary (background MC) - mode: {mode_label}  year: {year}\n")
        f.write("# MET_CUT = {:.1f}, DNN_CUT = {:.3f}\n\n".format(MET_CUT, DNN_CUT))

        f.write("Per-component counts (A,B,C,D,TOTAL)\n")
        f.write("{:<8s} {:<20s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s}\n".format("SVJ", "Component", "A", "B", "C", "D", "TOTAL"))
        f.write("-" * 86 + "\n")

        for svj in SVJ_ORDER:
            for comp in components:
                A,B,C,D = counts[(svj, comp)]
                total = A + B + C + D
                f.write("{:<8s} {:<20s} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f}\n".format(svj, comp[:20], A, B, C, D, total))

        f.write("\nAggregated totals (sum over backgrounds)\n")
        f.write("{:<8s} {:>14s} {:>14s} {:>14s} {:>14s}\n".format("SVJ", "A", "B", "C", "D"))
        f.write("-" * 70 + "\n")
        for svj in SVJ_ORDER:
            agg = aggregated[svj]
            f.write("{:<8s} {:14.6f} {:14.6f} {:14.6f} {:14.6f}\n".format(svj, agg["A"], agg["B"], agg["C"], agg["D"]))

        total_A = sum(aggregated[svj]["A"] for svj in SVJ_ORDER)
        total_B = sum(aggregated[svj]["B"] for svj in SVJ_ORDER)
        total_C = sum(aggregated[svj]["C"] for svj in SVJ_ORDER)
        total_D = sum(aggregated[svj]["D"] for svj in SVJ_ORDER)
        f.write("-" * 70 + "\n")
        f.write("{:<8s} {:14.6f} {:14.6f} {:14.6f} {:14.6f}\n".format("TOTAL", total_A, total_B, total_C, total_D))
    print(f"[OK] Written ABCD summary to: {outpath}")

# ---------------------------
# Plotting 4-panel ABCD
# ---------------------------
def plot_abcd_panels(counts, aggregated, components, outpng, title_suffix="", data_counts=None):
    svj_list = SVJ_ORDER
    x_labels = SVJ_XLABELS
    n_svj = len(svj_list)
    n_comp = len(components)
    region_names = ["A", "B", "C", "D"]

    clean_components = [clean_legend_label(c) for c in components]
    comp_colors = [color_for_label(c) for c in components]

    fig, axes = plt.subplots(1, 4, figsize=(18, 5), sharey=True)

    # CMS label (bold) left above plot
    fig.text(0.02, 0.98, "CMS Preliminary", fontsize=16, fontweight="bold", ha="left", va="top")
    fig.text(0.22, 0.98, "41.5 fb$^{-1}$ (13 TeV)", fontsize=14, ha="left", va="top")

    ymax = 0.0
    for ax_idx, region in enumerate(region_names):
        ax = axes[ax_idx]
        bottoms = np.zeros(n_svj)
        # stacked bars
        for i_comp in range(n_comp):
            vals = np.array([counts[(svj, components[i_comp])][ "ABCD".index(region) ] if False else counts[(svj, components[i_comp])][{"A":0,"B":1,"C":2,"D":3}[region]] for svj in svj_list])
            ax.bar(np.arange(n_svj), vals, bottom=bottoms, color=comp_colors[i_comp], label=clean_components[i_comp] if ax_idx==0 else None, edgecolor="none")
            bottoms += vals
        ymax = max(ymax, np.max(bottoms) if bottoms.size > 0 else 0.0)

        # optional data points for B,C,D
        if (data_counts is not None) and (region != "A"):
            data_vals = []
            for svj in svj_list:
                d = data_counts.get((svj, "Data"), (0.0, 0.0, 0.0, 0.0))
                # d is (A,B,C,D)
                if region == "B": data_vals.append(d[1])
                elif region == "C": data_vals.append(d[2])
                elif region == "D": data_vals.append(d[3])
            ax.errorbar(np.arange(n_svj), data_vals, fmt="ko", markersize=5, label="Data" if ax_idx==1 else None, zorder=10)

        ax.set_xticks(np.arange(n_svj))
        ax.set_xticklabels(x_labels)
        ax.set_title(region, fontsize=12, fontweight="bold")
        if ax_idx == 0:
            ax.set_ylabel("Events (background MC)")
        ax.set_xlabel("nSVJ")
        ax.set_yscale("log")
        ax.grid(axis="y", linestyle="--", alpha=0.3)

        # thin separators
        for sep in np.arange(0.5, n_svj - 0.5 + 1, 1.0):
            ax.axvline(sep - 0.5, color='black', linestyle='-', linewidth=0.4, alpha=0.08)

    # shared y-limits
    ymin = 0.1
    ymax_safe = max(ymax * 1.6, 1.0)
    ytop = 10 ** np.ceil(np.log10(ymax_safe)) if ymax_safe > 0 else 1.0
    axes[0].set_ylim(ymin, ytop)

    # legend (clean)
    handles = [plt.Rectangle((0,0),1,1, color=c) for c in comp_colors]
    labels = clean_components.copy()
    if data_counts is not None:
        handles.append(plt.Line2D([0],[0], marker='o', color='k', linestyle='', markersize=5))
        labels.append("Data")
    ncol = min(6, len(labels))
    fig.legend(handles, labels, loc='upper center', ncol=ncol, bbox_to_anchor=(0.5, 0.915), frameon=False)

    fig.tight_layout(rect=[0, 0, 1, 0.90])
    plt.savefig(outpng, dpi=300)
    plt.close(fig)
    print(f"[OK] Saved -> {outpng}")

# ---------------------------
# Single-mode runner
# ---------------------------
def run_mode(mode_label, maincut, Data, sgData, bgData, outdir, year):
    mode_prefix = mode_label.replace(" ", "_")
    counts, aggregated, components = build_bg_counts_by_SVJ(bgData, maincut, SVJ_ORDER)

    # Write TXT
    txtpath = os.path.join(outdir, f"{mode_prefix}_ABCD_counts_{year}.txt")
    save_all_to_txt(counts, aggregated, components, txtpath, mode_label, year)

    # MC-only plot
    mc_pdf = os.path.join(outdir, f"{mode_prefix}_ABCD_MC_{year}.pdf")
    plot_abcd_panels(counts, aggregated, components, mc_pdf, title_suffix=f"({mode_label})", data_counts=None)

    # Data overlay for B,C,D if Data present
    data_counts = {}
    if Data:
        data_counts = build_data_counts_by_SVJ(Data, maincut, SVJ_ORDER)
    data_pdf = os.path.join(outdir, f"{mode_prefix}_ABCD_Data_{year}.pdf")
    plot_abcd_panels(counts, aggregated, components, data_pdf, title_suffix=f"({mode_label})", data_counts=data_counts if data_counts else None)

# ---------------------------
# Main
# ---------------------------
def main():
    parser = optparse.OptionParser()
    parser.add_option('-d', '--dataset', dest='dataset', default='output/Current_Model_wp90', help='dataset base dir (contains YEAR_*.root)')
    parser.add_option('-o', '--outdir', dest='outdir', default='ABCD_output', help='output dir')
    opts, args = parser.parse_args()

    base_dataset = opts.dataset.rstrip("/") + "/"
    outdir = opts.outdir
    os.makedirs(outdir, exist_ok=True)

    years = ["2016", "2017", "2018"]

    for yr in years:
        print(f"\n[INFO] Processing {yr}")
        yout = os.path.join(outdir, yr)
        os.makedirs(yout, exist_ok=True)

        if yr == "2016":
            Data, sgData, bgData = getMCOnly(base_dataset, 1.0, yr)  # no data for 2016
        else:
            Data, sgData, bgData = getData(base_dataset, 1.0, yr)

        # PNET
        run_mode("PNET", "_pre_", Data, sgData, bgData, yout, yr)

        # WNAE
        run_mode("WNAE", "_pre_WNAE_", Data, sgData, bgData, yout, yr)

    print("\nAll done.")

if __name__ == "__main__":
    main()
