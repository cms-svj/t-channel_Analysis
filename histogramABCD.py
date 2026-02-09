#!/usr/bin/env python3
"""
ABCD stacked-yield plotter from 2D (MET vs DNN) histograms using utils.DataSetInfo.

Produces (per year, per mode):
  - {outdir}/{year}/{mode}_ABCD_counts_{year}.txt
  - {outdir}/{year}/{mode}_ABCD_yields_{year}_CMS.pdf              (MC only)
  - {outdir}/{year}/{mode}_ABCD_yields_{year}_CMS_with_data.pdf    (if --with-data and data exists)

Layout matches the "continuous A|B|C|D" CMS-style plot:
  A: nSVJ 0 1 2 3+ | B: 0 1 2 3+ | C: 0 1 2 3+ | D: 0 1 2 3+
"""

import os
import math
import optparse
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

import utils.DataSetInfo as info

# ---------------------------
# CONFIG (global)
# ---------------------------
MET_CUT = 250.0
DNN_CUT = 0.85
ABCD_HIST_NAME = "h_METvsDNN"  # prefix inside ROOT files

SVJ_ORDER   = ["0SVJ", "1SVJ", "2SVJ", "3PSVJ"]
SVJ_XLABELS = ["0", "1", "2", "3+"]

PLOT_REGIONS = ["A", "B", "C", "D"]

LUMI_BY_YEAR = {"2016": "35.9", "2017": "41.5", "2018": "59.7"}  # fb^-1

# Fixed stacking order (bottom -> top is the *reverse* of draw order below)
BKG_ORDER = ["QCD", "TT", "WJets", "ZJets", "ST"]

PROC_LABEL = {
    "QCD":  "QCD multijet",
    "TT":   r"$t\bar{t}$ + jets",
    "ZJets": r"$Z \to \nu\nu$ + jets",
    "WJets": r"$W \to \ell\nu$ + jets",
    "ST":   "Single top",
}

# Colors consistent with your second-script palette choice (edit as you like)
PROC_COLOR = {
    "QCD":  "#f89c20",
    "TT":   "#e42536",
    "WJets":"#5790fc",
    "ZJets":"#9c9ca1",
    "ST":   "#7a21dd",
}
PROC_HATCH = {k: None for k in PROC_COLOR.keys()}

# Heuristic mapping from DataSetInfo labels/filenames -> proc key
LABEL_KEYWORDS = {
    "QCD":  ["qcd"],
    "TT":   ["tt", "ttjets", "t#bar{t}", "tbar", "ttbar"],
    "WJets":["wjets", "w+jets", "wjets", "wto", "w->"],
    "ZJets":["zjets", "z#rightarrow#nu#nu", "z->nu", "z→", "z to", "znunu"],
    "ST":   ["single top", "st", "singletop"],
}

# ---------------------------
# Dataset loading
# ---------------------------
def getData(path, scale=1.0, year="2018"):
    Data = [
        info.DataSetInfo(basedir=path, fileName=f"{year}_Data.root", sys=-1.0, label="Data", scale=scale),
    ]
    bgData = [
        info.DataSetInfo(basedir=path, fileName=f"{year}_ST.root",     label="Single top",              scale=scale),
        info.DataSetInfo(basedir=path, fileName=f"{year}_TTJets.root", label="t#bar{t}",                scale=scale),
        info.DataSetInfo(basedir=path, fileName=f"{year}_ZJets.root",  label="Z#rightarrow#nu#nu+jets", scale=scale),
        info.DataSetInfo(basedir=path, fileName=f"{year}_WJets.root",  label="W+jets",                  scale=scale),
        info.DataSetInfo(basedir=path, fileName=f"{year}_QCD.root",    label="QCD",                     scale=scale),
    ]
    sgData = []  # keep for later
    return Data, sgData, bgData

def _dataset_exists(basedir, filename):
    return os.path.exists(os.path.join(basedir, filename))

def getData_with_missing_data_ok(path, scale=1.0, year="2018"):
    Data, sgData, bgData = getData(path, scale, year)
    datafile = f"{year}_Data.root"
    if not _dataset_exists(path, datafile):
        print(f"[INFO] No data file found for {year} ({datafile}); will produce MC-only plots unless --with-data is off.")
        Data = []
    return Data, sgData, bgData

# ---------------------------
# Label helpers
# ---------------------------
def get_comp_label(comp):
    for attr in ("label", "label_", "name", "fileName", "file_name"):
        if hasattr(comp, attr):
            val = getattr(comp, attr)
            if attr in ("fileName", "file_name") and isinstance(val, str):
                return os.path.splitext(val)[0]
            return str(val)
    return repr(comp)

def simplify_to_proc_key(label: str) -> str:
    L = label.lower()
    for key, patterns in LABEL_KEYWORDS.items():
        for p in patterns:
            if p in L:
                return key
    # try filename-like tokens
    if "tt" in L:
        return "TT"
    return "QCD" if "qcd" in L else "ST" if "st" in L else "WJets" if "w" in L else "ZJets" if "z" in L else "QCD"

# ---------------------------
# Robust integral wrapper
# ---------------------------
def _safe_get2DHistoIntegral(data_obj, histName, xmin, xmax, ymin, ymax, showEvents=False):
    """
    Normalizes return shape of DataSetInfo.get2DHistoIntegral across possible signatures:
      - (hist, integral, error)
      - (integral, error)
      - integral
    Returns: (hist_or_None, integral, error)
    """
    try:
        result = data_obj.get2DHistoIntegral(
            histName, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, showEvents=showEvents
        )
    except TypeError:
        try:
            result = data_obj.get2DHistoIntegral(histName, xmin, xmax, ymin, ymax, showEvents)
        except Exception as e:
            print(f"[WARN] get2DHistoIntegral failed for '{histName}' on {get_comp_label(data_obj)}: {e}")
            return None, 0.0, 0.0
    except Exception as e:
        print(f"[WARN] get2DHistoIntegral exception for '{histName}' on {get_comp_label(data_obj)}: {e}")
        return None, 0.0, 0.0

    if result is None:
        return None, 0.0, 0.0

    # Normalize
    try:
        n = len(result)
    except Exception:
        try:
            return None, float(result), 0.0
        except Exception:
            return None, 0.0, 0.0

    if n == 3:
        h, integ, err = result
        return h, float(integ or 0.0), float(err or 0.0)
    if n == 2:
        integ, err = result
        return None, float(integ or 0.0), float(err or 0.0)

    print(f"[WARN] Unexpected return for '{histName}' on {get_comp_label(data_obj)}: {result}")
    return None, 0.0, 0.0

# def get_region_count_from_dataset(data_obj, svj_label, maincut, met_cut=MET_CUT, dnn_cut=DNN_CUT):
#     histName = f"{ABCD_HIST_NAME}{maincut}{svj_label}"

#     # A: MET > met_cut, DNN > dnn_cut
#     _, A_int, A_err = _safe_get2DHistoIntegral(data_obj, histName, xmin=met_cut, xmax=1e9, ymin=dnn_cut, ymax=1.0, showEvents=True)
#     # B: MET < met_cut, DNN > dnn_cut
#     _, B_int, B_err = _safe_get2DHistoIntegral(data_obj, histName, xmin=0.0,   xmax=met_cut, ymin=dnn_cut, ymax=1.0, showEvents=True)
#     # C: MET > met_cut, DNN < dnn_cut
#     _, C_int, C_err = _safe_get2DHistoIntegral(data_obj, histName, xmin=met_cut, xmax=1e9, ymin=0.0,     ymax=dnn_cut, showEvents=True)
#     # D: MET < met_cut, DNN < dnn_cut
#     _, D_int, D_err = _safe_get2DHistoIntegral(data_obj, histName, xmin=0.0,   xmax=met_cut, ymin=0.0,     ymax=dnn_cut, showEvents=True)

#     return {
#         "A": (A_int, A_err),
#         "B": (B_int, B_err),
#         "C": (C_int, C_err),
#         "D": (D_int, D_err),
#     }
def get_region_count_from_dataset(data_obj, svj_label, maincut, met_cut=MET_CUT, dnn_cut=DNN_CUT):
    histName = f"{ABCD_HIST_NAME}{maincut}{svj_label}"

    # A: MET > met_cut, DNN > dnn_cut
    _, A_int, A_err = _safe_get2DHistoIntegral(data_obj, histName,
        xmin=met_cut, xmax=1e9, ymin=dnn_cut, ymax=1.0, showEvents=True)

    # B: MET < met_cut, DNN > dnn_cut
    _, B_int, B_err = _safe_get2DHistoIntegral(data_obj, histName,
        xmin=0.0, xmax=met_cut, ymin=dnn_cut, ymax=1.0, showEvents=True)

    # C: MET > met_cut, DNN < dnn_cut
    _, C_int, C_err = _safe_get2DHistoIntegral(data_obj, histName,
        xmin=met_cut, xmax=1e9, ymin=0.0, ymax=dnn_cut, showEvents=True)

    # # B: MET > met_cut, DNN < dnn_cut
    # _, B_int, B_err = _safe_get2DHistoIntegral(data_obj, histName,
    #     xmin=met_cut, xmax=1e9, ymin=0.0, ymax=dnn_cut, showEvents=True)

    # # C: MET < met_cut, DNN > dnn_cut
    # _, C_int, C_err = _safe_get2DHistoIntegral(data_obj, histName,
    #     xmin=0.0, xmax=met_cut, ymin=dnn_cut, ymax=1.0, showEvents=True)

    # D: MET < met_cut, DNN < dnn_cut
    _, D_int, D_err = _safe_get2DHistoIntegral(data_obj, histName,
        xmin=0.0, xmax=met_cut, ymin=0.0, ymax=dnn_cut, showEvents=True)

    return {"A": (A_int, A_err), "B": (B_int, B_err), "C": (C_int, C_err), "D": (D_int, D_err)}


# ---------------------------
# Yield building (per-proc, per bin)
# ---------------------------
def build_x_layout():
    """
    Continuous layout:
      A: 0 1 2 3 | B: 4 5 6 7 | C: 8 9 10 11 | D: 12 13 14 15
    """
    nsvj = len(SVJ_ORDER)
    x = []
    xticks = []
    xticklabels = []
    region_centers = {}
    region_boundaries = []

    for r_i, reg in enumerate(PLOT_REGIONS):
        base = r_i * nsvj
        if r_i > 0:
            region_boundaries.append(base - 0.5)

        for s_i in range(nsvj):
            xpos = base + s_i
            x.append(xpos)
            xticks.append(xpos)
            xticklabels.append(SVJ_XLABELS[s_i])

        region_centers[reg] = base + (nsvj - 1) / 2.0

    return np.array(x, dtype=float), xticks, xticklabels, region_centers, region_boundaries

def step_band_from_bins(xpos, y, yerr, width=1.0, eps=1e-6):
    left = xpos - width / 2.0
    right = xpos + width / 2.0

    xe = np.empty(2 * len(xpos))
    yhi = np.empty(2 * len(xpos))
    ylo = np.empty(2 * len(xpos))

    for i in range(len(xpos)):
        xe[2*i]   = left[i]
        xe[2*i+1] = right[i]
        up = y[i] + yerr[i]
        lo = max(eps, y[i] - yerr[i])
        yhi[2*i] = up
        yhi[2*i+1] = up
        ylo[2*i] = lo
        ylo[2*i+1] = lo

    return xe, ylo, yhi

def collect_yields(bgData, dataList, maincut):
    """
    Returns:
      proc_y[proc][idx], proc_e[proc][idx]
      data_y[idx], data_e[idx] (nan if missing)
      Also returns a component-level dict for txt writing.
    """
    x, _, _, _, _ = build_x_layout()
    nbin = len(x)
    nsvj = len(SVJ_ORDER)

    proc_y = {p: np.zeros(nbin, dtype=float) for p in BKG_ORDER}
    proc_e = {p: np.zeros(nbin, dtype=float) for p in BKG_ORDER}

    data_y = np.full(nbin, np.nan, dtype=float)
    data_e = np.full(nbin, np.nan, dtype=float)

    # component-level bookkeeping for your txt output
    comp_counts = {}  # (svj, label) -> (A,B,C,D)
    comp_errs   = {}  # (svj, label) -> (Ae,Be,Ce,De)
    components  = []

    # backgrounds
    for comp in bgData:
        label = get_comp_label(comp)
        components.append(label)
        proc_key = simplify_to_proc_key(label)

        for r_i, reg in enumerate(PLOT_REGIONS):
            for s_i, svj in enumerate(SVJ_ORDER):
                idx = r_i * nsvj + s_i
                regmap = get_region_count_from_dataset(comp, svj, maincut)
                y, e = regmap[reg]
                proc_y[proc_key][idx] += y
                # treat errors as uncorrelated across datasets: add in quadrature
                proc_e[proc_key][idx] = math.sqrt(proc_e[proc_key][idx]**2 + e**2)

        # store per-component A/B/C/D (summed over nothing; just that dataset)
        for svj in SVJ_ORDER:
            regmap = get_region_count_from_dataset(comp, svj, maincut)
            comp_counts[(svj, label)] = (regmap["A"][0], regmap["B"][0], regmap["C"][0], regmap["D"][0])
            comp_errs[(svj, label)]   = (regmap["A"][1], regmap["B"][1], regmap["C"][1], regmap["D"][1])

    # data (assume a single DataSetInfo entry, but handle list)
    for comp in dataList:
        for r_i, reg in enumerate(PLOT_REGIONS):
            for s_i, svj in enumerate(SVJ_ORDER):
                idx = r_i * nsvj + s_i
                regmap = get_region_count_from_dataset(comp, svj, maincut)
                y, e = regmap[reg]
                data_y[idx] = y
                data_e[idx] = e

    return proc_y, proc_e, data_y, data_e, comp_counts, comp_errs, components

# ---------------------------
# TXT writer
# ---------------------------
def save_all_to_txt(comp_counts, components, outpath, mode_label, year):
    # aggregate totals over datasets per svj
    aggregated = {svj: {"A": 0.0, "B": 0.0, "C": 0.0, "D": 0.0} for svj in SVJ_ORDER}
    for svj in SVJ_ORDER:
        for comp in components:
            A, B, C, D = comp_counts[(svj, comp)]
            aggregated[svj]["A"] += A
            aggregated[svj]["B"] += B
            aggregated[svj]["C"] += C
            aggregated[svj]["D"] += D

    with open(outpath, "w") as f:
        f.write(f"# ABCD summary (background MC) - mode: {mode_label}  year: {year}\n")
        f.write("# MET_CUT = {:.1f}, DNN_CUT = {:.3f}\n\n".format(MET_CUT, DNN_CUT))

        f.write("Per-component counts (A,B,C,D,TOTAL)\n")
        f.write("{:<8s} {:<30s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s}\n"
                .format("SVJ", "Component", "A", "B", "C", "D", "TOTAL"))
        f.write("-" * 102 + "\n")

        for svj in SVJ_ORDER:
            for comp in components:
                A, B, C, D = comp_counts[(svj, comp)]
                total = A + B + C + D
                f.write("{:<8s} {:<30s} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f}\n"
                        .format(svj, comp[:30], A, B, C, D, total))

        f.write("\nAggregated totals (sum over backgrounds)\n")
        f.write("{:<8s} {:>14s} {:>14s} {:>14s} {:>14s}\n".format("SVJ", "A", "B", "C", "D"))
        f.write("-" * 72 + "\n")
        for svj in SVJ_ORDER:
            agg = aggregated[svj]
            f.write("{:<8s} {:14.6f} {:14.6f} {:14.6f} {:14.6f}\n"
                    .format(svj, agg["A"], agg["B"], agg["C"], agg["D"]))

    print(f"[OK] Written ABCD summary to: {outpath}")

# ---------------------------
# Plotting (CMS-style continuous layout)
# ---------------------------
def plot_abcd_continuous(proc_y, proc_e, year, outpdf, with_data=False, data_y=None, data_e=None,
                         lumi_text=None, com_text="13", prelim=True):
    hep.style.use("CMS")
    plt.rcParams["figure.dpi"] = 150

    x, xticks, xticklabels, region_centers, region_boundaries = build_x_layout()

    # main + ratio panel if data requested and present
    has_data = with_data and (data_y is not None) and np.any(np.isfinite(data_y))

    if has_data:
        fig, (ax, rax) = plt.subplots(
            2, 1, figsize=(14.0, 9.5),
            gridspec_kw={"height_ratios": (3.4, 1.0), "hspace": 0.03},
            sharex=True,
            constrained_layout=True,
        )
    else:
        fig, ax = plt.subplots(1, 1, figsize=(14.0, 7.5), constrained_layout=True)
        rax = None

    # Stack: draw reverse so first in BKG_ORDER ends up on bottom
    bottoms = np.zeros_like(x, dtype=float)
    for p in reversed(BKG_ORDER):
        y = proc_y.get(p, np.zeros_like(x))
        ax.bar(
            x, y,
            bottom=bottoms,
            width=1.0,
            label=PROC_LABEL.get(p, p),
            color=PROC_COLOR.get(p, None),
            edgecolor="black",
            linewidth=0.2,
            hatch=PROC_HATCH.get(p, None),
            zorder=1,
        )
        bottoms += y

    # total MC + uncertainty (quadrature over proc_e)
    mc_tot = np.zeros_like(x, dtype=float)
    mc_var = np.zeros_like(x, dtype=float)
    for p in BKG_ORDER:
        mc_tot += proc_y[p]
        mc_var += proc_e[p] ** 2
    mc_err = np.sqrt(mc_var)

    xe, ylo, yhi = step_band_from_bins(x, mc_tot, mc_err, width=1.0, eps=1e-6)
    ax.fill_between(
        xe, ylo, yhi,
        step="pre",
        facecolor="none",
        edgecolor="0.6",
        linewidth=0.0,
        hatch="////",
        label="MC unc.",
        zorder=2,
    )

    # data overlay + ratio
    if has_data:
        m = np.isfinite(data_y)
        ax.errorbar(
            x[m], data_y[m], yerr=data_e[m],
            fmt="o", color="black", markersize=4,
            linewidth=1.0, capsize=0,
            label="Data",
            zorder=5,
        )

        denom = mc_tot.copy()
        ok = m & (denom > 0)
        ratio = np.full_like(x, np.nan, dtype=float)
        ratio_err = np.full_like(x, np.nan, dtype=float)
        ratio[ok] = data_y[ok] / denom[ok]
        ratio_err[ok] = data_e[ok] / denom[ok]  # data stat only

        rax.errorbar(
            x[ok], ratio[ok], yerr=ratio_err[ok],
            fmt="o", color="black", markersize=4,
            linewidth=1.0, capsize=0,
            zorder=5,
        )

        # unity line + MC unc band (relative)
        rax.axhline(1.0, color="black", linewidth=1.0, zorder=1)
        rel_lo = np.ones_like(mc_tot)
        rel_hi = np.ones_like(mc_tot)
        pos = mc_tot > 0
        rel_lo[pos] = np.maximum(0.0, (mc_tot[pos] - mc_err[pos]) / mc_tot[pos])
        rel_hi[pos] = (mc_tot[pos] + mc_err[pos]) / mc_tot[pos]

        # build step band for ratio
        left = x - 0.5
        right = x + 0.5
        xe_band = np.empty(2 * len(x))
        rlo_band = np.empty(2 * len(x))
        rhi_band = np.empty(2 * len(x))
        for i in range(len(x)):
            xe_band[2*i] = left[i]
            xe_band[2*i+1] = right[i]
            rlo_band[2*i] = rel_lo[i]
            rlo_band[2*i+1] = rel_lo[i]
            rhi_band[2*i] = rel_hi[i]
            rhi_band[2*i+1] = rel_hi[i]

        rax.fill_between(
            xe_band, rlo_band, rhi_band,
            step="pre",
            facecolor="none",
            edgecolor="gray",
            linewidth=0.0,
            hatch="////",
            zorder=2,
        )
        rax.set_ylabel("Data/MC")
        rax.set_ylim(0.0, 2.0)
        rax.grid(True, which="both", axis="y", linestyle=":", linewidth=0.8)

    # region separators
    for xb in region_boundaries:
        ax.axvline(xb, color="black", linewidth=1.6)
        if has_data:
            rax.axvline(xb, color="black", linewidth=1.6)

    # region labels inside plot
    ymax = float(np.nanmax(mc_tot)) if np.nanmax(mc_tot) > 0 else 1.0
    for reg, xc in region_centers.items():
        ax.text(xc, ymax * 0.35, reg, ha="center", va="bottom", fontsize=15, fontweight="bold")

    # axes formatting
    ax.set_yscale("log")
    ax.set_ylabel("Events")
    ax.set_xlim(-0.5, len(x) - 0.5)
    ax.set_ylim(1e-2, ymax * 50.0 if ymax > 0 else 10.0)

    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    if not has_data:
        ax.set_xlabel("nSVJ")
    else:
        rax.set_xlabel("nSVJ")

    ax.grid(True, which="both", axis="y", linestyle=":", linewidth=0.8)

    # legend
    ax.legend(
        loc="upper left",
        ncol=2,
        frameon=False,
        fontsize=10,
        handlelength=1.6,
        columnspacing=1.0,
    )

    if lumi_text is None:
        lumi_text = LUMI_BY_YEAR.get(year, "")
    hep.cms.label(
        ax=ax,
        label="Preliminary" if prelim else "",
        data=has_data,
        lumi=lumi_text,
        com=com_text,
    )

    fig.savefig(outpdf)
    plt.close(fig)
    print(f"[OK] wrote {outpdf}")

# ---------------------------
# Runner
# ---------------------------
def run_mode(mode_label, maincut, Data, bgData, outdir, year, with_data=False):
    mode_prefix = mode_label.replace(" ", "_")

    proc_y, proc_e, data_y, data_e, comp_counts, comp_errs, components = collect_yields(bgData, Data, maincut)

    txtpath = os.path.join(outdir, f"{mode_prefix}_ABCD_counts_{year}.txt")
    save_all_to_txt(comp_counts, components, txtpath, mode_label, year)

    # MC-only
    outpdf_mc = os.path.join(outdir, f"{mode_prefix}_ABCD_yields_{year}_CMS.pdf")
    plot_abcd_continuous(proc_y, proc_e, year, outpdf_mc, with_data=False)

    # With data (+ratio), if requested
    if with_data:
        outpdf_data = os.path.join(outdir, f"{mode_prefix}_ABCD_yields_{year}_CMS_with_data.pdf")
        plot_abcd_continuous(proc_y, proc_e, year, outpdf_data, with_data=True, data_y=data_y, data_e=data_e)

def main():
    parser = optparse.OptionParser()
    parser.add_option("-d", "--dataset", dest="dataset", default="output/Current_Model_wp90",
                      help="dataset base dir (contains YEAR_*.root)")
    parser.add_option("-o", "--outdir", dest="outdir", default="ABCD_output",
                      help="output dir")
    parser.add_option("--with-data", action="store_true", dest="with_data", default=False,
                      help="overlay data and draw ratio panel (requires YEAR_Data.root)")
    parser.add_option("--years", dest="years", default="2016,2017,2018",
                      help="comma-separated years, e.g. 2017,2018")
    opts, _ = parser.parse_args()

    base_dataset = opts.dataset.rstrip("/") + "/"
    outdir = opts.outdir
    os.makedirs(outdir, exist_ok=True)

    years = [y.strip() for y in opts.years.split(",") if y.strip()]
    for yr in years:
        print(f"\n[INFO] Processing {yr}")
        yout = os.path.join(outdir, yr)
        os.makedirs(yout, exist_ok=True)

        Data, _, bgData = getData_with_missing_data_ok(base_dataset, 1.0, yr)
    
        run_mode("PNET", "_pre_",      Data, bgData, yout, yr, with_data=opts.with_data)
        run_mode("WNAE",  "_pre_WNAE_", Data, bgData, yout, yr, with_data=opts.with_data)

    print("\nAll done.")

if __name__ == "__main__":
    main()
