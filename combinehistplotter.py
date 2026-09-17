#!/usr/bin/env python3
import os
import re
import math
import argparse
import numpy as np
import uproot

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import mplhep as hep

import Figure2_makerusingskims as f2

# --------------------------
# CONFIG
# --------------------------

# Plot/print region names (what you show)
PLOT_REGIONS = ["A", "B", "C", "D"]

# ROOT region mapping (what you read)
# (A<->B swap stays; C<->D swap added)
#ROOT_REGION_FOR = {"A": "B", "B": "A", "C": "D", "D": "C"}

# If you use this anywhere (optional / legacy), keep it consistent:
#REGIONS = ["B", "A", "D", "C"]
SVJ_ORDER = ["0SVJ", "1SVJ", "2SVJ", "3PSVJ"]


# No remapping: read the same region name that you plot
ROOT_REGION_FOR = {r: r for r in PLOT_REGIONS}
# REGIONS not needed; keep only if other code references it
REGIONS = PLOT_REGIONS

SVJ_XLABELS = ["0", "1", "2", "3+"]
# Background stacking order (bottom → top)
BKG_ORDER = ["QCD", "TTJets", "WJets", "ZJets", "ST"]

PROC_LABEL = {
    "QCD": "QCD",
    "TTJets": r"$t\bar{t}$+jets",
    "WJets": r"$W\rightarrow l\nu$+jets",
    "ZJets": r"$Z\rightarrow\nu\nu$+jets",
    "ST": "Single top",
}

# Explicit colors so QCD is visually obvious (bottom big block)
# PROC_COLOR = {
#     "QCD": "#7a21dd",
#     "TTJets": "#9c9ca1",
#     "WJetsToLNu": "#5790fc",
#     "ZJetsToNuNu": "#e42536",
#     "ST": "#f89c20",
# }

PROC_COLOR = {
    "QCD": "#f89c20",        # was ST
    "TTJets": "#e42536",     # was Z
    "WJets": "#5790fc",      # stays middle
    "ZJets": "#9c9ca1",      # was TT
    "ST": "#7a21dd",         # was QCD
}

# Optional: hatches (CMS-style often uses hatches for readability in B/W)
PROC_HATCH = {
    "QCD": None,
    "TTJets": None,
    "WJets": None,
    "ZJets": None,
    "ST": None,
}


# Signals to overlay (must match ROOT histogram names exactly)
SIGNAL_PROCS = [
    "mMed1000_mDark20_rinv0p3_yukawa1",
    "mMed1500_mDark20_rinv0p3_yukawa1",
    "mMed2000_mDark20_rinv0p3_yukawa1",
    "mMed4000_mDark20_rinv0p3_yukawa1",
]

SIG_COLORS = {
    "mMed1000_mDark20_rinv0p3_yukawa1": "red",
    "mMed1500_mDark20_rinv0p3_yukawa1": "blue",
    "mMed2000_mDark20_rinv0p3_yukawa1": "darkgreen",
    "mMed4000_mDark20_rinv0p3_yukawa1": "purple",
}

SIG_LINESTYLE = "--"
SIG_LINEWIDTH = 2.5
SIG_SCALE = 1.0   # change if you want ×10, ×50, etc.



DATA_NAME = "data_obs"

# --------------------------
# Helpers
# --------------------------
def strip_cycle(name: str) -> str:
    return name.split(";")[0]

def write_abcd_compact_table_txt(year, yields, outdir):
    """
    Write compact ABCD yield table (summed over backgrounds),
    with clean rounding and no scientific notation.
    """
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, f"ABCD_yields_compact_{year}.txt")

    def fmt(x):
        if abs(x) < 1e-6:
            return "0.000"
        if x >= 100:
            return f"{x:.0f}"
        if x >= 10:
            return f"{x:.1f}"
        if x >= 1:
            return f"{x:.2f}"
        return f"{x:.3f}"

    def fmt_ratio(num, den):
        if abs(den) < 1e-12:
            return "N/A"
        return fmt(num / den)

    # Sum over backgrounds
    total = {svj: {reg: 0.0 for reg in PLOT_REGIONS} for svj in SVJ_ORDER}
    for proc in BKG_ORDER:
        for reg in PLOT_REGIONS:
            for svj in SVJ_ORDER:
                total[svj][reg] += yields.get(proc, {}).get(reg, {}).get(svj, 0.0)

    # Column totals
    col_totals = {reg: sum(total[svj][reg] for svj in SVJ_ORDER) for reg in PLOT_REGIONS}

    with open(outpath, "w") as f:
        # Header
        f.write(
            f"{'SVJ':<6}"
            f"{'A':>14}{'B':>14}{'C':>14}{'D':>14}"
            f"{'B/D':>14}{'C/D':>14}\n"
        )
        f.write("-" * 90 + "\n")

        # Rows
        for svj in SVJ_ORDER:
            f.write(
                f"{svj:<6}"
                f"{fmt(total[svj]['A']):>14}"
                f"{fmt(total[svj]['B']):>14}"
                f"{fmt(total[svj]['C']):>14}"
                f"{fmt(total[svj]['D']):>14}"
                f"{fmt_ratio(total[svj]['B'], total[svj]['D']):>14}"
                f"{fmt_ratio(total[svj]['C'], total[svj]['D']):>14}\n"
            )

        # Totals row
        f.write("-" * 90 + "\n")
        f.write(
            f"{'TOTAL':<6}"
            f"{fmt(col_totals['A']):>14}"
            f"{fmt(col_totals['B']):>14}"
            f"{fmt(col_totals['C']):>14}"
            f"{fmt(col_totals['D']):>14}"
            f"{fmt_ratio(col_totals['B'], col_totals['D']):>14}"
            f"{fmt_ratio(col_totals['C'], col_totals['D']):>14}\n"
        )

    print(f"[OK] wrote compact ABCD table {outpath}")



def print_yield_summary(year, yields):
    """
    yields[proc][region][svj] = yield
    """
    print("\n" + "="*90)
    print(f"BACKGROUND YIELDS SUMMARY — {year} (Run2)")
    print("="*90)

    total_year = 0.0

    for region in PLOT_REGIONS:
        print(f"\nRegion {region}:")
        print("-"*70)

        region_total = 0.0

        for svj in SVJ_ORDER:
            svj_total = 0.0
            line = f"  nSVJ={svj:>4} : "

            for proc in BKG_ORDER:
                y = yields.get(proc, {}).get(region, {}).get(svj, 0.0)
                svj_total += y
                line += f"{proc}={y:8.2f}  "

            region_total += svj_total
            line += f"|  TOTAL={svj_total:9.2f}"
            print(line)

        print(f"  --> Region {region} total = {region_total:10.2f}")
        total_year += region_total

    print("\n" + "-"*70)
    print(f"TOTAL BACKGROUND (all regions, all nSVJ) = {total_year:12.2f}")
    print("="*90 + "\n")

def detect_years(file_):
    years = set()
    for k in file_.keys():
        k0 = strip_cycle(k)
        m = re.search(r"Y(2016|2017|2018)_Run2$", k0)
        if m:
            years.add(m.group(1))
    return sorted(years)

def svj_dir_name(svj: str, year: str):
    return f"{svj}Y{year}_Run2"

def get_dir(file_, path: str):
    try:
        return file_[path]
    except Exception:
        return None

def list_dir_items(dir_):
    return [strip_cycle(k) for k in dir_.keys()]

def hist_integral_and_err(h, include_flow=True):
    """
    Integrate a TH1-like object read by uproot.
    Error is sqrt(sum(variances)) if available, else Poisson sqrt(N) fallback.
    """
    vals, _ = h.to_numpy(flow=include_flow)
    y = float(np.sum(vals))

    try:
        v = h.variances(flow=include_flow)
        if v is not None:
            e2 = float(np.sum(v))
            return y, math.sqrt(e2) if e2 >= 0 else 0.0
    except Exception:
        pass

    return y, math.sqrt(y) if y >= 0 else 0.0

def read_region_bkgs(file_, svj: str, year: str, region: str, include_flow=True, with_data=False):
    """
    Returns:
      bkg_yields: dict(proc -> (yield, err)) or None if missing
      data: (yield, err) or None (only when with_data=True)
    """
    base = svj_dir_name(svj, year)
    dreg = get_dir(file_, f"{base}/{region}")
    if dreg is None:
        return None, None

    keys = list_dir_items(dreg)

    bkg = {}
    for proc in BKG_ORDER:
        if proc in keys:
            y, e = hist_integral_and_err(dreg[proc], include_flow=include_flow)
            bkg[proc] = (y, e)
        else:
            bkg[proc] = (0.0, 0.0)

    data = None
    if with_data and (DATA_NAME in keys):
        y, e = hist_integral_and_err(dreg[DATA_NAME], include_flow=include_flow)
        data = (y, e)

    return bkg, data

def build_x_layout():
    """
    Continuous layout:
      A: 0 1 2 3 | B: 4 5 6 7 | C: 8 9 10 11 | D: 12 13 14 15
    No gaps between bars.
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

    #return np.array(x), xticks, xticklabels, region_centers, region_boundaries
    return np.array(x, dtype=float), xticks, xticklabels, region_centers, region_boundaries

def step_band_from_bins(xpos, y, yerr, width=0.85, eps=1e-6):
    left = xpos - width / 2.0
    right = xpos + width / 2.0

    xe = np.empty(2 * len(xpos))
    yeu = np.empty(2 * len(xpos))
    yel = np.empty(2 * len(xpos))

    for i in range(len(xpos)):
        xe[2*i] = left[i]
        xe[2*i + 1] = right[i]

        up = y[i] + yerr[i]
        lo = y[i] - yerr[i]

        # IMPORTANT: never hit 0 on log axis
        lo = max(eps, lo)

        yeu[2*i] = up
        yeu[2*i + 1] = up
        yel[2*i] = lo
        yel[2*i + 1] = lo

    return xe, yel, yeu
    


# --------------------------
# Plotting
# --------------------------
def plot_year(file_path, year, outdir=".", include_flow=True, with_data=False,
              lumi_text=None, com_text="13", prelim=True, years_to_sum=None,
              output_label=None):
    os.makedirs(outdir, exist_ok=True)
    years_to_read = years_to_sum if years_to_sum is not None else [year]
    plot_label = output_label if output_label is not None else year

    # CMS style
    hep.style.use("CMS")
    plt.rcParams["figure.dpi"] = 150

    x, xticks, xticklabels, region_centers, region_boundaries = build_x_layout()
    nsvj = len(SVJ_ORDER)

    # 2-panel layout like CMS: main + ratio
    # Figure: 2-panel only if with_data, otherwise single panel
    if with_data:
        fig, (ax, rax) = plt.subplots(
            2, 1, figsize=(14.0, 9.5),
            gridspec_kw={"height_ratios": (3.4, 1.0), "hspace": 0.03},
            sharex=True,
            constrained_layout=True
        )
    else:
        fig, ax = plt.subplots(
            1, 1, figsize=(14.0, 7.5),
            constrained_layout=True
        )
        rax = None

 
    # Collect yields
    proc_y = {p: np.zeros_like(x, dtype=float) for p in BKG_ORDER}
    proc_var = {p: np.zeros_like(x, dtype=float) for p in BKG_ORDER}
    sig_y = {s: np.zeros_like(x, dtype=float) for s in SIGNAL_PROCS}

    data_y = np.zeros_like(x, dtype=float)
    data_var = np.zeros_like(x, dtype=float)
    data_seen = np.zeros_like(x, dtype=bool)


    with uproot.open(file_path) as f:
        for year_to_read in years_to_read:
            for r_i, reg_plot in enumerate(PLOT_REGIONS):
                reg_root = ROOT_REGION_FOR[reg_plot]

                for s_i, svj in enumerate(SVJ_ORDER):
                    idx = r_i * nsvj + s_i

                    # --- backgrounds ---
                    bkg, data = read_region_bkgs(
                        f, svj=svj, year=year_to_read, region=reg_root,
                        include_flow=include_flow, with_data=with_data
                    )
                    if bkg is None:
                        continue

                    for p in BKG_ORDER:
                        proc_y[p][idx] += bkg[p][0]
                        proc_var[p][idx] += bkg[p][1] ** 2

                    # --- signals ---
                    for sig in SIGNAL_PROCS:
                        sig_path = f"{svj}Y{year_to_read}_Run2/{reg_root}/{sig}"
                        if sig_path in f:
                            hsig = f[sig_path]
                            vals, _ = hsig.to_numpy(flow=include_flow)
                            sig_y[sig][idx] += vals.sum()

                    # --- data ---

                    if with_data and (data is not None):
                        data_y[idx] += data[0]
                        data_var[idx] += data[1] ** 2
                        data_seen[idx] = True

    proc_e = {p: np.sqrt(proc_var[p]) for p in BKG_ORDER}
    data_e = np.sqrt(data_var)
    data_y[~data_seen] = np.nan
    data_e[~data_seen] = np.nan


    # ---------------------------------------
    # DEBUG: check stacking numerically
    # First bin = Region A, nSVJ = 0
    # ---------------------------------------
    print(
        f"DEBUG {plot_label} (Region A, nSVJ=0):",
        {p: proc_y[p][0] for p in reversed(BKG_ORDER)},
        "TOTAL =", sum(proc_y[p][0] for p in BKG_ORDER)
    )


    # Stack (QCD first => bottom)
    # Stack (QCD first => bottom)
    bottoms = np.zeros_like(x)
    # Bars: draw on top
    bottoms = np.zeros_like(x)
    for p in reversed(BKG_ORDER):
        ax.bar(
            x, proc_y[p],
            bottom=bottoms,
            width=1.0,
            label=PROC_LABEL.get(p, p),
            color=PROC_COLOR.get(p, None),
            edgecolor="black",
            linewidth=0.2,
            hatch=PROC_HATCH.get(p, None),
            zorder=1,          # <-- was 1
        )
        bottoms += proc_y[p]

    

    # -------------------------------------------------
    # DEBUG #2: verify the stacked total equals sum(MC)
    # bottoms is now the total MC after stacking
    # First bin = Region A, nSVJ=0 (index 0)
    # -------------------------------------------------
    print(
        f"DEBUG {plot_label} stack check (Region A, nSVJ=0): "
        f"QCD={proc_y['QCD'][0]:.2f}  "
        f"MC_total_from_stack={bottoms[0]:.2f}  "
        f"MC_total_from_sum={sum(proc_y[p][0] for p in BKG_ORDER):.2f}"
    )


    # Total MC uncertainty band (from sum of per-process variances approximated as sum of per-process errors^2)
    mc_tot = np.zeros_like(x)
    mc_var = np.zeros_like(x)
    for p in reversed(BKG_ORDER):
        mc_tot += proc_y[p]
        mc_var += proc_e[p] ** 2
    mc_err = np.sqrt(mc_var)

    #xe, ylo, yhi = step_band_from_bins(x, mc_tot, mc_err, width=0.85)
    xe, ylo, yhi = step_band_from_bins(x, mc_tot, mc_err, width=1.0)

    # ax.fill_between(
    #     xe, ylo, yhi,
    #     step="pre",
    #     alpha=0.3,
    #     label="MC unc."
    # )
    # Draw MC uncertainty as a hatched band WITHOUT filling, so it doesn't wash out QCD
    # MC unc band: draw behind
    ax.fill_between(
        xe, ylo, yhi,
        step="pre",
        facecolor="none",
        edgecolor="0.6",
        linewidth=0.0,
        hatch="////",
        label="MC unc.",
        zorder=1              # <-- keep low
    )

    # -------------------------------------------------
    # Signal overlays (dashed lines)  <-- PUT IT HERE
    # -------------------------------------------------
    for sig in SIGNAL_PROCS:
        ysig = SIG_SCALE * sig_y[sig]
        ax.step(
            np.r_[x - 0.5, x[-1] + 0.5],
            np.r_[ysig, ysig[-1]],
            where="post",
            color=SIG_COLORS.get(sig, "red"),
            linestyle=SIG_LINESTYLE,
            linewidth=SIG_LINEWIDTH,
            label=sig if SIG_SCALE == 1.0 else f"{sig} × {SIG_SCALE:g}",
            zorder=5,
        )


    # Optional data overlay (OFF by default)
    if with_data:
        m = np.isfinite(data_y)
        ax.errorbar(
            x[m], data_y[m], yerr=data_e[m],
            fmt="o", color="black", markersize=4,
            linewidth=1.0, capsize=0,
            label="Data"
        )

        # Ratio: Data/Sim
        ratio = np.full_like(x, np.nan, dtype=float)
        ratio_err = np.full_like(x, np.nan, dtype=float)

        denom = mc_tot
        m2 = m & (denom > 0)
        ratio[m2] = data_y[m2] / denom[m2]
        # propagate data stat only for now (you can add mc_unc later)
        ratio_err[m2] = data_e[m2] / denom[m2]

        rax.errorbar(
            x[m2], ratio[m2], yerr=ratio_err[m2],
            fmt="o", color="black", markersize=4,
            linewidth=1.0, capsize=0
        )
        rax.set_ylabel("Data/Sim")
    # else:
    #     # Keep the panel for CMS-like layout, but no data points
    #     rax.set_ylabel("Data/Sim")

    # Ratio panel: unity line + MC unc band (as relative)
    if with_data:
        rax.axhline(1.0, color="black", linewidth=1.0)
        mct = mc_tot.copy()
        rel_lo = np.ones_like(mct)
        rel_hi = np.ones_like(mct)
        ok = mct > 0
        rel_lo[ok] = np.maximum(0.0, (mct[ok] - mc_err[ok]) / mct[ok])
        rel_hi[ok] = (mct[ok] + mc_err[ok]) / mct[ok]

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
            zorder=2
        )

        rax.set_ylim(0.0, 2.0)


    # Region separators and labels
    for xb in region_boundaries:
        ax.axvline(xb, color="black", linewidth=1.6)
        if with_data:
            rax.axvline(xb, color="black", linewidth=1.6)


    # Place A/B/C/D labels inside the main axis near top
    ymax = np.nanmax(mc_tot) if np.nanmax(mc_tot) > 0 else 1.0
    for reg, xc in region_centers.items():
        ax.text(xc, ymax * 0.35, reg, ha="center", va="bottom", fontsize=15, fontweight="bold")

    # Axes formatting
    ax.set_yscale("log")
    ax.set_ylabel("Events")

    ax.set_ylim(10e-3, ymax * 50.0 if ymax > 0 else 10.0)
    ax.set_xlim(-0.5, len(x) - 0.5)


    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    ax.set_xlabel("nSVJ")

    if with_data:
        rax.set_xlabel("nSVJ")


    # Grid similar to CMS plots
    ax.grid(True, which="both", axis="y", linestyle=":", linewidth=0.8)
    if with_data:
        rax.grid(True, which="both", axis="y", linestyle=":", linewidth=0.8)

    # Legend: keep compact and high, like screenshot
    ax.legend(
        loc="upper left",
        ncol=2,
        frameon=False,
        fontsize=10,
        handlelength=1.6,
        columnspacing=1.0
    )

    # CMS label (mplhep)
    # If you want the exact "CMS Preliminary 41.5 fb^{-1} (13 TeV)" formatting, set lumi_text accordingly.
    if lumi_text is None:
        # reasonable defaults if you don't pass lumi
        lumi_text = {"2016": "36.31", "2017": "42.07", "2018": "59.56", "Run2_Combined": "137.94"}.get(plot_label, "")
        #lumi_text = {"2016": "35.9 ", "2017": "41.5 ", "2018": "59.7 "}.get(year, "")
        lumi_text = {"2016": "36.31", "2017": "42.07", "2018": "59.56", "Run2_Combined": "137.94"}.get(plot_label, "")
    hep.cms.label(
        ax=ax,
        label="Preliminary" if prelim else "",
        data=with_data,
        lumi=lumi_text,
        com=com_text
    )

    #fig.tight_layout()

    suffix = "_with_data" if with_data else ""
    outpath = os.path.join(outdir, f"ABCD_yields_{plot_label}_CMS{suffix}.pdf")
    fig.savefig(outpath)
    plt.close(fig)
    print(f"[OK] wrote {outpath}")


def hist_arrays_no_flow(h):
    """
    Return bin values/edges/variances for shape plots.
    The DNN score x-axis is fixed to visible bins, so under/overflow are not
    folded into the plotted distribution.
    """
    vals, edges = h.to_numpy(flow=False)
    vals = np.asarray(vals, dtype=float)
    edges = np.asarray(edges, dtype=float)

    try:
        variances = h.variances(flow=False)
    except Exception:
        variances = None
    if variances is None:
        variances = np.abs(vals)
    else:
        variances = np.asarray(variances, dtype=float)

    return vals, edges, variances


def add_shape(target_vals, target_vars, target_edges, vals, variances, edges, label):
    if target_edges is None:
        target_edges = edges.copy()
        target_vals = np.zeros_like(vals, dtype=float)
        target_vars = np.zeros_like(vals, dtype=float)
    elif len(target_edges) != len(edges) or not np.allclose(target_edges, edges):
        raise RuntimeError(f"DNN score binning mismatch while adding {label}")

    target_vals += vals
    target_vars += variances
    return target_vals, target_vars, target_edges


def read_dnn_score_shapes(file_path, years_to_read, with_data=True):
    """
    Sum the combine-file TH1 shapes over ABCD regions and nSVJ categories.
    The same category/process histograms used for the ABCD yield integrals are
    treated here as the DNN-score distribution.
    """
    proc_vals = {p: None for p in BKG_ORDER}
    proc_vars = {p: None for p in BKG_ORDER}
    signal_vals = {s: None for s in SIGNAL_PROCS}
    signal_vars = {s: None for s in SIGNAL_PROCS}
    data_vals = None
    data_vars = None
    edges = None
    n_added = 0

    with uproot.open(file_path) as f:
        for year in years_to_read:
            for reg_plot in PLOT_REGIONS:
                reg_root = ROOT_REGION_FOR[reg_plot]
                for svj in SVJ_ORDER:
                    base = f"{svj}Y{year}_Run2/{reg_root}"

                    for proc in BKG_ORDER:
                        hist_path = f"{base}/{proc}"
                        if hist_path not in f:
                            continue
                        vals, hist_edges, variances = hist_arrays_no_flow(f[hist_path])
                        proc_vals[proc], proc_vars[proc], edges = add_shape(
                            proc_vals[proc], proc_vars[proc], edges,
                            vals, variances, hist_edges, hist_path
                        )
                        n_added += 1

                    for sig in SIGNAL_PROCS:
                        hist_path = f"{base}/{sig}"
                        if hist_path not in f:
                            continue
                        vals, hist_edges, variances = hist_arrays_no_flow(f[hist_path])
                        signal_vals[sig], signal_vars[sig], edges = add_shape(
                            signal_vals[sig], signal_vars[sig], edges,
                            vals, variances, hist_edges, hist_path
                        )

                    if with_data:
                        hist_path = f"{base}/{DATA_NAME}"
                        if hist_path in f:
                            vals, hist_edges, variances = hist_arrays_no_flow(f[hist_path])
                            data_vals, data_vars, edges = add_shape(
                                data_vals, data_vars, edges,
                                vals, variances, hist_edges, hist_path
                            )

    if edges is None or n_added == 0:
        raise RuntimeError(
            "No DNN score histograms found. Expected paths like "
            "0SVJY2018_Run2/A/QCD in the combine file."
        )

    zero = np.zeros(len(edges) - 1, dtype=float)
    for proc in BKG_ORDER:
        if proc_vals[proc] is None:
            proc_vals[proc] = zero.copy()
            proc_vars[proc] = zero.copy()
    for sig in SIGNAL_PROCS:
        if signal_vals[sig] is None:
            signal_vals[sig] = zero.copy()
            signal_vars[sig] = zero.copy()

    return edges, proc_vals, proc_vars, signal_vals, signal_vars, data_vals, data_vars


def step_values(edges, values):
    return edges, np.r_[values, values[-1]]


def add_unique_legend(ax, handles, labels, **kwargs):
    seen = set()
    uniq_handles = []
    uniq_labels = []
    for handle, label in zip(handles, labels):
        if label in seen:
            continue
        seen.add(label)
        uniq_handles.append(handle)
        uniq_labels.append(label)
    ax.legend(uniq_handles, uniq_labels, **kwargs)


def signal_display_label(sig):
    match = re.search(r"mMed(\d+)_mDark(\d+)_rinv([0-9p]+)_yukawa([0-9p]+)", sig)
    if not match:
        return sig
    m_med, m_dark, rinv, yukawa = match.groups()
    rinv = rinv.replace("p", ".")
    yukawa = yukawa.replace("p", ".")
    return (
        rf"$m_{{\Phi}}={m_med}$ GeV, "
        rf"$m_{{dark}}={m_dark}$ GeV, "
        rf"$r_{{inv}}={rinv}$, "
        rf"$\lambda={yukawa}$"
    )


def plot_dnn_score(file_path, year, outdir=".", with_data=False, lumi_text=None,
                   com_text="13", prelim=True, years_to_sum=None,
                   output_label=None):
    years_to_read = years_to_sum if years_to_sum is not None else [year]
    plot_label = output_label if output_label is not None else year
    subdir = "with_data" if with_data else "mc_only"
    plot_dir = os.path.join(outdir, "DNN_score", subdir)
    os.makedirs(plot_dir, exist_ok=True)

    edges, proc_vals, proc_vars, signal_vals, _, data_vals, data_vars = read_dnn_score_shapes(
        file_path, years_to_read, with_data=with_data
    )

    total_mc = sum(proc_vals[p] for p in BKG_ORDER)
    total_mc_integral = float(np.sum(total_mc))
    if total_mc_integral <= 0:
        raise RuntimeError(f"No positive MC yield available for DNN score plot {plot_label}")

    proc_plot = {p: proc_vals[p] / total_mc_integral for p in BKG_ORDER}
    mc_total_plot = total_mc / total_mc_integral

    data_plot = None
    data_err = None
    if with_data and data_vals is not None:
        data_integral = float(np.sum(data_vals))
        if data_integral > 0:
            data_plot = data_vals / data_integral
            data_err = np.sqrt(np.maximum(data_vars, 0.0)) / data_integral

    signal_plot = {}
    for sig in SIGNAL_PROCS:
        integral = float(np.sum(signal_vals[sig]))
        if integral > 0:
            signal_plot[sig] = signal_vals[sig] / integral

    hep.style.use("CMS")
    plt.rcParams["figure.dpi"] = 150
    fig, ax = plt.subplots(1, 1, figsize=(7.2, 7.0), constrained_layout=True)

    bottoms = np.zeros_like(mc_total_plot)
    stack_handles = []
    stack_labels = []
    for proc in reversed(BKG_ORDER):
        widths = np.diff(edges)
        bars = ax.bar(
            edges[:-1], proc_plot[proc],
            width=widths,
            bottom=bottoms,
            align="edge",
            color=PROC_COLOR.get(proc, None),
            edgecolor="black",
            linewidth=0.2,
            label=PROC_LABEL.get(proc, proc),
            zorder=1,
        )
        bottoms += proc_plot[proc]
        stack_handles.append(bars[0])
        stack_labels.append(PROC_LABEL.get(proc, proc))

    signal_handles = []
    signal_labels = []
    for sig in SIGNAL_PROCS:
        if sig not in signal_plot:
            continue
        x_step, y_step = step_values(edges, signal_plot[sig])
        line = ax.step(
            x_step, y_step,
            where="post",
            color=SIG_COLORS.get(sig, "black"),
            linestyle=SIG_LINESTYLE,
            linewidth=2.0,
            label=signal_display_label(sig),
            zorder=4,
        )[0]
        signal_handles.append(line)
        signal_labels.append(signal_display_label(sig))

    data_handle = None
    data_label = None
    if with_data and data_plot is not None:
        centers = 0.5 * (edges[:-1] + edges[1:])
        data_handle = ax.errorbar(
            centers, data_plot, yerr=data_err,
            fmt="o", color="black", markersize=4,
            linewidth=1.0, capsize=0,
            label="Data",
            zorder=5,
        )
        data_label = "Data"

    ax.set_yscale("log")
    positive = [np.max(mc_total_plot) if np.any(mc_total_plot > 0) else 0.0]
    positive.extend(np.max(v) for v in signal_plot.values() if np.any(v > 0))
    if data_plot is not None and np.any(data_plot > 0):
        positive.append(np.max(data_plot))
    ymax = max(positive) if positive else 1.0
    ax.set_ylim(1e-5, max(1.0, ymax * 250.0))
    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel("DNN score")
    ax.set_ylabel("Arbitrary units")
    ax.grid(True, which="both", axis="y", linestyle=":", linewidth=0.8)

    if lumi_text is None:
        lumi_text = {
            "2016": "36.31",
            "2017": "42.07",
            "2018": "59.56",
            "Run2_Combined": "137.94",
        }.get(plot_label, "")
    hep.cms.label(
        ax=ax,
        label="Preliminary" if prelim else "",
        data=with_data,
        lumi=lumi_text,
        com=com_text
    )

    handles = []
    labels = []
    if data_handle is not None:
        handles.append(data_handle)
        labels.append(data_label)
    handles.extend(stack_handles)
    labels.extend(stack_labels)
    handles.extend(signal_handles)
    labels.extend(signal_labels)
    add_unique_legend(
        ax, handles, labels,
        loc="upper left",
        bbox_to_anchor=(0.03, 0.91),
        ncol=2,
        frameon=False,
        fontsize=9,
        handlelength=1.7,
        columnspacing=0.9,
    )

    suffix = "_with_data" if with_data else "_mc_only"
    outbase = os.path.join(plot_dir, f"DNN_score_{plot_label}_CMS{suffix}")
    fig.savefig(f"{outbase}.pdf")
    fig.savefig(f"{outbase}.png")
    plt.close(fig)
    print(f"[OK] wrote {outbase}.pdf")
    print(f"[OK] wrote {outbase}.png")


def write_detailed_yield_table_txt(year, yields, outdir):
    """
    Write detailed ABCD yield table separated by background process,
    matching the requested ROOT-style text format.
    """
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, f"ABCD_yields_detailed_{year}.txt")

    # Labels formatted exactly like your example
    table_proc_labels = {
        "ST": "Single top",
        "TTJets": "t#bar{t}",
        "ZJets": "Z+jets",
        "WJets": "W+jets",
        "QCD": "QCD"
    }
    
    # Process order matching your example
    table_proc_order = ["ST", "TTJets", "ZJets", "WJets", "QCD"]

    def fmt_ratio(num, den):
        if abs(den) < 1e-12:
            return "N/A"
        return f"{num / den:.6f}"

    with open(outpath, "w") as f:
        # Header
        f.write(
            f"{'SVJ':<8} {'Component':<26} "
            f"{'A':>12} {'B':>12} {'C':>12} {'D':>12} "
            f"{'TOTAL':>12} {'B/D':>12} {'C/D':>12}\n"
        )
        f.write("-" * 116 + "\n")

        # Rows
        for svj in SVJ_ORDER:
            for proc in table_proc_order:
                # Safely get yields for each region, defaulting to 0.0
                yA = yields.get(proc, {}).get("A", {}).get(svj, 0.0)
                yB = yields.get(proc, {}).get("B", {}).get(svj, 0.0)
                yC = yields.get(proc, {}).get("C", {}).get(svj, 0.0)
                yD = yields.get(proc, {}).get("D", {}).get(svj, 0.0)
                ytot = yA + yB + yC + yD
                
                label = table_proc_labels.get(proc, proc)
                
                f.write(
                    f"{svj:<8} {label:<26} "
                    f"{yA:>12.6f} {yB:>12.6f} {yC:>12.6f} {yD:>12.6f} "
                    f"{ytot:>12.6f} {fmt_ratio(yB, yD):>12} {fmt_ratio(yC, yD):>12}\n"
                )

    print(f"[OK] wrote detailed ABCD table {outpath}")

def compute_background_fraction_summary(yields):
    """
    Sum ABCD yields into the same background-fraction inputs used by the
    LaTeX table and the per-nSVJ pie charts.
    """
    total_proc_svj = {p: {s: 0.0 for s in SVJ_ORDER} for p in BKG_ORDER}
    for p in BKG_ORDER:
        for s in SVJ_ORDER:
            total_proc_svj[p][s] = sum(
                yields.get(p, {}).get(r, {}).get(s, 0.0) for r in PLOT_REGIONS
            )

    total_svj = {
        s: sum(total_proc_svj[p][s] for p in BKG_ORDER) for s in SVJ_ORDER
    }
    total_proc_incl = {
        p: sum(total_proc_svj[p][s] for s in SVJ_ORDER) for p in BKG_ORDER
    }
    total_bkg_incl = sum(total_svj[s] for s in SVJ_ORDER)

    return total_proc_svj, total_svj, total_proc_incl, total_bkg_incl

# combinehistplotter's yield dict uses "WJets"/"ZJets"; Figure2 uses the
# longer "WJetsToLNu"/"ZJetsToNuNu" keys for the same processes. This maps
# between them so pie colors are pulled from the single shared palette.
PIE_TO_FIGURE2_PROC = {
    "QCD": "QCD",
    "TTJets": "TTJets",
    "WJets": "WJetsToLNu",
    "ZJets": "ZJetsToNuNu",
    "ST": "ST",
}


def write_fraction_pie_charts(year_label, yields, outdir, lumi_text=None,
                              com_text="13", prelim=True,
                              draw_cms_heading=True,
                              pie_text_fontsize=17,
                              tagger_label="PN"):
    """
    Write one combined figure with a background-fraction pie chart for each
    nSVJ bin (2x2 grid), sharing a single Figure2-style legend.
    """
    os.makedirs(outdir, exist_ok=True)

    total_proc_svj, total_svj, _, _ = compute_background_fraction_summary(yields)

    # Ordered to match Figure2_makerusingskims.LEGEND_ORDER so the same
    # process always gets the same color/position across every supplementary
    # plot, not just within this figure.
    pie_proc_order = ["QCD", "TTJets", "WJets", "ZJets", "ST"]
    tagger_latex = r"\mathrm{" + str(tagger_label) + "}"
    svj_titles = [
        rf"$n_{{\mathrm{{SVJ}}}}^{{{tagger_latex}}} = 0$",
        rf"$n_{{\mathrm{{SVJ}}}}^{{{tagger_latex}}} = 1$",
        rf"$n_{{\mathrm{{SVJ}}}}^{{{tagger_latex}}} = 2$",
        rf"$n_{{\mathrm{{SVJ}}}}^{{{tagger_latex}}} \geq 3$",
    ]

    hep.style.use("CMS")
    plt.rcParams["figure.dpi"] = 150

    def autopct_fmt(pct):
        return f"{pct:.0f}%" if pct >= 3.0 else ""

    if lumi_text is None:
        lumi_text = {
            "2016": "36.31",
            "2017": "42.07",
            "2018": "59.56",
            "Run2_Combined": "137.94",
        }.get(year_label, "")
    lumi_label = f"{lumi_text} fb$^{{-1}}$ ({com_text.strip()} TeV)" if lumi_text else f"({com_text.strip()} TeV)"
    cms_modifier = "Simulation Preliminary" if prelim else "Simulation"

    colors = [f2.PROC_COLOR_HEX[PIE_TO_FIGURE2_PROC[p]] for p in pie_proc_order]
    legend_labels = [PROC_LABEL.get(p, p) for p in pie_proc_order]

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 10.8))
    fig.subplots_adjust(left=0.03, right=0.80, top=0.84, bottom=0.03, wspace=0.05, hspace=0.25)

    for ax, svj, title in zip(axes.flat, SVJ_ORDER, svj_titles):
        values = [total_proc_svj[p][svj] for p in pie_proc_order]

        if total_svj[svj] > 0:
            _, _, autotexts = ax.pie(
                values,
                colors=colors,
                startangle=90,
                counterclock=False,
                autopct=autopct_fmt,
                pctdistance=1.14,
                radius=1.16,
                wedgeprops={"edgecolor": "black", "linewidth": 0.6},
                textprops={"fontsize": pie_text_fontsize, "color": "black"},
            )
            for txt in autotexts:
                txt.set_fontsize(pie_text_fontsize)
                txt.set_fontweight("bold")
        else:
            ax.text(0.5, 0.5, "No background", ha="center", va="center",
                    transform=ax.transAxes)

        ax.set_aspect("equal")
        ax.set_title(title, fontsize=22, pad=14)

    # One shared legend for the whole figure, styled like Figure2's ROOT
    # legends: single column, solid white box, no border.
    legend_handles = [
        mpatches.Patch(facecolor=color, edgecolor="black", linewidth=0.6, label=label)
        for color, label in zip(colors, legend_labels)
    ]
    legend = fig.legend(
        handles=legend_handles,
        loc="center left",
        ncol=1,
        frameon=True,
        fontsize=17,
        bbox_to_anchor=(0.80, 0.5),
        handlelength=1.5,
        handletextpad=0.8,
        labelspacing=1.2,
    )
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_edgecolor("none")

    if draw_cms_heading:
        fig.text(0.03, 0.975, "CMS", ha="left", va="top",
                 fontsize=32, fontweight="bold")
        fig.text(0.145, 0.975, cms_modifier, ha="left", va="top",
                 fontsize=22, fontstyle="italic")
        fig.text(0.80, 0.975, lumi_label, ha="right", va="top",
                 fontsize=18)

    outbase = os.path.join(outdir, f"ABCD_bkg_fraction_pie_{year_label}_CMS")
    fig.savefig(f"{outbase}.pdf", bbox_inches="tight")
    fig.savefig(f"{outbase}.png", bbox_inches="tight")
    plt.close(fig)

    print("[OK] wrote combined background fraction pie chart:")
    print(f"     {outbase}.pdf")
    return f"{outbase}.pdf"

def write_latex_fraction_tables(year_label, yields, outdir):
    """
    Calculates background fractions and efficiencies across nSVJ bins
    and writes them to a LaTeX file.
    """
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, f"ABCD_latex_tables_{year_label}.tex")
    display_year_label = year_label.replace("_", " ")

    # LaTeX mapping for process names
    tex_proc_map = {
        "QCD": "QCD multijet",
        "TTJets": "\\ttjets",
        "ZJets": "\\zjets",
        "WJets": "\\wjets",
        "ST": "Single top"
    }
    # Match your desired table order
    tex_proc_order = ["QCD", "TTJets", "ZJets", "WJets", "ST"]
    
    total_proc_svj, total_svj, total_proc_incl, total_bkg_incl = (
        compute_background_fraction_summary(yields)
    )

    with open(outpath, "w") as f:
        # --- TABLE 1: Background Fractions ---
        f.write("\\begin{table}[htbp]\n")
        f.write("\\centering\n")
        f.write("\\begin{tabular}{|c|c|c|c|c|c|} \n")
        f.write("    \\hline\n")
        f.write("    \\multirow{2}{*}{Background} & \\multicolumn{5}{c|}{Fraction [\\%]} \\\\ \n")
        f.write("    \\cline{2-6}\n")
        f.write("     & inclusive & $\\nsvjpn = 0$ & $\\nsvjpn = 1$ & $\\nsvjpn = 2$ & $\\nsvjpn \\geq 3$ \\\\ \n")
        f.write("    \\hline\n")

        for proc in tex_proc_order:
            name = tex_proc_map.get(proc, proc)
            
            # Inclusive fraction for this proc
            incl_frac = (total_proc_incl[proc] / total_bkg_incl * 100.0) if total_bkg_incl > 0 else 0.0
            
            # Fractions per SVJ bin
            svj_fracs = []
            for s in SVJ_ORDER:
                num = total_proc_svj[proc][s]
                den = total_svj[s]
                frac = (num / den * 100.0) if den > 0 else 0.0
                svj_fracs.append(frac)
            
            # Write row
            row_str = f"    {name:<15} & {incl_frac:>6.3f} & {svj_fracs[0]:>6.3f} & {svj_fracs[1]:>6.3f} & {svj_fracs[2]:>6.3f} & {svj_fracs[3]:>6.3f} \\\\\n"
            f.write(row_str)

        f.write("    \\hline\n")
        f.write("\\end{tabular}\n")
        f.write(f"\\caption{{Fraction of the different backgrounds for all events (first column) and in bins of number of ParticleNet-tagged SVJs (last columns) for the DNN introduced in Section~\\ref{{sec:closure_nsvjs}}. ({display_year_label})}}\n")
        f.write("\\label{table:abcd_bkg_fraction_nsvjpn}\n")
        f.write("\\end{table}\n\n")

        # --- FIGURE: Background Fractions as Pie Charts (one combined figure) ---
        pie_file = f"ABCD_bkg_fraction_pie_{year_label}_CMS.pdf"
        f.write("\\begin{figure}[htbp]\n")
        f.write("\\centering\n")
        f.write(f"\\includegraphics[width=0.85\\textwidth]{{{pie_file}}}\n")
        f.write(f"\\caption{{Background composition in bins of $\\nsvjpn$ for {display_year_label}.}}\n")
        f.write("\\label{fig:abcd_bkg_fraction_pie_nsvjpn}\n")
        f.write("\\end{figure}\n\n")

        # --- TABLE 2: Background Efficiency ---
        f.write("\\begin{table}[htbp]\n")
        f.write("\\centering\n")
        f.write("\\begin{tabular}{|c|c|} \n")
        f.write("    \\hline\n")
        f.write("    Category           & Background efficiency [\\%] \\\\ \n")
        f.write("    \\hline\n")

        svj_labels = ["$\\nsvjpn = 0$", "$\\nsvjpn = 1$", "$\\nsvjpn = 2$", "$\\nsvjpn \\geq  3$"]
        for s, label in zip(SVJ_ORDER, svj_labels):
            eff = (total_svj[s] / total_bkg_incl * 100.0) if total_bkg_incl > 0 else 0.0
            f.write(f"    {label:<18} & {eff:>6.3f} \\\\\n")

        f.write("    \\hline\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")

    print(f"[OK] wrote LaTeX tables to {outpath}")


# =============================================================================
# Signal region-A acceptance (fraction of each signal landing in region A,
# i.e. the ABCD signal region), inclusive and split by nSVJ bin.
# =============================================================================

# The 10-point mMed/rinv reference grid already used for the DNN/PNet AUC
# scans elsewhere in this analysis (supplementary_material_plotter.py's
# DNN_TRAINING_SAMPLES / PNET_TRAINING_SAMPLES): a mass scan at rinv=0.3 plus
# a rinv scan at mMed=2000, all at mDark=20, yukawa=1.
REFERENCE_SIGNALS = [
    ("mMed600_mDark20_rinv0p3_yukawa1", r"$m_{\Phi}=600$ GeV, $r_{\mathrm{inv}}=0.3$"),
    ("mMed800_mDark20_rinv0p3_yukawa1", r"$m_{\Phi}=800$ GeV, $r_{\mathrm{inv}}=0.3$"),
    ("mMed1000_mDark20_rinv0p3_yukawa1", r"$m_{\Phi}=1000$ GeV, $r_{\mathrm{inv}}=0.3$"),
    ("mMed1500_mDark20_rinv0p3_yukawa1", r"$m_{\Phi}=1500$ GeV, $r_{\mathrm{inv}}=0.3$"),
    ("mMed2000_mDark20_rinv0p1_yukawa1", r"$m_{\Phi}=2000$ GeV, $r_{\mathrm{inv}}=0.1$"),
    ("mMed2000_mDark20_rinv0p3_yukawa1", r"$m_{\Phi}=2000$ GeV, $r_{\mathrm{inv}}=0.3$"),
    ("mMed2000_mDark20_rinv0p5_yukawa1", r"$m_{\Phi}=2000$ GeV, $r_{\mathrm{inv}}=0.5$"),
    ("mMed2000_mDark20_rinv0p7_yukawa1", r"$m_{\Phi}=2000$ GeV, $r_{\mathrm{inv}}=0.7$"),
    ("mMed3000_mDark20_rinv0p3_yukawa1", r"$m_{\Phi}=3000$ GeV, $r_{\mathrm{inv}}=0.3$"),
    ("mMed4000_mDark20_rinv0p3_yukawa1", r"$m_{\Phi}=4000$ GeV, $r_{\mathrm{inv}}=0.3$"),
]

SIGNAL_REGION = "A"


def read_region_signal(file_, svj: str, year: str, region: str, signal_name: str) -> float:
    """Yield of one signal process in one (svj, year, region) cell, or 0.0 if missing."""
    base = svj_dir_name(svj, year)
    dreg = get_dir(file_, f"{base}/{region}")
    if dreg is None:
        return 0.0
    keys = list_dir_items(dreg)
    if signal_name not in keys:
        return 0.0
    y, _ = hist_integral_and_err(dreg[signal_name])
    return y


def compute_signal_region_a_acceptance(file_path: str, years, signals=None):
    """
    For each signal, for each nSVJ bin (and inclusive), the fraction of
    Run2-combined events landing in the ABCD signal region A:
        N_A / (N_A + N_B + N_C + N_D)
    All preselected events fall into exactly one of the four ABCD regions,
    so this denominator is the total selected signal yield for that
    (nSVJ bin / inclusive) slice.

    Returns: {signal_name: {"0SVJ": frac, ..., "Inclusive": frac}}
    """
    if signals is None:
        signals = [name for name, _ in REFERENCE_SIGNALS]

    results = {}
    with uproot.open(file_path) as f:
        for signal_name in signals:
            per_svj = {}
            region_totals_incl = {r: 0.0 for r in PLOT_REGIONS}
            for svj in SVJ_ORDER:
                region_totals = {r: 0.0 for r in PLOT_REGIONS}
                for year in years:
                    for region in PLOT_REGIONS:
                        y = read_region_signal(f, svj, year, region, signal_name)
                        region_totals[region] += y
                        region_totals_incl[region] += y
                denom = sum(region_totals.values())
                per_svj[svj] = (region_totals[SIGNAL_REGION] / denom) if denom > 0 else 0.0
            denom_incl = sum(region_totals_incl.values())
            per_svj["Inclusive"] = (region_totals_incl[SIGNAL_REGION] / denom_incl) if denom_incl > 0 else 0.0
            results[signal_name] = per_svj
    return results


def compute_signal_region_yields(file_path: str, years, signals=None, region: str = SIGNAL_REGION):
    """
    Raw (Run2-combined, already lumi-scaled) yield of each signal landing in
    one ABCD region, per nSVJ bin and inclusive. Unlike
    compute_signal_region_a_acceptance, this does not divide by anything --
    it's the numerator other denominators (generated yield, preselection
    yield, ...) can be compared against.

    Returns: {signal_name: {"0SVJ": yield, ..., "Inclusive": yield}}
    """
    if signals is None:
        signals = [name for name, _ in REFERENCE_SIGNALS]

    results = {}
    with uproot.open(file_path) as f:
        for signal_name in signals:
            per_svj = {}
            total = 0.0
            for svj in SVJ_ORDER:
                y = 0.0
                for year in years:
                    y += read_region_signal(f, svj, year, region, signal_name)
                per_svj[svj] = y
                total += y
            per_svj["Inclusive"] = total
            results[signal_name] = per_svj
    return results


def write_signal_region_a_table(acceptance: dict, outpath: str) -> None:
    """Write a standalone LaTeX table (no \\begin{table} wrapper needed by the
    caller -- this already includes one) of region-A signal acceptance."""
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    signal_labels = dict(REFERENCE_SIGNALS)
    columns = SVJ_ORDER + ["Inclusive"]
    # Matches the notation already used in the background-composition pie
    # chart labels (write_fraction_pie_charts) for the same quantity.
    column_headers = [
        r"$n_{\mathrm{SVJ}}^{\mathrm{PN}} = 0$",
        r"$n_{\mathrm{SVJ}}^{\mathrm{PN}} = 1$",
        r"$n_{\mathrm{SVJ}}^{\mathrm{PN}} = 2$",
        r"$n_{\mathrm{SVJ}}^{\mathrm{PN}} \geq 3$",
        "Inclusive",
    ]

    with open(outpath, "w") as f:
        f.write("\\begin{table}[htbp]\n")
        f.write("\\centering\n")
        f.write("\\begin{tabular}{|l|c|c|c|c|c|}\n")
        f.write("    \\hline\n")
        f.write(
            "    Signal & " + " & ".join(column_headers) + " \\\\ \n"
        )
        f.write("    \\hline\n")
        for signal_name, _ in REFERENCE_SIGNALS:
            label = signal_labels.get(signal_name, signal_name)
            row = acceptance.get(signal_name, {})
            cells = [f"{row.get(c, 0.0) * 100.0:.1f}\\%" for c in columns]
            f.write(f"    {label} & " + " & ".join(cells) + " \\\\ \n")
        f.write("    \\hline\n")
        f.write("\\end{tabular}\n")
        f.write(
            "\\caption{Fraction of Run 2 signal events landing in the ABCD signal "
            "region A ($N_A / (N_A+N_B+N_C+N_D)$), inclusively and in bins of "
            "$n_{\\mathrm{SVJ}}^{\\mathrm{PN}}$, for the reference $m_{\\Phi}$/$r_{\\mathrm{inv}}$ "
            "signal grid used for the AUC scans elsewhere in this document.}\n"
        )
        f.write("\\label{table:signal_region_a_acceptance}\n")
        f.write("\\end{table}\n")
    print(f"[OK] wrote signal region-A acceptance table: {outpath}")


def write_signal_region_a_vs_generated_table(region_a_yields: dict, generated_yields: dict, outpath: str) -> None:
    """
    Same layout as write_signal_region_a_table, but every nSVJ column shares
    one denominator per signal: the total generated (pre-selection) Run2
    yield from CutFlow/Initial, read the same way NMinusOne_maker_RA2.py
    normalizes its Lund-plane correction. Because nSVJ categorization happens
    downstream of generation, these columns are directly additive: they sum
    to the Inclusive column, unlike the region-total-normalized version above.
    """
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    signal_labels = dict(REFERENCE_SIGNALS)
    columns = SVJ_ORDER + ["Inclusive"]
    column_headers = [
        r"$n_{\mathrm{SVJ}}^{\mathrm{PN}} = 0$",
        r"$n_{\mathrm{SVJ}}^{\mathrm{PN}} = 1$",
        r"$n_{\mathrm{SVJ}}^{\mathrm{PN}} = 2$",
        r"$n_{\mathrm{SVJ}}^{\mathrm{PN}} \geq 3$",
        "Inclusive",
    ]

    with open(outpath, "w") as f:
        f.write("\\begin{table}[htbp]\n")
        f.write("\\centering\n")
        f.write("\\begin{tabular}{|l|c|c|c|c|c|}\n")
        f.write("    \\hline\n")
        f.write("    Signal & " + " & ".join(column_headers) + " \\\\ \n")
        f.write("    \\hline\n")
        for signal_name, _ in REFERENCE_SIGNALS:
            label = signal_labels.get(signal_name, signal_name)
            numer = region_a_yields.get(signal_name, {})
            denom = generated_yields.get(signal_name, 0.0)
            if denom > 0:
                cells = [f"{numer.get(c, 0.0) / denom * 100.0:.2f}\\%" for c in columns]
            else:
                # No skim directory on EOS for this sample (e.g. the two
                # DISABLED_SIGNAL_TOKENS entries in NMinusOne_maker_RA2.py) --
                # a generated normalization isn't available, so say so rather
                # than print a misleading 0.00%.
                cells = ["N/A" for _ in columns]
            f.write(f"    {label} & " + " & ".join(cells) + " \\\\ \n")
        f.write("    \\hline\n")
        f.write("\\end{tabular}\n")
        f.write(
            "\\caption{Fraction of \\emph{generated} Run 2 signal events landing "
            "in the ABCD signal region A ($N_A / N_{\\mathrm{generated}}$), "
            "inclusively and in bins of $n_{\\mathrm{SVJ}}^{\\mathrm{PN}}$. "
            "$N_{\\mathrm{generated}}$ is the Run2-combined, luminosity-scaled "
            "yield before any selection, read from each signal sample's "
            "CutFlow/Initial entry -- the same normalization "
            "NMinusOne\\_maker\\_RA2.py uses for its Lund-plane correction. "
            "Unlike the region-A-only table above, these nSVJ columns share a "
            "common denominator and sum to the Inclusive column.}\n"
        )
        f.write("\\label{table:signal_region_a_vs_generated}\n")
        f.write("\\end{table}\n")
    print(f"[OK] wrote signal region-A vs. generated table: {outpath}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", default='/uscms/home/nparmar/nobackup/SVJ/stat_inference_unblind_gapJetveto/stat_histograms/MET_mMed-fullScan_run2_pNet_V5_DNN85_WP90_MET250_BD_lWN_allsyst_min2p_distortion_alternative_trigger_unblind_data_gapJetveto_0to3PSVJ.root', help="Input ROOT file (combine hist file)")
    ap.add_argument("--outdir", default="/uscms/home/ashrivas/nobackup/Dark_Sector/t-channel_Analysis/Yield_plots/wGapveto_and_Qcdtrim", help="Output directory")
    ap.add_argument("--include-flow", action="store_true", help="Include under/overflow in integrals")
    ap.add_argument("--with-data", action="store_true", help="Overlay data_obs if present (OFF by default)")
    ap.add_argument("--years", nargs="*", default=None, help="Years to plot (e.g. 2016 2017 2018). Default: auto-detect")
    ap.add_argument("--lumi", default=None, help="Override lumi text (e.g. '41.5 fb$^{-1}$')")
    ap.add_argument("--com", default="13 ", help="Center-of-mass energy label (default: 13 TeV)")
    ap.add_argument("--final", action="store_true", help="Use 'CMS' instead of 'CMS Preliminary'")
    ap.add_argument("--skip-dnn-score", action="store_true", help="Do not write the DNN score stack plots")
    ap.add_argument("--pie-no-cms-heading", action="store_true", help="Do not draw the CMS/lumi heading on background-fraction pie charts")
    ap.add_argument("--pie-text-fontsize", type=int, default=17, help="Font size for pie-chart percentage labels")
    args = ap.parse_args()


    # with uproot.open(args.file) as f:
    #     year = years[0]

    with uproot.open(args.file) as f:
        years = args.years if args.years else detect_years(f)
    if not years:
        raise RuntimeError("No years detected. Check directory naming (e.g. 0SVJY2016_Run2).")


    # Initialize a master dictionary to hold the sums for ALL years combined
    total_yields = {
        proc: {
            region: {svj: 0.0 for svj in SVJ_ORDER} for region in PLOT_REGIONS
        } for proc in BKG_ORDER
    }
    for y in years:
        with uproot.open(args.file) as f:
            # --- in main(), inside: for y in years:  and inside: with uproot.open(args.file) as f:
            yields = {}

            for proc in BKG_ORDER:
                yields[proc] = {}

                for region_plot in PLOT_REGIONS:
                    region_root = ROOT_REGION_FOR[region_plot]   # <-- MUST be inside this loop
                    yields[proc][region_plot] = {}

                    for svj in SVJ_ORDER:
                        hist_path = f"{svj}Y{y}_Run2/{region_root}/{proc}"

                        if hist_path not in f:
                            yields[proc][region_plot][svj] = 0.0
                            continue

                        h = f[hist_path]
                        values, _ = h.to_numpy()
                        yields[proc][region_plot][svj] = float(values.sum())
                        # 1. First, calculate the value and assign it to 'val'
                        val = float(values.sum())
                        
                        # 2. Save it for the current year's tables/plots
                        yields[proc][region_plot][svj] = val
                        
                        # 3. Add it to the running total for the combined LaTeX tables
                        total_yields[proc][region_plot][svj] += val

        print_yield_summary(y, yields)
        write_abcd_compact_table_txt(y, yields, args.outdir)
        write_detailed_yield_table_txt(y, yields, args.outdir)  
        

        plot_year(
            args.file, y, outdir=args.outdir,
            include_flow=args.include_flow,
            with_data=args.with_data,
            lumi_text=args.lumi,
            com_text=args.com,
            prelim=(not args.final)
        )
        if not args.skip_dnn_score:
            plot_dnn_score(
                args.file, y, outdir=args.outdir,
                with_data=False,
                lumi_text=args.lumi,
                com_text=args.com,
                prelim=(not args.final)
            )
            plot_dnn_score(
                args.file, y, outdir=args.outdir,
                with_data=True,
                lumi_text=args.lumi,
                com_text=args.com,
                prelim=(not args.final)
            )
    plot_year(
        args.file, "Run2_Combined", outdir=args.outdir,
        include_flow=args.include_flow,
        with_data=args.with_data,
        lumi_text=args.lumi,
        com_text=args.com,
        prelim=(not args.final),
        years_to_sum=years,
        output_label="Run2_Combined"
    )
    if not args.skip_dnn_score:
        plot_dnn_score(
            args.file, "Run2_Combined", outdir=args.outdir,
            with_data=False,
            lumi_text=args.lumi,
            com_text=args.com,
            prelim=(not args.final),
            years_to_sum=years,
            output_label="Run2_Combined"
        )
        plot_dnn_score(
            args.file, "Run2_Combined", outdir=args.outdir,
            with_data=True,
            lumi_text=args.lumi,
            com_text=args.com,
            prelim=(not args.final),
            years_to_sum=years,
            output_label="Run2_Combined"
        )
    write_fraction_pie_charts(
        "Run2_Combined", total_yields, args.outdir,
        lumi_text=args.lumi,
        com_text=args.com,
        prelim=(not args.final),
        draw_cms_heading=(not args.pie_no_cms_heading),
        pie_text_fontsize=args.pie_text_fontsize,
    )
    write_latex_fraction_tables("Run2_Combined", total_yields, args.outdir)


if __name__ == "__main__":
    main()
