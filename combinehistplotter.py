#!/usr/bin/env python3
import os
import re
import math
import argparse
import numpy as np
import uproot

import matplotlib.pyplot as plt
import mplhep as hep

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
BKG_ORDER = ["QCD", "TTJets", "WJetsToLNu", "ZJetsToNuNu", "ST"]

PROC_LABEL = {
    "QCD": "QCD multijet",
    "TTJets": r"$t\bar{t}$ + jets",
    "ZJetsToNuNu": r"$Z \to \nu\nu$ + jets",
    "WJetsToLNu": r"$W \to \ell\nu$ + jets",
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
    "WJetsToLNu": "#5790fc", # stays middle
    "ZJetsToNuNu": "#9c9ca1",# was TT
    "ST": "#7a21dd",         # was QCD
}

# Optional: hatches (CMS-style often uses hatches for readability in B/W)
PROC_HATCH = {
    "QCD": None,
    "TTJets": None,
    "WJetsToLNu": None,
    "ZJetsToNuNu": None,
    "ST": None,
}


# Signals to overlay (must match ROOT histogram names exactly)
SIGNAL_PROCS = [
    "mMed500_rinv0p3",
    #"mMed700_rinv0p3",
    "mMed1000_rinv0p3",
    "mMed1500_rinv0p3",
    "mMed2000_rinv0p3",
    "mMed4000_rinv0p3",
]

SIG_COLORS = {
    "mMed500_rinv0p3": "orange",
    #"mMed700_rinv0p3": "green",
    "mMed1000_rinv0p3": "red",
    "mMed1500_rinv0p3": "blue",
    "mMed2000_rinv0p3": "darkgreen",
    "mMed4000_rinv0p3": "purple",
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
            f"{'A':>14}{'B':>14}{'C':>14}{'D':>14}\n"
        )
        f.write("-" * 62 + "\n")

        # Rows
        for svj in SVJ_ORDER:
            f.write(
                f"{svj:<6}"
                f"{fmt(total[svj]['A']):>14}"
                f"{fmt(total[svj]['B']):>14}"
                f"{fmt(total[svj]['C']):>14}"
                f"{fmt(total[svj]['D']):>14}\n"
            )

        # Totals row
        f.write("-" * 62 + "\n")
        f.write(
            f"{'TOTAL':<6}"
            f"{fmt(col_totals['A']):>14}"
            f"{fmt(col_totals['B']):>14}"
            f"{fmt(col_totals['C']):>14}"
            f"{fmt(col_totals['D']):>14}\n"
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
              lumi_text=None, com_text="13", prelim=True):
    os.makedirs(outdir, exist_ok=True)

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
    proc_e = {p: np.zeros_like(x, dtype=float) for p in BKG_ORDER}
    sig_y = {s: np.zeros_like(x, dtype=float) for s in SIGNAL_PROCS}

    data_y = np.full_like(x, np.nan, dtype=float)
    data_e = np.full_like(x, np.nan, dtype=float)


    with uproot.open(file_path) as f:
        for r_i, reg_plot in enumerate(PLOT_REGIONS):
            reg_root = ROOT_REGION_FOR[reg_plot]

            for s_i, svj in enumerate(SVJ_ORDER):
                idx = r_i * nsvj + s_i

                # --- backgrounds ---
                bkg, data = read_region_bkgs(
                    f, svj=svj, year=year, region=reg_root,
                    include_flow=include_flow, with_data=with_data
                )
                if bkg is None:
                    continue

                for p in BKG_ORDER:
                    proc_y[p][idx] = bkg[p][0]
                    proc_e[p][idx] = bkg[p][1]

                # --- signals ---
                for sig in SIGNAL_PROCS:
                    sig_path = f"{svj}Y{year}_Run2/{reg_root}/{sig}"
                    if sig_path in f:
                        hsig = f[sig_path]
                        vals, _ = hsig.to_numpy(flow=include_flow)
                        sig_y[sig][idx] = vals.sum()

                # --- data ---

                if with_data and (data is not None):
                    data_y[idx] = data[0]
                    data_e[idx] = data[1]


    # ---------------------------------------
    # DEBUG: check stacking numerically
    # First bin = Region A, nSVJ = 0
    # ---------------------------------------
    print(
        f"DEBUG {year} (Region A, nSVJ=0):",
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
        f"DEBUG {year} stack check (Region A, nSVJ=0): "
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

        # Ratio: Data/MC
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
        rax.set_ylabel("Data/MC")
    # else:
    #     # Keep the panel for CMS-like layout, but no data points
    #     rax.set_ylabel("Data/MC")

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
        lumi_text = {"2016": "35.9 ", "2017": "41.5 ", "2018": "59.7 "}.get(year, "")
    hep.cms.label(
        ax=ax,
        label="Preliminary" if prelim else "",
        data=with_data,
        lumi=lumi_text,
        com=com_text
    )

    #fig.tight_layout()

    suffix = "_with_data" if with_data else ""
    outpath = os.path.join(outdir, f"ABCD_yields_{year}_CMS{suffix}.pdf")
    fig.savefig(outpath)
    plt.close(fig)
    print(f"[OK] wrote {outpath}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", required=True, help="Input ROOT file (combine hist file)")
    ap.add_argument("--outdir", default="Yield_plots", help="Output directory")
    ap.add_argument("--include-flow", action="store_true", help="Include under/overflow in integrals")
    ap.add_argument("--with-data", action="store_true", help="Overlay data_obs if present (OFF by default)")
    ap.add_argument("--years", nargs="*", default=None, help="Years to plot (e.g. 2016 2017 2018). Default: auto-detect")
    ap.add_argument("--lumi", default=None, help="Override lumi text (e.g. '41.5 fb$^{-1}$')")
    ap.add_argument("--com", default="13 ", help="Center-of-mass energy label (default: 13 TeV)")
    ap.add_argument("--final", action="store_true", help="Use 'CMS' instead of 'CMS Preliminary'")
    args = ap.parse_args()

    with uproot.open(args.file) as f:
        years = args.years if args.years else detect_years(f)
    if not years:
        raise RuntimeError("No years detected. Check directory naming (e.g. 0SVJY2016_Run2).")

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

        print_yield_summary(y, yields)
        write_abcd_compact_table_txt(y, yields, args.outdir)




        plot_year(
            args.file, y, outdir=args.outdir,
            include_flow=args.include_flow,
            with_data=args.with_data,
            lumi_text=args.lumi,
            com_text=args.com,
            prelim=(not args.final)
        )


if __name__ == "__main__":
    main()