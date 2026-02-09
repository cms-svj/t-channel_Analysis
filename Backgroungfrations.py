import ROOT
import ROOTplotutils as pltutils
import optparse
import utils.DataSetInfo as info
import utils as TF
#import plotStack
import utils.CMS_lumi as CMS_lumi
import numpy as np
import os
from collections import defaultdict

# ============================================================
# GLOBAL accumulator for background yields (ALL YEARS)
# ============================================================

GLOBAL_BKG_YIELDS = defaultdict(lambda: {
    "0SVJ": 0.0,
    "1SVJ": 0.0,
    "2SVJ": 0.0,
    "2PSVJ": 0.0,
    "3PSVJ": 0.0,
})

# ============================================================
# SVJ bins
# ============================================================

def GetSVJbins(DNN_inner_edge, MET_inner_edge):
    return {
        "0SVJ": [DNN_inner_edge, MET_inner_edge],
        "1SVJ": [DNN_inner_edge, MET_inner_edge],
        "2SVJ": [DNN_inner_edge, MET_inner_edge],
        "2PSVJ": [DNN_inner_edge, MET_inner_edge],
        "3PSVJ": [DNN_inner_edge, MET_inner_edge],
    }

# ============================================================
# Dataset loading
# ============================================================

def getData(path, scale=1.0, year = "2018"):
    '''Uncomment the files that are to be included in the plots'''
    Data = [
        info.DataSetInfo(basedir=path, fileName=year+"_Data.root",        sys= -1.0, label="Data",        scale=scale),
    ]
    bgData = [
        info.DataSetInfo(basedir=path, fileName=year+"_ST.root",              label="Single top",              scale=scale, color=ROOT.TColor.GetColor("#5790fc")),
        info.DataSetInfo(basedir=path, fileName=year+"_TTJets.root",          label="t#bar{t}",                scale=scale, color=ROOT.TColor.GetColor("#f89c20")),
        info.DataSetInfo(basedir=path, fileName=year+"_ZJets.root",             label="Z#rightarrow#nu#nu+jets",    scale=scale, color=ROOT.TColor.GetColor("#e42536")),
        #info.DataSetInfo(basedir=path, fileName=year+"_WJets.root",              label="W+jets",                    scale=scale, color=ROOT.TColor.GetColor("#964a8b")),
        info.DataSetInfo(basedir=path, fileName=year+"_QCD.root",               label="QCD",                        scale=scale, color=ROOT.TColor.GetColor("#9c9ca1")),
        #info.DataSetInfo(basedir=path, fileName="2018_allBkg.root",               label="all bckg",                        scale=scale, color=ROOT.TColor.GetColor("#9c9ca1")),
    ]
    # Only include WJets if not 2016
    if year != "2016":
        bgData.insert(
            3,
            info.DataSetInfo(basedir=path, fileName=year+"_WJets.root", label="W+jets",
                             scale=scale, color=ROOT.TColor.GetColor("#964a8b"))
        )

    
    sgData = [
       
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="baseline", scale=scale, color=ROOT.kBlack),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",     label="mMed_600",  scale=scale, color=ROOT.kViolet+2),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-800_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",     label="mMed_800",  scale=scale, color=ROOT.kRed),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-1000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="mMed_1000", scale=scale, color=ROOT.kMagenta + 1),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-1500_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="mMed_1500", scale=scale, color=ROOT.kGray+4),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="mMed_3000", scale=scale, color=ROOT.kCyan),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="mMed_4000", scale=scale, color=ROOT.kBlue + 1),
        info.DataSetInfo(basedir=path, fileName=year+"_m2000_d20_r0p3_y1_N-1_M0_.root",    label="baseline", scale=scale, color=ROOT.kBlue + 1),
        # info.DataSetInfo(basedir=path, fileName=year+"_m600_d20_r0p3_y1_N-1_M0_.root",     label="mMed_600",  scale=scale, color=ROOT.kViolet+2),
        # info.DataSetInfo(basedir=path, fileName=year+"_m800_d20_r0p3_y1_N-1_M0_.root",     label="mMed_800",  scale=scale, color=ROOT.kRed),
        # info.DataSetInfo(basedir=path, fileName=year+"_m1000_d20_r0p3_y1_N-1_M0_.root",    label="mMed_1000", scale=scale, color=ROOT.kMagenta + 1),
        # info.DataSetInfo(basedir=path, fileName=year+"_m1500_d20_r0p3_y1_N-1_M0_.root",    label="mMed_1500", scale=scale, color=ROOT.kGray+4),
        # #info.DataSetInfo(basedir=path, fileName=year+"_m3000_d20_r0p3_y1_N-1_M0_.root",    label="mMed_3000", scale=scale, color=ROOT.kCyan),
        # info.DataSetInfo(basedir=path, fileName=year+"_m4000_d20_r0p3_y1_N-1_M0_.root",    label="mMed_4000", scale=scale, color=ROOT.kBlue + 1),
    ]
    return Data, sgData, bgData

# ============================================================
# ABCD regions
# ============================================================

def GetSubABCDregions(dnn_i, dnn_o, met_i, met_o):
    return [
        ("A", met_i, met_o, dnn_i, dnn_o),
        ("B", met_i, met_o, 0.0, dnn_i),
        ("C", 0.0, met_i, dnn_i, dnn_o),
        ("D", 0.0, met_i, 0.0, dnn_i),
    ]

# ============================================================
# Extract ABCD counts PER PROCESS
# ============================================================

def get_ABCD_counts_by_SVJ(
    Data, sgData, bgData,
    ABCDHistoVar, maincut,
    SVJBins,
    DNN_inner_edge, DNN_outer_edge,
    MET_inner_edge, MET_outer_edge,
):
    counts_by_SVJ = {}

    for data in bgData:
        process = data.label_

        for svj, _ in SVJBins.items():
            histName = ABCDHistoVar + maincut + svj

            NA = NB = NC = ND = 0.0

            for region, xmin, xmax, ymin, ymax in GetSubABCDregions(
                DNN_inner_edge, DNN_outer_edge,
                MET_inner_edge, MET_outer_edge
            ):
                _, integral, _ = data.get2DHistoIntegral(
                    histName,
                    xmin=xmin, xmax=xmax,
                    ymin=ymin, ymax=ymax,
                    showEvents=False
                )

                if region == "A": NA += integral
                if region == "B": NB += integral
                if region == "C": NC += integral
                if region == "D": ND += integral

            counts_by_SVJ[(svj, process)] = (NA, NB, NC, ND)

    return counts_by_SVJ

# ============================================================
# ABCD driver (UNCHANGED except accumulation)
# ============================================================

def compute_ABCD_prediction(
    Data, sgData, bgData,
    ABCDHistoVar, maincut, Year,
    DNN_inner_edges, DNN_outer_edges,
    MET_inner_edges, MET_outer_edges,
    VR
):
    outer_edge_results = {}

    for dnn_i, dnn_o, met_i, met_o in zip(
        DNN_inner_edges, DNN_outer_edges,
        MET_inner_edges, MET_outer_edges
    ):
        SVJBins = GetSVJbins(dnn_i, met_i)

        counts = get_ABCD_counts_by_SVJ(
            Data, sgData, bgData,
            ABCDHistoVar, maincut,
            SVJBins,
            dnn_i, dnn_o, met_i, met_o
        )

        outer_edge_results[(dnn_o, met_o)] = counts

    return outer_edge_results

# ============================================================
# WRITE FINAL TABLE (ALL YEARS)
# ============================================================

def write_background_fractions_txt(global_yields, output_txt):
    """
    Writes background fractions in a paper-style table.
    Percentages are formatted and aligned.
    """

    svj_bins = ["0SVJ", "1SVJ", "2SVJ", "2PSVJ", "3PSVJ"]

    # Paper-style ordering and labels
    process_order = [
        ("QCD", "QCD multijet"),
        ("t#bar{t}", "t̄t + jets"),
        ("Z#rightarrow#nu#nu+jets", "Z → νν + jets"),
        ("W+jets", "W + jets"),
        ("Single top", "Single top"),
    ]

    with open(output_txt, "w") as f:
        f.write("Table X: Background composition [%] using A+B+C+D (2016–2018)\n\n")
        f.write("-" * 63 + "\n")
        f.write("{:<22s} {:>7s} {:>7s} {:>7s} {:>8s} {:>8s}\n".format(
            "Process", "0 SVJ", "1 SVJ", "2 SVJ", "2P SVJ", "3P SVJ"
        ))
        f.write("-" * 63 + "\n")

        for key, label in process_order:
            if key not in global_yields:
                continue

            row = [label]
            for svj in svj_bins:
                total = sum(global_yields[p][svj] for p in global_yields)
                frac = 100.0 * global_yields[key][svj] / total if total > 0 else 0.0
                row.append(f"{frac:5.3f}")

            f.write("{:<22s} {:>7s} {:>7s} {:>7s} {:>8s} {:>8s}\n".format(*row))

        f.write("-" * 63 + "\n")

    print(f"[OK] Wrote formatted background table → {output_txt}")


def write_svj_background_efficiency_txt(global_yields, output_txt):
    """
    Writes SVJ background efficiency table (Table 18 style).
    """

    svj_bins = ["0SVJ", "1SVJ", "2SVJ", "3PSVJ"]

    # Total background across all SVJ bins
    total_bkg = sum(
        sum(global_yields[p][svj] for p in global_yields)
        for svj in svj_bins
    )

    with open(output_txt, "w") as f:
        f.write("Table Y: Background efficiency per SVJ category [%] (2016–2018)\n\n")
        f.write("-" * 45 + "\n")
        f.write("{:<12s} {:>15s}\n".format("SVJ category", "Efficiency [%]"))
        f.write("-" * 45 + "\n")

        for svj in svj_bins:
            svj_sum = sum(global_yields[p][svj] for p in global_yields)
            eff = 100.0 * svj_sum / total_bkg if total_bkg > 0 else 0.0
            f.write("{:<12s} {:>15.3f}\n".format(svj.replace("SVJ", " SVJ"), eff))

        f.write("-" * 45 + "\n")

    print(f"[OK] Wrote SVJ efficiency table → {output_txt}")

# ============================================================
# MAIN
# ============================================================

def main():
    parser = optparse.OptionParser()
    parser.add_option("-d", "--dataset", dest="dataset", default="testHadd")
    parser.add_option("-o", "--outputdir", dest="outputdir", default="ABCD_Output")
    options, _ = parser.parse_args()

    Years = ["2016", "2017", "2018"]

    for Year in Years:
        Data, sgData, bgData = getData(options.dataset + "/", 1.0, Year)

        DNN_outer_edges = np.linspace(1.0, 1.0, 1)
        DNN_inner_edges = np.linspace(0.6, 0.6, 1)
        MET_outer_edges = np.linspace(250, 250, 1)
        MET_inner_edges = np.linspace(210, 210, 1)

        outer_edge_results = compute_ABCD_prediction(
            Data, sgData, bgData,
            "h_METvsDNN", "_pre_", Year,
            DNN_inner_edges, DNN_outer_edges,
            MET_inner_edges, MET_outer_edges,
            "VR3"
        )

        # Accumulate background yields
        counts = list(outer_edge_results.values())[0]
        for (svj, process), (NA, NB, NC, ND) in counts.items():
            GLOBAL_BKG_YIELDS[process][svj] += (NA + NB + NC + ND)

    os.makedirs(options.outputdir, exist_ok=True)
    write_background_fractions_txt(
        GLOBAL_BKG_YIELDS,
        os.path.join(options.outputdir, "BackgroundFractions_AllYears.txt")
    )
    write_svj_background_efficiency_txt(
    GLOBAL_BKG_YIELDS,
    os.path.join(options.outputdir, "SVJ_BackgroundEfficiency_AllYears.txt")
    )



if __name__ == "__main__":
    main()
