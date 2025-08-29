import os
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from collections import defaultdict
import optparse
import ROOT

import ROOT
import ROOTplotutils as pltutils
import optparse
import pandas as pd
import utils.DataSetInfo as info
import MeTutils as TF
import plotStack
import utils.CMS_lumi as CMS_lumi
from MetAnalysis import GetSVJbins 
import numpy as np
import os 
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator

plt.style.use(hep.style.CMS)

process_colors = {
    "TTJets": "#f89c20",
    "WJets": "#964a8b",
    "ZJets": "#e42536",
    "QCD": "#9c9ca1",
    "ST": "#5790fc",
}

signal_style = {'linestyle': '--', 'linewidth': 2, 'color': 'black', 'label': 'Signal MC'}
data_marker = {'marker': 'o', 'linestyle': 'None', 'color': 'black', 'label': 'Data'}

def getData(path, scale=1.0, year=2016):
    '''Uncomment the files that are to be included in the plots'''
    Data = [
        info.DataSetInfo(basedir=path, fileName=year+"_Data.root",        sys= -1.0, label="Data",        scale=scale),
    ]
    bgData = [
        info.DataSetInfo(basedir=path, fileName=year+"_ST.root",              label="Single top",              scale=scale, color=ROOT.TColor.GetColor("#5790fc")),
        info.DataSetInfo(basedir=path, fileName=year+"_TTJets.root",          label="t#bar{t}",                scale=scale, color=ROOT.TColor.GetColor("#f89c20")),
        info.DataSetInfo(basedir=path, fileName=year+"_ZJets.root",             label="Z#rightarrow#nu#nu+jets",    scale=scale, color=ROOT.TColor.GetColor("#e42536")),
        info.DataSetInfo(basedir=path, fileName=year+"_WJets.root",              label="W+jets",                    scale=scale, color=ROOT.TColor.GetColor("#964a8b")),
        info.DataSetInfo(basedir=path, fileName=year+"_QCD.root",               label="QCD",                        scale=scale, color=ROOT.TColor.GetColor("#9c9ca1")),
    ]
    sgData = [
        # TODO: check if there is a cms guidline for line color
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="baseline", scale=scale, color=ROOT.kBlack),
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",     label="mMed_600",  scale=scale, color=ROOT.kViolet+2),
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-800_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",     label="mMed_800",  scale=scale, color=ROOT.kRed),
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-1000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="mMed_1000", scale=scale, color=ROOT.kMagenta + 1),
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-1500_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="mMed_1500", scale=scale, color=ROOT.kGray+4),
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="mMed_3000", scale=scale, color=ROOT.kCyan),
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="mMed_4000", scale=scale, color=ROOT.kBlue + 1),
        
    ]
    return Data, sgData, bgData
def plot_MET_stacked(histoname, Data, sgData, bgData, outdir, year):
    bin_edges = None
    bkg_contrib = defaultdict(list)
    signal_vals = []
    data_vals = []
    data_errs = []

    # Get background histograms
    for proc in bgData:
        hist = proc.get1DHisto(histoname)
        if bin_edges is None:
            bin_edges = [hist.GetBinLowEdge(i+1) for i in range(hist.GetNbinsX())] + [hist.GetBinLowEdge(hist.GetNbinsX()) + hist.GetBinWidth(hist.GetNbinsX())]
        bkg_contrib[proc.label] = [hist.GetBinContent(i+1) for i in range(hist.GetNbinsX())]

    # Get signal histogram
    for sig in sgData:
        hist = sig.get1DHisto(histoname)
        signal_vals = [hist.GetBinContent(i+1) for i in range(hist.GetNbinsX())]

    # Get data histogram
    for dat in Data:
        hist = dat.get1DHisto(histoname)
        data_vals = [hist.GetBinContent(i+1) for i in range(hist.GetNbinsX())]
        data_errs = [hist.GetBinError(i+1) for i in range(hist.GetNbinsX())]

    bin_centers = [(bin_edges[i] + bin_edges[i+1])/2 for i in range(len(bin_edges)-1)]
    plt.figure(figsize=(10, 6))
    bottoms = np.zeros(len(bin_centers))

    for label, vals in bkg_contrib.items():
        plt.bar(bin_centers, vals, width=np.diff(bin_edges), label=label, color=process_colors.get(label, None), edgecolor='black', bottom=bottoms)
        bottoms += vals

    plt.plot(bin_centers, signal_vals, **signal_style)
    plt.errorbar(bin_centers, data_vals, yerr=data_errs, **data_marker)

    hep.cms.label(rlabel="", year=year)
    plt.xlabel("MET [GeV]")
    plt.ylabel("Events")
    plt.title(f"Stacked MET Plot: {histoname}")
    plt.legend()
    plt.tight_layout()
    os.makedirs(outdir, exist_ok=True)
    plt.savefig(os.path.join(outdir, f"{histoname}.pdf"))
    plt.close()

def main():
    parser = optparse.OptionParser("usage: %prog [options]\n")
    parser.add_option('-d', '--dataset', dest='dataset', default='output/', help='dataset directory')
    parser.add_option('-y', dest='year', type='string', default='2018', help='Year')
    parser.add_option('-o', dest='outputdir', type='string', help='Output folder name')
    options, args = parser.parse_args()

    year = options.year
    dataset = options.dataset
    output_dir = options.outputdir or "stacked_MET_bySVJ"
    histnames = [
        "h_MET_pre_0SVJ", "h_MET_pre_1SVJ", "h_MET_pre_2J_1PSVJ", "h_MET_pre_2PSVJ",
        "h_MET_pre_2SVJ", "h_MET_pre_3J_1PSVJ", "h_MET_pre_3SVJ", "h_MET_pre_4J_1PSVJ",
        "h_MET_pre_4PSVJ"
    ]

    Data, sgData, bgData = getData(dataset + "/", 1.0, year)

    for histoname in histnames:
        plot_MET_stacked(histoname, Data, sgData, bgData, output_dir, year)

if __name__ == '__main__':
    main()
