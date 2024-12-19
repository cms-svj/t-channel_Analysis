# use coffea environment to do this
import uproot as up 
import mplhep as hep
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import coffea
import os
from hist import Hist
from scipy.stats import linregress
import pickle

mpl.use('Agg')

mpl.rc("font", family="serif", size=20)

plotFormat = "png"

def getYield(inputFile,cut):
    f = up.open(inputFile)
    hist = f["h_evtw" + cut]
    npHist, bins = hist.to_numpy() 
    coffeaHist = hist.to_hist() #*(1./np.sum(npHist)) # useful for ratioPlot
    return np.sum(npHist)

def sEff(si,sf):
    return sf/si

def fom(B,S):
    return np.sqrt( 2*( (S+B) * np.log(1+S/B) - S ) )

def plotFOMBySigPara(metrics,sigParaName):
    plt.figure()
    sigmaVals = metrics["sigmaVals"]
    bkgYields = metrics["QCDYields"]
    sigYieldsDict = metrics["sigYieldsDict"]
    preQCDYield = metrics["preQCDYield"]
    preSigYieldDict = metrics["preSigYieldDict"]
    for sigLab in sigYieldsDict:
        if sigParaName not in sigLab:
            continue
        sigYields = sigYieldsDict[sigLab]
        fomVals = fom(bkgYields,sigYields)
        plt.plot(sigmaVals,fomVals/np.max(fomVals),marker="o",label=f"{sigLab}",color=sigDict[sigLab][1])
        plt.axhline(fom(preQCDYield,preSigYieldDict[sigLab])/np.max(fomVals),linestyle="-",color=sigDict[sigLab][1])
    plt.xlabel("Number of Sigmas")
    plt.ylabel("Normalized FOM")
    plt.legend(loc=(1.04,0))
    plt.savefig(f"FOM_{year}_{sigParaName}.pdf",bbox_inches="tight") 

cutDict = {
"_pre_psFilterSig2":    2,
"_pre_psFilterSig2p5":  2.5,
"_pre_psFilterSig3":    3,
"_pre_psFilterSig3p5":  3.5,
"_pre_psFilterSig4":    4,
"_pre_psFilterSig4p5":  4.5,
"_pre_psFilterSig5":    5,
"_pre_psFilterSig5p5":  5.5,
"_pre_psFilterSig6":    6,
"_pre_psFilterSig6p5":  6.5,
"_pre_psFilterSig7":    7,
"_pre_psFilterSig7p5":  7.5,
# "_pre_psFilterSig8":    8,
# "_pre_psFilterSig8p5":  8.5,
}

sigDict = {
"mMed-600":   ["m600_d20_r0p3_y1_N-1_M0_.root",     "#ff7f0e"],
"mMed-800":   ["m800_d20_r0p3_y1_N-1_M0_.root",     "#2ca02c"],
"mMed-1000":  ["m1000_d20_r0p3_y1_N-1_M0_.root",    "#bcbd22"],
"mMed-1500":  ["m1500_d20_r0p3_y1_N-1_M0_.root",    "#17becf"],
"mMed-2000":  ["m2000_d20_r0p3_y1_N-1_M0_.root",    "#7f7f7f"],
"mMed-3000":  ["m3000_d20_r0p3_y1_N-1_M0_.root",    "#8c564b"],
"mMed-4000":  ["m4000_d20_r0p3_y1_N-1_M0_.root",    "#1f77b4"],
"mDark-1":    ["m2000_d1_r0p3_y1_N-1_M0_.root",     "#d62728"],
"mDark-20":   ["m2000_d20_r0p3_y1_N-1_M0_.root",    "#7f7f7f"],
"mDark-50":   ["m2000_d50_r0p3_y1_N-1_M0_.root",    "#9467bd"],
"mDark-100":  ["m2000_d100_r0p3_y1_N-1_M0_.root",   "#e377c2"],
"rinv-0.1":   ["m2000_d20_r0p1_y1_N-1_M0_.root",    "#000000"],
"rinv-0.3":   ["m2000_d20_r0p3_y1_N-1_M0_.root",    "#7f7f7f"],
"rinv-0.5":   ["m2000_d20_r0p5_y1_N-1_M0_.root",    "#ff7f0e"],
"rinv-0.7":   ["m2000_d20_r0p7_y1_N-1_M0_.root",    "#bcbd22"],
}

preYieldFolder = "/srv/skim_phiSpike_restudy_fixedTrgBug_clean"
defaultCut = "_pre_psFilterSig4" # this is the cut before any additional phi spike cut is applied
inputFolder = "/srv/skim_phiSpike_restudy_fixedTrgBug_pre_psFilterSig4_clean"
# years = [2016,2017,2018]
years = [
"2016",
"2017",
"2018PreHEM",
"2018PostHEM"
]

metricDict = {}
print("Extracting useful quantities for calculating metrics...")
for year in years:
    print(f"{year}")
    # FOM vs nSigmas
    bkgFile = f"{year}_QCD_all.root"
    preQCDYield = getYield(f"{preYieldFolder}/{bkgFile}",defaultCut)
    sigmaVals = []
    bkgYields = []
    for cutLab, cutVal in cutDict.items():
        bkgYield = getYield(f"{inputFolder}/{bkgFile}",cutLab)
        sigmaVals.append(cutVal)
        bkgYields.append(bkgYield)
    preSigYieldDict = {}
    sigYieldsDict = {}
    for sigLab, sigDetails in sigDict.items():
        sigFilePost = sigDetails[0]
        print(f"Processing {sigLab}")
        sigFile = f"{year}_{sigFilePost}"
        preSigYield = getYield(f"{preYieldFolder}/{sigFile}",defaultCut)
        preSigYieldDict[sigLab] = preSigYield
        sigYields = []
        for cutLab, cutVal in cutDict.items():
            sigYield = getYield(f"{inputFolder}/{sigFile}",cutLab)
            sigYields.append(sigYield)
        sigYieldsDict[sigLab] = np.array(sigYields)
    metricDict[year] = {"sigmaVals":np.array(sigmaVals),"preQCDYield":preQCDYield,"QCDYields":np.array(bkgYields),"preSigYieldDict":preSigYieldDict,"sigYieldsDict":sigYieldsDict}

for year, metrics in metricDict.items():
    print(year)
    plotFOMBySigPara(metrics,"mMed")
    plotFOMBySigPara(metrics,"mDark")
    plotFOMBySigPara(metrics,"rinv")

