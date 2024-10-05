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

def fom(B,S):
    return np.sqrt( 2*( (S+B) * np.log(1+S/B) - S ) )

cutDict = {
"_pre_psFilterSig4":    4,
"_pre_psFilterSig4p5":  4.5,
"_pre_psFilterSig5":    5,
"_pre_psFilterSig5p5":  5.5,
"_pre_psFilterSig6":    6,
}

inputFolder = "/srv/TM_phiSpike_opt_clean"
years = [2016,2017]

for year in years:
    bkgFile = f"{year}_Q_all.root"
    sigFile = f"{year}_m2000_d20_r0p3_y1_N-1_M0_.root"
    fomVals = []
    sigmaVals = []
    for cutLab, cutVal in cutDict.items():
        sigYield = getYield(f"{inputFolder}/{sigFile}",cutLab)
        bkgYield = getYield(f"{inputFolder}/{bkgFile}",cutLab)
        fomVals.append(fom(bkgYield,sigYield))
        sigmaVals.append(cutVal)
    plt.plot(sigmaVals,fomVals,label=f"{year}")
    plt.xlabel("Number of sigmas")
    plt.ylabel("FOM (QCD vs baseline signal)")
    plt.legend()
    plt.savefig("test.pdf",bbox_inches="tight")