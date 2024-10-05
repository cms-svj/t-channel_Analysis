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

def getHisto(inputFolder,fileName,varName,cut,dim=1):
    fullInputFilePath = "{}{}".format(inputFolder,fileName)
    f = up.open(fullInputFilePath)
    hist = f[varName + cut]
    if dim == 1:
        npHist, bins = hist.to_numpy() 
        histHist = hist.to_hist() #*(1./np.sum(npHist)) # normalization, however, ratioplot will not calculate the error properly
        return npHist, histHist, bins
    elif dim == 2:
        npHist = hist.to_numpy() 
        histHist = hist.to_hist() #*(1./np.sum(npHist)) # normalization, however, ratioplot will not calculate the error properly
        return npHist, histHist

def getTotal2DHist(inputFolder,inputFile,varName,cut):
    npHist, histHist = getHisto(inputFolder,inputFile,varName,cut,dim=2)
    npZ = npHist[0]
    npX = npHist[1]
    npY = npHist[2]
    return npZ,npX,npY

def getTicks(fullTicks, fullTickLabels, desiredLabels):
    linFit = linregress(fullTickLabels, fullTicks)
    m = linFit[0]
    b = linFit[1]
    desiredTicks = []
    for lab in desiredLabels:
        desiredTicks.append(m*lab+b)
    return desiredTicks, linFit

def actualToBinVal(val, linFit):
    m = linFit[0]
    b = linFit[1]
    return m*val+b

def bintoActualVal(binVal, linFit):
    m = linFit[0]
    b = linFit[1]
    return (binVal-b)/m

def widenHEMMask(linFitX, linFitY, zVals, extraMask=2):
    # extraMask = number of extra bins to widen the HEM veto by
    etaVetoMin = int(actualToBinVal(-3.05, linFitX) + 0.5)
    etaVetoMax = int(actualToBinVal(-1.35, linFitX) + 0.5)
    phiVetoMax = int(actualToBinVal(-0.82, linFitY) + 0.5) - extraMask
    phiVetoMin = int(actualToBinVal(-1.62, linFitY) + 0.5) + extraMask
    zVals[phiVetoMax:phiVetoMin,etaVetoMin:etaVetoMax] = 0

def uniqueCombine(a1,a2,b1,b2):
    c1 = list(a1) + list(b1)
    c2 = list(a2) + list(b2)
    s = []
    for i in range(len(c1)):
        s.append(f"{c1[i]}_{c2[i]}")
    ulist, uInds = np.unique(s,return_index=True)
    u1 = []
    u2 = []
    for i in uInds:
        u1.append(c1[i])
        u2.append(c2[i])
    return u1, u2

etaMin = -2.4
etaMax = 2.4
phiMin = -3.2
phiMax = 3.2
normalization = False ##
nSigma = 5
cut = "_pre_psFilterSig5" ##
nSigmaLab = str(np.round(nSigma,2)).replace(".","p")
plotFolderPre = f"plots/{cut[1:]}/sigma{nSigmaLab}"
inputFileDict = {
                    "2016": "/srv/TM_phiSpike_opt_clean/",
                    "2017": "/srv/TM_phiSpike_opt_clean/",
                    # "2018": "/srv/TM_phiSpike_opt_clean/",
}
samples = ["QCD"]
variables2D = {
                'j1Phivsj1Eta': [r"$\eta(j_1)$",r"$\phi(j_1)$"],
                'j2Phivsj2Eta': [r"$\eta(j_2)$",r"$\phi(j_2)$"],
                'j3Phivsj3Eta': [r"$\eta(j_3)$",r"$\phi(j_3)$"],
                'j4Phivsj4Eta': [r"$\eta(j_4)$",r"$\phi(j_4)$"],
}

flatNessDict = {}
phiSpikeHotSpots = {}
for year, inputFolder in inputFileDict.items():
    hsDict = {}
    for var,labels in variables2D.items():
        jiEtaHotSpots = []
        jiPhiHotSpots = []
        for sample in samples:
            plotFolder = f"{plotFolderPre}/{year}/{sample}"
            inputFile = f"{year}_{sample}_all.root"
            os.makedirs(plotFolder,exist_ok=True)
            varName = "h_{}".format(var)
            zVals, xVals, yVals = getTotal2DHist(inputFolder,inputFile,varName,cut)
            zVals = np.rot90(zVals,k=1)
            zValsUnnormed = zVals
            desiredYLabels = [-3,-2,-1,0,1,2,3]
            desiredYTicks, linFitY = getTicks(np.arange(zVals.shape[1] - 0.5,-1,-1), yVals, desiredYLabels)
            desiredXLabels = [-2,-1,0,1,2]
            desiredXTicks, linFitX = getTicks(np.arange(-.5,zVals.shape[0]), xVals, desiredXLabels)
            xmin = actualToBinVal(etaMin, linFitX)
            xmax = actualToBinVal(etaMax, linFitX)
            ymin = actualToBinVal(phiMin, linFitY)
            ymax = actualToBinVal(phiMax, linFitY)
            if normalization:
                highestEta = np.amax(zVals,axis=0)
                zVals = zVals/highestEta
                zVals -= np.average(zVals,axis=0)
                if "PostHEM" in plotFolder:
                    widenHEMMask(linFitX, linFitY, zVals, extraMask=2)
            # 2D plots heat map that shows the hot spots
            fig = plt.figure(figsize=(10, 10))
            plt.imshow(zVals)
            plt.yticks(desiredYTicks,desiredYLabels)
            plt.xticks(desiredXTicks,np.flip(desiredXLabels))
            plt.xlabel(labels[0])
            plt.ylabel(labels[1])
            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)
            plt.colorbar()
            if normalization:
                plt.savefig("{}/normed_2d_{}.{}".format(plotFolder,var,plotFormat), bbox_inches="tight")  
            else:
                plt.savefig("{}/2d_{}.{}".format(plotFolder,var,plotFormat), bbox_inches="tight")