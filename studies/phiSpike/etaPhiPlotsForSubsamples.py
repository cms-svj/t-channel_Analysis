# use coffea environment to do this
# see if some of the hot spots are caused by high weight events
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
import utility as utl 

mpl.use('Agg')

mpl.rc("font", family="serif", size=20)

plotFormat = "png"

cut = "_pre" ##
etaMin = -2.4
etaMax = 2.4
phiMin = -3.2
phiMax = 3.2
normalization = False ##
nSigma = 5
nSigmaLab = str(np.round(nSigma,2)).replace(".","p")
plotFolderPre = f"plots/{cut[1:]}/sigma{nSigmaLab}"
inputFileDict = {
                    # "2016": "/srv/TM_phiSpike_opt/",
                    "2017": "/srv/TM_phiSpike_opt/",
                    "2018": "/srv/TM_phiSpike_opt/",
}

sampleDict = {
    "QCD": [
                "Q80",
                "Q120",
                "Q170",
                "Q300",
                "Q470",
                "Q600",
                "Q800",
                "Q1000",
                "Q1400",
                "Q1800",
                "Q2400",
                "Q3200",
            ]
}

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
        for sample, subList in sampleDict.items():
            plotFolder = f"{plotFolderPre}/{year}/{sample}/subsample/"
            for sub in subList:
                inputFile = f"{year}_{sub}_sub.root"
                os.makedirs(plotFolder,exist_ok=True)
                varName = "h_{}".format(var)
                zVals, xVals, yVals = utl.getTotal2DHist(inputFolder,inputFile,varName,cut)
                zVals = np.rot90(zVals,k=1)
                zValsUnnormed = zVals
                desiredYLabels = [-3,-2,-1,0,1,2,3]
                desiredYTicks, linFitY = utl.getTicks(np.arange(zVals.shape[1] - 0.5,-1,-1), yVals, desiredYLabels)
                desiredXLabels = [-2,-1,0,1,2]
                desiredXTicks, linFitX = utl.getTicks(np.arange(-.5,zVals.shape[0]), xVals, desiredXLabels)
                xmin = utl.actualToBinVal(etaMin, linFitX)
                xmax = utl.actualToBinVal(etaMax, linFitX)
                ymin = utl.actualToBinVal(phiMin, linFitY)
                ymax = utl.actualToBinVal(phiMax, linFitY)
                if normalization:
                    highestEta = np.amax(zVals,axis=0)
                    zVals = zVals/highestEta
                    zVals -= np.average(zVals,axis=0)
                    if "PostHEM" in plotFolder:
                        utl.widenHEMMask(linFitX, linFitY, zVals, extraMask=2)
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
                    plt.savefig(f"{plotFolder}/normed_2d_{sub}_{var}.{plotFormat}", bbox_inches="tight")  
                else:
                    plt.savefig(f"{plotFolder}/2d_{sub}_{var}.{plotFormat}", bbox_inches="tight")