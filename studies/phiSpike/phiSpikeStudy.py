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

cut = "_pre" ##
etaMin = -2.4
etaMax = 2.4
phiMin = -3.2
phiMax = 3.2
normalization = True ##
nSigma = 6
nSigmaLab = str(np.round(nSigma,2)).replace(".","p")
plotFolderPre = f"plots/{cut[1:]}/sigma{nSigmaLab}"
inputFileDict = {
                    "2016": "/srv/TM_phiSpike_2016_clean/",
                    "2017": "/srv/TM_phiSpike_2017_clean/",
                    "2018PreHEM": "/srv/TM_phiSpike_2018PreHEM_clean/",
                    "2018PostHEM": "/srv/TM_phiSpike_2018PostHEM_clean/",
}
samples = ["Data","QCD"]
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
            # 1D plots to decide where the n-sigma threshold is
            zValsZoomed = zVals[int(ymax+0.5):int(ymin+0.5),int(xmin+0.5):int(xmax+0.5)] 
            flatzVals = zValsZoomed.flatten()
            flatzVals = flatzVals[np.isfinite(flatzVals)]
            mean = np.std(flatzVals)
            sigma = np.std(flatzVals)
            threshold = sigma*nSigma
            plt.figure(figsize=(10, 10))
            hist,bins = np.histogram(flatzVals,bins=100)
            hep.histplot(hist,bins)
            plt.ylabel("Number of Bins")
            plt.yscale("log")
            if normalization:
                plt.xlabel("Normalized Events")
                plt.axvline(threshold,color="red",linestyle="--")
                plt.savefig("{}/normed_1d_{}.{}".format(plotFolder,var,plotFormat), bbox_inches="tight")
            else:
                plt.xlabel("Events")
                plt.savefig("{}/1d_{}.{}".format(plotFolder,var,plotFormat), bbox_inches="tight")
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
                hotSpots = np.where(zVals>threshold)
                hotEtaVals = hotSpots[1]
                hotPhiVals = hotSpots[0]
                plt.plot(hotEtaVals,hotPhiVals,color="white",marker="o",linestyle="")
                jiEtaHotSpot = bintoActualVal(hotEtaVals, linFitX)
                jiPhiHotSpot = bintoActualVal(hotPhiVals, linFitY)
                jiEtaHotSpots,jiPhiHotSpots = uniqueCombine(jiEtaHotSpots,jiPhiHotSpots,jiEtaHotSpot,jiPhiHotSpot)
                plt.savefig("{}/normed_2d_{}_hotSpots.{}".format(plotFolder,var,plotFormat), bbox_inches="tight")  
            else:
                plt.savefig("{}/2d_{}.{}".format(plotFolder,var,plotFormat), bbox_inches="tight")
            # flatness measure after removing hot spots
            zValsUnnormedZoomed = zValsUnnormed[int(ymax+0.5):int(ymin+0.5),int(xmin+0.5):int(xmax+0.5)] 
            zValsClean = np.where(zValsZoomed>threshold,0,zValsUnnormedZoomed)
            etaClean = np.sum(zValsClean,axis=0) # sum over to see eta distribution
            phiClean = np.sum(zValsClean,axis=1) # sum over to see phi distribution
            plt.figure()
            eta = np.sum(zValsUnnormedZoomed,axis=0)
            plt.plot(eta)
            plt.plot(etaClean)
            plt.savefig("{}/etaClean_{}.{}".format(plotFolder,var,plotFormat), bbox_inches="tight")
            plt.figure()
            phi = np.sum(zValsUnnormedZoomed,axis=1)
            plt.plot(phi)            
            plt.plot(phiClean)
            plt.savefig("{}/phiClean_{}.{}".format(plotFolder,var,plotFormat), bbox_inches="tight")
        hsDict[var] = [np.around(jiEtaHotSpots,2),np.around(jiPhiHotSpots,2)]
    phiSpikeHotSpots[year] = hsDict

# # combine 2018 PreHEM and 2018 PostHEM hotspots into one single year
preDict = phiSpikeHotSpots["2018PreHEM"]
postDict = phiSpikeHotSpots["2018PostHEM"]
phiSpikeHotSpots["2018"] = {}
for var, hotSpotList in preDict.items():
    jiEtaHotSpots, jiPhiHotSpots = uniqueCombine(preDict[var][0],preDict[var][1],postDict[var][0],postDict[var][1])
    phiSpikeHotSpots["2018"][var] = [np.around(jiEtaHotSpots,2),np.around(jiPhiHotSpots,2)]

if normalization:
    with open(f"{plotFolderPre}/phiSpikeHotSpots.pkl", "wb") as outfile:
        pickle.dump(phiSpikeHotSpots, outfile, protocol=pickle.HIGHEST_PROTOCOL)