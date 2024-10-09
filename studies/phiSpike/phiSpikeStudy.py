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
import utility as utl 

mpl.use('Agg')

mpl.rc("font", family="serif", size=20)

plotFormat = "pdf"

cut = "_pre_psFilterSig4p5" # this cut (if other than _pre) shows the effect of the phi spike filter for the particular nSigma
nSigma = 4.5 # this nSigma determines what the threshold for hot spot identification
etaMin = -2.4
etaMax = 2.4
phiMin = -3.2
phiMax = 3.2
qcdGt300 = True # only uses QCD samples with pT bin>300 GeV to avoid high weight events messing with the hot spot identification
nSigmaLab = str(np.round(nSigma,2)).replace(".","p")
if qcdGt300:
    plotFolderPre = f"plots/{cut[1:]}/sigma{nSigmaLab}/QCD_gt300"
else:
    plotFolderPre = f"plots/{cut[1:]}/sigma{nSigmaLab}/QCD_all"
inputFileDict = {
                    "2016": "/srv/skim_appliedPhiSpike_opt_DataMC_puWeight_clean/",
                    "2017": "/srv/skim_appliedPhiSpike_opt_DataMC_puWeight_clean/",
                    # "2018": "/srv/TM_phiSpike_opt_removeHighWeightQCD_clean/",
                    "2018PreHEM": "/srv/skim_appliedPhiSpike_opt_DataMC_puWeight_clean/",    # the QCD file used here is actually for all 2018, but shouldn't change the outcome
                    "2018PostHEM": "/srv/skim_appliedPhiSpike_opt_DataMC_puWeight_clean/",  # the QCD file used here is actually for all 2018, but shouldn't change the outcome
}
# samples = ["Data","QCD"]
samples = ["QCD"]
variables2D = {
                'j1Phivsj1Eta': [r"$\eta(j_1)$",r"$\phi(j_1)$"],
                'j2Phivsj2Eta': [r"$\eta(j_2)$",r"$\phi(j_2)$"],
                'j3Phivsj3Eta': [r"$\eta(j_3)$",r"$\phi(j_3)$"],
                'j4Phivsj4Eta': [r"$\eta(j_4)$",r"$\phi(j_4)$"],
}

for normalization in [True,False]:
    flatNessDict = {}
    phiSpikeHotSpots = {}
    for year, inputFolder in inputFileDict.items():
        hsDict = {}
        for var,labels in variables2D.items():
            jiEtaHotSpots = []
            jiPhiHotSpots = []
            for sample in samples:
                plotFolder = f"{plotFolderPre}/{year}/{sample}"
                if qcdGt300:
                    if sample == "QCD":
                        inputFile = f"{year}_{sample}_gt300.root"
                    else:
                        inputFile = f"{year}_{sample}_all.root"
                else:
                    inputFile = f"{year}_{sample}_all.root"
                os.makedirs(plotFolder,exist_ok=True)
                varName = "h_{}".format(var)
                zVals, xVals, yVals = utl.getTotal2DHist(inputFolder,inputFile,varName,cut)
                zVals = np.rot90(zVals,k=1)
                zVals,xVals,yVals = utl.rebin2D(zVals,xVals,yVals,2,2)
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
                    jiEtaHotSpot = utl.bintoActualVal(hotEtaVals, linFitX)
                    jiPhiHotSpot = utl.bintoActualVal(hotPhiVals, linFitY)
                    jiEtaHotSpots,jiPhiHotSpots = utl.uniqueCombine(jiEtaHotSpots,jiPhiHotSpots,jiEtaHotSpot,jiPhiHotSpot)
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

    # combine 2018 PreHEM and 2018 PostHEM hotspots into one single year
    if "2018PreHEM" in inputFileDict.keys() and "2018PostHEM" in inputFileDict.keys():
        preDict = phiSpikeHotSpots["2018PreHEM"]
        postDict = phiSpikeHotSpots["2018PostHEM"]
        phiSpikeHotSpots["2018"] = {}
        for var, hotSpotList in preDict.items():
            jiEtaHotSpots, jiPhiHotSpots = utl.uniqueCombine(preDict[var][0],preDict[var][1],postDict[var][0],postDict[var][1])
            phiSpikeHotSpots["2018"][var] = [np.around(jiEtaHotSpots,2),np.around(jiPhiHotSpots,2)]

    if normalization:
        with open(f"{plotFolderPre}/phiSpikeHotSpots.pkl", "wb") as outfile:
            pickle.dump(phiSpikeHotSpots, outfile, protocol=pickle.HIGHEST_PROTOCOL)