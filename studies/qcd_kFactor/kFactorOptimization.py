# use coffea environment to do this
import uproot as up 
import mplhep as hep
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import coffea.hist as hist
import os

mpl.use('Agg')

mpl.rc("font", family="serif", size=16)

plotFormat = "pdf"

def normalize(histVals):
    return histVals/np.sum(histVals)

def getHisto(inputFile,varName,cut, rebin=False):
    f = up.open(inputFile)
    hist = f[varName + cut]
    npHist, bins = hist.to_numpy() 
    coffeaHist = hist.to_hist() #*(1./np.sum(npHist)) # useful for ratioPlot
    if rebin:
        npHist, bins = rebin_histogram(npHist, bins, rebin)
        coffeaHist = rebin_coffea_hist(coffeaHist, new_bins=2)
    return npHist, coffeaHist, bins

def rebin_coffea_hist(histogram, new_bins):
    if isinstance(new_bins, int):
        return histogram.rebin(histogram.axis('x').name, new_bins)  # Assuming 'x' is the axis name
    elif isinstance(new_bins, (list, np.ndarray)):
        return histogram.rebin(histogram.axis('x').name, hist.Bin(histogram.axis('x').name, new_bins))
    else:
        raise ValueError("new_bins must be either an integer or a list of bin edges.")


def getBkgHisto(inputFile,varName,cut,rebin=False):
    npHist, coffeaHist, bins = getHisto(inputFile,varName,cut,rebin)
    return coffeaHist, npHist, bins

def getTotalEvents(inputFile,cut):
    npHist, coffeaHist, bins = getHisto(inputFile,"h_eCounter",cut)
    return np.sum(npHist)

def absMeanError(numData,denData,avgRatio):
    ratio = numData/denData
    finiteRatio =  ratio[np.isfinite(ratio)]
    return np.nanmean(abs(finiteRatio - avgRatio))

def rebin_histogram(np_hist, bins, rebin_factor):
    assert len(np_hist) % rebin_factor == 0, "Bins must be divisible by rebin factor"

    # Rebin the histogram
    new_bins = bins[::rebin_factor]
    new_hist = np.add.reduceat(np_hist, np.arange(0, len(np_hist), rebin_factor))
    return new_hist, new_bins

def meanSquaredError(dataHist, bkgHist):
    # Calculate the MSE between data histogram and background histogram
    diff = dataHist - bkgHist
    return np.mean(diff[diff >= 0] ** 2)  # Only consider non-negative differences


import numpy as np

# def rebin1D(hist_vals, bin_edges, rebin_factor):
#     num_bins = len(hist_vals)
#     if num_bins % rebin_factor != 0:
#         raise ValueError(f"Number of bins ({num_bins}) must be divisible by rebin factor ({rebin_factor}).")
    
#     rebinned_vals = np.add.reduceat(hist_vals, np.arange(0, num_bins, rebin_factor))

#     rebinned_edges = bin_edges[::rebin_factor]
#     if rebinned_edges[-1] != bin_edges[-1]:
#         rebinned_edges = np.append(rebinned_edges, bin_edges[-1])  # Ensure the last edge is included

#     return rebinned_vals, rebinned_edges


inputFolder = "/uscms/home/nparmar/nobackup/SVJTchannel/t-channel_19_8_2024/t-channel_Analysis/keane_skim_outputs/skim_phiSpike_restudy_fixedTrgBug_pre_psFilterSig4"
kFacs = np.arange(1,10,1)
years = ["2017"]#,"2017","2018PreHEM","2018PostHEM"]

samples = ["ST","ZJets","WJets","TTJets"]
cut = "_pre_psFilterSig4"
reBin = False
varDict = {
    'h_MET':    ["MET [GeV]",[200, 1500]],}
#     'h_j1PtAK8':  ["MET [GeV]",[200, 1500]],
#     'h_j1PhiAK8': ["MET [GeV]",[200, 1500]],
#     'h_j2PtAK8':  ["MET [GeV]",[200, 1500]],
#     'h_j2PhiAK8': ["MET [GeV]",[200, 1500]],
# }


for year in years:
    for varName, varDetail in varDict.items():
        outputFolder = f"./plots/{year}/{varName}"
        os.makedirs(outputFolder,exist_ok=True)
        inputDataFile = f"{inputFolder}/{year}_Data.root"
        coffeaHistData, npHistData, bins = getBkgHisto(inputDataFile,varName,cut,rebin=reBin)
        inputQCDFile = f"{inputFolder}/{year}_QCD.root"
        ameList = []
        for kFac in kFacs:
            histoStacks = []
            labels = []
            coffeaHistQCD, npHistQCD, bins = getBkgHisto(inputQCDFile,varName,cut,rebin=reBin)
            coffeaHistQCD *= kFac
            totalHist = coffeaHistQCD.copy()
            totalNpHist = npHistQCD
            npHistQCD *= kFac 
            for sample in samples:
                inputFile = f"{inputFolder}/{year}_{sample}.root"
                coffeaHist, npHist, bins = getBkgHisto(inputFile,varName,cut,rebin = reBin)
                histoStacks.append(coffeaHist)
                totalHist += coffeaHist
                totalNpHist += npHist
                labels.append(sample)
            histoStacks.append(coffeaHistQCD)
            labels.append("QCD")
            plt.figure()
            hep.histplot(histoStacks,histtype="fill", stack=True, label=labels)
            plt.yscale("log")
            plt.legend()
            plt.savefig(f"{outputFolder}/stacked_{np.around(kFac,2)}.pdf")

            fig = plt.figure(figsize=(10, 8))
            main, sub = coffeaHistData.plot_ratio(totalHist,    
                rp_num_label="Data",
                rp_denom_label="Background",
                rp_uncert_draw_type="line",
                fp_label="Events",
                rp_ylim=[0,2])
            ax1,ax2 = fig.get_axes()
            ax1.set_yscale("log")
            # ax1.set_ylim(0.1,500000000000)
            # ax2.set_xlim(xDetail[1][0],xDetail[1][1])
            # ax2.set_ylim(0,2)
            # ax2.set_xlabel(xDetail[0])
            # mean absolute error
            avgRatio = getTotalEvents(inputDataFile,cut)/getTotalEvents(inputQCDFile,cut)
            ame = absMeanError(totalNpHist,npHistQCD,avgRatio)
            mse = meanSquaredError(npHistData, totalNpHist)
            ameList.append(ame)
            # ax2.axhline(avgRatio,linestyle="--",label="(Data yield)/(Bkg yield) = {:.2f}".format(avgRatio))
            # ax2.legend()
            # ax2.axhline(avgRatio, label=f"Ratio = {avgRatio:.2f}")
            ax2.legend(title=f'MSE: {mse:.2f}\n AME: {ame:.2f}')
            # print(ameList)
            plt.savefig(f"{outputFolder}/ratio_{np.around(kFac,2)}.pdf")
            # plt.savefig("{}{}vs{}_{}_data{}.{}".format(outputFolder,numLabel,denLabel,varName,cut,plotFormat))

# ameDict = {}
# for varName, xDetail in variables.items():
#     ameList = []
#     for i in range(len(etaLows)):
#         etaCutLabel = str(np.round(etaLows[i],2)).replace(".","p").replace("-","")
#         if noHEM:
#             cut = f"_pre"
#         else:
#             cut = f"_pre_eta_{etaCutLabel}"
#         coffeaHist, npHist, bins = getBkgHisto(numFile,varName,cut)
#         coffeaHistAPV, npHistAPV, bins = getBkgHisto(denFile,varName,cut)
#         maxInd = np.amax(np.where(npHist>0)[0])+1
#         fig = plt.figure(figsize=(10, 8))
#         main, sub = coffeaHist.plot_ratio(coffeaHistAPV,    
#             rp_num_label=numLabel,
#             rp_denom_label=denLabel,
#             rp_uncert_draw_type="line",
#             fp_label="Events",
#             rp_ylim=[0,2])
#         ax1,ax2 = fig.get_axes()
#         ax1.set_yscale("log")
#         ax1.set_ylim(0.1,500000000000)
#         ax2.set_xlim(xDetail[1][0],xDetail[1][1])
#         ax2.set_ylim(0,2)
#         ax2.set_xlabel(xDetail[0])
#         # mean absolute error
#         avgRatio = getTotalEvents(numFile,cut)/getTotalEvents(denFile,cut)
#         ame = absMeanError(npHist,npHistAPV,avgRatio)
#         ameList.append(ame)
#         ax2.hlines(avgRatio,np.min(bins),np.max(bins),linestyle="--",label="(preHEM yield)/(postHEM yield) = {:.2f}".format(avgRatio))
#         ax2.legend()
#         plt.savefig("{}{}vs{}_{}_data{}.{}".format(outputFolder,numLabel,denLabel,varName,cut,plotFormat))
#     ameDict[varName] = ameList

# print("ameDict")
# print(ameDict)

# if not noHEM:
#     plt.figure(figsize=(12,8))
#     for varName, ameScores in ameDict.items():
#         plt.plot(etaLows,ameScores/np.amax(ameScores),label=varName)
#         plt.ylabel("Normalized Difference between Actual and Expected Ratios")
#         plt.xlabel("$\eta$, $phi$ veto")
#         plt.legend()
#         plt.xticks(np.arange(-3.0, -3.51, -0.1), 
#             ('$-3.0<\eta<-1.4$\n$-1.57<\phi<-0.87$', '$-3.1<\eta<-1.3$\n$-1.67<\phi<-0.77$', '$-3.2<\eta<-1.2$\n$-1.77<\phi<-0.67$',   
#              '$-3.3<\eta<-1.1$\n$-1.87<\phi<-0.57$', '$-3.4<\eta<-1$\n$-1.97<\phi<-0.47$',  '$-3.5<\eta<-0.9$\n$-2.07<\phi<-0.37$',) ,fontsize = 12 )
#     plt.savefig("hemStudy.pdf")