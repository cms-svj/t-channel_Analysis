# use coffea environment to do this
import uproot as up 
import mplhep as hep
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import coffea
import os
mpl.use('Agg')

mpl.rc("font", family="serif", size=16)

plotFormat = "pdf"

def normalize(histVals):
    return histVals/np.sum(histVals)

def getHisto(inputFile,varName,cut):
    f = up.open(inputFile)
    hist = f[varName + cut]
    npHist, bins = hist.to_numpy() 
    coffeaHist = hist.to_hist() #*(1./np.sum(npHist)) # useful for ratioPlot
    return npHist, coffeaHist, bins

def getBkgHisto(inputFile,varName,cut):
    npHist, coffeaHist, bins = getHisto(inputFile,varName,cut)
    return coffeaHist, npHist, bins

def getTotalEvents(inputFile,cut):
    npHist, coffeaHist, bins = getHisto(inputFile,"h_eCounter",cut)
    return np.sum(npHist)

def absMeanError(numData,denData,avgRatio):
    ratio = numData/denData
    finiteRatio =  ratio[np.isfinite(ratio)]
    return np.nanmean(abs(finiteRatio - avgRatio))


inputFolder = "/uscms/home/keanet/nobackup/SVJ/t-channel_skimBug/t-channel_Analysis/skim_all_preselection_DataMC_clean"
kFacs = np.arange(1,5.1,0.2)
years = ["2016","2017","2018PreHEM","2018PostHEM"]
samples = ["ST","ZJets","WJets","TTJets"]
cut = "_pre_psFilterSig3p5"

varDict = {
    'h_MET':    ["MET [GeV]",[200, 1500]],
    'h_j1PtAK8':  ["MET [GeV]",[200, 1500]],
    'h_j1PhiAK8': ["MET [GeV]",[200, 1500]],
    'h_j2PtAK8':  ["MET [GeV]",[200, 1500]],
    'h_j2PhiAK8': ["MET [GeV]",[200, 1500]],
}


for year in years:
    for varName, varDetail in varDict.items():
        outputFolder = f"./plots/{year}/{varName}"
        os.makedirs(outputFolder,exist_ok=True)
        inputDataFile = f"{inputFolder}/{year}_Data.root"
        coffeaHistData, npHistData, bins = getBkgHisto(inputDataFile,varName,cut)
        inputQCDFile = f"{inputFolder}/{year}_QCD.root"
        ameList = []
        for kFac in kFacs:
            histoStacks = []
            labels = []
            coffeaHistQCD, npHistQCD, bins = getBkgHisto(inputQCDFile,varName,cut)
            coffeaHistQCD *= kFac
            totalHist = coffeaHistQCD.copy()
            totalNpHist = npHistQCD
            npHistQCD *= kFac 
            for sample in samples:
                inputFile = f"{inputFolder}/{year}_{sample}.root"
                coffeaHist, npHist, bins = getBkgHisto(inputFile,varName,cut)
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
            ameList.append(ame)
            # ax2.axhline(avgRatio,linestyle="--",label="(Data yield)/(Bkg yield) = {:.2f}".format(avgRatio))
            ax2.legend()
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