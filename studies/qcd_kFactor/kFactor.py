import uproot as up 
import mplhep as hep
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

mpl.use('Agg')
mpl.rc("font", family="serif", size=16)
plotFormat = "pdf"



def getHisto(inputFile, varName, cut, rebin=False):
    f = up.open(inputFile)
    hist = f[varName + cut]
    npHist, bins = hist.to_numpy()
    if rebin:
        npHist, bins = rebin_histogram(npHist, bins, rebin)
    return npHist, bins

def rebin_histogram(np_hist, bins, rebin_factor):
    assert len(np_hist) % rebin_factor == 0, "Bins must be divisible by rebin factor"
    new_bins = bins[::rebin_factor]
    new_hist = np.add.reduceat(np_hist, np.arange(0, len(np_hist), rebin_factor))
    return new_hist, new_bins

def getBkgHisto(inputFile, varName, cut, rebin=False):
    npHist, bins = getHisto(inputFile, varName, cut, rebin)
    return npHist, bins

def getTotalEvents(inputFile, cut):
    # npHist, _ = getHisto(inputFile, "h_eCounter", cut)
    npHist, _ = getHisto(inputFile, "h_evtw", cut)
    return np.sum(npHist)

def absMeanError(numData, denData, avgRatio):
    ratio = numData / denData
    finiteRatio = ratio[np.isfinite(ratio)]
    # print(f"finiteRatio = {finiteRatio}, avgRatio = {avgRatio}, abs = {abs(finiteRatio - avgRatio)}")
    return np.nanmean(abs(finiteRatio - avgRatio))

def meanSquaredError(dataHist, bkgHist):
    diff = dataHist - bkgHist
    # print(f"datahist = {dataHist},  bkgHist = {bkgHist},      diff mse = {diff},  diff[diff >= 0] ** 2 = {diff[diff >= 0] ** 2}")
    return np.mean(diff ** 2)

def normalize(histVals):
    return histVals / np.sum(histVals)

def calculate_kFactor(dataFile, qcdFile, mcFiles, cut):
    # Get total event yields for Data, QCD, and other MC backgrounds
    data_yield = getTotalEvents(dataFile, cut)
    qcd_yield = getTotalEvents(qcdFile, cut)
    mc_yield = sum(getTotalEvents(mcFile, cut) for mcFile in mcFiles)

    if qcd_yield == 0:
        raise ValueError("QCD yield is zero, cannot calculate kFactor.")
    
    # Calculate kFactor directly
    kFactor = abs(data_yield - mc_yield) / qcd_yield

    # Error calculation
    data_error = np.sqrt(data_yield)
    qcd_error = np.sqrt(qcd_yield)
    mc_errors = [np.sqrt(getTotalEvents(mcFile, cut)) for mcFile in mcFiles]
    mc_yield_error = np.sqrt(sum(e**2 for e in mc_errors))  # Combined error for MC backgrounds

    # Propagate errors to get kFactor error
    kFactor_error = np.sqrt(
        (data_error / qcd_yield) ** 2 +
        (mc_yield_error / qcd_yield) ** 2 +
        ((data_yield - mc_yield) * qcd_error / qcd_yield**2) ** 2
    )

    return kFactor, kFactor_error


def generate_yield_table(inputFolder, years, samples, cut, rebin=None, output_file="yield_table.csv"):
    # Initialize the table dictionary
    table_data = {year: {} for year in years}

    for year in years:
        inputDataFile = f"{inputFolder}/{year}_Data.root"
        inputQCDFile = f"{inputFolder}/{year}_QCD.root"
        
        # Get Data yield
        data_yield = getTotalEvents(inputDataFile, cut)
        table_data[year]['Data'] = data_yield
        
        # Get yields for each MC background sample
        inputBkgSample = []
        for sample in samples:
            inputFile = f"{inputFolder}/{year}_{sample}.root"
            inputBkgSample.append(inputFile)
            sample_yield = getTotalEvents(inputFile, cut)
            table_data[year][sample] = sample_yield
        
        # Get QCD yield
        qcd_yield = getTotalEvents(inputQCDFile, cut)
        table_data[year]['QCD'] = qcd_yield
        
        # Calculate total MC yield
        total_mc_yield = sum(table_data[year][sample] for sample in samples) + qcd_yield
        table_data[year]['Total_MC'] = total_mc_yield
        
        # Calculate kFactor and kFactor error
        kFac, kFacError = calculate_kFactor(inputDataFile, inputQCDFile, inputBkgSample, cut)
        table_data[year]['kFactor'] = f"{kFac:.4f} ± {kFacError:.4f}"
    
    # Convert table data into a DataFrame for display
    row_labels = samples + ['QCD', 'Total_MC', 'Data', 'kFactor']
    df = pd.DataFrame({year: {key: table_data[year].get(key, '') for key in row_labels} for year in years})
    
    # Save to CSV
    df.to_csv(output_file, index_label="Category")
    print(f"Yield table saved as '{output_file}'")
    return df
    

inputFolder = "/uscms/home/nparmar/nobackup/SVJTchannel/t-channel_19_8_2024/t-channel_Analysis/keane_skim_outputs/skim_phiSpike_restudy_fixedTrgBug_pre_psFilterSig4"
kFacs = np.arange(1, 2, 0.1)
years = ["2016","2017","2018PreHEM","2018PostHEM"]
samples = ["ST", "ZJets", "WJets", "TTJets"]
cut = "_pre_psFilterSig4"
reBin = 10

varDict = {
    'h_MET': ["MET [GeV]", [200, 1500]],
    'h_j1PtAK8': ["j1 pT [GeV]", [200, 1500]],
    'h_j2PtAK8': ["j2 pT [GeV]", [200, 1500]],
    'h_j1PhiAK8': ["j1 Phi", [-np.pi, np.pi]],
    'h_j2PhiAK8': ["j2 Phi", [-np.pi, np.pi]],
}
kFactor_df = generate_yield_table(inputFolder, years,samples,cut)
for year in years:
    inputDataFile = f"{inputFolder}/{year}_Data.root"
    inputQCDFile = f"{inputFolder}/{year}_QCD.root"
    inputBkgSample = []
    for sample in samples:
        inputFile = f"{inputFolder}/{year}_{sample}.root"
        inputBkgSample.append(inputFile)
    kFac, kFacError = calculate_kFactor(inputDataFile, inputQCDFile, inputBkgSample, cut)
    kFac = 1.0
    print(f"For the year = {year} the optimized kFac = {kFac:0.4f} +- {kFacError:0.4f}")
    for varName, varDetail in varDict.items():
        outputFolder = f"./plots/{year}/{varName}"
        os.makedirs(outputFolder, exist_ok=True)    
        npHistData, bins = getBkgHisto(inputDataFile, varName, cut, rebin=reBin)
        ameList = []
        inputBkgSample = []
        histoStacks = []
        labels = []
        npHistQCD, bins = getBkgHisto(inputQCDFile, varName, cut, rebin=reBin)
        npHistQCD *= kFac
        totalNpHist = npHistQCD.copy()
        for sample in samples:
            inputFile = f"{inputFolder}/{year}_{sample}.root"
            inputBkgSample.append(inputFile)
            npHist, bins = getBkgHisto(inputFile, varName, cut, rebin=reBin)
            histoStacks.append(npHist)
            totalNpHist += npHist 
            labels.append(sample)
        histoStacks.append(npHistQCD)
        labels.append("QCD")

        # Stacked Plot
        # plt.figure()
        # hep.histplot(histoStacks, bins=bins, histtype="fill", stack=True, label=labels)
        # plt.yscale("log")
        # plt.ylim(1, 1e6)  # Set y-axis limit
        # plt.legend(ncol=2)  # Make 2 columns in the legend
        # plt.savefig(f"{outputFolder}/stacked_{np.around(kFac, 2)}.pdf")

        # Ratio Plot
        # Set CMS style
        hep.style.use("CMS")
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw=dict(height_ratios= [3, 1],hspace=0.1), sharex=True)
        hep.cms.text('Preliminary',loc=0,ax=ax1)
        hep.histplot([npHistData], bins=bins, ax=ax1, label="Data", histtype="errorbar", color="black")
        hep.histplot(histoStacks, bins=bins, ax=ax1, histtype="fill", stack=True, label=labels)
        ax1.set_yscale("log")
        ax1.set_ylim(1, 1e7)
        ax1.set_ylabel("Events")
        ax1.legend(ncol=2)  # Make 2 columns in the legend

        ratio = np.divide(npHistData, totalNpHist, out=np.zeros_like(npHistData), where=totalNpHist != 0)
        ax2.plot(bins[:-1], ratio, marker='o', linestyle='', color="black")
        ax2.set_ylim(0, 2)
        ax2.set_ylabel("Data / Bkg")
        ax2.set_xlabel(varDetail[0])
        ax2.axhline(1, linestyle="--", color="red")

        # Calculate AME and Normalized AME
        avgRatio = getTotalEvents(inputDataFile, cut) / getTotalEvents(inputQCDFile, cut)
        ame = absMeanError(totalNpHist, npHistQCD, avgRatio)
        normalizedTotalNpHist = normalize(totalNpHist)
        normalizedNpHistQCD = normalize(npHistQCD)
        normalized_ame = absMeanError(normalizedTotalNpHist,normalizedNpHistQCD,avgRatio)
        mse = meanSquaredError(npHistData, totalNpHist)
        ameList.append(ame)
        # print(f"getTotalEvents({inputDataFile}, {cut}) = {getTotalEvents(inputDataFile, cut)},  getTotalEvents({inputQCDFile}, {cut}) = { getTotalEvents(inputQCDFile, cut)}, avg ratio = {avgRatio}, totalNpHist = {totalNpHist}, npHistQCD = {npHistQCD}, ame = {ame}, mse = {mse}")

        # Add AME and Normalized AME to Legend
        # dummy_legend = [plt.Line2D([0], [0], color='none', label=f'kFactor used = {kFac:0.3f} +- {kFacError:0.3f}')]
        # ax2.legend(handles=[dummy_legend[0], ax2.lines[0]], loc="lower right", fontsize=15)

        plt.savefig(f"{outputFolder}/ratio_{np.around(kFac, 2)}.pdf")