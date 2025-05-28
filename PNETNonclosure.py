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
AddCMSText = True

def GetSVJbins(DNN_inner_edge,MET_inner_edge):
    '''Change the SVJbin edges here these are the ones for DNN trained on all back '''    
    SVJbins = {
                "0SVJ" : [DNN_inner_edge,MET_inner_edge],
                "1SVJ" : [DNN_inner_edge,MET_inner_edge],
                "2PSVJ" : [DNN_inner_edge,MET_inner_edge],
    }
    return SVJbins

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


def GetSubABCDregions(DNN_inner_edge,DNN_outer_edge,MET_inner_edge,MET_outer_edge):
    #Defining the A B C D regions
    print('MeT',MET_outer_edge)
    print('dnn',DNN_outer_edge)
    #print('dnn/3',dnn/2)
    regions = [ ("dA", MET_inner_edge, MET_outer_edge,DNN_inner_edge ,DNN_outer_edge),
                    ("dB", MET_inner_edge, MET_outer_edge, 0, DNN_inner_edge),
                    ("dC", 0, MET_inner_edge, DNN_inner_edge,DNN_outer_edge),
                    ("dD", 0, MET_inner_edge, 0, DNN_inner_edge)
                ]
    return regions

def adjustRegionBoundaries(region, xmin, xmax, ymin, ymax):
    """ Adjust regions to avoid bin overlap based on the detected overlap. """
    if region == 'dA':
        xmin = xmin + 0.1 
        ymin = ymin + 0.01
    if region == 'dB':
        xmin = xmin + 0.1
    if region == 'dC':
        ymin = ymin + 0.01    

    return xmin, xmax, ymin, ymax

def GetABCDhistPerSVJBin(data, ABCDhistoVar, maincut, SVJbin,DNN_inner_edge,DNN_outer_edge,MET_inner_edge,MET_outer_edge):
    """Return the dictionary of ABCD histogram for only one svj bin"""

    bkgname = data.fileName.split('_')[1].replace('.root','')
    hist_dict = {region: ROOT.TH1F(f"h_{bkgname}_{region}",f"h_{bkgname}_{region}", 1,0,1) for region in ['dA','dB','dC','dD']}
    SVJ, (dnn, met) = SVJbin
    histName = ABCDhistoVar + maincut + SVJ
    regions = GetSubABCDregions(DNN_inner_edge,DNN_outer_edge,MET_inner_edge,MET_outer_edge)


    for region, xmin, xmax, ymin, ymax in regions:   # TODO: this for condition should be written into a function work in all the cases
        xmin, xmax, ymin, ymax = adjustRegionBoundaries(region, xmin, xmax, ymin, ymax)
        hist, histIntegral, integral_error = data.get2DHistoIntegral(histName, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, showEvents=True)
        print(f"hist integral is - {hist.Integral()} in region - {region} and from the function is  - {histIntegral}, xmin - {xmin}, xmax - {xmax}, ymin - {ymin}, ymax - {ymax}")
        hist_dict[region].Sumw2()
        hist_dict[region].SetBinContent(1, histIntegral)
        hist_dict[region].SetBinError(1, integral_error)
    return hist_dict


def get_ABCD_counts_by_SVJ(Data,sgData, bgData, ABCDHistoVar, maincut, SVJBins,DNN_inner_edge,DNN_outer_edge,MET_inner_edge,MET_outer_edge):
    """
    Extracts the ABCD region counts split by the number of SVJ jets for both signal and background.
    """
    counts_by_SVJ = {}
    error_by_SVJ = {}
    for data_list, label in [(Data,"Data"),(sgData, "Signal"), (bgData, "Background")]:
        for data in data_list:  # Iterate over individual dataset objects
            for svj_bin in SVJBins.items():
                #hist_dict = GetABCDhistPerSVJBin(data, ABCDHistoVar, maincut, svj_bin)
                hist_dict = GetABCDhistPerSVJBin(data, ABCDHistoVar, maincut, svj_bin,DNN_inner_edge,DNN_outer_edge,MET_inner_edge,MET_outer_edge)
                print('hist_dict:',hist_dict)
                #put if statment with pre and if maincut is pre Data = Data SR SR cut s now pre
                #Data_SR = TF.GetABCDhistPerSVJBin(Data[0], ABCDHistoVar, SRcut, SVJbins)
                if hasattr(svj_bin, 'name'):
                    bin_label = svj_bin.name
                elif hasattr(svj_bin, 'label'):
                    bin_label = svj_bin.label
                else:
                    bin_label = str(svj_bin)  # Fallback
                
                N_A = hist_dict['dA'].GetBinContent(1)
                N_B = hist_dict['dB'].GetBinContent(1)
                N_C = hist_dict['dC'].GetBinContent(1)
                N_D = hist_dict['dD'].GetBinContent(1)
                
                N_A_error = hist_dict['dA'].GetBinError(1)
                N_B_error = hist_dict['dB'].GetBinError(1)
                N_C_error = hist_dict['dC'].GetBinError(1)
                N_D_error = hist_dict['dD'].GetBinError(1)
                
                if (bin_label, label) not in counts_by_SVJ:
                    counts_by_SVJ[(bin_label, label)] = (0, 0, 0, 0)

                if (bin_label, label) not in error_by_SVJ:
                    error_by_SVJ[(bin_label, label)] = (0, 0, 0, 0)
                
                counts_by_SVJ[(bin_label, label)] = (
                    counts_by_SVJ[(bin_label, label)][0] + N_A,
                    counts_by_SVJ[(bin_label, label)][1] + N_B,
                    counts_by_SVJ[(bin_label, label)][2] + N_C,
                    counts_by_SVJ[(bin_label, label)][3] + N_D,
                )

                error_by_SVJ[(bin_label, label)] = (
                    error_by_SVJ[(bin_label, label)][0] + N_A_error,
                    error_by_SVJ[(bin_label, label)][1] + N_B_error,
                    error_by_SVJ[(bin_label, label)][2] + N_C_error,
                    error_by_SVJ[(bin_label, label)][3] + N_D_error,
                )

    
    return counts_by_SVJ,error_by_SVJ

def compute_ABCD_prediction(Data, sgData, bgData, ABCDHistoVar,maincut, Year, DNN_inner_edges, DNN_outer_edges, MET_inner_edges, MET_outer_edges,VR):
    outer_edge_results = {}   
    error_results_by_outer_edge = ()
    # Separate lists for each SVJ category and each component (Data, Signal, Background)
    obs_0SVJ_data, obs_0SVJ_sig, obs_0SVJ_bg = [], [], []
    pred_0SVJ_data, pred_0SVJ_sig, pred_0SVJ_bg = [], [], []

    obs_1SVJ_data, obs_1SVJ_sig, obs_1SVJ_bg = [], [], []
    pred_1SVJ_data, pred_1SVJ_sig, pred_1SVJ_bg = [], [], []

    obs_2PSVJ_data, obs_2PSVJ_sig, obs_2PSVJ_bg = [], [], []
    pred_2PSVJ_data, pred_2PSVJ_sig, pred_2PSVJ_bg = [], [], []

    obserr_0SVJ_data, obserr_0SVJ_sig, obserr_0SVJ_bg = [], [], []
    prederr_0SVJ_data, prederr_0SVJ_sig, prederr_0SVJ_bg = [], [], []

    obserr_1SVJ_data, obserr_1SVJ_sig, obserr_1SVJ_bg = [], [], []
    prederr_1SVJ_data, prederr_1SVJ_sig, prederr_1SVJ_bg = [], [], []

    obserr_2PSVJ_data, obserr_2PSVJ_sig, obserr_2PSVJ_bg = [], [], []
    prederr_2PSVJ_data, prederr_2PSVJ_sig, prederr_2PSVJ_bg = [], [], []

    for DNN_inner_edge,DNN_outer_edge,MET_inner_edge,MET_outer_edge in zip(DNN_inner_edges,DNN_outer_edges,MET_inner_edges,MET_outer_edges):
        if VR == 'VR1':
            dnn = DNN_inner_edge
            met = MET_inner_edge
            print(f"\nProcessing outer_edge: {MET_outer_edge}")
        elif VR == 'VR2':
            dnn = DNN_inner_edge
            met = MET_inner_edge
            print(f"\nProcessing outer_edge: {DNN_outer_edge}")
        elif VR == 'VR3':
            dnn = DNN_inner_edge
            met = MET_inner_edge
            print(f"\nProcessing outer_edge: {DNN_outer_edge} and {MET_outer_edge} ")

        SVJBins = GetSVJbins(dnn,met)
        counts_by_SVJ, error_by_SVJ = get_ABCD_counts_by_SVJ(Data,sgData, bgData, ABCDHistoVar, maincut, SVJBins, DNN_inner_edge, DNN_outer_edge, MET_inner_edge, MET_outer_edge)
        outer_edge_results[(DNN_outer_edge, MET_outer_edge)] = counts_by_SVJ
            
        for ((svj_bin, label), (N_A, N_B, N_C, N_D)), ((_, _), (NA_err, NB_err, NC_err, ND_err)) in zip(counts_by_SVJ.items(), error_by_SVJ.items()):
            print(svj_bin)

            # Determine which SVJ category this bin belongs to
            if "0SVJ" in svj_bin:
                obs_data, obs_sig, obs_bg = obs_0SVJ_data, obs_0SVJ_sig, obs_0SVJ_bg
                pred_data, pred_sig, pred_bg = pred_0SVJ_data, pred_0SVJ_sig, pred_0SVJ_bg
                obs_dataerr, obs_sigerr, obs_bgerr = obserr_0SVJ_data, obserr_0SVJ_sig, obserr_0SVJ_bg
                pred_dataerr, pred_sigerr, pred_bgerr = prederr_0SVJ_data, prederr_0SVJ_sig, prederr_0SVJ_bg
            elif "1SVJ" in svj_bin:
                obs_data, obs_sig, obs_bg = obs_1SVJ_data, obs_1SVJ_sig, obs_1SVJ_bg
                pred_data, pred_sig, pred_bg = pred_1SVJ_data, pred_1SVJ_sig, pred_1SVJ_bg
                obs_dataerr, obs_sigerr, obs_bgerr = obserr_1SVJ_data, obserr_1SVJ_sig, obserr_1SVJ_bg
                pred_dataerr, pred_sigerr, pred_bgerr = prederr_1SVJ_data, prederr_1SVJ_sig, prederr_1SVJ_bg
            elif "2PSVJ" in svj_bin:
                obs_data, obs_sig, obs_bg = obs_2PSVJ_data, obs_2PSVJ_sig, obs_2PSVJ_bg
                pred_data, pred_sig, pred_bg = pred_2PSVJ_data, pred_2PSVJ_sig, pred_2PSVJ_bg
                obs_dataerr, obs_sigerr, obs_bgerr = obserr_2PSVJ_data, obserr_2PSVJ_sig, obserr_2PSVJ_bg
                pred_dataerr, pred_sigerr, pred_bgerr = prederr_2PSVJ_data, prederr_2PSVJ_sig, prederr_2PSVJ_bg
            else:
                print(f"Warning: No valid SVJ label found for {label}")
                continue  # Skip if no valid SVJ label is found

            # Compute predicted A
            N_A_pred = (N_B * N_C) / N_D if N_D > 0 else float('nan')
            relative_error = np.sqrt((NB_err / N_B if N_B > 0 else float('nan') ) ** 2 + (NC_err / N_C if N_C > 0 else float('nan')) ** 2 + (ND_err / N_D  if N_D > 0 else float('nan')) ** 2 )
            NA_pred_err = N_A_pred * relative_error
            
            # Identify whether this entry belongs to Data, Signal, or Background
            if "Data" in label:
                obs_data.append(N_A)
                pred_data.append(N_A_pred)
                obs_dataerr.append(NA_err)
                pred_dataerr.append(NA_pred_err)
            elif "Signal" in label:
                obs_sig.append(N_A)
                pred_sig.append(N_A_pred)
                obs_sigerr.append(NA_err)
                pred_sigerr.append(NA_pred_err)
            elif "Background" in label:
                obs_bg.append(N_A)
                pred_bg.append(N_A_pred)
                obs_bgerr.append(NA_err)
                pred_bgerr.append(NA_pred_err)

            print(f"[{label}] outer_edge: {dnn,met}, SVJ Bin: {svj_bin} | Observed A: {N_A},Observed ERR A: {NA_err},Predicted A: {N_A_pred}, Predicted Err {NA_pred_err}")
            #print(f"0SVJ obs len: {len(obs_0SVJ_data)}, 1SVJ obs len: {len(obs_1SVJ_data)}, 2PSVJ obs len: {len(obs_2PSVJ_data)}")

    return (
        outer_edge_results,
        obs_0SVJ_data, pred_0SVJ_data, obs_0SVJ_sig, pred_0SVJ_sig, obs_0SVJ_bg, pred_0SVJ_bg,
        obs_1SVJ_data, pred_1SVJ_data, obs_1SVJ_sig, pred_1SVJ_sig, obs_1SVJ_bg, pred_1SVJ_bg,
        obs_2PSVJ_data, pred_2PSVJ_data, obs_2PSVJ_sig, pred_2PSVJ_sig, obs_2PSVJ_bg, pred_2PSVJ_bg,
        obserr_0SVJ_data, obserr_0SVJ_sig, obserr_0SVJ_bg,prederr_0SVJ_data, prederr_0SVJ_sig, prederr_0SVJ_bg,
        obserr_1SVJ_data, obserr_1SVJ_sig, obserr_1SVJ_bg,prederr_1SVJ_data, prederr_1SVJ_sig, prederr_1SVJ_bg,
        obserr_2PSVJ_data, obserr_2PSVJ_sig, obserr_2PSVJ_bg,prederr_2PSVJ_data, prederr_2PSVJ_sig, prederr_2PSVJ_bg


    )

def plot_ABCD_ratios(
        year,
        outer_edge_results, 
        obs_0SVJ_data, pred_0SVJ_data, obs_0SVJ_sig, pred_0SVJ_sig, obs_0SVJ_bg, pred_0SVJ_bg,
        obs_1SVJ_data, pred_1SVJ_data, obs_1SVJ_sig, pred_1SVJ_sig, obs_1SVJ_bg, pred_1SVJ_bg,
        obs_2PSVJ_data, pred_2PSVJ_data, obs_2PSVJ_sig, pred_2PSVJ_sig, obs_2PSVJ_bg, pred_2PSVJ_bg,
        obserr_0SVJ_data, obserr_0SVJ_sig, obserr_0SVJ_bg,prederr_0SVJ_data, prederr_0SVJ_sig, prederr_0SVJ_bg,
        obserr_1SVJ_data, obserr_1SVJ_sig, obserr_1SVJ_bg,prederr_1SVJ_data, prederr_1SVJ_sig, prederr_1SVJ_bg,
        obserr_2PSVJ_data, obserr_2PSVJ_sig, obserr_2PSVJ_bg,prederr_2PSVJ_data, prederr_2PSVJ_sig, prederr_2PSVJ_bg,
        Boundary_vals, output_dir,title):
    """
    Plots the ratio of observed to predicted A values for each SVJ type and SCJ category separately,
    and adds a subplot showing the difference between Data Ratio and Background Ratio with a rectangular bottom plot.
    """

    # Define x-axis values
    values = np.array(Boundary_vals)

    # Compute ratios safely
    def compute_ratio(obs,obs_err,pred,pred_err):
        obs = np.array(obs) 
        pred = np.array(pred)  
        obs_err = np.array(obs_err)
        pred_err = np.array(pred_err)

        """Avoid division by zero errors"""
        ratio =  np.where(pred > 0, obs / pred, float('nan'))
        nonclosure = np.ones_like(ratio)-np.reciprocal(ratio)

        ratio_err = np.where(pred > 0, ratio * np.sqrt((obs_err / obs) ** 2 + (pred_err / pred) ** 2), float('nan'))
        nonclosure_err = (1 / ratio ** 2) * ratio_err
        #nonclosure_err = ratio_err
        return nonclosure,nonclosure_err

    def signalcontamination(sig_obs,sig_obs_err,bckg_obs,bckg_err):
        sig_obs = np.array(sig_obs) 
        bckg_obs = np.array(bckg_obs)  
        sig_obs_err = np.array(sig_obs_err)
        bckg_err = np.array(bckg_err)
        MC = sig_obs+bckg_obs
        MCerr = np.sqrt((sig_obs_err**2)+(bckg_err**2))
        sigcontam =  np.where(bckg_obs > 0, sig_obs / MC, float('nan'))
        sigcont_ratio_err = np.where(MC > 0, sigcontam * np.sqrt((sig_obs_err / sig_obs) ** 2 + (MCerr / MC) ** 2), float('nan'))
        return sigcontam,sigcont_ratio_err
        
    def plot_signal_contamination(x_values, sig_contam, sig_contam_err, label, filename):
        hep.cms.label(rlabel="")
        plt.figure(figsize=(10, 6))
        plt.errorbar(x_values, sig_contam, yerr=sig_contam_err, fmt='o', color='purple', label=label, capsize=5)
        plt.axhline(0, linestyle='dashed', color='black', linewidth=2)
        hep.cms.label(rlabel="")
        plt.xlabel("Boundary Value",fontsize=16)
        plt.ylabel("Signal Contamination",fontsize=16)
        plt.title(f'Signal Contamination {label} {year}', fontsize=15.5, loc = 'right')
        plt.grid(True)
        plt.legend()
        plt.savefig(os.path.join(output_dir, filename), dpi=300)
        plt.close()

    # Compute non-closure and errors for each category
    ratio_0SVJ_data, errbars_0SVJ_data = compute_ratio(obs_0SVJ_data,obserr_0SVJ_data, pred_0SVJ_data,prederr_0SVJ_data)
    ratio_0SVJ_bg, errbars_0SVJ_bg = compute_ratio(obs_0SVJ_bg,obserr_0SVJ_bg, pred_0SVJ_bg,prederr_0SVJ_bg)
    ratio_0SVJ_sig, errbars_0SVJ_sig = compute_ratio(obs_0SVJ_sig,obserr_0SVJ_sig,prederr_0SVJ_sig ,pred_0SVJ_sig)

    ratio_1SVJ_data, errbars_1SVJ_data = compute_ratio(obs_1SVJ_data,obserr_1SVJ_data, pred_1SVJ_data,prederr_1SVJ_data)
    ratio_1SVJ_bg, errbars_1SVJ_bg = compute_ratio(obs_1SVJ_bg,obserr_1SVJ_bg, pred_1SVJ_bg,prederr_1SVJ_bg)
    ratio_1SVJ_sig, errbars_1SVJ_sig = compute_ratio(obs_1SVJ_sig,obserr_1SVJ_sig,prederr_1SVJ_sig ,pred_1SVJ_sig)

    ratio_2PSVJ_data, errbars_2PSVJ_data = compute_ratio(obs_2PSVJ_data,obserr_2PSVJ_data, pred_2PSVJ_data,prederr_2PSVJ_data)
    ratio_2PSVJ_bg, errbars_2PSVJ_bg = compute_ratio(obs_2PSVJ_bg,obserr_2PSVJ_bg, pred_2PSVJ_bg,prederr_2PSVJ_bg)
    ratio_2PSVJ_sig, errbars_2PSVJ_sig = compute_ratio(obs_2PSVJ_sig,obserr_2PSVJ_sig,prederr_2PSVJ_sig ,pred_2PSVJ_sig)

    # Compute Data Ratio - Background Ratio
    diff_0SVJ_data_bg = ratio_0SVJ_data - ratio_0SVJ_bg
    diff_1SVJ_data_bg = ratio_1SVJ_data - ratio_1SVJ_bg
    diff_2PSVJ_data_bg = ratio_2PSVJ_data - ratio_2PSVJ_bg
    
    errbars_0SVJ_diff = np.sqrt(errbars_0SVJ_data**2 + errbars_0SVJ_bg**2)
    errbars_1SVJ_diff = np.sqrt(errbars_1SVJ_data**2 + errbars_1SVJ_bg**2)
    errbars_2PSVJ_diff = np.sqrt(errbars_2PSVJ_data**2 + errbars_2PSVJ_bg**2)

    #compute Signal contamination 
    sig_contam_0SVJ, sig_contam_err_0SVJ = signalcontamination(obs_0SVJ_sig, obserr_0SVJ_sig, obs_0SVJ_bg, obserr_0SVJ_bg)
    sig_contam_1SVJ, sig_contam_err_1SVJ = signalcontamination(obs_1SVJ_sig, obserr_1SVJ_sig, obs_1SVJ_bg, obserr_1SVJ_bg)
    sig_contam_2PSVJ, sig_contam_err_2PSVJ = signalcontamination(obs_2PSVJ_sig, obserr_2PSVJ_sig, obs_2PSVJ_bg, obserr_2PSVJ_bg)

    # Function to plot and save the main plot and the difference plot
    def plot_and_save(x_values,non_closures,nonclosure_err, diff_ratios,diff_error, labels, colors,title, filename):
        plt.figure(figsize=(15, 12))  # Increase figure size to make room for the rectangular bottom plot
        
        # Create a gridspec for controlling subplot layout
        gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])  # Top plot takes 2/3 of the space, bottom plot 1/3

        # Plot the main ratio plot (top subplot)
        ax0 = plt.subplot(gs[0])  
        for non_closure, error, label, color in zip(non_closures, nonclosure_err, labels, colors):
            if len(non_closure) > 0:  
                ax0.errorbar(x_values, non_closure, yerr=error, fmt='o', color=color, label=label, capsize=5)
                ax0.plot(x_values, np.zeros_like(x_values), linestyle='dashed', color='black', linewidth=2)
        

        hep.cms.label(rlabel="")

        ax0.set_xlabel("Boundary Value")
        ax0.set_ylabel("Non-Closure")
        ax0.set_title(f'{title}', fontsize=18)
        ax0.legend()
        ax0.set_ylim(-0.40, 0.40)
        ax0.grid(True)

        # Plot the difference plot (bottom subplot)
        ax1 = plt.subplot(gs[1])  # Bottom subplot
        for diff_ratio,error, label, color in zip(diff_ratios,diff_error,labels, colors):
            if len(diff_ratio) > 0:  # Avoid empty plots
                ax1.errorbar(x_values,diff_ratio, yerr=error, fmt='o', color='black', label='data-background',capsize=5)
                ax1.plot(x_values, np.zeros_like(x_values), linestyle='dashed', color='black', linewidth=2)
        
        ax1.set_xlabel("Boundary Value")
        ax1.set_ylabel("Data - Background Sim ")
        ax1.set_title(f'Difference in Data and Background Ratios  {year}', fontsize=16)
        ax1.set_ylim(-0.60, 0.60)
       #ax1.legend()
        ax1.grid(True)

        for ax in [ax0, ax1]:
            for spine in ax.spines.values():
                spine.set_linewidth(2)  # Make borders bold

        ax0.yaxis.set_major_locator(MultipleLocator(0.1))  # Major ticks
        ax0.yaxis.set_minor_locator(MultipleLocator(0.05))  # Minor ticks

        #ax1.yaxis.set_major_locator(MultipleLocator(0.1))
        #ax1.yaxis.set_minor_locator(MultipleLocator(0.05))

        ax0.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
        # Save the plot
        plt.tight_layout()  # Adjust layout to prevent overlap
        plt.savefig(os.path.join(output_dir, filename), dpi=300)
        plt.close()

    plot_and_save(values, 
                [ratio_0SVJ_data, ratio_0SVJ_bg],
                [errbars_0SVJ_data, errbars_0SVJ_bg],
                [diff_0SVJ_data_bg],[errbars_0SVJ_diff], 
                ["Data 0SVJ", "Background MC 0SVJ"], 
                ["b", "r"],
                title, 
                "Ratio_Data_0SVJ.pdf")
                #"Ratio_Data_0SVJ.jpg")

    plot_and_save(values, 
                [ratio_1SVJ_data, ratio_1SVJ_bg],
                [errbars_1SVJ_data, errbars_1SVJ_bg],
                [diff_1SVJ_data_bg],[errbars_1SVJ_diff], 
                ["Data 1SVJ", "Background MC 1SVJ"], 
                ["b", "r"], 
                title,
                "Ratio_Data_1SVJ.pdf")
                #"Ratio_Data_1SVJ.jpg")
    print(f"data {ratio_1SVJ_data} err bars {errbars_1SVJ_data}")
    plot_and_save(values, 
                [ratio_2PSVJ_data, ratio_2PSVJ_bg],
                [errbars_2PSVJ_data, errbars_2PSVJ_bg],
                [diff_2PSVJ_data_bg],[errbars_2PSVJ_diff], 
                ["Data 2SVJ", "Background MC 2SVJ"], 
                ["b", "r"], 
                title,
                "Ratio_Data_2PSVJ.pdf")
                #"Ratio_Data_2PSVJ.jpg")
    
    #print(f"non closure data 2psvj {ratio_2PSVJ_data} ")
    plot_signal_contamination(values, sig_contam_0SVJ, sig_contam_err_0SVJ, "0SVJ", "signal_contamination_0SVJ.jpg")
    plot_signal_contamination(values, sig_contam_1SVJ, sig_contam_err_1SVJ, "1SVJ", "signal_contamination_1SVJ.jpg")
    plot_signal_contamination(values, sig_contam_2PSVJ, sig_contam_err_2PSVJ, "2PSVJ", "signal_contamination_2PSVJ.jpg")

    #print(f'Boundary Values {values[::-1]}')
    print(f'Boundary Values {values}')

def main():
    parser = optparse.OptionParser("usage: %prog [options]\n")
    parser.add_option('-d', '--dataset',    dest='dataset',                    default='testHadd_11242020',    help='dataset')
    parser.add_option('-y',                 dest='year',       type='string',  default='2018',                 help="Can pass in the run year")
    parser.add_option('-o',                 dest='outputdir',  type='string',                                  help="Output folder name")
    parser.add_option(    '--hemPeriod',  dest='hemPeriod', type=str, default=False,  help='HEM period (PreHEM or PostHEM), default includes entire sample',)
    options, args = parser.parse_args()
    print("Parsed options:", options)

    year = options.year
    hemPeriod = options.hemPeriod
    ABCDhistoVars = ["METvsDNN"]
    ABCDFolderName = "ABCD"
    CRCuts = "_cr_electron_"         
    maincuts = CRCuts
    
    SRCut = "_pre_"

    Years = ['2016','2017','2018']
    for Year in Years:
        Data, sgData, bgData = getData( options.dataset + "/", 1.0, Year)
        VR = 'VR1'
        output_dir = f'Nonclosure/PNET/VRI/{Year}/ControlRegion'
        os.makedirs(output_dir, exist_ok=True)
        DNN_outer_edges = np.linspace(1,1,30)
        DNN_inner_edges = np.linspace(0.6,0.6,30)
        MET_outer_edges = np.linspace(225,250,30)
        MET_inner_edges = np.linspace(210,225,30)
        outer_edge_results,obs_0SVJ_data, pred_0SVJ_data, obs_0SVJ_sig, pred_0SVJ_sig, obs_0SVJ_bg, pred_0SVJ_bg,obs_1SVJ_data, pred_1SVJ_data, obs_1SVJ_sig, pred_1SVJ_sig, obs_1SVJ_bg, pred_1SVJ_bg,obs_2PSVJ_data, pred_2PSVJ_data, obs_2PSVJ_sig, pred_2PSVJ_sig, obs_2PSVJ_bg, pred_2PSVJ_bg,obserr_0SVJ_data, obserr_0SVJ_sig, obserr_0SVJ_bg,prederr_0SVJ_data, prederr_0SVJ_sig, prederr_0SVJ_bg,obserr_1SVJ_data, obserr_1SVJ_sig, obserr_1SVJ_bg,prederr_1SVJ_data, prederr_1SVJ_sig, prederr_1SVJ_bg,obserr_2PSVJ_data, oberr_2PSVJ_sig, obserr_2PSVJ_bg,prederr_2PSVJ_data, prederr_2PSVJ_sig, prederr_2PSVJ_bg = compute_ABCD_prediction(Data, sgData, bgData, "h_METvsDNN", '_pre_',Year,DNN_inner_edges, DNN_outer_edges,MET_inner_edges,MET_outer_edges,VR)
        plot_ABCD_ratios(
            Year,
            outer_edge_results, 
            obs_0SVJ_data, pred_0SVJ_data, obs_0SVJ_sig, pred_0SVJ_sig, obs_0SVJ_bg, pred_0SVJ_bg,
            obs_1SVJ_data, pred_1SVJ_data, obs_1SVJ_sig, pred_1SVJ_sig, obs_1SVJ_bg, pred_1SVJ_bg,
            obs_2PSVJ_data, pred_2PSVJ_data, obs_2PSVJ_sig, pred_2PSVJ_sig, obs_2PSVJ_bg, pred_2PSVJ_bg,
            obserr_0SVJ_data, obserr_0SVJ_sig, obserr_0SVJ_bg,prederr_0SVJ_data, prederr_0SVJ_sig, prederr_0SVJ_bg,
            obserr_1SVJ_data, obserr_1SVJ_sig, obserr_1SVJ_bg,prederr_1SVJ_data, prederr_1SVJ_sig, prederr_1SVJ_bg,
            obserr_2PSVJ_data, oberr_2PSVJ_sig, obserr_2PSVJ_bg,prederr_2PSVJ_data, prederr_2PSVJ_sig, prederr_2PSVJ_bg,
            MET_outer_edges, output_dir,f'Control Region Non Closure split by MET (VR I) {Year}')

        VR = 'VR2'
        output_dir = f'Nonclosure/PNET/VRII/{Year}/ControlRegion'
        os.makedirs(output_dir, exist_ok=True)
        DNN_outer_edges = np.linspace(0.3,0.85,30)
        DNN_inner_edges = np.linspace(0.1,0.6,30)
        MET_outer_edges = np.linspace(1000,1000,30)
        MET_inner_edges = np.linspace(250,250,30)
        outer_edge_results,obs_0SVJ_data, pred_0SVJ_data, obs_0SVJ_sig, pred_0SVJ_sig, obs_0SVJ_bg, pred_0SVJ_bg,obs_1SVJ_data, pred_1SVJ_data, obs_1SVJ_sig, pred_1SVJ_sig, obs_1SVJ_bg, pred_1SVJ_bg,obs_2PSVJ_data, pred_2PSVJ_data, obs_2PSVJ_sig, pred_2PSVJ_sig, obs_2PSVJ_bg, pred_2PSVJ_bg,obserr_0SVJ_data, obserr_0SVJ_sig, obserr_0SVJ_bg,prederr_0SVJ_data, prederr_0SVJ_sig, prederr_0SVJ_bg,obserr_1SVJ_data, obserr_1SVJ_sig, obserr_1SVJ_bg,prederr_1SVJ_data, prederr_1SVJ_sig, prederr_1SVJ_bg,obserr_2PSVJ_data, oberr_2PSVJ_sig, obserr_2PSVJ_bg,prederr_2PSVJ_data, prederr_2PSVJ_sig, prederr_2PSVJ_bg = compute_ABCD_prediction(Data, sgData, bgData, "h_METvsDNN", '_pre_',Year,DNN_inner_edges, DNN_outer_edges,MET_inner_edges,MET_outer_edges,VR)
        plot_ABCD_ratios(
            Year,
            outer_edge_results, 
            obs_0SVJ_data, pred_0SVJ_data, obs_0SVJ_sig, pred_0SVJ_sig, obs_0SVJ_bg, pred_0SVJ_bg,
            obs_1SVJ_data, pred_1SVJ_data, obs_1SVJ_sig, pred_1SVJ_sig, obs_1SVJ_bg, pred_1SVJ_bg,
            obs_2PSVJ_data, pred_2PSVJ_data, obs_2PSVJ_sig, pred_2PSVJ_sig, obs_2PSVJ_bg, pred_2PSVJ_bg,
            obserr_0SVJ_data, obserr_0SVJ_sig, obserr_0SVJ_bg,prederr_0SVJ_data, prederr_0SVJ_sig, prederr_0SVJ_bg,
            obserr_1SVJ_data, obserr_1SVJ_sig, obserr_1SVJ_bg,prederr_1SVJ_data, prederr_1SVJ_sig, prederr_1SVJ_bg,
            obserr_2PSVJ_data, oberr_2PSVJ_sig, obserr_2PSVJ_bg,prederr_2PSVJ_data, prederr_2PSVJ_sig, prederr_2PSVJ_bg,
            DNN_outer_edges, output_dir,f'Control Region Non Closure split by DNN (VR II) {Year}')

        VR = 'VR3'
        output_dir = f'Nonclosure/PNET/VRIII/{Year}/ControlRegion'
        os.makedirs(output_dir, exist_ok=True)
        DNN_outer_edges = np.linspace(0.3,0.85,30)
        DNN_inner_edges = np.linspace(0.1,0.6,30)
        MET_outer_edges = np.linspace(225,250,30)
        MET_inner_edges = np.linspace(210,225,30)
        outer_edge_results,obs_0SVJ_data, pred_0SVJ_data, obs_0SVJ_sig, pred_0SVJ_sig, obs_0SVJ_bg, pred_0SVJ_bg,obs_1SVJ_data, pred_1SVJ_data, obs_1SVJ_sig, pred_1SVJ_sig, obs_1SVJ_bg, pred_1SVJ_bg,obs_2PSVJ_data, pred_2PSVJ_data, obs_2PSVJ_sig, pred_2PSVJ_sig, obs_2PSVJ_bg, pred_2PSVJ_bg,obserr_0SVJ_data, obserr_0SVJ_sig, obserr_0SVJ_bg,prederr_0SVJ_data, prederr_0SVJ_sig, prederr_0SVJ_bg,obserr_1SVJ_data, obserr_1SVJ_sig, obserr_1SVJ_bg,prederr_1SVJ_data, prederr_1SVJ_sig, prederr_1SVJ_bg,obserr_2PSVJ_data, oberr_2PSVJ_sig, obserr_2PSVJ_bg,prederr_2PSVJ_data, prederr_2PSVJ_sig, prederr_2PSVJ_bg = compute_ABCD_prediction(Data, sgData, bgData, "h_METvsDNN", '_pre_',Year,DNN_inner_edges, DNN_outer_edges,MET_inner_edges,MET_outer_edges,VR)
        plot_ABCD_ratios(
            Year,
            outer_edge_results, 
            obs_0SVJ_data, pred_0SVJ_data, obs_0SVJ_sig, pred_0SVJ_sig, obs_0SVJ_bg, pred_0SVJ_bg,
            obs_1SVJ_data, pred_1SVJ_data, obs_1SVJ_sig, pred_1SVJ_sig, obs_1SVJ_bg, pred_1SVJ_bg,
            obs_2PSVJ_data, pred_2PSVJ_data, obs_2PSVJ_sig, pred_2PSVJ_sig, obs_2PSVJ_bg, pred_2PSVJ_bg,
            obserr_0SVJ_data, obserr_0SVJ_sig, obserr_0SVJ_bg,prederr_0SVJ_data, prederr_0SVJ_sig, prederr_0SVJ_bg,
            obserr_1SVJ_data, obserr_1SVJ_sig, obserr_1SVJ_bg,prederr_1SVJ_data, prederr_1SVJ_sig, prederr_1SVJ_bg,
            obserr_2PSVJ_data, oberr_2PSVJ_sig, obserr_2PSVJ_bg,prederr_2PSVJ_data, prederr_2PSVJ_sig, prederr_2PSVJ_bg,
            DNN_outer_edges, output_dir,f'Control Region Non Closure split by DNN and MET(VR III) {Year}')

    # outer_edges = np.linspace(225,250,30)
    # inner_edges = np.linspace(210,225,30)

    # output_dir = 'Nonclosure/VRI-MET/ControlRegion'
    # outer_edge_results,obs_0SVJ_data, pred_0SVJ_data, obs_0SVJ_sig, pred_0SVJ_sig, obs_0SVJ_bg, pred_0SVJ_bg,obs_1SVJ_data, pred_1SVJ_data, obs_1SVJ_sig, pred_1SVJ_sig, obs_1SVJ_bg, pred_1SVJ_bg,obs_2PSVJ_data, pred_2PSVJ_data, obs_2PSVJ_sig, pred_2PSVJ_sig, obs_2PSVJ_bg, pred_2PSVJ_bg,obserr_0SVJ_data, obserr_0SVJ_sig, obserr_0SVJ_bg,prederr_0SVJ_data, prederr_0SVJ_sig, prederr_0SVJ_bg,obserr_1SVJ_data, obserr_1SVJ_sig, obserr_1SVJ_bg,prederr_1SVJ_data, prederr_1SVJ_sig, prederr_1SVJ_bg,obserr_2PSVJ_data, oberr_2PSVJ_sig, obserr_2PSVJ_bg,prederr_2PSVJ_data, prederr_2PSVJ_sig, prederr_2PSVJ_bg = compute_ABCD_prediction(Data, sgData, bgData, "h_METvsDNN", '_pre_',inner_edges, outer_edges)
    # #print('outer_edge_results ',outer_edge_results)
    # #print('observed_A',observed_A)
    # #print('ratio_A',ratio_A)
    # plot_ABCD_ratios(
    #     year,
    #     outer_edge_results, 
    #     obs_0SVJ_data, pred_0SVJ_data, obs_0SVJ_sig, pred_0SVJ_sig, obs_0SVJ_bg, pred_0SVJ_bg,
    #     obs_1SVJ_data, pred_1SVJ_data, obs_1SVJ_sig, pred_1SVJ_sig, obs_1SVJ_bg, pred_1SVJ_bg,
    #     obs_2PSVJ_data, pred_2PSVJ_data, obs_2PSVJ_sig, pred_2PSVJ_sig, obs_2PSVJ_bg, pred_2PSVJ_bg,
    #     obserr_0SVJ_data, obserr_0SVJ_sig, obserr_0SVJ_bg,prederr_0SVJ_data, prederr_0SVJ_sig, prederr_0SVJ_bg,
    #     obserr_1SVJ_data, obserr_1SVJ_sig, obserr_1SVJ_bg,prederr_1SVJ_data, prederr_1SVJ_sig, prederr_1SVJ_bg,
    #     obserr_2PSVJ_data, oberr_2PSVJ_sig, obserr_2PSVJ_bg,prederr_2PSVJ_data, prederr_2PSVJ_sig, prederr_2PSVJ_bg,
    #     outer_edges, output_dir)
 

if __name__ == '__main__':
    main()



