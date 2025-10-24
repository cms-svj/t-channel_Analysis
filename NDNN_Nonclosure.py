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
import matplotlib.colors as mcolors
from itertools import product
AddCMSText = True

def GetSVJbins(DNN_inner_edge,MET_inner_edge):
    '''Change the SVJbin edges here these are the ones for DNN trained on all back '''    
    SVJbins = {
                "0SVJ" : [DNN_inner_edge,MET_inner_edge],
                "1SVJ" : [DNN_inner_edge,MET_inner_edge],
                "2SVJ" : [DNN_inner_edge,MET_inner_edge],
                '2PSVJ': [DNN_inner_edge,MET_inner_edge],
                "3PSVJ": [DNN_inner_edge,MET_inner_edge],
    }
    return SVJbins

def getData(path, scale=1.0, year = "2018"):
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
        #info.DataSetInfo(basedir=path, fileName="2018_allBkg.root",               label="all bckg",                        scale=scale, color=ROOT.TColor.GetColor("#9c9ca1")),
    ]
    sgData = [
       
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="baseline", scale=scale, color=ROOT.kBlack),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",     label="mMed_600",  scale=scale, color=ROOT.kViolet+2),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-800_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",     label="mMed_800",  scale=scale, color=ROOT.kRed),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-1000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="mMed_1000", scale=scale, color=ROOT.kMagenta + 1),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-1500_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="mMed_1500", scale=scale, color=ROOT.kGray+4),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="mMed_3000", scale=scale, color=ROOT.kCyan),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="mMed_4000", scale=scale, color=ROOT.kBlue + 1),
        info.DataSetInfo(basedir=path, fileName=year+"_m2000_d20_r0p3_y1_N-1_M0_.root",    label="baseline", scale=scale, color=ROOT.kBlue + 1),
        # info.DataSetInfo(basedir=path, fileName=year+"_m600_d20_r0p3_y1_N-1_M0_.root",     label="mMed_600",  scale=scale, color=ROOT.kViolet+2),
        # info.DataSetInfo(basedir=path, fileName=year+"_m800_d20_r0p3_y1_N-1_M0_.root",     label="mMed_800",  scale=scale, color=ROOT.kRed),
        # info.DataSetInfo(basedir=path, fileName=year+"_m1000_d20_r0p3_y1_N-1_M0_.root",    label="mMed_1000", scale=scale, color=ROOT.kMagenta + 1),
        # info.DataSetInfo(basedir=path, fileName=year+"_m1500_d20_r0p3_y1_N-1_M0_.root",    label="mMed_1500", scale=scale, color=ROOT.kGray+4),
        # #info.DataSetInfo(basedir=path, fileName=year+"_m3000_d20_r0p3_y1_N-1_M0_.root",    label="mMed_3000", scale=scale, color=ROOT.kCyan),
        # info.DataSetInfo(basedir=path, fileName=year+"_m4000_d20_r0p3_y1_N-1_M0_.root",    label="mMed_4000", scale=scale, color=ROOT.kBlue + 1),
    ]
    return Data, sgData, bgData


def GetTFhistoAndPlot(ABCDhistDict_SR, ABCDhistDict_CR, SVJBins, xTitle="nSVJ", yTitle="Events", xmin=999.9, xmax=-999.9, year="2018", isLogY=False, saveName="TF.png",hemPeriod=False):
    ROOT.TH1.AddDirectory(False)
    skip_list = [f"{year}_QCD.root",f"{year}_ZJets.root",f"{year}_Data.root"]
    summed_SR = TF.SumHistograms(ABCDhistDict_SR, skip_list)
    summed_CR = TF.SumHistograms(ABCDhistDict_CR, skip_list)
    TFhist = pltutils.RatioHistogram(summed_SR,summed_CR)
    c1 = ROOT.TCanvas( "c", "c", 800, 700)
    c1, pad1, pad2 = pltutils.createCanvasPads(c1,isLogY)
    pad1.cd()
    ROOT.gStyle.SetOptStat("")
    leg = pltutils.SetupLegend(x1=0.4, NColumns=2)

    # make dummy to setup axis
    dummy = ROOT.TH1D("dummy", "dummy", len(SVJBins)*4, summed_SR.GetBinLowEdge(1), summed_SR.GetBinLowEdge(summed_SR.GetNbinsX()) + summed_SR.GetBinWidth(summed_SR.GetNbinsX()))
    # print(f"The dummy values are  - {summed_SR.GetBinLowEdge(1)}, {summed_SR.GetBinLowEdge(summed_SR.GetNbinsX())}, { summed_SR.GetBinWidth(summed_SR.GetNbinsX())}")
    ymax=10*9
    ymin=10
    lmax=10**9
    pltutils.setupDummy(dummy,leg,"", xTitle, yTitle, isLogY, xmin, xmax, ymin, ymax, lmax, isRatio=True)
    dummy.Draw("hist")
    
    pltutils.SetupLineHistStyle(summed_SR)
    pltutils.SetupLineHistStyle(summed_CR, color=ROOT.kRed)
    summed_SR.Draw("histe same")
    summed_CR.Draw("histe same")
    leg.AddEntry(summed_CR,"Control Region","L")
    leg.AddEntry(summed_SR,"Signal Region","L")

    lines_upperpad = pltutils.AddVerticalLine(summed_SR, SVJBins, ymax = 10**5)
    for line in lines_upperpad:
        line.Draw("same")
    
    pltutils.AddLabelsForABCD(summed_SR,SVJBins,yloc=5*10**4)
    leg.Draw("same")
    pad2.cd()
    TFhist = pltutils.SetupRatioStyle(TFhist, xTitle="nSVJ", yTitle="SR/CR", yTitleSize=0.13, ymax=2)
    TFhist.Draw("EX0P")
    lines_lowerpad = pltutils.AddVerticalLine(summed_SR, SVJBins, ymax = 2)
    for line in lines_lowerpad:
        line.Draw("same")
    if AddCMSText:
        pltutils.AddCMSLumiText(c1, year,isExtraText=True,hemPeriod=hemPeriod)
    
    # print(f"Summed SR Hist -  {pltutils.printBinContentAndError(summed_SR)}")
    # print(f"Summed CR Hist - {pltutils.printBinContentAndError(summed_CR)}")
    # print(f"TF Hist - {pltutils.printBinContentAndError(TFhist)}")

    c1.cd()
    c1.Update()
    # c1.RedrawAxis()
    ROOT.gPad.RedrawAxis()
    ROOT.gPad.RedrawAxis("G")
    c1.SaveAs(saveName+".png")
    c1.Close()
    del c1, leg

    return TFhist, summed_SR, summed_CR


def checks(DataList, ABCDHistoVar, SRcut, CRcut, SVJbins, plotOutputDir,year=2018, perSVJbin=False):
    Data, sgData, bgData = DataList

    for d in bgData:
        integralConsistency = TF.checkIntegralConsistency(d, ABCDHistoVar, CRcut, SVJbins)
        # print(f"Integral consistency  - {integralConsistency}")
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
                print(f'uning ')
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
    outer_edge_results = {}  

    # Separate lists for each SVJ category and each component (Data, Signal, Background)
    obs_0SVJ_data, obs_0SVJ_sig, obs_0SVJ_bg = [], [], []
    pred_0SVJ_data, pred_0SVJ_sig, pred_0SVJ_bg = [], [], []

    obs_1SVJ_data, obs_1SVJ_sig, obs_1SVJ_bg = [], [], []
    pred_1SVJ_data, pred_1SVJ_sig, pred_1SVJ_bg = [], [], []

    obs_2SVJ_data, obs_2SVJ_sig, obs_2SVJ_bg = [], [], []
    pred_2SVJ_data, pred_2SVJ_sig, pred_2SVJ_bg = [], [], []

    obs_2PSVJ_data, obs_2PSVJ_sig, obs_2PSVJ_bg = [], [], []
    pred_2PSVJ_data, pred_2PSVJ_sig, pred_2PSVJ_bg = [], [], []

    obs_3PSVJ_data, obs_3PSVJ_sig, obs_3PSVJ_bg = [], [], []
    pred_3PSVJ_data, pred_3PSVJ_sig, pred_3PSVJ_bg = [], [], []

    obserr_0SVJ_data, obserr_0SVJ_sig, obserr_0SVJ_bg = [], [], []
    prederr_0SVJ_data, prederr_0SVJ_sig, prederr_0SVJ_bg = [], [], []

    obserr_1SVJ_data, obserr_1SVJ_sig, obserr_1SVJ_bg = [], [], []
    prederr_1SVJ_data, prederr_1SVJ_sig, prederr_1SVJ_bg = [], [], []

    obserr_2SVJ_data, obserr_2SVJ_sig, obserr_2SVJ_bg = [], [], []
    prederr_2SVJ_data, prederr_2SVJ_sig, prederr_2SVJ_bg = [], [], []

    obserr_2PSVJ_data, obserr_2PSVJ_sig, obserr_2PSVJ_bg = [], [], []
    prederr_2PSVJ_data, prederr_2PSVJ_sig, prederr_2PSVJ_bg = [], [], []

    obserr_3PSVJ_data, obserr_3PSVJ_sig, obserr_3PSVJ_bg = [], [], []
    prederr_3PSVJ_data, prederr_3PSVJ_sig, prederr_3PSVJ_bg = [], [], []
    # **Nested loop over DNN and MET outer edges**
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
            elif "2SVJ" in svj_bin:
                obs_data, obs_sig, obs_bg = obs_2SVJ_data, obs_2SVJ_sig, obs_2SVJ_bg
                pred_data, pred_sig, pred_bg = pred_2SVJ_data, pred_2SVJ_sig, pred_2SVJ_bg
                obs_dataerr, obs_sigerr, obs_bgerr = obserr_2SVJ_data, obserr_2SVJ_sig, obserr_2SVJ_bg
                pred_dataerr, pred_sigerr, pred_bgerr = prederr_2SVJ_data, prederr_2SVJ_sig, prederr_2SVJ_bg
            elif "2PSVJ" in svj_bin:
                obs_data, obs_sig, obs_bg = obs_2PSVJ_data, obs_2PSVJ_sig, obs_2PSVJ_bg
                pred_data, pred_sig, pred_bg = pred_2PSVJ_data, pred_2PSVJ_sig, pred_2PSVJ_bg
                obs_dataerr, obs_sigerr, obs_bgerr = obserr_2PSVJ_data, obserr_2PSVJ_sig, obserr_2PSVJ_bg
                pred_dataerr, pred_sigerr, pred_bgerr = prederr_2PSVJ_data, prederr_2PSVJ_sig, prederr_2PSVJ_bg  
            elif "3PSVJ" in svj_bin:
                obs_data, obs_sig, obs_bg = obs_3PSVJ_data, obs_3PSVJ_sig, obs_3PSVJ_bg
                pred_data, pred_sig, pred_bg = pred_3PSVJ_data, pred_3PSVJ_sig, pred_3PSVJ_bg
                obs_dataerr, obs_sigerr, obs_bgerr = obserr_3PSVJ_data, obserr_3PSVJ_sig, obserr_3PSVJ_bg
                pred_dataerr, pred_sigerr, pred_bgerr = prederr_3PSVJ_data, prederr_3PSVJ_sig, prederr_3PSVJ_bg
            else:
                print(f"Warning: No valid SVJ label found for {label}")
                continue  

            # Compute predicted A
            N_A_pred = (N_B * N_C) / N_D if N_D > 0 else float('nan')
            relative_error = np.sqrt(
                (NB_err / N_B if N_B > 0 else 0) ** 2 +
                (NC_err / N_C if N_C > 0 else 0) ** 2 +
                (ND_err / N_D if N_D > 0 else 0) ** 2
            )
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

            #print(f"[{label}] DNN_outer={DNN_outer_edge}, MET_outer={MET_outer_edge}, SVJ Bin: {svj_bin} | Observed A: {N_A}, Observed Err A: {NA_err}, Predicted A: {N_A_pred}, Predicted Err: {NA_pred_err}")
            print(f"[{label}] outer_edge: {dnn,met}, SVJ Bin: {svj_bin} | Observed A: {N_A},Observed ERR A: {NA_err},Predicted A: {N_A_pred}, Predicted Err {NA_pred_err}")
    return (
        outer_edge_results,
        obs_0SVJ_data, pred_0SVJ_data, obs_0SVJ_sig, pred_0SVJ_sig, obs_0SVJ_bg, pred_0SVJ_bg,
        obs_1SVJ_data, pred_1SVJ_data, obs_1SVJ_sig, pred_1SVJ_sig, obs_1SVJ_bg, pred_1SVJ_bg,
        obs_2SVJ_data, pred_2SVJ_data, obs_2SVJ_sig, pred_2SVJ_sig, obs_2SVJ_bg, pred_2SVJ_bg,
        obs_2PSVJ_data, pred_2PSVJ_data, obs_2PSVJ_sig, pred_2PSVJ_sig, obs_2PSVJ_bg, pred_2PSVJ_bg,
        obs_3PSVJ_data, pred_3PSVJ_data, obs_3PSVJ_sig, pred_3PSVJ_sig, obs_3PSVJ_bg, pred_3PSVJ_bg,
        obserr_0SVJ_data, obserr_0SVJ_sig, obserr_0SVJ_bg, prederr_0SVJ_data, prederr_0SVJ_sig, prederr_0SVJ_bg,
        obserr_1SVJ_data, obserr_1SVJ_sig, obserr_1SVJ_bg, prederr_1SVJ_data, prederr_1SVJ_sig, prederr_1SVJ_bg,
        obserr_2SVJ_data, obserr_2SVJ_sig, obserr_2SVJ_bg, prederr_2SVJ_data, prederr_2SVJ_sig, prederr_2SVJ_bg,
        obserr_2PSVJ_data, obserr_2PSVJ_sig, obserr_2PSVJ_bg, prederr_2PSVJ_data, prederr_2PSVJ_sig, prederr_2PSVJ_bg,
        obserr_3PSVJ_data, obserr_3PSVJ_sig, obserr_3PSVJ_bg, prederr_3PSVJ_data, prederr_3PSVJ_sig, prederr_3PSVJ_bg,
    )

def compute_nonclosure(obs, pred,):
        obs = np.array(obs)
        pred = np.array(pred)
        #obs_err = np.array(obs_err)
        #pred_err = np.array(pred_err)

        ratio = np.where(pred > 0, obs / pred, float('nan'))
        nonclosure = np.ones_like(ratio) - np.reciprocal(ratio)

        #ratio_err = np.where(pred > 0, ratio * np.sqrt((obs_err / obs) ** 2 + (pred_err / pred) ** 2), float('nan'))
        #nonclosure_err = (1 / ratio ** 2) * ratio_err
        return nonclosure

def plot_nonclosure_histogram(nonclosure, x_edges, y_edges,max_nonclosure, title, filename):
            plt.figure(figsize=(10, 8))
            X, Y = np.meshgrid(x_edges, y_edges)
            
            # Transpose nonclosure if needed to match meshgrid orientation
            Z = np.abs(nonclosure).reshape(len(x_edges),len(y_edges))
            #print(f'RESHAPED nonclosure {Z}')
            cmap = plt.cm.viridis
            pcm = plt.pcolormesh(X, Y, Z, cmap=cmap, vmin=0, vmax=max_nonclosure, shading='auto')
            plt.colorbar(pcm, label='|Non-Closure|')
            hep.cms.label(rlabel="")
            plt.xlabel('MET [GeV]')
            plt.ylabel('DNN Score')
            plt.title(title,loc='right')
            plt.savefig(filename)
            plt.close()

def plot_ABCD_ratios(
        year,
        outer_edge_results, 
        obs_0SVJ_data, pred_0SVJ_data, obs_0SVJ_sig, pred_0SVJ_sig, obs_0SVJ_bg, pred_0SVJ_bg,
        obs_1SVJ_data, pred_1SVJ_data, obs_1SVJ_sig, pred_1SVJ_sig, obs_1SVJ_bg, pred_1SVJ_bg,
        obs_2SVJ_data, pred_2SVJ_data, obs_2SVJ_sig, pred_2SVJ_sig, obs_2SVJ_bg, pred_2SVJ_bg,
        obs_2PSVJ_data, pred_2PSVJ_data, obs_2PSVJ_sig, pred_2PSVJ_sig, obs_2PSVJ_bg, pred_2PSVJ_bg,
        obs_3PSVJ_data, pred_3PSVJ_data, obs_3PSVJ_sig, pred_3PSVJ_sig, obs_3PSVJ_bg, pred_3PSVJ_bg,
        obserr_0SVJ_data, obserr_0SVJ_sig, obserr_0SVJ_bg,prederr_0SVJ_data, prederr_0SVJ_sig, prederr_0SVJ_bg,
        obserr_1SVJ_data, obserr_1SVJ_sig, obserr_1SVJ_bg,prederr_1SVJ_data, prederr_1SVJ_sig, prederr_1SVJ_bg,
        obserr_2SVJ_data, obserr_2SVJ_sig, obserr_2SVJ_bg,prederr_2SVJ_data, prederr_2SVJ_sig, prederr_2SVJ_bg,
        obserr_2PSVJ_data, obserr_2PSVJ_sig, obserr_2PSVJ_bg,prederr_2PSVJ_data, prederr_2PSVJ_sig, prederr_2PSVJ_bg,
        obserr_3PSVJ_data, obserr_3PSVJ_sig, obserr_3PSVJ_bg, prederr_3PSVJ_data, prederr_3PSVJ_sig, prederr_3PSVJ_bg,
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

    ratio_2SVJ_data, errbars_2SVJ_data = compute_ratio(obs_2SVJ_data,obserr_2SVJ_data, pred_2SVJ_data,prederr_2SVJ_data)
    ratio_2SVJ_bg, errbars_2SVJ_bg = compute_ratio(obs_2SVJ_bg,obserr_2SVJ_bg, pred_2SVJ_bg,prederr_2SVJ_bg)
    ratio_2SVJ_sig, errbars_2SVJ_sig = compute_ratio(obs_2SVJ_sig,obserr_2SVJ_sig,prederr_2SVJ_sig ,pred_2SVJ_sig)

    ratio_2PSVJ_data, errbars_2PSVJ_data = compute_ratio(obs_2PSVJ_data,obserr_2PSVJ_data, pred_2PSVJ_data,prederr_2PSVJ_data)
    ratio_2PSVJ_bg, errbars_2PSVJ_bg = compute_ratio(obs_2PSVJ_bg,obserr_2PSVJ_bg, pred_2PSVJ_bg,prederr_2PSVJ_bg)
    ratio_2PSVJ_sig, errbars_2PSVJ_sig = compute_ratio(obs_2PSVJ_sig,obserr_2PSVJ_sig,prederr_2PSVJ_sig ,pred_2PSVJ_sig)

    ratio_3PSVJ_data, errbars_3PSVJ_data = compute_ratio(obs_3PSVJ_data,obserr_3PSVJ_data, pred_3PSVJ_data,prederr_3PSVJ_data)
    ratio_3PSVJ_bg, errbars_3PSVJ_bg = compute_ratio(obs_3PSVJ_bg,obserr_3PSVJ_bg, pred_3PSVJ_bg,prederr_3PSVJ_bg)
    ratio_3PSVJ_sig, errbars_3PSVJ_sig = compute_ratio(obs_3PSVJ_sig,obserr_3PSVJ_sig,prederr_3PSVJ_sig ,pred_3PSVJ_sig)

    # Compute Data Ratio - Background Ratio
    diff_0SVJ_data_bg = ratio_0SVJ_data - ratio_0SVJ_bg
    diff_1SVJ_data_bg = ratio_1SVJ_data - ratio_1SVJ_bg
    diff_2SVJ_data_bg = ratio_2SVJ_data - ratio_2SVJ_bg
    diff_2PSVJ_data_bg = ratio_2PSVJ_data - ratio_2PSVJ_bg
    diff_3PSVJ_data_bg = ratio_3PSVJ_data - ratio_3PSVJ_bg

    errbars_0SVJ_diff = np.sqrt(errbars_0SVJ_data**2 + errbars_0SVJ_bg**2)
    errbars_1SVJ_diff = np.sqrt(errbars_1SVJ_data**2 + errbars_1SVJ_bg**2)
    errbars_2SVJ_diff = np.sqrt(errbars_2SVJ_data**2 + errbars_2SVJ_bg**2)
    errbars_2PSVJ_diff = np.sqrt(errbars_2PSVJ_data**2 + errbars_2PSVJ_bg**2)
    errbars_3PSVJ_diff = np.sqrt(errbars_3PSVJ_data**2 + errbars_3PSVJ_bg**2)

    #compute Signal contamination 
    sig_contam_0SVJ, sig_contam_err_0SVJ = signalcontamination(obs_0SVJ_sig, obserr_0SVJ_sig, obs_0SVJ_bg, obserr_0SVJ_bg)
    sig_contam_1SVJ, sig_contam_err_1SVJ = signalcontamination(obs_1SVJ_sig, obserr_1SVJ_sig, obs_1SVJ_bg, obserr_1SVJ_bg)
    sig_contam_2SVJ, sig_contam_err_2SVJ = signalcontamination(obs_2SVJ_sig, obserr_2SVJ_sig, obs_2SVJ_bg, obserr_2SVJ_bg)
    sig_contam_2PSVJ, sig_contam_err_2PSVJ = signalcontamination(obs_2PSVJ_sig, obserr_2PSVJ_sig, obs_2PSVJ_bg, obserr_2PSVJ_bg)
    sig_contam_3PSVJ, sig_contam_err_3PSVJ = signalcontamination(obs_3PSVJ_sig, obserr_3PSVJ_sig, obs_3PSVJ_bg, obserr_3PSVJ_bg)

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

        txt_filename = filename.replace('.pdf', '.txt')
        txt_path = os.path.join(output_dir, txt_filename)
        with open(txt_path, 'w') as f:
            f.write("Boundary\t" + "\t".join([f"{label}_NonClosure\t{label}_NC_err" for label in labels]) + "\tDiff(Data-MC)\tDiff_er(Data-MC)r\n")
            for i in range(len(x_values)):
                row = [f"{x_values[i]:.5f}"]
                for nc, err in zip(non_closures, nonclosure_err):
                    row.append(f"{nc[i]:.5f}\t{err[i]:.5f}")
                row.append(f"{diff_ratios[0][i]:.5f}\t{diff_error[0][i]:.5f}")
                f.write("\t".join(row) + "\n")
                
    plot_and_save(values, 
                [ratio_0SVJ_data, ratio_0SVJ_bg],
                [errbars_0SVJ_data, errbars_0SVJ_bg],
                [diff_0SVJ_data_bg],[errbars_0SVJ_diff], 
                ["Data 0SVJ", "Background MC 0SVJ"], 
                ["b", "r"],
                title, 
                "Ratio_Data_0SVJ.pdf")

    plot_and_save(values, 
                [ratio_1SVJ_data, ratio_1SVJ_bg],
                [errbars_1SVJ_data, errbars_1SVJ_bg],
                [diff_1SVJ_data_bg],[errbars_1SVJ_diff], 
                ["Data 1SVJ", "Background MC 1SVJ"], 
                ["b", "r"], 
                title,
                "Ratio_Data_1SVJ.pdf")
    print(f"data {ratio_1SVJ_data} err bars {errbars_1SVJ_data}")

    plot_and_save(values, 
                [ratio_2SVJ_data, ratio_2SVJ_bg],
                [errbars_2SVJ_data, errbars_2SVJ_bg],
                [diff_2SVJ_data_bg],[errbars_2SVJ_diff], 
                ["Data 2SVJ", "Background MC 2SVJ"], 
                ["b", "r"], 
                title,
                "Ratio_Data_2SVJ.pdf")

    plot_and_save(values, 
                [ratio_2PSVJ_data, ratio_2PSVJ_bg],
                [errbars_2PSVJ_data, errbars_2PSVJ_bg],
                [diff_2PSVJ_data_bg],[errbars_2PSVJ_diff], 
                ["Data 2PSVJ", "Background MC 2PSVJ"], 
                ["b", "r"], 
                title,
                "Ratio_Data_2PSVJ.pdf")

    plot_and_save(values, 
                [ratio_3PSVJ_data, ratio_3PSVJ_bg],
                [errbars_3PSVJ_data, errbars_3PSVJ_bg],
                [diff_3PSVJ_data_bg],[errbars_3PSVJ_diff], 
                ["Data 3PSVJ", "Background MC 3PSVJ"], 
                ["b", "r"], 
                title,
                "Ratio_Data_3PSVJ.pdf")
    
    #print(f"non closure data 2psvj {ratio_2SVJ_data} ")
    # plot_signal_contamination(values, sig_contam_0SVJ, sig_contam_err_0SVJ, "0SVJ", "signal_contamination_0SVJ.jpg")
    # plot_signal_contamination(values, sig_contam_1SVJ, sig_contam_err_1SVJ, "1SVJ", "signal_contamination_1SVJ.jpg")
    # plot_signal_contamination(values, sig_contam_2SVJ, sig_contam_err_2SVJ, "2SVJ", "signal_contamination_2SVJ.jpg")
    # plot_signal_contamination(values, sig_contam_3PSVJ, sig_contam_err_3PSVJ, "3PSVJ", "signal_contamination_3PSVJ.jpg")
    

    #print(f'Boundary Values {values[::-1]}')
    print(f'Boundary Values {values}')

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

def plot_signal_contamination(x_values, sig_contam, sig_contam_err, label,year,output_dir, filename):
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

#!/usr/bin/env python3
"""
Fixed main() for NDNN_Nonclosure-style script:
 - Single-year (no loop)
 - Computes PNET (_pre_) and WNAE (_pre_WNAE_) ABCD predictions
 - Computes background non-closure for 0/1/2/3 SVJ
 - Saves 2D non-closure heatmaps via plot_nonclosure_histogram
"""

import os
import sys
import optparse
import traceback

import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep

# Use CMS style for any plotting done inside helper functions
hep.style.use("CMS")


def main():
    parser = optparse.OptionParser("usage: %prog [options]\n")
    parser.add_option('-d', '--dataset', dest='dataset',
                      default='testHadd_11242020', help='dataset')
    parser.add_option('-y', dest='year', type='string', default='2018',
                      help="Run year (e.g. 2018)")
    parser.add_option('-o', '--outputdir', dest='outputdir', type='string',
                      default='Hyperparamtrained_NC',
                      help="Top-level output folder name (e.g. Hyperparamtrained_NC)")
    parser.add_option('--hemPeriod', dest='hemPeriod', type=str, default=False,
                      help='HEM period (PreHEM or PostHEM), default includes entire sample')

    options, args = parser.parse_args()
    print("Parsed options:", options)

    base_output = options.outputdir
    year = options.year
    hemPeriod = options.hemPeriod

    # --- Binning setup (adjust bin counts/edges as you need) ---
    met_bins = 50
    dnn_bins = 50

    # Inner/outer edges used by compute_ABCD_prediction (keep same orientation as your original)
    MET_inner_edges = np.linspace(205, 1000, met_bins)
    MET_outer_edges = np.linspace(20000, 20000, met_bins)  # placeholder 'outer' value as before

    DNN_inner_edges = np.linspace(0.1, 0.9, dnn_bins)
    DNN_outer_edges = np.linspace(1, 1, dnn_bins)  # placeholder 'outer' value as before

    # --- Load data for provided year ---
    try:
        Data, sgData, bgData = getData(options.dataset + "/", 1.0, year)
    except Exception as e:
        print("[ERROR] getData failed:", e)
        traceback.print_exc()
        sys.exit(1)

    # ------------------------
    # PNET (_pre_) processing
    # ------------------------
    output_dir_pnet = f'{base_output}/PNET/2D-NonC/{year}/ControlRegion'
    os.makedirs(output_dir_pnet, exist_ok=True)
    print(f"[INFO] Output directory (PNET): {output_dir_pnet}")

    try:
        (
            outer_edge_results,
            obs_0SVJ_data, pred_0SVJ_data, obs_0SVJ_sig, pred_0SVJ_sig, obs_0SVJ_bg, pred_0SVJ_bg,
            obs_1SVJ_data, pred_1SVJ_data, obs_1SVJ_sig, pred_1SVJ_sig, obs_1SVJ_bg, pred_1SVJ_bg,
            obs_2SVJ_data, pred_2SVJ_data, obs_2SVJ_sig, pred_2SVJ_sig, obs_2SVJ_bg, pred_2SVJ_bg,
            obs_3SVJ_data, pred_3SVJ_data, obs_3SVJ_sig, pred_3SVJ_sig, obs_3SVJ_bg, pred_3SVJ_bg,
            obserr_0SVJ_data, obserr_0SVJ_sig, obserr_0SVJ_bg, prederr_0SVJ_data, prederr_0SVJ_sig, prederr_0SVJ_bg,
            obserr_1SVJ_data, obserr_1SVJ_sig, obserr_1SVJ_bg, prederr_1SVJ_data, prederr_1SVJ_sig, prederr_1SVJ_bg,
            obserr_2SVJ_data, obserr_2SVJ_sig, obserr_2SVJ_bg, prederr_2SVJ_data, prederr_2SVJ_sig, prederr_2SVJ_bg,
            obserr_3SVJ_data, obserr_3SVJ_sig, obserr_3SVJ_bg, prederr_3SVJ_data, prederr_3SVJ_sig, prederr_3SVJ_bg
        ) = compute_ABCD_prediction(
            Data, sgData, bgData,
            "h_METvsDNN",          # histogram name
            '_pre_',               # cut tag for PNET
            DNN_inner_edges,
            DNN_outer_edges,
            MET_inner_edges,
            MET_outer_edges
        )
    except Exception as e:
        print("[ERROR] compute_ABCD_prediction (PNET) failed:", e)
        traceback.print_exc()
        sys.exit(1)

    # Compute background non-closure arrays (expecting shapes suitable for plotting)
    try:
        nonclosure_0SVJ = compute_nonclosure(obs_0SVJ_bg, pred_0SVJ_bg)
        nonclosure_1SVJ = compute_nonclosure(obs_1SVJ_bg, pred_1SVJ_bg)
        nonclosure_2SVJ = compute_nonclosure(obs_2SVJ_bg, pred_2SVJ_bg)
        nonclosure_3SVJ = compute_nonclosure(obs_3SVJ_bg, pred_3SVJ_bg)
    except Exception as e:
        print("[ERROR] compute_nonclosure (PNET) failed:", e)
        traceback.print_exc()
        sys.exit(1)

    print(f"[INFO] PNET: non-closure shapes: 0SVJ={getattr(nonclosure_0SVJ,'shape', 'unk')}, 1SVJ={getattr(nonclosure_1SVJ,'shape','unk')}")

    # Plot and save (MET edges X, DNN edges Y). Adjust edge arrays to match plotting function expectation:
    MET_edges_plot = np.linspace(200, 1000, met_bins)
    DNN_edges_plot = np.linspace(0.0, 1.0, dnn_bins)

    try:
        plot_nonclosure_histogram(nonclosure_0SVJ, MET_edges_plot, DNN_edges_plot, 1, '0SVJ PNET', f"{output_dir_pnet}/nonclosure_0SVJ_PNET.png")
        plot_nonclosure_histogram(nonclosure_1SVJ, MET_edges_plot, DNN_edges_plot, 1, '1SVJ PNET', f"{output_dir_pnet}/nonclosure_1SVJ_PNET.png")
        plot_nonclosure_histogram(nonclosure_2SVJ, MET_edges_plot, DNN_edges_plot, 1, '2SVJ PNET', f"{output_dir_pnet}/nonclosure_2SVJ_PNET.png")
        plot_nonclosure_histogram(nonclosure_3SVJ, MET_edges_plot, DNN_edges_plot, 1, '3SVJ PNET', f"{output_dir_pnet}/nonclosure_3SVJ_PNET.png")
    except Exception as e:
        print("[ERROR] plot_nonclosure_histogram (PNET) failed:", e)
        traceback.print_exc()
        # continue to WNAE section even if plotting PNET failed

    # ------------------------
    # WNAE (_pre_WNAE_) processing
    # ------------------------
    output_dir_wnae = f'{base_output}/WNAE/2D-NonC/{year}/ControlRegion'
    os.makedirs(output_dir_wnae, exist_ok=True)
    print(f"[INFO] Output directory (WNAE): {output_dir_wnae}")

    try:
        (
            outer_edge_results_wnae,
            obs_0SVJ_data, pred_0SVJ_data, obs_0SVJ_sig, pred_0SVJ_sig, obs_0SVJ_bg, pred_0SVJ_bg,
            obs_1SVJ_data, pred_1SVJ_data, obs_1SVJ_sig, pred_1SVJ_sig, obs_1SVJ_bg, pred_1SVJ_bg,
            obs_2SVJ_data, pred_2SVJ_data, obs_2SVJ_sig, pred_2SVJ_sig, obs_2SVJ_bg, pred_2SVJ_bg,
            obs_3SVJ_data, pred_3SVJ_data, obs_3SVJ_sig, pred_3SVJ_sig, obs_3SVJ_bg, pred_3SVJ_bg,
            obserr_0SVJ_data, obserr_0SVJ_sig, obserr_0SVJ_bg, prederr_0SVJ_data, prederr_0SVJ_sig, prederr_0SVJ_bg,
            obserr_1SVJ_data, obserr_1SVJ_sig, obserr_1SVJ_bg, prederr_1SVJ_data, prederr_1SVJ_sig, prederr_1SVJ_bg,
            obserr_2SVJ_data, obserr_2SVJ_sig, obserr_2SVJ_bg, prederr_2SVJ_data, prederr_2SVJ_sig, prederr_2SVJ_bg,
            obserr_3SVJ_data, obserr_3SVJ_sig, obserr_3SVJ_bg, prederr_3SVJ_data, prederr_3SVJ_sig, prederr_3SVJ_bg
        ) = compute_ABCD_prediction(
            Data, sgData, bgData,
            "h_METvsDNN",          # histogram name
            '_pre_WNAE_',          # cut tag for WNAE
            DNN_inner_edges,
            DNN_outer_edges,
            MET_inner_edges,
            MET_outer_edges
        )
    except Exception as e:
        print("[ERROR] compute_ABCD_prediction (WNAE) failed:", e)
        traceback.print_exc()
        sys.exit(1)

    # Compute background non-closure arrays for WNAE
    try:
        nonclosure_0SVJ_w = compute_nonclosure(obs_0SVJ_bg, pred_0SVJ_bg)
        nonclosure_1SVJ_w = compute_nonclosure(obs_1SVJ_bg, pred_1SVJ_bg)
        nonclosure_2SVJ_w = compute_nonclosure(obs_2SVJ_bg, pred_2SVJ_bg)
        nonclosure_3SVJ_w = compute_nonclosure(obs_3SVJ_bg, pred_3SVJ_bg)
    except Exception as e:
        print("[ERROR] compute_nonclosure (WNAE) failed:", e)
        traceback.print_exc()
        sys.exit(1)

    print(f"[INFO] WNAE: non-closure shapes: 0SVJ={getattr(nonclosure_0SVJ_w,'shape','unk')}")

    # Plot and save WNAE non-closure
    try:
        plot_nonclosure_histogram(nonclosure_0SVJ_w, MET_edges_plot, DNN_edges_plot, 1, '0SVJ WNAE', f"{output_dir_wnae}/nonclosure_0SVJ_WNAE.png")
        plot_nonclosure_histogram(nonclosure_1SVJ_w, MET_edges_plot, DNN_edges_plot, 1, '1SVJ WNAE', f"{output_dir_wnae}/nonclosure_1SVJ_WNAE.png")
        plot_nonclosure_histogram(nonclosure_2SVJ_w, MET_edges_plot, DNN_edges_plot, 1, '2SVJ WNAE', f"{output_dir_wnae}/nonclosure_2SVJ_WNAE.png")
        plot_nonclosure_histogram(nonclosure_3SVJ_w, MET_edges_plot, DNN_edges_plot, 1, '3SVJ WNAE', f"{output_dir_wnae}/nonclosure_3SVJ_WNAE.png")
    except Exception as e:
        print("[ERROR] plot_nonclosure_histogram (WNAE) failed:", e)
        traceback.print_exc()

    print("[INFO] Completed PNET and WNAE non-closure computation and plotting.")


if __name__ == '__main__':
    main()

