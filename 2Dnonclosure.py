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
# ROOT.TH1.SetDefaultSumw2()
# ROOT.TH2.SetDefaultSumw2()

# def GetSVJbins(outer_edge):
#     '''Change the SVJbin edges here these are the ones for DNN trained on only QCD '''    
#     SVJbins = {
#                 "0SVJ" : [0.5,350.0/outer_edge],
#                 "1SVJ" : [0.56,450.0/outer_edge],
#                 "2PSVJ" : [0.56,350.0/outer_edge],
#     }
#    return SVJbins


def GetSVJbins(DNN_inner_edge,MET_inner_edge):
    '''Change the SVJbin edges here these are the ones for DNN trained on all back '''    
    SVJbins = {
                "0SVJ" : [DNN_inner_edge,MET_inner_edge],
                "1SVJ" : [DNN_inner_edge,MET_inner_edge],
                "2PSVJ" : [DNN_inner_edge,MET_inner_edge],
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

def plotABCD(dataList, ABCDHistoVar, maincut, SVJBins, plotOutputDir, xTitle="nSVJ", yTitle="Events", xmin=999.9, xmax=-999.9, isLogY=False, year="2018", isRatio=False, hemPeriod=False):
    '''
    Plots either a stacked histogram for background with or without the signal (when data is empty)
    or a ratio plot between tData and background histograms. 
    '''
    ROOT.TH1.AddDirectory(False)
    firstpass = True
    Data, sgData, bgData = dataList
    if maincut == "_pre_":
        Data = None
    if isRatio and Data is None: 
        print("Not passed data cannot make ratio plot, making normal plot instead") 
        isRatio = False

    # Plot 1: Histogram (Background + Signal + Data)
    c1 = ROOT.TCanvas("c1", "c1", 800, 800)
    c1.cd()
    pltutils.SetupGPad(logY=isLogY)
    ROOT.gStyle.SetOptStat("")

    # TLegend 
    leg = pltutils.SetupLegend(NColumns=2, textSize=0.024)

    if bgData:
        bgDataMergedABCDDict = TF.GetABCDhistDict(bgData, ABCDHistoVar, maincut, SVJBins, isStack=True, merge=True)
        bgStackedHist, bgSummedHist = pltutils.StackedHistogram(bgDataMergedABCDDict)
        for bgDataHist in bgDataMergedABCDDict.keys():
            leg.AddEntry(bgDataMergedABCDDict[bgDataHist], bgDataHist.label_, "F")
        if firstpass:
            firstpass = False
            dummy = ROOT.TH1D("dummy", "dummy", len(SVJBins)*4, bgSummedHist.GetBinLowEdge(1), bgSummedHist.GetBinLowEdge(bgSummedHist.GetNbinsX()) + bgSummedHist.GetBinWidth(bgSummedHist.GetNbinsX()))
            ymax = 10**9
            ymin = 10
            lmax = 10**9
            for bin_idx in range(1, bgSummedHist.GetNbinsX() + 1):
                bin_label = bgSummedHist.GetXaxis().GetBinLabel(bin_idx)
                dummy.GetXaxis().SetBinLabel(bin_idx, bin_label)
            pltutils.setupDummy(dummy, leg, "", xTitle, yTitle, isLogY, xmin, xmax, ymin, ymax, lmax, isRatio=False)
            dummy.Draw("hist")
            
        bgStackedHist.Draw("hist F same")
        lines_upperpad = pltutils.AddVerticalLine(dummy, SVJBins, ymax=10**5)
        for line in lines_upperpad:
            line.Draw("same")
    
    if sgData: 
        sgDataMergedABCDDict = TF.GetABCDhistDict(sgData, ABCDHistoVar, maincut, SVJBins, merge=True)
        linestylenumber = 0
        linestyle = [ROOT.kSolid, ROOT.kDashed, ROOT.kDotted]
        for sgDataHist in sgDataMergedABCDDict.keys():
            leg.AddEntry(sgDataMergedABCDDict[sgDataHist], sgDataHist.label_, "L")
            sgDataMergedABCDDict[sgDataHist].SetLineStyle(linestyle[linestylenumber % 3])
            linestylenumber += 1
            sgDataMergedABCDDict[sgDataHist].SetLineWidth(3)
            sgDataMergedABCDDict[sgDataHist].Draw("hist same") 
            
    if Data is not None:
        for d in Data:
            DataHist = TF.GetABCDhist(d, ABCDHistoVar, maincut, SVJBins, merge=True)
            ROOT.gStyle.SetErrorX(0.)
            DataHist = pltutils.SetupDataStyle(DataHist)
            leg.AddEntry(DataHist, d.label_)
            DataHist.Draw("P same")
    
    leg.Draw("same")
    #pltutils.AddLabelsForABCD(dummy, SVJBins, yloc=5*10**4, labels=["dA", "dB", "dC", "dD"])
    dummy.Draw("AXIS same")
    
    if AddCMSText:
        pltutils.AddCMSLumiText(c1, year, isExtraText=True, hemPeriod=hemPeriod)
    
    c1.cd()
    c1.Update()
    ROOT.gPad.RedrawAxis()
    SaveName = plotOutputDir + "/ABCDPlot_" + maincut + "_Hist"
    c1.SaveAs(SaveName + ".png")
    c1.Close()
    del c1, leg

    # Plot 2: Data/MC Ratio
    if isRatio and Data is not None:
        c2 = ROOT.TCanvas("c2", "c2", 800, 700)
        c2, pad1, pad2 = pltutils.createCanvasPads(c2, isLogY)
        pad1.cd()
        ROOT.gStyle.SetOptStat("")
        leg = pltutils.SetupLegend(NColumns=2)

        # Make dummy to setup axis
        dummy = ROOT.TH1D("dummy", "dummy", len(SVJBins)*4, bgSummedHist.GetBinLowEdge(1), bgSummedHist.GetBinLowEdge(bgSummedHist.GetNbinsX()) + bgSummedHist.GetBinWidth(bgSummedHist.GetNbinsX()))
        ymax = 10**9
        ymin = 10
        lmax = 10**9
        pltutils.setupDummy(dummy, leg, "", xTitle, yTitle, isLogY, xmin, xmax, ymin, ymax, lmax, isRatio=True)
        dummy.Draw("hist")

        # Draw the ratio plot
        pad2.cd()
        ratioHist = pltutils.RatioHistogram(DataHist, bgSummedHist)
        ratioHist = pltutils.SetupRatioStyle(ratioHist, xTitle, yTitle="Data/MC", yTitleSize=0.13)
        ratioHist.Draw("EX0P")
        lines_lowerpad = pltutils.AddVerticalLine(DataHist, SVJBins, ymax=2)
        for line in lines_lowerpad:
            line.Draw("same")
        
        if AddCMSText:
            pltutils.AddCMSLumiText(c2, year, isExtraText=True, hemPeriod=hemPeriod)
        
        c2.cd()
        c2.Update()
        ROOT.gPad.RedrawAxis()
        SaveName = plotOutputDir + "/ABCDPlot_" + maincut + "_Ratio"
        c2.SaveAs(SaveName + ".png")
        c2.Close()
        del c2, leg

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
    #outer_edges = [1.1,1.2,1.3,1.5,1.6,1.7,1.8,1.9,2]
    #NA_Predicted = [] 
    #for outer_edge in outer_edges:
    #    regions = GetSubABCDregions(met,dnn,outer_edge=2)    
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

def compute_ABCD_prediction(Data, sgData, bgData, ABCDHistoVar, maincut, DNN_inner_edges, DNN_outer_edges, MET_inner_edges, MET_outer_edges):
    outer_edge_results = {}  

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

    # **Nested loop over DNN and MET outer edges**
    for (DNN_inner_edge, DNN_outer_edge), (MET_inner_edge, MET_outer_edge) in product(zip(DNN_inner_edges, DNN_outer_edges),zip(MET_inner_edges, MET_outer_edges)):
        print(f"\nProcessing outer_edge: DNN={DNN_outer_edge}, MET={MET_outer_edge}")

        SVJBins = GetSVJbins(DNN_inner_edge, MET_inner_edge)  # Ensure correct binning
        
        # **FIX: Pass both MET_inner_edge and MET_outer_edge**
        counts_by_SVJ, error_by_SVJ = get_ABCD_counts_by_SVJ(
            Data, sgData, bgData, ABCDHistoVar, maincut, SVJBins, 
            DNN_inner_edge, DNN_outer_edge, MET_inner_edge, MET_outer_edge
        )
            
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

    return (
        outer_edge_results,
        obs_0SVJ_data, pred_0SVJ_data, obs_0SVJ_sig, pred_0SVJ_sig, obs_0SVJ_bg, pred_0SVJ_bg,
        obs_1SVJ_data, pred_1SVJ_data, obs_1SVJ_sig, pred_1SVJ_sig, obs_1SVJ_bg, pred_1SVJ_bg,
        obs_2PSVJ_data, pred_2PSVJ_data, obs_2PSVJ_sig, pred_2PSVJ_sig, obs_2PSVJ_bg, pred_2PSVJ_bg,
        obserr_0SVJ_data, obserr_0SVJ_sig, obserr_0SVJ_bg, prederr_0SVJ_data, prederr_0SVJ_sig, prederr_0SVJ_bg,
        obserr_1SVJ_data, obserr_1SVJ_sig, obserr_1SVJ_bg, prederr_1SVJ_data, prederr_1SVJ_sig, prederr_1SVJ_bg,
        obserr_2PSVJ_data, obserr_2PSVJ_sig, obserr_2PSVJ_bg, prederr_2PSVJ_data, prederr_2PSVJ_sig, prederr_2PSVJ_bg
    )




def main():
    parser = optparse.OptionParser("usage: %prog [options]\n")
    parser.add_option('-d', '--dataset',    dest='dataset',                    default='testHadd_11242020',    help='dataset')
    parser.add_option('-y',                 dest='year',       type='string',  default='2018',                 help="Can pass in the run year")
    parser.add_option('-o',                 dest='outputdir',  type='string',                                  help="Output folder name")
    parser.add_option(    '--hemPeriod',  dest='hemPeriod', type=str, default=False,  help='HEM period (PreHEM or PostHEM), default includes entire sample',)
    options, args = parser.parse_args()
    print("Parsed options:", options)

    # outROOTfileName = options.outfile
    year = options.year
    hemPeriod = options.hemPeriod
    ABCDhistoVars = ["METvsDNN"]
    ABCDFolderName = "ABCD"
    SRCut = "_pre_"
    #CRCuts = ["_lcr_pre_"]#,"_cr_muon_","_cr_electron_"]
    
    #DataCut = "_pre_"                                                               
    CRCuts = ["_lcr_pre_"]        
    maincuts = CRCuts
    Data, sgData, bgData = getData( options.dataset + "/", 1.0, year)

    bins = 10
    MET_outer_edges = np.linspace(20000, 20000, bins)
    MET_inner_edges = np.linspace(205, 1000, bins)

    DNN_outer_edges = np.linspace(1,1,bins)
    DNN_inner_edges = np.linspace(0.1, 0.90,bins)
    
    Boundary_vals = []
    for MET_outer_edge in MET_outer_edges:
        for DNN_outer_edge in DNN_outer_edges:
            boundary_val = [DNN_outer_edge, MET_outer_edge]
            Boundary_vals.append(boundary_val)

    # Define output directory
    output_dir = f'Nonclosure/2D-NonC/{year}'
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Compute ABCD predictions
    outer_edge_results, obs_0SVJ_data, pred_0SVJ_data, obs_0SVJ_sig, pred_0SVJ_sig, obs_0SVJ_bg, pred_0SVJ_bg, \
    obs_1SVJ_data, pred_1SVJ_data, obs_1SVJ_sig, pred_1SVJ_sig, obs_1SVJ_bg, pred_1SVJ_bg, \
    obs_2PSVJ_data, pred_2PSVJ_data, obs_2PSVJ_sig, pred_2PSVJ_sig, obs_2PSVJ_bg, pred_2PSVJ_bg, \
    obserr_0SVJ_data, obserr_0SVJ_sig, obserr_0SVJ_bg, prederr_0SVJ_data, prederr_0SVJ_sig, prederr_0SVJ_bg, \
    obserr_1SVJ_data, obserr_1SVJ_sig, obserr_1SVJ_bg, prederr_1SVJ_data, prederr_1SVJ_sig, prederr_1SVJ_bg, \
    obserr_2PSVJ_data, oberr_2PSVJ_sig, obserr_2PSVJ_bg, prederr_2PSVJ_data, prederr_2PSVJ_sig, prederr_2PSVJ_bg = \
        compute_ABCD_prediction(Data, sgData, bgData, "h_METvsDNN", '_pre_',DNN_inner_edges, DNN_outer_edges, MET_inner_edges, MET_outer_edges)

    # Function to compute non-closure
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

    # Compute non-closure for each SVJ category
    nonclosure_0SVJ = compute_nonclosure(obs_0SVJ_bg, pred_0SVJ_bg)
    print(f'observed NA: {len(obs_0SVJ_bg)}')
    print(f'predicted NA: {len(pred_0SVJ_bg)}')
    print(f'non closure: {len(nonclosure_0SVJ)}')
    print(f'output dir {output_dir}')
    nonclosure_1SVJ = compute_nonclosure(obs_1SVJ_bg, pred_1SVJ_bg)
    nonclosure_2PSVJ = compute_nonclosure(obs_2PSVJ_bg, pred_2PSVJ_bg)
    # Function to plot 2D histograms for non-closure
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

    DNN_edges = np.linspace(0,1,bins)
    MET_edges = np.linspace(200,1000,bins)
    # Plot Non-Closure Histograms
    plot_nonclosure_histogram(nonclosure_0SVJ, MET_edges,DNN_edges,1, '0SVJ', f"{output_dir}/nonclosure_0SVJ_inneredgeonly.png")
    plot_nonclosure_histogram(nonclosure_1SVJ, MET_edges,DNN_edges,1, '1SVJ', f"{output_dir}/nonclosure_1SVJ_inneredgeonly.png")
    plot_nonclosure_histogram(nonclosure_2PSVJ, MET_edges,DNN_edges,1,'2PSVJ', f"{output_dir}/nonclosure_2PSVJ_inneredgeonly.png")

 
if __name__ == '__main__':
    main()


#ADD (0.85,250) as the WP 