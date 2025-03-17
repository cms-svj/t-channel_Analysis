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
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator
AddCMSText = True
# ROOT.TH1.SetDefaultSumw2()
# ROOT.TH2.SetDefaultSumw2()

def GetSVJbins(factor):
    '''Change the SVJbin edges here'''
    SVJbins = {
                "0SVJ" : [0.5/factor,350.0],
                "1SVJ" : [0.56/factor,450.0],
                "2PSVJ" : [0.56/factor,350.0],
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

def PlotValidation(TFhist, hist, expected_SR, MC_QCD_ZJets_SR, SVJBins,xTitle="nSVJ", yTitle="Events", xmin=999.9, xmax = -999.9, isLogY = False, year="2018", saveName="Validation", isDataCR=False,hemPeriod=False):
    '''
    validation is done by taking the ratio between the SR (Expected) and the Predicted SR.
    '''
    ROOT.TH1.AddDirectory(False)
    predicted_Data_SR = TF.Validation(hist,TFhist, MC_QCD_ZJets_SR)
    RatioHist = pltutils.RatioHistogram(expected_SR,predicted_Data_SR)
    c1 = ROOT.TCanvas( "c", "c", 800, 700)
    c1, pad1, pad2 = pltutils.createCanvasPads(c1,isLogY)
    pad1.cd()
    ROOT.gStyle.SetOptStat("")
    leg = pltutils.SetupLegend(x1=0.2, textSize=0.04)
    pltutils.SetupLineHistStyle(expected_SR)
    pltutils.SetupLineHistStyle(predicted_Data_SR, color=ROOT.kRed)
   
    # make dummy plot for axis
    dummy = ROOT.TH1D("dummy", "dummy", len(SVJBins)*4, expected_SR.GetBinLowEdge(1), expected_SR.GetBinLowEdge(expected_SR.GetNbinsX()) + expected_SR.GetBinWidth(expected_SR.GetNbinsX()))
    # print(f"The dummy values are  - {expected_SR.GetBinLowEdge(1)}, {expected_SR.GetBinLowEdge(expected_SR.GetNbinsX())}, { expected_SR.GetBinWidth(expected_SR.GetNbinsX())}")
    ymax=10**9
    ymin=10
    lmax=10**9
    pltutils.setupDummy(dummy,leg,"", xTitle, yTitle, isLogY, xmin, xmax, ymin, ymax, lmax, isRatio=True)
    dummy.Draw("hist")

    expected_SR.Draw("histe same")
    predicted_Data_SR.Draw("histe same")
    leg.AddEntry(expected_SR,"Expected SR (MC)","L")
    if isDataCR:
        leg.AddEntry(predicted_Data_SR,"Predicted SR ((Data_{CR} - MC^{QCD,ZJets}_{CR}) #times TF + MC^{QCD,ZJets}_{SR})","L")
    else:
        leg.AddEntry(predicted_Data_SR,"Predicted SR (CR #times TF + SR MC_{QCD,ZJets})","L")
    lines_upperpad = pltutils.AddVerticalLine(expected_SR, SVJBins, ymax = 10**5)
    for line in lines_upperpad:
        line.Draw("same")
    
    pltutils.AddLabelsForABCD(predicted_Data_SR,SVJBins,yloc=5*10**4)
    leg.Draw("same")
    pad2.cd()
    
    RatioHist = pltutils.SetupRatioStyle(RatioHist, xTitle="nSVJ", yTitle="Expected/Predicted", yTitleSize=0.1, ymin=0.6, ymax=1.4)
    RatioHist.Draw("EX0P")
    lines_lowerpad = pltutils.AddVerticalLine(expected_SR, SVJBins, ymin=0.5, ymax = 1.5)
    for line in lines_lowerpad:
        line.Draw("same")
    if AddCMSText:
        pltutils.AddCMSLumiText(c1, year,isExtraText=True,hemPeriod=hemPeriod)
    
    
    # print(f"Expected Hist -  {pltutils.printBinContentAndError(expected_SR)}")
    # print(f"Predicted Hist - {pltutils.printBinContentAndError(predicted_Data_SR)}")
    # print(f"Ratio Hist - {pltutils.printBinContentAndError(RatioHist)}")
    
    
    c1.cd()
    c1.Update()
    # c1.RedrawAxis()
    ROOT.gPad.RedrawAxis()
    ROOT.gPad.RedrawAxis("G")
    c1.SaveAs(saveName+".png")
    c1.Close()
    del c1, leg


def plotDataMCRatio(dataList, ABCDHistoVar, maincut, SVJBins, plotOutputDir, xTitle="nSVJ", year="2018"):
    """
    Creates and saves a separate Data/MC ratio plot.
    """
    ROOT.TH1.AddDirectory(False)
    Data, _, bgData = dataList


    print(f"+++++++++++++++++++++++Generating Data/MC ratio plot for {maincut}, saving to {plotOutputDir}/ABCD/")
    if not Data:
        print("No data provided, cannot make Data/MC ratio plot.")
        return
    
    c1 = ROOT.TCanvas("c_ratio", "c_ratio", 800, 700)
    pads = pltutils.createCanvasPads(c1, isLogY=False)

    c1, pad1, pad2 = pads  # Expecting three elements: canvas, upper pad, lower pad
    pad2.cd()


    pad2.cd()
    
    # Get background summed histogram
    bgDataMergedABCDDict = TF.GetABCDhistDict(bgData, ABCDHistoVar, maincut, SVJBins, isStack=True, merge=True)
    _, bgSummedHist = pltutils.StackedHistogram(bgDataMergedABCDDict)
    
    for d in Data:
        DataHist = TF.GetABCDhist(d, ABCDHistoVar, maincut, SVJBins, merge=True)
        ROOT.gStyle.SetErrorX(0.)
        ratioHist = pltutils.RatioHistogram(DataHist, bgSummedHist)
        ratioHist = pltutils.SetupRatioStyle(ratioHist, xTitle, yTitle="Data/MC", yTitleSize=0.13)
        ratioHist.Draw("EX0P")
        
        # Add vertical lines to indicate bin separation
        lines_lowerpad = pltutils.AddVerticalLine(DataHist, SVJBins, ymax=2)
        for line in lines_lowerpad:
            line.Draw("same")
    
    c1.cd()
    c1.Update()
    ROOT.gPad.RedrawAxis()
    SaveName = f"{plotOutputDir}/ABCD/Data_MC_Ratio_{maincut}.png"
    c1.SaveAs(SaveName)
    c1.Close()
    print(f"Saved Data/MC ratio plot: {SaveName}")



def checks(DataList, ABCDHistoVar, SRcut, CRcut, SVJbins, plotOutputDir,year=2018, perSVJbin=False):
    Data, sgData, bgData = DataList

    for d in bgData:
        integralConsistency = TF.checkIntegralConsistency(d, ABCDHistoVar, CRcut, SVJbins)
        # print(f"Integral consistency  - {integralConsistency}")
def GetSubABCDregions(met,dnn,factor):
    #Defining the A B C D regions
    print('MeT',met)
    print('dnn',dnn/factor)
    #print('dnn/3',dnn/2)
    regions = [ ("dA", met, 20000, dnn/factor, dnn),
                    ("dB", met, 20000, 0, dnn/factor),
                    ("dC", 0, met, dnn/factor,dnn),
                    ("dD", 0, met, 0, dnn/factor)
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

def GetABCDhistPerSVJBin(data, ABCDhistoVar, maincut, SVJbin,factor):
    """Return the dictionary of ABCD histogram for only one svj bin"""

    bkgname = data.fileName.split('_')[1].replace('.root','')
    hist_dict = {region: ROOT.TH1F(f"h_{bkgname}_{region}",f"h_{bkgname}_{region}", 1,0,1) for region in ['dA','dB','dC','dD']}
    SVJ, (dnn, met) = SVJbin
    histName = ABCDhistoVar + maincut + SVJ
    #factors = [1.1,1.2,1.3,1.5,1.6,1.7,1.8,1.9,2]
    #NA_Predicted = [] 
    #for factor in factors:
    #    regions = GetSubABCDregions(met,dnn,factor=2)    
    regions = GetSubABCDregions(met,dnn,factor)

    for region, xmin, xmax, ymin, ymax in regions:   # TODO: this for condition should be written into a function work in all the cases
        xmin, xmax, ymin, ymax = adjustRegionBoundaries(region, xmin, xmax, ymin, ymax)
        hist, histIntegral, integral_error = data.get2DHistoIntegral(histName, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, showEvents=True)
        print(f"hist integral is - {hist.Integral()} in region - {region} and from the function is  - {histIntegral}, xmin - {xmin}, xmax - {xmax}, ymin - {ymin}, ymax - {ymax}")
        hist_dict[region].Sumw2()
        hist_dict[region].SetBinContent(1, histIntegral)
        hist_dict[region].SetBinError(1, integral_error)
    return hist_dict


def get_ABCD_counts_by_SVJ(Data,sgData, bgData, ABCDHistoVar, maincut, SVJBins,factor):
    """
    Extracts the ABCD region counts split by the number of SVJ jets for both signal and background.
    """
    counts_by_SVJ = {}
    error_by_SVJ = {}
    for data_list, label in [(Data,"Data"),(sgData, "Signal"), (bgData, "Background")]:
        for data in data_list:  # Iterate over individual dataset objects
            for svj_bin in SVJBins.items():
                #hist_dict = GetABCDhistPerSVJBin(data, ABCDHistoVar, maincut, svj_bin)
                hist_dict = GetABCDhistPerSVJBin(data, ABCDHistoVar, maincut, svj_bin,factor)
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

def compute_ABCD_prediction(Data, sgData, bgData, ABCDHistoVar, maincut, factors):
    factor_results = {}  
    error_results_by_factor = ()
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

    for factor in factors:
        print(f"\nProcessing factor: {factor}")
        SVJBins = GetSVJbins(factor)
        counts_by_SVJ,error_by_SVJ = get_ABCD_counts_by_SVJ(Data, sgData, bgData, ABCDHistoVar, maincut, SVJBins, factor)
        factor_results[factor] = counts_by_SVJ
            
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
            relative_error = np.sqrt((NB_err / N_B) ** 2 + (NC_err / N_C) ** 2 + (ND_err / N_D) ** 2)
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

            print(f"[{label}] Factor: {factor}, SVJ Bin: {svj_bin} | Observed A: {N_A},Observed ERR A: {NA_err},Predicted A: {N_A_pred}, Predicted Err {NA_pred_err}")
            #print(f"0SVJ obs len: {len(obs_0SVJ_data)}, 1SVJ obs len: {len(obs_1SVJ_data)}, 2PSVJ obs len: {len(obs_2PSVJ_data)}")

    return (
        factor_results,
        obs_0SVJ_data, pred_0SVJ_data, obs_0SVJ_sig, pred_0SVJ_sig, obs_0SVJ_bg, pred_0SVJ_bg,
        obs_1SVJ_data, pred_1SVJ_data, obs_1SVJ_sig, pred_1SVJ_sig, obs_1SVJ_bg, pred_1SVJ_bg,
        obs_2PSVJ_data, pred_2PSVJ_data, obs_2PSVJ_sig, pred_2PSVJ_sig, obs_2PSVJ_bg, pred_2PSVJ_bg,
        obserr_0SVJ_data, obserr_0SVJ_sig, obserr_0SVJ_bg,prederr_0SVJ_data, prederr_0SVJ_sig, prederr_0SVJ_bg,
        obserr_1SVJ_data, obserr_1SVJ_sig, obserr_1SVJ_bg,prederr_1SVJ_data, prederr_1SVJ_sig, prederr_1SVJ_bg,
        obserr_2PSVJ_data, obserr_2PSVJ_sig, obserr_2PSVJ_bg,prederr_2PSVJ_data, prederr_2PSVJ_sig, prederr_2PSVJ_bg


    )

def plot_ABCD_ratios(
        year,
        factor_results, 
        obs_0SVJ_data, pred_0SVJ_data, obs_0SVJ_sig, pred_0SVJ_sig, obs_0SVJ_bg, pred_0SVJ_bg,
        obs_1SVJ_data, pred_1SVJ_data, obs_1SVJ_sig, pred_1SVJ_sig, obs_1SVJ_bg, pred_1SVJ_bg,
        obs_2PSVJ_data, pred_2PSVJ_data, obs_2PSVJ_sig, pred_2PSVJ_sig, obs_2PSVJ_bg, pred_2PSVJ_bg,
        obserr_0SVJ_data, obserr_0SVJ_sig, obserr_0SVJ_bg,prederr_0SVJ_data, prederr_0SVJ_sig, prederr_0SVJ_bg,
        obserr_1SVJ_data, obserr_1SVJ_sig, obserr_1SVJ_bg,prederr_1SVJ_data, prederr_1SVJ_sig, prederr_1SVJ_bg,
        obserr_2PSVJ_data, obserr_2PSVJ_sig, obserr_2PSVJ_bg,prederr_2PSVJ_data, prederr_2PSVJ_sig, prederr_2PSVJ_bg,
        factors, output_dir):
    """
    Plots the ratio of observed to predicted A values for each SVJ type and SCJ category separately,
    and adds a subplot showing the difference between Data Ratio and Background Ratio with a rectangular bottom plot.
    """

    # Define x-axis values
    values = 0.5 / np.array(factors)
    #values = 350 / np.array(factors)

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

    # Function to plot and save the main plot and the difference plot
    def plot_and_save(x_values,non_closures,nonclosure_err, diff_ratios,diff_error, labels, colors, filename):
        plt.figure(figsize=(15, 12))  # Increase figure size to make room for the rectangular bottom plot
        
        # Create a gridspec for controlling subplot layout
        gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])  # Top plot takes 2/3 of the space, bottom plot 1/3

        # Plot the main ratio plot (top subplot)
        ax0 = plt.subplot(gs[0])  # Top subplot
        for non_closure, error, label, color in zip(non_closures, nonclosure_err, labels, colors):
            if len(non_closure) > 0:  # Avoid empty plots
                ax0.errorbar(x_values, non_closure, yerr=error, fmt='o', color=color, label=label, capsize=5)
                ax0.plot(x_values, np.zeros_like(x_values), linestyle='dashed', color='black', linewidth=2)

        hep.cms.label(rlabel="")

        ax0.set_xlabel("Boundary Value")
        ax0.set_ylabel("Non-Closure")
        ax0.set_title(f'Control Region Non Closure split by DNN (VR II) {year}', fontsize=15.5)
        ax0.set_ylim(-0.6, 0.6)
        ax0.legend()
        ax0.grid(True)

        # Plot the difference plot (bottom subplot)
        ax1 = plt.subplot(gs[1])  # Bottom subplot
        for diff_ratio,error, label, color in zip(diff_ratios,diff_error,labels, colors):
            if len(diff_ratio) > 0:  # Avoid empty plots
                #ax1.scatter(x_values, diff_ratio, marker='o', color=color, label=label)
                ax1.errorbar(x_values,diff_ratio, yerr=error, fmt='o', color='black', label=label,capsize=5)
                ax1.plot(x_values, np.zeros_like(x_values), linestyle='dashed', color='black', linewidth=2)
        
        ax1.set_xlabel("Boundary Value")
        ax1.set_ylabel("Data - Background Sim ")
        ax1.set_title(f'Difference in Data and Background Ratios  {year}', fontsize=16)
        ax1.set_ylim(-0.7, 0.7)
        #ax1.legend()
        ax1.grid(True)

        for ax in [ax0, ax1]:
            for spine in ax.spines.values():
                spine.set_linewidth(2)  # Make borders bold

        ax0.yaxis.set_major_locator(MultipleLocator(0.1))  # Major ticks
        ax0.yaxis.set_minor_locator(MultipleLocator(0.05))  # Minor ticks

        ax1.yaxis.set_major_locator(MultipleLocator(0.1))
        ax1.yaxis.set_minor_locator(MultipleLocator(0.05))

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
                "Ratio_Data_0SVJ.jpg")

    plot_and_save(values, 
                [ratio_1SVJ_data, ratio_1SVJ_bg],
                [errbars_1SVJ_data, errbars_1SVJ_bg],
                [diff_1SVJ_data_bg],[errbars_1SVJ_diff], 
                ["Data 1SVJ", "Background MC 1SVJ"], 
                ["b", "r"], 
                "Ratio_Data_1SVJ.jpg")
    print(f"data {ratio_1SVJ_data} err bars {errbars_1SVJ_data}")
    plot_and_save(values, 
                [ratio_2PSVJ_data, ratio_2PSVJ_bg],
                [errbars_2PSVJ_data, errbars_2PSVJ_bg],
                [diff_2PSVJ_data_bg],[errbars_2PSVJ_diff], 
                ["Data 2SVJ", "Background MC 2SVJ"], 
                ["b", "r"], 
                "Ratio_Data_2PSVJ.jpg")
    

    print(f'Boundary Values {values[::-1]}')


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
    CRCuts = ["data/MC"]         
    maincuts = CRCuts
    Data, sgData, bgData = getData( options.dataset + "/", 1.0, year)

    signal_factors = np.linspace(1/2, 1, 20)
    controlregion_factors = [1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2]


    factors = controlregion_factors
    '''
    if (factors == signal_factors).all():
        output_dir = 'Noncloser/SignalRegion'
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    elif factors == controlregion_factors:
        output_dir = 'Noncloser/ControlRegion'
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    '''
    output_dir = 'Nonclosure/VRII-DNN/ControlRegion/'
    factor_results,obs_0SVJ_data, pred_0SVJ_data, obs_0SVJ_sig, pred_0SVJ_sig, obs_0SVJ_bg, pred_0SVJ_bg,obs_1SVJ_data, pred_1SVJ_data, obs_1SVJ_sig, pred_1SVJ_sig, obs_1SVJ_bg, pred_1SVJ_bg,obs_2PSVJ_data, pred_2PSVJ_data, obs_2PSVJ_sig, pred_2PSVJ_sig, obs_2PSVJ_bg, pred_2PSVJ_bg,obserr_0SVJ_data, obserr_0SVJ_sig, obserr_0SVJ_bg,prederr_0SVJ_data, prederr_0SVJ_sig, prederr_0SVJ_bg,obserr_1SVJ_data, obserr_1SVJ_sig, obserr_1SVJ_bg,prederr_1SVJ_data, prederr_1SVJ_sig, prederr_1SVJ_bg,obserr_2PSVJ_data, oberr_2PSVJ_sig, obserr_2PSVJ_bg,prederr_2PSVJ_data, prederr_2PSVJ_sig, prederr_2PSVJ_bg = compute_ABCD_prediction(Data, sgData, bgData, "h_METvsDNN", "_pre_", factors)
    #print('factor_results ',factor_results)
    #print('observed_A',observed_A)
    #print('ratio_A',ratio_A)
    plot_ABCD_ratios(
        year,
        factor_results, 
        obs_0SVJ_data, pred_0SVJ_data, obs_0SVJ_sig, pred_0SVJ_sig, obs_0SVJ_bg, pred_0SVJ_bg,
        obs_1SVJ_data, pred_1SVJ_data, obs_1SVJ_sig, pred_1SVJ_sig, obs_1SVJ_bg, pred_1SVJ_bg,
        obs_2PSVJ_data, pred_2PSVJ_data, obs_2PSVJ_sig, pred_2PSVJ_sig, obs_2PSVJ_bg, pred_2PSVJ_bg,
        obserr_0SVJ_data, obserr_0SVJ_sig, obserr_0SVJ_bg,prederr_0SVJ_data, prederr_0SVJ_sig, prederr_0SVJ_bg,
        obserr_1SVJ_data, obserr_1SVJ_sig, obserr_1SVJ_bg,prederr_1SVJ_data, prederr_1SVJ_sig, prederr_1SVJ_bg,
        obserr_2PSVJ_data, oberr_2PSVJ_sig, obserr_2PSVJ_bg,prederr_2PSVJ_data, prederr_2PSVJ_sig, prederr_2PSVJ_bg,
        factors, output_dir)
 

if __name__ == '__main__':
    main()



 # TODO: Add a function that reads the SVJbin edges
# TODO: Add a function to make table of number of events in the ABCD region, TF, Validation
# TODO: Add a function to make a table in terms of SVJbins vs bkgs.
# TODO: Remove the 0 and 2 from the ratio plot pad 2
# TODO: Redraw axis properly for the dummy.
# TODO: Remove the grid from the pad1. 