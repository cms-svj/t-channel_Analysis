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
        hist_dict[region].SetBinContent(1, histIntegral)
        hist_dict[region].SetBinError(1, integral_error)
        hist_dict[region].Sumw2()
    return hist_dict


def get_ABCD_counts_by_SVJ(Data,sgData, bgData, ABCDHistoVar, maincut, SVJBins,factor):
    """
    Extracts the ABCD region counts split by the number of SVJ jets for both signal and background.
    """
    counts_by_SVJ = {}
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
                
                if (bin_label, label) not in counts_by_SVJ:
                    counts_by_SVJ[(bin_label, label)] = (0, 0, 0, 0)
                
                counts_by_SVJ[(bin_label, label)] = (
                    counts_by_SVJ[(bin_label, label)][0] + N_A,
                    counts_by_SVJ[(bin_label, label)][1] + N_B,
                    counts_by_SVJ[(bin_label, label)][2] + N_C,
                    counts_by_SVJ[(bin_label, label)][3] + N_D,
                )
    
    return counts_by_SVJ

def compute_ABCD_prediction_by_SVJ(Data,sgData, bgData, ABCDHistoVar, maincut,factor):
    
    SVJBins = GetSVJbins(factor)
    counts_by_SVJ = get_ABCD_counts_by_SVJ(Data,sgData, bgData, ABCDHistoVar, maincut, SVJBins,factor)
    #plot_ABCD_results(counts_by_SVJ)

    for (svj_bin, label), (N_A, N_B, N_C, N_D) in counts_by_SVJ.items():
        print(f'[{label}] A:{N_A} B:{N_B} C:{N_C} D:{N_D}')  # Print with dataset label
        
        if N_D is None or N_D <= 0 or np.isnan(N_D):
            print(f"Warning: N_D is zero or undefined for {label} in SVJ Bin {svj_bin}")
            N_A_pred = float('nan')
        else:
            N_A_pred = (N_B * N_C) / N_D
        
        ratio = N_A / N_A_pred if N_A_pred > 0 else float('nan')
        print(f"[{label}] SVJ Bin: {svj_bin} | Observed A: {N_A}, Predicted A: {N_A_pred}, Ratio: {ratio}")
        
    return counts_by_SVJ,N_A,N_A_pred


def compute_ABCD_prediction(Data, sgData, bgData, ABCDHistoVar, maincut, factors):
    factor_results = {}  # Stores ABCD counts for each factor
    SVJ_labels = ['0SVJ', '1SVJ', '2PSVJ']

    obsA0SVJ, obsA1SVJ, obsA2PSVJ = [], [], []
    predA0SVJ, predA1SVJ, predA2PSVJ = [], [], []
    
    DataobsA0SVJ, DataobsA1SVJ, DataobsA2PSVJ = [], [], []
    DatapredA0SVJ, DatapredA1SVJ, DatapredA2PSVJ = [], [], []
    
    SignalobsA0SVJ, SignalobsA1SVJ, SignalobsA2PSVJ = [], [], []
    SignalpredA0SVJ, SignalpredA1SVJ, SignalpredA2PSVJ = [], [], []
    
    BackgroundobsA0SVJ, BackgroundobsA1SVJ, BackgroundobsA2PSVJ = [], [], []
    BackgroundpredA0SVJ, BackgroundpredA1SVJ, BackgroundpredA2PSVJ = [], [], []

    for factor in factors:
        print(f"\nProcessing factor: {factor}")
        SVJBins = GetSVJbins(factor)
        counts_by_SVJ = get_ABCD_counts_by_SVJ(Data, sgData, bgData, ABCDHistoVar, maincut, SVJBins, factor)
        factor_results[factor] = counts_by_SVJ

        for (svj_bin, label), (N_A, N_B, N_C, N_D) in counts_by_SVJ.items():
            N_A_pred = (N_B * N_C) / N_D if N_D > 0 else float('nan')
            
            print(factor)
            #print("l******abel",label)
            #print(type(label))
            print(svj_bin)
            print(N_A)
            #print(type(svj_bin))
            if svj_bin == "('0SVJ', [0.25, 350.0])":  # 0SVJ bin
                if label == 'Data':
                    DataobsA0SVJ.append(N_A)
                    DatapredA0SVJ.append(N_A_pred)
                elif label == 'Signal':
                    SignalobsA0SVJ.append(N_A)
                    SignalpredA0SVJ.append(N_A_pred)
                elif label == 'Background':
                    BackgroundobsA0SVJ.append(N_A)
                    BackgroundpredA0SVJ.append(N_A_pred)
            elif svj_bin == "('1SVJ', [0.28, 450.0])":  # 1SVJ bin
                if label == 'Data':
                    DataobsA1SVJ.append(N_A)
                    DatapredA1SVJ.append(N_A_pred)
                elif label == 'Signal':
                    SignalobsA1SVJ.append(N_A)
                    SignalpredA1SVJ.append(N_A_pred)
                elif label == 'Background':
                    BackgroundobsA1SVJ.append(N_A)
                    BackgroundpredA1SVJ.append(N_A_pred)
            elif svj_bin == "('2PSVJ', [0.28, 350.0])" :  # 2PSVJ bin
                if label == 'Data':
                    DataobsA2PSVJ.append(N_A)
                    DatapredA2PSVJ.append(N_A_pred)
                elif label == 'Signal':
                    SignalobsA2PSVJ.append(N_A)
                    SignalpredA2PSVJ.append(N_A_pred)
                elif label == 'Background':
                    BackgroundobsA2PSVJ.append(N_A)
                    BackgroundpredA2PSVJ.append(N_A_pred)
            
            print(f"[{label}] Factor: {factor}, SVJ Bin: {svj_bin} | Observed A: {N_A}, Predicted A: {N_A_pred}")
            
    return factor_results, DataobsA0SVJ, DataobsA1SVJ, DataobsA2PSVJ, SignalobsA0SVJ, SignalobsA1SVJ, SignalobsA2PSVJ, BackgroundobsA0SVJ, BackgroundobsA1SVJ, BackgroundobsA2PSVJ, DatapredA0SVJ, DatapredA1SVJ, DatapredA2PSVJ, SignalpredA0SVJ, SignalpredA1SVJ, SignalpredA2PSVJ, BackgroundpredA0SVJ, BackgroundpredA1SVJ, BackgroundpredA2PSVJ


def plot_all_SVJ_ratios(values, DataobsA0SVJ, SignalobsA0SVJ, BackgroundobsA0SVJ, 
                        DataobsA1SVJ, SignalobsA1SVJ, BackgroundobsA1SVJ, 
                        DataobsA2PSVJ, SignalobsA2PSVJ, BackgroundobsA2PSVJ,
                        DatapredA0SVJ, SignalpredA0SVJ, BackgroundpredA0SVJ, 
                        DatapredA1SVJ, SignalpredA1SVJ, BackgroundpredA1SVJ, 
                        DatapredA2PSVJ, SignalpredA2PSVJ, BackgroundpredA2PSVJ,
                        output_dir):
    """
    Function to plot observed vs predicted ratios for all SVJ types: 0SVJ, 1SVJ, 2PSVJ.
    """
    # Define a helper function to calculate ratios
    def calculate_ratios(obs, pred):
        return np.array(obs) / np.array(pred) if len(pred) > 0 else np.zeros(len(obs))

    print(f"  Data OBS SVJ (0SVJ): {DataobsA0SVJ}")
    # Calculate ratios for all three SVJ types
    Data_ratios_0SVJ = calculate_ratios(DataobsA0SVJ, DatapredA0SVJ)
    Signal_ratios_0SVJ = calculate_ratios(SignalobsA0SVJ, SignalpredA0SVJ)
    Background_ratios_0SVJ = calculate_ratios(BackgroundobsA0SVJ, BackgroundpredA0SVJ)

    Data_ratios_1SVJ = calculate_ratios(DataobsA1SVJ, DatapredA1SVJ)
    Signal_ratios_1SVJ = calculate_ratios(SignalobsA1SVJ, SignalpredA1SVJ)
    Background_ratios_1SVJ = calculate_ratios(BackgroundobsA1SVJ, BackgroundpredA1SVJ)

    Data_ratios_2PSVJ = calculate_ratios(DataobsA2PSVJ, DatapredA2PSVJ)
    Signal_ratios_2PSVJ = calculate_ratios(SignalobsA2PSVJ, SignalpredA2PSVJ)
    Background_ratios_2PSVJ = calculate_ratios(BackgroundobsA2PSVJ, BackgroundpredA2PSVJ)

    print("\nRatios for 0SVJ:")
    print(f"  Data Ratios (0SVJ): {Data_ratios_0SVJ}")
    print(f"  Signal Ratios (0SVJ): {Signal_ratios_0SVJ}")
    print(f"  Background Ratios (0SVJ): {Background_ratios_0SVJ}")

    print("\nRatios for 1SVJ:")
    print(f"  Data Ratios (1SVJ): {Data_ratios_1SVJ}")
    print(f"  Signal Ratios (1SVJ): {Signal_ratios_1SVJ}")
    print(f"  Background Ratios (1SVJ): {Background_ratios_1SVJ}")

    print("\nRatios for 2PSVJ:")
    print(f"  Data Ratios (2PSVJ): {Data_ratios_2PSVJ}")
    print(f"  Signal Ratios (2PSVJ): {Signal_ratios_2PSVJ}")
    print(f"  Background Ratios (2PSVJ): {Background_ratios_2PSVJ}")

    print("Ratios calculated for 0SVJ, 1SVJ, and 2PSVJ.")

    # Define the color palette for the plots
    colors = ['blue', 'red', 'green', 'purple', 'orange', 'brown']

    # Plot for 0SVJ
    plt.figure(figsize=(10, 7))
    plt.scatter(values, Data_ratios_0SVJ, marker='o', color=colors[0], label="Data 0SVJ")
    plt.scatter(values, Background_ratios_0SVJ, marker='x', color=colors[1], label="Background Sim 0SVJ")
    plt.xlabel("Boundary Value")
    plt.ylabel("Observed A / Predicted A")
    plt.legend()
    plt.grid(True)
    plt.title("0SVJ")
    plt.savefig(os.path.join(output_dir, "NonCloser_0SVJ.jpg"), dpi=300)
    plt.close()

    # Plot for 1SVJ
    plt.figure(figsize=(10, 7))
    plt.scatter(values, Data_ratios_1SVJ, marker='o', color=colors[0], label="Data 1SVJ")
    plt.scatter(values, Background_ratios_1SVJ, marker='x', color=colors[1], label="Background Sim 1SVJ")
    plt.xlabel("Boundary Value")
    plt.ylabel("Observed A / Predicted A")
    plt.legend()
    plt.grid(True)
    plt.title("1SVJ")
    plt.savefig(os.path.join(output_dir, "NonCloser_1SVJ.jpg"), dpi=300)
    plt.close()

    # Plot for 2PSVJ
    plt.figure(figsize=(10, 7))
    plt.scatter(values, Data_ratios_2PSVJ, marker='o', color=colors[0], label="Data 2PSVJ")
    plt.scatter(values, Background_ratios_2PSVJ, marker='x', color=colors[1], label="Background Sim 2PSVJ")
    plt.xlabel("Boundary Value")
    plt.ylabel("Observed A / Predicted A")
    plt.legend()
    plt.grid(True)
    plt.title("2PSVJ")
    plt.savefig(os.path.join(output_dir, "NonCloser_2PSVJ.jpg"), dpi=300)
    plt.close()

    # Plot the difference between Data and Background for each SVJ type
    def plot_ratio_diff(Data_ratios, Background_ratios, svj_type):
        Ratio_Diff = np.array(Data_ratios) - np.array(Background_ratios)
        plt.figure(figsize=(10, 5))
        plt.scatter(values, Ratio_Diff, marker='o', color=colors[2], label=f"Data - Sim {svj_type}")
        plt.plot(values, np.array([0] * len(values)), linestyle='--', color='r')
        plt.xlabel("Boundary Value")
        plt.ylabel("Observed A / Predicted A")
        plt.legend()
        plt.grid(True)
        plt.title(f"{svj_type} Data - Background")
        plt.savefig(os.path.join(output_dir, f"Data-Sim_{svj_type}.jpg"), dpi=300)
        plt.close()

    plot_ratio_diff(Data_ratios_0SVJ, Background_ratios_0SVJ, "0SVJ")
    plot_ratio_diff(Data_ratios_1SVJ, Background_ratios_1SVJ, "1SVJ")
    plot_ratio_diff(Data_ratios_2PSVJ, Background_ratios_2PSVJ, "2PSVJ")

    print("Saved all plots for 0SVJ, 1SVJ, and 2PSVJ.")







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
    CRCuts = ["_lcr_pre_"]#,"_cr_muon_","_cr_electron_"]
        
    maincuts = [SRCut] + CRCuts
    Data, sgData, bgData = getData( options.dataset + "/", 1.0, year)

    signal_factors = np.linspace(1/2, 1, 20)
    controlregion_factors = [1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2]
    factors = controlregion_factors
    values = values = 0.5 / np.array(factors)
    output_dir = 'Noncloser/ControlRegion'
    #factor_results, observed_A, predicted_A, ratio_A = compute_ABCD_prediction(Data, sgData, bgData, "h_METvsDNN", "_pre_", factors)
    #plot_ABCD_ratios(factor_results,observed_A, predicted_A, factors, GetSVJbins(factors[0]),output_dir)

    factor_results, DataobsA0SVJ, DataobsA1SVJ, DataobsA2PSVJ, SignalobsA0SVJ, SignalobsA1SVJ, SignalobsA2PSVJ, BackgroundobsA0SVJ, BackgroundobsA1SVJ, BackgroundobsA2PSVJ, DatapredA0SVJ, DatapredA1SVJ, DatapredA2PSVJ, SignalpredA0SVJ, SignalpredA1SVJ, SignalpredA2PSVJ, BackgroundpredA0SVJ, BackgroundpredA1SVJ, BackgroundpredA2PSVJ = compute_ABCD_prediction(Data, sgData, bgData, "h_METvsDNN", "_pre_", factors)
    #print(factor_results)
    plot_all_SVJ_ratios(values, DataobsA0SVJ, SignalobsA0SVJ, BackgroundobsA0SVJ, 
                     DataobsA1SVJ, SignalobsA1SVJ, BackgroundobsA1SVJ, 
                     DataobsA2PSVJ, SignalobsA2PSVJ, BackgroundobsA2PSVJ,
                     DatapredA0SVJ, SignalpredA0SVJ, BackgroundpredA0SVJ, 
                     DatapredA1SVJ, SignalpredA1SVJ, BackgroundpredA1SVJ, 
                     DatapredA2PSVJ, SignalpredA2PSVJ, BackgroundpredA2PSVJ,
                     output_dir)
    
 

if __name__ == '__main__':
    main()



 # TODO: Add a function that reads the SVJbin edges
# TODO: Add a function to make table of number of events in the ABCD region, TF, Validation
# TODO: Add a function to make a table in terms of SVJbins vs bkgs.
# TODO: Remove the 0 and 2 from the ratio plot pad 2
# TODO: Redraw axis properly for the dummy.
# TODO: Remove the grid from the pad1. 