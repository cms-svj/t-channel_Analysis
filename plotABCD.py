import ROOT
ROOT.gROOT.SetBatch(True)
import math
import utils.DataSetInfo as info
import optparse
import copy
import os
from array import array
import numpy as np
import pandas as pd
import utils.CMS_lumi as CMS_lumi
from utils.var import var as vars
import plotStack as plotutils
import plotHist as plthist

# Global variable to write the CMS text 
AddCMSText = True

def dataframe_to_latex(df):
    '''Convert the DataFrame to LaTeX'''
    return df.to_latex(index=True)
   

def predictedNA(B,C,D,E=0,F=0,extended=True):
    '''Calculate the Prediction for region A'''
    if extended:
        if D*E*B == 0:
            return 0
        else:
            return ((C*B/D)**2)*(F/(E*B))
    else:
        if D == 0:
            return 0
        else:
            return B*C/D

def closureError(predicted,true):
    '''Calculte Closure error'''
    if true == 0: 
        return 0
    else:
        return abs(predicted-true)/true

def AddCMSLumiText(canvas, year, isExtraText= False, extraText="Preliminary"):
    '''Add CMS Lumi Text on the Plots'''
    if year == "2017":
            lumi = "41.5"
    elif year == "2016":
        lumi = "59.7"
    else:
        lumi = "59.7"
    CMS_lumi.writeExtraText = isExtraText
    CMS_lumi.extraText = extraText
    CMS_lumi.lumi_sqrtS = lumi + " fb^{-1} (13 TeV)"
    iPeriod = 0
    iPos = 0
    CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)
    canvas.cd()
    canvas.Update()
    canvas.RedrawAxis()


def getData(path, scale=1.0, year = "2018"):
    '''Uncomment the files that are to be included in the plots'''
    Data = [
        # info.DataSetInfo(basedir=path, fileName="2017_Data.root",        sys= -1.0, label="Data",        scale=scale),
    ]
    bgData = [
        info.DataSetInfo(basedir=path, fileName=year+"_ST.root",              label="Single top",              scale=scale, color=(ROOT.kRed + 1)),
        info.DataSetInfo(basedir=path, fileName=year+"_TTJets.root",          label="t#bar{t}",                scale=scale, color=(ROOT.kBlue - 6)),
        info.DataSetInfo(basedir=path, fileName=year+"_ZJets.root",             label="Z#rightarrow#nu#nu+jets",    scale=scale, color=(ROOT.kGray + 1)),
        info.DataSetInfo(basedir=path, fileName=year+"_WJets.root",              label="W+jets",                    scale=scale, color=(ROOT.kYellow + 1)),
        info.DataSetInfo(basedir=path, fileName=year+"_QCD.root",               label="QCD",                        scale=scale, color=(ROOT.kGreen + 1)),
    ]
    sgData = [

        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="baseline", scale=scale, color=ROOT.kBlack),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-600.root",     label="mMed_600",  scale=scale, color=ROOT.kViolet+2),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-800.root",     label="mMed_800",  scale=scale, color=ROOT.kRed),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-1000.root",    label="mMed_1000", scale=scale, color=ROOT.kMagenta + 1),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-1500.root",    label="mMed_1500", scale=scale, color=ROOT.kGray+4),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-3000.root",    label="mMed_3000", scale=scale, color=ROOT.kCyan),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-4000.root",    label="mMed_4000", scale=scale, color=ROOT.kBlue + 1),
        
    ]
    return Data, sgData, bgData

def GetABCDregions(met, dnn):
    '''Defining the A B C D regions'''
    regions = [ ("A", met, 20000, dnn, 1.0),
                    ("B", met, 20000, 0, dnn),
                    ("C", 0, met, dnn, 1.0),
                    ("D", 0, met, 0, dnn)
                ]
    return regions

def GetABCDhist(data, ABCDhistoVar, maincut, SVJbins, isStack=False):
    """Returns the individual A,B,C,D histogram for a given data( No data list ), SVJbins, maincut and histovar. It loops over SVJbins and returns 4 histograms."""
    xmin, binwidth = int(list(SVJbins.keys())[0][0]),len(SVJbins)
    xmax = xmin + binwidth
    print(f"xmin = {xmin}, xmax = {xmax} , binwidth =  {binwidth}")
    hist_dict = {region: ROOT.TH1F(f"{region}_{data.label_}_{maincut}", f"{region}_{data.label_}_{maincut}", binwidth, xmin, xmax) for region in ['A', 'B', 'C', 'D']}

    for i, SVJ in enumerate(SVJbins.keys()):
        histName = ABCDhistoVar + maincut + SVJ
        dnn, met = SVJbins[SVJ]
        regions = GetABCDregions(met, dnn)
        for region, xmin, xmax, ymin, ymax in regions:
            histIntegral = data.get2DHistoIntegral(histName, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, showEvents=True)[1]
            # hist_dict[region].SetBinContent(int(SVJ[0])+1, histIntegral ) #TODO: Add integral error and set it as bin error
            hist_dict[region].SetBinContent(i+1,histIntegral) # Bin number starts from 1.

    if isStack:
        for hist in hist_dict.values():
            hist.SetFillColor(data.getColor())
            hist.SetFillStyle(3001)

    return hist_dict

def GetABCDhistDict(dataList, ABCDhistoVar, maincut, SVJbins,isStack=False):
    '''Returns a dictionary containing A,B,C,D histograms for a each files passed in the data list'''
    ABCDhistDict = {}
    for data in dataList:
        histograms = GetABCDhist(data, ABCDhistoVar, maincut, SVJbins,isStack=isStack)
        ABCDhistDict.update({data.fileName : [histograms['A'], histograms['B'], histograms['C'], histograms['D']]})
    return ABCDhistDict

def SumABCDhistList(ABCDhistList, skipList=None):
    '''Returns A,B,C,D histograms which are summed over all the files in the data list. It skips the files provided in the skipList'''
    firstPass = True
    for fileName, histograms in ABCDhistList.items():
        if skipList != None and fileName in skipList:
            continue
        if firstPass:
            sumHistograms = [hist.Clone() for hist in histograms]
            for sumHist in sumHistograms:
                sumHist.Reset()
            firstPass = False
        print(f"fileName = {fileName}")
        for sumHist, hist in zip(sumHistograms, histograms):
            sumHist.Add(hist)
    return sumHistograms

def MergeSingleABCDHist(A_hist, B_hist, C_hist, D_hist, SVJbins, filename):
    """
    Merges A, B, C, D histograms into a single histogram for a single file.
    
    Parameters:
        A_hist: Histogram for region A.
        B_hist: Histogram for region B.
        C_hist: Histogram for region C.
        D_hist: Histogram for region D.
        SVJbins: Dictionary containing the SVJ bins information for labeling.
        filename: Name of the file used for naming the histogram.
    
    Returns:
        merged_hist: The merged histogram for the given file.
    """
    # Get the number of bins from A histogram (assuming A, B, C, D have the same number of bins)
    nBins_A = A_hist.GetNbinsX()
    bin_width = A_hist.GetBinWidth(1)
    
    # The total number of bins in the new histogram is 4 * number of bins in A (since A, B, C, D have same size)
    nBins_total = nBins_A * 4
    
    # Create the new merged histogram
    merged_hist = ROOT.TH1F(f"merged_{filename}", f"merged_{filename}", nBins_total, 0, nBins_total * bin_width)
    
    # Loop over all bins and fill the new histogram with the content of A, B, C, and D histograms
    for i in range(1, nBins_A + 1):
        # Set bin content for A region
        merged_hist.SetBinContent(i, A_hist.GetBinContent(i))
        merged_hist.GetXaxis().SetBinLabel(i, f"SVJ{i-1}")
        
        # Set bin content for B region (offset by nBins_A)
        merged_hist.SetBinContent(i + nBins_A, B_hist.GetBinContent(i))
        merged_hist.GetXaxis().SetBinLabel(i + nBins_A, f"SVJ{i-1}")
        
        # Set bin content for C region (offset by 2 * nBins_A)
        merged_hist.SetBinContent(i + 2 * nBins_A, C_hist.GetBinContent(i))
        merged_hist.GetXaxis().SetBinLabel(i + 2 * nBins_A, f"SVJ{i-1}")
        
        # Set bin content for D region (offset by 3 * nBins_A)
        merged_hist.SetBinContent(i + 3 * nBins_A, D_hist.GetBinContent(i))
        merged_hist.GetXaxis().SetBinLabel(i + 3 * nBins_A, f"SVJ{i-1}")
    
    return merged_hist


def MergeABCDHistogramsForAllFiles(ABCDhistDict, SVJbins):
    """
    Loops over the ABCDhistDict and merges histograms for all files.
    
    Parameters:
        ABCDhistDict: Dictionary with the structure {filename: [A_hist, B_hist, C_hist, D_hist]}
        SVJbins: Dictionary containing the SVJ bins information.
    
    Returns:
        merged_hist_dict: Dictionary containing the merged histograms for each file.
    """
    merged_hist_dict = {}

    # Loop over each entry in the ABCDhistDict (each file)
    for filename, histograms in ABCDhistDict.items():
        A_hist, B_hist, C_hist, D_hist = histograms
        
        # Use the MergeSingleABCDHist function to merge the histograms for this file
        merged_hist = MergeSingleABCDHist(A_hist, B_hist, C_hist, D_hist, SVJbins, filename)
        
        # Store the merged histogram in the dictionary
        merged_hist_dict[filename] = merged_hist
    
    return merged_hist_dict


def MergeLastTwoBins(hist): # TODO: Generalize the function to work with any bins
    num_bins = hist.GetNbinsX()
    print(f"In merge - num_bins = {num_bins}")
    new_hist = ROOT.TH1D(hist.GetName() + "_merged", hist.GetTitle(), num_bins - 1, 0, num_bins - 1)
    for i in range(1, num_bins):
        new_hist.SetBinContent(i, hist.GetBinContent(i))
        new_hist.SetBinError(i, hist.GetBinError(i))
    # Merge the last two bins
    last_bin_content = hist.GetBinContent(num_bins) + hist.GetBinContent(num_bins - 1)
    last_bin_error = (hist.GetBinError(num_bins)**2 + hist.GetBinError(num_bins - 1)**2)**0.5
    new_hist.SetBinContent(num_bins - 1, last_bin_content)
    new_hist.SetBinError(num_bins - 1, last_bin_error)
    return new_hist

def MergeAllLastTwoBins(histlist):
    '''Merge last two bins for all the histograms provided'''
    merged_histogram_list = []
    for hist in histlist:
        merged_hist = MergeLastTwoBins(hist)
        merged_histogram_list.append(merged_hist)
    return merged_histogram_list


def rebinHistogramstoABCD(histList):
    """Takes a list of histograms and returns a new list of histograms rebinned according to the bins of A, B, C, and D."""# TODO: Figure out the issue, not rebinning properly as desired
    newHistList = []
    for hist in histList:
        newHist = hist.Clone()
        newHist.Reset()
        for bin in range(1, hist.GetNbinsX() + 1):
            newHist.Fill(hist.GetBinContent(bin))
        newHistList.append(newHist)
    return newHistList

def SetupGPad(leftMargin=0.15, rightMargin=0.05, topMargin=0.08, bottomMargin=0.12, logY = False):
    ROOT.gPad.Clear()
    ROOT.gPad.SetLeftMargin(leftMargin)
    ROOT.gPad.SetRightMargin(rightMargin)
    ROOT.gPad.SetTopMargin(topMargin)
    ROOT.gPad.SetBottomMargin(bottomMargin)
    ROOT.gPad.SetTicks(1,1)
    ROOT.gPad.SetLogy(logY)
    ROOT.gStyle.SetOptStat("")
    ROOT.gPad.SetGrid()

def SetupLegend(x1=0.65,y1=0.55,x2=0.85,y2=0.85, fillStyle=0, borderSize=0, lineWidth=1, NColumns=1, textFont=42, textSize=0.06):
    leg = ROOT.TLegend(x1, y1, x2, y2)
    leg.SetFillStyle(fillStyle)
    leg.SetBorderSize(borderSize)
    leg.SetLineWidth(lineWidth)
    leg.SetNColumns(NColumns)
    leg.SetTextFont(textFont)
    ROOT.gStyle.SetLegendTextSize(textSize)
    return leg    

def SetupBeforeRatioPlot(hist1, nDivision=3, xtitle = "nSVJ" ,title = "", xtitleSize=1.5, xtitleOffset=1,isStack=False):
    hist1.SetTitle(title)
    if isStack:
        if hist1.GetNhists() > 0:
            hist1.GetHists().At(0).GetXaxis().SetTitle(xtitle)
            hist1.GetHists().At(0).GetXaxis().SetNdivisions(nDivision)
            hist1.GetHists().At(0).GetXaxis().SetTitleSize(xtitleSize)
            hist1.GetHists().At(0).GetXaxis().SetTitleOffset(xtitleOffset)
    else:
        hist1.GetXaxis().SetTitle(xtitle)
        hist1.GetXaxis().SetNdivisions(nDivision)
        hist1.GetXaxis().SetTitleSize(xtitleSize)
        hist1.GetXaxis().SetTitleOffset(xtitleOffset)
    return hist1


def SetupRatioPlot(ratio, yRange=[10**-1,10**5], yUppertitle="Events",yLowertitle="Predict/MC",lowerYMax=2, xLabelSize=0.06, isCenterLabel=True, separationMargin=0.0, splitFraction=0.3, markerStyle=8,is4plots=False):
    """Makes the Ratio Plot and beautifys it"""
    if is4plots:
        # ratio.SetLeftMargin()
        ratio.SetRightMargin(0)
        ratio.SetUpTopMargin(0)
        ratio.SetUpBottomMargin(0)
        ratio.SetLowTopMargin(0)
    
    ratio.GetUpperRefYaxis().SetRangeUser(yRange[0],yRange[1])
    ratio.GetUpperRefYaxis().SetTitle(yUppertitle)
    # ratio.GetUpperRefYaxis().SetTitleOffset(1.0)
    # ratio.GetUpperRefYaxis().SetTitleSize(0.08)
    # ratio.GetUpperRefYaxis().SetLabelSize(0.06)
    
    ratio.GetLowerRefYaxis().SetTitle(yLowertitle)
    # ratio.GetLowerRefYaxis().SetTitleOffset(1.4)
    # ratio.GetLowerRefYaxis().SetTitleSize(0.06)
    # ratio.GetLowerRefYaxis().SetLabelSize(0.06)
    
    ratio.GetLowerRefGraph().SetMarkerStyle(markerStyle)
    # ratio.GetLowerRefGraph().SetMaximum(lowerYMax)
    # ratio.GetLowerRefGraph().SetMinimum(0)
    ratio.GetLowerRefXaxis().CenterLabels(isCenterLabel)
    ratio.GetLowerRefXaxis().SetLabelSize(xLabelSize)
    ratio.SetSplitFraction(splitFraction)
    ratio.SetSeparationMargin(separationMargin)
    



def plotABCD(data,ABCDhistoVar,maincut, SVJbins,outputPath="./",stList=None,SVJbinContent=None,isLogY=False,rebinSVJtoABCD=False,year=2017,scenario = "d0_wp7_p0i0"):
    # TODO: update the code to work with the ABCD functions and update it for the 
    # TODO: Add a option to make ratio plots when data file is passed
    
    ROOT.TH2.AddDirectory(False)
    canvas = ROOT.TCanvas("canvas", "ABCD plot", 1200, 600)
    canvas.Divide(len(SVJbins),1,0.0,0.0)

    stacks,signalhist,dummy,leg = [], [], [], []
    bkghist = [[] for _ in range(len(SVJbins.keys()))]
    
    for iCanvas,SVJ in enumerate(SVJbins.keys()):
        dnn, met = SVJbins[SVJ][0], SVJbins[SVJ][1]
        print(f"iCanvas = {iCanvas}, key = {SVJ}, DNN = {dnn}, met = {met}")

        ABCDregions = GetABCDregions(met, dnn)
        canvas.cd(iCanvas+1)
        
        ROOT.gPad.Clear()
        if(iCanvas == 0):
          ROOT.gPad.SetLeftMargin(0.16)
          if(iCanvas==4 and len(SVJbins)==5):
            ROOT.gPad.SetRightMargin(0.05)
        if(iCanvas!=0):
            ROOT.gPad.SetLeftMargin(0.001)
        ROOT.gPad.SetTopMargin(0.0)
        ROOT.gPad.SetBottomMargin(0.15)
        ROOT.gPad.SetLogy(isLogY)
        ROOT.gPad.SetGrid()

        if iCanvas == 0:
            leg.append(SetupLegend(0.16, 0.78, 0.98, 0.98))
        else:
            leg.append(SetupLegend(0.02, 0.78, 0.98, 0.98))

        stacks.append(ROOT.THStack(f"{iCanvas}",f"{iCanvas}"))
        print(f"***** Working on region {iCanvas} *****")
        firstPass = True
        hMC = None
        labels, bkgcontrib = [], []
        xmin, xmax, binwidth = 1, 5, 4
        regionsum = []
        histName = ABCDhistoVar + maincut + SVJ
        for d in data[1]:
            bkghist[iCanvas].append(ROOT.TH1F(f"h_{d.label_}",f"{d.label_}",binwidth,xmin,xmax))
            regionIntegral = 0
            for i, region in enumerate(ABCDregions, start=1):
                hist, regionIntegral = d.get2DHistoIntegral(histName, xmin=region[1], xmax=region[2], ymin=region[3], ymax=region[4], showEvents=True)
                print(f"region - {region[0]} --- Integral -- {regionIntegral}")
                bkghist[iCanvas][-1].SetBinContent(i, regionIntegral)
                if (stList!=None):
                    # regionstring = "('{}', {}, {}, '{}', {})".format(region)
                    newEntry = stList + [d.label_,region[0],SVJ,regionIntegral]
                    SVJbinContent.loc[len(SVJbinContent.index)]=newEntry
            bkghist[iCanvas][-1].SetFillColor(d.getColor())
            bkghist[iCanvas][-1].SetFillStyle(3001)
            stacks[-1].Add(copy.deepcopy(bkghist[iCanvas][-1])) # Adding the histogram to the Stack
            leg[-1].AddEntry(bkghist[iCanvas][-1], f"{d.label_}", "F")
            if firstPass:
                hMC = copy.deepcopy(bkghist[iCanvas][-1])
                firstPass = False
            else:
                hMC.Add(copy.deepcopy(bkghist[iCanvas][-1]))
        
        ymax=10**10
        ymin=10**-4
        lmax=10**11
        if(iCanvas==len(SVJbins)-1):
            xTitle = f"ABCDregions"
        else:
            xTitle = ""
        yTitle = "Events"    
        dummy.append(ROOT.TH1D(f"dummy_{iCanvas}", f"dummy{iCanvas}", binwidth, xmin, xmax))
        plotutils.setupDummy(dummy[-1], leg[-1], "", xTitle, yTitle, isLogY, xmin, xmax, ymin, ymax, lmax, title=SVJ,isABCD=True)
        dummy[-1].GetXaxis().SetBinLabel(1,"A")
        dummy[-1].GetXaxis().SetBinLabel(2,"B")
        dummy[-1].GetXaxis().SetBinLabel(3,"C")
        dummy[-1].GetXaxis().SetBinLabel(4,"D")
        dummy[-1].Draw("hist")
        if(iCanvas!=0):
            dummy[-1].GetYaxis().SetLabelSize(0)
            dummy[-1].GetYaxis().SetTitleSize(0)
        stacks[-1].Draw("hist F same")

        # Setup signal histo
        linestylenumber = 0
        linestyle = [ROOT.kSolid,ROOT.kDashed,ROOT.kDotted]
        signaltoPlot = ["baseline","mMed_600","mMed_4000"]
        if(data[2]):
            for d in data[2]:
                signalhist.append(ROOT.TH1F(f"h_{d.label_}",f"{d.label_}",binwidth,xmin,xmax) )  
                for i, region in enumerate(ABCDregions, start=1):
                    hist, regionIntegral = d.get2DHistoIntegral(histName, xmin=region[1], xmax=region[2], ymin=region[3], ymax=region[4], showEvents=True)
                    signalhist[-1].SetBinContent(i, regionIntegral)
                    print(f"region - {region[0]} --- Integral -- {regionIntegral}")
                    if (stList!=None):
                        newEntry = stList + [d.label_,region,SVJ,regionIntegral]
                        SVJbinContent.loc[len(SVJbinContent.index)]=newEntry
            if d.label_ in signaltoPlot:
                signalhist[-1].SetLineStyle(linestyle[linestylenumber%3] )
                signalhist[-1].SetLineWidth(3)
                linestylenumber+=1
                signalhist[-1].SetLineColor(d.getColor())
                leg[-1].AddEntry(signalhist[-1], d.label_, "L")
                signalhist[-1].Draw("hist same")    
        leg[-1].Draw("same")
        leg[-1].SetHeader(f"{SVJ} [{met},{dnn}]")
        canvas.Update()
        ROOT.gPad.RedrawAxis()
        ROOT.gPad.RedrawAxis("G") 

    AddCMSLumiText(canvas, year, isExtraText=True)
    canvas.cd()
    canvas.Update()
  
    print("saving the histogram")
    savestring = ABCDhistoVar+maincut
    canvas.SaveAs(outputPath+"/"+savestring+"ABCD.png")
    print("CLosing the canvas")
    canvas.Close()
    del canvas
    # del leg
    del hMC


def PredictPlotABCD(data,ABCDhistoVar,maincut, SVJbins,outputPath="./",stList=None,ABCDPredList=None,isLogY=False,isNoSum=False,year=2017,scenario = "d0_wp7_p0i0"): #TODO: Make this function compatible with the functions
     #TODO: Update this function with the GetABCDhistDict Functions
    '''Takes in multiple data and make the Prediction plots, if there is only one bkg passed then makes prediction only for that bkg, else it sums all the bkg'''
    ROOT.TH2.AddDirectory(False)
    Ahistlist, Bhistlist, Chistlist, Dhistlist = [],[],[],[]
    # xmin, xmax, binwidth = 1, 5, 4
    xmin, binwidth = int(list(SVJbins.keys())[0][0]),len(SVJbins)
    xmax = xmin + binwidth
    for d in data:
        Ahistlist.append(ROOT.TH1F(f"A_{d.label_}",f"A_{d.label_}",binwidth,xmin,xmax))
        Bhistlist.append(ROOT.TH1F(f"B_{d.label_}",f"B_{d.label_}",binwidth,xmin,xmax))
        Chistlist.append(ROOT.TH1F(f"C_{d.label_}",f"C_{d.label_}",binwidth,xmin,xmax))
        Dhistlist.append(ROOT.TH1F(f"D_{d.label_}",f"D_{d.label_}",binwidth,xmin,xmax))

        for SVJ in SVJbins.keys():
            histName = ABCDhistoVar + maincut + SVJ
            dnn, met = SVJbins[SVJ][0], SVJbins[SVJ][1]
            A,B,C,D = 0,0,0,0
            A = d.get2DHistoIntegral(histName, xmin=met, xmax=10000, ymin=dnn, ymax=1.0, showEvents=True)[1]
            B = d.get2DHistoIntegral(histName, xmin=met, xmax=10000, ymin=0, ymax=dnn, showEvents=True)[1]
            C = d.get2DHistoIntegral(histName, xmin=0, xmax=met, ymin=dnn, ymax=1.0, showEvents=True)[1]
            D = d.get2DHistoIntegral(histName, xmin=0, xmax=met, ymin=0, ymax=dnn, showEvents=True)[1]
            Ahistlist[-1].SetBinContent(int(SVJ[0]),A)
            Bhistlist[-1].SetBinContent(int(SVJ[0]),B)
            Chistlist[-1].SetBinContent(int(SVJ[0]),C)
            Dhistlist[-1].SetBinContent(int(SVJ[0]),D)

    sumA_histogram = ROOT.TH1F("sumA", "Sum of A region", binwidth, xmin, xmax)
    sumB_histogram = ROOT.TH1F("sumB", "Sum of B region", binwidth, xmin, xmax)
    sumC_histogram = ROOT.TH1F("sumC", "Sum of C region", binwidth, xmin, xmax)
    sumD_histogram = ROOT.TH1F("sumD", "Sum of D region", binwidth, xmin, xmax)

    for histA, histB, histC, histD in zip(Ahistlist, Bhistlist, Chistlist, Dhistlist):
        sumA_histogram.Add(histA)
        sumB_histogram.Add(histB)
        sumC_histogram.Add(histC)
        sumD_histogram.Add(histD)
    
    ratio_histogram = sumC_histogram.Clone("ratio_histogram")
    ratio_histogram.Divide(sumD_histogram)
    ymax = 10**4
    ymin = 10**-1
    xTitle = "SVJ bins"
    yTitle = "Events"
    for d,Ahist,Bhist,Chist,Dhist in zip(data,Ahistlist,Bhistlist,Chistlist,Dhistlist):
        PredictionAhist = Bhist.Clone(f"predictionA_{d.label_}")
        if isNoSum:
            CbyD = Chist.Clone(f"CbyD{d.label_}")
            CbyD.Divide(Dhist)
            PredictionAhist.Multiply(CbyD)
        else:
            print("Errors - this is invoked")
            PredictionAhist.Multiply(ratio_histogram)

        for SVJ in SVJbins.keys():
            A = Ahist.GetBinContent(int(SVJ[0]))
            PredA = PredictionAhist.GetBinContent(int(SVJ[0]))
            error = closureError(PredA,A)
            if stList!= None:
                newEntry = stList + [d.label_,SVJ,A, PredA, error, Bhist.GetBinContent(int(SVJ[0])), Chist.GetBinContent(int(SVJ[0])),Dhist.GetBinContent(int(SVJ[0]))]
                ABCDPredList.loc[len(ABCDPredList.index)] = newEntry 
        
        c1 = ROOT.TCanvas(f"c_{d.label_}",f"c_{d.label_}",800,700)
        c1.cd()
        SetupGPad(logY=isLogY)
        leg = SetupLegend()
        PredictionAhist.SetTitle("")
        PredictionAhist.GetXaxis().SetTitle(xTitle)
        PredictionAhist.GetXaxis().SetNdivisions(4)
       
        ratio = ROOT.TRatioPlot(PredictionAhist,Ahist)
        ratio.SetH1DrawOpt("histe")
        ratio.SetH2DrawOpt("histe same")
        ratio.Draw()
        Ahist.SetLineColor(ROOT.kRed)
        SetupRatioPlot(ratio)
        leg.AddEntry(PredictionAhist,f"Pred {d.label_}","L")
        leg.AddEntry(Ahist,f"MC {d.label_}","L")
        ratio.GetUpperPad().cd()
        leg.Draw()
        c1.Update()
        ROOT.gPad.RedrawAxis()
        ROOT.gPad.RedrawAxis("G")
        
        if isNoSum:
            savestring = ABCDhistoVar+maincut+f"PredictionPlotRatio_CbyDnoSum_{d.label_}"
        else:
            savestring = ABCDhistoVar+maincut+f"PredictionPlotRatio{d.label_}"
        c1.SaveAs(outputPath+"/"+savestring+".png")
        c1.Close()
        del c1
        
 

def PredictABCD(data,ABCDhistoVar,maincut, SVJbins,outputPath="./",stList=None,ABCDPredList=None,isStack=False,ifSavename=None,isOnlyCbyD=False,isLogY=False,year=2017,scenario = "d0_wp7_p0i0"):
    #TODO: Update this function with the GetABCDhistDict Functions
    ROOT.TH2.AddDirectory(False)
    c1 = ROOT.TCanvas( "c", "c", 800, 700)
    c1.cd()
    ROOT.gPad.Clear()
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.08)
    ROOT.gPad.SetBottomMargin(0.12)
    ROOT.gPad.SetTicks(1,1)
    ROOT.gPad.SetLogy(isLogY)
    ROOT.gStyle.SetOptStat("")
    leg = ROOT.TLegend(0.6, 0.55, 0.85, 0.85)
    #nColumns = 3 if(len(data[1]) >= 3) else 1
    nColumns = 1
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetLineWidth(1)
    leg.SetNColumns(nColumns)
    leg.SetTextFont(42)
    ROOT.gStyle.SetLegendTextSize(0.06)
    # xmin, xmax, binwidth = 1, 5, 4
    xmin, binwidth = int(list(SVJbins.keys())[0][0]),len(SVJbins)
    xmax = xmin + binwidth
    predictHist = ROOT.TH1F(f"predA","predA",binwidth,xmin,xmax)
    AHist = ROOT.TH1F(f"A","A",binwidth,xmin,xmax)
    # stack = ROOT.THStack("BkgPred","BkgPred")
    # bkghist = 
    
    Ahistlist, Bhistlist, Chistlist, Dhistlist = [],[],[],[]
    for d in data:
        Ahistlist.append(ROOT.TH1F(f"A_{d.label_}",f"A_{d.label_}",binwidth,xmin,xmax))
        Bhistlist.append(ROOT.TH1F(f"B_{d.label_}",f"B_{d.label_}",binwidth,xmin,xmax))
        Chistlist.append(ROOT.TH1F(f"C_{d.label_}",f"C_{d.label_}",binwidth,xmin,xmax))
        Dhistlist.append(ROOT.TH1F(f"D_{d.label_}",f"D_{d.label_}",binwidth,xmin,xmax))

        for SVJ in SVJbins.keys():
            histName = ABCDhistoVar + maincut + SVJ
            dnn, met = SVJbins[SVJ][0], SVJbins[SVJ][1]
            
            A,B,C,D = 0,0,0,0
            A = d.get2DHistoIntegral(histName, xmin=met, xmax=10000, ymin=dnn, ymax=1.0, showEvents=True)[1]
            B = d.get2DHistoIntegral(histName, xmin=met, xmax=10000, ymin=0, ymax=dnn, showEvents=True)[1]
            C = d.get2DHistoIntegral(histName, xmin=0, xmax=met, ymin=dnn, ymax=1.0, showEvents=True)[1]
            D = d.get2DHistoIntegral(histName, xmin=0, xmax=met, ymin=0, ymax=dnn, showEvents=True)[1]
            Ahistlist[-1].SetBinContent(int(SVJ[0]),A)
            Bhistlist[-1].SetBinContent(int(SVJ[0]),B)
            Chistlist[-1].SetBinContent(int(SVJ[0]),C)
            Dhistlist[-1].SetBinContent(int(SVJ[0]),D)

    sumA_histogram = ROOT.TH1F("sumA", "Sum of A region", binwidth, xmin, xmax)
    sumB_histogram = ROOT.TH1F("sumB", "Sum of B region", binwidth, xmin, xmax)
    sumC_histogram = ROOT.TH1F("sumC", "Sum of C region", binwidth, xmin, xmax)
    sumD_histogram = ROOT.TH1F("sumD", "Sum of D region", binwidth, xmin, xmax)

    # Add histograms to respective sum histograms
    for histA, histB, histC, histD in zip(Ahistlist, Bhistlist, Chistlist, Dhistlist):
        sumA_histogram.Add(histA)
        sumB_histogram.Add(histB)
        sumC_histogram.Add(histC)
        sumD_histogram.Add(histD)
    
    ratio_histogram = sumC_histogram.Clone("ratio_histogram")
    ratio_histogram.Divide(sumD_histogram)
    if isStack:
        stack = ROOT.THStack("prediction_stack", "Stack of Prediction Histograms")
    PredictionAhistlist = []
    for d,Ahist,Bhist,Chist,Dhist in zip(data,Ahistlist,Bhistlist,Chistlist,Dhistlist):
        PredictionAhist = Bhist.Clone(f"predictionA_{d.label_}")
        if isOnlyCbyD:
            CbyD = Chist.Clone(f"CbyD_{d.label_}")
            CbyD.Divide(Dhist)
            PredictionAhist.Multiply(CbyD)
            # stack.Add(CbyD)
            # PredictionAhist.SetFillColor(d.getColor())
            # PredictionAhist.SetFillStyle(3001)
        else:
            PredictionAhist.Multiply(ratio_histogram)
        
        if isStack:
            stack.Add(PredictionAhist)
            PredictionAhist.SetFillColor(d.getColor())
            PredictionAhist.SetFillStyle(3001)
        else:
            PredictionAhist.SetLineColor(ROOT.kBlue)
            PredictionAhist.SetLineStyle(ROOT.kDashed)
        leg.AddEntry(PredictionAhist,f"Pred {d.label_}","F")
        PredictionAhistlist.append(PredictionAhist)
        

    # Adding the All the data prediction of Region A
    SumPredA = ROOT.TH1F("sumPredA","Sum of PredA",binwidth,xmin,xmax)
    for hist in PredictionAhistlist:
        SumPredA.Add(hist)

    ymax = 10**4
    ymin = 10**-1
    xTitle = "SVJ bins"
    yTitle = "Events"

    if isStack:
        stack.SetTitle("")
        if stack.GetNhists() > 0:
            stack.GetHists().At(0).GetXaxis().SetTitle("nSVJ")
            stack.GetHists().At(0).GetXaxis().SetNdivisions(4)
        stack.SetMaximum(ymax)
        stack.SetMinimum(ymin)
        ratio = ROOT.TRatioPlot(stack,sumA_histogram)
        ratio.SetH1DrawOpt("hist F")
    else:
        SumPredA.SetTitle("")
        SumPredA.GetXaxis().SetTitle(xTitle)
        SumPredA.GetXaxis().SetNdivisions(4)
        ratio = ROOT.TRatioPlot(SumPredA,sumA_histogram)
        ratio.SetH1DrawOpt("histe")
        SumPredA.SetMaximum(ymax)
        SumPredA.SetMinimum(ymin)

    ratio.SetH2DrawOpt("histe same")
    ratio.Draw()
    ratio.GetLowerRefGraph().SetMarkerStyle(8)
    ratio.GetLowerRefYaxis().SetTitle("Predict/MC")
    ratio.GetUpperRefYaxis().SetTitle("Events")
    # ratio.GetUpperRefYaxis().SetRangeUser(ymin,ymax)
    # lower_hist = ratio.GetLowerRefGraph().GetHistogram()
    # lower_hist.GetXaxis().SetTitle("nSVJ")
    sumA_histogram.SetLineColor(ROOT.kRed)
    ratio.SetSplitFraction(0.3)
    ratio.GetLowerRefGraph().SetMaximum(1.5)
    ratio.GetLowerRefXaxis().CenterLabels()  
    ratio.SetSeparationMargin(0.0)  
    leg.AddEntry(sumA_histogram,"MC A","L")
    ratio.GetUpperPad().cd()
    leg.Draw()
    c1.Update()
    ROOT.gPad.RedrawAxis()
    ROOT.gPad.RedrawAxis("G")


    # Writing the ABCDPredList dataframe 
    for d,Ahist,Bhist,Chist,Dhist,PredAhist in zip(data,Ahistlist,Bhistlist,Chistlist,Dhistlist,PredictionAhistlist):
        for SVJ in SVJbins.keys():
            A = Ahist.GetBinContent(int(SVJ[0]))
            PredA = PredAhist.GetBinContent(int(SVJ[0]))
            error = closureError(PredA,A)
            if stList!= None:
                newEntry = stList + [d.label_,SVJ,A, PredA, error, Bhist.GetBinContent(int(SVJ[0])), Chist.GetBinContent(int(SVJ[0])),Dhist.GetBinContent(int(SVJ[0]))]
                ABCDPredList.loc[len(ABCDPredList.index)] = newEntry 

    if isStack and isOnlyCbyD:
        savestring = ABCDhistoVar+maincut+"PredictionPlotRatioStack_CbyD"
    elif not isOnlyCbyD:
        savestring = ABCDhistoVar+maincut+"PredictionPlotRatioStack"
    elif ifSavename!= None:   
        savestring = ABCDhistoVar+maincut+f"PredictionPlotRatio{ifSavename}"
    else:
        savestring = ABCDhistoVar+maincut+f"PredictionPlotRatio"
    c1.SaveAs(outputPath+"/"+savestring+".png")
    c1.Close()
    del c1

            
def formatBinContent(Table,cutName,rowName,columnName,plotOutDir,ABCDFolderName):
    ABCDregion_values = Table['ABCDregion'].unique()
    for region in ABCDregion_values:
        print(region)
        cutTable = Table[(Table['cut'] == cutName) & (Table['ABCDregion'] == region)]
        # print(cutTable)
        pivot_df = cutTable.pivot_table(values='content', index=rowName, columns=columnName, aggfunc='sum')
        # Calculate the percentage contribution for each SVJ bin
        pivot_df = pivot_df.round(2)
        percentage_df = pivot_df.apply(lambda x: x / x.sum() * 100, axis=1)
        # Create a new dataframe with content and percentage
        new_df = pd.DataFrame()
        percentage_df = percentage_df.round(2)
        for col in pivot_df.columns:
            new_df[col] = pivot_df[col].map(str) + ' (' + percentage_df[col].map(str) + '%)'
            # new_df[col] =  percentage_df[col].map(str) + '%'
        # Add a row and a column that contains the sum of the content
        new_df.loc['Total'] = pivot_df.sum()
        new_df['Total'] = pivot_df.sum(axis=1)
        formatted_df = new_df.round(2)
        formatted_df.to_csv("{}/ABCDContentFormatted_{}_{}.csv".format(plotOutDir,cutName,region))
        # formatted_df = formatted_df.loc[len(formatted_df)] = pd.Series(dtype='float64')
        del formatted_df, new_df,percentage_df,pivot_df,cutTable
                    
    
def PlotTF(data, ABCDhistoVar, SRcut, CRcuts, SVJbins, outputPath="./",stList=None, TFContent=None, isLogY=True, year=2017, rebinSVJtoABCD = False, isSumAllCR = False):
    skipData = [f"{year}_QCD.root", f"{year}_ZJets.root"] # These backgrounds do not contribute to Lost Lepton
    SR_ABCDdict = GetABCDhistDict(data, ABCDhistoVar, SRcut, SVJbins)
    CR_dict, sumCR_dict = {}, {}
    for CRcut in CRcuts:
        crABCDdict = GetABCDhistDict(data, ABCDhistoVar, CRcut, SVJbins)
        CR_dict.update(crABCDdict)
        sumCR_dict.update({f"{CRcut}": SumABCDhistList(crABCDdict,skipData)})  

    sumCR_ABCDlist = SumABCDhistList(sumCR_dict)
    sumSR_ABCDlist = SumABCDhistList(SR_ABCDdict,skipData)
    ABCDregion = ["A","B","C","D"]
    leg, ratioTF = [], []
    ROOT.TH2.AddDirectory(False)
    
    canvas = ROOT.TCanvas("TFPlot", "TF plot", 1200,600)
    canvas.Divide(len(sumSR_ABCDlist),1,0,0)
    
    for iCanvas, (SRhist, CRhist, region) in enumerate(zip(sumSR_ABCDlist,sumCR_ABCDlist,ABCDregion)):
        print(f"Working on the Region - SR - {SRhist}, CR - {CRhist}")
        canvas.cd(iCanvas+1)
        # if iCanvas == 0:
        #     SetupGPad(leftMargin=0.16,rightMargin=0, topMargin=0, bottomMargin=0.15, logY=isLogY)
        # elif iCanvas == len(sumSR_ABCDlist)-1:
        #     SetupGPad(leftMargin=0.01, rightMargin=0.05, topMargin=0,bottomMargin=0.15, logY=isLogY)
        # else:
        #     SetupGPad(leftMargin=0.01,rightMargin=0.0,topMargin=0.0, bottomMargin=0.15,logY=isLogY)
    
        leg.append(SetupLegend(0.8,0.8,0.98,0.98))
        # SetupGPad(logY=isLogY)
        ROOT.gPad.SetLogy(isLogY)
        SRhist.SetLineColor(ROOT.kBlue)
        CRhist.SetLineColor(ROOT.kRed)
        leg[-1].AddEntry(CRhist,"CR","l")
        leg[-1].AddEntry(SRhist,"SR","l")
        print(f"Region {iCanvas+1} ------ SRhist bin content - {SRhist.GetBinContent(1)}, CRhist bin content - {CRhist.GetBinContent(1)}")
        if iCanvas == 3:
            SRhist = SetupBeforeRatioPlot(SRhist,nDivision=4, xtitle="nSVJ",xtitleOffset=1.5)
            # SRhist.GetXaxis().SetTitleSize(10)
        else:
            SRhist = SetupBeforeRatioPlot(SRhist, nDivision=4,xtitle="")
        ratioTF.append(ROOT.TRatioPlot(SRhist, CRhist))
        ratioTF[-1].Draw()
        ratioTF[-1].SetH1DrawOpt("histe")
        ratioTF[-1].SetH2DrawOpt("histe")
        
        ylowmax = 4
        yRange = [10**-1,10**6]
        SetupRatioPlot(ratioTF[-1],yRange=yRange, yLowertitle="SR/CR",lowerYMax=ylowmax,is4plots=True)
        
        ratioTF[-1].SetLowBottomMargin(0.4)
        ratioTF[-1].GetUpperPad().cd()
        ratioTF[-1].SetGridlines(0)
        # ratioTF[-1].GetLowerRefGraph().SetMinimum(0)
        # ratioTF[-1].GetLow
        leg[-1].Draw()
        canvas.Update()
        ROOT.gPad.RedrawAxis()

        # Writing the TF content dictionary
        for i in range(0,SRhist.GetNbinsX()):
            if(stList != None):
                if CRhist.GetBinContent(i+1) == 0:
                    tf = 0
                else:
                    tf = SRhist.GetBinContent(i+1)/CRhist.GetBinContent(i+1)
                newEntry  = stList + [region,i,SRhist.GetBinContent(i+1), CRhist.GetBinContent(i+1), tf]
                TFContent.loc[len(TFContent.index)] = newEntry

    if AddCMSText:
        AddCMSLumiText(canvas, year, isExtraText= True)
    if rebinSVJtoABCD:
        savestring = ABCDhistoVar+SRcut+"_"+CRcut+"_TFPlot_inABCDbins"
    else:
        savestring = ABCDhistoVar+SRcut+"_"+CRcut+"_TFPlot_Sum_merged"

    canvas.SaveAs(outputPath+"/"+savestring+".png")
    canvas.Close()
    del canvas
    del leg

def PlotMETDNNRatio(data, histVar2D, maincut, SVJBins, isNorm = False, outputPath="./", year="2017"):
    ROOT.TH2.AddDirectory(False)
    for SVJ in SVJBins.keys():
        histName = histVar2D + maincut + SVJ
        dnn, met = SVJBins[SVJ]
        for d in data:
            # Get the histogram projections
            METhistLowDnn = d.getXProjection(histName, dnn, '<', rebinx=10)
            METhistHighDnn = d.getXProjection(histName, dnn, '>', rebinx=10) 
            METhistLowDnn.SetLineColor(ROOT.kBlue)
            METhistHighDnn.SetLineColor(ROOT.kRed)
            METhistHighDnn.SetLineWidth(2)
            METhistLowDnn.SetLineWidth(2)
            # Normalize the histograms
            plotutils.normHisto(METhistHighDnn,isNorm)
            plotutils.normHisto(METhistLowDnn,isNorm)
            # Working with the canvas
            canvas = ROOT.TCanvas(f"METDNN_{SVJ}_{d.fileName}",f"METDNN {SVJ} {d.fileName}", 800, 700)
            SetupGPad(rightMargin=0.005,topMargin=0.08,bottomMargin=0.12,logY=True)
            ratio = ROOT.TRatioPlot(METhistHighDnn,METhistLowDnn)
            ratio.Draw()
            ratio.SetH1DrawOpt("histe")
            ratio.SetH2DrawOpt("histe")
            SetupBeforeRatioPlot(METhistHighDnn, nDivision=False, xtitle="MET")
            SetupRatioPlot(ratio,yRange=[10**-4,10**1] ,yLowertitle="High/Low",xLabelSize=0.03,isCenterLabel=False)
            leg = SetupLegend(x1=0.55)
            leg.AddEntry(METhistHighDnn, f"dnn_score > {dnn}", "L")
            leg.AddEntry(METhistLowDnn, f"dnn_score < {dnn}", "L")
            legName = d.fileName.replace(".root","")
            leg.SetHeader(f"{legName} {SVJ}")
            ratio.GetUpperPad().cd()
            leg.Draw()
            canvas.Update()
            ROOT.gPad.RedrawAxis()
            # AddCMSLumiText(canvas,year,isExtraText=True)
            savestring = histVar2D+maincut+SVJ+f"{legName}_{dnn}"
            canvas.SaveAs(outputPath+"/"+savestring+".png")
            canvas.Close()
            del canvas, METhistHighDnn, METhistLowDnn, ratio, leg
            

def PlotABCDRatio(data,ABCDhistoVar,maincut, SVJbins,outputPath="./",stList=None,ABCDPredList=None,isStack=False,ifSavename=None,isRatio=False,isLogY=False,year=2018,scenario = "d0_wp7_p0i0"):
    ROOT.TH2.AddDirectory(False)
    if isRatio:
        c1 = ROOT.TCanvas( "c", "c", 800, 700)
        c1, pad1, pad2 = plotutils.createCanvasPads(c1,isLogY)
        pad1.cd()
    else:
        c1 = ROOT.TCanvas( "c", "c", 800, 800)
        c1.cd()
        SetupGPad(leftMargin=0.15,rightMargin=0.05,topMargin=0.08,bottomMargin=0.12,,logY=isLogY)

    if isRatio:
        leg = SetupLegend(0.17,0.7,0.95,0.88,NColumns=2,textSize=0.05)
    else:
        leg = SetupLegend(0.17,0.7,0.95,0.88,NColumns=2,textSize=0.024)

    cut_ABCDdict = GetABCDhistDict(data, ABCDhistoVar, maincut, SVJbins)
    merged_ABCDdist = MergeABCDHistogramsForAllFiles(cut_ABCDdict, SVJbins)

        

 
            
    

    
def main():
    parser = optparse.OptionParser("usage: %prog [options]\n")
    parser.add_option('-b',                 dest='isNormBkg',  action="store_true",                            help="Normalized Background and Signal plots")
    parser.add_option('-d', '--dataset',    dest='dataset',                    default='testHadd_11242020',    help='dataset')
    parser.add_option('-m',                 dest='manySigs',   action="store_true",                            help="Plot ROC curves with many signals vs. QCD")
    parser.add_option('-n',                 dest='isNorm',     action="store_true",                            help="Normalize stack plots")
    parser.add_option('-s',                 dest='onlySig',    action="store_true",                            help="Plot only signals")
    parser.add_option('-y',                 dest='year',       type='string',  default='2018',                 help="Can pass in the run year")
    parser.add_option('-o',                 dest='outputdir',  type='string',                                  help="Output folder name")
    parser.add_option('-w',                 dest='scenario',  type='string',   default='d0_w7p0i0',                               help="Scenario")
    options, args = parser.parse_args()
    scenario = options.scenario
    year = options.year
    
    # change the values for different scenarios
    SVJbins = {
                "0SVJ" : [0.5,150.0],
                "1SVJ" : [0.56,450.0],
                "2PSVJ" : [0.56,350.0],
    }
    ABCDhistoVars = ["METvsDNN"]
    ABCDFolderName = "ABCD"
    SRCut = ["_pre_"]
    CRCuts = ["_cr_muon_","_cr_electron_","_lcr_pre_"]
    maincuts = SRCut #+ CRCuts
    
    # Create the ABCDregions list using a list comprehension
    Data, sgData, bgData = getData( options.dataset + "/", 1.0, year)
    SVJbinContent = pd.DataFrame(columns=["cut","var","source","ABCDregion","SVJbin","content"])
    TFContent = pd.DataFrame(columns=["CRCut","ABCDregion","SVJbin","SR","CR","TF"])
    ABCDPredList = pd.DataFrame(columns = ["cut","var","source","SVJbin","A","Pred A","error" ,"B", "C","D"])
    if options.outputdir:
        plotOutDir = "outputPlots/{}".format(options.outputdir)
    else: 
        plotOutDir = "outputPlots/{}".format(options.dataset)

    
    # myvars = key : ["xlabel", no. of bins, xmin,xmax, npzinfo, flattenInfo, weightName]

    for histName in ABCDhistoVars:
        for maincut in maincuts:
            plotutils.makeDirs(plotOutDir,maincut,ABCDFolderName)
            plotABCDdir = plotOutDir+'/'+ABCDFolderName+'/'+maincut[1:]
            stList = [maincut,histName]
            # PredictPlotABCD(bgData,"h_"+histName,maincut,SVJbins,outputPath=plotABCDdir,stList=stList,ABCDPredList=ABCDPredList,isLogY=True,isNoSum=True,year=year,scenario=scenario)
            # PlotMETDNNRatio(bgData, "h_"+histName,maincut,SVJbins, outputPath=plotABCDdir, year = year, isNorm=True)
            plotABCD((Data,bgData,sgData),"h_"+histName,maincut,SVJbins,outputPath=plotABCDdir,stList=stList,SVJbinContent=SVJbinContent,isLogY=True,year=year,scenario=scenario)
            # PredictABCD(bgData,"h_"+histName,maincut,SVJbins,outputPath=plotABCDdir,stList=stList,ABCDPredList=ABCDPredList,isStack=True, isOnlyCbyD=True, isLogY=True,year=year,scenario=scenario)
            # PredictABCD(bgData,"h_"+histName,maincut,SVJbins,outputPath=plotABCDdir,stList=stList,ABCDPredList=ABCDPredList,isStack=True, isOnlyCbyD=False, isLogY=True,year=year,scenario=scenario)
        
        # for CRCut in CRCuts:
        #     # plotutils.makeDirs(plotOutDir,SRCut, ABCDFolderName)
        # plotTFDir = plotOutDir+'/'+ABCDFolderName
        # stList = ["Sum of CR"]
        # PlotTF(bgData,"h_"+histName, SRCut[0], CRCuts, SVJbins, outputPath=plotTFDir, stList=stList, TFContent=TFContent, isLogY=True, year=year, rebinSVJtoABCD=False )

    # TFContent = TFContent.round(3)
    # TFContent.to_csv("{}/{}/TFContent_Sum_merged.csv".format(plotOutDir,ABCDFolderName))
    # SVJbinContent.to_csv("{}/{}/allSVJBinContents.csv".format(plotOutDir,ABCDFolderName))
    
    # plotting transfer Factor
    
    
    # ABCDPredList = ABCDPredList.round(2)
    # ABCDPredList.to_csv("{}/{}/ABCDPredList_bkgCbyD.csv".format(plotOutDir,ABCDFolderName))
    # for maincut in maincuts:
    #     plotABCDdir = plotOutDir+'/'+ABCDFolderName+'/'+maincut[1:]
    #     formatBinContent(SVJbinContent,maincut,rowName=["SVJbin"],columnName=["source"],plotOutDir=plotABCDdir,ABCDFolderName=ABCDFolderName)

    
    # vars2D = ["METvsnNMedEvent","nNMedEventvsnjetsAK8","nNMedEventvsnTagSVJ","nNMedEventvsnTruSVJ","nTagSVJvsnjetsAK8","nTruSVJvsnjetsAK8","METvsDNN"]
    # for var in vars2D:
    #     for cut in cutsImportant:
    #         plotutils.makeDirs(plotOutDir,cut,"Histo2D")
    #         histo




if __name__ == '__main__':
    main()
