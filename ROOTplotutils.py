import ROOT
import numpy as np
import utils.CMS_lumi as CMS_lumi
import math
import os

def makeDirs(plotOutDir,cut,plotType):
    if not os.path.exists(plotOutDir+"/"+plotType+"/"+cut[1:]):
        os.makedirs(plotOutDir+"/"+plotType+"/"+cut[1:])

def StackedHistogram(hist_dict):
    '''Creates a stacked histogram from a dictionary of histograms and returns both the stacked histogram and the summed histogram.'''
    
    # Create a THStack to hold the histograms
    stacked_hist = ROOT.THStack("stacked_hist", "Stacked Histogram")
    
    # Get the first histogram's binning to create an empty histogram for summing
    first_key = next(iter(hist_dict))
    summed_hist = hist_dict[first_key].Clone("summed_hist")
    summed_hist.Reset()  # Clear the content, keep the structure

    # Iterate over the histograms in the dictionary, add them to the stack, and sum them
    for key in hist_dict:
        hist = hist_dict[key]
        stacked_hist.Add(hist)
        summed_hist.Add(hist)  # Sum the histograms together
    
    # Return both the stacked histogram and the summed histogram
    return stacked_hist, summed_hist


def RatioHistogram(numerator, denominator, yTitle="num/den", xTitle = "xTitle", Title = ""):
    '''Creates a ratio histogram by dividing the numerator histogram by the denominator histogram.'''
    
    # Clone the numerator histogram to avoid modifying the original
    ratio_hist = numerator.Clone("ratio_hist")
    
    # Set the title and axis labels (optional)
    ratio_hist.SetTitle(Title)
    ratio_hist.GetYaxis().SetTitle(yTitle)
    ratio_hist.GetXaxis().SetTitle(xTitle)  # Use the same x-axis label as the numerator
    
    # Perform the division
    ratio_hist.Divide(denominator)
    # ratio_hist.Sumw2()
    
    return ratio_hist

def NonclosureHistogram(numerator, denominator, yTitle="|num-den|/den", xTitle="xTitle", Title=""):
    '''Creates a nonclosure histogram: abs(numerator - denominator) / denominator.'''
    nonclosure_hist = numerator.Clone("nonclosure_hist")
    nonclosure_hist.SetTitle(Title)
    nonclosure_hist.GetYaxis().SetTitle(yTitle)
    nonclosure_hist.GetXaxis().SetTitle(xTitle)

    for i in range(1, nonclosure_hist.GetNbinsX() + 1):
        num_val = numerator.GetBinContent(i)
        den_val = denominator.GetBinContent(i)
        if den_val != 0:
            nonclosure = abs(num_val - den_val) / abs(den_val)
        else:
            nonclosure = 0.0
        nonclosure_hist.SetBinContent(i, nonclosure)
        # Error propagation (optional, can be improved for your use case)
        num_err = numerator.GetBinError(i)
        den_err = denominator.GetBinError(i)
        if den_val != 0:
            err = math.sqrt(
                (num_err / abs(den_val)) ** 2 +
                ((num_val - den_val) * den_err / (den_val ** 2)) ** 2
            )
        else:
            err = 0.0
        nonclosure_hist.SetBinError(i, err)
    return nonclosure_hist

def SetupRatioStyle(ratioHist, xTitle, yTitle, yTitleSize=0.15, ymin=0,ymax=2):
    ratioHist.SetLineColor(ROOT.kBlack)
    ratioHist.SetMarkerStyle(20)
    ratioHist.SetMarkerSize(1)
    ratioHist.SetMarkerColor(ROOT.kBlack)
    ratioHist.SetTitle("")
    
    ratioHist.SetMinimum(ymin)
    ratioHist.SetMaximum(ymax)
    ratioHist.SetStats(0)
    # ymax = ratioHist.GetMaximumBin()
    # print("ymax in create Ratio = ",ymax)
    # ratioHist.SetMaximum(ymax)
	# Adjust y-axis settings
    x = ratioHist.GetXaxis()
    y = ratioHist.GetYaxis()
    
    x.SetTitleOffset(0.65)
    y.SetTitleOffset(0.4)
    x.SetTitleSize(0.2)
    y.SetTitleSize(yTitleSize)
    x.SetLabelSize(0.2)
    y.SetLabelSize(0.13)
    x.SetTitle(xTitle)

    y.SetTitle(yTitle)

    # if(x.GetNdivisions() % 100 > 5): x.SetNdivisions(6, 5, 0)

    y.SetNdivisions(505)
    return ratioHist

def AddCMSLumiText(canvas, year, isExtraText= False, extraText="Preliminary", hemPeriod = False):
    '''Add CMS Lumi Text on the Plots'''
    if year == "2017":
        lumi = "41.5"
    elif year == "2016":
        lumi = "35.9"
    elif hemPeriod == "PostHEM" and year == "2018":
        lumi = "38.7"
    elif hemPeriod == "PreHEM" and year == "2018":
        lumi = "21.1"
    elif year == "2018" and hemPeriod == False:
        lumi = "59.7"
    CMS_lumi.writeExtraText = isExtraText
    CMS_lumi.extraText = extraText
    CMS_lumi.lumi_sqrtS = lumi + " fb^{-1} (13 TeV)"
    iPeriod = 0
    iPos = 10
    CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)
    canvas.cd()
    canvas.Update()
    canvas.RedrawAxis()


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

def SetupLegend(x1=0.17,y1=0.7,x2=0.95,y2=0.88, fillStyle=0, borderSize=0, lineWidth=1, NColumns=1, textFont=42, textSize=0.05):
    leg = ROOT.TLegend(x1, y1, x2, y2)
    leg.SetFillStyle(fillStyle)
    leg.SetBorderSize(borderSize)
    leg.SetLineWidth(lineWidth)
    leg.SetNColumns(NColumns)
    leg.SetTextFont(textFont)
    ROOT.gStyle.SetLegendTextSize(textSize)
    return leg

def createCanvasPads(c,isLogY):
	# Upper histogram plot is pad1
    # eps = 0.005
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1.0, 0.97)
    pad1.SetBottomMargin(0.01)  # joins upper and lower plot
    pad1.SetLeftMargin(0.10)
    pad1.SetRightMargin(0.05)
    pad1.SetTopMargin(0.1)
    # pad1.SetBottomMargin(0.12)
    pad1.SetTicks(1,1)
    pad1.SetLogy(isLogY)
    pad1.SetGrid(1,1)
    pad1.Draw()
    # Lower ratio plot is pad2
    c.cd()  # returns to main canvas before defining pad2
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.0, 1.0, 0.3)
    pad2.SetTopMargin(0)  # joins upper and lower plot
    # pad2.SetBottomMargin(0.3)
    pad2.SetLeftMargin(0.10)
    pad2.SetRightMargin(0.05)
    # pad2.SetTopMargin(0.08)
    pad2.SetBottomMargin(0.35)
    pad2.SetTicks(1,1)
    pad2.SetGrid()
    pad2.Draw()

    return c, pad1, pad2

def setupAxes(dummy, xOffset, yOffset, xTitle, yTitle, xLabel, yLabel,title=""):
    dummy.SetStats(0)
    dummy.SetTitle(title)
    dummy.GetXaxis().SetTitleOffset(xOffset)
    dummy.GetYaxis().SetTitleOffset(yOffset)
    dummy.GetXaxis().SetTitleSize(xTitle)
    dummy.GetYaxis().SetTitleSize(yTitle)
    dummy.GetXaxis().SetLabelSize(xLabel)
    dummy.GetYaxis().SetLabelSize(yLabel)
    if(dummy.GetXaxis().GetNdivisions() % 100 > 5): dummy.GetXaxis().SetNdivisions(6,5,0)

def setupDummy(dummy, leg, histName, xAxisLabel, yAxisLabel, isLogY, xmin, xmax, ymin, ymax, lmax, title="", norm=False, normBkg=False,isRatio=False,isABCD=False):
    
    if isRatio:
        setupAxes(dummy, 0, 1.05, 0.0, 0.05, 0.0, 0.05)
        dummy.GetXaxis().SetTitle("")
    elif isABCD:
        setupAxes(dummy, 0.6, 0.8, 0.1, 0.1, 0.08, 0.06,title)
        dummy.GetXaxis().SetTitle(xAxisLabel)
        dummy.GetXaxis().CenterLabels()
        
    else:    
        setupAxes(dummy, 1.2, 1.6, 0.045, 0.045, 0.045, 0.045)
        dummy.GetXaxis().SetTitle(xAxisLabel)
    
    dummy.GetYaxis().SetTitle(yAxisLabel)
    dummy.SetTitle(histName)
    #Set the y-range of the histogram
    if(isLogY):
        if norm:
            default = 0.00001
        else:
            default = 0.02
        locMin = min(default, max(default, 0.05 * ymin))
        legSpan = (math.log10(3*ymax) - math.log10(locMin)) * (leg.GetY1() - ROOT.gPad.GetBottomMargin()) / ((1 - ROOT.gPad.GetTopMargin()) - ROOT.gPad.GetBottomMargin())
        legMin = legSpan + math.log10(locMin)
        if(math.log10(lmax) > legMin):
            scale = (math.log10(lmax) - math.log10(locMin)) / (legMin - math.log10(locMin))
            if norm:
                ymax = 2.
            else:
                ymax = pow(ymax/locMin, scale)*locMin
                # ymax = 10**8
        # dummy.GetYaxis().SetRangeUser(locMin, 10*ymax)
        dummy.GetYaxis().SetRangeUser(locMin, ymax*10e-5)
    else:
        locMin = 0.0
        legMin = (1.2*ymax - locMin) * (leg.GetY1() - ROOT.gPad.GetBottomMargin()) / ((1 - ROOT.gPad.GetTopMargin()) - ROOT.gPad.GetBottomMargin())
        if(lmax > legMin): ymax *= (lmax - locMin)/(legMin - locMin)
        # dummy.GetYaxis().SetRangeUser(0.0, ymax*1.2)
        dummy.GetYaxis().SetRangeUser(0.0, ymax)
    #set x-axis range
    if not isABCD:
        if(xmin < xmax): dummy.GetXaxis().SetRangeUser(xmin, xmax)


def SetupDataStyle(datahist):
    ROOT.gStyle.SetErrorX(0.)
    datahist.SetMarkerStyle(20)
    datahist.SetMarkerSize(1)
    datahist.SetLineColor(ROOT.kBlack)
    return datahist

def SetupLineHistStyle(hist, color=ROOT.kBlue, width=2, style = ROOT.kSolid):
    hist.SetLineColor(color)
    hist.SetLineWidth(width)
    hist.SetLineStyle(style)

def AddVerticalLineAtBinEnd(histogram, bin_number, ymin= 0, ymax = None, color=ROOT.kBlack, line_width=2, line_style = ROOT.kSolid):
    '''
    Adds a vertical line at the end of a specified bin in the given histogram.

    Parameters:
        histogram: The histogram to which the vertical line will be added.
        bin_number: The bin number where the line will be drawn (1-based index).
        color: Color of the vertical line (default: black).
        line_width: Line width of the vertical line (default: 2).
    '''
    if bin_number < 1 or bin_number > histogram.GetNbinsX():
        print("Bin number out of range.")
        return
    
    x_end = histogram.GetBinLowEdge(bin_number + 1)  # End of the bin
    # ymin = histogram.GetMinimum()
    if ymax is None:
        ymax = histogram.GetMaximum()*1.1
    
    # Create a TLine for the vertical line
    line = ROOT.TLine(x_end, ymin, x_end, ymax)
    print(f"In the Line function --- x_end = {x_end}, ymin = {ymin}, ymax={ymax}")
    line.SetLineColor(color)
    line.SetLineWidth(line_width)
    line.SetLineStyle(line_style)
    return line
    # line.Draw("same")

def AddVerticalLineForABCD(histogram, SVJbins, ymin=0 , ymax=None, color=ROOT.kBlack, line_width=2):
    '''
    Adds vertical lines to separate A, B, C, D regions in the given histogram.

    Parameters:
        histogram: The merged histogram for which the lines will be drawn.
        SVJbins: Dictionary containing the SVJ bins information for determining bin boundaries.
        color: Color of the vertical lines (default: black).
        line_width: Line width of the vertical lines (default: 2).
    '''
    # Get bin boundaries from SVJbins
    binwidth = histogram.GetNbinsX()
    region_width = binwidth/4
    bin_boundaries = [(i+1)*region_width for i in range(len(SVJbins))]  # List of bin edges
    print(f"For the Vertical line bin boundaries - {bin_boundaries}, binwidth = {binwidth}, region_width = {region_width}")
    # Draw vertical lines at the end of each bin
    for bin_number in bin_boundaries:
        AddVerticalLineAtBinEnd(histogram, int(bin_number), ymin=ymin, ymax=ymax, color=color, line_width=line_width)

def AddVerticalLine(hist, SVJBins, ymin = 0, ymax = None, color = ROOT.kBlack, line_width=3):
    if ymax is None:
        ymax = hist.GetMaximum()*1.1
    num_bins_per_region = len(SVJBins)
    line = []
    positions = [num_bins_per_region * (i + 1) for i in range(3)]

    for x in positions:
        line.append(ROOT.TLine(x, ymin, x, ymax))
        line[-1].SetLineColor(color)
        line[-1].SetLineWidth(line_width)
        # line[-1].Draw("same")
    return line
        
def AddLabelsForABCD(hist, SVJbins, yloc=None):
    """
    Add A, B, C, D labels to the histogram at specified positions based on SVJbins.
    
    Parameters:
        hist (ROOT.TH1F): The histogram to add labels to.
        SVJbins (dict): Dictionary of SVJ bins.
    """
    num_bins_per_region = len(SVJbins)
    labels = ["A", "B", "C", "D"] # TODO : read labels directy for the regions
    
    # Calculate positions for the labels
    positions = [num_bins_per_region * (i + 0.5) + 1 for i in range(4)]

    # Create TLatex object for adding labels
    latex = ROOT.TLatex()
    latex.SetTextSize(0.06)  # Set text size
    latex.SetTextAlign(22)   # Center text alignment
    if yloc is None:
        yloc = hist.GetMaximum()
    # Add labels to the histogram at the calculated positions
    for label, pos in zip(labels, positions):
        x = hist.GetBinCenter(int(pos))
        latex.DrawLatex(x, yloc * 1.1, label)



def printBinContentAndError(hist):
    # Loop over all bins of the histogram
    for binx in range(1, hist.GetNbinsX() + 1):
        if isinstance(hist, ROOT.TH2):  # For 2D histograms
            for biny in range(1, hist.GetNbinsY() + 1):
                content = hist.GetBinContent(binx, biny)
                error = hist.GetBinError(binx, biny)
                print(f"Bin({binx}, {biny}): Content = {content}, Error = {error}")
        else:  # For 1D histograms
            content = hist.GetBinContent(binx)
            error = hist.GetBinError(binx)
            print(f"Bin({binx}): Content = {content}, Error = {error}")
