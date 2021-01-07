import ROOT
ROOT.gROOT.SetBatch(True)
import math
import utils.DataSetInfo as info
import optparse
import copy
import math
import os
from array import array
import numpy as np
import utils.CMS_lumi as CMS_lumi

def normHisto(hist, doNorm=False):
    if doNorm:
        hist.Scale(1.0/hist.Integral())

def simpleSig(hSig, hBg):
    sig = 0.0
    for i in range(0, hSig.GetNbinsX()):
        totBG = hBg.GetBinContent(i)
        nSig = hSig.GetBinContent(i)
        if(totBG > 1.0 and nSig > 1.0):
            s = nSig / math.sqrt( totBG + (0.3*totBG)**2 )
            sig = math.sqrt(sig**2 + s**2)
    return sig

def getBGHistos(data, histoName, rebinx, xmin, xmax):
    hs = ROOT.THStack()
    hMC = None
    hList = []
    firstPass = True
    for d in data[1]:
        h = d.getHisto(histoName, rebinx=rebinx, xmin=xmin, xmax=xmax, fill=True, showEvents=True)
        hist = copy.deepcopy(h)
        hs.Add(hist)
        hList.append((hist, d.legEntry()))
        if(firstPass):
            hMC = hist
            firstPass = False
        else:
            hMC.Add(hist)
    return hs, hMC, hList

def getData(path, scale=1.0, year = "2018"):
    Data = [
        #info.DataSetInfo(basedir=path, fileName=year+"_Data.root",        sys= -1.0, label="Data",        scale=scale),
    ]

    # Normal
    bgData = [
        # info.DataSetInfo(basedir=path, fileName=year+"_Triboson.root",        label="VVV",                     scale=scale, color=(ROOT.kGray)),
        # info.DataSetInfo(basedir=path, fileName=year+"_Diboson.root",         label="VV",                      scale=scale, color=(ROOT.kMagenta + 1)),
        # info.DataSetInfo(basedir=path, fileName=year+"_DYJetsToLL_M-50.root", label="Z#gamma*+jets",           scale=scale, color=(ROOT.kOrange + 2)),
        # info.DataSetInfo(basedir=path, fileName=year+"_TTX.root",             label="ttX",                     scale=scale, color=(ROOT.kCyan + 1)),
        # info.DataSetInfo(basedir=path, fileName=year+"_ST.root",              label="Single top",              scale=scale, color=(ROOT.kRed + 1)),
        info.DataSetInfo(basedir=path, fileName=year+"_ZJets.root",           label="Z#rightarrow#nu#nu+jets", scale=scale, color=(ROOT.kGray + 1)),
        info.DataSetInfo(basedir=path, fileName=year+"_TTJets.root",          label="t#bar{t}",                scale=scale, color=(ROOT.kBlue - 6)),
        info.DataSetInfo(basedir=path, fileName=year+"_WJets.root",           label="W+jets",                  scale=scale, color=(ROOT.kYellow + 1)),
        info.DataSetInfo(basedir=path, fileName=year+"_QCD.root",             label="QCD",                     scale=scale, color=(ROOT.kGreen + 1)),
    ]
    #
    sgData = [
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="t-ch 3000", scale=scale, color=ROOT.kMagenta+1),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",     label="t-ch 600",  scale=scale, color=ROOT.kBlack),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-800_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",     label="t-ch 800",  scale=scale, color=ROOT.kGreen),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="t-ch 2000", scale=scale, color=ROOT.kBlue),
        # info.DataSetInfo(basedir=path, fileName="2017_mZprime-3000_mDark-20_rinv-0p3_alpha-peak.root",           label="s-ch baseline", scale=scale, color=ROOT.kRed),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-6000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="t-ch 6000", scale=scale, color=ROOT.kCyan,)
        ## varying mMed
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-400_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",     label="mMed 400",  scale=scale, color=ROOT.kMagenta + 1),
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-800_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",     label="mMed 800",  scale=scale, color=ROOT.kOrange + 2),
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="mMed 2000", scale=scale, color=ROOT.kBlue),
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="mMed 3000", scale=scale, color=ROOT.kGreen+2),
        info.DataSetInfo(basedir=path, fileName="2017_mZprime-3000_mDark-20_rinv-0p3_alpha-peak.root",          label="s-ch 2100", scale=scale, color=ROOT.kRed),
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-6000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="mMed 6000", scale=scale, color=ROOT.kCyan),
        ## varying mDark
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-3000_mDark-1_rinv-0p3_alpha-peak_yukawa-1.root",    label="M-3000_mD-1",scale=scale, color=ROOT.kBlue),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="M-3000_mD-20",scale=scale, color=ROOT.kGreen+2),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-3000_mDark-50_rinv-0p3_alpha-peak_yukawa-1.root",    label="M-3000_mD-50",scale=scale, color=ROOT.kRed),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-3000_mDark-100_rinv-0p3_alpha-peak_yukawa-1.root",    label="M-3000_mD-100",scale=scale, color=ROOT.kCyan),
        ## varying rinv
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-3000_mDark-20_rinv-0p1_alpha-peak_yukawa-1.root",    label="M-3000_r-0p1",scale=scale, color=ROOT.kBlue),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="M-3000_r-0p3",scale=scale, color=ROOT.kGreen+2),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-3000_mDark-20_rinv-0p5_alpha-peak_yukawa-1.root",    label="M-3000_r-0p5",scale=scale, color=ROOT.kRed),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-3000_mDark-20_rinv-0p7_alpha-peak_yukawa-1.root",    label="M-3000_r-0p7",scale=scale, color=ROOT.kCyan),
        ## varying rinv at mMed 800
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-800_mDark-20_rinv-0p1_alpha-peak_yukawa-1.root",     label="M-800_r-0p1", scale=scale, color=ROOT.kOrange + 2),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-800_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",     label="M-800_r-0p3", scale=scale, color=ROOT.kMagenta + 1),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-800_mDark-20_rinv-0p5_alpha-peak_yukawa-1.root",     label="M-800_r-0p5", scale=scale, color=ROOT.kBlack),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-800_mDark-20_rinv-0p7_alpha-peak_yukawa-1.root",     label="M-800_r-0p7", scale=scale, color=ROOT.kGreen),
    ]


    return Data, sgData, bgData

def setupAxes(dummy, xOffset, yOffset, xTitle, yTitle, xLabel, yLabel):
    dummy.SetStats(0)
    dummy.SetTitle("")
    dummy.GetXaxis().SetTitleOffset(xOffset)
    dummy.GetYaxis().SetTitleOffset(yOffset)
    dummy.GetXaxis().SetTitleSize(xTitle)
    dummy.GetYaxis().SetTitleSize(yTitle)
    dummy.GetXaxis().SetLabelSize(xLabel)
    dummy.GetYaxis().SetLabelSize(yLabel)
    if(dummy.GetXaxis().GetNdivisions() % 100 > 5): dummy.GetXaxis().SetNdivisions(6, 5, 0)

def setupDummy(dummy, leg, histName, xAxisLabel, yAxisLabel, isLogY, xmin, xmax, ymin, ymax, lmax, norm=False, normBkg=False):
    setupAxes(dummy, 1.2, 1.6, 0.045, 0.045, 0.045, 0.045)
    dummy.GetYaxis().SetTitle(yAxisLabel)
    dummy.GetXaxis().SetTitle(xAxisLabel)
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
        dummy.GetYaxis().SetRangeUser(locMin, 10*ymax)
    else:
        locMin = 0.0
        legMin = (1.2*ymax - locMin) * (leg.GetY1() - ROOT.gPad.GetBottomMargin()) / ((1 - ROOT.gPad.GetTopMargin()) - ROOT.gPad.GetBottomMargin())
        if(lmax > legMin): ymax *= (lmax - locMin)/(legMin - locMin)
        dummy.GetYaxis().SetRangeUser(0.0, ymax*1.2)
    #set x-axis range
    if(xmin < xmax): dummy.GetXaxis().SetRangeUser(xmin, xmax)

def makeRocVec(h):
    h.Scale( 1.0 / h.Integral() );
    v, cuts = [], []
    for i in range(0, h.GetNbinsX()+1):
        val = h.Integral(i, h.GetNbinsX())
        # if round(val,2) > 1.0:
        #     print (i)
        #     print (h.GetNbinsX())
        #     print (h.Integral())
        #     print (h.Integral(i, h.GetNbinsX()))
        v.append(val)
        cuts.append(h.GetBinLowEdge(i)+h.GetBinWidth(i))
    return v, cuts

def drawRocCurve(fType, rocBgVec, rocSigVec, leg, rebinx,manySigs=False):
    h = []

    # deciding whether to flip the ROC curves based on the baseline vs. QCD curve
    flip = False
    baselineNames = ["t-ch 3000","mMed 3000","M-800_r-0p3","M-3000_mD-20"]
    for mBg, cutBg, lBg, cBg in rocBgVec:
        for mSig, cutSig, lSig, cSig in rocSigVec:
            n = len(mBg)
            if ("QCD" in lBg) and (lSig in baselineNames):
                mBgAr = [1] + mBg + [0]
                mSigAr = [0] + mSig + [0]
                gAr = ROOT.TGraph(n, array("d", mBgAr), array("d", mSigAr))
                gArea = gAr.Integral()
                print (gArea)
                if gArea < 0.5:
                    flip = True
                    break

    if manySigs:
        for rbv in rocBgVec:
            if rbv[2] == "QCD":
                QCDVec = rbv
        rocBgVec = [QCDVec]
    else:
        for rsv in rocSigVec:
            if rsv[2] in baselineNames:
                baseVec = rsv
        rocSigVec = [baseVec]

    for mBg, cutBg, lBg, cBg in rocBgVec:
        for mSig, cutSig, lSig, cSig in rocSigVec:
            flip = True
            n = len(mBg)
            rv = ">cut"
            if flip:
                mBg_f = 1 - np.array(mBg)
                mSig_f = 1 - np.array(mSig)
                rv = "<cut"

            if manySigs:
                col = cSig
            else:
                col = cBg

            g = ROOT.TGraph(n, array("d", mBg_f), array("d", mSig_f))
            for i in range(0,n):
                if ((i % rebinx == 0) and (rebinx != -1)):
                    latex = ROOT.TLatex(g.GetX()[i], g.GetY()[i],str(round(cutSig[i],2)))
                    latex.SetTextSize(0.02)
                    latex.SetTextColor(ROOT.kRed)
                    g.GetListOfFunctions().Add(latex) # add cut values
            g.SetLineWidth(2)
            g.SetLineColor(col)
            g.SetMarkerSize(0.7)
            g.SetMarkerStyle(ROOT.kFullSquare)
            g.SetMarkerColor(col)
            g.Draw("same LP text")
            leg.AddEntry(g, fType + " " + lBg + " vs " + lSig + "_" + rv, "LP")
            h.append(g)
    return h

def plotROC(data, histoName, outputPath="./", isLogY=False, rebinx=-1.0, xmin=999.9, xmax=-999.9, norm=False, manySigs=False):
    #This is a magic incantation to disassociate opened histograms from their files so the files can be closed
    ROOT.TH1.AddDirectory(False)

    #create the canvas for the plot
    c1 = ROOT.TCanvas( "c", "c", 800, 800)
    c1.cd()
    ROOT.gPad.Clear()
    ROOT.gStyle.SetOptStat("")
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.08)
    ROOT.gPad.SetBottomMargin(0.12)
    ROOT.gPad.SetTicks(1,1)
    ROOT.gPad.SetLogy(isLogY)

    #Create TLegend
    leg = ROOT.TLegend(0.17, 0.72, 0.95, 0.88)
    #nColumns = 3 if(len(data[1]) >= 3) else 1
    nColumns = 2
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetLineWidth(1)
    leg.SetNColumns(nColumns)
    leg.SetTextFont(42)
    ROOT.gStyle.SetLegendTextSize(0.024)

    rocBgVec = []
    for d in data[1]:
        h = d.getHisto(histoName, rebinx=-1, xmin=xmin, xmax=xmax, fill=True, showEvents=False)
        rocBgVec.append(makeRocVec(h) + ( d.legEntry(), d.getColor()))

    rocSigVec = []
    for d in data[2]:
        h = d.getHisto(histoName, rebinx=-1, xmin=xmin, xmax=xmax, fill=True, showEvents=False)
        rocSigVec.append(makeRocVec(h) + (d.legEntry(), d.getColor()))

    #create a dummy histogram to act as the axes
    ymax=1.0
    ymin=10**-4
    lmax=1.0
    dummy = ROOT.TH1D("dummy", "dummy", 1000, 0.0, 1.0)
    setupDummy(dummy, leg, "", "#epsilon_{ bg}", "#epsilon_{ sg}", isLogY, xmin, xmax, ymin, ymax, lmax)
    dummy.Draw("hist")
    leg.Draw("same")
    history = drawRocCurve("", rocBgVec, rocSigVec, leg, rebinx, manySigs)

    line1 = ROOT.TF1( "line1","1",0,1)
    line1.SetLineColor(ROOT.kBlack)
    line1.Draw("same")
    line2 = ROOT.TF1( "line2","x",0,1)
    line2.SetLineColor(ROOT.kBlack)
    line2.SetLineStyle(ROOT.kDotted)
    line2.Draw("same")

    dummy.Draw("AXIS same")

    # CMS label
    CMS_lumi.writeExtraText = 1
    lumi = "59.7"

    CMS_lumi.lumi_sqrtS = lumi + " fb^{-1} (13 TeV)"

    iPeriod = 0
    iPos = 0

    CMS_lumi.CMS_lumi(c1, iPeriod, iPos)
    c1.cd()
    c1.Update();
    c1.RedrawAxis()

    c1.SaveAs(outputPath+"/roc/"+histoName+"_ROC.png")
    c1.Close()
    del c1
    del leg

def plotStack(data, histoName, outputPath="./", xTitle="", yTitle="", isLogY=False, rebinx=-1.0, xmin=999.9, xmax=-999.9, norm=False, normBkg=False, onlySig=False):
    #This is a magic incantation to disassociate opened histograms from their files so the files can be closed
    ROOT.TH1.AddDirectory(False)

    #create the canvas for the plot
    c1 = ROOT.TCanvas( "c", "c", 800, 800)
    c1.cd()
    ROOT.gPad.Clear()
    ROOT.gStyle.SetOptStat("")
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.08)
    ROOT.gPad.SetBottomMargin(0.12)
    ROOT.gPad.SetTicks(1,1)
    ROOT.gPad.SetLogy(isLogY)

    #Create TLegend
    leg = ROOT.TLegend(0.17, 0.72, 0.95, 0.88)
    #nColumns = 3 if(len(data[1]) >= 3) else 1
    nColumns = 2
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetLineWidth(1)
    leg.SetNColumns(nColumns)
    leg.SetTextFont(42)
    ROOT.gStyle.SetLegendTextSize(0.024)

    #Setup background histos
    hs = ROOT.THStack()
    hMC = None
    firstPass = True
    for d in data[1]:
        h = d.getHisto(histoName, rebinx=rebinx, xmin=xmin, xmax=xmax, fill=True, showEvents=True)
        if normBkg:
            normHisto(h, True)
            h.SetLineWidth(3)
            h.SetFillStyle(3955)
        hs.Add(copy.deepcopy(h))
        leg.AddEntry(h, d.legEntry(), "F")
        if(firstPass):
            hMC = copy.deepcopy(h)
            firstPass = False
        else:
            hMC.Add(copy.deepcopy(h))
    # there is a bug with getBGHistos. Once fixed, can delete lines 294-305, and uncomment
    # the line below and lines 313-314
    # hs, hMC, hList = getBGHistos(data, histoName, rebinx, xmin, xmax)
    if norm:
        normHisto(hMC, True)
    #Fill background legend
    # for h in hList:
    #     leg.AddEntry(h[0], h[1], "F")

    #create a dummy histogram to act as the axes
    if norm:
        ymax=10**1
        ymin=10**-12
        lmax=10**1
    else:
        ymax=10**11
        ymin=10**-4
        lmax=10**12
    dummy = ROOT.TH1D("dummy", "dummy", 1000, hMC.GetBinLowEdge(1), hMC.GetBinLowEdge(hMC.GetNbinsX()) + hMC.GetBinWidth(hMC.GetNbinsX()))
    setupDummy(dummy, leg, "", xTitle, yTitle, isLogY, xmin, xmax, ymin, ymax, lmax, norm, normBkg)
    # setupDummy(dummy, leg, histoName, xTitle, yTitle, isLogY, xmin, xmax, ymin, ymax, lmax, norm, normBkg)
    if normBkg:
        dummy.SetMaximum(100)
        dummy.SetMinimum(0.00001)
    dummy.Draw("hist")
    if norm:
        hMC.Draw("hist same")
        leg.Clear()
        leg.AddEntry(hMC, "Total Background", "L")
    elif normBkg:
        hs.Draw("nostackHIST same")
        hs.SetMaximum(100)
        hs.SetMinimum(0.00001)
    elif onlySig:
        leg.Clear()
    else:
        hs.Draw("hist F same")
    leg.Draw("same")

    #Setup signal histos
    history = []
    sig = 0.0
    if(data[2]):
        #firstPass=True
        for d in data[2]:
            h = d.getHisto(histoName, rebinx=rebinx, xmin=xmin, xmax=xmax, showEvents=True)
            #if(firstPass):
            sig = round(simpleSig(h, hMC),2)
            #firstPass=False
            #print(d.legEntry(), round(simpleSig(h, hMC),2))
            h.SetLineStyle(ROOT.kDashed)
            h.SetLineWidth(3)
            leg.AddEntry(h, d.legEntry()+", {}".format(sig), "L")
            if norm or normBkg:
                normHisto(h, True)
            h.Draw("hist same")
            history.append(h)

    #Draw significance
    significance = ROOT.TLatex()
    significance.SetNDC(True)
    significance.SetTextAlign(11)
    significance.SetTextFont(52)
    significance.SetTextSize(0.030)
    #significance.DrawLatex(0.45, 0.72, ("Significance = #frac{N_{s}}{#sqrt{N_{b}+#left(0.3N_{b}#right)^{2}}} = "+str(sig)))
    #significance.DrawLatex(0.45, 0.72, ("Significance = #frac{N_{s}}{#sqrt{N_{b}+#left(0.3N_{b}#right)^{2}}}"))

    if onlySig:
        dummy.SetMaximum(10**8)
    dummy.Draw("AXIS same")

    # CMS label
    CMS_lumi.writeExtraText = 1
    lumi = "59.7"

    CMS_lumi.lumi_sqrtS = lumi + " fb^{-1} (13 TeV)"

    iPeriod = 0
    iPos = 0

    CMS_lumi.CMS_lumi(c1, iPeriod, iPos)
    c1.cd()
    c1.Update();
    c1.RedrawAxis()

    if norm:
        c1.SaveAs(outputPath+"/"+histoName+"_norm.png")
    elif normBkg:
        c1.SaveAs(outputPath+"/"+histoName+"_normBkg.png")
    else:
        c1.SaveAs(outputPath+"/"+histoName+".png")

    c1.Close()
    del c1
    del leg
    del hMC

def main():
    parser = optparse.OptionParser("usage: %prog [options]\n")
    parser.add_option('-y',                 dest='year',       type='string',  default='2018',                 help="Can pass in the run year")
    parser.add_option('-d', '--dataset',    dest='dataset',                    default='testHadd_11242020',    help='dataset')
    parser.add_option('-n',                 dest='isNorm',     action="store_true",                            help="Normalize stack plots")
    parser.add_option('-b',                 dest='isNormBkg',  action="store_true",                            help="Normalized Bakground and Signal plots")
    parser.add_option('-s',                 dest='onlySig',    action="store_true",                            help="Plot only signals")
    parser.add_option('-m',                 dest='manySigs',   action="store_true",                            help="Plot ROC curves with many signals vs. QCD")
    options, args = parser.parse_args()

    year = options.year
    # cuts = ["", "_ge2AK8j", "_ge2AK8j_lp6METrST", "_ge2AK8j_l1p5dEta12", "_baseline"]
    #cuts = ["_ge2AK8j"]
    cuts = [""]
    Data, sgData, bgData = getData("condor/" + options.dataset + "/", 1.0, year)
    #Data, sgData, bgData = getData("condor/MakeNJetsDists_"+year+"/", 1.0, year)

    plotOutDir = "plots"

    if not os.path.exists(plotOutDir):
        os.makedirs(plotOutDir)

    plotDict = {
    # "plotname":               [xlabel,                                            ylabel,     xmin,   xmax,   rebinx_stack,   rebinx_roc, cuts]
    "h_njets":                  ["Number of AK4 Jets",                              "Events",   0,      20,     -1,             1,          ["","_ge2AK8j","_ge2AK4j"]],
    "h_njetsAK8":               ["Number of AK8 Jets",                              "Events",   0,      12,     -1,             1,          ["","_ge2AK8j","_ge2AK4j"]],
    "h_nb":                     ["Number of B Jets",                                "Events",   0,      20,     -1,             1,          [""]],
    "h_nl":                     ["Number of Leptons",                               "Events",   0,      10,     -1,             1,          [""]],
    "h_ht":                     ["H_{T} [GeV]",                                     "Events",   0,      5000,   20,             10,         ["","_ge2AK8j","_ge2AK4j"]],
    "h_st":                     ["S_{T} [GeV]",                                     "Events",   0,      5000,   20,             10,         ["","_ge2AK8j","_ge2AK4j"]],
    "h_met":                    ["MET [GeV]",                                       "Events",   0,      2000,   20,             10,          ["","_ge2AK8j","_ge2AK4j"]],
    "h_jPt":                    ["p_{T}(j) [GeV]",                                  "Events",   0,      2000,   10,             5,          [""]],
    "h_jEta":                   ["#eta(j)",                                         "Events",   -6,     6,      10,             10,         [""]],
    "h_jPhi":                   ["#phi(j)",                                         "Events",   -4,     4,      10,             10,         [""]],
    "h_jAxismajor":             ["#sigma_{major}(j)",                               "Events",   0,      0.5,    -1,             5,          [""]],
    "h_jAxisminor":             ["#sigma_{minor}(j)",                               "Events",   0,      0.3,    -1,             5,          [""]],
    "h_jPtD":                   ["ptD(j)",                                          "Events",   0,      1.2,    -1,             2,          [""]],
    "h_jPtAK8":                 ["p_{T}(J) [GeV]",                                  "Events",   0,      2000,   10,             5,          [""]],
    "h_jEtaAK8":                ["#eta(J)",                                         "Events",   -6,     6,      10,             10,         [""]],
    "h_jPhiAK8":                ["#phi(J)",                                         "Events",   -4,     4,      10,             10,         [""]],
    "h_jAxismajorAK8":          ["#sigma_{major}(J)",                               "Events",   0,      0.5,    -1,             5,          [""]],
    "h_jAxisminorAK8":          ["#sigma_{minor}(J)",                               "Events",   0,      0.3,    -1,             5,          [""]],
    "h_jGirthAK8":              ["girth(J)",                                        "Events",   0,      0.5,    -1,             5,          [""]],
    "h_jPtDAK8":                ["ptD(J)",                                          "Events",   0,      1.2,    -1,             2,          [""]],
    "h_jTau1AK8":               ["#tau_{1}(J)",                                     "Events",   0,      0.8,    -1,             3,          [""]],
    "h_jTau2AK8":               ["#tau_{2}(J)",                                     "Events",   0,      0.65,   -1,             3,          [""]],
    "h_jTau3AK8":               ["#tau_{3}(J)",                                     "Events",   0,      0.35,   -1,             3,          [""]],
    # "h_jTau21AK8":              ["#tau_{2}/#tau_{1}(J)",                            "Events",   0,      1.3,    -1,             3,          [""]],
    # "h_jTau32AK8":              ["#tau_{3}/#tau_{2}(J)",                            "Events",   0,      1.3,    -1,             3,          [""]],
    "h_jSoftDropMassAK8":       ["m_{SD}(J)",                                       "Events",   0,      200,    -1,             3,          [""]],
    "h_mT":                     ["m_{T}",                                           "Events",   0,      5000,   20,             20,         ["","_ge2AK8j","_ge2AK4j"]],
    "h_METrHT_pt30":            ["MET/H_{T}",                                       "Events",   0,      10,     2,              2,          ["","_ge2AK8j","_ge2AK4j"]],
    "h_METrST_pt30":            ["MET/S_{T}",                                       "Events",   0,      1,      5,              5,          ["","_ge2AK8j","_ge2AK4j"]],
    "h_dEtaJ12":                ["#Delta#eta(J_{1},J_{2})",                         "Events",   0,      6,      5,              5,          ["_ge2AK8j"]],
    "h_dRJ12":                  ["#Delta R(J_{1},J_{2})",                           "Events",   0,      6,      2,              3,          ["_ge2AK8j"]],
    "h_dPhiJ1MET":              ["#Delta#phi(J_{1},MET)",                           "Events",   0,      3.15,   2,              5,          ["_ge2AK8j"]],
    "h_dPhiJ2MET":              ["#Delta#phi(J_{2},MET)",                           "Events",   0,      3.15,   2,              5,          ["_ge2AK8j"]],
    "h_dPhiJ1METrdPhiJ2MET":    ["#Delta#phi(J_{1},MET)/#Delta#phi(J_{2},MET)",     "Events",   0,      100,    2,              2,          ["_ge2AK8j"]],
    "h_dPhiMinJMET":            ["#Delta#phi_{min}(J,MET)",                         "Events",   0,      3.15,   2,              4,          ["_ge2AK8j"]],
    "h_j1PtAK8":                ["p_{T}(J_{1}) [GeV]",                              "Events",   0,      2000,   10,             5,          ["_ge2AK8j"]],
    "h_j1EtaAK8":               ["#eta(J_{1})",                                     "Events",   -6,     6,      10,             10,         ["_ge2AK8j"]],
    "h_j1PhiAK8":               ["#phi(J_{1})",                                     "Events",   -4,     4,      10,             10,         ["_ge2AK8j"]],
    "h_j1AxismajorAK8":         ["#sigma_{major}(J_{1})",                           "Events",   0,      0.5,    -1,             5,          ["_ge2AK8j"]],
    "h_j1AxisminorAK8":         ["#sigma_{minor}(J_{1})",                           "Events",   0,      0.3,    -1,             5,          ["_ge2AK8j"]],
    "h_j1GirthAK8":             ["girth (J_{1})",                                   "Events",   0,      0.5,    -1,             5,          ["_ge2AK8j"]],
    "h_j1PtDAK8":               ["ptD (J_{1})",                                     "Events",   0,      1.2,    -1,             2,          ["_ge2AK8j"]],
    "h_j1Tau1AK8":              ["#tau_{1}(J_{1})",                                 "Events",   0,      0.8,    -1,             3,          ["_ge2AK8j"]],
    "h_j1Tau2AK8":              ["#tau_{2}(J_{1})",                                 "Events",   0,      0.65,   -1,             3,          ["_ge2AK8j"]],
    "h_j1Tau3AK8":              ["#tau_{3}(J_{1})",                                 "Events",   0,      0.35,   -1,             3,          ["_ge2AK8j"]],
    "h_j1Tau21AK8":             ["#tau_{21}(J_{1})",                                "Events",   0,      1.3,    -1,             3,          ["_ge2AK8j"]],
    "h_j1Tau32AK8":             ["#tau_{32}(J_{1})",                                "Events",   0,      1.3,    -1,             3,          ["_ge2AK8j"]],
    "h_j1SoftDropMassAK8":      ["m_{SD}(J_{1})",                                   "Events",   0,      200,    -1,             3,          ["_ge2AK8j"]],
    "h_j2PtAK8":                ["p_{T}(J_{2}) [GeV]",                              "Events",   0,      2000,   10,             5,          ["_ge2AK8j"]],
    "h_j2EtaAK8":               ["#eta(J_{2})",                                     "Events",   -6,     6,      10,             10,         ["_ge2AK8j"]],
    "h_j2PhiAK8":               ["#phi(J_{2})",                                     "Events",   -4,     4,      10,             10,         ["_ge2AK8j"]],
    "h_j2AxismajorAK8":         ["#sigma_{major}(J_{2})",                           "Events",   0,      0.5,    -1,             5,          ["_ge2AK8j"]],
    "h_j2AxisminorAK8":         ["#sigma_{minor}(J_{2})",                           "Events",   0,      0.3,    -1,             5,          ["_ge2AK8j"]],
    "h_j2GirthAK8":             ["girth (J_{2})",                                   "Events",   0,      0.5,    -1,             5,          ["_ge2AK8j"]],
    "h_j2PtDAK8":               ["ptD (J_{2})",                                     "Events",   0,      1.2,    -1,             2,          ["_ge2AK8j"]],
    "h_j2Tau1AK8":              ["#tau_{1}(J_{2})",                                 "Events",   0,      0.8,    -1,             3,          ["_ge2AK8j"]],
    "h_j2Tau2AK8":              ["#tau_{2}(J_{2})",                                 "Events",   0,      0.65,   -1,             3,          ["_ge2AK8j"]],
    "h_j2Tau3AK8":              ["#tau_{3}(J_{2})",                                 "Events",   0,      0.35,   -1,             3,          ["_ge2AK8j"]],
    "h_j2Tau21AK8":             ["#tau_{21}(J_{2})",                                "Events",   0,      1.3,    -1,             3,          ["_ge2AK8j"]],
    "h_j2Tau32AK8":             ["#tau_{32}(J_{2})",                                "Events",   0,      1.3,    -1,             3,          ["_ge2AK8j"]],
    "h_j2SoftDropMassAK8":      ["m_{SD}(J_{2})",                                   "Events",   0,      200,    -1,             3,          ["_ge2AK8j"]],
    "h_dEtaj12":                ["#Delta#eta(j_{1},j_{2})",                         "Events",   0,      6,      5,              5,          ["_ge2AK4j"]],
    "h_dRj12":                  ["#Delta R(j_{1},j_{2})",                           "Events",   0,      6,      2,              3,          ["_ge2AK4j"]],
    "h_dPhij1MET":              ["#Delta#phi(j_{1},MET)",                           "Events",   0,      3.15,   2,              5,          ["_ge2AK4j"]],
    "h_dPhij2MET":              ["#Delta#phi(j_{2},MET)",                           "Events",   0,      3.15,   2,              5,          ["_ge2AK4j"]],
    "h_dPhij1METrdPhij2MET":    ["#Delta#phi(j_{1},MET)/#Delta#phi(j_{2},MET)",     "Events",   0,      100,    2,              2,          ["_ge2AK4j"]],
    "h_dPhiMinjMET":            ["#Delta#phi_{min}(j,MET)",                         "Events",   0,      3.15,   2,              4,          ["_ge2AK4j"]],
    "h_j1Pt":                   ["p_{T}(j_{1}) [GeV]",                              "Events",   0,      2000,   10,             5,          ["_ge2AK4j"]],
    "h_j1Eta":                  ["#eta(j_{1})",                                     "Events",   -6,     6,      10,             10,         ["_ge2AK4j"]],
    "h_j1Phi":                  ["#phi(j_{1})",                                     "Events",   -4,     4,      10,             10,         ["_ge2AK4j"]],
    "h_j1Axismajor":            ["#sigma_{major}(j_{1})",                           "Events",   0,      0.5,    -1,             5,          ["_ge2AK4j"]],
    "h_j1Axisminor":            ["#sigma_{minor}(j_{1})",                           "Events",   0,      0.3,    -1,             5,          ["_ge2AK4j"]],
    "h_j1PtD":                  ["ptD (j_{1})",                                     "Events",   0,      1.2,    -1,             2,          ["_ge2AK4j"]],
    "h_j2Pt":                   ["p_{T}(j_{2}) [GeV]",                              "Events",   0,      2000,   10,             5,          ["_ge2AK4j"]],
    "h_j2Eta":                  ["#eta(j_{2})",                                     "Events",   -6,     6,      10,             10,         ["_ge2AK4j"]],
    "h_j2Phi":                  ["#phi(j_{2})",                                     "Events",   -4,     4,      10,             10,         ["_ge2AK4j"]],
    "h_j2Axismajor":            ["#sigma_{major}(j_{2})",                           "Events",   0,      0.5,    -1,             5,          ["_ge2AK4j"]],
    "h_j2Axisminor":            ["#sigma_{minor}(j_{2})",                           "Events",   0,      0.3,    -1,             5,          ["_ge2AK4j"]],
    "h_j2PtD":                  ["ptD (j_{2})",                                     "Events",   0,      1.2,    -1,             2,          ["_ge2AK4j"]],
    }

    for histName,details in plotDict.items():
        isNorm = options.isNorm
        isNormBkg = options.isNormBkg
        onlySig = options.onlySig
        manySigs = options.manySigs

        for cut in details[6]:
            plotROC(  (Data, bgData, sgData), histName+cut, plotOutDir,                         isLogY=False,   rebinx=details[5], manySigs=manySigs)
            # plotStack((Data, bgData, sgData), histName+cut, plotOutDir, details[0], details[1], isLogY=True,    rebinx=details[4], norm=isNorm, xmin=details[2], xmax=details[3], normBkg=isNormBkg, onlySig=onlySig)

if __name__ == '__main__':
    main()
