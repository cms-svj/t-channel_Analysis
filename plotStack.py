import ROOT
ROOT.gROOT.SetBatch(True)
import math
import utils.DataSetInfo as info
import optparse
import copy
import math
import os

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

def getData(path, scale=1.0, year = "2018"):
    Data = [
        #info.DataSetInfo(basedir=path, fileName=year+"_Data.root",        sys= -1.0, label="Data",        scale=scale),
    ]

    bgData = [
        #info.DataSetInfo(basedir=path, fileName=year+"_Triboson.root",        label="VVV",                     scale=scale, color=(ROOT.kGray)),
        #info.DataSetInfo(basedir=path, fileName=year+"_Diboson.root",         label="VV",                      scale=scale, color=(ROOT.kMagenta + 1)),
        info.DataSetInfo(basedir=path, fileName=year+"_ZJets.root",           label="Z#rightarrow#nu#nu+jets", scale=scale, color=(ROOT.kGray + 1)),
        info.DataSetInfo(basedir=path, fileName=year+"_TTJets.root",          label="t#bar{t}",                scale=scale, color=(ROOT.kBlue - 6)),
        # info.DataSetInfo(basedir=path, fileName=year+"_DYJetsToLL_M-50.root", label="Z#gamma*+jets",           scale=scale, color=(ROOT.kOrange + 2)),
        #info.DataSetInfo(basedir=path, fileName=year+"_TTX.root",             label="ttX",                     scale=scale, color=(ROOT.kCyan + 1)),
        #info.DataSetInfo(basedir=path, fileName=year+"_ST.root",              label="Single top",              scale=scale, color=(ROOT.kRed + 1)),
        info.DataSetInfo(basedir=path, fileName=year+"_WJets.root",           label="W+jets",                  scale=scale, color=(ROOT.kYellow + 1)),
        info.DataSetInfo(basedir=path, fileName=year+"_QCD.root",             label="QCD",                     scale=scale, color=(ROOT.kGreen + 1)),
    ]

    sgData = [
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-200_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",     label="t-ch 200",  scale=scale, color=ROOT.kOrange + 2),
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-400_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",     label="t-ch 400",  scale=scale, color=ROOT.kMagenta + 1),
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",     label="t-ch 600",  scale=scale, color=ROOT.kBlack),
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-800_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",     label="t-ch 800",  scale=scale, color=ROOT.kGreen),
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="t-ch 2000", scale=scale, color=ROOT.kBlue),
        # info.DataSetInfo(basedir=path, fileName=year+"_mZprime-2100_mDark-20_rinv-0p3_alpha-peak.root",          label="s-ch 2100", scale=scale, color=ROOT.kRed),
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="t-ch 3000", scale=scale, color=ROOT.kGreen+2),
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-6000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="t-ch 6000", scale=scale, color=ROOT.kCyan,)
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

def setupDummy(dummy, leg, histName, xAxisLabel, yAxisLabel, isLogY, xmin, xmax, ymin, ymax, lmax):
    setupAxes(dummy, 1.2, 1.6, 0.045, 0.045, 0.045, 0.045)
    dummy.GetYaxis().SetTitle(yAxisLabel)
    dummy.GetXaxis().SetTitle(xAxisLabel)
    dummy.SetTitle(histName)
    #Set the y-range of the histogram
    if(isLogY):
        #default = 0.02
        default = 0.000001
        locMin = min(default, max(default, 0.05 * ymin))
        legSpan = (math.log10(3*ymax) - math.log10(locMin)) * (leg.GetY1() - ROOT.gPad.GetBottomMargin()) / ((1 - ROOT.gPad.GetTopMargin()) - ROOT.gPad.GetBottomMargin())
        legMin = legSpan + math.log10(locMin)
        if(math.log10(lmax) > legMin):
            scale = (math.log10(lmax) - math.log10(locMin)) / (legMin - math.log10(locMin))
            ymax = pow(ymax/locMin, scale)*locMin
        dummy.GetYaxis().SetRangeUser(locMin, 10*ymax)
    else:
        locMin = 0.0
        legMin = (1.2*ymax - locMin) * (leg.GetY1() - ROOT.gPad.GetBottomMargin()) / ((1 - ROOT.gPad.GetTopMargin()) - ROOT.gPad.GetBottomMargin())
        if(lmax > legMin): ymax *= (lmax - locMin)/(legMin - locMin)
        dummy.GetYaxis().SetRangeUser(0.0, ymax*1.2)
    #set x-axis range
    if(xmin < xmax): dummy.GetXaxis().SetRangeUser(xmin, xmax)

#def smartMax():


def plotStack(data, histoName, outputPath="./", xTitle="", yTitle="", isLogY=False, rebinx=-1.0, xmin=999.9, xmax=-999.9):
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
        hs.Add(copy.deepcopy(h))
        leg.AddEntry(h, d.legEntry(), "F")
        if(firstPass):
            hMC = copy.deepcopy(h)
            firstPass = False
        else:
            hMC.Add(copy.deepcopy(h))
    # normHisto(hMC, False)

    #create a dummy histogram to act as the axes
    ymax=10**11
    # ymax=10**1
    ymin=10**-4
    # ymin=10**-12
    lmax=10**12
    # lmax=10**1
    dummy = ROOT.TH1D("dummy", "dummy", 1000, hMC.GetBinLowEdge(1), hMC.GetBinLowEdge(hMC.GetNbinsX()) + hMC.GetBinWidth(hMC.GetNbinsX()))
    setupDummy(dummy, leg, histoName, xTitle, yTitle, isLogY, xmin, xmax, ymin, ymax, lmax)
    dummy.Draw("hist")
    # hMC.Draw("hist same")
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
            # normHisto(h, False)
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

    dummy.Draw("AXIS same")

    c1.SaveAs(outputPath+"/"+histoName+".pdf")
    c1.Close()
    del c1
    del leg
    del hMC

def main():
    parser = optparse.OptionParser("usage: %prog [options]\n")
    parser.add_option('-y', dest='year', type='string', default='2016', help="Can pass in the run year")
    options, args = parser.parse_args()

    #year = options.year
    year = "2018"
    # cuts = ["", "_ge2AK8j", "_ge2AK8j_lp6METrST", "_ge2AK8j_l1p5dEta12", "_baseline"]
    #cuts = ["_ge2AK8j"]
    cuts = [""]
    Data, sgData, bgData = getData("condor/testHadd/", 1.0, year)
    #Data, sgData, bgData = getData("condor/MakeNJetsDists_"+year+"/", 1.0, year)

    plotOutDir = "plots"

    if not os.path.exists(plotOutDir):
        os.makedirs(plotOutDir)

    for cut in cuts:
        plotStack((Data, bgData, sgData), "h_njets"+cut,                plotOutDir, "N_{j}",                           "A.U.", isLogY=True, rebinx=-1, xmin=0, xmax=20)
        plotStack((Data, bgData, sgData), "h_njetsAK8"+cut,             plotOutDir, "N_{J}",                           "A.U.", isLogY=True, rebinx=-1, xmin=0, xmax=12)
        # plotStack((Data, bgData, sgData), "h_ntops"+cut,       "./", "N_{t}",                           "A.U.", isLogY=True, rebinx=-1, xmin=0, xmax=6)
        # plotStack((Data, bgData, sgData), "h_nb"+cut,          "./", "N_{b}",                           "A.U.", isLogY=True, rebinx=-1)
        #plotStack((Data, bgData, sgData), "h_nl"+cut,          "./", "N_{lep}",                         "A.U.", isLogY=True, rebinx=-1)
        #plotStack((Data, bgData, sgData), "h_ne"+cut,          "./", "N_{el}",                          "A.U.", isLogY=True, rebinx=-1)
        #plotStack((Data, bgData, sgData), "h_nm"+cut,          "./", "N_{mu}",                          "A.U.", isLogY=True, rebinx=-1)
        plotStack((Data, bgData, sgData), "h_ht"+cut,                   plotOutDir, "H_{T}",                           "A.U.", isLogY=True, rebinx=20)
        plotStack((Data, bgData, sgData), "h_st"+cut,                   plotOutDir, "S_{T}",                           "A.U.", isLogY=True, rebinx=20)
        plotStack((Data, bgData, sgData), "h_met"+cut,                  plotOutDir, "MET",                             "A.U.", isLogY=True, rebinx=20, xmin=0, xmax=2000)
        plotStack((Data, bgData, sgData), "h_jPt"+cut,                  plotOutDir, "pT_{j}",                          "A.U.", isLogY=True, rebinx=10)
        plotStack((Data, bgData, sgData), "h_jEta"+cut,                 plotOutDir, "#eta_{j}",                        "A.U.", isLogY=True, rebinx=10)
        plotStack((Data, bgData, sgData), "h_jPhi"+cut,                 plotOutDir, "#phi_{j}",                        "A.U.", isLogY=True, rebinx=10)
        plotStack((Data, bgData, sgData), "h_jPtAK8"+cut,               plotOutDir, "pT_{J}",                          "A.U.", isLogY=True, rebinx=10)
        plotStack((Data, bgData, sgData), "h_jEtaAK8"+cut,              plotOutDir, "#eta_{J}",                        "A.U.", isLogY=True, rebinx=10)
        plotStack((Data, bgData, sgData), "h_jPhiAK8"+cut,              plotOutDir, "#phi_{J}",                        "A.U.", isLogY=True, rebinx=10)
        plotStack((Data, bgData, sgData), "h_dEtaJ12"+cut,              plotOutDir, "#Delta#eta(J_{1},J_{2})",         "A.U.", isLogY=True, rebinx=5, xmin=0, xmax=6)
        plotStack((Data, bgData, sgData), "h_dRJ12"+cut,                plotOutDir, "#DeltaR(J_{1},J_{2})",            "A.U.", isLogY=True, rebinx=2, xmin=0, xmax=6)
        plotStack((Data, bgData, sgData), "h_dPhiJ1MET"+cut,            plotOutDir, "#Delta#phi(J_{1},MET)",           "A.U.", isLogY=True, rebinx=2, xmin=0, xmax=3.15)
        plotStack((Data, bgData, sgData), "h_dPhiJ2MET"+cut,            plotOutDir, "#Delta#phi(J_{2},MET)",           "A.U.", isLogY=True, rebinx=2, xmin=0, xmax=3.15)
        plotStack((Data, bgData, sgData), "h_dPhiMinJMET"+cut,          plotOutDir, "#Delta#phi_{min}(J_{1,2},MET)",   "A.U.", isLogY=True, rebinx=2, xmin=0, xmax=3.15)
        plotStack((Data, bgData, sgData), "h_dPhiJ1METrdPhiJ2MET"+cut,  plotOutDir, "#Delta#phi Ratio",        "A.U.", isLogY=True, rebinx=2)
        # plotStack((Data, bgData, sgData), "h_mjjM"+cut,        "./", "m_{JJ}",                          "A.U.", isLogY=True, rebinx=10)
        # plotStack((Data, bgData, sgData), "h_mjjPt"+cut,       "./", "pT_{JJ}",                         "A.U.", isLogY=True, rebinx=10)
        # plotStack((Data, bgData, sgData), "h_mjjEta"+cut,      "./", "#eta_{JJ}",                       "A.U.", isLogY=True, rebinx=10)
        # plotStack((Data, bgData, sgData), "h_mT"+cut,          "./", "mT",                              "A.U.", isLogY=True, rebinx=20)
        plotStack((Data, bgData, sgData), "h_METrHT_pt30"+cut,          plotOutDir, "MET/HT",                          "A.U.", isLogY=True, rebinx=2, xmin=0, xmax=10)
        plotStack((Data, bgData, sgData), "h_METrST_pt30"+cut,          plotOutDir, "MET/ST",                          "A.U.", isLogY=True, rebinx=5)

if __name__ == '__main__':
    main()
