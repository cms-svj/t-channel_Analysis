import ROOT
ROOT.gROOT.SetBatch(True)
import math
import DataSetInfo as info
import optparse
import copy
import math

def getData(path, scale=1.0, year = "2016"):
    Data = [
        #info.DataSetInfo(basedir=path, fileName=year+"_Data.root",        sys= -1.0, label="Data",        scale=scale),
    ]
    
    bgData = [
        info.DataSetInfo(basedir=path, fileName=year+"_Triboson.root",        label="Triboson",        scale=scale, color=(ROOT.kGray)),
        info.DataSetInfo(basedir=path, fileName=year+"_Diboson.root",         label="Diboson",         scale=scale, color=(ROOT.kMagenta + 1)),
        info.DataSetInfo(basedir=path, fileName=year+"_DYJetsToLL_M-50.root", label="DYJetsToLL_M-50", scale=scale, color=(ROOT.kOrange + 2)),        
        info.DataSetInfo(basedir=path, fileName=year+"_TTX.root",             label="TTX",             scale=scale, color=(ROOT.kCyan + 1)),
        info.DataSetInfo(basedir=path, fileName=year+"_WJets.root",           label="WJets",           scale=scale, color=(ROOT.kYellow + 1)),
        info.DataSetInfo(basedir=path, fileName=year+"_ST.root",              label="ST",              scale=scale, color=(ROOT.kRed + 1)),
        info.DataSetInfo(basedir=path, fileName=year+"_TT.root",              label="T#bar{T}",        scale=scale, color=(ROOT.kBlue - 6)),
        info.DataSetInfo(basedir=path, fileName=year+"_QCD.root",             label="QCD",             scale=scale, color=(ROOT.kGreen + 1)),
    ]

    sgData = [
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-20_rinv-0p3_alpha-peak.root",    label="t-ch 2000", scale=scale, color=ROOT.kCyan),
        info.DataSetInfo(basedir=path, fileName=year+"_mZprime-2100_mDark-20_rinv-0p3_alpha-peak.root", label="s-ch 2100", scale=scale, color=ROOT.kRed),
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-3000_mDark-20_rinv-0p3_alpha-peak.root",    label="t-ch 3000", scale=scale, color=ROOT.kGreen+2),
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
        locMin = min(0.02, max(0.02, 0.05 * ymin))
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
    

def plotStack(data, histoName, outputPath="./", xTitle="", yTitle="", isLogY=False, rebinx=-1.0, xmin=None, xmax=None):
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
    leg = ROOT.TLegend(0.20, 0.75, 0.95, 0.88)
    nColumns = 3 if(len(data[1]) >= 3) else 1
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetLineWidth(1)
    leg.SetNColumns(nColumns)
    leg.SetTextFont(42)

    #Setup background histos
    hs = ROOT.THStack()
    hMC = None
    firstPass = True
    for d in data[1]:
        h = d.getHisto(histoName, rebinx=rebinx, xmin=xmin, xmax=xmax, fill=True)
        hs.Add(copy.deepcopy(h))        
        leg.AddEntry(h, d.legEntry(), "F")
        if(firstPass):
            hMC = copy.deepcopy(h)
            firstPass = False
        else:
            hMC.Add(copy.deepcopy(h))

    #create a dummy histogram to act as the axes
    ymax=10**11
    ymin=10**-4
    lmax=10**12
    dummy = ROOT.TH1D("dummy", "dummy", 1000, hMC.GetBinLowEdge(1), hMC.GetBinLowEdge(hMC.GetNbinsX()) + hMC.GetBinWidth(hMC.GetNbinsX()))
    setupDummy(dummy, leg, histoName, xTitle, yTitle, isLogY, xmin, xmax, ymin, ymax, lmax)
    dummy.Draw("hist")    
    hs.Draw("hist F same")
    leg.Draw("same")

    #Setup signal histos
    history = []
    if(data[2]):
        for d in data[2]:
            h = d.getHisto(histoName, rebinx=rebinx, xmin=xmin, xmax=xmax)
            h.SetLineStyle(ROOT.kDashed)
            h.SetLineWidth(3)
            leg.AddEntry(h, d.legEntry(), "L")
            h.Draw("hist same")
            history.append(h)
    
    dummy.Draw("AXIS same")

    c1.SaveAs(outputPath+"/"+histoName+".png")
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
    cuts = ["", "_ge2AK8j"]
    #cuts = ["_ge2AK8j"]
    #cuts = [""]
    Data, sgData, bgData = getData("condor/MakeNJetsDists_"+year+"/", 1.0, year)
    
    for cut in cuts:
        plotStack((Data, bgData, sgData), "h_njets"+cut,       "./", "N_{j}",                   "Events / bin", isLogY=True, rebinx=-1, xmin=0, xmax=20)
        plotStack((Data, bgData, sgData), "h_njetsAK8"+cut,    "./", "N_{j}",                   "Events / bin", isLogY=True, rebinx=-1, xmin=0, xmax=12)
        plotStack((Data, bgData, sgData), "h_ntops"+cut,       "./", "N_{t}",                   "Events / bin", isLogY=True, rebinx=-1, xmin=0, xmax=6)
        plotStack((Data, bgData, sgData), "h_nb"+cut,          "./", "N_{b}",                   "Events / bin", isLogY=True, rebinx=-1)
        plotStack((Data, bgData, sgData), "h_ht"+cut,          "./", "H_{T}",                   "Events / bin", isLogY=True, rebinx=10)
        plotStack((Data, bgData, sgData), "h_st"+cut,          "./", "S_{T}",                   "Events / bin", isLogY=True, rebinx=10)
        plotStack((Data, bgData, sgData), "h_met"+cut,         "./", "MET",                     "Events / bin", isLogY=True, rebinx=10, xmin=0, xmax=2000)
        plotStack((Data, bgData, sgData), "h_jPt"+cut,         "./", "pT_{j}",                  "Events / bin", isLogY=True, rebinx=10)
        plotStack((Data, bgData, sgData), "h_jEta"+cut,        "./", "#eta_{j}",                "Events / bin", isLogY=True, rebinx=10)
        plotStack((Data, bgData, sgData), "h_jPhi"+cut,        "./", "#phi_{j}",                "Events / bin", isLogY=True, rebinx=10)
        plotStack((Data, bgData, sgData), "h_jPtAK8"+cut,      "./", "pT_{j8}",                 "Events / bin", isLogY=True, rebinx=10)
        plotStack((Data, bgData, sgData), "h_jEtaAK8"+cut,     "./", "#eta_{j8}",               "Events / bin", isLogY=True, rebinx=10)
        plotStack((Data, bgData, sgData), "h_jPhiAK8"+cut,     "./", "#phi_{j8}",               "Events / bin", isLogY=True, rebinx=10)
        plotStack((Data, bgData, sgData), "h_dEtaJ12"+cut,     "./", "#Delta#eta(j_{1},j_{2})", "Events / bin", isLogY=True, rebinx=5, xmin=0, xmax=6)
        plotStack((Data, bgData, sgData), "h_dRJ12"+cut,       "./", "#DeltaR(j_{1},j_{2})",    "Events / bin", isLogY=True, rebinx=2, xmin=0, xmax=6)
        plotStack((Data, bgData, sgData), "h_dPhiJ1MET"+cut,   "./", "#Delta#phi(j_{i},MET)",   "Events / bin", isLogY=True, rebinx=2, xmin=0, xmax=3.15)
        plotStack((Data, bgData, sgData), "h_dPhiJ2MET"+cut,   "./", "#Delta#phi(j_{2},MET)",   "Events / bin", isLogY=True, rebinx=2, xmin=0, xmax=3.15)
        plotStack((Data, bgData, sgData), "h_mjjM"+cut,        "./", "m_{jj}",                  "Events / bin", isLogY=True, rebinx=10) 
        plotStack((Data, bgData, sgData), "h_mjjPt"+cut,       "./", "pT_{jj}",                 "Events / bin", isLogY=True, rebinx=10)
        plotStack((Data, bgData, sgData), "h_mjjEta"+cut,      "./", "#eta_{jj}",               "Events / bin", isLogY=True, rebinx=10)
        plotStack((Data, bgData, sgData), "h_mT"+cut,          "./", "mT",                      "Events / bin", isLogY=True, rebinx=10) 
        plotStack((Data, bgData, sgData), "h_METrHT_pt30"+cut, "./", "MET/HT",                  "Events / bin", isLogY=True, rebinx=2, xmin=0, xmax=10)
        plotStack((Data, bgData, sgData), "h_METrST_pt30"+cut, "./", "MET/ST",                  "Events / bin", isLogY=True, rebinx=5)

if __name__ == '__main__':
    main()
