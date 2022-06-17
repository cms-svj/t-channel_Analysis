import ROOT

ROOT.gROOT.SetBatch( True )
ROOT.gStyle.SetOptStat("")
ROOT.gStyle.SetPaintTextFormat("3.2f")
ROOT.gStyle.SetFrameLineWidth(2)
ROOT.gStyle.SetEndErrorSize(0)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()

# Create a canvas and determine if it should be split for a ratio plot
# Margins are scaled on-the-fly so that distances are the same in either
# scenario.
def makeCanvas(doLogY):

    canvas = ROOT.TCanvas("cv", "cv", 1000, 800)

    # Split the canvas 70 / 30 by default if doing ratio
    # scale parameter keeps text sizes in ratio panel the
    # same as in the upper panel
    split      = 0.3
    upperSplit = 1.0
    lowerSplit = 1.0
    scale      = 1.0
    TopMargin    = 0.06
    BottomMargin = 0.12
    RightMargin  = 0.04
    LeftMargin   = 0.16
    upperSplit = 1.0-split
    lowerSplit = split
    scale = upperSplit / lowerSplit

    canvas.Divide(1,2)

    canvas.cd(1)
    ROOT.gPad.SetPad(0.0, split, 1.0, 1.0)
    ROOT.gPad.SetTopMargin(TopMargin / upperSplit)
    ROOT.gPad.SetBottomMargin(0)
    ROOT.gPad.SetLeftMargin(LeftMargin)
    ROOT.gPad.SetRightMargin(RightMargin)
    if doLogY:
        ROOT.gPad.SetLogy()

    canvas.cd(2)
    ROOT.gPad.SetPad(0.0, 0.0, 1.0, split)
    ROOT.gPad.SetTopMargin(0)
    ROOT.gPad.SetBottomMargin(BottomMargin / lowerSplit)
    ROOT.gPad.SetLeftMargin(LeftMargin)
    ROOT.gPad.SetRightMargin(RightMargin)

    return canvas

def getHisto(f, name, color=None, lw=None, ls=None, ms=None, scale=None):
    h = f.Get(name)
    if color: h.SetLineColor(color)
    if lw:    h.SetLineWidth(lw)
    if ls:    h.SetLineStyle(ls)
    if ms:    h.SetMarkerStyle(ms)
    if scale: h.Scale(scale)
    h.GetYaxis().SetTitle("Events")
    h.GetYaxis().SetLabelSize(0.05)
    return h

def drawPlots(f, suffix1, suffix2, scale=1.0, drawName=""):
    nsvjJetsAK8      = getHisto(f, "h_nsvjJetsAK8_qual_trg_st_{}".format(suffix1),      ROOT.kRed,   3, None, None,               1.0)
    nsvjJetsAK8Plus1 = getHisto(f, "h_nsvjJetsAK8Plus1_qual_trg_st_{}".format(suffix1+suffix2), ROOT.kBlack, 3, 1,    ROOT.kFullDotLarge, scale)
    nsvjJetsAK8Plus2 = getHisto(f, "h_nsvjJetsAK8Plus2_qual_trg_st_{}".format(suffix1+suffix2), ROOT.kBlack, 3, 1,    ROOT.kFullDotLarge, scale)
    nsvjJetsAK8Plus3 = getHisto(f, "h_nsvjJetsAK8Plus3_qual_trg_st_{}".format(suffix1+suffix2), ROOT.kBlack, 3, 1,    ROOT.kFullDotLarge, scale)
    nsvjJetsAK8Plus4 = getHisto(f, "h_nsvjJetsAK8Plus4_qual_trg_st_{}".format(suffix1+suffix2), ROOT.kBlack, 3, 1,    ROOT.kFullDotLarge, scale)
    nsvjJetsAK8Plus1.Add(nsvjJetsAK8Plus2)
    nsvjJetsAK8Plus1.Add(nsvjJetsAK8Plus3)
    nsvjJetsAK8Plus1.Add(nsvjJetsAK8Plus4)
    nsvjJetsAK8Pred = nsvjJetsAK8Plus1

    ratio = nsvjJetsAK8.Clone("ratio")
    ratio.SetTitle("")
    ratio.GetYaxis().SetTitle("True/Pred.")
    ratio.Divide(nsvjJetsAK8Pred)
    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetLineWidth(3)
    ratio.SetLineStyle(1)
    ratio.SetMarkerStyle(ROOT.kFullDotLarge)
    ratio.GetXaxis().SetLabelSize(0.08)
    ratio.GetYaxis().SetLabelSize(0.08)
    ratio.GetXaxis().SetTitleSize(0.08)
    ratio.GetYaxis().SetTitleSize(0.08)
    ratio.SetMaximum(21.0)

    #fit = ROOT.TF1('fit','[0]*TMath::Poisson(x-[1],[2])',1+2,5+2)
    #fit.SetLineColor(ROOT.kRed)
    #fit.SetParameter(0,100.0)
    #fit.SetParameter(1, -1)
    #fit.SetParameter(2, 10.0)
    #ratio.Fit(fit,"Q", "", 1+2, 5+2)

    canvas = makeCanvas(doLogY=True)
    ROOT.gStyle.SetOptStat(0)
    nsvjJetsAK8.SetMinimum(1.5e-2)    
    nsvjJetsAK8.SetMaximum(1e8)

    canvas.cd(1)
    nsvjJetsAK8.Draw("hist E")
    nsvjJetsAK8Pred.Draw("PEX0 same")

    canvas.cd(2)
    ROOT.gPad.SetGridy()
    ratio.Draw("E0P")
    #fit.Draw("same")

    canvas.SaveAs("fakerateBGEstimation{}{}.gif".format(suffix2, drawName))
    canvas.SaveAs("fakerateBGEstimation{}{}.pdf".format(suffix2, drawName))

def main():    
    f = ROOT.TFile("2018_QCD.root")
    drawPlots(f, "0nim", "_0J", 1.3, "SR")
    drawPlots(f, "0nim", "_1J", 1.1, "SR")
    drawPlots(f, "0nim", "_2J", 1.1, "SR")

    drawPlots(f, "ge1nim", "_0J", 2.5,"CR")
    drawPlots(f, "ge1nim", "_1J", 2.0 ,"CR")
    drawPlots(f, "ge1nim", "_2J", 1.5, "CR")

    #f = ROOT.TFile("test.root")
    #drawPlots(f, "0nim", "_0J", 5.5, "SR")
    #drawPlots(f, "0nim", "_1J", 2.1, "SR")
    #drawPlots(f, "0nim", "_2J", 1.15, "SR")

    #drawPlots(f, "ge1nim", "_0J", 17.0,"CR")
    #drawPlots(f, "ge1nim", "_1J", 4.2 ,"CR")
    #drawPlots(f, "ge1nim", "_2J", 1.2, "CR")

if __name__ == '__main__':
    main()
