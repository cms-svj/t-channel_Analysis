import ROOT

ROOT.gROOT.SetBatch( True )

def getHisto(f, name, color=None, lw=None, ls=None, ms=None, scale=None):
    h = f.Get(name)
    if color: h.SetLineColor(color)
    if lw:    h.SetLineWidth(lw)
    if ls:    h.SetLineStyle(ls)
    if ms:    h.SetMarkerStyle(ms)
    if scale: h.Scale(scale)
    return h

def drawPlots(f, suffix1, suffix2, scale=1.0, drawName=""):
    nsvjJetsAK8      = getHisto(f, "h_nsvjJetsAK8_qual_trg_st_{}".format(suffix1),      ROOT.kRed,   3, None, None,               1.0)
    nsvjJetsAK8Plus1 = getHisto(f, "h_nsvjJetsAK8Plus1_qual_trg_st_{}".format(suffix1+suffix2), ROOT.kBlack, 3, 1,    ROOT.kFullDotLarge, scale)
    nsvjJetsAK8Plus2 = getHisto(f, "h_nsvjJetsAK8Plus2_qual_trg_st_{}".format(suffix1+suffix2), ROOT.kBlack, 3, 1,    ROOT.kFullDotLarge, scale)
    nsvjJetsAK8Plus3 = getHisto(f, "h_nsvjJetsAK8Plus3_qual_trg_st_{}".format(suffix1+suffix2), ROOT.kBlack, 3, 1,    ROOT.kFullDotLarge, scale)
    nsvjJetsAK8Plus4 = getHisto(f, "h_nsvjJetsAK8Plus4_qual_trg_st_{}".format(suffix1+suffix2), ROOT.kBlack, 3, 1,    ROOT.kFullDotLarge, scale)

    canvas = ROOT.TCanvas("cv","cv",1000,800)
    #ROOT.gPad.SetLogy()
    nsvjJetsAK8.SetMinimum(1e-1)    
    nsvjJetsAK8.SetMaximum(1e4)


    nsvjJetsAK8.Draw("hist E")
    nsvjJetsAK8Plus1.Draw("PEX0 same")
    nsvjJetsAK8Plus2.Draw("PEX0 same")
    nsvjJetsAK8Plus3.Draw("PEX0 same")
    nsvjJetsAK8Plus4.Draw("PEX0 same")

    canvas.SaveAs("fakerateBGEstimation{}{}.gif".format(suffix2, drawName))
    canvas.SaveAs("fakerateBGEstimation{}{}.pdf".format(suffix2, drawName))

def main():    
    f = ROOT.TFile("test.root")
    drawPlots(f, "0nim", "_0J", 5.5, "SR")
    drawPlots(f, "0nim", "_1J", 2.1, "SR")
    drawPlots(f, "0nim", "_2J", 1.15, "SR")

    drawPlots(f, "ge1nim", "_0J", 17.0,"CR")
    drawPlots(f, "ge1nim", "_1J", 4.2 ,"CR")
    drawPlots(f, "ge1nim", "_2J", 1.2, "CR")

if __name__ == '__main__':
    main()
