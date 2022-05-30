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
    nsvjJetsAK8      = getHisto(f, "h_nsvjJetsAK8_qual_trg_st_0nim{}".format(suffix1),      ROOT.kRed,   3, None, None,               1.0)
    nsvjJetsAK8Plus1 = getHisto(f, "h_nsvjJetsAK8Plus1_qual_trg_st_0nim{}".format(suffix2), ROOT.kBlack, 3, 1,    ROOT.kFullDotLarge, scale)
    nsvjJetsAK8Plus2 = getHisto(f, "h_nsvjJetsAK8Plus2_qual_trg_st_0nim{}".format(suffix2), ROOT.kBlack, 3, 1,    ROOT.kFullDotLarge, scale)
    nsvjJetsAK8Plus3 = getHisto(f, "h_nsvjJetsAK8Plus3_qual_trg_st_0nim{}".format(suffix2), ROOT.kBlack, 3, 1,    ROOT.kFullDotLarge, scale)
    nsvjJetsAK8Plus4 = getHisto(f, "h_nsvjJetsAK8Plus4_qual_trg_st_0nim{}".format(suffix2), ROOT.kBlack, 3, 1,    ROOT.kFullDotLarge, scale)

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
    #f = ROOT.TFile("pbdCode.root")
    f = ROOT.TFile("test.root")
    drawPlots(f,      "", "_0J", 5.5, "pbdCode")
    drawPlots(f, "_ge1J", "_1J", 2.0, "pbdCode")
    drawPlots(f, "_ge2J", "_2J", 1.0, "pbdCode")

    #f = ROOT.TFile("myCode.test")
    #drawPlots(f,      "", "_0J", 1.85,"myCode")
    #drawPlots(f, "_ge1J", "_1J", 1.2 ,"myCode")
    #drawPlots(f, "_ge2J", "_2J", 1.0, "myCode")

if __name__ == '__main__':
    main()
