import ROOT

ROOT.gROOT.SetBatch( True )

def setStyle(h, color=None, lw=None, ls=None, ms=None, scale=None, fill=None):
    if color: 
        h.SetLineColor(color)
        h.SetMarkerColor(color)
    if lw:    h.SetLineWidth(lw)
    if ls:    h.SetLineStyle(ls)
    if ms:    h.SetMarkerStyle(ms)
    if scale: h.Scale(scale)
    if fill:  h.SetFillColor(fill)
    return h

def getHisto(f, name, color=None, lw=None, ls=None, ms=None, scale=None, fill=None):
    h = f.Get(name)
    h = setStyle(h, color, lw, ls, ms, scale, fill)
    return h

def cleanHistos(num, den):
    nXBins = den.GetXaxis().GetNbins()
    for i in range(0, nXBins+1):
        n = num.GetBinContent(i)
        d = den.GetBinContent(i)
        if(d == 0.0 and n > 0.0):
            print(n, d)
            num.SetBinContent(i, 0.0)
            den.SetBinContent(i, 0.0)
    return num, den

def drawPlots(f, suffix1, suffix2, scale=1.0, drawName="", rebin=None):
    num = getHisto(f, "h_{}".format(suffix1), ROOT.kBlack,  3, ROOT.kDotted, None, None, ROOT.kAzure+1)
    den = getHisto(f, "h_{}".format(suffix2), ROOT.kBlack , 3, 1)
    #num, den = cleanHistos(num,den)
    if rebin:
        num.Rebin(rebin)
        den.Rebin(rebin)

    #eff  = ROOT.TGraphAsymmErrors(num,den, "cp")
    eff = ROOT.TGraphAsymmErrors(num,den, "n")
    eff = setStyle(eff,  ROOT.kRed,   3, 1, ROOT.kFullDotLarge)
    eff.SetName(drawName)

    effT = ROOT.TEfficiency(num,den)
    #effT.SetStatisticOption(ROOT.TEfficiency.kFCP)
    effT.SetStatisticOption(ROOT.TEfficiency.kFNormal)
    effT = setStyle(effT, ROOT.kBlack, 3, 1, ROOT.kFullDotLarge)
    effT.SetName(drawName)

    scale = 1.0 / den.GetMaximum()
    num.Scale(scale)
    den.Scale(scale)

    canvas = ROOT.TCanvas("cv","cv",1000,800)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gPad.SetLogy()
    den.SetMinimum(1e-2)    
    den.SetMaximum(1.2)

    line = ROOT.TF1("line","1.0",-5000.0,5000.0)
    line.SetLineColor(ROOT.kBlack)

    den.Draw("hist")
    num.Draw("hist same")
    eff.Draw("P0 same")
    #effT.Draw("P0 same")
    line.Draw("same")

    canvas.SaveAs("fakerate{}.gif".format(drawName))
    canvas.SaveAs("fakerate{}.pdf".format(drawName))

    return eff, effT

def drawSummaryPlot(name, hList):
    canvas = ROOT.TCanvas("cv","cv",1000,800)
    colors = [ROOT.kRed, ROOT.kBlack, ROOT.kBlue, ROOT.kGreen+2]

    hBase = hList[0][0].GetHistogram()
    hBase.SetMaximum(1.5)
    hBase.Draw()

    legend = ROOT.TLegend(0.15,0.65,0.45,0.85)
    line = ROOT.TF1("line","1.0",-5000.0,5000.0)
    line.SetLineColor(ROOT.kBlack)
    
    for i, t in enumerate(hList):
        h , hT = t
        h.SetLineColor(colors[i])
        hT.SetLineColor(colors[i])
        h.SetMarkerColor(colors[i])
        hT.SetMarkerColor(colors[i])
        legend.AddEntry(h, h.GetName(),"l")
        h.Draw("same")
        #hT.Draw("same")

    line.Draw("same")
    legend.Draw()

    canvas.SaveAs("fakerate{}_summary.gif".format(name))
    canvas.SaveAs("fakerate{}_summary.pdf".format(name))


def main():    
    f = ROOT.TFile("test.root")
    hList = {"Pt_":[],"Eta":[]}
    hList["Pt_"].append(drawPlots(f, "svjPtAK8_qual_trg_st",  "jPtAK8_qual_trg_st",  1.0, "jPt_Fakerate",  10))
    hList["Eta"].append(drawPlots(f, "svjEtaAK8_qual_trg_st", "jEtaAK8_qual_trg_st", 1.0, "jEta_Fakerate", 10))

    hList["Pt_"].append(drawPlots(f, "svjPtAK8_qual_trg_st_0nim",  "jPtAK8_qual_trg_st_0nim",  1.0, "jPt_Fakerate_SR",  10))
    hList["Eta"].append(drawPlots(f, "svjEtaAK8_qual_trg_st_0nim", "jEtaAK8_qual_trg_st_0nim", 1.0, "jEta_Fakerate_SR", 10))

    hList["Pt_"].append(drawPlots(f, "svjPtAK8_qual_trg_st_ge1nim",  "jPtAK8_qual_trg_st_ge1nim",  1.0, "jPt_Fakerate_CR",  10))
    hList["Eta"].append(drawPlots(f, "svjEtaAK8_qual_trg_st_ge1nim", "jEtaAK8_qual_trg_st_ge1nim", 1.0, "jEta_Fakerate_CR", 10))

    outfile = ROOT.TFile("fakerate.root","RECREATE")
    for n,l in hList.items():
        drawSummaryPlot(n, l)
        for h, hT in l:
            h.Write()
    outfile.Close()

if __name__ == '__main__':
    main()
