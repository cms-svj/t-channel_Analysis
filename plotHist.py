# General code in ROOT to plot 1D and 2D histograms with rebining and other options.

import ROOT

class Histogram:
    def __init__(self, hist):
        self.hist = hist

    def beautify(self, title="", xtitle="", ytitle="", binlabels=None, rebin=1, underflow=False, overflow=False):
        # Set the histogram title and axis labels
        self.hist.SetTitle(title)
        self.hist.GetXaxis().SetTitle(xtitle)
        self.hist.GetYaxis().SetTitle(ytitle)

        # Set the bin labels if provided
        if binlabels is not None:
            for i, label in enumerate(binlabels):
                self.hist.GetXaxis().SetBinLabel(i+1, label)

        # Rebin the histogram if requested
        if rebin > 1:
            self.hist.Rebin(rebin)

        # Include underflow and overflow bins if requested
        if underflow:
            self.hist.SetBinContent(1, self.hist.GetBinContent(0) + self.hist.GetBinContent(1))
            self.hist.SetBinContent(0, 0)
        if overflow:
            self.hist.SetBinContent(self.hist.GetNbinsX(), self.hist.GetBinContent(self.hist.GetNbinsX()) + self.hist.GetBinContent(self.hist.GetNbinsX()+1))
            self.hist.SetBinContent(self.hist.GetNbinsX()+1, 0)

        # Set the font size
        ROOT.gStyle.SetLabelSize(0.04, "XYZ")
        ROOT.gStyle.SetTitleSize(0.05, "XYZ")

        # Set the pad margins
        ROOT.gPad.SetLeftMargin(0.15)
        ROOT.gPad.SetRightMargin(0.05)
        ROOT.gPad.SetBottomMargin(0.15)
        ROOT.gPad.SetTopMargin(0.05)

    def plot_1d(self, title="", xtitle="", ytitle="", binlabels=None, rebin=1, underflow=False, overflow=False, draw_options=""):
        # Beautify the histogram
        self.beautify(title=title, xtitle=xtitle, ytitle=ytitle, binlabels=binlabels, rebin=rebin, underflow=underflow, overflow=overflow)

        # Draw the 1D histogram
        c = ROOT.TCanvas()
        self.hist.Draw(draw_options)
        c.Draw()

    def plot_2d(self, title="", xtitle="", ytitle="", xmin=None, xmax=None, ymin=None, ymax=None, binlabels=None, rebin=1, underflow=False, overflow=False, draw_options="colz"):
        # Beautify the histogram
        self.beautify(title=title, xtitle=xtitle, ytitle=ytitle, binlabels=binlabels, rebin=rebin, underflow=underflow, overflow=overflow)
        # set range 
        if xmin is not None and xmax is not None:
            self.hist.GetXaxis().SetRangeUser(xmin, xmax)
        if ymin is not None and ymax is not None:
            self.hist.GetYaxis().SetRangeUser(ymin, ymax)
        # Draw the 2D histogram
        c = ROOT.TCanvas()
        self.hist.Draw(draw_options)
        c.Draw()
        return c

    def plot_ratio(self, hist2, title="", xtitle="", ytitle="", binlabels=None, rebin=1, underflow=False, overflow=False, draw_options=""):
        # Beautify the histograms
        self.beautify(title=title, xtitle=xtitle, ytitle=ytitle, binlabels=binlabels, rebin=rebin, underflow=underflow, overflow=overflow)
        hist2.beautify(title=title, xtitle=xtitle, ytitle=ytitle, binlabels=binlabels, rebin=rebin, underflow=underflow, overflow=overflow)

        # Create the ratio plot
        ratio = self.hist.Clone()
        ratio.Divide(hist2.hist)

        # Beautify the ratio plot
        ratio.SetTitle("")
        ratio.GetXaxis().SetTitle(self.hist.GetXaxis().GetTitle())
        ratio.GetYaxis().SetTitle("Ratio")
        ratio.GetYaxis().SetRangeUser(0.0, 2.0)

        # Draw the ratio plot
        c = ROOT.TCanvas()
        ratio.Draw(draw_options)
        c.Draw()
