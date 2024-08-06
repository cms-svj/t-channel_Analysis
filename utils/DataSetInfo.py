import ROOT

class DataSetInfo:
    def __init__(self, basedir, fileName, label, color=None, sys=None, processName=None, process=None, rate=None, lumiSys=None, scale=-1.0):
        self.basedir = basedir
        self.fileName = fileName
        self.label_ = label
        self.color = color
        try:
            self.file = ROOT.TFile.Open(basedir+fileName)
        except:
            pass
        self.sys = sys
        self.processName = processName
        self.process = process
        self.rate = rate
        self.lumiSys = lumiSys
        self.scale = scale
        self.numEvents = None

    def getHisto(self, name, rebinx=-1.0, rebiny=-1.0, xmin=None, xmax=None, fill=None, showEvents=False, overflow=True):
        if not self.file.GetListOfKeys().Contains(str(name)):
            print ("\33[31m"+"Error: Histo \""+name+"\" not in file \""+self.fileName+"\""+"\033[0m")
        histo = self.file.Get(name)
        if overflow:
            lastbin = histo.GetNbinsX()
            histo.Fill(histo.GetBinCenter(lastbin), histo.GetBinContent(lastbin + 1))
        if(self.scale != -1.0): histo.Scale(self.scale)
        if(rebinx != -1.0): histo.RebinX(rebinx)
        if(rebiny != -1.0): histo.RebinY(rebiny)
        if(xmin != None): histo.GetXaxis().SetRangeUser(xmin, xmax)
        if(self.color):
            histo.SetLineColor(self.color)
            if(fill): histo.SetFillColor(self.color)
        if(showEvents):
            self.numEvents = histo.Integral()
            self.label = "{0} ({1})".format(self.label_, int(round(self.numEvents)))
        else:
            self.label = self.label_
        return histo
    
    def get2DHistoIntegral(self, name, rebinx=-1.0, rebiny=-1.0, xmin=None, xmax=None, ymin=None, ymax=None, fill=None, showEvents=False, overflow=True):
        if not self.file.GetListOfKeys().Contains(str(name)):
            print("\33[31m" + "Error: Histo \"" + name + "\" not in file \"" + self.fileName + "\"" + "\033[0m")
        # print(f"xmin = {xmin}, xmax = {xmax}, ymin = {ymin}, ymax = {ymax}")
        histo = self.file.Get(name)
        if overflow:
            lastbinx = histo.GetNbinsX() # Should this be +1 ?
            lastbiny = histo.GetNbinsY()
            histo.Fill(histo.GetXaxis().GetBinCenter(lastbinx), histo.GetYaxis().GetBinCenter(lastbiny), histo.GetBinContent(lastbinx+1, lastbiny+1))
            
        if self.scale != -1.0:
            histo.Scale(self.scale)
        if rebinx != -1.0:
            histo.RebinX(rebinx)
        if rebiny != -1.0:
            histo.RebinY(rebiny)
        if xmin is not None and xmax is not None:
            histo.GetXaxis().SetRangeUser(xmin, xmax)
        if ymin is not None and ymax is not None:
            histo.GetYaxis().SetRangeUser(ymin, ymax)
        # Normal integral    
        self.numEvents = histo.Integral()
        # if withError:
        #     error = Double(0)

        if showEvents:
            self.label = "{0} ({1})".format(self.label_, int(round(self.numEvents)))
        else:
            self.label = self.label_
        return histo, self.numEvents

    def getXProjection(self, name, yValue, condition, rebinx=-1.0, rebiny=-1.0, xmin=None, xmax=None):
        histo = self.file.Get(name)
        nYbins = histo.GetNbinsY()
        yValueBin = histo.GetYaxis().FindBin(yValue)
        if condition == '>':
            print("This is envoked with the condition >")
            xProjectionhist = histo.ProjectionX(f"Projection_greater_{yValue}", yValueBin, nYbins+1) # Overflow bin is considered
        elif condition == '<':
            print("This is envoked with the condition <")
            xProjectionhist = histo.ProjectionX(f"Projection_less_{yValue}", 0, yValueBin) # Underflow bin is considered 
        
        if self.scale != -1.0: xProjectionhist.Scale(self.scale)
        if rebinx != -1.0: xProjectionhist.RebinX(rebinx)
        if rebiny != -1.0: xProjectionhist.RebinY(rebiny)
        if xmin is not None and xmax is not None:
            xProjectionhist.GetXaxis().SetRangeUser(xmin,xmax)
        return xProjectionhist

    def getFile(self):
        return self.file

    def legEntry(self):
        return self.label

    def getColor(self):
        return self.color

    def __del__(self):
        self.file.Close()
