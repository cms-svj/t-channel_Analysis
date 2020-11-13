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

    def getHisto(self, name, rebinx=-1.0, rebiny=-1.0, xmin=None, xmax=None, fill=None, showEvents=False):
        if not self.file.GetListOfKeys().Contains(str(name)):
            print ("\33[31m"+"Error: Histo \""+name+"\" not in file \""+self.fileName+"\""+"\033[0m")
        histo = self.file.Get(name)
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

    def getFile(self):
        return self.file

    def legEntry(self):
        return self.label

    def __del__(self):
        self.file.Close()
