import ROOT
import ctypes
import array

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
        # Check if the histogram exists in the file
        if not self.file.GetListOfKeys().Contains(str(name)):
            print("\33[31m" + "Error: Histo \"" + name + "\" not in file \"" + self.fileName + "\"" + "\033[0m")
            return None, None
        
        # Get the histogram
        histo = self.file.Get(name)
        histo.Sumw2()
        # Handle overflow by adding the overflow bin content to the last bin if specified
        if overflow:
            lastbinx = histo.GetNbinsX()
            lastbiny = histo.GetNbinsY()
            histo.Fill(histo.GetXaxis().GetBinCenter(lastbinx), histo.GetYaxis().GetBinCenter(lastbiny), histo.GetBinContent(lastbinx + 1, lastbiny + 1))

        # Apply scaling if specified
        if self.scale != -1.0:
            histo.Scale(self.scale)
        
        # Rebin histograms if needed
        if rebinx != -1.0:
            histo.RebinX(rebinx)
        if rebiny != -1.0:
            histo.RebinY(rebiny)
        
        # Apply range restrictions on the x and y axes
        if xmin is not None and xmax is not None:
            histo.GetXaxis().SetRangeUser(xmin, xmax)
        if ymin is not None and ymax is not None:
            histo.GetYaxis().SetRangeUser(ymin, ymax)

        # Find bin numbers for the specified range
        binx_min = histo.GetXaxis().FindBin(xmin) if xmin is not None else 1
        last_bin_x = histo.GetNbinsX()

        # Check if xmax is out of range, and include the overflow bin if overflow is True
        if xmax is not None and xmax > histo.GetXaxis().GetXmax() and overflow:
            binx_max = last_bin_x + 1  # Include the overflow bin
        else:
            binx_max = histo.GetXaxis().FindBin(xmax) if xmax is not None else last_bin_x

        biny_min = histo.GetYaxis().FindBin(ymin) if ymin is not None else 1
        biny_max = histo.GetYaxis().FindBin(ymax) if ymax is not None else histo.GetNbinsY()

        # Calculate the integral and its error within the correct bin range
        integral_error = ctypes.c_double(0.0)  # To store the error in the integral
        integral_normal = histo.Integral()
        integral = histo.IntegralAndError(binx_min, binx_max, biny_min, biny_max, integral_error)

        # print(f"in the datasetinfo.py, the integral and  error values given are - xmin - {binx_min}, xmax - {binx_max}, ymin = {biny_min}, ymax = {biny_max} ")
        # print(f"integral using histo.Integral - {integral_normal} and the integral using IntegralAndError - {integral}")
        # Save the number of events (integral)
        self.numEvents = integral

        # Update the label to include event count if specified
        if showEvents:
            self.label = "{0} ({1})".format(self.label_, int(round(self.numEvents)))
        else:
            self.label = self.label_
        
        return histo, integral, integral_error.value
    
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
        if isinstance(rebinx, list):  # If rebinx is a list of bin edges
            bin_array = array('d', rebinx)  # Convert list to array of doubles
            xProjectionhist = xProjectionhist.Rebin(len(bin_array) - 1, f"{xProjectionhist.GetName()}_rebinned", bin_array)
        elif rebinx != -1.0:  # If rebinx is a single integer
            xProjectionhist.Rebin(int(rebinx))
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
