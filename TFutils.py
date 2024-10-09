import ROOT
import numpy as np

def GetABCDregions(met,dnn):
    '''Defining the A B C D regions'''
    regions = [ ("A", met, 20000, dnn, 1.0),
                    ("B", met, 20000, 0, dnn),
                    ("C", 0, met, dnn, 1.0),
                    ("D", 0, met, 0, dnn)
                ]
    return regions


def GetABCDhist(data, ABCDhistoVar, maincut, SVJbins, isStack=False, merge=False):
    """
    Computes individual A, B, C, D histograms and optionally merges them into a single histogram.
    
    Parameters:
        data: Data object used to retrieve histograms and perform calculations.
        ABCDhistoVar: Variable name for histograms.
        maincut: Main cut used in histogram names.
        SVJbins: Dictionary containing the SVJ bins information.
        isStack: Boolean flag to indicate if histograms should be styled for stacking.
        merge: Boolean flag to indicate if histograms should be merged.
    
    Returns:
        hist_dict: Dictionary of histograms for regions A, B, C, and D.
        merged_hist: Merged histogram if merge=True; otherwise, None.
    """
    
    # Determine histogram bin details
    xmin, binwidth = int(list(SVJbins.keys())[0][0]), len(SVJbins)
    xmax = xmin + binwidth
    nBins_total = binwidth * 4  # Total bins for merged histogram
    xmax_merged = nBins_total  # Assuming bin width is the same for all histograms
    # print(f"xmin = {xmin}, xmax = {xmax}, nBins_total = {nBins_total}, binwidht = {binwidth}, bin_width = {xmax_merged}")
    # Initialize histograms for regions A, B, C, D
    hist_dict = {region: ROOT.TH1F(f"{region}_{data.label_}_{maincut}", f"{region}_{data.label_}_{maincut}", binwidth, xmin, xmax) for region in ['A', 'B', 'C', 'D']}
    
    # Initialize merged histogram if needed
    if merge:
        merged_hist = ROOT.TH1F(f"merged_{data.label_}_{maincut}", f"merged_{data.label_}_{maincut}", nBins_total, 0, xmax_merged)
    else:
        merged_hist = None
    error = None
    # Process each SVJ bin
    for i, SVJ in enumerate(SVJbins.keys()):
        histName = ABCDhistoVar + maincut + SVJ
        dnn, met = SVJbins[SVJ]
        regions = GetABCDregions(met, dnn)
        
        for region, xmin, xmax, ymin, ymax in regions:
            hist, histIntegral, integral_error = data.get2DHistoIntegral(histName, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, showEvents=True)
            hist_dict[region].SetBinContent(i + 1, histIntegral)
            hist_dict[region].SetBinError(i+1, integral_error)
            # print(f"the hist_dict for region - {region}, bin - {i+1} the integral is  - {histIntegral} and error - {integral_error}, xmin = {xmin}, xmax = {xmax}, ymin = {ymin}, ymax = {ymax}")
            hist_dict[region].Sumw2()
        
        # Add bin contents to the merged histogram if merging is enabled
    if merge:
        for region_index, region in enumerate(['A', 'B', 'C', 'D']):
            offset = region_index * binwidth
            # print(f"region index = {region_index}, binwidth = {binwidth}, offset = {offset}, SVJ = {SVJ}, region = {region}, hist_dict = {hist_dict[region]}")
            for bin, SVJ in enumerate(SVJbins.keys()):
                merged_hist.SetBinContent(bin + offset+1, hist_dict[region].GetBinContent(bin+1))
                merged_hist.SetBinError(bin+offset+1, hist_dict[region].GetBinError(bin+1))
                merged_hist.Sumw2()
                if SVJ == "2PSVJ":
                    merged_hist.GetXaxis().SetBinLabel(bin + offset+1, f"{SVJ[0]}+")
                else:
                    merged_hist.GetXaxis().SetBinLabel(bin + offset+1, SVJ[0])
                # print(f"bin = {bin}, offset = {offset}")

    # Apply stack styling if needed
    if isStack:
        for hist in hist_dict.values():
            hist.SetFillColor(data.getColor())
            hist.SetFillStyle(2001)
    if merge:
        if isStack:
            merged_hist.SetFillColor(data.getColor())
            merged_hist.SetFillStyle(3001)
        return merged_hist
    else:
        return hist_dict


def GetABCDhistDict(dataList, ABCDhistoVar, maincut, SVJbins, isStack=False, merge=False):
    """
    Returns a dictionary containing A, B, C, D histograms or only merged histograms based on the merge flag.
    
    Parameters:
        dataList: List of data objects.
        ABCDhistoVar: Variable name for histograms.
        maincut: Main cut used in histogram names.
        SVJbins: Dictionary containing the SVJ bins information.
        isStack: Boolean flag to indicate if histograms should be styled for stacking.
        merge: Boolean flag to indicate if histograms should be merged.
    
    Returns:
        ABCDhistDict: Dictionary where keys are file names and values are either:
                      - List of histograms [A, B, C, D] if merge=False,
                      - Merged histogram if merge=True.
    """
    ABCDhistDict = {}
    
    for data in dataList:
        if merge:
            merged_hist = GetABCDhist(data, ABCDhistoVar, maincut, SVJbins, isStack=isStack, merge=merge)
            ABCDhistDict.update({data: merged_hist})
        else:
            histograms = GetABCDhist(data, ABCDhistoVar, maincut, SVJbins, isStack=isStack, merge=False)
            ABCDhistDict.update({data: [histograms['A'], histograms['B'], histograms['C'], histograms['D']]})
    
    return ABCDhistDict


def SumABCDhistList(ABCDhistList, skipList=None):
    '''Returns A,B,C,D histograms which are summed over all the files in the data list. It skips the files provided in the skipList'''
    firstPass = True
    for d, histograms in ABCDhistList.items():
        if skipList != None and d.fileName in skipList:
            continue
        if firstPass:
            sumHistograms = [hist.Clone() for hist in histograms]
            for sumHist in sumHistograms:
                sumHist.Reset()
            firstPass = False
        print(f"fileName = {d.fileName}")
        for sumHist, hist in zip(sumHistograms, histograms):
            sumHist.Add(hist)
    return sumHistograms

def SumHistograms(hist_dict, skipList=None):
    """
    Sums histograms from a dictionary of histograms, with an option to skip specific histograms.
    
    Parameters:
        hist_dict (dict): Dictionary of histograms to sum.
        skipList (list, optional): List of keys to skip while summing histograms.
    
    Returns:
        ROOT.TH1: The summed histogram.
    """
    first_key = next(iter(hist_dict))
    summed_hist = hist_dict[first_key].Clone("summed_hist")
    summed_hist.Reset()  # Clear the content, keep the structure

    # Iterate over the histograms and sum them
    for d, hist in hist_dict.items():
        if skipList != None and d.fileName in skipList:
            continue  # Skip this histogram if its key is in the skipList
        summed_hist.Add(hist)  # Sum the histograms

    return summed_hist

def rebinHistogramstoABCD(histList):
    """Takes a list of histograms and returns a new list of histograms rebinned according to the bins of A, B, C, and D."""# TODO: Figure out the issue, not rebinning properly as desired
    newHistList = []
    for hist in histList:
        newHist = hist.Clone()
        newHist.Reset()
        for bin in range(1, hist.GetNbinsX() + 1):
            newHist.Fill(hist.GetBinContent(bin))
        newHistList.append(newHist)
    return newHistList


def Validation(Data_CR, TFhist, SR_nonLL):
    predicted_Data_SR = Data_CR.Clone("predicted_Data_SR")
    predicted_Data_SR.Multiply(TFhist)  # TFhist * Data_CR
    predicted_Data_SR.Add(SR_nonLL)  # Add MC QCD + ZJets in SR
    return predicted_Data_SR
    
