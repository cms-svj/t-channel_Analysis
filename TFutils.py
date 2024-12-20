import ROOT
import numpy as np
import os

def GetABCDregions(met,dnn):
    '''Defining the A B C D regions'''
    regions = [ ("A", met, 20000, dnn, 1.0),
                    ("B", met, 20000, 0, dnn),
                    ("C", 0, met, dnn, 1.0),
                    ("D", 0, met, 0, dnn)
                ]
    return regions

# def  GetABCDhistinMETProjections(data, data_type, ABCDhistoVar, maincut, SVJbin, binning):
#     """Return the A,B,C,D histograms for the given SVJbin using the 2D variable provided for the ABCDhistVar"""
#     if data_type == "bkg":
#         filename = data.fileName.split()[1].replace('.root','')
#     elif data_type == "sgn":
#         filename = "sgn"
#     elif data_type == "data":
#         filename = "obs"
    
#     histName =ABCDhistoVar + maincut + SVJbin
#     dnn, met = SVJbin

def GetABCDhistPerSVJBin(data, ABCDhistoVar, maincut, SVJbin):
    """Return the dictionary of ABCD histogram for only one svj bin"""

    bkgname = data.fileName.split('_')[1].replace('.root','')
    hist_dict = {region: ROOT.TH1F(f"h_{bkgname}_{region}",f"h_{bkgname}_{region}", 1,0,1) for region in ['A','B','C','D']}
    SVJ, (dnn, met) = SVJbin
    histName = ABCDhistoVar + maincut + SVJ
    regions = GetABCDregions(met,dnn)

    for region, xmin, xmax, ymin, ymax in regions:
        hist, histIntegral, integral_error = data.get2DHistoIntegral(histName, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, showEvents=True)
        hist_dict[region].SetBinContent(1, histIntegral)
        hist_dict[region].SetBinError(1, integral_error)
        hist_dict[region].Sumw2()
    return hist_dict



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
    bkgname = data.fileName.split('_')[1].replace('.root', '')
    # Determine histogram bin details
    xmin, binwidth = int(list(SVJbins.keys())[0][0]), len(SVJbins)
    xmax = xmin + binwidth
    nBins_total = binwidth * 4  # Total bins for merged histogram
    xmax_merged = nBins_total  # Assuming bin width is the same for all histograms
    # print(f"xmin = {xmin}, xmax = {xmax}, nBins_total = {nBins_total}, binwidht = {binwidth}, bin_width = {xmax_merged}")
    # Initialize histograms for regions A, B, C, D
    hist_dict = {region: ROOT.TH1F(f"h_{bkgname}_{region}", f"h_{bkgname}_{region}_{maincut}", binwidth, xmin, xmax) for region in ['A', 'B', 'C', 'D']}
    
    # Initialize merged histogram if needed
    if merge:
        merged_hist = ROOT.TH1F(f"merged_{bkgname}_{maincut}", f"merged_{bkgname}_{maincut}", nBins_total, 0, xmax_merged)
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


def GetABCDhistDict(dataList, ABCDhistoVar, maincut, SVJbins, isStack=False, merge=False,perSVJbin=False):
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
        elif perSVJbin:
            histograms = GetABCDhistPerSVJBin(data, ABCDhistoVar, maincut, SVJbins)
            ABCDhistDict.update({data: [histograms['A'], histograms['B'], histograms['C'], histograms['D']]})
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
            sumHistograms = [hist.Clone("Summed") for hist in histograms]
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
    

def SaveHistDictToFile(hist_data, rootfileName):
    """
    Saves the histograms in hist_data to a ROOT file.
    Supports both dictionaries (region: histogram) and lists of histograms.
    
    For lists, the histograms are saved in folders "A", "B", "C", "D" in that order.
    
    Parameters:
        hist_data: Either a dictionary {region: histogram} or a list of histograms (interpreted as A, B, C, D).
        rootfileName: Name of the output ROOT file.
    """
    # Ensure the directory exists
    output_dir = os.path.dirname(rootfileName)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)  # Create directories if they don't exist
    
    # Create the ROOT file
    output_file = ROOT.TFile(rootfileName, "RECREATE")
    
    if isinstance(hist_data, ROOT.TH1):  # Single histogram
        folder = output_file.mkdir('A')
        folder.cd()
        hist_data.Write()  # Write the histogram directly
    elif isinstance(hist_data, dict):  # If it's a dictionary
        for region, hist in hist_data.items():
            # Create or navigate to the folder corresponding to the region
            folder = output_file.mkdir(region)
            folder.cd()  # Change directory to the folder
            
            # Write the histogram into the folder
            hist.Write()
    elif isinstance(hist_data, list):  # If it's a list
        if len(hist_data) != 4:
            raise ValueError("For list input, exactly 4 histograms are expected (for A, B, C, D).")
        
        # Define folder names for regions
        regions = ["A", "B", "C", "D"]
        for region, hist in zip(regions, hist_data):
            # Create or navigate to the folder corresponding to the region
            folder = output_file.mkdir(region)
            folder.cd()  # Change directory to the folder
            
            # Write the histogram into the folder
            hist.Write()
    else:
        raise ValueError("Unsupported data type for hist_data. Must be dict or list.")
    
    # Close the ROOT file
    output_file.Close()
    print(f"Histograms saved successfully to {rootfileName}.")