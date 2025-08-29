import ROOT
import numpy as np
import os


def GetABCDregions(met,dnn):
    #Defining the A B C D regions
    print('MeT',met)
    print('dnn',dnn/2)
    #print('dnn/3',dnn/2)
    regions = [ ("dA", met, 20000, dnn/2, dnn),
                    ("dB", met, 20000, 0, dnn/2),
                    ("dC", 0, met, dnn/2,dnn),
                    ("dD", 0, met, 0, dnn/2)
                ]
    return regions
'''
def GetABCDregions(met, dnn, dnnFractions=None):
    ##Defining multiple ABCD regions using different DNN fraction divisions
    if dnnFractions is None:
        dnnFractions = [2]  # Default to dnn/2 if no list is provided

    all_regions = {}
    for frac in dnnFractions:
        dnn_cut = dnn / frac
        regions = [
            ("A", met, 20000, dnn_cut, dnn),
            ("B", met, 20000, 0, dnn_cut),
            ("C", 0, met, dnn_cut, dnn),
            ("D", 0, met, 0, dnn_cut)
        ]
        all_regions[f"dnn_div_{frac}"] = regions  # Store regions with a unique key

    return all_regions
'''

def adjustRegionBoundaries(region, xmin, xmax, ymin, ymax):
    """ Adjust regions to avoid bin overlap based on the detected overlap. """
    if region == 'dA':
        xmin = xmin + 0.1 
        ymin = ymin + 0.01
    if region == 'dB':
        xmin = xmin + 0.1
    if region == 'dC':
        ymin = ymin + 0.01    

    return xmin, xmax, ymin, ymax

def GetABCDhistPerSVJBin(data, ABCDhistoVar, maincut, SVJbin):
    """Return the dictionary of ABCD histogram for only one svj bin"""

    bkgname = data.fileName.split('_')[1].replace('.root','')
    hist_dict = {region: ROOT.TH1F(f"h_{bkgname}_{region}",f"h_{bkgname}_{region}", 1,0,1) for region in ['dA','dB','dC','dD']}
    SVJ, (dnn, met) = SVJbin
    histName = ABCDhistoVar + maincut + SVJ
    regions = GetABCDregions(met,dnn)

    for region, xmin, xmax, ymin, ymax in regions:   # TODO: this for condition should be written into a function work in all the cases
        xmin, xmax, ymin, ymax = adjustRegionBoundaries(region, xmin, xmax, ymin, ymax)
        hist, histIntegral, integral_error = data.get2DHistoIntegral(histName, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, showEvents=True)
        print(f"hist integral is - {hist.Integral()} in region - {region} and from the function is  - {histIntegral}, xmin - {xmin}, xmax - {xmax}, ymin - {ymin}, ymax - {ymax}")
        hist_dict[region].SetBinContent(1, histIntegral)
        hist_dict[region].SetBinError(1, integral_error)
        hist_dict[region].Sumw2()
    return hist_dict

def checkIntegralConsistency(data, ABCDhistoVar, maincut, SVJbin):
    """Check if the sum of the integrals of all regions equals the original integral"""

    bkgname = data.fileName.split('_')[1].replace('.root','')
    SVJ, (dnn, met) = SVJbin
    histName = ABCDhistoVar + maincut + SVJ
    print(f"The file is - {data.fileName} and the SVJ region is - {SVJbin}")
    # Get the original integral without considering the regions
    original_hist, original_integral, original_integral_error = data.get2DHistoIntegral(histName, showEvents=False)

    # Initialize the sum of the integrals of all regions
    total_integral = 0.0

    # Define the regions
    regions = GetABCDregions(met, dnn)

    # Sum the integrals over the regions
    for region, xmin, xmax, ymin, ymax in regions:
        xmin, xmax, ymin, ymax = adjustRegionBoundaries(region, xmin, xmax, ymin, ymax)
        hist, histIntegral, integral_error = data.get2DHistoIntegral(histName, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, showEvents=False)
        # print(f"the ABCD region is - {region}, and integral is - {histIntegral}")
        print(f"hist integral is - {hist.Integral()} in region - {region} and from the function is  - {histIntegral}, xmin - {xmin}, xmax - {xmax}, ymin - {ymin}, ymax - {ymax}")
        
        total_integral += histIntegral

    # Compare the sum of the integrals of the regions with the original integral
    print(f"Original integral: {original_integral}")
    print(f"Total integral of all regions: {total_integral}")

    # Check if they are approximately equal
    if abs(original_integral - total_integral) < 1e-5:  # Tolerance can be adjusted
        print("The sum of the region integrals matches the original integral.")
        return True
    else:
        print("The sum of the region integrals does NOT match the original integral.")
        return False

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
    hist_dict = {region: ROOT.TH1F(f"h_{bkgname}_{region}", f"h_{bkgname}_{region}_{maincut}", binwidth, xmin, xmax) for region in ['dA', 'dB', 'dC', 'dD']}
    
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
            xmin, xmax, ymin, ymax = adjustRegionBoundaries(region, xmin, xmax, ymin, ymax)
            hist, histIntegral, integral_error = data.get2DHistoIntegral(histName, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, showEvents=True)
            hist_dict[region].SetBinContent(i + 1, histIntegral)
            hist_dict[region].SetBinError(i+1, integral_error)
            # print(f"the hist_dict for region - {region}, bin - {i+1} the integral is  - {histIntegral} and error - {integral_error}, xmin = {xmin}, xmax = {xmax}, ymin = {ymin}, ymax = {ymax}")
            hist_dict[region].Sumw2()
        
        # Add bin contents to the merged histogram if merging is enabled
    if merge:
        for region_index, region in enumerate(['dA', 'dB', 'dC', 'dD']):
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
            ABCDhistDict.update({data: [histograms['dA'], histograms['dB'], histograms['dC'], histograms['dD']]})
        else:
            histograms = GetABCDhist(data, ABCDhistoVar, maincut, SVJbins, isStack=isStack, merge=False)
            ABCDhistDict.update({data: [histograms['dA'], histograms['dB'], histograms['dC'], histograms['dD']]})
    
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

    Parameters:
        hist_data: Dictionary {region: histogram} or list of 4 histograms (A, B, C, D).
        rootfileName: Name of the output ROOT file.
    """
    # Ensure the output directory exists
    output_dir = os.path.dirname(rootfileName)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Create the ROOT file
    output_file = ROOT.TFile(rootfileName, "RECREATE")

    try:
        if isinstance(hist_data, ROOT.TH1):  # Single histogram case
            folder = output_file.mkdir('A')
            folder.cd()
            hist_clone = hist_data.Clone()
            hist_clone.SetDirectory(0)  # Prevent ROOT from automatically managing this
            hist_clone.Write(hist_clone.GetName(), ROOT.TObject.kOverwrite)
        
        elif isinstance(hist_data, dict):  # Dictionary of histograms
            for region, hist in hist_data.items():
                folder = output_file.GetDirectory(region)
                if not folder:
                    folder = output_file.mkdir(region)
                folder.cd()
                
                hist_clone = hist.Clone()
                hist_clone.SetDirectory(0)  # Prevent memory issues
                hist_clone.Write(hist_clone.GetName(), ROOT.TObject.kOverwrite)
        
        elif isinstance(hist_data, list):  # List of histograms (A, B, C, D)
            if len(hist_data) != 4:
                raise ValueError("For list input, exactly 4 histograms are expected (for A, B, C, D).")
            
            regions = ["dA", "dB", "dC", "dD"]
            for region, hist in zip(regions, hist_data):
                folder = output_file.GetDirectory(region)
                if not folder:
                    folder = output_file.mkdir(region)
                folder.cd()
                
                hist_clone = hist.Clone()
                hist_clone.SetDirectory(0)  # Prevent memory issues
                hist_clone.Write(hist_clone.GetName(), ROOT.TObject.kOverwrite)

    finally:
        output_file.Write("", ROOT.TObject.kOverwrite)  # Ensure all histograms are saved
        output_file.Close()
        print(f"Histograms saved successfully to {rootfileName}.")
