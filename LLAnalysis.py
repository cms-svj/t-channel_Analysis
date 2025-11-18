import ROOT
import ROOTplotutils as pltutils
import optparse
import pandas as pd
import utils.DataSetInfo as info
import TFutils as TF
import plotStack
import utils.CMS_lumi as CMS_lumi

AddCMSText = True
# ROOT.TH1.SetDefaultSumw2()
# ROOT.TH2.SetDefaultSumw2()

def GetSVJbins():
    '''Change the SVJbin edges here'''
    SVJbins = {
                "0SVJ" : [0.5,350.0],
                "1SVJ" : [0.56,450.0],
                "2PSVJ" : [0.56,350.0],
    }
    return SVJbins

def METbinningDict():
    '''Get a dictonary for met bins'''
    return {
        "2016": {
            "1SVJ": [200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 5000.0],
            "2PSVJ": [200.0, 250.0, 300.0, 350.0, 400.0, 5000.0],
        },
        "2017": {
            "1SVJ": [200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 850.0, 5000.0],
            "2PSVJ": [200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 5000.0],
        },
        "2018": {
            "1SVJ": [200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 850.0, 5000.0],
            "2PSVJ": [200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 5000.0],
        },
    }

def getData(path, scale=1.0, year = "2018"):
    '''Uncomment the files that are to be included in the plots'''
    Data = [
        info.DataSetInfo(basedir=path, fileName=year+"_Data.root",        sys= -1.0, label="Data",        scale=scale),
    ]
    bgData = [
        info.DataSetInfo(basedir=path, fileName=year+"_ST.root",              label="Single top",              scale=scale, color=ROOT.TColor.GetColor("#5790fc")),
        info.DataSetInfo(basedir=path, fileName=year+"_TTJets.root",          label="t#bar{t}",                scale=scale, color=ROOT.TColor.GetColor("#f89c20")),
        info.DataSetInfo(basedir=path, fileName=year+"_ZJets.root",             label="Z#rightarrow#nu#nu+jets",    scale=scale, color=ROOT.TColor.GetColor("#e42536")),
        info.DataSetInfo(basedir=path, fileName=year+"_WJets.root",              label="W+jets",                    scale=scale, color=ROOT.TColor.GetColor("#964a8b")),
        info.DataSetInfo(basedir=path, fileName=year+"_QCD.root",               label="QCD",                        scale=scale, color=ROOT.TColor.GetColor("#9c9ca1")),
    ]
    sgData = [
        # TODO: check if there is a cms guidline for line color
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="baseline", scale=scale, color=ROOT.kBlack),
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",     label="mMed_600",  scale=scale, color=ROOT.kViolet+2),
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-800_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",     label="mMed_800",  scale=scale, color=ROOT.kRed),
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-1000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="mMed_1000", scale=scale, color=ROOT.kMagenta + 1),
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-1500_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="mMed_1500", scale=scale, color=ROOT.kGray+4),
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="mMed_3000", scale=scale, color=ROOT.kCyan),
        info.DataSetInfo(basedir=path, fileName=year+"_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="mMed_4000", scale=scale, color=ROOT.kBlue + 1),
        
    ]
    return Data, sgData, bgData

def plotABCDwithMETbins(dataList, ABCDHistoVar, maincut, SVJBins, plotOutputDir, xTitle="MET [GeV]",yTitle="Events",xmin=999.9, xmax = -999.9, isLogY = False, year="2018", isRatio=False,hemPeriod=False ):
    '''
    Plots either a stacked histogram for background with or without the signal (when data is empty) for MET binnings
    or a ratio plot between Data and background in case for Control region 
    '''
    ROOT.TH1.AddDirectory(False)
    firstpass = True
    Data, sgData, bgData = dataList
    if maincut == "_pre_":
        Data = None
    if isRatio and Data is None: 
        print("\033[31mWarning: Not passed data cannot make ratio plot, making normal plot instead .\033[0m") 
        isRatio = False
    
    bgHists = TF.GetABCDhistDict(bgData, )

    canvas = ROOT.TCanvas("canvas","Canvas with ratio pads", 1000,800)
    canvas.Divide(2,2)

    for i, region in enumerate(['A','B','C','D'], start=1):
        canvas.cd(i)


    


def plotABCD(dataList, ABCDHistoVar, maincut, SVJBins, plotOutputDir, xTitle="nSVJ",yTitle="Events",xmin=999.9, xmax = -999.9, isLogY = False, year="2018", isRatio=False,hemPeriod=False):
    '''
    Plots either a stacked histogram for background with or without the signal (when data is empty)
    or a ratio plot between tData and background histograms. 
    '''
    ROOT.TH1.AddDirectory(False)
    firstpass = True
    Data, sgData, bgData = dataList
    if maincut == "_pre_":
        Data = None
    # if maincut == "_lcr_pre_":
    #     sgData = None
    if isRatio and Data is None: 
        print("Not passed data cannot make ratio plot, making normal plot instead") 
        isRatio = False
    if isRatio: 
        c1 = ROOT.TCanvas( "c", "c", 800, 700)
        c1, pad1, pad2 = pltutils.createCanvasPads(c1,isLogY)
        pad1.cd()
    else: 
        c1 = ROOT.TCanvas( "c", "c", 800, 800)
        c1.cd()
        pltutils.SetupGPad(logY=isLogY)
    
    ROOT.gStyle.SetOptStat("")

    #TLegend 
    if isRatio:
        leg = pltutils.SetupLegend(NColumns=2)
    else:
        leg = pltutils.SetupLegend(NColumns=2, textSize=0.024)

    if bgData:
        bgDataMergedABCDDict = TF.GetABCDhistDict(bgData,ABCDHistoVar, maincut, SVJBins, isStack=True, merge=True)
        bgStackedHist, bgSummedHist = pltutils.StackedHistogram(bgDataMergedABCDDict)
        for bgDataHist in bgDataMergedABCDDict.keys():
            leg.AddEntry(bgDataMergedABCDDict[bgDataHist], bgDataHist.label_,"F")
        if firstpass:
            firstpass = False
            dummy = ROOT.TH1D("dummy", "dummy", len(SVJBins)*4, bgSummedHist.GetBinLowEdge(1), bgSummedHist.GetBinLowEdge(bgSummedHist.GetNbinsX()) + bgSummedHist.GetBinWidth(bgSummedHist.GetNbinsX()))
            # print(f"The dummy values are  - {bgSummedHist.GetBinLowEdge(1)}, {bgSummedHist.GetBinLowEdge(bgSummedHist.GetNbinsX())}, { bgSummedHist.GetBinWidth(bgSummedHist.GetNbinsX())}")
            ymax=10**9
            ymin=10
            lmax=10**9
            # Set bin labels for the dummy histogram to match bgSummedHist
            for bin_idx in range(1, bgSummedHist.GetNbinsX() + 1):
                bin_label = bgSummedHist.GetXaxis().GetBinLabel(bin_idx)
                dummy.GetXaxis().SetBinLabel(bin_idx, bin_label)
            pltutils.setupDummy(dummy,leg,"", xTitle, yTitle, isLogY, xmin, xmax, ymin, ymax, lmax, isRatio=isRatio)
            dummy.Draw("hist")
            
        bgStackedHist.Draw("hist F same")
        lines_upperpad = pltutils.AddVerticalLine(dummy, SVJBins, ymax = 10**5)
        for line in lines_upperpad:
            line.Draw("same")
    
    if sgData: 
        sgDataMergedABCDDict = TF.GetABCDhistDict(sgData, ABCDHistoVar, maincut, SVJBins, merge=True)
        linestylenumber = 0
        linestyle = [ROOT.kSolid,ROOT.kDashed,ROOT.kDotted]
        for sgDataHist in sgDataMergedABCDDict.keys():
            leg.AddEntry(sgDataMergedABCDDict[sgDataHist], sgDataHist.label_,"L")
            sgDataMergedABCDDict[sgDataHist].SetLineStyle(linestyle[linestylenumber%3])
            sgDataMergedABCDDict[sgDataHist].SetLineColor(sgDataHist.color)
            linestylenumber+=1
            sgDataMergedABCDDict[sgDataHist].SetLineWidth(3)
            sgDataMergedABCDDict[sgDataHist].Draw("hist same") 
            
    if Data is not None:
        for d in Data:
            DataHist = TF.GetABCDhist(d, ABCDHistoVar, maincut, SVJBins, merge=True)
            ROOT.gStyle.SetErrorX(0.)
            DataHist = pltutils.SetupDataStyle(DataHist)
            leg.AddEntry(DataHist, d.label_)
            DataHist.Draw("P same")
            
            if isRatio:
                pad2.cd()
                ratioHist = pltutils.RatioHistogram(DataHist,bgSummedHist)
                ratioHist = pltutils.SetupRatioStyle(ratioHist, xTitle, yTitle="Data/MC", yTitleSize=0.13)
                ratioHist.Draw("EX0P")
                lines_lowerpad = pltutils.AddVerticalLine(DataHist, SVJBins, ymax = 2)
                for line in lines_lowerpad:
                    line.Draw("same")
        pad1.cd()
    leg.Draw("same")
    
    # pltutils.AddVerticalLineForABCD(dummy, SVJBins)
    pltutils.AddLabelsForABCD(dummy, SVJBins,yloc=5*10**4)
    dummy.Draw("AXIS same")
    # pltutils.AddVerticalLineAtBinEnd(dummy, 3)
    if AddCMSText:
        pltutils.AddCMSLumiText(c1, year, isExtraText=True,hemPeriod=hemPeriod)
    c1.cd()
    c1.Update()
    # c1.RedrawAxis()
    ROOT.gPad.RedrawAxis()
    ROOT.gPad.RedrawAxis("G")
    SaveName = plotOutputDir+"/ABCDPlot_"+maincut
    c1.SaveAs(SaveName+".png")
    c1.Close()
    del c1, leg

def GetTFhistoAndPlot(ABCDhistDict_SR, ABCDhistDict_CR, SVJBins, xTitle="nSVJ", yTitle="Events", xmin=999.9, xmax=-999.9, year="2018", isLogY=False, saveName="TF.png",hemPeriod=False):
    ROOT.TH1.AddDirectory(False)
    skip_list = [f"{year}_QCD.root",f"{year}_ZJets.root",f"{year}_Data.root"]
    summed_SR = TF.SumHistograms(ABCDhistDict_SR, skip_list)
    summed_CR = TF.SumHistograms(ABCDhistDict_CR, skip_list)
    TFhist = pltutils.RatioHistogram(summed_SR,summed_CR)
    c1 = ROOT.TCanvas( "c", "c", 800, 700)
    c1, pad1, pad2 = pltutils.createCanvasPads(c1,isLogY)
    pad1.cd()
    ROOT.gStyle.SetOptStat("")
    leg = pltutils.SetupLegend(x1=0.4, NColumns=2)

    # make dummy to setup axis
    dummy = ROOT.TH1D("dummy", "dummy", len(SVJBins)*4, summed_SR.GetBinLowEdge(1), summed_SR.GetBinLowEdge(summed_SR.GetNbinsX()) + summed_SR.GetBinWidth(summed_SR.GetNbinsX()))
    # print(f"The dummy values are  - {summed_SR.GetBinLowEdge(1)}, {summed_SR.GetBinLowEdge(summed_SR.GetNbinsX())}, { summed_SR.GetBinWidth(summed_SR.GetNbinsX())}")
    ymax=10*9
    ymin=10
    lmax=10**9
    pltutils.setupDummy(dummy,leg,"", xTitle, yTitle, isLogY, xmin, xmax, ymin, ymax, lmax, isRatio=True)
    dummy.Draw("hist")
    
    pltutils.SetupLineHistStyle(summed_SR)
    pltutils.SetupLineHistStyle(summed_CR, color=ROOT.kRed)
    summed_SR.Draw("histe same")
    summed_CR.Draw("histe same")
    leg.AddEntry(summed_CR,"Control Region","L")
    leg.AddEntry(summed_SR,"Signal Region","L")

    lines_upperpad = pltutils.AddVerticalLine(summed_SR, SVJBins, ymax = 10**5)
    for line in lines_upperpad:
        line.Draw("same")
    
    pltutils.AddLabelsForABCD(summed_SR,SVJBins,yloc=5*10**4)
    leg.Draw("same")
    pad2.cd()
    TFhist = pltutils.SetupRatioStyle(TFhist, xTitle="nSVJ", yTitle="SR/CR", yTitleSize=0.13, ymax=2)
    TFhist.Draw("EX0P")
    lines_lowerpad = pltutils.AddVerticalLine(summed_SR, SVJBins, ymax = 2)
    for line in lines_lowerpad:
        line.Draw("same")
    if AddCMSText:
        pltutils.AddCMSLumiText(c1, year,isExtraText=True,hemPeriod=hemPeriod)
    
    # print(f"Summed SR Hist -  {pltutils.printBinContentAndError(summed_SR)}")
    # print(f"Summed CR Hist - {pltutils.printBinContentAndError(summed_CR)}")
    # print(f"TF Hist - {pltutils.printBinContentAndError(TFhist)}")

    c1.cd()
    c1.Update()
    # c1.RedrawAxis()
    ROOT.gPad.RedrawAxis()
    ROOT.gPad.RedrawAxis("G")
    c1.SaveAs(saveName+".png")
    c1.Close()
    del c1, leg

    return TFhist, summed_SR, summed_CR

def PlotValidation(TFhist, hist, expected_SR, MC_QCD_ZJets_SR, SVJBins,xTitle="nSVJ", yTitle="Events", xmin=999.9, xmax = -999.9, isLogY = False, year="2018", saveName="Validation", isDataCR=False,hemPeriod=False):
    '''
    validation is done by taking the ratio between the SR (Expected) and the Predicted SR.
    '''
    ROOT.TH1.AddDirectory(False)
    predicted_Data_SR = TF.Validation(hist,TFhist, MC_QCD_ZJets_SR)
    RatioHist = pltutils.RatioHistogram(expected_SR,predicted_Data_SR)
    c1 = ROOT.TCanvas( "c", "c", 800, 700)
    c1, pad1, pad2 = pltutils.createCanvasPads(c1,isLogY)
    pad1.cd()
    ROOT.gStyle.SetOptStat("")
    leg = pltutils.SetupLegend(x1=0.2, textSize=0.04)
    pltutils.SetupLineHistStyle(expected_SR)
    pltutils.SetupLineHistStyle(predicted_Data_SR, color=ROOT.kRed)
   
    # make dummy plot for axis
    dummy = ROOT.TH1D("dummy", "dummy", len(SVJBins)*4, expected_SR.GetBinLowEdge(1), expected_SR.GetBinLowEdge(expected_SR.GetNbinsX()) + expected_SR.GetBinWidth(expected_SR.GetNbinsX()))
    # print(f"The dummy values are  - {expected_SR.GetBinLowEdge(1)}, {expected_SR.GetBinLowEdge(expected_SR.GetNbinsX())}, { expected_SR.GetBinWidth(expected_SR.GetNbinsX())}")
    ymax=10**9
    ymin=10
    lmax=10**9
    pltutils.setupDummy(dummy,leg,"", xTitle, yTitle, isLogY, xmin, xmax, ymin, ymax, lmax, isRatio=True)
    dummy.Draw("hist")

    expected_SR.Draw("histe same")
    predicted_Data_SR.Draw("histe same")
    leg.AddEntry(expected_SR,"Expected SR (MC)","L")
    if isDataCR:
        leg.AddEntry(predicted_Data_SR,"Predicted SR ((Data_{CR} - MC^{QCD,ZJets}_{CR}) #times TF + MC^{QCD,ZJets}_{SR})","L")
    else:
        leg.AddEntry(predicted_Data_SR,"Predicted SR (CR #times TF + SR MC_{QCD,ZJets})","L")
    lines_upperpad = pltutils.AddVerticalLine(expected_SR, SVJBins, ymax = 10**5)
    for line in lines_upperpad:
        line.Draw("same")
    
    pltutils.AddLabelsForABCD(predicted_Data_SR,SVJBins,yloc=5*10**4)
    leg.Draw("same")
    pad2.cd()
    
    RatioHist = pltutils.SetupRatioStyle(RatioHist, xTitle="nSVJ", yTitle="Expected/Predicted", yTitleSize=0.1, ymin=0.6, ymax=1.4)
    RatioHist.Draw("EX0P")
    lines_lowerpad = pltutils.AddVerticalLine(expected_SR, SVJBins, ymin=0.5, ymax = 1.5)
    for line in lines_lowerpad:
        line.Draw("same")
    if AddCMSText:
        pltutils.AddCMSLumiText(c1, year,isExtraText=True,hemPeriod=hemPeriod)
    
    
    # print(f"Expected Hist -  {pltutils.printBinContentAndError(expected_SR)}")
    # print(f"Predicted Hist - {pltutils.printBinContentAndError(predicted_Data_SR)}")
    # print(f"Ratio Hist - {pltutils.printBinContentAndError(RatioHist)}")
    
    
    c1.cd()
    c1.Update()
    # c1.RedrawAxis()
    ROOT.gPad.RedrawAxis()
    ROOT.gPad.RedrawAxis("G")
    c1.SaveAs(saveName+".png")
    c1.Close()
    del c1, leg



def TFAnalysis(dataList, ABCDHistoVar, SRcut,CRcut, SVJBins,  plotOutputDir, xTitle="nSVJ",yTitle="Events",xmin=999.9, xmax = -999.9, isLogY = False, year="2018",hemPeriod=False):
    '''
    Calculate the transfer factor and make the plots for the Trandfer factor and if Data is available then make validation plots
    '''
    ROOT.TH1.AddDirectory(False)
    Data, bgData = dataList

    bgDataMergedABCDDict_SR = TF.GetABCDhistDict(bgData, ABCDHistoVar, SRcut, SVJBins, merge=True)
    bgDataMergedABCDDict_CR = TF.GetABCDhistDict(bgData, ABCDHistoVar, CRcut, SVJBins, merge=True)
    print("***************** Working on TF   ****************")
    TFhisto, LLinSR, LLinCR = GetTFhistoAndPlot(bgDataMergedABCDDict_SR, bgDataMergedABCDDict_CR, SVJBins,year = year, isLogY=isLogY, saveName=plotOutputDir+f"/TFplot_{SRcut}_{CRcut}",hemPeriod=hemPeriod)
    histogram_list = [TFhisto, LLinSR, LLinCR]
    for hist in bgDataMergedABCDDict_SR.values():
        histogram_list.append(hist)
    for hist in bgDataMergedABCDDict_CR.values():
        histogram_list.append(hist)
    for data in Data:
        DataHist_CR = TF.GetABCDhist(data, ABCDHistoVar, CRcut, SVJBins, merge=True)
        histogram_list.append(DataHist_CR)
        # print(f" ***********  Working on MC_QCD_ZJets_SR  ****************    skip list is - {year}_TTJets,{year}_WJets,{year}_ST") 
        MC_QCD_ZJets_SR = TF.SumHistograms(bgDataMergedABCDDict_SR,skipList=[f"{year}_TTJets.root",f"{year}_WJets.root",f"{year}_ST.root"])
        MC_QCD_ZJets_CR = TF.SumHistograms(bgDataMergedABCDDict_CR,skipList=[f"{year}_TTJets.root",f"{year}_WJets.root",f"{year}_ST.root"])
        DataHist_Sub_QCD_ZJets = DataHist_CR.Clone("Data CR - MC QCD ZJets")
        DataHist_Sub_QCD_ZJets.Add(MC_QCD_ZJets_CR, -1)
        histogram_list.append(MC_QCD_ZJets_SR)
        expected_SR = LLinSR.Clone("expected SR")
        expected_SR.Add(MC_QCD_ZJets_SR)
        # PlotValidation(TFhisto, DataHist_CR, expected_SR, MC_QCD_ZJets_SR, SVJBins, xTitle="Bin index", yTitle="Events", isLogY=isLogY, year=year, saveName=plotOutputDir+f"/Validation_{SRcut}_{CRcut}_usingDataCR", isDataCR=True,hemPeriod=hemPeriod)
        PlotValidation(TFhisto, LLinCR, expected_SR, MC_QCD_ZJets_SR, SVJBins, xTitle="Bin index", yTitle="Events", isLogY=isLogY, year=year, saveName=plotOutputDir+f"/Validation_{SRcut}_{CRcut}_usingCR",hemPeriod=hemPeriod)
        PlotValidation(TFhisto, DataHist_Sub_QCD_ZJets, expected_SR, MC_QCD_ZJets_SR, SVJBins, xTitle="Bin index", yTitle="Events", isLogY=isLogY, year=year, saveName=plotOutputDir+f"/Validation_{SRcut}_{CRcut}_usingDataCR_Sub_MC_QCD_ZJets", isDataCR=True,hemPeriod=hemPeriod)

    return histogram_list     

def SaveROOTFiles(DataList, ABCDHistoVar, SRcut, CRcut, SVJbin, plotOutputDir,year=2018, perSVJbin=False,isMETbinsCollapsed = True ,Metbinning=None):
    """
    Generate and save various ROOT files for analysis, including A, B, C, D histograms,
    summed backgrounds, transfer factors, LL estimates, and QCD estimates compatible with the StatInference Framework.
    
    Parameters:
        DataList: Tuple of (Data, sgData, bgData).
        ABCDHistoVar: Variable name for histograms.
        SRcut: String defining the SR cut.
        CRcut: String defining the CR cut.
        SVJbins: Dictionary containing SVJ bin definitions.
        plotOutputDir: Directory where the ROOT files will be saved.
    """
    Data, sgData, bgData = DataList
    print(SVJbin[0])
    print(f"year - {year}, SVJbin - {SVJbin[0]}")
    Metbinning = METbinningDict().get(year, {}).get(SVJbin[0], [])

    # Generate and save individual A, B, C, D histograms for Data, Signal, and Background
    for dataset, cutname in [(Data, "Data"), (sgData, "Signal"), (bgData, "Background")]:  # TODO: cutname is not required here, fix this
        for data in dataset:
            hist_dict = TF.GetABCDhistDict([data], ABCDHistoVar, SRcut, SVJbin, isStack=False, merge=False,perSVJbin=perSVJbin,isMETbinsCollapsed=isMETbinsCollapsed, Metbinning=Metbinning)
            print(hist_dict)
            if cutname != "Data":
                TF.SaveHistogramsToFile(hist_dict[data], f"{plotOutputDir}/{data.fileName.replace('.root', f'_SR.root')}")
            hist_dict = TF.GetABCDhistDict([data], ABCDHistoVar, CRcut, SVJbin, isStack=False, merge=False,perSVJbin=perSVJbin,isMETbinsCollapsed=isMETbinsCollapsed, Metbinning=Metbinning)
            if cutname != "Signal":
                TF.SaveHistogramsToFile(hist_dict[data], f"{plotOutputDir}/{data.fileName.replace('.root', f'_CR.root')}")
    #         del hist_dict

    # # Sum all bgData histograms for total background
    # total_bg_SR = TF.SumABCDhistList(TF.GetABCDhistDict(bgData, ABCDHistoVar, SRcut, SVJbin, isStack=False, merge=False,perSVJbin=perSVJbin,isMETbinsCollapsed=isMETbinsCollapsed, Metbinning=Metbinning))
    # TF.SaveHistogramsToFile(total_bg_SR, f"{plotOutputDir}/{year}_total_bkg_SR.root")
    # total_bg_CR = TF.SumABCDhistList(TF.GetABCDhistDict(bgData, ABCDHistoVar, CRcut, SVJbin, isStack=False, merge=False,perSVJbin=perSVJbin,isMETbinsCollapsed=isMETbinsCollapsed, Metbinning=Metbinning))
    # TF.SaveHistogramsToFile(total_bg_CR, f"{plotOutputDir}/{year}_total_bkg_CR.root")

    # # Sum only lost lepton backgrounds (WJets, TTJets, ST) for mearged bins
    # skip_list = [f"{year}_QCD.root",f"{year}_ZJets.root",f"{year}_Data.root"]
    # ll_bkg_SR_merged = TF.SumABCDhistList(TF.GetABCDhistDict(bgData, ABCDHistoVar, SRcut, SVJbin, isStack=False, merge=False,perSVJbin=perSVJbin,isMETbinsCollapsed=True), skipList=skip_list)
    # # ll_bkg_SR.SetTitle(f"LL bkg in SR")
    # TF.SaveHistogramsToFile(ll_bkg_SR_merged, f"{plotOutputDir}/{year}_LL_bkg_SR_merged.root")
    # ll_bkg_CR_merged = TF.SumABCDhistList(TF.GetABCDhistDict(bgData, ABCDHistoVar, CRcut, SVJbin, isStack=False, merge=False,perSVJbin=perSVJbin,isMETbinsCollapsed=True), skipList=skip_list)
    # # ll_bkg_SR.SetTitle(f"LL bkg in CR")
    # TF.SaveHistogramsToFile(ll_bkg_CR_merged, f"{plotOutputDir}/{year}_LL_bkg_CR_merged.root")

    # lost lepton backgrounds without the mearged bins
    # ll_bkg_SR = TF.SumABCDhistList(TF.GetABCDhistDict(bgData, ABCDHistoVar, SRcut, SVJbin, isStack=False, merge=False,perSVJbin=perSVJbin,isMETbinsCollapsed=isMETbinsCollapsed, Metbinning=Metbinning), skipList=skip_list)
    # ll_bkg_SR.SetTitle(f"LL bkg in SR")
    # print(ll_bkg_SR)
    # TF.SaveHistogramsToFile(ll_bkg_SR, f"{plotOutputDir}/{year}_LL_bkg_SR.root")
    # ll_bkg_CR = TF.SumABCDhistList(TF.GetABCDhistDict(bgData, ABCDHistoVar, CRcut, SVJbin, isStack=False, merge=False,perSVJbin=perSVJbin,isMETbinsCollapsed=isMETbinsCollapsed, Metbinning=Metbinning), skipList=skip_list)
    # # ll_bkg_SR.SetTitle(f"LL bkg in CR")
    # TF.SaveHistogramsToFile(ll_bkg_CR, f"{plotOutputDir}/{year}_LL_bkg_CR.root")


    # # print(f"ll_bkg_SR - {ll_bkg_SR} ")
    # # Compute Transfer Factors
    # tf_hist = []
    # for i, region in enumerate(['A', 'B', 'C', 'D']):
    #     tf_hist.append(ll_bkg_SR_merged[i].Clone(f"h_TF_{region}"))
    #     tf_hist[i].SetTitle(f"Transfer Factor in {region}")
    #     tf_hist[i].Divide(ll_bkg_CR_merged[i])
    # TF.SaveHistogramsToFile(tf_hist, f"{plotOutputDir}/{year}_LL_transfer_factors.root")

    # # Removing MC QCD and Zjets from the Data CR
    # skip_list_qcd_zjets = [f"{year}_TTJets.root",f"{year}_WJets.root",f"{year}_ST.root",f"{year}_Data.root"]
    # MC_QCD_ZJets_CR = TF.SumABCDhistList(TF.GetABCDhistDict(bgData, ABCDHistoVar, CRcut, SVJbin, isStack=False, merge=False,perSVJbin=perSVJbin,isMETbinsCollapsed=isMETbinsCollapsed, Metbinning=Metbinning), skipList=skip_list_qcd_zjets)
    # TF.SaveHistogramsToFile(MC_QCD_ZJets_CR, f"{plotOutputDir}/{year}_MC_QCD_ZJets_CR.root")

    # Data_CR = TF.GetABCDhistPerSVJBin(Data[0], ABCDHistoVar, CRcut, SVJbin,isMETbinsCollapsed=isMETbinsCollapsed, Metbinning=Metbinning)
    
    # data_sub_qcd_ZJets = []
    # for i,region in enumerate(['A', 'B', 'C', 'D']):
    #     data_sub_qcd_ZJets.append(Data_CR[region].Clone(f"h_Data_sub_QCD_ZJets_CR_{region}"))
    #     data_sub_qcd_ZJets[i].SetTitle(f"Data CR - MC QCD ZJets in CR for {region}")
    #     data_sub_qcd_ZJets[i].Add(MC_QCD_ZJets_CR[i],-1)
    # TF.SaveHistogramsToFile(data_sub_qcd_ZJets, f"{plotOutputDir}/{year}_Data_CR_subQCDZJets.root")


    # # Compute LL Estimate in SR
    # ll_estimate_in_sr = []
    
    # for i, region in enumerate(['A', 'B', 'C', 'D']):
    #     ll_estimate_in_sr.append(data_sub_qcd_ZJets[i].Clone(f"h_LL_Estimate_{region}"))
    #     ll_estimate_in_sr[i].SetTitle(f"LL estimate in SR for {region}")
    #     ll_estimate_in_sr[i].Multiply(tf_hist[i])
    # TF.SaveHistogramsToFile(ll_estimate_in_sr, f"{plotOutputDir}/{year}_LL_Estimate_in_SR.root")


    # ll_estimate_in_sr_data_cr_no_qcdZjets_subtraction = []
    
    # for i, region in enumerate(['A', 'B', 'C', 'D']):
    #     ll_estimate_in_sr_data_cr_no_qcdZjets_subtraction.append(Data_CR[region].Clone(f"h_LL_Estimate_in_sr_data_cr_no_qcdzjets_subtraction{region}"))
    #     ll_estimate_in_sr_data_cr_no_qcdZjets_subtraction[i].SetTitle(f"LL estimate in SR for {region} where QCD and ZJets is not subtracted from Data CR")
    #     ll_estimate_in_sr_data_cr_no_qcdZjets_subtraction[i].Multiply(tf_hist[i])
    # TF.SaveHistogramsToFile(ll_estimate_in_sr_data_cr_no_qcdZjets_subtraction, f"{plotOutputDir}/{year}_LL_Estimate_in_SR_DataCR_without_qcdzjets_subtraction.root")

    # # Subtract LL estimate from total background in SR
    # total_bkg_minus_ll = []
    # for i, region in enumerate(['A', 'B', 'C', 'D']):
    #     total_bkg_minus_ll.append(total_bg_SR[i].Clone(f"h_total_minus_LL_{region}"))
    #     total_bkg_minus_ll[i].SetTitle(f"Total bkg in SR - LL estimate in SR for {region}")
    #     total_bkg_minus_ll[i].Add(ll_estimate_in_sr[i], -1)
    # TF.SaveHistogramsToFile(total_bkg_minus_ll, f"{plotOutputDir}/{year}_LL_estimate_subtracted_total_bkg_in_SR.root")

    # # Subtract LL estimate and ZJets from total background in SR for QCD Estimate
    # qcd_estimate_in_sr = []
    # zdata = [data for data in bgData if "ZJets.root" in data.fileName]
    # # print(f"z data - {zdata[0].fileName}")
    # # Get ZJets histograms (using GetABCDhist)
    # if perSVJbin:
    #     zjets_hist = TF.GetABCDhistPerSVJBin(zdata[0],ABCDHistoVar,SRcut,SVJbin,isMETbinsCollapsed=isMETbinsCollapsed, Metbinning=Metbinning)
    # else:
    #     zjets_hist = TF.GetABCDhist(zdata[0],ABCDHistoVar,SRcut,SVJbin,isStack=False)

    # # Debugging: Check the structure of zjets_hist
    # # print(f"ZJets histograms: {zjets_hist}")

    # # Compute QCD Estimate in SR
    # for i, region in enumerate(['A', 'B', 'C', 'D']):
    #     qcd_estimate_in_sr.append(total_bkg_minus_ll[i].Clone(f"h_QCD_Estimate_SR_{region}"))
    #     qcd_estimate_in_sr[i].SetTitle(f"QCD estimate in SR for region {region}")
    #     qcd_estimate_in_sr[i].Add(zjets_hist[region], -1)  # Indexing zjets_hist by i
    # TF.SaveHistogramsToFile(qcd_estimate_in_sr, f"{plotOutputDir}/{year}_QCD_Estimate_in_SR.root")

    # # Create a histogram for C/D
    # cd_ratio_hist = qcd_estimate_in_sr[2].Clone("h_QCD_C_over_D")  # Clone C histogram
    # cd_ratio_hist.Divide(qcd_estimate_in_sr[3])  # Divide by D histogram
    # cd_ratio_hist.SetTitle("C/D ratio for QCD estimate in SR")

    # # Save the C/D histogram
    # TF.SaveHistogramsToFile(cd_ratio_hist, f"{plotOutputDir}/{year}_QCD_C_over_D.root")

    # # Compute QCD estimate in region A: C/D * B
    # qcd_in_A_hist = qcd_estimate_in_sr[1].Clone("h_QCD_in_A")  # Clone B histogram
    # qcd_in_A_hist.Multiply(cd_ratio_hist)  # Multiply by C/D
    # qcd_in_A_hist.SetTitle("QCD estimate in region A using ABCD method")

    # # Save the QCD in A histogram
    # TF.SaveHistogramsToFile(qcd_in_A_hist, f"{plotOutputDir}/{year}_QCD_in_A.root")

def checks(DataList, ABCDHistoVar, SRcut, CRcut, SVJbins, plotOutputDir,year=2018, perSVJbin=False):
    Data, sgData, bgData = DataList

    for d in bgData:
        integralConsistency = TF.checkIntegralConsistency(d, ABCDHistoVar, CRcut, SVJbins)
        # print(f"Integral consistency  - {integralConsistency}")


def metbinningChecks(DataList, ABCDHistoVar, SRcut, CRcut, SVJbin, plotOutputDir,year=2018, perSVJbin=False):
    Data, sgData, bgData = DataList
    print(SVJbin[0])
    print(f"year - {year}, SVJbin - {SVJbin[0]}")
    metbins = METbinningDict().get(year, {}).get(SVJbin[0], [])
    
    # Generate and save individual A, B, C, D histograms for Data, Signal, and Background
    for dataset, cutname in [(Data, "Data"), (sgData, "Signal"), (bgData, "Background")]:  # TODO: cutname is not required here, fix this
        for data in dataset:
            hist_dict = TF.GetABCDhistDict([data], ABCDHistoVar, SRcut, SVJbin, isStack=False, merge=False,perSVJbin=perSVJbin,isMETbinsCollapsed = False,Metbinning=metbins)
            print(hist_dict)
            TF.SaveHistogramsToFile(hist_dict[data], f"{plotOutputDir}/{data.fileName.replace('.root', f'_SR.root')}")
            # hist_dict = TF.GetABCDhistDict([data], ABCDHistoVar, CRcut, SVJbins, isStack=False, merge=False,perSVJbin=perSVJbin)
            # TF.SaveHistogramsToFile(hist_dict[data], f"{plotOutputDir}/{data.fileName.replace('.root', f'_CR.root')}")
            del hist_dict    

def main():
    parser = optparse.OptionParser("usage: %prog [options]\n")
    parser.add_option('-d', '--dataset',    dest='dataset',                    default='testHadd_11242020',    help='dataset')
    parser.add_option('-y',                 dest='year',       type='string',  default='2018',                 help="Can pass in the run year")
    parser.add_option('-o',                 dest='outputdir',  type='string',                                  help="Output folder name")
    parser.add_option(    '--hemPeriod',  dest='hemPeriod', type=str, default=False,  help='HEM period (PreHEM or PostHEM), default includes entire sample',)
    options, args = parser.parse_args()
    print("Parsed options:", options)

    # outROOTfileName = options.outfile
    year = options.year
    hemPeriod = options.hemPeriod
    SVJbins = GetSVJbins()
    ABCDhistoVars = ["METvsDNN"]
    ABCDFolderName = "ABCD"
    SRCut = "_pre_"
    CRCuts = ["_lcr_pre_"]#,"_cr_muon_","_cr_electron_"]
    maincuts = [SRCut] + CRCuts
    Data, sgData, bgData = getData( options.dataset + "/", 1.0, year)
    # print(f"Data = {Data[0].fileName} \n bgData = {bgData[0].fileName}")
    
    histograms = []

    if options.outputdir:
        plotOutDir = "outputPlots/{}".format(options.outputdir)
    else: 
        plotOutDir = "outputPlots/{}".format(options.dataset)

    for histName in ABCDhistoVars:
        for maincut in maincuts:
            pltutils.makeDirs(plotOutDir,maincut,ABCDFolderName)
            plotABCDdir = plotOutDir+'/'+ABCDFolderName+'/'+maincut[1:]
            plotABCD((Data, sgData, bgData), "h_"+histName, maincut, SVJbins, plotOutputDir=plotABCDdir, isLogY = True, year = year, isRatio= True, hemPeriod=hemPeriod)


        for CRcut in CRCuts:
            pltutils.makeDirs(plotOutDir,CRcut, "TransferFactors")
            plotTFdir= plotOutDir+'/'+"TransferFactors"+'/'+CRcut[1:]
            TFhistos = TFAnalysis((Data, bgData), "h_"+histName, SRCut, CRcut, SVJbins, plotOutputDir=plotTFdir, isLogY = True, year=year, hemPeriod=hemPeriod)

    # for SVJ in SVJbins.items():
    #     key, values = SVJ
    #     for histName in ABCDhistoVars:
    #         for CRCut in CRCuts:
    #             pltutils.makeDirs(plotOutDir,CRCut, f"{SVJ[0]}")
    #             Rootfiles = plotOutDir+f"/{SVJ[0]}"
    #             # checks((Data,sgData,bgData), "h_"+histName, SRCut,CRCut, SVJ, plotOutputDir=Rootfiles, year = year, perSVJbin=True)
    #             # metbinningChecks((Data,sgData,bgData), "h_"+histName, SRCut,CRCut, SVJ, plotOutputDir=Rootfiles, year = year, perSVJbin=True)

    #             SaveROOTFiles((Data,sgData,bgData), "h_"+histName, SRCut,CRCut, SVJ, plotOutputDir=Rootfiles, year = year, perSVJbin=True,isMETbinsCollapsed=False)



    

if __name__ == '__main__':
    main()


# TODO: Add a function that reads the SVJbin edges
# TODO: Add a function to make table of number of events in the ABCD region, TF, Validation
# TODO: Add a function to make a table in terms of SVJbins vs bkgs.
# TODO: Remove the 0 and 2 from the ratio plot pad 2
# TODO: Redraw axis properly for the dummy.
# TODO: Remove the grid from the pad1. 
