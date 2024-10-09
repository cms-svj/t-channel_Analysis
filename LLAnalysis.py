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

def getData(path, scale=1.0, year = "2018"):
    '''Uncomment the files that are to be included in the plots'''
    Data = [
        info.DataSetInfo(basedir=path, fileName=year+"_Data.root",        sys= -1.0, label="Data",        scale=scale),
    ]
    bgData = [
        info.DataSetInfo(basedir=path, fileName=year+"_ST.root",              label="Single top",              scale=scale, color=(ROOT.kRed + 1)),
        info.DataSetInfo(basedir=path, fileName=year+"_TTJets.root",          label="t#bar{t}",                scale=scale, color=(ROOT.kBlue - 6)),
        info.DataSetInfo(basedir=path, fileName=year+"_ZJets.root",             label="Z#rightarrow#nu#nu+jets",    scale=scale, color=(ROOT.kGray + 1)),
        info.DataSetInfo(basedir=path, fileName=year+"_WJets.root",              label="W+jets",                    scale=scale, color=(ROOT.kYellow + 1)),
        info.DataSetInfo(basedir=path, fileName=year+"_QCD.root",               label="QCD",                        scale=scale, color=(ROOT.kGreen + 1)),
    ]
    sgData = [

        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="baseline", scale=scale, color=ROOT.kBlack),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-600.root",     label="mMed_600",  scale=scale, color=ROOT.kViolet+2),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-800.root",     label="mMed_800",  scale=scale, color=ROOT.kRed),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-1000.root",    label="mMed_1000", scale=scale, color=ROOT.kMagenta + 1),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-1500.root",    label="mMed_1500", scale=scale, color=ROOT.kGray+4),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-3000.root",    label="mMed_3000", scale=scale, color=ROOT.kCyan),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-4000.root",    label="mMed_4000", scale=scale, color=ROOT.kBlue + 1),
        
    ]
    return Data, sgData, bgData

def plotABCD(dataList, ABCDHistoVar, maincut, SVJBins, plotOutputDir, xTitle="nSVJ",yTitle="Events",xmin=999.9, xmax = -999.9, isLogY = False, year="2018", isRatio=False):
    '''
    Plots either a stacked histogram for background with or without the signal (when data is empty)
    or a ratio plot between tData and background histograms. 
    '''
    ROOT.TH1.AddDirectory(False)
    firstpass = True
    Data, sgData, bgData = dataList
    if maincut == "_pre_":
        Data = None
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
            print(f"The dummy values are  - {bgSummedHist.GetBinLowEdge(1)}, {bgSummedHist.GetBinLowEdge(bgSummedHist.GetNbinsX())}, { bgSummedHist.GetBinWidth(bgSummedHist.GetNbinsX())}")
            ymax=10**11
            ymin=10
            lmax=10**11
            pltutils.setupDummy(dummy,leg,"", xTitle, yTitle, isLogY, xmin, xmax, ymin, ymax, lmax, isRatio=isRatio)
            dummy.Draw("hist")
            
        bgStackedHist.Draw("hist F same")
        lines_upperpad = pltutils.AddVerticalLine(dummy, SVJBins, ymax = 10**7)
        for line in lines_upperpad:
            line.Draw("same")
    
    if sgData: 
        sgDataMergedABCDDict = TF.GetABCDhistDict(sgData, ABCDHistoVar, maincut, SVJBins, merge=True)
        linestylenumber = 0
        linestyle = [ROOT.kSolid,ROOT.kDashed,ROOT.kDotted]
        for sgDataHist in sgDataMergedABCDDict.keys():
            leg.AddEntry(sgDataMergedABCDDict[sgDataHist], sgDataHist.label_,"L")
            sgDataHist.SetLineStyle(linestyle[linestylenumber%3])
            linestylenumber+=1
            sgDataHist.SetLineWidth(3)
            sgDataHist.Draw("hist same") 
            
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
    pltutils.AddLabelsForABCD(dummy, SVJBins,yloc=10**6)
    dummy.Draw("AXIS same")
    # pltutils.AddVerticalLineAtBinEnd(dummy, 3)
    if AddCMSText:
        pltutils.AddCMSLumiText(c1, year, isExtraText=True)
    c1.cd()
    c1.Update()
    # c1.RedrawAxis()
    ROOT.gPad.RedrawAxis()
    ROOT.gPad.RedrawAxis("G")
    SaveName = plotOutputDir+"/ABCDPlot_"+maincut
    c1.SaveAs(SaveName+".png")
    c1.Close()
    del c1, leg

def GetTFhistoAndPlot(ABCDhistDict_SR, ABCDhistDict_CR, SVJBins, xTitle="nSVJ", yTitle="Events", xmin=999.9, xmax=-999.9, year="2018", isLogY=False, saveName="TF.png"):
    ROOT.TH1.AddDirectory(False)
    skip_list = [f"{year}_QCD.root",f"{year}_ZJets.root",f"{year}_Data.root"]
    summed_SR = TF.SumHistograms(ABCDhistDict_SR, skip_list)
    summed_CR = TF.SumHistograms(ABCDhistDict_CR, skip_list)
    TFhist = pltutils.RatioHistogram(summed_SR,summed_CR)
    c1 = ROOT.TCanvas( "c", "c", 800, 700)
    c1, pad1, pad2 = pltutils.createCanvasPads(c1,isLogY)
    pad1.cd()
    ROOT.gStyle.SetOptStat("")
    leg = pltutils.SetupLegend(x1=0.7)

    # make dummy to setup axis
    dummy = ROOT.TH1D("dummy", "dummy", len(SVJBins)*4, summed_SR.GetBinLowEdge(1), summed_SR.GetBinLowEdge(summed_SR.GetNbinsX()) + summed_SR.GetBinWidth(summed_SR.GetNbinsX()))
    print(f"The dummy values are  - {summed_SR.GetBinLowEdge(1)}, {summed_SR.GetBinLowEdge(summed_SR.GetNbinsX())}, { summed_SR.GetBinWidth(summed_SR.GetNbinsX())}")
    ymax=10**10
    ymin=10
    lmax=10**11
    pltutils.setupDummy(dummy,leg,"", xTitle, yTitle, isLogY, xmin, xmax, ymin, ymax, lmax, isRatio=True)
    dummy.Draw("hist")
    
    pltutils.SetupLineHistStyle(summed_SR)
    pltutils.SetupLineHistStyle(summed_CR, color=ROOT.kRed)
    summed_SR.Draw("histe same")
    summed_CR.Draw("histe same")
    leg.AddEntry(summed_CR,"Control Region","L")
    leg.AddEntry(summed_SR,"Signal Region","L")

    lines_upperpad = pltutils.AddVerticalLine(summed_SR, SVJBins, ymax = 10**7)
    for line in lines_upperpad:
        line.Draw("same")
    
    pltutils.AddLabelsForABCD(summed_SR,SVJBins,yloc=10**6)
    leg.Draw("same")
    pad2.cd()
    TFhist = pltutils.SetupRatioStyle(TFhist, xTitle="nSVJ", yTitle="SR/CR", yTitleSize=0.13, ymax=2)
    TFhist.Draw("EX0P")
    lines_lowerpad = pltutils.AddVerticalLine(summed_SR, SVJBins, ymax = 2)
    for line in lines_lowerpad:
        line.Draw("same")
    if AddCMSText:
        pltutils.AddCMSLumiText(c1, year,isExtraText=True)
    
    print(f"Summed SR Hist -  {pltutils.printBinContentAndError(summed_SR)}")
    print(f"Summed CR Hist - {pltutils.printBinContentAndError(summed_CR)}")
    print(f"TF Hist - {pltutils.printBinContentAndError(TFhist)}")

    c1.cd()
    c1.Update()
    # c1.RedrawAxis()
    ROOT.gPad.RedrawAxis()
    ROOT.gPad.RedrawAxis("G")
    c1.SaveAs(saveName+".png")
    c1.Close()
    del c1, leg

    return TFhist, summed_SR, summed_CR

def PlotValidation(TFhist, Data_CR, expected_SR, MC_QCD_ZJets_SR, SVJBins,xTitle="nSVJ", yTitle="Events", xmin=999.9, xmax = -999.9, isLogY = False, year="2018", saveName="Validation.png"):
    '''
    validation is done by taking the ratio between the SR (Expected) and the Predicted SR.
    '''
    ROOT.TH1.AddDirectory(False)
    predicted_Data_SR = TF.Validation(Data_CR,TFhist, MC_QCD_ZJets_SR)
    RatioHist = pltutils.RatioHistogram(expected_SR,predicted_Data_SR)
    c1 = ROOT.TCanvas( "c", "c", 800, 700)
    c1, pad1, pad2 = pltutils.createCanvasPads(c1,isLogY)
    pad1.cd()
    ROOT.gStyle.SetOptStat("")
    leg = pltutils.SetupLegend(x1=0.2)
    pltutils.SetupLineHistStyle(expected_SR)
    pltutils.SetupLineHistStyle(predicted_Data_SR, color=ROOT.kRed)
   
    # make dummy plot for axis
    dummy = ROOT.TH1D("dummy", "dummy", len(SVJBins)*4, expected_SR.GetBinLowEdge(1), expected_SR.GetBinLowEdge(expected_SR.GetNbinsX()) + expected_SR.GetBinWidth(expected_SR.GetNbinsX()))
    print(f"The dummy values are  - {expected_SR.GetBinLowEdge(1)}, {expected_SR.GetBinLowEdge(expected_SR.GetNbinsX())}, { expected_SR.GetBinWidth(expected_SR.GetNbinsX())}")
    ymax=10**12
    ymin=10
    lmax=10**11
    pltutils.setupDummy(dummy,leg,"", xTitle, yTitle, isLogY, xmin, xmax, ymin, ymax, lmax, isRatio=True)
    dummy.Draw("hist")

    expected_SR.Draw("histe same")
    predicted_Data_SR.Draw("histe same")
    leg.AddEntry(expected_SR,"Expected SR (MC)","L")
    leg.AddEntry(predicted_Data_SR,"Predicted SR (CR Data #times TF + SR MC_{QCD,ZJets})","L")

    lines_upperpad = pltutils.AddVerticalLine(expected_SR, SVJBins, ymax = 10**7)
    for line in lines_upperpad:
        line.Draw("same")
    
    pltutils.AddLabelsForABCD(predicted_Data_SR,SVJBins,yloc=10**6)
    leg.Draw("same")
    pad2.cd()
    
    RatioHist = pltutils.SetupRatioStyle(RatioHist, xTitle="nSVJ", yTitle="Expected/Predicted", yTitleSize=0.1)
    RatioHist.Draw("EX0P")
    lines_lowerpad = pltutils.AddVerticalLine(expected_SR, SVJBins, ymax = 2)
    for line in lines_lowerpad:
        line.Draw("same")
    if AddCMSText:
        pltutils.AddCMSLumiText(c1, year,isExtraText=True)
    
    
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



def TFAnalysis(dataList, ABCDHistoVar, SRcut,CRcut, SVJBins,  plotOutputDir, xTitle="nSVJ",yTitle="Events",xmin=999.9, xmax = -999.9, isLogY = False, year="2018"):
    '''
    Calculate the transfer factor and make the plots for the Trandfer factor and if Data is available then make validation plots
    '''
    ROOT.TH1.AddDirectory(False)
    Data, bgData = dataList

    bgDataMergedABCDDict_SR = TF.GetABCDhistDict(bgData, ABCDHistoVar, SRcut, SVJBins, merge=True)
    bgDataMergedABCDDict_CR = TF.GetABCDhistDict(bgData, ABCDHistoVar, CRcut, SVJBins, merge=True)

    TFhisto, SR, CR = GetTFhistoAndPlot(bgDataMergedABCDDict_SR, bgDataMergedABCDDict_CR, SVJBins, year, isLogY=isLogY, saveName=plotOutputDir+f"/TFplot_{SRcut}_{CRcut}")

    for data in Data:
        DataHist_CR = TF.GetABCDhist(data, ABCDHistoVar, CRcut, SVJBins, merge=True)
        MC_QCD_ZJets_SR = TF.SumHistograms(bgDataMergedABCDDict_SR,skipList=[f"{year}_TTJets",f"{year}_WJets",f"{year}_ST"])
        PlotValidation(TFhisto, DataHist_CR, SR, MC_QCD_ZJets_SR, SVJBins, xTitle="Bin index", yTitle="Events", isLogY=isLogY, year=year, saveName=plotOutputDir+f"/Validation_{SRcut}_{CRcut}")
        


def main():
    parser = optparse.OptionParser("usage: %prog [options]\n")
    parser.add_option('-b',                 dest='isNormBkg',  action="store_true",                            help="Normalized Background and Signal plots")
    parser.add_option('-d', '--dataset',    dest='dataset',                    default='testHadd_11242020',    help='dataset')
    parser.add_option('-s',                 dest='onlySig',    action="store_true",                            help="Plot only signals")
    parser.add_option('-y',                 dest='year',       type='string',  default='2018',                 help="Can pass in the run year")
    parser.add_option('-o',                 dest='outputdir',  type='string',                                  help="Output folder name")
    parser.add_option('-w',                 dest='scenario',  type='string',   default='d0_w7p0i0',                               help="Scenario")
    parser.add_option(    '--hemPeriod',  dest='hemPeriod', type=str, default=False,  help='HEM period (PreHEM or PostHEM), default includes entire sample',)
    options, args = parser.parse_args()
    print("Parsed options:", options)

    year = options.year
    SVJbins = GetSVJbins()
    ABCDhistoVars = ["METvsDNN"]
    ABCDFolderName = "ABCD"
    SRCut = "_pre_"
    CRCuts = ["_lcr_pre_"]#,"_cr_muon_","_cr_electron_"]
    maincuts = [SRCut] + CRCuts
    Data, sgData, bgData = getData( options.dataset + "/", 1.0, year)
    print(f"Data = {Data[0].fileName} \n bgData = {bgData[0].fileName}")
    
    if options.outputdir:
        plotOutDir = "outputPlots/{}".format(options.outputdir)
    else: 
        plotOutDir = "outputPlots/{}".format(options.dataset)

    for histName in ABCDhistoVars:
        for maincut in maincuts:
            pltutils.makeDirs(plotOutDir,maincut,ABCDFolderName)
            plotABCDdir = plotOutDir+'/'+ABCDFolderName+'/'+maincut[1:]
            plotABCD((Data, sgData, bgData), "h_"+histName, maincut, SVJbins, plotOutputDir=plotABCDdir, isLogY = True, year = year, isRatio= True)


        for CRcut in CRCuts:
            pltutils.makeDirs(plotOutDir,CRcut, "TransferFactors")
            plotTFdir= plotOutDir+'/'+"TransferFactors"+'/'+CRcut[1:]
            TFAnalysis((Data, bgData), "h_"+histName, SRCut, CRcut, SVJbins, plotOutputDir=plotTFdir, isLogY = True, year=year)

if __name__ == '__main__':
    main()


    # TODO: Add a function that reads the SVJbin edges
# TODO: Add a function to make table of number of events in the ABCD region, TF, Validation
# TODO: Add a function to make a table in terms of SVJbins vs bkgs.
# TODO: Check the working with signal for all the plots.
# TODO: Remove the 0 and 2 from the ratio plot pad 2
# TODO: Redraw axis properly for the dummy.
    # TODO: Remove the grid from the pad1. 
    