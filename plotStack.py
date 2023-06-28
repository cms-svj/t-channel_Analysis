import ROOT
ROOT.gROOT.SetBatch(True)
import math
import utils.DataSetInfo as info
import optparse
import copy
import os
from array import array
import numpy as np
import utils.CMS_lumi as CMS_lumi
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
from utils.var import var as vars
# from utils.variables import variables as vars

mpl.rc("font", family="serif", size=15)

# color definition, needed for python plots for consistent color scheme
colorDict = {
    (ROOT.kGray + 1): "#999999",
    (ROOT.kYellow + 1): "#cccc19",
    (ROOT.kBlue - 6): "#6666cc",
    (ROOT.kGreen + 1): "#4bd42d",
    (ROOT.kMagenta + 1): "#cd2bcc",
    (ROOT.kOrange+2): "#cc660d",
    (ROOT.kRed): "#f2231b",
    (ROOT.kViolet-1): "#9a27cc",
    (ROOT.kOrange+3): "#663303",
    (ROOT.kPink+1): "#f69acc",
    (ROOT.kGray): "#cccccc",
    (ROOT.kYellow+2): "#999910",
    (ROOT.kCyan): "#5efdff",
}

# signal vs. background figure of merit
def fom(S,B):
    return np.sqrt(2 * ( (S+B) * np.log(1+S/B) - S) )

def signif(S,B):
    return S/(np.sqrt( B + (0.3*B)**2 ))

def getLabel(label):
    si = label.find("(")
    sourceLabel = label[:si-1]
    return sourceLabel

def makeDirs(plotOutDir,cut,plotType):
    if not os.path.exists(plotOutDir+"/"+plotType+"/"+cut[1:]):
        os.makedirs(plotOutDir+"/"+plotType+"/"+cut[1:])

def find_nearest(trgeff, target):
    trgeff = np.asarray(trgeff)
    idx = (np.abs(trgeff - target)).argmin()
    return int(trgeff[idx])

def divisorGenerator(n):
    large_divisors = []
    for i in range(1, int(np.sqrt(n) + 1)):
        if n % i == 0:
            yield i
            if i*i != n:
                large_divisors.append(n / i)
    for divisor in reversed(large_divisors):
        yield divisor

def rebinCalc(nBins,target):
    rebinFloat = nBins/float(target)
    allDivs = list(divisorGenerator(nBins))
    return find_nearest(allDivs, rebinFloat)

def normHisto(hist, doNorm=False):
    if doNorm:
        if hist.Integral() > 0:
            hist.Scale(1.0/hist.Integral())

def simpleSig(hSig, hBg):
    sig = 0.0
    for i in range(0, hSig.GetNbinsX()):
        totBG = hBg.GetBinContent(i)
        nSig = hSig.GetBinContent(i)
        if(totBG > 1.0 and nSig > 1.0):
            s = nSig / math.sqrt( totBG + (0.3*totBG)**2 )
            sig = math.sqrt(sig**2 + s**2)
    return sig

def getBGHistos(data, histoName, rebinx, xmin, xmax):
    hs = ROOT.THStack()
    hMC = None
    hList = []
    firstPass = True
    for d in data[1]:
        h = d.getHisto(histoName, rebinx=rebinx, xmin=xmin, xmax=xmax, fill=True, showEvents=True)
        hist = copy.deepcopy(h)
        hs.Add(hist)
        hList.append((hist, d.legEntry()))
        if(firstPass):
            hMC = hist
            firstPass = False
        else:
            hMC.Add(hist)
    return hs, hMC, hList

def getData(path, scale=1.0, year = "2018"):
    Data = [
        # info.DataSetInfo(basedir=path, fileName=year+"_DataSR.root",        sys= -1.0, label="Data",        scale=scale),
        info.DataSetInfo(basedir=path, fileName=year+"_DataCR.root",        sys= -1.0, label="Data",        scale=scale),
        # info.DataSetInfo(basedir=path, fileName="2018_Data.root",        sys= -1.0, label="Data",        scale=scale),
        # info.DataSetInfo(basedir=path, fileName="2017_Data.root",        sys= -1.0, label="Data",        scale=scale),
        # info.DataSetInfo(basedir=path, fileName="2016_Data.root",        sys= -1.0, label="Data",        scale=scale),
    ]
    # print("Data = ",Data)

    # qdm_qsmDir = "condor/testHadd_main_01062022_noEtaCut_pT170_withJetCat_pairProduction"
    # Normal
    bgData = [
        # info.DataSetInfo(basedir=path, fileName=year+"_Triboson.root",        label="VVV",                     scale=scale, color=(ROOT.kGray)),
        # info.DataSetInfo(basedir=path, fileName=year+"_Diboson.root",         label="VV",                      scale=scale, color=(ROOT.kMagenta + 1)),
        # info.DataSetInfo(basedir=path, fileName=year+"_DYJetsToLL_M-50.root", label="Z#gamma*+jets",           scale=scale, color=(ROOT.kOrange + 2)),
        # info.DataSetInfo(basedir=path, fileName=year+"_TTX.root",             label="ttX",                     scale=scale, color=(ROOT.kCyan + 1)),
        info.DataSetInfo(basedir=path, fileName=year+"_ST.root",              label="Single top",              scale=scale, color=(ROOT.kRed + 1)),
        # info.DataSetInfo(basedir=path, fileName=year+"_ZJets.root",           label="Z#rightarrow#nu#nu+jets", scale=scale, color=(ROOT.kGray + 1)),
    
        # info.DataSetInfo(basedir=path, fileName=year+"_ST_tZq.root",          label="ST tZq",                scale=scale, color=(ROOT.kRed)),
        # info.DataSetInfo(basedir=path, fileName=year+"_ST_s-channel.root",          label="ST s-channel",                scale=scale, color=(ROOT.kGreen + 3)),
        # info.DataSetInfo(basedir=path, fileName=year+"_ST_tW.root",          label="ST tW",                scale=scale, color=(ROOT.kPink + 7)),
        # info.DataSetInfo(basedir=path, fileName=year+"_ST_t-channel.root",          label="ST t-channel",                scale=scale, color=(ROOT.kTeal)),
        

        info.DataSetInfo(basedir=path, fileName=year+"_TTJets.root",          label="t#bar{t}",                scale=scale, color=(ROOT.kBlue - 6)),
        # info.DataSetInfo(basedir=path, fileName=year+"_WJets.root",           label="W+jets",                  scale=scale, color=(ROOT.kYellow + 1)),
        # info.DataSetInfo(basedir=path, fileName=year+"_QCD.root",             label="QCD",                     scale=scale, color=(ROOT.kGreen + 1)),
        # full bkg sample
        # info.DataSetInfo(basedir=path, fileName=year+"_mTTJetsmini_Inc_noEtaCut_pT50.root",     label="t#bar{t}",                scale=scale, color=(ROOT.kBlue - 6)),
        info.DataSetInfo(basedir=path, fileName=year+"_ZJets.root",             label="Z#rightarrow#nu#nu+jets",    scale=scale, color=(ROOT.kGray + 1)),
        info.DataSetInfo(basedir=path, fileName=year+"_WJets.root",              label="W+jets",                    scale=scale, color=(ROOT.kYellow + 1)),
        # info.DataSetInfo(basedir=path, fileName=year+"_TT.root",              label="t#bar{t} (pow)",             scale=scale, color=(ROOT.kBlue - 6)),
        info.DataSetInfo(basedir=path, fileName=year+"_QCD.root",               label="QCD",                        scale=scale, color=(ROOT.kGreen + 1)),
    ]
    #
    sgData = [

        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",     label="mMed 600",  scale=scale, color=ROOT.kMagenta + 1),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-1_rinv-0p3_alpha-peak_yukawa-1.root",     label="M-2000_mD-1",scale=scale, color=ROOT.kCyan),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-20_rinv-0p1_alpha-peak_yukawa-1.root",    label="M-2000_r-0p1",scale=scale, color=ROOT.kGray+3),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="baseline", scale=scale, color=ROOT.kViolet+2),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-50_rinv-0p3_alpha-peak_yukawa-1.root",    label="M-2000_mD-50",scale=scale, color=ROOT.kViolet+5),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-20_rinv-0p5_alpha-peak_yukawa-1.root",    label="M-2000_r-0p5",scale=scale, color=ROOT.kYellow),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="mMed 4000", scale=scale, color=ROOT.kRed),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-100_rinv-0p3_alpha-peak_yukawa-1.root",   label="M-2000_mD-100",scale=scale, color=ROOT.kPink+9),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-20_rinv-0p7_alpha-peak_yukawa-1.root",    label="M-2000_r-0p7",scale=scale, color=ROOT.kCyan+4),

        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="t-ch 3000", scale=scale, color=ROOT.kMagenta+1),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",     label="t-ch 600",  scale=scale, color=ROOT.kBlack),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-800_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",     label="t-ch 800",  scale=scale, color=ROOT.kGreen),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="t-ch 2000", scale=scale, color=ROOT.kBlue),
        # info.DataSetInfo(basedir=path, fileName="2017_mZprime-3000_mDark-20_rinv-0p3_alpha-peak.root",           label="s-ch baseline", scale=scale, color=ROOT.kRed),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-6000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="t-ch 6000", scale=scale, color=ROOT.kCyan,)
        ## varying mMed

        # # info.DataSetInfo(basedir=path, fileName=year+"_mMed-500_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",     label="mMed 500",  scale=scale, color=ROOT.kMagenta + 1),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",     label="mMed 600",  scale=scale, color=ROOT.kMagenta + 1),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-800_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",     label="mMed 800",  scale=scale, color=ROOT.kRed),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-1000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="mMed 1000", scale=scale, color=ROOT.kBlack),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-1500_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="mMed 1500", scale=scale, color=ROOT.kGray+4),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="baseline", scale=scale, color=ROOT.kOrange+2),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="mMed 3000", scale=scale, color=ROOT.kCyan),
        
        # # info.DataSetInfo(basedir=path, fileName="2017_mZprime-3000_mDark-20_rinv-0p3_alpha-peak.root",          label="s-ch 3000", scale=scale, color=ROOT.kRed),
        # # info.DataSetInfo(basedir=path, fileName=year+"_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1_noEtaCut_pT170.root",    label="mMed 4000", scale=scale, color=ROOT.kRed),
        
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="mMed 4000", scale=scale, color=ROOT.kRed+2),
        
        # ## varying mDark

        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-1_rinv-0p3_alpha-peak_yukawa-1.root",    label="M-2000_mD-1",scale=scale, color=ROOT.kViolet-5),
        
        # # info.DataSetInfo(basedir=path, fileName=year+"_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1_noEtaCut_pT170.root",    label="M-3000_mD-20",scale=scale, color=ROOT.kGreen+2),
        # # info.DataSetInfo(basedir=path, fileName=year+"_mMed-3000_mDark-50_rinv-0p3_alpha-peak_yukawa-1_noEtaCut_pT170.root",    label="M-3000_mD-50",scale=scale, color=ROOT.kRed),
        
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-50_rinv-0p3_alpha-peak_yukawa-1.root",    label="M-2000_mD-50",scale=scale, color=ROOT.kViolet+3),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-100_rinv-0p3_alpha-peak_yukawa-1.root",    label="M-2000_mD-100",scale=scale, color=ROOT.kViolet+6),
        
        # ## varying rinv
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-20_rinv-0p1_alpha-peak_yukawa-1.root",    label="M-2000_r-0p1",scale=scale, color=ROOT.kOrange+3),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-20_rinv-0p5_alpha-peak_yukawa-1.root",    label="M-2000_r-0p5",scale=scale, color=ROOT.kOrange+9),
        
        # # info.DataSetInfo(basedir=path, fileName=year+"_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1_noEtaCut_pT170.root",    label="M-3000_r-0p3",scale=scale, color=ROOT.kGreen+2),
        # # info.DataSetInfo(basedir=path, fileName=year+"_mMed-3000_mDark-20_rinv-0p5_alpha-peak_yukawa-1_noEtaCut_pT170.root",    label="M-3000_r-0p5",scale=scale, color=ROOT.kRed),
        
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-20_rinv-0p7_alpha-peak_yukawa-1.root",    label="M-2000_r-0p7",scale=scale, color=ROOT.kOrange-9),
        
        # ## varying alpha
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-20_rinv-0p3_alpha-low_yukawa-1.root",    label="M-2000_a-low",scale=scale, color=ROOT.kGray),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1_noEtaCut_pT170.root",    label="M-3000_a-peak",scale=scale, color=ROOT.kGreen+2),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-20_rinv-0p3_alpha-high_yukawa-1.root",    label="M-2000_a-high",scale=scale, color=ROOT.kYellow+2),
        ## varying rinv at mMed 800
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-800_mDark-20_rinv-0p1_alpha-peak_yukawa-1.root",     label="M-800_r-0p1", scale=scale, color=ROOT.kOrange + 2),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-800_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",     label="M-800_r-0p3", scale=scale, color=ROOT.kMagenta + 1),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-800_mDark-20_rinv-0p5_alpha-peak_yukawa-1.root",     label="M-800_r-0p5", scale=scale, color=ROOT.kBlack),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-800_mDark-20_rinv-0p7_alpha-peak_yukawa-1.root",     label="M-800_r-0p7", scale=scale, color=ROOT.kGreen),
        ## comparing QdM and QsM jets
        # info.DataSetInfo(basedir="{}/QdM/".format(qdm_qsmDir), fileName=year+"_mMed400.root",     label="mMed 400 QdM",  scale=scale, color=ROOT.kMagenta + 1),
        # info.DataSetInfo(basedir="{}/QsM/".format(qdm_qsmDir), fileName=year+"_mMed400.root",     label="mMed 400 QsM",  scale=scale, color=ROOT.kBlack),
        # info.DataSetInfo(basedir="{}/QdM/".format(qdm_qsmDir), fileName=year+"_mMed600.root",     label="mMed 600 QdM",  scale=scale, color=ROOT.kMagenta + 1),
        # info.DataSetInfo(basedir="{}/QsM/".format(qdm_qsmDir), fileName=year+"_mMed600.root",     label="mMed 600 QsM",  scale=scale, color=ROOT.kBlack),
        # info.DataSetInfo(basedir="{}/QdM/".format(qdm_qsmDir), fileName=year+"_mMed800.root",     label="mMed 800 QdM",  scale=scale, color=ROOT.kMagenta + 1),
        # info.DataSetInfo(basedir="{}/QsM/".format(qdm_qsmDir), fileName=year+"_mMed800.root",     label="mMed 800 QsM",  scale=scale, color=ROOT.kBlack),
        # info.DataSetInfo(basedir="{}/QdM/".format(qdm_qsmDir), fileName=year+"_mMed1000.root",    label="mMed 1000 QdM", scale=scale, color=ROOT.kMagenta + 1),
        # info.DataSetInfo(basedir="{}/QsM/".format(qdm_qsmDir), fileName=year+"_mMed1000.root",    label="mMed 1000 QsM", scale=scale, color=ROOT.kBlack),

    ]
    # print(sgData, bgData)
    return Data, sgData, bgData

def setupAxes(dummy, xOffset, yOffset, xTitle, yTitle, xLabel, yLabel):
    dummy.SetStats(0)
    dummy.SetTitle("")
    dummy.GetXaxis().SetTitleOffset(xOffset)
    dummy.GetYaxis().SetTitleOffset(yOffset)
    dummy.GetXaxis().SetTitleSize(xTitle)
    dummy.GetYaxis().SetTitleSize(yTitle)
    dummy.GetXaxis().SetLabelSize(xLabel)
    dummy.GetYaxis().SetLabelSize(yLabel)
    if(dummy.GetXaxis().GetNdivisions() % 100 > 5): dummy.GetXaxis().SetNdivisions(6, 5, 0)

def setupDummy(dummy, leg, histName, xAxisLabel, yAxisLabel, isLogY, xmin, xmax, ymin, ymax, lmax, norm=False, normBkg=False,isRatio=False):
    
    if isRatio:
        setupAxes(dummy, 0, 1.05, 0.0, 0.05, 0.0, 0.05)
        dummy.GetXaxis().SetTitle("")
    
    else:    
        setupAxes(dummy, 1.2, 1.6, 0.045, 0.045, 0.045, 0.045)
        dummy.GetXaxis().SetTitle(xAxisLabel)
    
    dummy.GetYaxis().SetTitle(yAxisLabel)
    dummy.SetTitle(histName)
    #Set the y-range of the histogram
    if(isLogY):
        if norm:
            default = 0.00001
        else:
            default = 0.02
        locMin = min(default, max(default, 0.05 * ymin))
        legSpan = (math.log10(3*ymax) - math.log10(locMin)) * (leg.GetY1() - ROOT.gPad.GetBottomMargin()) / ((1 - ROOT.gPad.GetTopMargin()) - ROOT.gPad.GetBottomMargin())
        legMin = legSpan + math.log10(locMin)
        if(math.log10(lmax) > legMin):
            scale = (math.log10(lmax) - math.log10(locMin)) / (legMin - math.log10(locMin))
            if norm:
                ymax = 2.
            else:
                ymax = pow(ymax/locMin, scale)*locMin
                # ymax = 10**8
        # dummy.GetYaxis().SetRangeUser(locMin, 10*ymax)
        dummy.GetYaxis().SetRangeUser(locMin, ymax*10e-5)
    else:
        locMin = 0.0
        legMin = (1.2*ymax - locMin) * (leg.GetY1() - ROOT.gPad.GetBottomMargin()) / ((1 - ROOT.gPad.GetTopMargin()) - ROOT.gPad.GetBottomMargin())
        if(lmax > legMin): ymax *= (lmax - locMin)/(legMin - locMin)
        # dummy.GetYaxis().SetRangeUser(0.0, ymax*1.2)
        dummy.GetYaxis().SetRangeUser(0.0, ymax)
    #set x-axis range
    if(xmin < xmax): dummy.GetXaxis().SetRangeUser(xmin, xmax)
    # print("ymax in dummy = ",ymax)

def makeRocVec(h,reverse=False,ignoreUnderflow=False):
    if h.Integral() > 0.0:
        h.Scale( 1.0 / h.Integral() );
    v, cuts = [], []
    si = 1
    if ignoreUnderflow == True:
        si = -1
    for i in range(si, h.GetNbinsX()+1):
        if reverse:
            val = h.Integral(si, i)
        else:
            val = h.Integral(i, h.GetNbinsX())
        v.append(val)
        cuts.append(h.GetBinLowEdge(i))
    return v, cuts

def ROCArea(n,mBg,mSig):
    mBgAr = [1] + mBg + [0]
    mSigAr = [0] + mSig + [0]
    gAr = ROOT.TGraph(n, array("d", mBgAr), array("d", mSigAr))
    gArea = round(gAr.Integral(),2) # original way
    # gArea = round(np.trapz(mBgAr,mSigAr),2)
    return gArea

def drawRocCurve(fType, rocBgVec, rocSigVec, leg, manySigs=False, stList=None, allRocValues=None, baseline="baseline", mainBkg = "QCD"):
    # saving all the ROC scores for all signals and backgrounds
    rocValues = pd.DataFrame(columns=["cut","var","sig","bkg","roc_auc","cutDir","varCut","cBg","cSig","mBg","mSig"])
    for mBg, cutBg, lBg, cBg in rocBgVec:
        for mSig, cutSig, lSig, cSig in rocSigVec:
            n = len(mBg)
            gArea = ROCArea(n,mBg,mSig)
            rv = ">=cut"
            if gArea < 0.5:
                mBg_f = 1 - np.array(mBg)
                mSig_f = 1 - np.array(mSig)
                rv = "<=cut"
                gArea = 1 - gArea
            else:
                mBg_f = mBg
                mSig_f = mSig
            rocValues.loc[len(rocValues.index)] = stList + [lSig,lBg,round(gArea,3),rv,cutSig,cBg,cSig,mBg_f,mSig_f]
            allRocValues.loc[len(allRocValues.index)] = stList + [lSig,lBg,round(gArea,3),rv,cutSig,colorDict[cBg],colorDict[cSig],mBg_f,mSig_f]
    if manySigs:
        rocValues = rocValues[rocValues["bkg"] == mainBkg]
        colLabel = "cSig"
        varMCLabel = "sig"
        mainMC = mainBkg
    else:
        rocValues = rocValues[rocValues["sig"] == baseline]
        colLabel = "cBg"
        varMCLabel = "bkg"
        mainMC = baseline
    h = []
    for varMC in list(rocValues[varMCLabel]):
        datai = rocValues[rocValues[varMCLabel] == varMC].iloc[0]
        n = len(datai["mBg"])
        g = ROOT.TGraph(n, array("d", datai["mBg"]), array("d", datai["mSig"]))
        rebinx = rebinCalc(n,20)
        for i in range(0,n):
            if i % rebinx == 0:
                latex = ROOT.TLatex(g.GetX()[i], g.GetY()[i],str(round(datai["varCut"][i],2)))
                latex.SetTextSize(0.02)
                latex.SetTextColor(ROOT.kRed)
                g.GetListOfFunctions().Add(latex) # add cut values
        g.SetLineWidth(2)
        g.SetLineColor(datai[colLabel])
        g.SetMarkerSize(0.7)
        g.SetMarkerStyle(ROOT.kFullSquare)
        g.SetMarkerColor(datai[colLabel])
        g.Draw("same LP text")
        leg.AddEntry(g, "#splitline{" + fType + " " + mainMC + " vs " + varMC + "_" + datai["cutDir"] + "}{("+"{:.2f}".format(datai["roc_auc"])+")}", "LP")
        h.append(g)
    return h

def plotSignificance(data, histName, totalBin, xlab, plotOutDir, cut, isLogY=False, rebinx=-1.0, xmin=999.9, xmax=-999.9, reverseCut=False, signifValues=None):
    rocBgVec = []
    histoName = histName + cut
    outputPath = plotOutDir+"/FOM/"+cut[1:]
    rebinValue = rebinx # how many bins to merge into 1 bin

    # background
    print("histoName",histoName)
    for d in data[1]:
        h = d.getHisto(histoName, rebinx=-1, xmin=xmin, xmax=xmax, fill=True, showEvents=False, overflow=True)
        h.Rebin(rebinValue)
        hIn = h.Integral()
        effList = np.array(makeRocVec(h,reverseCut)[0]) * hIn
        rocBgVec.append([effList])

    sigLabelList = []
    # signal
    rocSigVec = []
    for d in data[2]:
        h = d.getHisto(histoName, rebinx=-1, xmin=xmin, xmax=xmax, fill=True, showEvents=False, overflow=True)
        h.Rebin(rebinValue)
        hIn = h.Integral()
        eff = np.array(makeRocVec(h,reverseCut)[0])
        effList = eff * hIn
        rocSigVec.append([effList,d.legEntry()])
        sigLabelList.append(d.legEntry())
        cutValues = np.array(makeRocVec(h)[1])

    B = np.zeros(len(cutValues))
    for rbv in rocBgVec:
        B += rbv[0]

    fomList = []
    normedfomList = []
    normedcutList = []
    cutList = []

    colorList = ['b','g','r','c','m','y']
    lineStyles = ["solid","dashed","dotted"]
    # comparing the locations of maximum FOM for different signals
    fig = plt.figure(figsize=(12,8))
    ax = plt.subplot(111)
    for i in range(len(rocSigVec)):
        rsv = rocSigVec[i]
        fo = fom(rsv[0],B)
        foReal = np.ma.masked_invalid(fo)
        fomList.append(foReal)
        maxFOM = foReal.max()
        normedfomList.append(foReal/maxFOM)
        cutList.append(cutValues)
        newEntry = [cut,histName.replace("h_",""),rsv[1],maxFOM]
        signifValues.loc[len(signifValues.index)] = newEntry
        lstyle = "solid"
        mstyle = "o"
        if i > len(colorList) - 1:
            lstyle = "dashed"
            mstyle = "D"
        if i > len(colorList)*2 - 1:
            lstyle = "^"
        ax.plot(cutValues[:-1],foReal[1:]/maxFOM,label=rsv[1] + " ({:.1e})".format(maxFOM), marker=mstyle, linestyle=lstyle, color=colorList[i%len(colorList)])

    pltTitle = ">= cut"
    if reverseCut:
        pltTitle = "<= cut"

    ax.plot(cutList[0][:-1],np.ma.masked_invalid(normedfomList).mean(axis=0)[1:],label="Average",linewidth=5,alpha=0.5,color="black")
    ax.set_title(pltTitle)
    ax.legend(loc='upper right', fontsize=12, ncol=3)
    ax.set_ylabel("Normalized FOM ( sqrt(2((S+B)*log(1+S/B)-S)) )")
    ax.set_xlabel(xlab)
    ax.set_ylim(0,1.5)
    plt.savefig(outputPath + "/FOM_" + histoName+".png")

    # # comparing signal efficiency for different signals
    # plt.figure(figsize=(12,8))
    # for i in range(len(sEffList)):
    #     seff = sEffList[i]
    #     siglab = sigLabelList[i]
    #     plt.step(cutList[i],seff,label=siglab)
    # plt.legend()
    # plt.ylabel("Signal Efficiency")
    # plt.xlabel(xlab)
    # plt.grid()
    # plt.savefig(outputPath + "/sigEff_" + histoName+".png")

    # # comparing FOM and signal efficiency for each signal
    # for i in range(len(sigLabelList)):
    #     slabel = sigLabelList[i]
    #     foms = fomList[i]
    #     cutValues = cutList[i]
    #     seffs = sEffList[i]
    #
    #     fig, ax1 = plt.subplots()
    #
    #     color = 'red'
    #     ax1.set_xlabel(xlab)
    #     ax1.set_ylabel("FOM: sqrt(2((S+B)*log(1+S/B)-S))")
    #     ax1.step(cutValues,foms,color=color)
    #     ax1.tick_params(axis='y',labelcolor=color)
    #
    #     ax2 = ax1.twinx()
    #
    #     color = 'blue'
    #     ax2.set_ylabel('Signal Efficiency', color=color)
    #     ax2.step(cutValues,seffs,color=color)
    #     ax2.tick_params(axis='y',labelcolor=color)
    #
    #     fig.tight_layout()
    #     # for some reason, the x grid line is just not showing right, that's why we have the following code
    #     for cut in np.arange(min(cutValues),max(cutValues),50):
    #         plt.vlines(cut,0,max(seffs),color="silver",linewidth=0.5)
    #     plt.grid()
    #     plt.savefig(outputPath + "/FOMSEff_" + histoName + "_" +  slabel + ".png")

    plt.close()




def plotROC(data, histoName, outputPath="./", isLogY=False, xmin=999.9, xmax=-999.9, norm=False, manySigs=False, stList=None, allRocValues=None):
    #This is a magic incantation to disassociate opened histograms from their files so the files can be closed
    ROOT.TH1.AddDirectory(False)

    #create the canvas for the plot
    c1 = ROOT.TCanvas( "c", "c", 800, 800)

    c1.cd()

    ROOT.gPad.Clear()
    ROOT.gStyle.SetOptStat("")
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.08)
    ROOT.gPad.SetBottomMargin(0.12)
    ROOT.gPad.SetTicks(1,1)
    ROOT.gPad.SetLogy(isLogY)

    #Create TLegend
    leg = ROOT.TLegend(0.17, 0.72, 0.95, 0.88)
    #nColumns = 3 if(len(data[1]) >= 3) else 1
    nColumns = 2
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetLineWidth(1)
    leg.SetNColumns(nColumns)
    leg.SetTextFont(42)
    ROOT.gStyle.SetLegendTextSize(0.024)

    rocBgVec = []
    for d in data[1]:
        h = d.getHisto(histoName, rebinx=-1, xmin=xmin, xmax=xmax, fill=True, showEvents=False, overflow=False)
        rocBgVec.append(makeRocVec(h) + ( d.legEntry(), d.getColor()))

    rocSigVec = []
    for d in data[2]:
        h = d.getHisto(histoName, rebinx=-1, xmin=xmin, xmax=xmax, fill=True, showEvents=False, overflow=False)
        rocSigVec.append(makeRocVec(h) + (d.legEntry(), d.getColor()))

    #create a dummy histogram to act as the axes
    ymax=1.0
    ymin=10**-4
    lmax=1.0
    dummy = ROOT.TH1D("dummy", "dummy", 1000, 0.0, 1.0)
    setupDummy(dummy, leg, "", "#epsilon_{ bg}", "#epsilon_{ sg}", isLogY, xmin, xmax, ymin, ymax, lmax)
    # print("ymax = ",ymax)
    dummy.Draw("hist")
    leg.Draw("same")
    print(histoName)
    history = drawRocCurve("", rocBgVec, rocSigVec, leg, manySigs, stList, allRocValues)

    line1 = ROOT.TF1( "line1","1",0,1)
    line1.SetLineColor(ROOT.kBlack)
    line1.Draw("same")
    line2 = ROOT.TF1( "line2","x",0,1)
    line2.SetLineColor(ROOT.kBlack)
    line2.SetLineStyle(ROOT.kDotted)
    line2.Draw("same")

    dummy.Draw("AXIS same")
    # dummy.GetXaxis().SetRangeUser(0,0.1)
    # dummy.SetMinimum(0.8)
    # dummy.SetMaximum(1.1)
    # CMS label
    CMS_lumi.writeExtraText = 1
    lumi = "59.7"

    CMS_lumi.lumi_sqrtS = lumi + " fb^{-1} (13 TeV)"

    iPeriod = 0
    iPos = 0

    CMS_lumi.CMS_lumi(c1, iPeriod, iPos)
    c1.cd()
    c1.Update();
    c1.RedrawAxis()

    c1.SaveAs(outputPath+"/"+histoName+"_ROC.png")
    c1.Close()
    del c1
    del leg

def createRatio(h1, h2, xtitle):
    h3 = h1.Clone("h3"+xtitle)
    h3.SetLineColor(ROOT.kBlack)
    h3.SetMarkerStyle(20)
    h3.SetMarkerSize(1)
    h3.SetMarkerColor(ROOT.kBlack)
    h3.SetTitle("")
    
    h3.SetMinimum(0)
    h3.SetMaximum(2)
	# Set up plot for markers and errors
	#h3.Sumw2()
    h3.SetStats(0)
    h3.Divide(h2)
    # ymax = h3.GetMaximumBin()
    # print("ymax in create Ratio = ",ymax)
    # h3.SetMaximum(ymax)
	# Adjust y-axis settings
    x = h3.GetXaxis()
    y = h3.GetYaxis()
    
    x.SetTitleOffset(0.65)
    y.SetTitleOffset(0.27)
    x.SetTitleSize(0.2)
    y.SetTitleSize(0.15)
    x.SetLabelSize(0.13)
    y.SetLabelSize(0.13)
    x.SetTitle(xtitle)
    # print("xtitle",xtitle)
    y.SetTitle("Data/MC")

    if(x.GetNdivisions() % 100 > 5): x.SetNdivisions(6, 5, 0)

    y.SetNdivisions(505)

    # y.SetTitleFont(43)
    
    # y.SetLabelFont(43)
 
    # x.SetTitleFont(43)
    # x.SetLabelFont(43)
    # x.SetLabelOffset(0.05)

    return h3

def createCanvasPads(c,isLogY):
	# Upper histogram plot is pad1
    # eps = 0.005
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1.0, 0.97)
    pad1.SetBottomMargin(0.01)  # joins upper and lower plot
    pad1.SetLeftMargin(0.10)
    pad1.SetRightMargin(0.05)
    pad1.SetTopMargin(0.1)
    # pad1.SetBottomMargin(0.12)
    pad1.SetTicks(1,1)
    pad1.SetLogy(isLogY)
    pad1.Draw()
    # Lower ratio plot is pad2
    c.cd()  # returns to main canvas before defining pad2
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.0, 1.0, 0.3)
    pad2.SetTopMargin(0)  # joins upper and lower plot
    # pad2.SetBottomMargin(0.3)
    pad2.SetLeftMargin(0.10)
    pad2.SetRightMargin(0.05)
    # pad2.SetTopMargin(0.08)
    pad2.SetBottomMargin(0.35)
    pad2.SetTicks(1,1)
    pad2.SetGrid()
    pad2.Draw()

    return c, pad1, pad2


def plotStack(data, histoName, totalBin, outputPath="./", xTitle="", yTitle="", isLogY=False, xmin=999.9, xmax=-999.9, norm=False, normBkg=False, onlySig=False, stList=None, yieldValues=None, isRatio=False):
    #This is a magic incantation to disassociate opened histograms from their files so the files can be closed
    ROOT.TH1.AddDirectory(False)
    # print("Data in plot stack = ",data)
    #create the canvas for the plot
    
    
    if isRatio:
        c1 = ROOT.TCanvas( "c", "c", 800, 700)
        c1, pad1, pad2 = createCanvasPads(c1,isLogY)
        pad1.cd()
    else:
        c1 = ROOT.TCanvas( "c", "c", 800, 800)
        c1.cd()
        ROOT.gPad.Clear()
        ROOT.gPad.SetLeftMargin(0.15)
        ROOT.gPad.SetRightMargin(0.05)
        ROOT.gPad.SetTopMargin(0.08)
        ROOT.gPad.SetBottomMargin(0.12)
        ROOT.gPad.SetTicks(1,1)
        ROOT.gPad.SetLogy(isLogY)
    
    ROOT.gStyle.SetOptStat("")

    #Create TLegend
    leg = ROOT.TLegend(0.17, 0.7, 0.95, 0.88)
    #nColumns = 3 if(len(data[1]) >= 3) else 1
    nColumns = 2
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetLineWidth(1)
    leg.SetNColumns(nColumns)
    leg.SetTextFont(42)
    if isRatio:
        ROOT.gStyle.SetLegendTextSize(0.05)
    else:    
        ROOT.gStyle.SetLegendTextSize(0.024)
            
    #Setup background histos
    rebinx = rebinCalc(totalBin,40)
    # print("rebinx = ",rebinx)
    hs = ROOT.THStack()
    hMC = None
    firstPass = True
    bkghist = None
    history = []
    
    # setup background histos
    # print("before setup bkg histos detail key = {}, {}, {}".format(histoName,xmin, xmax))
    for d in data[1]:
        h = d.getHisto(histoName, rebinx=rebinx, xmin=xmin, xmax=xmax, fill=True, showEvents=True)
        if (stList != None) and (not normBkg):
            newEntry = stList + [getLabel(d.legEntry()),round(h.Integral())]
            yieldValues.loc[len(yieldValues.index)] = newEntry
            # print(newEntry)
        if normBkg:
            normHisto(h, True)
            h.SetLineWidth(3)
            h.SetFillStyle(3955) 
        hs.Add(copy.deepcopy(h))
        leg.AddEntry(h, d.legEntry(), "F")
        # bkghist += h
        if(firstPass):
            hMC = copy.deepcopy(h)
            bkghist = h
            firstPass = False
        else:
            hMC.Add(copy.deepcopy(h))
            bkghist.Add(h)

    

    print("hs = ", hs)
    print("hMC = ",hMC)
    # there is a bug with getBGHistos. Once fixed, can delete lines 294-305, and uncomment
    # the line below and lines 313-314
    # hs, hMC, hList = getBGHistos(data, histoName, rebinx, xmin, xmax)
    if norm:
        normHisto(hMC, True)
    #Fill background legend
    # for h in hList:
    #     leg.AddEntry(h[0], h[1], "F")

    #create a dummy histogram to act as the axes
    if norm:
        ymax=10**1
        ymin=10**-12
        lmax=10**1
    else:
        ymax=10**11
        ymin=10**-4
        lmax=10**12
    dummy = ROOT.TH1D("dummy", "dummy", 1000, hMC.GetBinLowEdge(1), hMC.GetBinLowEdge(hMC.GetNbinsX()) + hMC.GetBinWidth(hMC.GetNbinsX()))
    setupDummy(dummy, leg, "", xTitle, yTitle, isLogY, xmin, xmax, ymin, ymax, lmax, norm, normBkg,isRatio)
    # print("near setup dummy detail key = {}, {}, {}".format(histoName,xmin, xmax))
        
    # setupDummy(dummy, leg, histoName, xTitle, yTitle, isLogY, xmin, xmax, ymin, ymax, lmax, norm, normBkg)
    if normBkg:
        dummy.SetMaximum(100)
        dummy.SetMinimum(0.00001)
    dummy.Draw("hist")
    if norm:
        hMC.Draw("hist same")
        leg.Clear()
        leg.AddEntry(hMC, "Total Background", "L")
    elif normBkg:
        hs.Draw("nostackHIST same")
        hs.SetMaximum(100)
        hs.SetMinimum(0.00001)
    elif onlySig:
        leg.Clear()
    else:
        hs.Draw("hist F same")
    leg.Draw("same")

    # print("detail key = {}, {}, {}".format(histoName,xmin, xmax))
    
    #Setup signal histos
    
    sig = 0.0
    linestylenumber = 0
    linestyle = [ROOT.kSolid,ROOT.kDashed,ROOT.kDotted]
    if(data[2]):
        #firstPass=True
        for d in data[2]:
            h = d.getHisto(histoName, rebinx=rebinx, xmin=xmin, xmax=xmax, showEvents=True)
            if (stList != None) and (not normBkg):
                newEntry = stList + [getLabel(d.legEntry()),round(h.Integral())]
                yieldValues.loc[len(yieldValues.index)] = newEntry
                # print("Signal = ",newEntry)
            #if(firstPass):
            sig = round(simpleSig(h, hMC),2)
            #firstPass=False
            #print(d.legEntry(), round(simpleSig(h, hMC),2))
            # h.SetLineStyle(ROOT.kDashed)
            h.SetLineStyle(linestyle[linestylenumber%3] )
            linestylenumber+=1
            h.SetLineWidth(3)
            leg.AddEntry(h, d.legEntry()+", {}".format(sig), "L")
            if norm or normBkg:
                normHisto(h, True)
            h.Draw("hist same")
            history.append(h)
        
    # Setup data histogram
    if(data[0]):
        for d in data[0]:
            datahist = d.getHisto(histoName, rebinx=rebinx,xmin=xmin,xmax=xmax,showEvents=True)
            if (stList != None) and (not normBkg):
                newEntry = stList + [getLabel(d.legEntry()),round(datahist.Integral())] # For the yield value calculation
                yieldValues.loc[len(yieldValues.index)] = newEntry
                # print("Data = ",newEntry)
                #firstPass=False
            ROOT.gStyle.SetErrorX(0.)
            datahist.SetMarkerStyle(20)
            datahist.SetMarkerSize(1)
            datahist.SetLineColor(ROOT.kBlack)
            leg.AddEntry(datahist, d.legEntry(), "P")
            # if norm or normBkg:
            #     normHisto(datahist, True)
            datahist.Draw("P same")
            dhist = datahist
            history.append(datahist)
            # Print ratio plot
            if isRatio:
                pad2.cd()
                # print("bkghist = ",bkghist)
                # setupAxes(dummy, 1.2, 1.6, 0.045, 0.045, 0.045, 0.045)
                ratio = createRatio(datahist,bkghist,xTitle)
                ratio.Draw("EX0P")



    

    #Draw significance
    significance = ROOT.TLatex()
    significance.SetNDC(True)
    significance.SetTextAlign(11)
    significance.SetTextFont(52)
    significance.SetTextSize(0.030)
    #significance.DrawLatex(0.45, 0.72, ("Significance = #frac{N_{s}}{#sqrt{N_{b}+#left(0.3N_{b}#right)^{2}}} = "+str(sig)))
    #significance.DrawLatex(0.45, 0.72, ("Significance = #frac{N_{s}}{#sqrt{N_{b}+#left(0.3N_{b}#right)^{2}}}"))

    if onlySig:
        dummy.SetMaximum(10**8)
    dummy.Draw("AXIS same")

    # ran = [870, 2385.0]
    # vl = ROOT.TLine(ran[0],0,ran[0],10**6)
    # vl.SetLineWidth(2)
    # vl.SetLineColor(ROOT.kRed)
    # vl.Draw("same")
    # vl2 = ROOT.TLine(ran[1],0,ran[1],10**6)
    # vl2.SetLineWidth(2)
    # vl2.SetLineColor(ROOT.kRed)
    # vl2.Draw("same")
    # CMS label
    CMS_lumi.writeExtraText = 1
    lumi = "59.7"

    CMS_lumi.lumi_sqrtS = lumi + " fb^{-1} (13 TeV)"

    iPeriod = 0
    iPos = 0

    CMS_lumi.CMS_lumi(c1, iPeriod, iPos)
    c1.cd()
    c1.Update();
    c1.RedrawAxis()

    if norm:
        c1.SaveAs(outputPath+"/"+histoName+"_norm.png")
    elif normBkg:
        c1.SaveAs(outputPath+"/"+histoName+"_normBkg.png")
    else:
        c1.SaveAs(outputPath+"/"+histoName+".png")

    c1.Close()
    del c1
    del leg
    del hMC

def main():
    parser = optparse.OptionParser("usage: %prog [options]\n")
    parser.add_option('-b',                 dest='isNormBkg',  action="store_true",                            help="Normalized Background and Signal plots")
    parser.add_option('-d', '--dataset',    dest='dataset',                    default='testHadd_11242020',    help='dataset')
    parser.add_option('-j', '--jNVar',      help='make histograms for nth jet variables', dest='jNVar', default=False, action='store_true')
    parser.add_option('-m',                 dest='manySigs',   action="store_true",                            help="Plot ROC curves with many signals vs. QCD")
    parser.add_option('-n',                 dest='isNorm',     action="store_true",                            help="Normalize stack plots")
    parser.add_option('-s',                 dest='onlySig',    action="store_true",                            help="Plot only signals")
    parser.add_option('-y',                 dest='year',       type='string',  default='2018',                 help="Can pass in the run year")
    parser.add_option('-o',                 dest='outputdir',  type='string',                                  help="Output folder name")
    # parser.add_optio
    options, args = parser.parse_args()

    year = options.year
    # cuts = ["", "_ge2AK8j", "_ge2AK8j_lp6METrST", "_ge2AK8j_l1p5dEta12", "_baseline"]
    #cuts = ["_ge2AK8j"]
    # cutsImportant = ["_qual_trg_st","_qual_trg_st_0nim","_qual_trg_st_ge1nim"]
    # cutsImportant = ["","_2PJ","_2PJ_nl","_qual_trg_2PJ", "_qual_trg_st_2PJ"]
    # cutsImportant = ["","_2PJ","_2PJ_nl","_qual_trg_2PJ","_qual_trg_st_2PJ","_qual_trg_st_ht_2PJ_dphimin"]
    # cutsImportant = ["_cr"]
    # cutsImportant = ["issue_ht","_issue_met"]
    cutsImportant = ["_cr_muon","_cr_electron"]
    


    Data, sgData, bgData = getData("condor/" + options.dataset + "/", 1.0, year)
   
    # print(sgData)
    #Data, sgData, bgData = getData("condor/MakeNJetsDists_"+year+"/", 1.0, year)
    allRocValues = pd.DataFrame(columns=["cut","var","sig","bkg","roc_auc","cutDir","cutSig","cBg","cSig","mBg_f","mSig_f"])
    yieldValues = pd.DataFrame(columns=["cut","var","source","yield"])
    signifValues = pd.DataFrame(columns=["cut","var","source","max signif."])
    if options.outputdir:
        plotOutDir = "output/{}".format(options.outputdir)
    else: 
        plotOutDir = "output/{}".format(options.dataset)

    preVars = {
        "h_njets":False,
        "h_njetsAK8":False,
        "h_nb":False,
        "h_ht":False,
        "h_st":False,
        "h_met":False,
        # "h_mT":True,
        # "h_METrHT_pt30":False,
        # "h_METrST_pt30":False,
        # "h_dEtaj12AK8":True,
        # "h_dRJ12AK8":True,
        # "h_dPhij1METAK8":False,
        # "h_dPhij2METAK8":False,
        # "h_dPhij1rdPhij2AK8":False,
        # "h_dPhiMinjMETAK8":False,
        # "h_dEtaj12":True,
        # "h_dRJ12":True,
        # "h_dPhij1MET":False,
        # "h_dPhij2MET":False,
        # "h_dPhij1rdPhij2":True,
        # "h_dPhiMinjMET":True,
        # "h_mT2_f4_msm":False,
        # "h_mT2_f4_msm_dEta":False,
        # "h_mT2_f4_msm_dPhi":False,
        # "h_mT2_f4_msm_dR":False,
    }
    varsSkip = [
    "eCounter",
    "evtw",
    "jw",
    "fjw"
    ]

# myvars = key : ["xlabel", no. of bins, xmin,xmax, npzinfo, flattenInfo, weightName]
    
    for histName,details in vars(options.jNVar).items():
    # for histName, details in myVars.items():
        # print(histName)
        # print(details)
        isNorm = options.isNorm
        isNormBkg = options.isNormBkg
        onlySig = options.onlySig
        manySigs = options.manySigs
        if histName in varsSkip:
            continue
        #if details[6] != "evtw":
        #    continue
        for cut in cutsImportant:
            makeDirs(plotOutDir,cut,"Stacked")
            makeDirs(plotOutDir,cut,"roc")
            makeDirs(plotOutDir,cut,"FOM")
            makeDirs(plotOutDir,cut,"NormedStacked")
            stList = [cut,histName]
            # print("Data = ",Data)
            # plotROC(  (Data, bgData, sgData), "h_"+histName+cut, plotOutDir+"/roc/"+cut[1:], isLogY=False,   manySigs=manySigs, stList=stList, allRocValues=allRocValues)
            # plotStack((Data, bgData, sgData), "h_"+histName+cut, details[1], plotOutDir+"/Stacked/"+cut[1:], details[0], "Events", isLogY=True, norm=isNorm, xmin=details[2], xmax=details[3], normBkg=False, onlySig=onlySig, stList=stList, yieldValues=yieldValues,isRatio=False)
            plotStack((Data, bgData, sgData), "h_"+histName+cut, details[1], plotOutDir+"/Stacked/"+cut[1:], details[0], "Events", isLogY=True,norm=isNorm, xmin=details[2], xmax=details[3], normBkg=False, onlySig=onlySig, stList=stList, yieldValues=yieldValues,isRatio=True)
            # plotStack((Data, bgData, sgData), "h_"+histName+cut, details[1], plotOutDir+"/NormedStacked/"+cut[1:], details[0], "Events", isLogY=True, norm=isNorm, xmin=details[2], xmax=details[3], normBkg=True, onlySig=onlySig, stList=stList, yieldValues=yieldValues)
            if histName in preVars.keys():
                plotSignificance((Data, bgData, sgData), "h_"+histName, details[1], details[0], plotOutDir, cut,                    isLogY=False, reverseCut=preVars[histName], signifValues=signifValues)
    yieldValues.to_csv("{}/yieldValues.csv".format(plotOutDir))
    signifValues.to_csv("{}/signifValues.csv".format(plotOutDir))
    allRocValues.to_csv("{}/allRocValues.csv".format(plotOutDir))
    # print(yieldValues)

if __name__ == '__main__':
    main()
