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
import itertools
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
        info.DataSetInfo(basedir=path, fileName=f"{year}_Data.root",        sys= -1.0, label="Data",        scale=scale),
        # info.DataSetInfo(basedir=path, fileName="2018_Data.root",        sys= -1.0, label="Data",        scale=scale),
        # info.DataSetInfo(basedir=path, fileName="2018_Data_notSkims.root",        sys= -1.0, label="Data",        scale=scale),
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

        # # info.DataSetInfo(basedir=path, fileName=year+"_mMed-500_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",     label="mMed_500",  scale=scale, color=ROOT.kMagenta + 1),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1.root",    label="baseline", scale=scale, color=ROOT.kBlack),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-600.root",     label="mMed_600",  scale=scale, color=ROOT.kViolet+2),
        # # info.DataSetInfo(basedir=path, fileName=year+"_mMed-800.root",     label="mMed_800",  scale=scale, color=ROOT.kRed),
        # # info.DataSetInfo(basedir=path, fileName=year+"_mMed-1000.root",    label="mMed_1000", scale=scale, color=ROOT.kMagenta + 1),
        # # info.DataSetInfo(basedir=path, fileName=year+"_mMed-1500.root",    label="mMed_1500", scale=scale, color=ROOT.kGray+4),
        # # info.DataSetInfo(basedir=path, fileName=year+"_mMed-3000.root",    label="mMed_3000", scale=scale, color=ROOT.kCyan),
        # info.DataSetInfo(basedir=path, fileName=year+"_mMed-4000.root",    label="mMed_4000", scale=scale, color=ROOT.kBlue + 1),
        
        # # info.DataSetInfo(basedir=path, fileName="2017_mZprime-3000_mDark-20_rinv-0p3_alpha-peak.root",          label="s-ch 3000", scale=scale, color=ROOT.kRed),
        # # info.DataSetInfo(basedir=path, fileName=year+"_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1_noEtaCut_pT170.root",    label="mMed 4000", scale=scale, color=ROOT.kRed),
        
        
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

def setupAxes(dummy, xOffset, yOffset, xTitle, yTitle, xLabel, yLabel,title=""):
    dummy.SetStats(0)
    dummy.SetTitle(title)
    dummy.GetXaxis().SetTitleOffset(xOffset)
    dummy.GetYaxis().SetTitleOffset(yOffset)
    dummy.GetXaxis().SetTitleSize(xTitle)
    dummy.GetYaxis().SetTitleSize(yTitle)
    dummy.GetXaxis().SetLabelSize(xLabel)
    dummy.GetYaxis().SetLabelSize(yLabel)
    if(dummy.GetXaxis().GetNdivisions() % 100 > 5): dummy.GetXaxis().SetNdivisions(6, 5, 0)

def setupDummy(dummy, leg, histName, xAxisLabel, yAxisLabel, isLogY, xmin, xmax, ymin, ymax, lmax, title="", norm=False, normBkg=False,isRatio=False,isABCD=False):
    
    if isRatio:
        setupAxes(dummy, 0, 1.05, 0.0, 0.05, 0.0, 0.05)
        dummy.GetXaxis().SetTitle("")
    elif isABCD:
        setupAxes(dummy, 0.6, 0.8, 0.1, 0.1, 0.08, 0.06,title)
        dummy.GetXaxis().SetTitle(xAxisLabel)
        dummy.GetXaxis().CenterLabels()
        
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
    if not isABCD:
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
    print(f"xTitle = {xtitle}")
    print(f"Type of h2: {type(h2)}")
    if h1.GetEntries() == 0 or h2.GetEntries() == 0:
        print("One of the histograms is empty.")
    h3 = h1.Clone("h3")
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


def plotStack(data, histoName, totalBin, outputPath="./", xTitle="", yTitle="", isLogY=False, xmin=999.9, xmax=-999.9, norm=False, normBkg=False, onlySig=False, stList=None, yieldValues=None, year="2017", isRatio=False, hemPeriod=False):
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
    # bkghist = None
    bkghist = []
    history = []
    labels= []
    bkgcontrib = []
    
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
            # h.SetLineColor(d.color) 
        elif norm:
            normHisto(h,True)
            h.SetLineColor(d.Color)
            # h.Draw("hist same")
            leg.AddEntry(h, d.legEntry(), "L")
        else:    
            h.SetFillStyle(3001)
        hs.Add(copy.deepcopy(h))
        bkgcontrib.append(round(h.Integral()))
        labels.append(d.legEntry())
        
        # leg.AddEntry(h, d.legEntry(), "F")
        # bkghist += h
        bkghist.append(copy.deepcopy(h))
        if(firstPass):
            hMC = copy.deepcopy(h)
            # bkghist = h
            firstPass = False
        else:
            hMC.Add(copy.deepcopy(h))
            # bkghist.Add(h)

    # Print legends in the reverse order along with the bkg contribution 
    totalOfAllhisto = sum(bkgcontrib)
    percentOfAllhisto = [(x/totalOfAllhisto)*100 if totalOfAllhisto != 0 else x for x in bkgcontrib  ]
    print(f"bkghist type {type(bkghist)} \n hMC type - {type(hMC)} ")
    for histo, label, bkg, percent in zip(reversed(bkghist), reversed(labels), reversed(bkgcontrib), reversed(percentOfAllhisto)):
            # print(f"working on the legend values - histo - {histo}, labels - {label}, bkgcontrib - {bkg}, percentofAllhisto - {percent:.2f}")
            if normBkg:
                leg.AddEntry(histo, f"{label}({percent:.2f}%)", "L")
            else:
                leg.AddEntry(histo, f"{label}({percent:.2f}%)", "F")


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
        for bkg in bkghist:
            bkg.Draw("hist same")
            
    #     # hMC.Draw("hist same")
    #     # leg.Clear()
        # leg.AddEntry(hMC, "Total Background", "L")
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
                print("bkghist = ",)
                # setupAxes(dummy, 1.2, 1.6, 0.045, 0.045, 0.045, 0.045)
                ratio = createRatio(datahist,hMC,xTitle)
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
    print(f"hemPeriod = {hemPeriod},  year = {year}")
    # lumi = "59.7"
    if year == "2017":
        lumi = "41.5"
    elif year == "2016":
        lumi = "35.9"
    elif hemPeriod == "PostHEM" and year == "2018":
        lumi = "38.7"
    elif hemPeriod == "PreHEM" and year == "2018":
        lumi = "21.1"
    elif year == "2018" and hemPeriod == False:
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


def plotABCDSingle(data, histoName, cuts, region, outputPath="./", title="", xTitle="", yTitle="", met_min=0, met_max=2000, tagger_min=0.0, tagger_max=1.0, isLogY=False, norm=False, normBkg=False, onlySig=False, stList=None, yieldValues=None, isRatio=False):
    ROOT.TH2.AddDirectory(False)
    print(f"region working on is - {title} and the cuts are - {cuts}")
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

    leg = ROOT.TLegend(0.17, 0.7, 0.95, 0.88)
    #nColumns = 3 if(len(data[1]) >= 3) else 1
    nColumns = 2
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetLineWidth(1)
    leg.SetNColumns(nColumns)
    leg.SetTextFont(42)
    ROOT.gStyle.SetLegendTextSize(0.024)

    hs = ROOT.THStack()
    hMC = None
    firstPass = True
    bkghist = []
    bkgcontrib = []
    labels= []
    # setup background histos
    # print("before setup bkg histos detail key = {}, {}, {}".format(histoName,xmin, xmax))
    for d in data[1]:
        # if (stList != None) and (not normBkg):
        #     newEntry = stList + [getLabel(d.legEntry()),round(h.Integral())]
        #     yieldValues.loc[len(yieldValues.index)] = newEntry
            # print(newEntry)
        print(f"the dataset is {d}, and label is {d.label_}, the len of the cuts = {len(cuts)}")
        # print(d.legEntry())
        # h = ROOT.TH1F(f"h_{d.label_}",f"{d.label_}",len(cuts),0,len(cuts))
        bkghist.append(ROOT.TH1F(f"h_{d.label_}",f"{d.label_}",len(cuts),0,len(cuts)) ) 
        for cutname in cuts: # loop over the the SVJbins 
            histo2DName = histoName + cutname 

            numberOfEvents = d.get2DHistoIntegral(histo2DName, xmin=met_min, xmax=met_max, ymin=tagger_min, ymax=tagger_max, showEvents=True) # using the Integral method to find the number of events
            # h.Fill(cutname,numberOfEvents)
            bkghist[-1].Fill(cutname,numberOfEvents)
            # print(f"histo name = {histo2DName}, num of Events calculated = {numberOfEvents}")
        # if (stList != None) and (not normBkg):
        #     newEntry = stList + [getLabel(d.legEntry()),round(numberOfEvents)]
        #     yieldValues.loc[len(yieldValues.index)] = newEntry
        bkghist[-1].SetFillColor(d.getColor())
        hs.Add(copy.deepcopy(bkghist[-1]))
        labels.append(d.label_)
        bkgcontrib.append(round(bkghist[-1].Integral()))

        # leg.AddEntry(bkghist[-1], d.label_+f"({round(bkghist[-1].Integral())})", "F")
        print(f"The value of legend entry here is - {d.legEntry()}")
        # bkghist += h
        if(firstPass):
            hMC = copy.deepcopy(bkghist[-1])
            # bkghist = hMC[-1]
            firstPass = False
        else:
            hMC.Add(copy.deepcopy(bkghist[-1]))
            # bkghist.Add(hMC[-1])
    # for hist in hMC:
    #     leg.AddEntry(h,)
    ymax=10**11
    ymin=10**-4
    lmax=10**12
    xmin = 0
    xmax =len(cuts)
    dummy = ROOT.TH1D("dummy", "dummy", 1000, hMC.GetBinLowEdge(1), hMC.GetBinLowEdge(hMC.GetNbinsX()) + hMC.GetBinWidth(hMC.GetNbinsX()))
    setupDummy(dummy, leg, "", xTitle, yTitle, isLogY, xmin, xmax, ymin, ymax, lmax, norm, normBkg,isRatio)
    dummy.Draw("hist")
    dummy.SetTitle(title)
    totalOfAllhisto = sum(bkgcontrib)
    percentOfAllhisto = [(x/totalOfAllhisto)*100 if totalOfAllhisto != 0 else x for x in bkgcontrib  ]
    for histo, label, bkg, percent in zip(reversed(bkghist), reversed(labels), reversed(bkgcontrib), reversed(percentOfAllhisto)):
            print(f"working on the legend values - histo - {histo}, labels - {label}, bkgcontrib - {bkg}, percentofAllhisto - {percent:.2f}")
            leg.AddEntry(histo, f"{label}({bkg}, {percent:.2f}%)", "F")
    hs.Draw("hist F same")
    leg.Draw("same")
    linestylenumber = 0
    linestyle = [ROOT.kSolid,ROOT.kDashed,ROOT.kDotted]
    if(data[2]):
        #firstPass=True
        for d in data[2]:
            print(f"the dataset is {d}, and label is {d.label_}, the len of the cuts = {len(cuts)}")
            h = ROOT.TH1F(f"h_{d.label_}",f"{d.label_}",len(cuts),0,len(cuts))  
            for cutname in cuts: # loop over the the SVJbins 
                histo2DName = histoName + cutname 
                numberOfEvents = d.get2DHistoIntegral(histo2DName, xmin=met_min, xmax=met_max, ymin=tagger_min, ymax=tagger_max, showEvents=True) # using the Integral method to find the number of events
                h.Fill(cutname,numberOfEvents)
                print(f"histo name = {histo2DName}, num of Events calculated = {numberOfEvents}")
        
            if (stList != None) and (not normBkg):
                newEntry = stList + [getLabel(d.legEntry()),round(h.Integral())]
                yieldValues.loc[len(yieldValues.index)] = newEntry
                # print("Signal = ",newEntry)
            # h.SetLineStyle(ROOT.kDashed)
            sig = round(simpleSig(h, hMC),2)
            h.SetLineStyle(linestyle[linestylenumber%3] )
            linestylenumber+=1
            h.SetLineWidth(3)
            leg.AddEntry(h, d.label_+" ({}), {}".format(round(h.Integral()),sig), "L")
            print(f"The value of legend entry here is - {d.legEntry()}")
        
            h.Draw("hist same")
    #         # history.append(h)
    # print("legend - ",leg)
    # leg.Draw()
    regionName=title.split(":")   
    c1.Update()
    c1.RedrawAxis() 
    c1.SaveAs(outputPath+"/"+histoName+region+regionName[0]+".png")
    c1.Close()
    del c1
    del leg
    del hMC

## TODO make a function to define the 4regions using the TPad for better positioning.


def plotABCD(data, histoName, maincut, cuts, ABCDregions, outputPath="./", xTitle="", yTitle="", isLogY=False, norm=False, normBkg=False, onlySig=False, stList=None, SVJbinContent=None, isRatio=False,year=2018,scenario = "d0_wp7_p0i0",optionLL=False):
    ROOT.TH2.AddDirectory(False)
    canvas = ROOT.TCanvas("canvas", "Stack Plot", 1200, 600)
    canvas.Divide(len(ABCDregions),1,0.0,0.0)
    # TODO : Make it compatible with ratio plots
    met = ABCDregions[0][1]
    dnn_score = ABCDregions[0][3]
    print(f"**********    Working on the ABCD region with the met - {met} and dnn score - {dnn_score}    ***************")

    # need to use these list because otherwise it only prints on the last histogram, some memory issue with python?
    stacks = []
    signalhist = []
    contributions = []   
    dummy = []
    leg = []
    bkghist = [[] for _ in range(len(ABCDregions))]
    regionbkgsum = []
    for iRegion,(region, met_min, met_max, tagger_min, tagger_max) in enumerate(ABCDregions):
        
        # print(iRegion)
        canvas.cd(iRegion+1)
        ROOT.gPad.Clear()
        if(iRegion==0):
            ROOT.gPad.SetLeftMargin(0.16)
        if(iRegion==3):
            ROOT.gPad.SetRightMargin(0.05)
        if(iRegion!=0):
            ROOT.gPad.SetLeftMargin(0.001)

        ROOT.gPad.SetTopMargin(0.0)
        ROOT.gPad.SetBottomMargin(0.15)
        # ROOT.gPad.SetTicks(0)
        ROOT.gPad.SetLogy(isLogY)
        # ROOT.gPad.SetFrameLineWidth(2)
        ROOT.gPad.SetGrid()
        # ROOT.SetTitleSize(0.1,"XY")
        
        
        # leg = ROOT.TLegend(0.1, 0.78, 0.98, 0.98)
        # leg.append(ROOT.TLegend(0.16, 0.78, 0.98, 0.98))
        if (iRegion==0):
            leg.append(ROOT.TLegend(0.16, 0.78, 0.98, 0.98))
        else:
            leg.append(ROOT.TLegend(0.02, 0.78, 0.98, 0.98))
        # leg.append(ROOT.TLegend(0.17, 0.7, 0.95, 0.88)) 
        #nColumns = 3 if(len(data[1]) >= 3) else 1
        nColumns = 1
        leg[-1].SetFillStyle(0)
        leg[-1].SetBorderSize(0)
        leg[-1].SetLineWidth(1)
        leg[-1].SetNColumns(nColumns)
        leg[-1].SetTextFont(42)
        ROOT.gStyle.SetLegendTextSize(0.06)

        
        stacks.append(ROOT.THStack(f"{iRegion}",f"{iRegion}"))
        print(f"***** Working on region {iRegion} *****")
        hMC = None
        # bkghist = []
        firstPass = True
        labels = []
        bkgcontrib = []
        # xmin = 0
        # xmax = len(cuts)
        # binwidth = len(cuts)
        xmin = 1
        xmax = 5
        binwidth = 4
        regionsum=[]
        # Setup background histo
        for d in data[1]:
            
        # if ('QCD' not in d.fileName) and ('ZJets' not in d.fileName):
        #     # print(f"the dataset is {d}, and label is {d.label_}, the len of the cuts = {len(cuts)}")
            # print(d.legEntry())
            # h = ROOT.TH1F(f"h_{d.label_}",f"{d.label_}",len(cuts),0,len(cuts))  
            bkghist[iRegion].append(ROOT.TH1F(f"h_{d.label_}",f"{d.label_}",binwidth,xmin,xmax))
            for cutname in cuts: # loop over the the SVJbins 
                histo2DName = histoName + cutname 

                numberOfEvents = d.get2DHistoIntegral(histo2DName, xmin=met_min, xmax=met_max, ymin=tagger_min, ymax=tagger_max, showEvents=True) # using the Integral method to find the number of events
                bkghist[iRegion][-1].Fill(cutname,numberOfEvents)
                if (stList!=None):
                    newEntry = stList + [region,d.label_,cutname,numberOfEvents]
                    SVJbinContent.loc[len(SVJbinContent.index)]=newEntry
                # print(f"histo name = {histo2DName}, num of Events calculated = {numberOfEvents}")
            
            bkghist[iRegion][-1].SetFillColor(d.getColor())
            bkghist[iRegion][-1].SetFillStyle(3001)
            labels.append(d.label_)
            regionsum.append(bkghist[iRegion][-1].Integral())
            bkgcontrib.append(round(bkghist[iRegion][-1].Integral()))
            stacks[-1].Add(copy.deepcopy(bkghist[iRegion][-1])) # Adding the histogram to the Stack
            # leg[-1].AddEntry(bkghist[-1], d.label_+f"({round(bkghist[-1].Integral())})", "F")
            # contributions.append((histoName,region,d.label_,round(bkghist[-1].Integral())))
            # bkgcontrib.append(round(bkghist[-1].Integral()))
            if firstPass:
                hMC = copy.deepcopy(bkghist[iRegion][-1])
                # bkghist = h
                firstPass = False
            else:
                # bkghist.Add(h)
                hMC.Add(copy.deepcopy(bkghist[iRegion][-1]))
        # print(f"Bkghist list = {bkghist}")
        regionbkgsum.append(regionsum)
        # sumContrib = sum(bkgcontrib)
        totalOfAllhisto = sum(bkgcontrib)
        # numbkgContrib = np.array(bkgcontrib)
        # sumContrib = np.sum(bkgcontrib)
        percentOfAllhisto = [(x/totalOfAllhisto)*100 if totalOfAllhisto != 0  else x for x in bkgcontrib  ]
        # np.set_printoptions(precision=3)
        
        for histo, label, bkg, percent in zip(reversed(bkghist[iRegion]), reversed(labels), reversed(bkgcontrib), reversed(percentOfAllhisto)):
            print(f"working on the legend values - histo - {histo}, labels - {label}, bkgcontrib - {bkg}, percentofAllhisto - {percent:.2f}")
            leg[-1].AddEntry(histo, f"{label}({bkg}, {percent:.2f}%)", "F")

        ymax=10**10
        ymin=10**-4
        lmax=10**11
        # xmin = 0
        # xmax = len(cuts)
        # xmin = 1
        # xmax =len(cuts)
        if(iRegion==3):
            xTitle = f"nSVJ({met},{dnn_score})"
        else:
            xTitle = ""
        # dummy.append(ROOT.TH1D(f"dummy_{iRegion}", f"dummy{iRegion}", 1000, hMC.GetBinLowEdge(1), hMC.GetBinLowEdge(hMC.GetNbinsX()) + hMC.GetBinWidth(hMC.GetNbinsX())))
        dummy.append(ROOT.TH1D(f"dummy_{iRegion}", f"dummy{iRegion}", binwidth, xmin, xmax))
        setupDummy(dummy[-1], leg[-1], "", xTitle, yTitle, isLogY, xmin, xmax, ymin, ymax, lmax, title=region,isABCD=True)
        # for i,cutname in enumerate(cuts):
        #     dummy[-1].GetXaxis().SetBinLabel(i+1,cutname)
        dummy[-1].Draw("hist")
        # dummy[-1].SetTitle(region)
        # dummy[-1].GetXaxis().SetTickSize(0.1)
        # dummy[-1].GetYaxis().SetTickSize(0.1)
        
        # dummy[-1].SetTitleOffset(0.02)
        # dummy[-1].SetTitleSize(0.4)
        if(iRegion!=0):
            dummy[-1].GetYaxis().SetLabelSize(0)
            dummy[-1].GetYaxis().SetTitleSize(0)
        # dummy[-1].SetTitleSize()
        # ROOT.gStyle
        print("Drawing the Stack plot")
        # stacks[-1].SetTitle(region)
        stacks[-1].Draw("hist F same text")
        

        # canvas.Update()
        # Setup signal histo
        linestylenumber = 0
        linestyle = [ROOT.kSolid,ROOT.kDashed,ROOT.kDotted]
        signaltoPlot = ["baseline","mMed_600","mMed_4000"]
        if(data[2]):
            for d in data[2]:
                if ('QCD' not in d.fileName) and ('ZJets' not in d.fileName):
                # print(f"the dataset is {d}, and label is {d.label_}, the len of the cuts = {len(cuts)}")
                    signalhist.append(ROOT.TH1F(f"h_{d.label_}",f"{d.label_}",binwidth,xmin,xmax) )  
                
                
                    for cutname in cuts: # loop over the the SVJbins 
                        histo2DName = histoName + cutname 
                        numberOfEvents = d.get2DHistoIntegral(histo2DName, xmin=met_min, xmax=met_max, ymin=tagger_min, ymax=tagger_max, showEvents=True) # using the Integral method to find the number of events
                        signalhist[-1].Fill(cutname,numberOfEvents)
                        if (stList!=None):
                            newEntry = stList + [region,d.label_,cutname,numberOfEvents]
                            SVJbinContent.loc[len(SVJbinContent.index)]=newEntry
                        print(f"histo name = {histo2DName}, num of Events calculated = {numberOfEvents}")
                
                    
                    
                    
                    sig = round(simpleSig(signalhist[-1], hMC),2)
                    # print(f"The value of legend entry here is - {d.legEntry()}")
            
                    print("Drawing the signal plot")
                    # plot only few plots 
                    if d.label_ in signaltoPlot:
                        signalhist[-1].SetLineStyle(linestyle[linestylenumber%3] )
                        signalhist[-1].SetLineWidth(3)
                        linestylenumber+=1
                        signalhist[-1].SetLineColor(d.getColor())
                        leg[-1].AddEntry(signalhist[-1], d.label_+" ({}), {}".format(round(signalhist[-1].Integral()),sig), "L")
                        signalhist[-1].Draw("hist same")
                    # signalhist.append(copy.deepcopy(h))
                    # history.append(h)
            
        leg[-1].Draw("same")
        canvas.Update()
        ROOT.gPad.RedrawAxis()
        ROOT.gPad.RedrawAxis("G")
        # print("contribution value = {}",contributions)
        # print(f"The total in the {region} is {}")
        print(f"The region = {region}  --- the bkgcontrib = {bkgcontrib} and the sum = {totalOfAllhisto} and percentage = {percentOfAllhisto}")
        contributions.append((histoName,cutname,region,percentOfAllhisto))
    print("contribution value = {}",contributions)
    # print(f"the bkgcontrib variable is - {bkgcontrib} \n and the bkghist goes as - {bkghist}")
        


    # canvas.cd(1)
    CMS_lumi.cmsText = f"Scenario {scenario}"
    CMS_lumi.writeExtraText = False
    # CMS_lumi.extraText = "Simulation"
    # if year == "2017":
    #     lumi = "41.5"
    # elif year == "2016":
    #     lumi = "59.7"
    # else:
    #     lumi = "59.7"

    # CMS_lumi.lumi_sqrtS = lumi + " fb^{-1} (13 TeV)"

    iPeriod = 0
    iPos = 0

    CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)
    canvas.cd()
    canvas.Update()
    # canvas.RedrawAxis("G")

    print("saving the histogram")
    # cutSaveName = cuts.split("_")
    savestring = histoName+maincut+f"_MET{met}_DNN{dnn_score}"
    canvas.SaveAs(outputPath+"/"+savestring+"_forLL.png")
    # PredictionABCD(regionbkgsum)
    print("CLosing the canvas")
    canvas.Close()
    del canvas
    # del leg
    del hMC

def plotTransferFactors(data, histoName, SRcut, CRcuts,SVJbins, ABCDregions, outputPath="./", xTitle="", yTitle="", isLogY=False, norm=False, normBkg=False, onlySig=False, stList=None, TFContent=None, year=2018,scenario = "d0_wp7_p0i0",optionAVG=False,optionSum=False):
    ROOT.TH2.AddDirectory(False)
    canvas = ROOT.TCanvas(f"canvas", f"ABCD TF Plot", 1200, 600)
    canvas.Divide(len(ABCDregions),1,0.0,0.0)
    # TODO : Make it compatible with ratio plots
    met = ABCDregions[0][1]
    dnn_score = ABCDregions[0][3]
    print(f"**********    Working on the ABCD region with the met - {met} and dnn score - {dnn_score}    ***************")

    # need to use these list because otherwise it only prints on the last histogram, some memory issue with python?
    histSR = []
    histCR = []
    histTF = []
    ratioTF = []
    leg = []
    bkghist = [[] for _ in range(len(ABCDregions))]
    
    for iRegion,(region, met_min, met_max, tagger_min, tagger_max) in enumerate(ABCDregions):
        # print(iRegion)
        canvas.cd(iRegion+1)
        ROOT.gPad.Clear()
        ROOT.gPad.SetRightMargin(0)
        ROOT.gPad.SetLeftMargin(0)
        ROOT.gPad.SetTopMargin(0.0)
        ROOT.gPad.SetBottomMargin(0.15)
        ROOT.gPad.SetLogy(isLogY)
        ROOT.gPad.SetGrid()
        ROOT.gStyle.SetOptStat(0)
        
        leg.append(ROOT.TLegend(0.8, 0.8, 0.98, 0.98))
        nColumns = 1
        leg[-1].SetFillStyle(0)
        leg[-1].SetBorderSize(0)
        leg[-1].SetLineWidth(1)
        leg[-1].SetNColumns(nColumns)
        leg[-1].SetTextFont(42)
        ROOT.gStyle.SetLegendTextSize(0.06)

        print(f"***** Working on region {iRegion} *****")
        xmin = 1
        # xmax = len(SVJbins)
        binwidth = len(SVJbins)
        # xmin = int(SVJbins[0][0])
        xmax = 5
        print(f"xmin = {int(SVJbins[0][0])}")
        print(f"binwidth == {binwidth}")
        
        histSR.append(ROOT.TH1F(f"h_SR{iRegion}",f"SR{iRegion}",binwidth,xmin,xmax))
        histCR.append(ROOT.TH1F(f"h_CR{iRegion}",f"CR{iRegion}",binwidth,xmin,xmax))
        histTF.append(ROOT.TH1F(f"h_TF{iRegion}",f"TF{iRegion}",binwidth,xmin,xmax))
        
        print(f"cr = {CRcuts} , sr = {SRcut}")
        for SVJ in SVJbins:
            noEventsCR = []
            for CRcut in CRcuts:
                print(f"CR cut in the loop -- {CRcut}, and CRcuts = {CRcuts}")
                numberOfEvents_cr,numberOfEvents_sr = 0,0
                histo2DName_sr = histoName + SRcut + SVJ
                histo2DName_cr = histoName + CRcut + SVJ
                print(" =================================",histo2DName_cr)
                for d in data[1]:
                    if ('QCD' not in d.fileName) and ('ZJets' not in d.fileName): 
                        print("files used = ",d.fileName)
                        numberOfEvents_cr += d.get2DHistoIntegral(histo2DName_cr, xmin=met_min, xmax=met_max, ymin=tagger_min, ymax=tagger_max, showEvents=True) # using the Integral method to find the number of events
                        numberOfEvents_sr += d.get2DHistoIntegral(histo2DName_sr, xmin=met_min, xmax=met_max, ymin=tagger_min, ymax=tagger_max, showEvents=True) # using the Integral method to find the number of events
                noEventsCR.append(numberOfEvents_cr)
            histSR[-1].SetBinContent(int(SVJ[0]),numberOfEvents_sr)
            if optionAVG:
                averageCR = (noEventsCR[0]+noEventsCR[1])/2
                if averageCR!=0:
                    TFvalues = numberOfEvents_sr/averageCR
                else:
                    TFvalues = 0 
                histCR[-1].SetBinContent(int(SVJ[0]),averageCR) 
                newEntry = [region,SRcut,CRcut+"",SVJ,TFvalues]
            elif optionSum:
                sumCR = sum(noEventsCR)
                if sumCR!= 0:
                    TFvalues = numberOfEvents_sr/sumCR
                else:
                    TFvalues = 0
                histCR[-1].SetBinContent(int(SVJ[0]),sumCR) 
                newEntry = [region,SRcut,CRcut+"sum",SVJ,TFvalues]
            else:
                if numberOfEvents_cr!= 0:
                    TFvalues = numberOfEvents_sr/numberOfEvents_cr
                else:
                    TFvalues = 0
                histCR[-1].SetBinContent(int(SVJ[0]),numberOfEvents_cr)
                newEntry = [region,SRcut,CRcut,SVJ,TFvalues]
             
            
            TFContent.loc[len(TFContent.index)] = newEntry

        histSR[-1].SetLineColor(ROOT.kRed)
        histCR[-1].SetLineColor(ROOT.kBlue)
        histSR[-1].SetStats(0)
        histSR[-1].SetTitle("")
        print(histCR[-1].GetBinContent(1))  ## TODO : The bin content value is used as the bin error, need to fix this.
        print(histCR[-1].GetBinError(1))

        if(iRegion==3):
            histSR[-1].GetXaxis().SetTitle(f"nSVJ({met},{dnn_score})")
            histSR[-1].GetXaxis().SetTitleSize(1.5)
            histSR[-1].GetXaxis().SetTitleOffset(1)
        histSR[-1].GetXaxis().SetNdivisions(4)
        # setupAxes(histSR[-1],0.6, 0.8, 0.1, 0.1, 0.08, 0.06) 
        ratioTF.append(ROOT.TRatioPlot(histSR[-1],histCR[-1]))   
        leg[-1].AddEntry(histCR[-1],"CR","l")
        leg[-1].AddEntry(histSR[-1],"SR","l")
        ratioTF[-1].Draw()
        ratioTF[-1].SetH1DrawOpt("histe")
        ratioTF[-1].SetH2DrawOpt("histe")
        ratioTF[-1].GetLowerRefGraph().SetMarkerStyle(8)
        ymax=10**6
        ymin=10**-6
        ratioTF[-1].SetSplitFraction(0.4)
        # Ratio plot margin setup and size of the plots
        if iRegion == 0:
            ratioTF[-1].SetLeftMargin(0.18)
            ratioTF[-1].SetRightMargin(0)
            ratioTF[-1].GetUpperRefYaxis().SetTitle("Events")
            ratioTF[-1].GetUpperRefYaxis().SetTitleOffset(1.0)
            ratioTF[-1].GetUpperRefYaxis().SetTitleSize(0.08)
            ratioTF[-1].GetUpperRefYaxis().SetLabelSize(0.06)
            
            ratioTF[-1].GetLowerRefYaxis().SetTitle("SR/CR")
            ratioTF[-1].GetLowerRefYaxis().SetTitleOffset(1.4)
            ratioTF[-1].GetLowerRefYaxis().SetTitleSize(0.06)
            ratioTF[-1].GetLowerRefYaxis().SetLabelSize(0.06)
        
        else:
            ratioTF[-1].SetLeftMargin(0)
            ratioTF[-1].SetRightMargin(0)        
            ratioTF[-1].GetLowerRefYaxis().SetLabelSize(0)
            ratioTF[-1].GetUpperRefYaxis().SetLabelSize(0)
            ratioTF[-1].GetLowerRefYaxis().SetTitleSize(0)
            ratioTF[-1].GetUpperRefYaxis().SetTitleSize(0)
        ratioTF[-1].SetUpTopMargin(0)
        ratioTF[-1].SetUpBottomMargin(0)
        ratioTF[-1].SetLowTopMargin(0)
        ratioTF[-1].SetLowBottomMargin(0.4)
        # ratioTF[-1].GetXaxis().SetTitle(f"nSVJ({met},{dnn_score})")
        # if(iRegion==3):
        #     ratioTF[-1].GetLowerRefGraph().GetXaxis().SetTitleSize(0.1)
        #     ratioTF[-1].GetLowerRefGraph().GetXaxis().SetTitleOffset(0.7)
        # else:
        #     ratioTF[-1].GetLowerRefXaxis().SetTitle("")
        
        ratioTF[-1].GetLowerRefXaxis().SetLabelSize(0.06)    
        ratioTF[-1].SetSeparationMargin(0.0)
        
        # 
        ratioTF[-1].GetLowerRefXaxis().CenterLabels()
        # ratioTF[-1].GetLowerRefYaxis().SetRangeUser(0.0,20.0)
        if optionSum:
            ratioTF[-1].GetLowerRefGraph().SetMaximum(8.0)
        else:
            ratioTF[-1].GetLowerRefGraph().SetMaximum(16.0)
        
        ratioTF[-1].GetUpperRefYaxis().SetRangeUser(ymin,ymax)
        
        ratioTF[-1].GetUpperPad().cd()
        leg[-1].Draw()
        canvas.Update()
        ROOT.gPad.RedrawAxis()
        # ROOT.gPad.RedrawAxis("G")

    
    iPeriod = 0
    iPos = 0
    CMS_lumi.cmsText = f"Scenario {scenario}"
    CMS_lumi.writeExtraText = False
    CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)
    canvas.cd()
    canvas.Update()
    print("saving the histogram")
    # cutSaveName = cuts.split("_")
    if optionAVG:
        savestring = histoName+f"avgCR_{SRcut}_MET{met}_DNN{dnn_score}"
    elif optionSum:
        savestring = histoName+f"SumCR_{SRcut}_MET{met}_DNN{dnn_score}"
    else:
        savestring = histoName+f"{CRcut}_{SRcut}_MET{met}_DNN{dnn_score}"
    # savestring = histoName+f"avgCR_{SRCut}_MET{met}_DNN{dnn_score}"
    canvas.SaveAs(outputPath+"/"+savestring+".png")
    canvas.Close()
    del canvas
    
def PredictionABCD(data,histoName,cutname,SVJbins,ABCDregions,outputPath,scenario="d0_wp7_p0i0",isLogY=False):
    ROOT.TH2.AddDirectory(False)
    c1 = ROOT.TCanvas( "c", "c", 800, 700)
    c1.cd()
    ROOT.gPad.Clear()
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.08)
    ROOT.gPad.SetBottomMargin(0.12)
    ROOT.gPad.SetTicks(1,1)
    ROOT.gPad.SetLogy(isLogY)
    ROOT.gStyle.SetOptStat("")
    # ROOT.gStyle.SetErrorX(0)

    met = 150
    dnn_score = 0.7
    leg = ROOT.TLegend(0.7, 0.6, 0.75, 0.75)
    #nColumns = 3 if(len(data[1]) >= 3) else 1
    nColumns = 1
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetLineWidth(1)
    leg.SetNColumns(nColumns)
    leg.SetTextFont(42)
    ROOT.gStyle.SetLegendTextSize(0.06)
    xmin = 1
    xmax = 5
    binwidth = 4
    predictHist = ROOT.TH1F(f"predA","predA",binwidth,xmin,xmax)
    AHist = ROOT.TH1F(f"A","A",binwidth,xmin,xmax)
    for SVJ in SVJbins:
        a_totalEvents,b_totalEvents,c_totalEvents,d_totalEvents = 0,0,0,0
        histo2DName = histoName + cutname + SVJ
        for d in data[1]:
            # if 'QCD' in d.fileName:
            a_totalEvents += d.get2DHistoIntegral(histo2DName, xmin=met, xmax=5000, ymin=dnn_score, ymax=1.0, showEvents=True)
            b_totalEvents += d.get2DHistoIntegral(histo2DName, xmin=met, xmax=5000, ymin=0, ymax=dnn_score, showEvents=True)
            c_totalEvents += d.get2DHistoIntegral(histo2DName, xmin=0, xmax=met, ymin=dnn_score, ymax=1.0, showEvents=True)
            d_totalEvents += d.get2DHistoIntegral(histo2DName, xmin=0, xmax=met, ymin=0, ymax=dnn_score, showEvents=True)
        predicted_A_events = b_totalEvents*(c_totalEvents/d_totalEvents)
        print(f"Ahist = {a_totalEvents}, Bhist = {b_totalEvents}, Chist = {c_totalEvents}, dhist = {d_totalEvents}")
        print(f"predicted_A_events = {predicted_A_events}   Ahist = {a_totalEvents}")
        predictHist.SetBinContent(int(SVJ[0]),predicted_A_events)
        AHist.SetBinContent(int(SVJ[0]),a_totalEvents)
        closure = np.abs((a_totalEvents*d_totalEvents)-(b_totalEvents*c_totalEvents))/((a_totalEvents*d_totalEvents)+(b_totalEvents*c_totalEvents))
        print(f"------ Closure ==  {closure}  ------")

    ymax = 10**6
    ymin = 10**-4
    lmax = 10**6
    xTitle = "SVJ bins"
    yTitle = "Events"
    # dummy = ROOT.TH1D(f"dummy", f"dummy", binwidth, xmin, xmax)
    # setupDummy(dummy, leg, "", xTitle, yTitle, isLogY, xmin, xmax, ymin, ymax, lmax,title="")

    # dummy.Draw("hist")
    # predictHist.Draw("hist same")
    # AHist.Draw("hist same")
    predictHist.SetTitle("")
    predictHist.GetXaxis().SetTitle(f"nSVJ({met},{dnn_score})")
    predictHist.GetXaxis().SetNdivisions(4)
    ratio = ROOT.TRatioPlot(predictHist,AHist)
    ratio.Draw()
    ratio.SetH1DrawOpt("histe text")
    ratio.SetH2DrawOpt("same histe text")
    ratio.GetLowerRefGraph().SetMarkerStyle(8)
    ratio.GetLowerRefYaxis().SetTitle("Predict/MC")
    ratio.GetUpperRefYaxis().SetTitle("Events")
    ratio.GetUpperRefYaxis().SetRangeUser(ymin,ymax)
    # dummy.GetXaxis().SetNdivisions(4)
    predictHist.SetLineStyle(ROOT.kSolid)
    AHist.SetLineStyle(ROOT.kDashed)
    predictHist.SetLineColor(ROOT.kViolet)
    AHist.SetLineColor(ROOT.kRed)
    ratio.SetSplitFraction(0.3)
    ratio.GetLowerRefGraph().SetMaximum(2.5)
    ratio.GetLowerRefXaxis().CenterLabels()  
    ratio.SetSeparationMargin(0.0)  
    leg.AddEntry(predictHist,"Predicted A","L")
    leg.AddEntry(AHist,"MC A","L")
    ratio.GetUpperPad().cd()
    leg.Draw()
    c1.Update()
    ROOT.gPad.RedrawAxis()
    ROOT.gPad.RedrawAxis("G")

    # # iPeriod = 0
    # # iPos = 0
    # # CMS_lumi.cmsText = f"Scenario {scenario}"
    # # CMS_lumi.writeExtraText = False
    # # CMS_lumi.CMS_lumi(c1, iPeriod, iPos)
    # c1.cd()
    # c1.Update()
    savestring = histoName+cutname+"PredictionPlotRatio"
    c1.SaveAs(outputPath+"/"+savestring+".png")
    c1.Close()
    del c1



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
    parser.add_option('-w',                 dest='scenario',  type='string',   default='d0_w7p0i0',                               help="Scenario")
    parser.add_option(    '--hemPeriod',  dest='hemPeriod', type=str, default=False,  help='HEM period (PreHEM or PostHEM), default includes entire sample')
    options, args = parser.parse_args()
    scenario = options.scenario
    year = options.year
    hemPeriod = options.hemPeriod
    
    SVJbins = ["0SVJ","1SVJ","2SVJ","3SVJ","4PSVJ"]
    regions = ["_cr_muon_","_cr_electron_","_pre_dphimin_"]

    ABCDregions = [ ("A: MET > 150, event tagger > 0.7", 150, 2000, 0.7, 1.0),
                ("B: MET > 150, event tagger < 0.7", 150, 2000, 0, 0.7),
                ("C: MET < 150, event tagger > 0.7", 0, 150, 0.7, 1.0),
                ("D: MET < 150, event tagger < 0.7", 0, 150, 0, 0.7)
            ]

    # Define the MET range
    # met_values = list(range(150, 501, 50))
    met_values = [150]
    dnn_value = 0.7
    ABCDregions = []
    # ABCDFolderName = "ABCD"
    ABCDFolderName = "ABCD_withoutSVJ0"
    # ABCDFolderName = "ABCD_withoutSVJ0and1"
    # Create the ABCDregions list using a list comprehension
    for met in met_values:
        ABCD = [ (f"A: MET > {met}, event tagger > {dnn_value}", met, 2000, dnn_value, 1.0),
                (f"B: MET > {met}, event tagger < {dnn_value}", met, 2000, 0, dnn_value),
                (f"C: MET < {met}, event tagger > {dnn_value}", 0, met, dnn_value, 1.0),
                (f"D: MET < {met}, event tagger < {dnn_value}", 0, met, 0, 0.7)
        ]
        ABCDregions.append(ABCD)
        
    print("ABCD region is - ",ABCDregions)
    cutsImportant = ["_pre","_lcr_pre"]#,"_qual_2PJ_st_dphimin_nl","_qual_2PJ_st_dphimin_ll"]
    
    # for region in regions:
    #     cutsImportant.append(region)
    #     for bin in SVJbins:
    #         cutsImportant.append(region + bin)
    
    print(f"cutsImportant - {cutsImportant}")
    # Data, sgData, bgData = getData("condor/" + options.dataset + "/", 1.0, year)
    # Data, sgData, bgData = getData("DaskHadd/" + options.dataset + "/", 1.0, year)
    Data, sgData, bgData = getData( options.dataset + "/", 1.0, year)
   
    # print(sgData)
    #Data, sgData, bgData = getData("condor/MakeNJetsDists_"+year+"/", 1.0, year)
    allRocValues = pd.DataFrame(columns=["cut","var","sig","bkg","roc_auc","cutDir","cutSig","cBg","cSig","mBg_f","mSig_f"])
    yieldValues = pd.DataFrame(columns=["cut","var","source","yield"])
    signifValues = pd.DataFrame(columns=["cut","var","source","max signif."])
    SVJbinContent = pd.DataFrame(columns=["cut","var","ABCDregion","source","SVJbin","content"])
    TFContent = pd.DataFrame(columns=["ABCDregion","SR","CR","SVJbin","content"])
    if options.outputdir:
        plotOutDir = "outputPlots/{}".format(options.outputdir)
    else: 
        plotOutDir = "outputPlots/{}".format(options.dataset)

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

    ABCDhistoVars = ["METvsDNN"]
    
    # myvars = key : ["xlabel", no. of bins, xmin,xmax, npzinfo, flattenInfo, weightName]
    
    for histName,details in vars(options.jNVar).items():
        # print("HistNames - {}".format(histName))
    # # for histName, details in myVars.items():
        print(histName)
        print(details)
        isNorm = options.isNorm
        isNormBkg = options.isNormBkg
        onlySig = options.onlySig
        manySigs = options.manySigs
        if histName in varsSkip:
            continue
        # if details[6] != "evtw":
        #    continue
        for cut in cutsImportant:
            print(f"cut - {cut}")
            makeDirs(plotOutDir,cut,"Stacked")
            makeDirs(plotOutDir,cut,"roc")
            makeDirs(plotOutDir,cut,"FOM")
            makeDirs(plotOutDir,cut,"NormedStacked")
            stList = [cut,histName]
            print("Data = ",Data)
            # plotROC(  (Data, bgData, sgData), "h_"+histName+cut, plotOutDir+"/roc/"+cut[1:], isLogY=False,   manySigs=manySigs, stList=stList, allRocValues=allRocValues)
            # plotStack((Data, bgData, sgData), "h_"+histName+cut, details[1], plotOutDir+"/Stacked/"+cut[1:], details[0], "Events", isLogY=True, norm=isNorm, xmin=details[2], xmax=details[3], normBkg=False, onlySig=onlySig, stList=stList, yieldValues=yieldValues,isRatio=False)
            plotStack((Data, bgData, sgData), "h_"+histName+cut, details[1], plotOutDir+"/Stacked/"+cut[1:], details[0], "Events", isLogY=True,norm=isNorm, xmin=details[2], xmax=details[3], normBkg=False, onlySig=onlySig, stList=stList, yieldValues=yieldValues,year= year, isRatio=True, hemPeriod = hemPeriod)
            # plotStack((Data, bgData, sgData), "h_"+histName+cut, details[1], plotOutDir+"/NormedStacked/"+cut[1:], details[0], "Events", isLogY=True, norm=isNorm, xmin=details[2], xmax=details[3], normBkg=True, onlySig=onlySig, stList=stList, yieldValues=yieldValues)
    # #         # if histName in preVars.keys():
    # # #             plotSignificance((Data, bgData, sgData), "h_"+histName, details[1], details[0], plotOutDir, cut,                    isLogY=False, reverseCut=preVars[histName], signifValues=signifValues)
    # for i,met in enumerate(met_values):    
    #     for histName in ABCDhistoVars:
    #         for region in regions:
    #             makeDirs(plotOutDir,region,ABCDFolderName)
    #             cuts = [region + SVJ for SVJ in SVJbins]
    #             current_abcd_regions = ABCDregions[i]
    #             print(f"The met value is - {met}      current ABCDregion - {current_abcd_regions}")
    #             stList = [region,histName]
    #             PredictionABCD((Data,bgData,sgData),"h_"+histName,region,SVJbins,ABCDregions,plotOutDir+"/"+ABCDFolderName,scenario=scenario,isLogY=True)
    #             # plotABCD((Data,bgData,sgData),"h_"+histName,region,cuts,current_abcd_regions,outputPath=plotOutDir+"/"+ABCDFolderName+"/"+region[1:],xTitle="nSVJ",yTitle="Events",isLogY=True,stList=stList,SVJbinContent=SVJbinContent,year=year,scenario=scenario)
                
    #         # for title, met_min, met_max, tagger_min, tagger_max in ABCDregions:
    #             # plotABCDSingle((Data,bgData,sgData),"h_"+histName,cuts,region,outputPath=plotOutDir+"/ABCD/", title=title,xTitle="nSVJ",yTitle="number of Events",met_min=met_min,met_max=met_max, tagger_min=tagger_min, tagger_max=tagger_max, isLogY=True)
    yieldValues.to_csv("{}/yieldValues.csv".format(plotOutDir))
    # # signifValues.to_csv("{}/signifValues.csv".format(plotOutDir))
    # # allRocValues.to_csv("{}/allRocValues.csv".format(plotOutDir))
    # SVJbinContent.to_csv("{}/{}/allSVJBinContents.csv".format(plotOutDir,ABCDFolderName))
    # # print(yieldValues)
    # SRCut = "_pre_dphimin_"
    # CRcuts = ["_cr_muon_","_cr_electron_"]
                
    # for histName in ABCDhistoVars:
    #     for CRcut in CRcuts:
    #         CRcut = [CRcut]
    #         plotTransferFactors((Data,bgData,sgData),"h_"+histName,SRCut,CRcut,SVJbins,ABCDregions[0],outputPath=plotOutDir+"/"+ABCDFolderName,xTitle="nSVJ",yTitle="Transfer Factors",isLogY=True,TFContent=TFContent,year=year,scenario=scenario)
    #     plotTransferFactors((Data,bgData,sgData),"h_"+histName,SRCut,CRcuts,SVJbins,ABCDregions[0],outputPath=plotOutDir+"/"+ABCDFolderName,xTitle="nSVJ",yTitle="Transfer Factors",isLogY=True,TFContent=TFContent,year=year,scenario=scenario,optionAVG=True)
    #     plotTransferFactors((Data,bgData,sgData),"h_"+histName,SRCut,CRcuts,SVJbins,ABCDregions[0],outputPath=plotOutDir+"/"+ABCDFolderName,xTitle="nSVJ",yTitle="Transfer Factors",isLogY=True,TFContent=TFContent,year=year,scenario=scenario,optionSum=True)
    # TFContent.to_csv("{}/{}/TFIndividualContent.csv".format(plotOutDir,ABCDFolderName))

if __name__ == '__main__':
    main()
