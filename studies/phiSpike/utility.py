import uproot as up 
import mplhep as hep
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import coffea
import os
from hist import Hist
from scipy.stats import linregress
import pickle

def getHisto(inputFolder,fileName,varName,cut,dim=1):
    fullInputFilePath = "{}{}".format(inputFolder,fileName)
    f = up.open(fullInputFilePath)
    hist = f[varName + cut]
    if dim == 1:
        npHist, bins = hist.to_numpy() 
        histHist = hist.to_hist() #*(1./np.sum(npHist)) # normalization, however, ratioplot will not calculate the error properly
        return npHist, histHist, bins
    elif dim == 2:
        npHist = hist.to_numpy()
        histHist = hist.to_hist() #*(1./np.sum(npHist)) # normalization, however, ratioplot will not calculate the error properly
        return npHist, histHist

def getTotal2DHist(inputFolder,inputFile,varName,cut):
    npHist, histHist = getHisto(inputFolder,inputFile,varName,cut,dim=2)
    npZ = npHist[0]
    npX = npHist[1]
    npY = npHist[2]
    npZList = []
    for row in npZ:
        npZList.append(list(row))
    return npZ,npX,npY

def getTicks(fullTicks, fullTickLabels, desiredLabels):
    linFit = linregress(fullTickLabels, fullTicks)
    m = linFit[0]
    b = linFit[1]
    desiredTicks = []
    for lab in desiredLabels:
        desiredTicks.append(m*lab+b)
    return desiredTicks, linFit

def actualToBinVal(val, linFit):
    m = linFit[0]
    b = linFit[1]
    return m*val+b

def bintoActualVal(binVal, linFit):
    m = linFit[0]
    b = linFit[1]
    return (binVal-b)/m

def widenHEMMask(linFitX, linFitY, zVals, extraMask=2):
    # extraMask = number of extra bins to widen the HEM veto by
    etaVetoMin = int(actualToBinVal(-3.05, linFitX) + 0.5)
    etaVetoMax = int(actualToBinVal(-1.35, linFitX) + 0.5)
    phiVetoMax = int(actualToBinVal(-0.82, linFitY) + 0.5) - extraMask
    phiVetoMin = int(actualToBinVal(-1.62, linFitY) + 0.5) + extraMask
    zVals[phiVetoMax:phiVetoMin,etaVetoMin:etaVetoMax] = 0

def uniqueCombine(a1,a2,b1,b2):
    c1 = list(a1) + list(b1)
    c2 = list(a2) + list(b2)
    s = []
    for i in range(len(c1)):
        s.append(f"{c1[i]}_{c2[i]}")
    ulist, uInds = np.unique(s,return_index=True)
    u1 = []
    u2 = []
    for i in uInds:
        u1.append(c1[i])
        u2.append(c2[i])
    return u1, u2

def rebin2D(zVals,xVals,yVals,xRebin,yRebin):
    # rebin in y-direction
    yReducedZVals = np.add.reduceat(zVals,np.arange(0,zVals.shape[0],yRebin),axis=0)
    # rebin in x-direction
    xyReducedZVals = np.add.reduceat(yReducedZVals,np.arange(0,zVals.shape[1],xRebin),axis=1)
    rebinNpx = []
    for i in range(0,len(xVals),xRebin):
        rebinNpx.append(xVals[i])
    rebinNpy = []
    for i in range(0,len(yVals),yRebin):
        rebinNpy.append(yVals[i])  
    return xyReducedZVals, np.array(rebinNpx), np.array(rebinNpy)