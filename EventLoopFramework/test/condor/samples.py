from ctypes import cdll
from ctypes import c_char_p
from ctypes import c_double
from ctypes import POINTER 
from ctypes import c_int

from os import environ

class SampleCollection:
    def __init__(self, ssfile, scfile):
        self.lib = cdll.LoadLibrary(environ["CMSSW_BASE"] + "../obj/samplesModule.so")
        self.ss = self.lib.SS_new(ssfile)
        self.obj = self.lib.SC_new(self.ss, scfile)
        self.lib.SC_samples.restype = POINTER(c_char_p)
        self.lib.SC_samples_names.restype = POINTER(c_char_p)
        self.lib.SC_samples_nEvts.restype = POINTER(c_int)
        self.lib.SS_samples.restype = POINTER(c_char_p)
        self.lib.SS_samples_names.restype = POINTER(c_char_p)
        self.lib.SC_samplecollection_names.restype = POINTER(c_char_p)
        self.lib.SC_samplecollection_lumis.restype = POINTER(c_double)
        self.lib.SC_fixed_lumi.restype = c_double        

    def nSamples(self, name):
        return self.lib.SC_samples_size(self.obj, name)

    def getFixedLumi(self):
        return self.lib.SC_fixed_lumi()

    def sampleList(self, name):
        files = self.lib.SC_samples(self.obj, name)
        names = self.lib.SC_samples_names(self.obj, name)
        nEvts = self.lib.SC_samples_nEvts(self.obj, name)
        list = [(files[i],names[i],nEvts[i]) for i in xrange(self.lib.SC_samples_size(self.obj, name))]
        return list

    def sampleCollectionList(self):
        names = self.lib.SC_samplecollection_names(self.obj)
        list = [names[i] for i in xrange(self.lib.SC_samplecollection_size(self.obj))]
        return list

    def sampleCollectionLumiList(self):
        names = self.lib.SC_samplecollection_names(self.obj)
        lumis = self.lib.SC_samplecollection_lumis(self.obj)
        list = [(names[i],lumis[i]) for i in xrange(self.lib.SC_samplecollection_size(self.obj))]
        return dict(list)

    def sampleSetList(self):
        names = self.lib.SS_samples_names(self.ss)
        files = self.lib.SS_samples(self.ss)
        list = [(names[i],files[i]) for i in xrange(self.lib.SS_samples_size(self.ss))]
        return list

class SampleSet:
    def __init__(self, ssfile):
        self.lib = cdll.LoadLibrary('../obj/samplesModule.so')
        self.ss = self.lib.SS_new(ssfile)
        self.lib.SS_samples.restype = POINTER(c_char_p)
        self.lib.SS_samples_names.restype = POINTER(c_char_p)

    def sampleSetList(self):
        names = self.lib.SS_samples_names(self.ss)
        files = self.lib.SS_samples(self.ss)
        list = [(names[i],files[i]) for i in xrange(self.lib.SS_samples_size(self.ss))]
        return list
