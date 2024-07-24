#!/usr/bin/env python

from coffea import processor
import uproot as up
import sys,os
from utils import samples as s
import time
from optparse import OptionParser
from glob import glob
import numpy as np
from utils.python.lpc_dask import *

def main():
    # start run time clock
    tstart = time.time()

    # get options from command line
    parser = OptionParser()
    parser.add_option('-d', '--dataset',   help='dataset',           dest='dataset')
    parser.add_option('-j', '--jNVar',     help='make histograms for nth jet variables', dest='jNVar', default=False, action='store_true')
    parser.add_option('-N', '--nFiles',    help='nFiles',            dest='nFiles',    type=int, default=-1)
    parser.add_option('-M', '--startFile', help='startFile',         dest='startFile', type=int, default=0)
    parser.add_option(      '--condor',    help='running on condor', dest='condor',              default=False, action='store_true')
    parser.add_option(      '--dask',      help='run w/ dask', dest='dask',              default=False, action='store_true')
    parser.add_option(      '--port',      help='port for dask status dashboard (localhost:port)', dest='port', type=int, default=8787)
    parser.add_option(      '--mincores',  help='dask waits for min # cores', dest='mincores', type=int, default=4)
    parser.add_option(      '--quiet',     help='suppress status printouts', dest='quiet',              default=False, action='store_true')
    parser.add_option('-w', '--workers',   help='Number of workers to use for multi-worker executors (e.g. futures or condor)', dest='workers', type=int, default=8)
    parser.add_option('-s', '--chunksize', help='Chunk size',        dest='chunksize', type=int, default=10000)
    parser.add_option('-m', '--maxchunks', help='Max number of chunks (for testing)',        dest='maxchunks', type=int, default=None)
    parser.add_option('-b', '--jobs',      help='Number of workers to use for condor dask', dest='jobs', type=int, default=1)
    parser.add_option(      '--outHistF',  help='Output directory for histogram files',      dest='outHistF', type=str, default="./")
    parser.add_option(      '--hemPeriod', help='HEM period (PreHEM or PostHEM), default includes entire sample',            dest='hemPeriod', type=str, default="")
    parser.add_option(      '--tcut',      help='Cut for training files: _pre, _pre_1PSVJ',  dest='tcut', type=str, default="_pre")    
    parser.add_option('-z', '--eth',       help='Use trained model from eth', dest='eth',  default=False, action='store_true')
    parser.add_option('-f', '--sFactor',   help='Scale factor', dest='sFactor',  default=False, action='store_true')
    parser.add_option(      '--training',  help='For which training should the files be made: NN, PN',           dest='training', default="PN")
    parser.add_option(      '--skimSource',help='Use skim files instead of TreeMaker ntuples ', dest='skimSource', default=False, action='store_true')
    options, args = parser.parse_args()

    # set output root file
    sample = options.dataset
    # training kind
    trainingKind = options.training
    if trainingKind == "NN":
        from processors.rootProcessor_NN import MainProcessor
    elif trainingKind == "PN":
        # from processors.rootProcessor_varModule import MainProcessor
        from processors.trainFileProcessor import MainProcessor
    # getting dictionary of files from a sample collection e.g. "2016_QCD, 2016_WJets, 2016_TTJets, 2016_ZJets"
    fileset = s.getFileset(sample, options.startFile, options.nFiles, options.skimSource, True)

    # run processor
    MainExecutor = processor.futures_executor
    if options.dask:
        runProcessWithErrorHandling(fileset,sample,MainExecutor,MainProcessor,options,trainingKind=trainingKind,trainFileProduction=True)
    else:
        exe_args = {
            'workers': options.workers, 
            'schema': processor.TreeMakerSchema
        }
        run_processor(fileset,sample,MainExecutor,MainProcessor,options,exe_args,trainingKind=trainingKind,trainFileProduction=True)

    dt = time.time() - tstart
    print("run time: %.2f [sec]" % (dt))

if __name__ == "__main__":
    main()

#python analyze_npz.py -d [filename] -w 2 -s 10000
