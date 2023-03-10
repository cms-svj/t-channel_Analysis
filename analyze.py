#!/usr/bin/env python

from coffea import processor
from processors.mainProcessor import MainProcessor
from utils import samples as s
import time
from optparse import OptionParser
import numpy as np
from utils.python.lpc_dask import *

def main():
    ###########################################################################################################
    # Removing numpy warnings; might not be a good idea
    ###########################################################################################################
    np.seterr(all='ignore')

    # start run time clock
    tstart = time.time()

    ###########################################################################################################
    # get options from command line
    ###########################################################################################################
    parser = OptionParser()
    parser.add_option('-d', '--dataset',   help='dataset', dest='dataset', type=str, default="2018_mMed-1000_mDark-20_rinv-0p3_alpha-peak_yukawa-1")
    parser.add_option('-j', '--jNVar',     help='make histograms for nth jet variables', dest='jNVar', default=False, action='store_true')
    parser.add_option('-N', '--nFiles',    help='nFiles',            dest='nFiles',    type=int, default=-1)
    parser.add_option('-M', '--startFile', help='startFile',         dest='startFile', type=int, default=0)
    parser.add_option(      '--condor',    help='running on condor', dest='condor',              default=False, action='store_true')
    parser.add_option(      '--dask',      help='run w/ dask', dest='dask',              default=False, action='store_true')
    parser.add_option(      '--port',      help='port for dask status dashboard (localhost:port)', dest='port', type=int, default=8787)
    parser.add_option(      '--mincores',  help='dask waits for min # cores', dest='mincores', type=int, default=4)
    parser.add_option(      '--quiet',     help='suppress status printouts', dest='quiet',              default=False, action='store_true')
    parser.add_option('-w', '--workers',   help='Number of workers to use for multi-worker executors (e.g. futures or condor)', dest='workers', type=int, default=8)
    parser.add_option('-b', '--jobs',      help='Number of workers to use for condor dask', dest='jobs', type=int, default=100)
    parser.add_option('-s', '--chunksize', help='Chunk size',        dest='chunksize', type=int, default=10000)
    parser.add_option('-m', '--maxchunks', help='Max number of chunks (for testing)',        dest='maxchunks', type=int, default=None)
    options, args = parser.parse_args()

    ###########################################################################################################
    # set output root file
    ###########################################################################################################
    sample = options.dataset

    # getting dictionary of files from a sample collection e.g. "2016_QCD, 2016_WJets, 2016_TTJets, 2016_ZJets"
    fileset = s.getFileset(sample, True, options.startFile, options.nFiles)
    sf = s.sfGetter(sample,False)
    print("scaleFactor = {}".format(sf))

    ###########################################################################################################
    # get executor/processor args
    ###########################################################################################################
    MainExecutor = processor.futures_executor(workers=options.workers)

    if options.dask:
        runProcessWithErrorHandling(fileset,sample,sf,MainExecutor,MainProcessor,options)
    else:
        run_processor(fileset,sample,sf,MainExecutor,MainProcessor,options)

    ###########################################################################################################
    # print run time in seconds
    ###########################################################################################################
    dt = time.time() - tstart
    print("run time: %.2f [sec]" % (dt))

if __name__ == "__main__":
    main()
