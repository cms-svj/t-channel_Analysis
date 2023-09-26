#!/usr/bin/env python

from coffea import processor
from utils import samples as s
import time
import numpy as np
from utils.python.lpc_dask import *
from magiconfig import ArgumentParser, MagiConfigOptions, ArgumentDefaultsRawHelpFormatter
from utils.data.DNNEventClassifier import configs as c

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
    parser = ArgumentParser(config_options=MagiConfigOptions(strict = False, default="utils/data/DNNEventClassifier/pre_nnOutput_v2/config_out.py"),formatter_class=ArgumentDefaultsRawHelpFormatter)
    parser.add_argument('-d', '--dataset',   help='dataset', dest='dataset', type=str, default="2018_mMed-1000_mDark-20_rinv-0p3_alpha-peak_yukawa-1")
    parser.add_argument('-j', '--jNVar',     help='make histograms for nth jet variables', dest='jNVar', default=False, action='store_true')
    parser.add_argument('-N', '--nFiles',    help='nFiles',            dest='nFiles',    type=int, default=-1)
    parser.add_argument('-M', '--startFile', help='startFile',         dest='startFile', type=int, default=0)
    parser.add_argument(      '--condor',    help='running on condor', dest='condor',              default=False, action='store_true')
    parser.add_argument(      '--dask',      help='run w/ dask', dest='dask',              default=False, action='store_true')
    parser.add_argument(      '--port',      help='port for dask status dashboard (localhost:port)', dest='port', type=int, default=8787)
    parser.add_argument(      '--mincores',  help='dask waits for min # cores', dest='mincores', type=int, default=4)
    parser.add_argument(      '--quiet',     help='suppress status printouts', dest='quiet',              default=False, action='store_true')
    parser.add_argument('-w', '--workers',   help='Number of workers to use for multi-worker executors (e.g. futures or condor)', dest='workers', type=int, default=8)
    parser.add_argument('-b', '--jobs',      help='Number of workers to use for condor dask', dest='jobs', type=int, default=1)
    parser.add_argument('-s', '--chunksize', help='Chunk size',        dest='chunksize', type=int, default=10000)
    parser.add_argument('-m', '--maxchunks', help='Max number of chunks (for testing)',        dest='maxchunks', type=int, default=None)
    parser.add_argument('-t', '--eTagLoc',   help='Location of the event tagger model',        dest='eTagLoc', type=str, default="utils/data/DNNEventClassifier/pre_nnOutput_v2")
    parser.add_argument(      '--outHistF',  help='Output directory for histogram files',      dest='outHistF', type=str, default="./")
    parser.add_argument(      '--hemPeriod', help='HEM period (PreHEM or PostHEM), default includes entire sample',            dest='hemPeriod', type=str, default="")
    parser.add_argument(      '--hemStudy',  help='HEM study',         dest='hemStudy',             default=False, action='store_true')
    parser.add_argument(      '--slimProc',  help='Slimmed processor for fasting processing',       dest='slimProc',             default=False, action='store_true')  
    parser.add_argument(      '--tcut',      help='Cut for training files: _pre, _pre_1PSVJ',  dest='tcut', type=str, default="_pre")    
    for arg in c.config_schema:
        parser.add_config_argument(arg)
    options = parser.parse_args()
    ###########################################################################################################
    # set output root file
    ###########################################################################################################
    sample = options.dataset

    # getting dictionary of files from a sample collection e.g. "2016_QCD, 2016_WJets, 2016_TTJets, 2016_ZJets"
    fileset = s.getFileset(sample, True, options.startFile, options.nFiles)
    sf = s.sfGetter(sample,True)
    print("scaleFactor = {}".format(sf))

    ###########################################################################################################
    # get event level NN information
    ###########################################################################################################
    hyper = options.hyper
    features = options.features
    evtTaggerDict = {
        "hyper":hyper,
        "features":features,
        "evtTaggerLocation":options.eTagLoc
    }
    ###########################################################################################################
    # get executor/processor args
    ###########################################################################################################
    if options.hemStudy:
        from processors.hemPhiSpikeProcessor import MainProcessor
    elif options.slimProc:
        from processors.slimProcessor import MainProcessor
    else:
        from processors.mainProcessor import MainProcessor

    MainExecutor = processor.futures_executor

    if options.dask:
        runProcessWithErrorHandling(fileset,sample,sf,MainExecutor,MainProcessor,options,evtTaggerDict)
    else:
        exe_args = {
            'workers': options.workers, 
            'schema': processor.TreeMakerSchema
        }
        run_processor(fileset,sample,sf,MainExecutor,MainProcessor,options,exe_args,evtTaggerDict)

    ###########################################################################################################
    # print run time in seconds
    ###########################################################################################################
    dt = time.time() - tstart
    print("run time: %.2f [sec]" % (dt))

if __name__ == "__main__":
    main()
