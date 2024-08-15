#!/usr/bin/env python

from coffea import processor
from utils import samples as s
import time
import numpy as np
from utils.python.lpc_dask import *
from magiconfig import ArgumentParser
from utils.data.DNNEventClassifier import configs as c
from utils.coffea.n_tree_maker_schema import NTreeMakerSchema
import joblib
import json
import torch
from pydoc import locate

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
    parser = ArgumentParser()
    parser.add_argument('-d', '--dataset',      help='dataset', dest='dataset', nargs="+", default=["2018_mMed-1000_mDark-20_rinv-0p3_alpha-peak_yukawa-1"])
    parser.add_argument('-j', '--jNVar',        help='make histograms for nth jet variables', dest='jNVar', default=False, action='store_true')
    parser.add_argument('-N', '--nFiles',       help='nFiles',            dest='nFiles',    type=int, nargs="+", default=[-1])
    parser.add_argument('-M', '--startFile',    help='startFile',         dest='startFile', type=int, nargs="+", default=[0])
    parser.add_argument(      '--condor',       help='running on condor', dest='condor',              default=False, action='store_true')
    parser.add_argument(      '--dask',         help='run w/ dask', dest='dask',              default=False, action='store_true')
    parser.add_argument(      '--port',         help='port for dask status dashboard (localhost:port)', dest='port', type=int, default=2542)
    parser.add_argument(      '--mincores',     help='dask waits for min # cores', dest='mincores', type=int, default=4)
    parser.add_argument(      '--quiet',        help='suppress status printouts', dest='quiet',              default=False, action='store_true')
    parser.add_argument('-w', '--workers',      help='Number of workers to use for multi-worker executors (e.g. futures or condor)', dest='workers', type=int, default=8)
    parser.add_argument('-b', '--jobs',         help='Number of workers to use for condor dask', dest='jobs', type=int, default=1)
    parser.add_argument('-s', '--chunksize',    help='Chunk size',        dest='chunksize', type=int, default=10000)
    parser.add_argument('-m', '--maxchunks',    help='Max number of chunks (for testing)',        dest='maxchunks', type=int, default=None)
    parser.add_argument('-t', '--eTagLoc',      help='Location of the event tagger model',        dest='eTagLoc', type=str, default="utils/data/DNNEventClassifier/sdt_QCD_disco_0p001_closure_0p02_damp_1_net_64_32_16_8_1Evt_pn")
    parser.add_argument(      '--outHistF',     help='Output directory for histogram files',      dest='outHistF', type=str, default="./")
    parser.add_argument(      '--hemPeriod',    help='HEM period (PreHEM or PostHEM), default includes entire sample',            dest='hemPeriod', type=str, default="")
    parser.add_argument(      '--hemStudy',     help='HEM study',         dest='hemStudy',             default=False, action='store_true')
    parser.add_argument(      '--slimProc',     help='Slimmed processor for fasting processing',       dest='slimProc',             default=False, action='store_true')  
    parser.add_argument(      '--tcut',         help='Cut for training files: _pre, _pre_1PSVJ',  dest='tcut', type=str, default="_pre")    
    parser.add_argument('-i', '--issues',       help='Run the dataTestProcessor', dest='issue', default=False, action='store_true')
    parser.add_argument('-f', '--sFactor',      help='Scale factor', dest='sFactor',  default=False, action='store_true')
    parser.add_argument('-o', '--outputFile',   help='Output file name ', dest='outputFile', default=False, type=str)
    parser.add_argument(      '--skimSource',   help='Use skim files instead of TreeMaker ntuples ', dest='skimSource', default=False, action='store_true')
    parser.add_argument(      '--skimCut',      help='The selection of cuts that have been applied to the TM ntuples when making the skims: t_channel_pre_selection, t_channel_lost_lepton_control_region, etc.',  dest='skimCut', type=str, default="t_channel_pre_selection")
    parser.add_argument(      '--runJetTag',    help='Run jet tagger.', dest='runJetTag',  default=False, action='store_true')
    parser.add_argument(      '--runEvtClass',  help='Run event classifier.', dest='runEvtClass',  default=False, action='store_true')


    options = parser.parse_args()
    ###########################################################################################################
    # set output root file
    ###########################################################################################################
    sample = options.dataset

    # getting dictionary of files from a sample collection e.g. "2016_QCD, 2016_WJets, 2016_TTJets, 2016_ZJets"
    fileset = s.getFilesetFromList(sample, options.startFile, options.nFiles, options.skimSource, options.skimCut, True)
    # sf = s.sfGetter(sample,True)
    # print("scaleFactor = {}".format(sf))

    ###########################################################################################################
    # get event level DNN information
    ###########################################################################################################
    eTagLoc = options.eTagLoc
    with open(f"{eTagLoc}/order_of_variables.py", "r") as fp:
        order_of_variables = json.load(fp)
    # load event classifier
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    checkpoint = torch.load(f"{eTagLoc}/model.pt", map_location=device )
    config_training_model = f"{eTagLoc}/training_model".replace("/",".")
    # utils.data.DNNEventClassifier.damp_1_0p001_closure_0p06_net_64_32_16_8_LJP_1EvtABCD_fixedSkimBugs.files_preparation_options
    opt_training_model = locate(config_training_model)
    model = opt_training_model.model 
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()
    evtTaggerDict = {
        "location":eTagLoc.replace("/","."),
        "scaler":joblib.load(f"{eTagLoc}/scaler.joblib"),
        "model": model,
    }
    ###########################################################################################################
    # get executor/processor args
    ###########################################################################################################
    if options.hemStudy:
        from processors.hemPhiSpikeProcessor import MainProcessor
    elif options.slimProc:
        from processors.slimProcessor import MainProcessor
    elif options.issue:
        from processors.dataTestProcessor import MainProcessor
    else:
        from processors.mainProcessor import MainProcessor

    MainExecutor = processor.futures_executor

    if options.dask:
        runProcessWithErrorHandling(fileset,sample,MainExecutor,MainProcessor,options,evtTaggerDict)
    else:
        schema = processor.TreeMakerSchema
        if options.skimSource:
            schema = NTreeMakerSchema
        exe_args = {
            'workers': options.workers, 
            'schema': schema
        }
        run_processor(fileset,sample,MainExecutor,MainProcessor,options,exe_args,evtTaggerDict)

    ###########################################################################################################
    # print run time in seconds
    ###########################################################################################################
    dt = time.time() - tstart
    print("run time: %.2f [sec]" % (dt))

if __name__ == "__main__":
    main()
