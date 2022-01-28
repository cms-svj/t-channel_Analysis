#!/usr/bin/env python

from coffea import hist, processor
from processors.mainProcessor import MainProcessor
import uproot
import sys,os
from utils import samples as s
import time
from optparse import OptionParser
from glob import glob
import numpy as np
from magiconfig import MagiConfig
from utils.models import DNN, DNN_GRF
import torch

def use_dask(condor,njobs,port):
    from dask.distributed import Client
    from lpc_dask.lpc_dask import HTCondorCluster
    import socket

    # make list of local package directories (absolute paths) that should be sent to jobs
    initpylist = [os.path.abspath(os.path.dirname(x)) for x in glob('*/__init__.py')]
    initpylist.append("patch.sh")
    job_extra = {'transfer_input_files': ','.join(initpylist)}

    extra = ['--worker-port 10002:10100']

    hostname = socket.gethostname()

    if condor:
        cluster = HTCondorCluster(
            scheduler_options = {'host': f'{hostname}:10000', 'dashboard_address': ':{}'.format(port)},
            cores=1,
            memory="2GB",
            disk="2GB",
            python='python',
            nanny=False,
            extra=extra,
            job_extra=job_extra,
        )

        cluster.scale(jobs=njobs)

        client = Client(cluster,
            timeout=100
        )
    else:
        client = Client()

    exe_args = {
        'client': client,
        'savemetrics': True,
        'schema': None,
        'nano': False,
        'align_clusters': True
    }

    return exe_args

def main():
    # start run time clock
    tstart = time.time()

    # get options from command line
    parser = OptionParser()
    parser.add_option('-d', '--dataset',   help='dataset',           dest='dataset')
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
    options, args = parser.parse_args()

    # set output root file
    sample = options.dataset
    outfile = "MyAnalysis_%s_%d.root" % (sample, options.startFile) if options.condor or options.dask else "test.root"

    # getting dictionary of files from a sample collection e.g. "2016_QCD, 2016_WJets, 2016_TTJets, 2016_ZJets"
    fileset = s.getFileset(sample, True, options.startFile, options.nFiles)

    # get processor args
    exe_args = {'workers': options.workers, 'flatten': False}
    if options.dask:
        exe_args = use_dask(options.condor,options.workers,options.port)
        if options.quiet: exe_args['status'] = False

        client = exe_args['client']
        while len(client.ncores()) < options.mincores:
            print('Waiting for more cores to spin up, currently there are {0} available...'.format(len(client.ncores())))
            print('Dask client info ->', client)
            time.sleep(10)

    sf = s.sfGetter(sample)
    print("scaleFactor = {}".format(sf))
    # open saved neural network
    device = torch.device('cpu')
    modelLocation = "."
    varSet = ['njets', 'njetsAK8', 'nb', 'dPhij1rdPhij2AK8', 'dPhiMinjMETAK8', 'dEtaj12AK8', 'dRJ12AK8', 'jGirthAK8', 'jTau1AK8', 'jTau2AK8', 'jTau3AK8', 'jTau21AK8', 'jTau32AK8', 'jSoftDropMassAK8', 'jAxisminorAK8', 'jAxismajorAK8', 'jPtDAK8', 'jecfN2b1AK8', 'jecfN3b1AK8', 'jEleEFractAK8', 'jMuEFractAK8', 'jNeuHadEFractAK8', 'jPhoEFractAK8', 'jPhoMultAK8', 'jNeuMultAK8', 'jNeuHadMultAK8', 'jMuMultAK8', 'jEleMultAK8', 'jChHadMultAK8', 'jChMultAK8', 'jNeuEmEFractAK8', 'jHfHadEFractAK8', 'jHfEMEFractAK8', 'jChEMEFractAK8', 'jMultAK8', 'jecfN3b2AK8', 'jecfN2b2AK8', 'jPhiAK8', 'jEtaAK8']
    hyper = MagiConfig(batchSize=2000, dropout=0.3, epochs=10, lambdaDC=0.0, lambdaGR=1.0, lambdaReg=0.0, lambdaTag=1.0, learning_rate=0.001, n_pTBins=35, num_of_layers_features=2, num_of_layers_pT=5, num_of_layers_tag=2, num_of_nodes=40, pTBins=[50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2500, 3000, 3500, 4000, 4500], rseed=30)
    model = DNN_GRF(n_var=len(varSet),  n_layers_features=hyper.num_of_layers_features, n_layers_tag=hyper.num_of_layers_tag, n_layers_pT=hyper.num_of_layers_pT, n_nodes=hyper.num_of_nodes, n_outputs=2, n_pTBins=hyper.n_pTBins, drop_out_p=hyper.dropout).to(device=device)
    print("Loading model from file {}/net.pth".format(modelLocation))
    model.load_state_dict(torch.load("{}/net.pth".format(modelLocation),map_location=device))
    model.eval()
    model.to('cpu')
    normMeanStd = np.load("{}/normMeanStd.npz".format(modelLocation))
    normMean = normMeanStd["normMean"]
    normStd = normMeanStd["normStd"]

    # run processor
    output = processor.run_uproot_job(
        fileset,
        treename='TreeMaker2/PreSelection',
        # processor_instance=MainProcessor(sf,model,varSet,normMean,normStd),
        processor_instance=MainProcessor(sf,model,varSet,normMean,normStd),
        executor=processor.dask_executor if options.dask else processor.futures_executor,
        executor_args=exe_args,
        chunksize=options.chunksize,
        maxchunks=options.maxchunks,
    )

    # export the histograms to root files
    ## the loop makes sure we are only saving the histograms that are filled
    fout = uproot.recreate(outfile)
    if isinstance(output,tuple): output = output[0]
    for key,H in output.items():
        if type(H) is hist.Hist and H._sumw2 is not None:
            fout[key] = hist.export1d(H)
    fout.close()

    # print run time in seconds
    dt = time.time() - tstart
    print("run time: %.2f [sec]" % (dt))

if __name__ == "__main__":
    main()
