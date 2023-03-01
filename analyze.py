#!/usr/bin/env python

from coffea import processor
import hist as h
from processors.mainProcessor import MainProcessor
import uproot
import sys,os
from utils import samples as s
import time
from optparse import OptionParser
from glob import glob
import numpy as np

def use_dask(condor,jobs,port):
    from dask.distributed import Client
    from lpcjobqueue import LPCCondorCluster
    import socket
    import dask

    if condor:
        proxyfile = '{0}/x509up_u{1}'.format(os.environ['HOME'], os.getuid())
        dask.config.set({
            'distributed.worker.profile.interval': '1d',
            'distributed.worker.profile.cycle': '2d',
            'distributed.worker.profile.low-level': False,
        })
        cluster = LPCCondorCluster(
            cores=1,
            memory="2GB",
            disk="2GB",
            transfer_input_files=[f'{os.getenv("TCHANNEL_BASE")}/utils',f'{os.getenv("TCHANNEL_BASE")}/processors'],
            ship_env=True,
            log_directory=None,
            job_extra_directives={'x509userproxy': proxyfile},
            death_timeout=180
        )
        if jobs < 20:
            cluster.scale(jobs=jobs)
        else:
            cluster.adapt(minimum=50,maximum=150)

        client = Client(cluster,
            timeout=100
        )
    else:
        client = Client()

    exe_args = {
        'client': client,
        'savemetrics': True,
        'schema': processor.TreeMakerSchema,
        #'nano': False,
        'align_clusters': True
    }

    return exe_args

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
    parser.add_option('-s', '--chunksize', help='Chunk size',        dest='chunksize', type=int, default=10000)
    parser.add_option('-m', '--maxchunks', help='Max number of chunks (for testing)',        dest='maxchunks', type=int, default=None)
    options, args = parser.parse_args()

    ###########################################################################################################
    # set output root file
    ###########################################################################################################
    sample = options.dataset
    outfile = "MyAnalysis_%s_%d.root" % (sample, options.startFile) if options.condor or options.dask else "test.root"

    # getting dictionary of files from a sample collection e.g. "2016_QCD, 2016_WJets, 2016_TTJets, 2016_ZJets"
    fileset = s.getFileset(sample, True, options.startFile, options.nFiles)
    sf = s.sfGetter(sample,False)
    print("scaleFactor = {}".format(sf))

    ###########################################################################################################
    # get executor/processor args
    ###########################################################################################################
    exe_args = {
        'workers': options.workers, 
        'schema': processor.TreeMakerSchema
    }
    MainExecutor = processor.futures_executor

    if options.dask:
        exe_args = use_dask(options.condor,options.workers,options.port)
        if options.quiet: exe_args['status'] = False

        client = exe_args['client']
        while len(client.ncores()) < options.mincores:
            print('Waiting for more cores to spin up, currently there are {0} available...'.format(len(client.ncores())))
            print('Dask client info ->', client)
            time.sleep(10)

        MainExecutor = processor.dask_executor

    ###########################################################################################################
    # run processor
    ###########################################################################################################
    output = processor.run_uproot_job(
        fileset,
        treename='TreeMaker2/PreSelection',
        processor_instance=MainProcessor(sample,sf,options.jNVar),
        executor=MainExecutor,
        executor_args=exe_args,
        chunksize=options.chunksize,
        maxchunks=options.maxchunks,
    )

    ###########################################################################################################
    # export the histograms to root files
    ## the loop makes sure we are only saving the histograms that are filled
    ###########################################################################################################
    fout = uproot.recreate(outfile)
    if isinstance(output,tuple): output = output[0]
    output = dict(sorted(output.items()))
    for key,H in output.items():
        if len(H.axes) == 0:
            print("Somethings wrong with this histogram \"{}\" skipping the write out".format(H.label))
            continue
        fout[key] = H
    fout.close()

    ###########################################################################################################
    # print run time in seconds
    ###########################################################################################################
    dt = time.time() - tstart
    print("run time: %.2f [sec]" % (dt))

if __name__ == "__main__":
    main()
