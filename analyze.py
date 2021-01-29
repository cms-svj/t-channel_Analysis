#!/usr/bin/env python

from coffea import hist, processor
from processors.mainProcessor import MainProcessor
import uproot
import sys,os
from utils import samples as s
import time
from optparse import OptionParser
from glob import glob

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

    # run processor
    output = processor.run_uproot_job(
        fileset,
        treename='TreeMaker2/PreSelection',
        processor_instance=MainProcessor(),
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
