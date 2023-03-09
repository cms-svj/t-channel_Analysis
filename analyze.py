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
import dask.distributed
from distributed import CancelledError
from distributed.scheduler import KilledWorker
from distributed.comm.core import CommClosedError

def use_dask(condor,jobs,port):
    from dask.distributed import Client
    from lpcjobqueue import LPCCondorCluster
    import socket
    import dask

    hostname = socket.gethostname()

    if condor:
        cluster = LPCCondorCluster(
            scheduler_options = {'host': f'{hostname}:10000', 'dashboard_address': ':{}'.format(port)},
            cores=1,
            memory="4GB",
            disk="2GB",
            transfer_input_files=[f'{os.getenv("TCHANNEL_BASE")}/utils',f'{os.getenv("TCHANNEL_BASE")}/processors'],
            log_directory=None,
            death_timeout=180
        )
        if jobs < 100:
            cluster.scale(jobs=jobs)
        else:    
            cluster.scale(jobs=100)
        # if jobs < 20:
        #     cluster.scale(jobs=jobs)
        # else:
        #     cluster.adapt(minimum=100,maximum=jobs)
        client = Client(cluster,
            timeout=100
        )
    else:
        client = Client()
    return client

def restart_client(client):
    try:
        client.restart(timeout=180, wait_for_workers=False)
        # Wait for workers in dedicated function
        n_workers = len(client.scheduler_info()['workers'])
        client.wait_for_workers(n_workers // 2)
        time.sleep(2)  # Short wait for 5 seconds.
    except:
        time.sleep(10)
        pass

def run_processor(fileset,sample,sf,MainExecutor,options):
    ###########################################################################################################
    # run processor
    ###########################################################################################################
    run = processor.Runner(
        executor=MainExecutor,
        schema=processor.TreeMakerSchema,
        chunksize=options.chunksize,
        maxchunks=options.maxchunks,
        align_clusters = True,
        savemetrics = True
    )

    output = run(
        fileset,
        "TreeMaker2/PreSelection",
        processor_instance=MainProcessor(sample,sf,options.jNVar),
    )

    ###########################################################################################################
    # export the histograms to root files
    ## the loop makes sure we are only saving the histograms that are filled
    ###########################################################################################################
    outfile = "MyAnalysis_%s_%d.root" % (sample, options.startFile) if options.condor or options.dask else "test.root"
    fout = uproot.recreate(outfile)
    if isinstance(output,tuple): output = output[0]
    output = dict(sorted(output.items()))
    for key,H in output.items():
        if len(H.axes) == 0:
            print("Somethings wrong with this histogram \"{}\" skipping the write out".format(H.label))
            continue
        fout[key] = H
    fout.close()

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
        import dask.distributed
        from distributed import CancelledError
        from distributed.scheduler import KilledWorker
        from distributed.comm.core import CommClosedError
        
        client = use_dask(options.condor,options.jobs,options.port)
        print("Waiting for at least one worker...")
        client.wait_for_workers(1)
        # while len(client.ncores()) < options.mincores:
        #     print('Waiting for more cores to spin up, currently there are {0} available...'.format(len(client.ncores())))
        #     print('Dask client info ->', client)
        #     time.sleep(10)
        MainExecutor = processor.dask_executor(client=client)
        # error handling
        try:
            run_processor(fileset,sample,sf,MainExecutor,options)
        except (CancelledError, CommClosedError) as err:
            print(
            err, "Cluster communication error, restart workers and trying again.")
            restart_client(client)
            run_processor(fileset,sample,sf,MainExecutor,options)
        except KilledWorker as err:
            print(
            err, """Some job consistently got worker nodes killed. Assuming dirty
            worker state, restarting client and and trying just once. If issue
            persists this should be a true error that needs to be looked at.""")
            restart_client(client)
            run_processor(fileset,sample,sf,MainExecutor,options)
        except TimeoutError as err:
            print(
            err, """Connection time out error, this is is typically a worker node
            start up issue, retrying indefinitely. Make sure all grid
            requirements are met if issue persists.""")
            restart_client(client)
            run_processor(fileset,sample,sf,MainExecutor,options)
        except OSError as err:
            err_str = str(err)
            if 'XRootD' in err_str and 'Operation expired' in err_str:
                print(
                err, """Time out error during xrootd, the offending file will be
                added to a block-listed if the same error for that file occurs more
                than 3 times.""")
                bad_file = err_str.split()[-1]  # Getting the bad file.
                print("bad_file",bad_file)
                run_processor(fileset,sample,sf,MainExecutor,options)

            if 'XRootD' in err_str:
            # Note, missing files would be FileNotFound error, not an OSError.
                print(
                err, """Error during xrootd, the job will be automatically
                restarted, but make sure all xrootd requirements are met if issue
                persists (grid certificates validity session).""")
                restart_client(client)
                run_processor(fileset,sample,sf,MainExecutor,options)
            else:
                raise err
        except Exception as err:
            print("Unrecognized error raised, this is actually a problem.")
            raise err
    else:
        run_processor(fileset,sample,sf,MainExecutor,options)

    ###########################################################################################################
    # print run time in seconds
    ###########################################################################################################
    dt = time.time() - tstart
    print("run time: %.2f [sec]" % (dt))

if __name__ == "__main__":
    main()
