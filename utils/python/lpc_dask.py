import dask.distributed
from coffea import processor
import uproot
import os 

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

def run_processor(fileset,sample,sf,MainExecutor,MainProcessor,options):
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

def runProcessWithErrorHandling(fileset,sample,sf,MainExecutor,MainProcessor,options):
    import dask.distributed
    from distributed import CancelledError
    from distributed.scheduler import KilledWorker
    from distributed.comm.core import CommClosedError

    client = use_dask(options.condor,options.jobs,options.port)
    print("Waiting for at least one worker...")
    client.wait_for_workers(1)

    MainExecutor = processor.dask_executor(client=client)
    # error handling
    try:
      run_processor(fileset,sample,sf,MainExecutor,MainProcessor,options)
    except (CancelledError, CommClosedError) as err:
      print(
      err, "Cluster communication error, restart workers and trying again.")
      restart_client(client)
      run_processor(fileset,sample,sf,MainExecutor,MainProcessor,options)
    except KilledWorker as err:
      print(
      err, """Some job consistently got worker nodes killed. Assuming dirty
      worker state, restarting client and and trying just once. If issue
      persists this should be a true error that needs to be looked at.""")
      restart_client(client)
      run_processor(fileset,sample,sf,MainExecutor,MainProcessor,options)
    except TimeoutError as err:
      print(
      err, """Connection time out error, this is is typically a worker node
      start up issue, retrying indefinitely. Make sure all grid
      requirements are met if issue persists.""")
      restart_client(client)
      run_processor(fileset,sample,sf,MainExecutor,MainProcessor,options)
    except OSError as err:
      err_str = str(err)
      if 'XRootD' in err_str and 'Operation expired' in err_str:
          print(
          err, """Time out error during xrootd, the offending file will be
          added to a block-listed if the same error for that file occurs more
          than 3 times.""")
          bad_file = err_str.split()[-1]  # Getting the bad file.
          print("bad_file",bad_file)
          run_processor(fileset,sample,sf,MainExecutor,MainProcessor,options)

      if 'XRootD' in err_str:
      # Note, missing files would be FileNotFound error, not an OSError.
          print(
          err, """Error during xrootd, the job will be automatically
          restarted, but make sure all xrootd requirements are met if issue
          persists (grid certificates validity session).""")
          restart_client(client)
          run_processor(fileset,sample,sf,MainExecutor,MainProcessor,options)
      else:
          raise err
    except Exception as err:
      print("Unrecognized error raised, this is actually a problem.")
      raise err