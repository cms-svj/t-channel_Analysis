import dask.distributed
from coffea import processor
import uproot
import os 
from utils.coffea.n_tree_maker_schema import NTreeMakerSchema
from utils.shortNameDict import shortNameDict

def use_dask(condor,jobs,port):
    from dask.distributed import Client, LocalCluster
    from lpcjobqueue import LPCCondorCluster
    import socket

    hostname = socket.gethostname()

    if condor:
        cluster = LPCCondorCluster(
            scheduler_options = {'host': f'{hostname}:10000', 'dashboard_address': ':{}'.format(port)},
            cores=1,
            memory="4GB",
            disk="2GB",
            transfer_input_files=[f'{os.getenv("TCHANNEL_BASE")}/input',f'{os.getenv("TCHANNEL_BASE")}/utils',f'{os.getenv("TCHANNEL_BASE")}/processors',f'{os.getenv("TCHANNEL_BASE")}/coffeaenv/lib64/python3.8/site-packages/magiconfig.py'], # for some reason magiconfig cannot be imported properly when running dask condor
            log_directory=None,
            death_timeout=180
        )
        cluster.scale(jobs=jobs)
        # if jobs < 20:
        #     cluster.scale(jobs=jobs)
        # else:
        #     cluster.adapt(minimum=100,maximum=jobs)
        client = Client(cluster,
            timeout=100
        )
    else:
        cluster = LocalCluster(memory_limit="4GiB",scheduler_port=0)
        cluster.scale(4)
        client = Client(cluster)
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

def out_file_name_creator(outHistF,dataset,nFiles,startFile,condor,dask,hemPeriod):
    details = ""
    sampleList = dataset
    nFilesList = nFiles
    startFileList = startFile
    year = sampleList[0][:sampleList[0].find("_")]
    for i in range(len(sampleList)):
        shortName = shortNameDict[sampleList[i].replace(f"{year}_","")]
        details += f"{shortName}_N{nFilesList[i]}_M{startFileList[i]}_"
    outfile = f"{outHistF}/{year}_{details}{hemPeriod}.root" 
    return outfile

def run_processor(fileset,sample,MainExecutor,MainProcessor,args,exe_args,evtTaggerDict={},trainingKind="",trainFileProduction=False):
    ###########################################################################################################
    # run processor
    ###########################################################################################################
    if args.skimSource:
        treename='Events'
    else:
        treename='TreeMaker2/PreSelection'
    output = processor.run_uproot_job(
        fileset,
        treename=treename,
        processor_instance=MainProcessor(jNVar=args.jNVar,hemStudy=args.hemStudy,hemPeriod=args.hemPeriod,evtTaggerDict=evtTaggerDict,tcut=args.tcut,sFactor=args.sFactor,skimSource=args.skimSource,skimCut=args.skimCut,runJetTag=args.runJetTag,runEvtClass=args.runEvtClass),
        executor=MainExecutor,
        executor_args=exe_args,
        chunksize=args.chunksize,
        maxchunks=args.maxchunks,
    )
    outHistF = args.outHistF
    os.makedirs(outHistF,exist_ok=True)

    if trainFileProduction:
        # saving input to root files.
        # When inputs are too large, saving the inputs could return an error. That's why we need to save them in batches.
        maxNumOfJets = 100000
        jVarName = "njetsAK8"
        if trainingKind == "PN":
            jVarName = "pT"
        if isinstance(output,tuple): output = output[0]
        output = dict(sorted(output.items()))
        numOfJets = len(output[jVarName].value)
        print("numOfJets",numOfJets)
        sampleName = sample.replace("_TuneCP5_13TeV-madgraphMLM-pythia8","")
        partNumber = 0
        for i in range(0,numOfJets,maxNumOfJets):
            outputNPZ = {}
            fileName = "{}/tree_{}_M{}_{}.root".format(outHistF,sampleName,args.startFile,partNumber)
            partNumber += 1
            file = uproot.recreate(fileName)
            for v in output.keys():
                outShape = (output[v].value).shape
                if (len(outShape) > 1) and (trainingKind == "NN"):
                    for j in range(outShape[1]):
                        try:
                            outputNPZ[f"{v}_{j}"] = output[v].value[:,j][i:i+maxNumOfJets]
                        except:
                            outputNPZ[f"{v}_{j}"] = output[v].value[:,j][i:]
                else:
                    try:
                        outputNPZ[v] = output[v].value[i:i+maxNumOfJets]
                    except:
                        outputNPZ[v] = output[v].value[i:]
            file["tree"] = outputNPZ
            file["tree"].show()
    else:
    ###########################################################################################################
    # export the histograms to root files
    ## the loop makes sure we are only saving the histograms that are filled
    ###########################################################################################################
        outfile = out_file_name_creator(outHistF,args.dataset,args.nFiles,args.startFile,args.condor,args.dask,args.hemPeriod)
        fout = uproot.recreate(outfile)
        if isinstance(output,tuple): output = output[0]
        output = dict(sorted(output.items()))
        for key,H in output.items():
            if len(H.axes) == 0:
                print("Somethings wrong with this histogram \"{}\" skipping the write out".format(H.label))
                continue
            fout[key] = H
        fout.close()

def runProcessWithErrorHandling(fileset,sample,MainExecutor,MainProcessor,args,evtTaggerDict={},trainingKind="",trainFileProduction=False):
    import dask.distributed
    from distributed import CancelledError
    from distributed.scheduler import KilledWorker
    from distributed.comm.core import CommClosedError
    client = use_dask(args.condor,args.jobs,args.port)
    print("Waiting for at least one worker...")
    client.wait_for_workers(1)
    schema = processor.TreeMakerSchema
    if args.skimSource:
        schema = NTreeMakerSchema
    exe_args = {
        'client': client,
        'savemetrics': True,
        'schema': schema,
        #'nano': False,
        'align_clusters': True
    }
    MainExecutor = processor.dask_executor
    # error handling
    try:
      run_processor(fileset,sample,MainExecutor,MainProcessor,args,exe_args,evtTaggerDict,trainingKind,trainFileProduction)
    except (CancelledError, CommClosedError) as err:
      print(
      err, "Cluster communication error, restart workers and trying again.")
      restart_client(client)
      run_processor(fileset,sample,MainExecutor,MainProcessor,args,exe_args,evtTaggerDict,trainingKind,trainFileProduction)
    except KilledWorker as err:
      print(
      err, """Some job consistently got worker nodes killed. Assuming dirty
      worker state, restarting client and and trying just once. If issue
      persists this should be a true error that needs to be looked at.""")
      restart_client(client)
      run_processor(fileset,sample,MainExecutor,MainProcessor,args,exe_args,evtTaggerDict,trainingKind,trainFileProduction)
    except TimeoutError as err:
      print(
      err, """Connection time out error, this is is typically a worker node
      start up issue, retrying indefinitely. Make sure all grid
      requirements are met if issue persists.""")
      restart_client(client)
      run_processor(fileset,sample,MainExecutor,MainProcessor,args,exe_args,evtTaggerDict,trainingKind,trainFileProduction)
    except OSError as err:
      err_str = str(err)
      if 'XRootD' in err_str and 'Operation expired' in err_str:
          print(
          err, """Time out error during xrootd, the offending file will be
          added to a block-listed if the same error for that file occurs more
          than 3 times.""")
          bad_file = err_str.split()[-1]  # Getting the bad file.
          print("bad_file",bad_file)
          run_processor(fileset,sample,MainExecutor,MainProcessor,args,exe_args,evtTaggerDict,trainingKind,trainFileProduction)

      if 'XRootD' in err_str:
      # Note, missing files would be FileNotFound error, not an OSError.
          print(
          err, """Error during xrootd, the job will be automatically
          restarted, but make sure all xrootd requirements are met if issue
          persists (grid certificates validity session).""")
          restart_client(client)
          run_processor(fileset,sample,MainExecutor,MainProcessor,args,exe_args,evtTaggerDict,trainingKind,trainFileProduction)
      else:
          raise err
    except Exception as err:
      print("Unrecognized error raised, this is actually a problem.")
      raise err