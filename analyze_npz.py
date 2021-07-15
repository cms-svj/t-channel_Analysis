#!/usr/bin/env python

from coffea import hist, processor
from processors.npzProcessor import MainProcessor
import uproot
import sys,os
from utils import samples as s
import time
from optparse import OptionParser
from glob import glob
import numpy as np

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

def dataStream(output,fname):
    np.savez(
    fname,
    evtw=output["evtw"].value,
    jw=output["jw"].value,
    fjw=output["fjw"].value,
    njets=output["njets"].value,
    njetsAK8=output["njetsAK8"].value,
    nb=output["nb"].value,
    nl=output["nl"].value,
    ht=output["ht"].value,
    st=output["st"].value,
    met=output["met"].value,
    madHT=output["madHT"].value,
    jPt=output["jPt"].value,
    jEta=output["jEta"].value,
    jPhi=output["jPhi"].value,
    jAxismajor=output["jAxismajor"].value,
    jAxisminor=output["jAxisminor"].value,
    jPtD=output["jPtD"].value,
    dPhiMinjMET=output["dPhiMinjMET"].value,
    dEtaj12=output["dEtaj12"].value,
    dRJ12=output["dRJ12"].value,
    jPtAK8=output["jPtAK8"].value,
    jEtaAK8=output["jEtaAK8"].value,
    jPhiAK8=output["jPhiAK8"].value,
    jAxismajorAK8=output["jAxismajorAK8"].value,
    jAxisminorAK8=output["jAxisminorAK8"].value,
    jGirthAK8=output["jGirthAK8"].value,
    jPtDAK8=output["jPtDAK8"].value,
    jTau1AK8=output["jTau1AK8"].value,
    jTau2AK8=output["jTau2AK8"].value,
    jTau3AK8=output["jTau3AK8"].value,
    jTau21AK8=output["jTau21AK8"].value,
    jTau32AK8=output["jTau32AK8"].value,
    jSoftDropMassAK8=output["jSoftDropMassAK8"].value,
    dPhiMinjMETAK8=output["dPhiMinjMETAK8"].value,
    dEtaj12AK8=output["dEtaj12AK8"].value,
    dRJ12AK8=output["dRJ12AK8"].value,
    mT=output["mT"].value,
    METrHT_pt30=output["METrHT_pt30"].value,
    METrST_pt30=output["METrST_pt30"].value,
    j1Pt=output["j1Pt"].value,
    j1Eta=output["j1Eta"].value,
    j1Phi=output["j1Phi"].value,
    j1Axismajor=output["j1Axismajor"].value,
    j1Axisminor=output["j1Axisminor"].value,
    j1PtD=output["j1PtD"].value,
    dPhij1MET=output["dPhij1MET"].value,
    j2Pt=output["j2Pt"].value,
    j2Eta=output["j2Eta"].value,
    j2Phi=output["j2Phi"].value,
    j2Axismajor=output["j2Axismajor"].value,
    j2Axisminor=output["j2Axisminor"].value,
    j2PtD=output["j2PtD"].value,
    dPhij2MET=output["dPhij2MET"].value,
    dPhij1rdPhij2=output["dPhij1rdPhij2"].value,
    j1PtAK8=output["j1PtAK8"].value,
    j1EtaAK8=output["j1EtaAK8"].value,
    j1PhiAK8=output["j1PhiAK8"].value,
    j1AxismajorAK8=output["j1AxismajorAK8"].value,
    j1AxisminorAK8=output["j1AxisminorAK8"].value,
    j1GirthAK8=output["j1GirthAK8"].value,
    j1PtDAK8=output["j1PtDAK8"].value,
    j1Tau1AK8=output["j1Tau1AK8"].value,
    j1Tau2AK8=output["j1Tau2AK8"].value,
    j1Tau3AK8=output["j1Tau3AK8"].value,
    j1Tau21AK8=output["j1Tau21AK8"].value,
    j1Tau32AK8=output["j1Tau32AK8"].value,
    j1SoftDropMassAK8=output["j1SoftDropMassAK8"].value,
    dPhij1METAK8=output["dPhij1METAK8"].value,
    j2PtAK8=output["j2PtAK8"].value,
    j2EtaAK8=output["j2EtaAK8"].value,
    j2PhiAK8=output["j2PhiAK8"].value,
    j2AxismajorAK8=output["j2AxismajorAK8"].value,
    j2AxisminorAK8=output["j2AxisminorAK8"].value,
    j2GirthAK8=output["j2GirthAK8"].value,
    j2PtDAK8=output["j2PtDAK8"].value,
    j2Tau1AK8=output["j2Tau1AK8"].value,
    j2Tau2AK8=output["j2Tau2AK8"].value,
    j2Tau3AK8=output["j2Tau3AK8"].value,
    j2Tau21AK8=output["j2Tau21AK8"].value,
    j2Tau32AK8=output["j2Tau32AK8"].value,
    j2SoftDropMassAK8=output["j2SoftDropMassAK8"].value,
    dPhij2METAK8=output["dPhij2METAK8"].value,
    dPhij1rdPhij2AK8=output["dPhij1rdPhij2AK8"].value
    )

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
    dataStream(output, sample)

    # print run time in seconds
    dt = time.time() - tstart
    print("run time: %.2f [sec]" % (dt))

if __name__ == "__main__":
    main()
