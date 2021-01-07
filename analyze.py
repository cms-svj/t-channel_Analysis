#!/usr/bin/env python
from coffea import hist, processor
from processors.mainProcessor import MainProcessor
import uproot
import sys
from utils import samples as s
import time
from optparse import OptionParser

def main():
    # start run time clock
    tstart = time.time()

    # get options from command line
    parser = OptionParser()
    parser.add_option('-d', '--dataset',   help='dataset',           dest='dataset')
    parser.add_option('-N', '--nFiles',    help='nFiles',            dest='nFiles',    type=int, default=-1)
    parser.add_option('-M', '--startFile', help='startFile',         dest='startFile', type=int, default=0)
    parser.add_option(      '--condor',    help='running on condor', dest='condor',              default=False, action='store_true')
    parser.add_option('-w', '--workers',   help='Number of workers to use for multi-worker executors (e.g. futures or condor)', dest='workers', type=int, default=8)
    parser.add_option('-s', '--chunksize', help='Chunk size',        dest='chunksize', type=int, default=10000)
    options, args = parser.parse_args()

    # set output root file
    sample = options.dataset
    outfile = "MyAnalysis_%s_%d.root" % (sample, options.startFile) if options.condor else "test.root"

    # getting dictionary of files from a sample collection e.g. "2016_QCD, 2016_WJets, 2016_TTJets, 2016_ZJets"
    fileset = s.getFileset(sample, True, options.startFile, options.nFiles)
    # run processor
    output = processor.run_uproot_job(
        fileset,
        treename='TreeMaker2/PreSelection',
        processor_instance=MainProcessor(),
        executor=processor.futures_executor,
        executor_args={'workers': options.workers, 'flatten': False},
        chunksize=options.chunksize,
    )

    # export the histograms to root files
    ## the loop makes sure we are only saving the histograms that are filled
    fout = uproot.recreate(outfile)
    for key,H in output.items():
        if type(H) is hist.Hist and H._sumw2 is not None:
            fout[key] = hist.export1d(H)
    fout.close()

    # print run time in seconds
    dt = time.time() - tstart
    print("run time: %.2f [sec]" % (dt))

if __name__ == "__main__":
    main()
