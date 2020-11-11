import sys, os
from os import system, environ
sys.path = [environ["TCHANNEL_BASE"],] + sys.path
from utils.python.samples import Samples as s
import optparse

def removeCopies(x):
  return list(dict.fromkeys(x))

def makeExeAndFriendsTarball(filestoTransfer, fname, path):
    system("mkdir -p %s" % fname)
    for fn in removeCopies(filestoTransfer):
        system("cd %s; ln -s %s" % (fname, fn))

    tarallinputs = "tar czf %s/%s.tar.gz %s --dereference"% (path, fname, fname)
    system(tarallinputs)
    system("rm -r %s" % fname)

def getDatasets(datasets):
    if datasets:
        return datasets.split(',')
    else:
        print("No dataset specified")
        exit(0)

def main():
    # parse command line arguments
    parser = optparse.OptionParser("usage: %prog [options]\n")
    parser.add_option ('-n',             dest='numfile',  type='int',                         default = 10,            help="number of files per job")
    parser.add_option ('-d',             dest='datasets', type='string',                      default = '',            help="List of datasets, comma separated")
    parser.add_option ('-c',             dest='noSubmit',                action='store_true', default = False,         help="Do not submit jobs.  Only create condor_submit.txt.")
    parser.add_option ('-w','--workers', dest='workers',  type='int',                         default = 2,             help='Number of workers to use for multi-worker executors (e.g. futures or condor)')
    parser.add_option ('--output',       dest='outPath',  type='string',                      default = '.',           help="Name of directory where output of each condor job goes")
    options, args = parser.parse_args()

    # prepare the list of hardcoded files to transfer
    filestoTransfer = []

    # add top of jdl file
    fileParts = []
    fileParts.append("Universe   = vanilla\n")
    fileParts.append("Executable = run_Analyzer_condor.sh\n")
    fileParts.append("Transfer_Input_Files = %s/%s.tar.gz, %s/exestuff.tar.gz\n" % (options.outPath,"tchannel",options.outPath))
    fileParts.append("Should_Transfer_Files = YES\n")
    fileParts.append("WhenToTransferOutput = ON_EXIT\n")
    fileParts.append("x509userproxy = $ENV(X509_USER_PROXY)\n\n")

    # loop over all sample collections in the dataset
    datasets = getDatasets(options.datasets)
    nFilesPerJob = options.numfile
    numberOfJobs = 0
    print("-"*50)
    for sc in datasets:
        print(sc)

        # create the directory
        if not os.path.isdir("%s/output-files/%s" % (options.outPath, sc)):
            os.makedirs("%s/output-files/%s" % (options.outPath, sc))

        # loop over all samples in the sample collection
        samples = s().getFileset(sc, False)
        for n, rFiles in samples.items():
            count = len(rFiles)
            print("    %-40s %d" % (n, count))

            # loop over the root files that will be in each job
            for startFileNum in range(0, count, nFilesPerJob):
                numberOfJobs+=1
                outputDir = "%s/output-files/%s" % (options.outPath, sc)

                # list the output files that will be transfered to output directory
                outputFiles = [
                   "MyAnalysis_%s_%s.root" % (n, startFileNum),
                ]
                transfer = "transfer_output_remaps = \""
                for f_ in outputFiles:
                    transfer += "%s = %s/%s" % (f_, outputDir, f_)
                    if f_ != outputFiles[-1]:
                        transfer += "; "
                transfer += "\"\n"

                # add each job to the jdl file
                fileParts.append(transfer)
                fileParts.append("Arguments = %s %i %i %i\n"%(n, nFilesPerJob, startFileNum, options.workers))
                fileParts.append("Output = %s/log-files/MyAnalysis_%s_%i.stdout\n"%(options.outPath, n, startFileNum))
                fileParts.append("Error = %s/log-files/MyAnalysis_%s_%i.stderr\n"%(options.outPath, n, startFileNum))
                fileParts.append("Log = %s/log-files/MyAnalysis_%s_%i.log\n"%(options.outPath, n, startFileNum))
                fileParts.append("Queue\n\n")
        print("-"*50)

    # write out the jdl file
    fout = open("condor_submit.jdl", "w")
    fout.write(''.join(fileParts))
    fout.close()

    # print number jobs to run
    print("Number of Jobs:", numberOfJobs)

    # only runs when you submit
    if not options.noSubmit:
        # tar up working area to send with each job
        makeExeAndFriendsTarball(filestoTransfer, "exestuff", options.outPath)
        #system("tar --exclude-caches-all --exclude-vcs -zcf %s/tchannel.tar.gz -C ${TCHANNEL_BASE}/.. tchannel --exclude=src --exclude=tmp" % options.outPath)
        system("tar czf %s/tchannel.tar.gz -C ${TCHANNEL_BASE} . --exclude=EventLoopFramework --exclude=condor --exclude=notebooks --exclude=root --exclude=.git --exclude=coffeaenv.tar.gz" % options.outPath)

        # submit the jobs to condor
        system('mkdir -p %s/log-files' % options.outPath)
        system("echo 'condor_submit condor_submit.jdl'")
        system('condor_submit condor_submit.jdl')

if __name__ == "__main__":
    main()
