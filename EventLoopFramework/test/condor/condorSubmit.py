#!/cvmfs/cms.cern.ch/slc6_amd64_gcc630/cms/cmssw/CMSSW_9_3_3/external/slc6_amd64_gcc630/bin/python
####!${SRT_CMSSW_RELEASE_BASE_SCRAMRTDEL}/external/${SCRAM_ARCH}/bin/python

import sys, os
from os import system, environ

from samples import SampleCollection
import optparse 
import subprocess

def removeCopies(x):
  return list(dict.fromkeys(x))

def makeExeAndFriendsTarball(filestoTransfer, fname, path):
    system("mkdir -p %s" % fname)
    for fn in removeCopies(filestoTransfer):
        system("cd %s; ln -s %s" % (fname, fn))
        
    tarallinputs = "tar czf %s/%s.tar.gz %s --dereference"% (path, fname, fname)
    system(tarallinputs)
    system("rm -r %s" % fname)

def getTopTaggerTrainingFile(topTaggerFile):
    name = ""
    with file(topTaggerFile) as meowttcfgFile:
        for line in meowttcfgFile:
            if "modelFile" in line:
                name = line.split("=")[1].strip().strip("\"")
                break
    return name

def main():
    repo = "Analyzer/Analyzer"    
    # Parse command line arguments
    parser = optparse.OptionParser("usage: %prog [options]\n")    
    parser.add_option ('-n',        dest='numfile',  type='int',                         default = 10,            help="number of files per job")
    parser.add_option ('-d',        dest='datasets', type='string',                      default = '',            help="List of datasets, comma separated")
    parser.add_option ('-l',        dest='dataCollections',         action='store_true', default = False,         help="List all datacollections")
    parser.add_option ('-L',        dest='dataCollectionslong',     action='store_true', default = False,         help="List all datacollections and sub collections")
    parser.add_option ('-c',        dest='noSubmit',                action='store_true', default = False,         help="Do not submit jobs.  Only create condor_submit.txt.")
    parser.add_option ('--output',  dest='outPath',  type='string',                      default = '.',           help="Name of directory where output of each condor job goes")
    parser.add_option ('--analyze', dest='analyze',                                      default = 'Analyze1Lep', help="AnalyzeBackground, AnalyzeEventSelection, Analyze0Lep, Analyze1Lep, MakeNJetDists")    
    options, args = parser.parse_args()
    
    # Prepare the list of files to transfer
    mvaFileName2016 = getTopTaggerTrainingFile(environ["CMSSW_BASE"] + "/src/%s/test/TopTaggerCfg_2016.cfg" % repo)
    mvaFileName2017 = getTopTaggerTrainingFile(environ["CMSSW_BASE"] + "/src/%s/test/TopTaggerCfg_2017.cfg" % repo)
    mvaFileName2018 = getTopTaggerTrainingFile(environ["CMSSW_BASE"] + "/src/%s/test/TopTaggerCfg_2018.cfg" % repo)

    filestoTransfer = [environ["CMSSW_BASE"] + "/src/%s/test/MyAnalysis" % repo, 
                       environ["CMSSW_BASE"] + "/src/%s/test/%s" % (repo,mvaFileName2016),
                       environ["CMSSW_BASE"] + "/src/%s/test/%s" % (repo,mvaFileName2017),
                       environ["CMSSW_BASE"] + "/src/%s/test/%s" % (repo,mvaFileName2018),
                       environ["CMSSW_BASE"] + "/src/%s/test/TopTaggerCfg_2016.cfg" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/TopTaggerCfg_2017.cfg" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/TopTaggerCfg_2018.cfg" % repo,
                       environ["CMSSW_BASE"] + "/src/TopTagger/TopTagger/test/libTopTagger.so",
                       environ["CMSSW_BASE"] + "/src/%s/test/sampleSets.cfg" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/sampleCollections.cfg" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/Mass_Regression.cfg" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/DeepEventShape_2016.cfg" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/DeepEventShape_2017.cfg" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/DeepEventShape_2018pre.cfg" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/DeepEventShape_2018post.cfg" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/DeepEventShape_NonIsoMuon_2016.cfg" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/DeepEventShape_NonIsoMuon_2017.cfg" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/DeepEventShape_NonIsoMuon_2018pre.cfg" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/DeepEventShape_NonIsoMuon_2018post.cfg" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/keras_frozen_Regression.pb" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/keras_frozen_2016.pb" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/keras_frozen_2017.pb" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/keras_frozen_2018pre.pb" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/keras_frozen_2018post.pb" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/allInOne_BTagEff.root" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/allInOne_SFMean.root" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/allInOne_leptonSF_2016.root" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/allInOne_leptonSF_2017.root" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/allInOne_leptonSF_2018.root" % repo, 
                       environ["CMSSW_BASE"] + "/src/%s/test/PileupHistograms_0121_69p2mb_pm4p6.root" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/pu_ratio.root" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/PileupHistograms_2018_69mb_pm5.root" % repo, 
                       environ["CMSSW_BASE"] + "/src/%s/test/CSVv2_Moriond17_B_H.csv" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/DeepCSV_102XSF_WP_V1.csv" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/DeepCSV_2016LegacySF_WP_V1.csv" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/DeepCSV_94XSF_WP_V4_B_F.csv" % repo,
                       environ["CMSSW_BASE"] + "/src/%s/test/L1prefiring_jetpt_2017BtoF.root" % repo,
                       ]
    
    print "--------------Files to Transfer-----------------"
    for i in filestoTransfer:    
        print i
    print "------------------------------------------------"
    
    sc = SampleCollection("../sampleSets.cfg", "../sampleCollections.cfg")
    if options.dataCollections or options.dataCollectionslong:
        scl = sc.sampleCollectionList()
        for sampleCollection in scl:
            sl = sc.sampleList(sampleCollection)
            print sampleCollection
            if options.dataCollectionslong:
                sys.stdout.write("\t")
                for sample in sl:
                    sys.stdout.write("%s  "%sample[1])
                print ""
                print ""
        exit(0)
    
    datasets = []
    if options.datasets:
        datasets = options.datasets.split(',')
    else:
        print "No dataset specified"
        exit(0)
    
    fileParts = []
    fileParts.append("Universe   = vanilla\n")
    fileParts.append("Executable = run_Analyzer_condor.sh\n")
    fileParts.append("Transfer_Input_Files = %s/%s.tar.gz, %s/exestuff.tar.gz\n" % (options.outPath,environ["CMSSW_VERSION"],options.outPath))
    fileParts.append("Should_Transfer_Files = YES\n")
    fileParts.append("WhenToTransferOutput = ON_EXIT\n")
    fileParts.append("x509userproxy = $ENV(X509_USER_PROXY)\n\n")

    nFilesPerJob = options.numfile
    numberOfJobs = 0
    for ds in datasets:
        ds = ds.strip()
        # create the directory
        if not os.path.isdir("%s/output-files/%s" % (options.outPath, ds)):
            os.makedirs("%s/output-files/%s" % (options.outPath, ds))
    
        for s, n, e in sc.sampleList(ds):
            print "SampleSet:", n, ", nEvents:", e
            f = open(s)
            if not f == None:
                count = 0
                for l in f:
                    if '.root' in l:
                        count = count + 1
                for startFileNum in xrange(0, count, nFilesPerJob):
                    numberOfJobs+=1
                    outputDir = "%s/output-files/%s" % (options.outPath, ds)
                    outputFiles = [
                        "MyAnalysis_%s_%s.root" % (n, startFileNum),
                        "MyAnalysis_%s_%s_Train.root" % (n, startFileNum),
                        "MyAnalysis_%s_%s_Test.root" % (n, startFileNum),
                        "MyAnalysis_%s_%s_Val.root" % (n, startFileNum),
                    ]
                    transfer = "transfer_output_remaps = \""
                    for f_ in outputFiles:
                        transfer += "%s = %s/%s" % (f_, outputDir, f_)
                        if f_ != outputFiles[-1]:
                            transfer += "; "
                    transfer += "\"\n"                    
                    fileParts.append(transfer)
                    fileParts.append("Arguments = %s %i %i %s %s %s\n"%(n, nFilesPerJob, startFileNum, s, options.analyze, environ["CMSSW_VERSION"]))
                    fileParts.append("Output = %s/log-files/MyAnalysis_%s_%i.stdout\n"%(options.outPath, n, startFileNum))
                    fileParts.append("Error = %s/log-files/MyAnalysis_%s_%i.stderr\n"%(options.outPath, n, startFileNum))
                    fileParts.append("Log = %s/log-files/MyAnalysis_%s_%i.log\n"%(options.outPath, n, startFileNum))
                    fileParts.append("Queue\n\n")
    
                f.close()
    
    fout = open("condor_submit.txt", "w")
    fout.write(''.join(fileParts))
    fout.close()

    if not options.dataCollections and not options.dataCollectionslong:
        makeExeAndFriendsTarball(filestoTransfer, "exestuff", options.outPath)
        system("tar --exclude-caches-all --exclude-vcs -zcf %s/${CMSSW_VERSION}.tar.gz -C ${CMSSW_BASE}/.. ${CMSSW_VERSION} --exclude=src --exclude=tmp" % options.outPath)
        
    if not options.noSubmit: 
        system('mkdir -p %s/log-files' % options.outPath)
        system("echo 'condor_submit condor_submit.txt'")
        system('condor_submit condor_submit.txt')
    else:
        print "------------------------------------------"
        print "Number of Jobs:", numberOfJobs
        print "------------------------------------------"

if __name__ == "__main__":
    main()
