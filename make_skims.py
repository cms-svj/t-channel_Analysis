import ROOT
import optparse
import numpy as np
import os
from utils import samples as s

def skimmer(fileset,sample,cutstring,outputDir,output,isCondor):
    for file in fileset[sample]:
        print("file = ", file) 
        inputFileScope = ROOT.TFile.Open(file)
        inputTreeScope = inputFileScope.Get("TreeMaker2/PreSelection")
        if not inputTreeScope:
            print(" Empty file - {}".format(file))
            inputFileScope.Close()
            continue
        outputfileName = output+file.split('/')[-1]
        print("filenames - ",outputfileName)
        outputFileScope = ROOT.TFile(outputfileName,"RECREATE")
        outputTreeScope = inputTreeScope.CopyTree(cutstring)                                                                                                                                       
        outputTreeScope.Write()
        outputFileScope.Close()
        inputFileScope.Close()
        # copy the file to the 
        # if isCondor:
            # eosStore = "/eos/uscms/store/user/lpcdarkqcd/tchannel_UL/skims/"
            # print("Output Dir - ",outputDir)
            # if not os.path.isdir(eosStore+outputDir):
            #     os.makedirs(eosStore+outputDir)
            # os.system("xrdcp -f "+outputfileName+" root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims/"+outputDir+".")
            # remove the file created
            # os.system("rm "+outputfileName)


def main():
    # Get options from the command line
    parser = optparse.OptionParser()
    parser.add_option('-d','--dataset', help='dataset',dest='dataset', type=str, default="2018_METData")
    parser.add_option('-n', '--nFiles',  help='nFiles', dest='nFiles',  type=int, default=-1)
    parser.add_option('-M', '--startFile', help='startFile',         dest='startFile', type=int, default=0)
    parser.add_option('-o', '--output',  help='output',dest='output',  type='string',default='Skim_try_')
    parser.add_option(      '--condor',    help='running on condor', dest='condor',              default=False, action='store_true')
    options, args = parser.parse_args()

    sample = options.dataset

    # Configure file name is using condor
    output = "Skim_" if options.condor else options.output

    # get dictionary of files from the input samples
    fileset = s.getFileset(sample, True, options.startFile, options.nFiles)
    year, dataset = options.dataset.split("_")
    outputDir = year+'/'+dataset+'/'
    eosStore = 'root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims/'

    # Define cutstrings here - 
    cutstring = "((HT+MET)>1200) && (@JetsAK8.size() > 1)"
    # cutstring = "(HT+MET)>800"

    # call the skimmer 
    skimmer(fileset,sample,cutstring,outputDir,output,options.condor)

if __name__ == "__main__":
    main()



