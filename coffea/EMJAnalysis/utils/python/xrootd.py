import importlib, json
from XRootD import client
from XRootD.client.flags import DirListFlags, StatInfoFlags, OpenFlags, MkDirFlags, QueryCode

def get_files_xrootd(endpoint,base_path,dicts=[],verbose=0,debug=None):
    """
    Get a dictionary of sample names and files from an XRootD accessible storage site.
    The sample names (strings) will be the keys and the values will be a list of filenames.
    The required inputs will be:
        endpoint (str): an XRootD endpoint in the format \"root://<endpoint>[:port]/\"
        base_path (str): the base path of the sample locations (i.e. \"/store/user/<username>\")
        dicts (list(str)): The TreeMaker condorSub dictionary of samples (i.e. \"EMJ2016\")
    The optional arguments are:
        verbose (int): the verbosity can be increased by specifying a number > 0. Acceptable values are:
            1. Print the final dictionary
            2. The above plus print the result of the dirlist of each folder
        debug (tuple(int,int)): If set, this will limit the number of results output in the dictionary
            - The first int specifies the number of samples (keys) to keep
            - The second int specifies the number of files per sample to keep
 
    Example::
        samples = get_files_xrootd( endpoint = "root://cmseos.fnal.gov/",
                                    base_path = "/store/user/<username>/",
                                    verbose=1,
                                    debug=(2,2)
                                  )
    """
    
    # Setup the XRootD client
    xrdfs = client.FileSystem(endpoint)

    # Get the sample data
    sample_dict = {}
    for d in dicts:
        dict_path=("TreeMaker/Production/test/condorSub/dict_%s" % d).replace('/','.')
        flist = importlib.import_module(dict_path).flist
        
        for sample in flist['samples']:
            for folder in sample:
                # Equivalent to 'xrdfs root://hepxrd01.colorado.edu:1094/ ls /store/user/aperloff/ExoEMJAnalysis2020/PrivateSamples/<dataset_name>/'
                sample_folder, sample_name = folder.split('.')
                path = "%s%s/%s" % (base_path,sample_folder,sample_name)
                path_with_endpoint = "%s%s" % (endpoint,path)
                status, listing = xrdfs.dirlist(path,DirListFlags.STAT)
                print 
                if status.status == 0:
                    sample_dict[sample_name]=[("%s%s" % (path_with_endpoint,entry.name)) for entry in listing]
                    if verbose > 1: print(str(sample_name) + " : " + str(len(sample_dict[sample_name])) + " files")
                else:
                    if verbose > 1: print("Status: %s  Unable to get files from %s" % (status.status,path_with_endpoint))

    # Limit the processes and files if debug is set to True
    if debug != None:
        sample_dict_debug = {}
        if debug[0]>len(sample_dict):
            raise IndexError("The first debug argument ("+str(debug[0])+") is larger than the number of samples found.")
        for i in range(0,debug[0]):
            key   = list(sample_dict.keys())[i]
            value = list(sample_dict.values())[i][0:debug[1]]
            sample_dict_debug[key] = value
        sample_dict = sample_dict_debug                    

    if verbose > 0: print(json.dumps(sample_dict,sort_keys=True, indent=4))
    return sample_dict