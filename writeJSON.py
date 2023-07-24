import os 
import json
import optparse

def create_json_file(dir_path, json_file_path,write_dir,dictkey):
    # List all files in the directory
    print("dir_path = {} \n json_file_path = {} \n writedir = {}".format(dir_path,json_file_path,write_dir))
    files = os.listdir(dir_path)
    file_paths = []
    # print(files)
    # Loop through all files and add them to the dictionary
    for file_name in files:
        file_path = os.path.join(write_dir, file_name)
        file_paths.append(file_path)
    file_dict = {dictkey : file_paths}
    # Write the dictionary to a JSON file
    with open(json_file_path, 'w') as json_file:
        json.dump(file_dict, json_file,indent=4)


parser = optparse.OptionParser("usage: %prog [options]\n")
parser.add_option('-p', '--path',   dest='path',       type='string',                               help="path to the dataset")
# parser.add_option('-o', '--output', dest='outputFile', type='string', default='try.JSON',           help="Name of the output JSON file")
parser.add_option(      '--dest',   dest='destDir',    type='string', default='input/sampleJSONs/backgrounds/',  help="Destination directory")
parser.add_option('-y',             dest='year',       type='string',    default='2018',               help="Year of the data type")
parser.add_option('-k',             dest='key',        type='string',    default='',                   help="key for the dictionary")
parser.add_option(      '--isSkim', dest='isskim',  action='store_true', default=False,                help="Skim option")
# parser.add_option('-d', '--dataset',dest='dataset',    type='string', default=)
options, args = parser.parse_args()

if options.year == "2018":
    location = 'Summer20UL18/'
elif options.year == "2017":
    location = 'Summer20UL17/'
elif options.year == "2016":
    location = 'Summer20UL16/'
else: 
    print("Error")


rooteos =  "root://cmseos.fnal.gov//store/user/"
susyhad = "lpcsusyhad/SusyRA2Analysis2015/Run2ProductionV20/"+location
eos = "/eos/uscms/store/user/"

skims = "/eos/uscms/store/user/lpcdarkqcd/tchannel_UL/skims/"

if options.isskim:
    jsonFilePath = options.destDir+options.year+'/'+options.year+'_Skim_'+options.path+'.json'
    writedir = rooteos+'lpcdarkqcd/tchannel_UL/skims/'+options.year+'/'+options.path #The dir location to be written in the JSON file, the file is stored in the lpcdarkqcd area
    dirPath = eos+'lpcdarkqcd/tchannel_UL/skims/'+options.year+'/'+options.path
else:
    jsonFilePath = options.destDir+options.year+'/'+options.year+'_'+options.path+'.json'
    writedir = rooteos+susyhad+options.path #The dir location to be written in the JSON file
    dirPath = eos+susyhad+options.path


dictkey = options.year+"_"+options.key

create_json_file(dirPath,jsonFilePath,writedir,dictkey)

