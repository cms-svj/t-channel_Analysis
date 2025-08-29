import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import uproot
import h5py
import yaml
import pickle
import mplhep as hep
import awkward as ak
import sklearn.metrics as metrics

plt.rcParams['figure.dpi'] = 100



filename = "root://cmseos.fnal.gov//store/user/lpcdarkqcd/tchannel_UL/skims/2017/t_channel_pre_selection/nominal//HTMHT/part-0.root"




#TT_filename = "/uscms/home/nbruhwil/nobackup/CMSDAS2025/StauLongExercise/CMSSW_13_0_10/src/final_mc_ntuples/TTTo2L2Nu_postEE.root"


#left_signal_file = uproot.open(left_signal_filename)
#right_signal_file = uproot.open(right_signal_filename)
#DY_file = uproot.open(DY_filename)
#TT_file = uproot.open(TT_filename)
#WW_file = uproot.open(WW_filename)
#files = [(left_signal_file, "left_stau"), (right_signal_file, "right_stau"), (DY_file, "DY"), (TT_file, "TT"), (WW_file, "WW")]

file = uproot.open(filename)
print("File keys: ", file.keys())
print()
print("Events: ", file["Events"].keys())
print()



NNscore = file["Events"]["JetsAK8_pNetJetTaggerScore"].array()
print(f" NN Score range:", min(NNscore),max(NNscore))                                        

plt.hist(NNscore, 50,(-1.1,1.1), histtype="step", density=True)
plt.title("NN Score")
plt.xlabel("SCore")
plt.ylabel("Events")
plt.legend()

print(NNscore)
plt.legend()                                
