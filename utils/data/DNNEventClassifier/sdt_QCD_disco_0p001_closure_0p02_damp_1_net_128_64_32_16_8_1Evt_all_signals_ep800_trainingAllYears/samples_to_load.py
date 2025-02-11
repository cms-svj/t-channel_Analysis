#Here need to define signal and background samples to process
import os
cwd = os.getcwd()

__years = [
    "2016",
    "2017",
    "2018",
]

__signal_names = [
    "training_t-channel_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
    "training_t-channel_mMed-800_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
    "training_t-channel_mMed-1000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
    "training_t-channel_mMed-1500_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
    "training_t-channel_mMed-2000_mDark-20_rinv-0p1_alpha-peak_yukawa-1",
    "training_t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
    "training_t-channel_mMed-2000_mDark-20_rinv-0p5_alpha-peak_yukawa-1",
    "training_t-channel_mMed-2000_mDark-20_rinv-0p7_alpha-peak_yukawa-1",
    "training_t-channel_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
    "training_t-channel_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
]

__bkgs = [
#    #"QCD_Pt_170to300",  # Not enough MC events passing selections - very high weights!
    "QCD_Pt_300to470",
    "QCD_Pt_470to600",
    "QCD_Pt_600to800",
    "QCD_Pt_800to1000",
    "QCD_Pt_1000to1400",
    "QCD_Pt_1400to1800",
    "QCD_Pt_1800to2400",
    "QCD_Pt_2400to3200",
    "QCD_Pt_3200toInf",
]

# Important note!!!
# The key in background list will be used to merge bins together simply
# by checking if the key is present in the name of the bin!
# E.g. the QCD_Pt_* bins will be merged into the QCD background
backgrounds_list = [
    "2016_QCD",
    "2017_QCD",
    "2018_QCD",
]  # Will be read in samples_nfiles.py

tt_stitch = True

lumi_per_year = {
    "2016": 36.31e3,
    "2017": 41.48e3,
    "2018": 59.83e3,
}

# If not in in dict, no k-factor
k_factors = {
    "2016_QCD": 1.279,
    "2017_QCD": 1.464,
    "2018_QCD": 1.510,
}

samples = {
    "signals_samples" : {
        f"{year}_{signal_name}": f"{cwd}/samples/t_channel/dataset7_allBackground/t_channel_pre_selection/training_set1/signal/{year}_{signal_name}_nominal_files.txt"
        for signal_name in __signal_names
        for year in __years
    },
    "background_samples": {
        f"{year}_{bin_name}": f"{cwd}/samples/t_channel/dataset7_allBackground/t_channel_pre_selection/training_set1/background/{year}_{bin_name}_files.txt"
        for bin_name in __bkgs
        for year in __years
    },
}

## Add systematic variations
#systematics = [
#    "jer_up",
#    "jer_down",
#    "jec_up",
#    "jec_down",
#]
#
#for systematic in systematics:
#    samples["signals_samples"].update({
#        f"{signal_name}_{systematic}": f"{cwd}/samples/t_channel/dataset8/training_set4/signal/{signal_name}_{systematic}_files.txt"
#        for signal_name in __signal_names
#    })
#


n_files = {
    "signals_samples" : {k: -1 for k in samples["signals_samples"].keys()},
    "background_samples" : {k: -1 for k in samples["background_samples"].keys()},
}

