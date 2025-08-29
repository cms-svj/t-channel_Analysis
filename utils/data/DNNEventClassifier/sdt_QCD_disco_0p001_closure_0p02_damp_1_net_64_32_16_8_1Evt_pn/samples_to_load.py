#Here need to define signal and background samples to process
import os
cwd = os.getcwd()

__signal_names = [
#    "t-channel_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
#    "t-channel_mMed-800_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
    "t-channel_mMed-1000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
#    "t-channel_mMed-1500_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
#    "t-channel_mMed-2000_mDark-1_rinv-0p3_alpha-peak_yukawa-1",
    "t-channel_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
#    "t-channel_mMed-2000_mDark-50_rinv-0p3_alpha-peak_yukawa-1",
#    "t-channel_mMed-2000_mDark-100_rinv-0p3_alpha-peak_yukawa-1",
    "t-channel_mMed-2000_mDark-20_rinv-0p1_alpha-peak_yukawa-1",
    "t-channel_mMed-2000_mDark-20_rinv-0p5_alpha-peak_yukawa-1",
    "t-channel_mMed-2000_mDark-20_rinv-0p7_alpha-peak_yukawa-1",
    "t-channel_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
    "t-channel_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
]

__qcd_bins = [
    "QCD_Pt_170to300",
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
__backgrounds_list = ["QCD"]  # Will be read in samples_nfiles.py


samples = {
    "signals_samples" : {
        signal_name: f"{cwd}/samples/t_channel/dataset6_QCD_pn/signal/{signal_name}_files.txt"
        for signal_name in __signal_names
    },
    "background_samples": {
        bin_name: f"{cwd}/samples/t_channel/dataset6_QCD_pn/background/{bin_name}_files.txt"
        for bin_name in __qcd_bins
    },
}

