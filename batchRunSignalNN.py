from utils import samples as s
import os
import sys

reRunDetail = [
"2018_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1",  
"2018_mMed-800_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
"2018_mMed-1000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
"2018_mMed-1500_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
"2018_mMed-2000_mDark-100_rinv-0p3_alpha-peak_yukawa-1", 
"2018_mMed-2000_mDark-1_rinv-0p3_alpha-peak_yukawa-1", 
"2018_mMed-2000_mDark-20_rinv-0p1_alpha-peak_yukawa-1", 
"2018_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1", 
"2018_mMed-2000_mDark-20_rinv-0p5_alpha-peak_yukawa-1", 
"2018_mMed-2000_mDark-20_rinv-0p7_alpha-peak_yukawa-1", 
"2018_mMed-2000_mDark-50_rinv-0p3_alpha-peak_yukawa-1", 
"2018_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1", 
"2018_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
]

for rerun in reRunDetail:
    sample = rerun
    command = "python analyze_root_varModule.py -d {} -w 4 -s 1000 --training NN".format(sample) # make training files for ML
    os.system(command)
    command = "mv trainFile.root tree_SVJ_{}_NN.root".format(sample.replace("2018_",""))
    os.system(command)