from utils import samples as s
import os
import sys

reRunDetail = [
# "2016_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1",  
# # "2016_mMed-800_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
# # "2016_mMed-1000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
# # "2016_mMed-1500_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
# "2016_mMed-2000_mDark-100_rinv-0p3_alpha-peak_yukawa-1", 
# "2016_mMed-2000_mDark-1_rinv-0p3_alpha-peak_yukawa-1", 
# "2016_mMed-2000_mDark-20_rinv-0p1_alpha-peak_yukawa-1", 
# "2016_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1", 
# # "2016_mMed-2000_mDark-20_rinv-0p5_alpha-peak_yukawa-1", 
# "2016_mMed-2000_mDark-20_rinv-0p7_alpha-peak_yukawa-1", 
# # "2016_mMed-2000_mDark-50_rinv-0p3_alpha-peak_yukawa-1", 
# # "2016_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1", 
# "2016_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
# "2017_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1",  
# # "2017_mMed-800_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
# # "2017_mMed-1000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
# # "2017_mMed-1500_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
# "2017_mMed-2000_mDark-100_rinv-0p3_alpha-peak_yukawa-1", 
# "2017_mMed-2000_mDark-1_rinv-0p3_alpha-peak_yukawa-1", 
# "2017_mMed-2000_mDark-20_rinv-0p1_alpha-peak_yukawa-1", 
# "2017_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1", 
# # "2017_mMed-2000_mDark-20_rinv-0p5_alpha-peak_yukawa-1", 
# "2017_mMed-2000_mDark-20_rinv-0p7_alpha-peak_yukawa-1", 
# # "2017_mMed-2000_mDark-50_rinv-0p3_alpha-peak_yukawa-1", 
# # "2017_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1", 
# "2017_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
"2018_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1",  
"2018_mMed-800_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
"2018_mMed-1000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
"2018_mMed-1500_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
"2018_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1", 
"2018_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1", 
"2018_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
"2018_mMed-2000_mDark-100_rinv-0p3_alpha-peak_yukawa-1", 
"2018_mMed-2000_mDark-1_rinv-0p3_alpha-peak_yukawa-1", 
"2018_mMed-2000_mDark-20_rinv-0p1_alpha-peak_yukawa-1", 
"2018_mMed-2000_mDark-20_rinv-0p5_alpha-peak_yukawa-1", 
"2018_mMed-2000_mDark-20_rinv-0p7_alpha-peak_yukawa-1", 
"2018_mMed-2000_mDark-50_rinv-0p3_alpha-peak_yukawa-1", 
]

# make neural network training files
# tcut = "_pre_1PSVJ" # _pre or _pre_1PSVJ
# outHistF = "output/eventTaggerOutput/trainingFiles/wpt0p9{}/".format(tcut)
# for rerun in reRunDetail:
#     sample = rerun
#     command = "python analyze_root_varModule.py -d {} -w 10 -s 1000 --training NN --outHistF {} --tcut {}".format(sample,outHistF,tcut)
#     print(command)
#     os.system(command)

# make histograms
eTagName = "d1_w7_p0i1"
evtTaggerLoc = f"utils/data/DNNEventClassifier/{eTagName}"
outHistF = f"output/eventTaggerOutput/ABCDHistograms/{eTagName}"
for rerun in reRunDetail:
    sample = rerun
    command = f"python analyze.py -d {sample} -w 15 -s 1000 --outHistF {outHistF} -C {evtTaggerLoc}/config_out.py -t {evtTaggerLoc}"
    print(command)
    os.system(command)
    command = f"mv {outHistF}/test.root {outHistF}/{sample}.root"
    print(command)
    os.system(command)

# "2016_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1",  
# "2016_mMed-800_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
# "2016_mMed-1000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
# "2016_mMed-1500_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
# "2016_mMed-2000_mDark-100_rinv-0p3_alpha-peak_yukawa-1", 
# "2016_mMed-2000_mDark-1_rinv-0p3_alpha-peak_yukawa-1", 
# "2016_mMed-2000_mDark-20_rinv-0p1_alpha-peak_yukawa-1", 
# "2016_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1", 
# "2016_mMed-2000_mDark-20_rinv-0p5_alpha-peak_yukawa-1", 
# "2016_mMed-2000_mDark-20_rinv-0p7_alpha-peak_yukawa-1", 
# "2016_mMed-2000_mDark-50_rinv-0p3_alpha-peak_yukawa-1", 
# "2016_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1", 
# "2016_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
# "2017_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1",  
# "2017_mMed-800_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
# "2017_mMed-1000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
# "2017_mMed-1500_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
# "2017_mMed-2000_mDark-100_rinv-0p3_alpha-peak_yukawa-1", 
# "2017_mMed-2000_mDark-1_rinv-0p3_alpha-peak_yukawa-1", 
# "2017_mMed-2000_mDark-20_rinv-0p1_alpha-peak_yukawa-1", 
# "2017_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1", 
# "2017_mMed-2000_mDark-20_rinv-0p5_alpha-peak_yukawa-1", 
# "2017_mMed-2000_mDark-20_rinv-0p7_alpha-peak_yukawa-1", 
# "2017_mMed-2000_mDark-50_rinv-0p3_alpha-peak_yukawa-1", 
# "2017_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1", 
# "2017_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
# "2018_mMed-600_mDark-20_rinv-0p3_alpha-peak_yukawa-1",  
# "2018_mMed-800_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
# "2018_mMed-1000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
# "2018_mMed-1500_mDark-20_rinv-0p3_alpha-peak_yukawa-1",
# "2018_mMed-2000_mDark-100_rinv-0p3_alpha-peak_yukawa-1", 
# "2018_mMed-2000_mDark-1_rinv-0p3_alpha-peak_yukawa-1", 
# "2018_mMed-2000_mDark-20_rinv-0p1_alpha-peak_yukawa-1", 
# "2018_mMed-2000_mDark-20_rinv-0p3_alpha-peak_yukawa-1", 
# "2018_mMed-2000_mDark-20_rinv-0p5_alpha-peak_yukawa-1", 
# "2018_mMed-2000_mDark-20_rinv-0p7_alpha-peak_yukawa-1", 
# "2018_mMed-2000_mDark-50_rinv-0p3_alpha-peak_yukawa-1", 
# "2018_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1", 
# "2018_mMed-4000_mDark-20_rinv-0p3_alpha-peak_yukawa-1",