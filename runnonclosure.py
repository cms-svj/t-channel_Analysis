import subprocess

commands = [
    "python3 MeTclosure.py -d output/Preselection_DNNallbckg -y 2016 -o 2016_preselection_DNNtrainedallbckgs",
    "python3 DNNclosure.py -d output/Preselection_DNNallbckg -y 2016 -o 2016_preselection_DNNtrainedallbckgs",
    "python3 VRIIInonclosure.py -d output/Preselection_DNNallbckg -y 2016 -o 2016_preselection_DNNtrainedallbckgs",
    "python3 SignalRegionCloser.py -d output/Preselection_DNNallbckg -y 2016 -o 2016_preselection_DNNtrainedallbckgs",
    "python3 SignalRegionCloser_MET.py -d output/Preselection_DNNallbckg -y 2016 -o 2016_preselection_DNNtrainedallbckgs",
    "python3 SignalRegionClosureVRIII.py -d output/Preselection_DNNallbckg -y 2016 -o 2016_preselection_DNNtrainedallbckgs",


    
    "python3 2Dnonclosure.py -d output/Preselection_DNNallbckg -y 2016 -o 2016_preselection_DNNtrainedallbckgs",
    "python3 2Dnonclosure.py -d output/Preselection_DNNallbckg -y 2017 -o 2017_preselection_DNNtrainedallbckgs"
    "python3 2Dnonclosure.py -d output/Preselection_DNNallbckg -y 2018 -o 2018_preselection_DNNtrainedallbckgs",



]

for cmd in commands:
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True)