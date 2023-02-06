# t-channel_Analysis

This folder houses the analysis code used to look for semi-visible jets in t-channel.

## Setup Coffea and Torch Environment

The Coffea+Torch environment can be setup in two ways. Each method adds a python virtual environment to install the python packages that are missing. The first option is based on a singularity container made for machine learning and the second option is built on a LCG environment. These instructions will assume you are using the cmslpc cluster.

### Logging into cluster

Begin by logging into `cmslpc`:
```bash
ssh -L localhost:8NNN:localhost:8NNN <username>@cmslpc-sl7.fnal.gov
```
Remember to replace `<username>` by your username. The `-L` option is needed if you'd like to work with Jupyter. Replace the `NNN` in `8NNN` with a unique number.
If you are not going to work with Jupyter, you can omit `-L localhost:8NNN:localhost:8NNN`.

### Python Virtual Environment

To begin the initial setup, run the following commands:
```bash
cd <working_directory>
git clone git@github.com:cms-svj/t-channel_Analysis.git
cd t-channel_Analysis
./setup.sh
```
Remember to replace `<working_directory>` with the directory where you want your files/folders to appear. You can change the name of the virtual environment by using the `-n` option and you can use the development version of coffea by using the `-d` option. These commands only need to be run during the initial setup. When doing your day-to-day tasks, you can skip these. The LCG method can be setup by skipping the Singularity launching step and by adding a `-l` flag to the setup.sh script.

To activate the `coffeaenv` environment and set the Jupyter paths, run the command (every time):
```bash
cd <working_directory>/t-channel_Analysis
source init.sh
```

If you're using LCG, the last line above changes:
```
source initLCG.sh
```

When you are done working and would like to ``de-activate'' the `coffeaenv` environment, run the command:
```bash
deactivate
```
This shell function was given to you by the virtual environment.

To remove the virtual environment and the associated files (i.e. inverse of the setup script), you can use the run the following command:
```bash
cd <working_directory>/t-channel_Analysis
./clean.sh
```
The `clean.sh` script has the same `-n` and `-d` options as in the `setup.sh` script.

### Using Jupyter

To launch a Jupyter server (replace `8NNN` with the unique number):
```bash
./jupy.sh 8NNN
```

To open the Jupyter interface in your local web browser, navigate to the "Server url: http://127.0.0.1:8NNN/?token=..." printed in the terminal.

### Common Good Practices

It's a good idea to get/renew a voms ticket if you're going to be working with XRootD/GFAL or other grid tools:
```bash
voms-proxy-init -voms cms --valid 192:00
```

### Uploading New ParticleNet Model
To upload a new particleNet model, you need to copy and paste the trained model (name the model as `particleNetModel.pth`) to the `t-channel_Analysis` directory. You may need to update `utils/data/GNNTagger/svj.yaml`.
1. Convert the particleNet model with jit:
```
python makeJitModel.py
```
The jitted model is called `model.pt`.

2. Upload `model.pt` to `triton/svj_tch_gnn/1` in [this container](https://test-burt-3.okddev.fnal.gov/buckets/triton-models/browse). You will need to connect to [Fermilab VPN](https://redtop.fnal.gov/guide-to-vpn-connections-to-fermilab/) before accessing the container.

You may also need to ask [Burt Holzman](https://computing.fnal.gov/burt-holzman/) for permission to access the container. 

### Running the Analysis
Currently the default setup does not work for the neural network. So instead of doing `source init.sh` while you are in `t-channel_Analysis`,
you should follow the steps in https://github.com/cms-svj/t-channel_Analysis/issues/22 to set up the environment.
Also need to add `python -m pip install --ignore-installed magiconfig` to the set up steps above.
After that, you will also need to make soft links to the neural network file and the npz file that contains the normalization information.
For now you can do the following:
```
ln -s /uscms/home/keanet/nobackup/SVJ/Tagger/SVJTaggerNN/logs/test_tch_normMeanStd/net.pth net.pth
ln -s /uscms/home/keanet/nobackup/SVJ/Tagger/SVJTaggerNN/logs/test_tch_normMeanStd/normMeanStd.npz normMeanStd.npz
```
To make histograms locally using the signal and background ntuples, make sure you are in `t-channel_Analysis`.
```bash
python analyze.py -d <sample label> -N <number of files>
```
To make neural network training files locally using the signal and background ntuples, make sure you are using the LCG coffea environment. The way the code is set up now, we need to `hadd` the training root files together, but the singularity environment does not support `hadd`
```bash
python analyze_root_varModule.py -d <sample label> -N <number of files>
```
The output files are called `test.root` and `trainFile.root` respectively.
<sample label> can be anything in the `input/sampleLabels.txt`, but be careful with the t-channel signals, because the JSON files contain both the pair production and full t-channel files. So for example, `2018_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1` alone will probably make the code run over the pair production sample. We need to update `utils/samples.py` to make it smarter, but now just use <sample label> that are more specific. For example, `2018_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1_13TeV-madgraphMLM` will grab the full t-channel samples, while `2018_mMed-3000_mDark-20_rinv-0p3_alpha-peak_yukawa-1_13TeV-pythia8` will grab the pair production samples.

Running things locally is usually done for debugging and testing purposes.
* `-d`: sample labels for list of input files to run over, which can be found in `input/sampleLabels.txt`.
* `-N`: number of files from the sample to run over. Default is -1.
* `-M`: index of the first file to run over.
* `-w`: number of workers.
* `-s`: chunksize; an input for the coffea processor. It determines how many events to process at the same time.  
Sometimes it is helpful to use `-w 1 -s 2` for debugging.

Dask can be used to run `analyze.py` in parallel, either locally or on Condor. The relevant options are:
* `--dask`: run w/ dask
* `--port`: port for dask status dashboard
* `--mincores`: dask waits for min # cores
* `--quiet`: suppress status printouts  
To view the status dashboard, specify `--port 8NNN` (using the forwarded port from the earlier ssh command)
and navigate to `localhost:8NNN` in a web browser.

To run jobs on condor, cd into the `condor` directory and run
```bash
source initCondor.sh
python condorSubmit.py -d 2018_QCD,2018_mMed,2018_TTJets,2018_WJets,2018_ZJets -n 5 -w 1 --output [output directory] -p --pout [eos output directory for storing the training files]
```
* actually the argument after `--output` does not do anything in this case.
* also be careful when using `2018_mMed` after the `-d` flag, because the JSON files contain both the pair production and full t-channel signals, so those files have very similar names, and `utils/samples.py` may not grab the desired inputs. Still need to work on `samples.py` to make it smarter.

This will run over all the backgrounds (QCD, TTJets, WJets, ZJets) and the t-channel signals (the s-channel signals are labeled as 2016_mZprime,2017_mZprime). -n 10 means each job will use 10 root files as input, while -w 1 means we are using 1 CPU per job. An higher number of CPU used will use too much memory causing the job to be held, while a higher number of input files can make the job run longer and may also cause memory issue. After the jobs have finished running, the output histogram root files should be in `condor/testDir` (set by the --output flag).
* `-d`: sample labels for list of input files to run over. Can use the labels found in `input/sampleLabels.txt` or more general labels such as 2018_QCD.
* `-n`: number of files from the sample to run over. Default is -1.
* `-c`: do not submit jobs to condor, but can still see the number of jobs that would have been submitted.
* `-w`: number of workers.
* `--output`: output directory.

To hadd all the output files produced from the previous step, make sure you are in the `condor` directory and run
```bash
python hadder.py -H testHadd2 -p testDir/output-files -o
```
* `-H`: output directory.
* `-p`: input directory.
* `-o`: overwrite existing files with the same names.

To plot the histograms, cd back into `t-channel_Analysis` and run
```bash
python plotStack.py
```
The plots are stored in `plots`.
