# t-channel_Analysis

This folder houses the analysis code used to look for semi-visible jets in t-channel.

## Setup Coffea Environment

The Coffea environment can be setup in two ways. One makes use of the python virtual environment mechanism. The other uses a Docker image maintained by the Coffea development team. Both methods can be used on one's private computer or on a remote computing cluster. The use of either Docker or Singularity will be determined by the virtualization software allowed/installed on the system. These instructions will assume you are using the cmslpc cluster.

### Logging into cluster

Begin by logging into `cmslpc`:
```bash
ssh -L localhost:8NNN:localhost:8NNN  <username>@cmslpc-sl7.fnal.gov
```
Remember to replace `<username>` by your username. The `-L` option is needed if you'd like to work with Jupyter. Replace the `NNN` in `8NNN` with a unique number.
If you are not going to work with Jupyter, you can just do:
```bash
ssh -L localhost:8NNN:localhost:8NNN  <username>@cmslpc-sl7.fnal.gov
```
### Python Virtual Environment

To begin the initial setup, run the following commands:
```bash
cd <working_directory>
git clone git@github.com:cms-svj/t-channel_Analysis.git
cd t-channel_Analysis
./setup.sh
```
Remember to replace `<working_directory>` with the directory where you want your files/folders to appear. You can change the name of the virtual environment by using the `-n` option and you can use the development version of coffea by using the `-d` option. These commands only need to be run during the initial setup. When doing your day-to-day tasks, you can skip these.

To activate the `coffeaenv` environment and set the Jupyter paths, run the command (every time):
```bash
cd <working_directory>/t-channel_Analysis
source init.sh
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

### Docker Image

The directions below will assume you are working on the `cmslpc` cluster, but can be easily adapted for a different computing cluster or local machine. If you're on a local machine, it's likely that all of the singularity commands can be replaced by Docker commands.

Set cache directory for singularity, go to that directory, and get the docker container for `coffea-dask` and the t-channel SVJ software:
```bash
export SINGULARITY_CACHEDIR=/uscms_data/d1/`whoami`/.singularity
cd /uscms_data/d1/`whoami`/<working_directory>/
singularity pull docker://coffeateam/coffea-dask:latest
singularity shell -B ${PWD}:/work docker://coffeateam/coffea-dask:latest
cd /work
git clone git@github.com:cms-svj/t-channel_Analysis.git
```
Remember to replace `<working_directory>`. You only need to do the `pull` command if updating the Docker image. Otherwise you can skip that command and go straight to the `shell` command.

If you'd like to work with Jupyter, perform the following commands from within the container to set the Jupyter paths:
```bash
export JUPYTER_PATH=/work/.jupyter
export JUPYTER_RUNTIME_DIR=/work/.local/share/jupyter/runtime
export JUPYTER_DATA_DIR=/work/.local/share/jupyter
export IPYTHONDIR=/work/.ipython
```

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
### Running the Analysis
To make histograms locally using the signal and background ntuples, make sure you are in `t-channel_Analysis`
```bash
python analyze.py -d <sample label> -N <number of files>
```
This is usually done for debugging and testing purposes. The list of sample labels can be found in `input/sampleLabels.txt`.
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

To make histograms on condor, cd into the `condor` directory and run
```bash
python condorSubmit.py -d 2018_QCD,2018_TTJets,2018_WJets,2018_ZJets,2018_mMed -n 10 -w 1 --output testDir
```
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
