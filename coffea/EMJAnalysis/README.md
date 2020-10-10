# EMJAnalysis

This folder houses the analysis code used to look for emerging jets.

## Setup Coffea Environment

The Coffea environment can be setup in two ways. One makes use of the python virtual environment mechanism. The other uses a Docker image maintained by the Coffea development team. Both methods can be used on one's private computer or on a remote computing cluster. The use of either Docker or Singularity will be determined by the virtualization software allowed/installed on the system. These instructions will assume you are using the cmslpc cluster.

### Logging into cluster

Begin by logging into `cmslpc`:
```bash
ssh -L localhost:8NNN:localhost:8NNN  <username>@cmslpc-sl7.fnal.gov
```
Remember to replace `<username>` by your username. The `-L` option is needed if you'd like to work with Jupyter. Replace the `NNN` in `8NNN` with a unique number.

### Python Virtual Environment

To begin the initial setup, run the following commands:
```bash
cd <working_directory>
git clone git@github.com:cms-svj/t-channel_Analysis.git
cd t-channel_Analysis/coffea/EMJAnalysis
./setup.sh
```
Remember to replace `<working_directory>` with the directory where you want your files/folders to appear. You can change the name of the virtual environment by using the `-n` option and you can use the development version of coffea by using the `-d` option. These commands only need to be run during the initial setup. When doing your day-to-day tasks, you can skip these.

To activate the `coffeaenv` environment and set the Jupyter paths, run the command (every time):
```bash
cd <working_directory>/EMJAnalysis
source init.sh
```

When you are done working and would like to ``de-activate'' the `coffeaenv` environment, run the command:
```bash
deactivate
```
This shell function was given to you by the virtual environment.

To remove the virtual environment and the associated files (i.e. inverse of the setup script), you can use the run the following command:
```bash
cd <working_directory>/EMJAnalysis
./clean.sh
```
The `clean.sh` script has the same `-n` and `-d` options as in the `setup.sh` script.

### Docker Image

The directions below will assume you are working on the `cmslpc` cluster, but can be easily adapted for a different computing cluster or local machine. If you're on a local machine, it's likely that all of the singularity commands can be replaced by Docker commands.

Set cache directory for singularity, go to that directory, and get the docker container for `coffea-dask` and the EMJ software:
```bash
export SINGULARITY_CACHEDIR=/uscms_data/d1/`whoami`/.singularity
cd /uscms_data/d1/`whoami`/<working_directory>/
singularity pull docker://coffeateam/coffea-dask:latest
singularity shell -B ${PWD}:/work docker://coffeateam/coffea-dask:latest
cd /work
git clone ssh://git@gitlab.cern.ch:7999/cms-emj/EMJAnalysis.git
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
