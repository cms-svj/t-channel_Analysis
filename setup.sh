#!/usr/bin/env bash

case `uname` in
  Linux) ECHO="echo -e" ;;
  *) ECHO="echo" ;;
esac

usage(){
	EXIT=$1
	$ECHO "setup.sh [options]"
	$ECHO
	$ECHO "Options:"
	$ECHO "-d              \tuse the developer branch of Coffea (default = 0)"
	$ECHO "-l              \tuse default LCG environment"
	$ECHO "-h              \tprint this message and exit"
	$ECHO "-n [NAME]       \toverride the name of the virtual environment (default = coffeaenv)"
	exit $EXIT
}

NAME=coffeaenv
LCG=/cvmfs/sft.cern.ch/lcg/views/LCG_101cuda/x86_64-centos7-gcc8-opt
SC=/cvmfs/unpacked.cern.ch/registry.hub.docker.com/fnallpc/fnallpc-docker:pytorch-1.9.0-cuda11.1-cudnn8-runtime-singularity
DEV=0
useLCG=0

# check arguments
while getopts "dlhn:" opt; do
	case "$opt" in
		d) DEV=1
		;;
                l) useLCG=1
                ;;
		h) usage 0
		;;
		n) NAME=$OPTARG
		;;
		:) printf "missing argument for -%s\n" "$OPTARG" >&2
		   usage -1
		;;
		\?) printf "illegal option: -%s\n" "$OPTARG" >&2
		    usage -2
		;;
	esac
done

# Setup the base environment
if [[ "$useLCG" == "1" ]]; then
        $ECHO "\nGetting the LCG environment ... "
        source $LCG/setup.sh
        pyenvflag=--copies
        NAME=${NAME}LCG
elif [[ "$SINGULARITY_CONTAINER" == "$SC" ]]; then
        $ECHO "\nBuilding env on top of Singularity container \"$SINGULARITY_CONTAINER\" ... "
        pyenvflag=--system-site-packages
else
        $ECHO "Error: Expected to be in Singularity container \"$SC\""
        $ECHO "Error: Either launch the correct Singularity container or specify that you are using LCG with the \"-l\" flag"
        exit
fi

# Finding path to env
pypath=`which python | sed 's/bin\/python//g'`

# Install most of the needed software in a virtual environment
# following https://aarongorka.com/blog/portable-virtualenv/, an alternative is https://github.com/pantsbuild/pex
$ECHO "\nMaking and activating the virtual environment ... "
python -m venv $pyenvflag $NAME
source $NAME/bin/activate

# Setting up Dask
$ECHO "\nSetup for Dask on LPC ... \n"
pyversion=$(python -c"import sys; print('{}.{}'.format(sys.version_info.major,sys.version_info.minor))")
pypackages=lib/python${pyversion}/site-packages/
siteprefix=${pypath}/${pypackages}
# need to remove python path from site to avoid dask conflicts
export PYTHONPATH=""
ln -sf ${siteprefix}/pyxrootd ${NAME}/${pypackages}/pyxrootd
ln -sf ${siteprefix}/XRootD   ${NAME}/${pypackages}/XRootD
git clone git@github.com:cms-svj/lpc_dask

# pip installing extra python packages
$ECHO "\nInstalling 'pip' packages ... \n"
python -m pip install --no-cache-dir pip --upgrade
python -m pip install --no-cache-dir dask[dataframe]==2020.12.0 distributed==2020.12.0 dask-jobqueue
python -m pip install --no-cache-dir magiconfig
if [[ "$useLCG" == "1" ]]; then
        python -m pip install --no-cache-dir torch==1.9 --upgrade        
fi
python -m pip install --no-cache-dir mt2
if [[ "$DEV" == "1" ]]; then
	$ECHO "\nInstalling the 'development' version of Coffea ... "
	python -m pip install --no-cache-dir flake8 pytest coverage
	git clone https://github.com/CoffeaTeam/coffea
	cd coffea
	python -m pip install --no-cache-dir --editable .[dask,spark,parsl]
	cd ..
else
        $ECHO "\nInstalling the 'production' version of Coffea ... "
	python -m pip install --no-cache-dir coffea[dask,spark,parsl]==0.7.14
fi

# apply patches
./patch.sh $NAME

# Clone TreeMaker for its lists of samples and files
$ECHO "\nCloning the TreeMaker repository ..."
git clone git@github.com:TreeMaker/TreeMaker.git ${NAME}/${pypackages}/TreeMaker/

# Setup the activation script for the virtual environment
$ECHO "\nSetting up the activation script for the virtual environment ... "
sed -i '40s/.*/VIRTUAL_ENV="$(cd "$(dirname "$(dirname "${BASH_SOURCE[0]}" )")" \&\& pwd)"/' $NAME/bin/activate
find $NAME/bin/ -type f -print0 | xargs -0 -P 4 sed -i '1s/#!.*python$/#!\/usr\/bin\/env python/'
sed -i "2a source ${pypath}/setup.sh"'\nexport PYTHONPATH=""' $NAME/bin/activate
sed -i "4a source ${pypath}/setup.csh"'\nsetenv PYTHONPATH ""' $NAME/bin/activate.csh

# Setting up jupyter
$ECHO "\nSetting up the ipython/jupyter kernel ... "
storage_dir=$(readlink -f $PWD)
ipython kernel install --prefix=${storage_dir}/.local --name=$NAME

# Finishing up
tar -zcf ${NAME}.tar.gz ${NAME}
deactivate
$ECHO "\nFINISHED"

