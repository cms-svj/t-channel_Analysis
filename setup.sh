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
	$ECHO "-h              \tprint this message and exit"
	$ECHO "-n [NAME]       \toverride the name of the virtual environment (default = coffeaenv)"
	exit $EXIT
}

NAME=coffeaenv
LCG=/cvmfs/sft.cern.ch/lcg/views/LCG_96python3/x86_64-centos7-gcc8-opt
DEV=0

# check arguments
while getopts "dhn:" opt; do
	case "$opt" in
		d) DEV=1
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

# Setup the LCG environment
$ECHO "Getting the LCG environment ... "
source $LCG/setup.sh

# Install most of the needed software in a virtual environment
# following https://aarongorka.com/blog/portable-virtualenv/, an alternative is https://github.com/pantsbuild/pex
$ECHO "\nMaking and activiating the virtual environment ... "
python -m venv --copies $NAME
source $NAME/bin/activate
$ECHO "\nInstalling 'pip' packages ... "
python -m pip install --no-cache-dir setuptools pip --upgrade
python -m pip install --no-cache-dir xxhash
python -m pip install --no-cache-dir uproot4
if [[ "$DEV" == "1" ]]; then
	$ECHO "\nInstalling the 'development' version of Coffea ... "
	python -m pip install --no-cache-dir flake8 pytest coverage
	git clone https://github.com/CoffeaTeam/coffea
	cd coffea
	python -m pip install --no-cache-dir --editable .[dask,spark,parsl]
	cd ..
else
	$ECHO "Installing the 'production' version of Coffea ... "
	python -m pip install --no-cache-dir coffea[dask,spark,parsl]
fi

# Clone TreeMaker for its lists of samples and files
$ECHO "\nCloning the TreeMaker repository ..."
git clone git@github.com:TreeMaker/TreeMaker.git ${NAME}/lib/python3.6/site-packages/TreeMaker/

# Setup the activation script for the virtual environment
$ECHO "\nSetting up the activation script for the virtual environment ... "
sed -i '40s/.*/VIRTUAL_ENV="$(cd "$(dirname "$(dirname "${BASH_SOURCE[0]}" )")" \&\& pwd)"/' $NAME/bin/activate
find coffeaenv/bin/ -type f -print0 | xargs -0 -P 4 sed -i '1s/#!.*python$/#!\/usr\/bin\/env python/'
sed -i "2a source ${LCG}/setup.sh" $NAME/bin/activate
sed -i "4a source ${LCG}/setup.csh" $NAME/bin/activate.csh

$ECHO "\nSetting up the ipython/jupyter kernel ... "
storage_dir=$(readlink -f $PWD)
ipython kernel install --prefix=${storage_dir}/.local --name=$NAME
tar -zcf ${NAME}.tar.gz ${NAME}

deactivate
$ECHO "\nFINISHED"
