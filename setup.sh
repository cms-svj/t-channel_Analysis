#!/usr/bin/env bash

source local.sh

case `uname` in
  Linux) ECHO="echo -e" ;;
  *) ECHO="echo" ;;
esac

NAME=coffeaenv
LCG=$TCHANNEL_LCG
SC=$TCHANNEL_SC
DEV=0
useLCG=0

usage(){
	EXIT=$1
	$ECHO "setup.sh [options]"
	$ECHO
	$ECHO "Options:"
	$ECHO "-d              \tuse the developer branch of Coffea"
	$ECHO "-l              \tuse LCG environment (appends LCG to venv name)"
	$ECHO "-h              \tprint this message and exit"
	$ECHO "-n [NAME]       \toverride the name of the virtual environment (default = $NAME)"
	exit $EXIT
}

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
if [[ "$useLCG" -eq 1 ]]; then
        $ECHO "\nGetting the LCG environment ... "
        source $LCG/setup.sh
        pyenvflag=--copies
        NAME=${NAME}LCG
elif [[ "$SINGULARITY_CONTAINER" == "" ]]; then
        ./launchSingularity.sh "$0 $@"
        exit 0
else
        $ECHO "\nBuilding env on top of Singularity container \"$SINGULARITY_CONTAINER\" ... "
        pyenvflag=--system-site-packages
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

# pip installing extra python packages
$ECHO "\nInstalling 'pip' packages ... \n"
python -m pip install --no-cache-dir pip --upgrade
# set up lpcjobqueue
python -m pip install --no-cache-dir git+https://github.com/CoffeaTeam/lpcjobqueue.git@v0.2.9
python -m pip install --no-cache-dir magiconfig
if [[ "$useLCG" -eq 1 ]]; then
        python -m pip install --no-cache-dir torch==1.9 --upgrade
	python -m pip install --no-cache-dir dask[dataframe]==2020.12.0 distributed==2020.12.0 dask-jobqueue
	python -m pip install --no-cache-dir tritonclient[grpc]
	python -m pip install --no-cache-dir tritonclient[http]
fi
python -m pip install --no-cache-dir mt2
if [[ "$DEV" -eq 1 ]]; then
	$ECHO "\nInstalling the 'development' version of Coffea ... "
	python -m pip install --no-cache-dir flake8 pytest coverage
	git clone https://github.com/CoffeaTeam/coffea
	cd coffea
	python -m pip install --no-cache-dir --editable .[dask,spark,parsl]
	cd ..
else
        $ECHO "\nInstalling the 'production' version of Coffea ... "
	# python -m pip install --no-cache-dir coffea[dask,spark,parsl]==0.7.17
fi
# python -m pip install --no-cache-dir numpy==1.24 # older numpy version can cause error
# Clone TreeMaker for its lists of samples and files
$ECHO "\nCloning the TreeMaker repository ..."
git clone https://github.com/TreeMaker/TreeMaker.git TreeMaker

# Setup the activation script for the virtual environment
$ECHO "\nSetting up the activation script for the virtual environment ... "
sed -i '40s/.*/VIRTUAL_ENV="$(cd "$(dirname "$(dirname "${BASH_SOURCE[0]}" )")" \&\& pwd)"/' $NAME/bin/activate
find $NAME/bin/ -type f -print0 | xargs -0 -P 4 sed -i '1s/#!.*python$/#!\/usr\/bin\/env python/'
if [[ "$useLCG" -eq 1 ]]; then
	sed -i "2a source ${pypath}/setup.sh"'\nexport PYTHONPATH=""' $NAME/bin/activate
	sed -i "4a source ${pypath}/setup.csh"'\nsetenv PYTHONPATH ""' $NAME/bin/activate.csh
fi

# Setting up jupyter
$ECHO "\nSetting up the ipython/jupyter kernel ... "
export TMPDIR=$(mktemp -d -p .)
python -m ipykernel install \
        --user                             \
        --name coffea-svj                  \
        --display-name "coffea for SVJ"    \
        --env PYTHONPATH /srv:$PWD  \
        --env TCHANNEL_BASE /srv               \
        --env PYTHONNOUSERSITE 1
rm -rf $TMPDIR && unset TMPDIR

# Finishing up
tar -zcf ${NAME}.tar.gz ${NAME}
deactivate
$ECHO "\nFINISHED"

