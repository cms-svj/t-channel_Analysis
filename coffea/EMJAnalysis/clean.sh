#!/usr/bin/env bash

case `uname` in
  Linux) ECHO="echo -e" ;;
  *) ECHO="echo" ;;
esac

usage(){
	EXIT=$1
	$ECHO "clean.sh [options]"
	$ECHO
	$ECHO "Options:"
	$ECHO "-d              \tuse the developer branch of Coffea (default = 0)"
	$ECHO "-h              \tprint this message and exit"
	$ECHO "-n [NAME]       \toverride the name of the virtual environment (default = coffeaenv)"
	exit $EXIT
}

NAME=coffeaenv
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

$ECHO "Removing the virtual environment ... "
rm -rf ${NAME} ${NAME}.tar.gz

if [[ "$DEV" == "1" ]]; then
	$ECHO "\nRemoving the 'development' version of Coffea ... "
	rm -rf coffea
fi

$ECHO "\nRemoving the ipython/jupyter kernel ... "
storage_dir=$(readlink -f $PWD)
rm -rf ${storage_dir}/.local/share/jupyter/kernels/${NAME}

$ECHO "\nFINISHED"
