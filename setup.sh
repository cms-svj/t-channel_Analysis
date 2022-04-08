#!/usr/bin/env bash

case `uname` in
  Linux) ECHO="echo -e" ;;
  *) ECHO="echo" ;;
esac

NAME=coffeaenv

$ECHO "\nMaking and activating the virtual environment ... "
python -m venv --system-site-packages $NAME
source $NAME/bin/activate

$ECHO "\nInstalling 'pip' packages ... "
python -m pip install --ignore-installed cython coffea==0.7.14
python -m pip install --ignore-installed mt2
python -m pip install --ignore-installed magiconfig

deactivate
$ECHO "\nFINISHED"

