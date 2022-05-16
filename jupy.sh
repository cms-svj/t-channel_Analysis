#!/bin/bash

PORT=$1
if [ -z "$PORT" ]; then
	PORT=8888
fi

# generate 32 bit random hex string
TOKEN=$(python -c 'import random; print("%030x" % random.randrange(16**30))')

echo "Server url: http://127.0.0.1:${PORT}/?token=${TOKEN}"
echo ""

jupyter notebook --ip 0.0.0.0 --no-browser --notebook-dir . --port $PORT --NotebookApp.token=$TOKEN
