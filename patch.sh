#!/bin/bash

VENV=$1

pyversion=$(python -c"import sys; print('{}.{}'.format(sys.version_info.major,sys.version_info.minor))")

sed -i 's/ if issubclass(schema, schemas.BaseSchema):/ if schema is not None and issubclass(schema, schemas.BaseSchema):/'         ${VENV}/lib/python${pyversion}/site-packages/coffea/processor/executor.py
sed -i 's/"xrootd_handler": uproot4.source.xrootd.XRootDSource/"xrootd_handler": uproot4.source.xrootd.MultithreadedXRootDSource/' ${VENV}/lib/python${pyversion}/site-packages/uproot4/reading.py
