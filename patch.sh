#!/bin/bash

VENV=$1

sed -i 's/ if issubclass(schema, schemas.BaseSchema):/ if schema is not None and issubclass(schema, schemas.BaseSchema):/' ${VENV}/lib/python3.6/site-packages/coffea/processor/executor.py

sed -i 's/"xrootd_handler": uproot4.source.xrootd.XRootDSource/"xrootd_handler": uproot4.source.xrootd.MultithreadedXRootDSource/' ${VENV}/lib/python3.6/site-packages/uproot4/reading.py
