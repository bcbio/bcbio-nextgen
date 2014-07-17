#!/usr/bin/env bash
# Run unit tests

BCBIO_DIR=$(dirname $(readlink -f `which bcbio_nextgen.py`))
unset PYTHONHOME
unset PYTHONPATH
git clone https://github.com/roryk/bcbio-nextgen-test-data.git
$BCBIO_DIR/nosetests -v -s -a unit
