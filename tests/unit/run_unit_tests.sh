#!/usr/bin/env bash
# Run unit tests

BCBIO_DIR=$(dirname $(readlink -f `which bcbio_nextgen.py`))
unset PYTHONHOME
unset PYTHONPATH
if [ -d "bcbio-nextgen-test-data" ]; then
    cd bcbio-nextgen-test-data; git pull; cd ../
else
    git clone https://github.com/roryk/bcbio-nextgen-test-data.git
fi
$BCBIO_DIR/nosetests -v -s -a unit
