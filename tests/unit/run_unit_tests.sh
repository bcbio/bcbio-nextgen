#!/usr/bin/env bash
# Run unit tests

BCBIO_DIR=$(dirname $(readlink -f `which bcbio_nextgen.py`))
unset PYTHONHOME
unset PYTHONPATH
$BCBIO_DIR/python manage_unit_test_data.py
$BCBIO_DIR/nosetests -v -s -a unit
