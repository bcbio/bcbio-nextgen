#!/usr/bin/env bash
# Simple command to run tests using installed bcbio_nextgen python/nose
# Can pass an optional argument with the type of tests to run
# ./run_tests.sh rnaseq
# ./run_tests.sh speed=1
# ./run_tests.sh devel

ATTR=${1:-speed=1}
BCBIO_DIR=$(dirname $(readlink -f `which bcbio_nextgen.py`))
$BCBIO_DIR/nosetests -v -s -a $ATTR
