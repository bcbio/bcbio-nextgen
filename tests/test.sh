#!/usr/bin/env bash
set -e

readlinkf(){ perl -MCwd -e 'print Cwd::abs_path shift' $1;}

BCBIO_DIR=$(dirname "$(readlinkf `which bcbio_nextgen.py`)")
TEST_CLASS="test_automated_analysis.py:AutomatedAnalysisTest"


echo 
cmd="$BCBIO_DIR/nosetests -v -s $TEST_CLASS.$1"
echo $cmd
$cmd
