#!/usr/bin/env bash
# Simple command to run tests using installed bcbio_nextgen python/nose
# Can pass an optional argument with the type of tests to run
# ./run_tests.sh rnaseq
# ./run_tests.sh speed=1
# ./run_tests.sh devel

# Portable resolution of symlinks http://stackoverflow.com/a/24572274/252589
# readlink -f does not work on Macs
readlinkf(){ perl -MCwd -e 'print Cwd::abs_path shift' $1;}

ATTR=${1:-speed=1}
if [ -n "$1" ]; then
	# Remaining args are passed unmodified to the test runner.
	shift
fi
if [[ "$ATTR" == docker* && "`which bcbio_vm.py`" != "" ]]; then
    BCBIO_DIR=$(dirname "$(readlinkf `which bcbio_vm.py`)")
else
    BCBIO_DIR=$(dirname "$(readlinkf `which bcbio_nextgen.py`)")
fi
unset PYTHONHOME
unset PYTHONPATH
export PYTHONNOUSERSITE=1
"$BCBIO_DIR/nosetests" -v -s -a $ATTR "$@"
