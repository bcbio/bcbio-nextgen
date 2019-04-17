#!/usr/bin/env bash
# Portable resolution of symlinks http://stackoverflow.com/a/24572274/252589
# readlink -f does not work on Macs
readlinkf(){ perl -MCwd -e 'print Cwd::abs_path shift' $1;}

MARK=${1:-speed1}
if [ -n "$1" ]; then
	# Remaining args are passed unmodified to the test runner.
	shift
fi

if [ "--help" = "$MARK" -o "-h" = "$MARK" ]; then
	cat <<EOHELP

USAGE

   run_tests.sh -h|--help
   run_tests.sh <selection>

DESCRIPTION

  Simple command to run tests using installed bcbio_nextgen python/nose. Pass optional argument to select type of tests to run:

    ./run_tests.sh rnaseq
    ./run_tests.sh speed=1
    ./run_tests.sh devel
    ./run_tests.sh docker
    ./run_tests.sh docker_ipython

ENVIRONMENT

  The script unsets values for PYTHONHOME and PYTHONPATH.

  Use PYTEST to specify an alternative path to the py.test executable.

EOHELP
	exit 1
fi

if [[ ${MARK} == docker* && "`which bcbio_vm.py`" != "" ]]; then
	BCBIO_DIR=$(dirname "$(readlinkf `which bcbio_vm.py`)")
else
	BCBIO_DIR=$(dirname "$(readlinkf `which bcbio_nextgen.py`)")
fi

unset PYTHONHOME
unset PYTHONPATH
export PYTHONNOUSERSITE=1

# Ensure version.py exists in raw cloned bcbio directory
if [ -d "../bcbio/pipeline" ]; then
	[ -f ../bcbio/pipeline/version.py ] || touch ../bcbio/pipeline/version.py
fi

if [ -z "$PYTEST" ]; then
	PYTEST="$BCBIO_DIR/py.test"
fi
if [ ! -x "$PYTEST" ]; then
	echo "E: py.test (at '$PYTEST') not found or not executable."
	exit 1
fi
"$PYTEST" -p no:cacheprovider -p no:stepwise -v -s -m ${MARK} "$@"
