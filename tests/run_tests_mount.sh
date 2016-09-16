#!/usr/bin/env bash
# Simple command to run tests using installed bcbio_nextgen python/nose
# Can pass an optional argument with the type of tests to run
# ./run_tests.sh rnaseq
# ./run_tests.sh speed=1
# ./run_tests.sh devel
# ./run_tests.sh docker
# ./run_tests.sh docker_ipython

# Portable resolution of symlinks http://stackoverflow.com/a/24572274/252589
# readlink -f does not work on Macs

BUCKET='testbucket'
readlinkf(){ perl -MCwd -e 'print Cwd::abs_path shift' $1;}
checkmounted(){ grep -qs "$1" "/proc/mounts";}
ATTR=${1:-speed=1}
set -e

if [ -n "$1" ]; then
	# Remaining args are passed unmodified to the test runner.
	shift
fi

if checkmounted "$BUCKET"; then
	echo "S3 bucket is mounted."
else 
	echo "Mounting the bucket..."
	"$GOPATH/bin/goofys" --endpoint "s3.eu-central-1.amazonaws.com" --sse testbcbio "/mnt/$BUCKET"
	if checkmounted "$BUCKET"; then
		echo "Successfully mounted S3 bucket."
	else
		echo "Something went wrong with the mount. Exiting now."
		exit
	fi
fi

BCBIO_DIR=$(dirname "$(readlinkf `which bcbio_nextgen.py`)")

unset PYTHONHOME
unset PYTHONPATH
export PYTHONNOUSERSITE=1

echo "$BCBIO_DIR/nosetests" -v -s -a $ATTR "$@"
echo
#"$BCBIO_DIR/nosetests" -v -s -a $ATTR "$@"
