#!/bin/bash

set -e -x

/usr/bin/apt-get -qq update
/usr/bin/apt-get install -q -y build-essential python-setuptools zlib1g-dev
/usr/bin/apt-get -y autoremove
/usr/bin/apt-get clean

NEW_PATH='${HOME}/local/share/bcbio/anaconda/bin:${HOME}/local/bin:${PATH}'
grep -qxF "PATH=${NEW_PATH}" ~vagrant/.profile || echo "PATH=${NEW_PATH}" >> ~vagrant/.profile
grep -qxF 'export BCBIO_TEST_DIR=/tmp/bcbio' ~vagrant/.profile || echo 'export BCBIO_TEST_DIR=/tmp/bcbio' >> ~vagrant/.profile

wget --no-verbose -P /opt https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
