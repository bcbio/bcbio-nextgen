#!/bin/bash

set -e -x

/usr/bin/apt-get -qq update
/usr/bin/apt-get install -q -y build-essential docker.io python-setuptools zlib1g-dev
/usr/bin/apt-get -y autoremove
/usr/bin/apt-get clean

/bin/systemctl start docker.service
/usr/sbin/usermod -aG docker vagrant

NEW_PATH='${HOME}/local/share/bcbio/anaconda/bin:${HOME}/local/bin:${PATH}'
/bin/grep -qxF "PATH=${NEW_PATH}" ~vagrant/.profile || /bin/echo "PATH=${NEW_PATH}" >> ~vagrant/.profile
/bin/grep -qxF 'export BCBIO_TEST_DIR=/tmp/bcbio' ~vagrant/.profile || /bin/echo 'export BCBIO_TEST_DIR=/tmp/bcbio' >> ~vagrant/.profile

/usr/bin/wget --no-verbose -P /opt https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
