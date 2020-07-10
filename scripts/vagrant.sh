#!/bin/bash

# Vagrant shell provisioning script https://www.vagrantup.com/docs/provisioning/shell.html

set -e -x

/usr/bin/apt-get -qq update
/usr/bin/apt-get install -q -y build-essential docker.io python-setuptools
/usr/bin/apt-get -y autoremove
/usr/bin/apt-get clean

/bin/systemctl start docker.service
/usr/sbin/usermod -aG docker vagrant

NEW_PATH='${HOME}/local/share/bcbio/anaconda/bin:${HOME}/local/bin:${PATH}'
/bin/grep -qxF "PATH=${NEW_PATH}" ~vagrant/.profile || /bin/echo "PATH=${NEW_PATH}" >> ~vagrant/.profile
# configure test directory within the internal VM file system because synced folders do not support certain features required by bcbio
# for example: mkfifo is required for CWL Cromwell runs, etc
/bin/grep -qxF 'export BCBIO_TEST_DIR=/tmp/bcbio' ~vagrant/.profile || /bin/echo 'export BCBIO_TEST_DIR=/tmp/bcbio' >> ~vagrant/.profile
