#!/bin/bash

set -e

wget --no-verbose https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/scripts/vm/Vagrantfile
wget --no-verbose https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/scripts/vm/bootstrap.sh

vagrant up

echo "Your VM is up and running. You can complete bcbio setup by doing the following:

vagrant ssh
python /opt/bcbio_nextgen_install.py ~/local/share/bcbio --tooldir=~/local --nodata

You might have to answer a few questions during the install process, but other than that it should be all set!"
