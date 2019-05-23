#!/bin/bash
set -e
wget https://raw.github.com/bcbio/bcbio-nextgen/master/scripts/vm/Vagrantfile
wget https://raw.github.com/bcbio/bcbio-nextgen/master/scripts/vm/bootstrap.sh
vagrant up

echo "Your VM is all set up. You can connect to the VM by typing 'vagrant ssh'."
echo "You can complete the setup by doing the following:
vagrant ssh
cd tmp
python bcbio_nextgen_install.py ~/local/share/bcbio-nextgen --tooldir=~/local --nodata"
echo "You might have to answer a few questions during the install process, but other"
echo "than that it should be all set!"
