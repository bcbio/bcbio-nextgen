#!/bin/bash
wget https://raw.github.com/chapmanb/bcbio-nextgen/master/scripts/vm/Vagrantfile
wget https://raw.github.com/chapmanb/bcbio-nextgen/master/scripts/vm/bootstrap.sh
vagrant up

echo "Your VM is all set up. You can connect to the VM by typing 'vagrant ssh'."
echo "You can complete the setup by doing the following:\n
vagrant ssh
cd tmp
python bcbio_nextgen_install.py /usr/local/share/bcbio-nextgen --tooldir=/usr/local\n"
echo "You might have to answer a few questions during the install process, but other\n"
echo "than that it should be all set!"
