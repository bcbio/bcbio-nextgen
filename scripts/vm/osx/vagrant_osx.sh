#!/bin/bash
mkdir bcbio
cd bcbio
vagrant box add ubuntu_nextgen http://cloud-images.ubuntu.com/precise/current/precise-server-cloudimg-vagrant-amd64-disk1.box
curl -O http://raw.github.com/chapmanb/bcbio-nextgen/blob/master/scripts/vm/Vagrantfile
curl -O http://raw.github.com/chapmanb/bcbio-nextgen/blob/master/scripts/vm/osx/bootstrap.sh
vagrant up

echo "Your VM is all set up. You can connect to the VM by typing 'vagrant ssh'."
echo "You can complete the setup by doing the following:\n
vagrant ssh
cd tmp
sudo python bcbio_nextgen_install.py /usr/local/share/bcbio-nextgen --tooldir=/usr/local\n"
echo "You might have to answer a few questions during the install process, but other\n"
echo "than that it should be all set!"
