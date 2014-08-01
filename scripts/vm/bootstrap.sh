#!/bin/bash
sudo apt-get update
sudo apt-get install -y build-essential
sudo apt-get install -y git
sudo apt-get install -y libz-dev
mkdir tmp
cd tmp
wget https://raw.github.com/chapmanb/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
## python bcbio_nextgen_install.py /usr/local/share/bcbio-nextgen --tooldir=/usr/local --sudo --genomes GRCh37 --aligners bwa
