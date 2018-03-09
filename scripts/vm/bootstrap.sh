#!/bin/bash
sudo apt-get update
sudo apt-get install -y build-essential zlib1g-dev wget curl python-setuptools git libz-dev
mkdir tmp
cd tmp
wget https://raw.github.com/bcbio/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
## python bcbio_nextgen_install.py /usr/local/share/bcbio-nextgen --tooldir=/usr/local --sudo --genomes GRCh37 --aligners bwa
