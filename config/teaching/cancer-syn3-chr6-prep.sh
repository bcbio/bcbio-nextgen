#!/bin/bash
# Setup directory for running a minimal DREAM synthetic 3 dataset
# The original DREAM data is subset to exomes on chromosome 6
set -eu -o pipefail

mkdir -p cancer-syn3-chr6
cd cancer-syn3-chr6
wget -c https://s3.amazonaws.com/bcbio_nextgen/dream/cancer-syn3-chr6-input.tar.gz
tar -xzvpf cancer-syn3-chr6-input.tar.gz
rm -f cancer-syn3-chr6-input.tar.gz

mkdir -p config
cd config
wget -c https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/teaching/cancer-syn3-chr6.yaml
cd ..

mkdir -p work
