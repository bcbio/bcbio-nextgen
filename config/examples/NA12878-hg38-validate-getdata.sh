#!/bin/bash
#
# Retrieve data for hg38/GRCh38 validation
# https://bcbio-nextgen.readthedocs.org/en/latest/contents/testing.html#example-pipelines

set -eu -o pipefail

mkdir -p config
cd config
wget -c https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/NA12878-hg38-validate.yaml
cd ..

mkdir -p input
cd input
# Genome data
wget -c -O NA12878_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz
wget -c -O NA12878_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_2.fastq.gz

mkdir -p work
