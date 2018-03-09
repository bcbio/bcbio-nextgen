#!/bin/bash
set -eu -o pipefail

# Data retrieval script for validation comparing alignment methods, preparation approaches
# and variant callers for an NA12878 exome dataset from EdgeBio.
#
# See the bcbio-nextgen documentation for full instructions to
# run this analysis:
# https://bcbio-nextgen.readthedocs.org/en/latest/contents/testing.html#example-pipelines

mkdir -p config
cd config
wget https://raw.github.com/bcbio/bcbio-nextgen/master/config/examples/NA12878-exome-methodcmp.yaml
cd ..

mkdir -p input
cd input
wget -c -O NA12878-NGv3-LAB1360-A_1.fastq.gz https://s3.amazonaws.com/bcbio_nextgen/NA12878-NGv3-LAB1360-A_1.fastq.gz
wget -c -O NA12878-NGv3-LAB1360-A_2.fastq.gz https://s3.amazonaws.com/bcbio_nextgen/NA12878-NGv3-LAB1360-A_2.fastq.gz
cd ..

mkdir -p work
