#!/bin/bash
set -eu -o pipefail

# We need about 100Gb for the input files. Confirm we have the space.
REQ_DISK_SPACE=100  
df --block-size=G --output='avail' . | sed s/G//g | awk -v req_disk_space=${REQ_DISK_SPACE} '{ if ($1 !~ /Avail/ && $1 < req_disk_space ) printf("Warning: Not enough disk space.\n Warning: Requires %sGb in total but only has %sGb.\n", req_disk_space, $1)  }'

mkdir -p seqc-test/input
cd seqc-test
wget --no-check-certificate https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/rnaseq-seqc.yaml

cd input
for SAMPLE in SRR950078 SRR950079 SRR950080 SRR950081 SRR950082 SRR950083
do 
   wget -O ${SAMPLE}_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR950/${SAMPLE}/${SAMPLE}_1.fastq.gz
   wget -O ${SAMPLE}_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR950/${SAMPLE}/${SAMPLE}_2.fastq.gz
done 

wget --no-check-certificate https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/seqc.csv

cd ../
bcbio_nextgen.py -w template rnaseq-seqc.yaml input/seqc.csv input
