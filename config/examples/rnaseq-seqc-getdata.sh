#!/bin/bash
mkdir -p seqc-test/input
cd seqc-test
wget --no-check-certificate https://raw.githubusercontent.com/chapmanb/bcbio-nextgen/master/config/examples/rnaseq-seqc.yaml

# These six input files are about 50Gb in total 
cd input
for SAMPLE in SRR950078 SRR950079 SRR950080 SRR950081 SRR950082 SRR950083
do 
   wget -O ${SAMPLE}_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR950/${SAMPLE}/${SAMPLE}_1.fastq.gz
   wget -O ${SAMPLE}_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR950/${SAMPLE}/${SAMPLE}_2.fastq.gz
done 

wget --no-check-certificate https://raw.githubusercontent.com/chapmanb/bcbio-nextgen/master/config/examples/seqc.csv

cd ../
bcbio_nextgen.py -w template rnaseq-seqc.yaml input/seqc.csv input
