#!/bin/bash
mkdir -p seqc-test seqc-test/input
cd seqc-test
wget --no-check-certificate https://raw.githubusercontent.com/chapmanb/bcbio-nextgen/master/config/examples/rnaseq-seqc.yaml
cd input
wget -O SRR950078_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR950/SRR950078/SRR950078_1.fastq.gz
wget -O SRR950078_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR950/SRR950078/SRR950078_2.fastq.gz
wget -O SRR950079_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR950/SRR950079/SRR950079_1.fastq.gz
wget -O SRR950079_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR950/SRR950079/SRR950079_2.fastq.gz
wget -O SRR950080_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR950/SRR950080/SRR950080_1.fastq.gz
wget -O SRR950080_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR950/SRR950080/SRR950080_2.fastq.gz
wget -O SRR950081_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR950/SRR950081/SRR950081_1.fastq.gz
wget -O SRR950081_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR950/SRR950081/SRR950081_2.fastq.gz
wget -O SRR950082_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR950/SRR950082/SRR950082_1.fastq.gz
wget -O SRR950082_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR950/SRR950082/SRR950082_1.fastq.gz
wget -O SRR950083_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR950/SRR950083/SRR950083_2.fastq.gz
wget -O SRR950083_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR950/SRR950083/SRR950083_2.fastq.gz
# wget -O SRR950084_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR950/SRR950084/SRR950084_2.fastq.gz
# wget -O SRR950084_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR950/SRR950084/SRR950084_2.fastq.gz
# wget -O SRR950085_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR950/SRR950085/SRR950085_2.fastq.gz
# wget -O SRR950085_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR950/SRR950085/SRR950085_2.fastq.gz
# wget -O SRR950086_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR950/SRR950086/SRR950086_2.fastq.gz
# wget -O SRR950086_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR950/SRR950086/SRR950086_2.fastq.gz
# wget -O SRR950087_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR950/SRR950087/SRR950087_2.fastq.gz
# wget -O SRR950087_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR950/SRR950087/SRR950087_2.fastq.gz

wget --no-check-certificate https://raw.githubusercontent.com/chapmanb/bcbio-nextgen/master/config/examples/seqc.csv

cd ../
bcbio_nextgen.py -w template rnaseq-seqc.yaml input/seqc.csv input
