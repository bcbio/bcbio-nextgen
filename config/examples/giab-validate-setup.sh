#!/bin/bash
# 
# Prepare a run directory for evaluations against
# Genome in a Bottle truth sets.
# 
# https://bcbio-nextgen.readthedocs.org/en/latest/contents/testing.html#example-pipelines
#
set -eu -o pipefail

mkdir -p config
cd config
wget -c https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/giab-validate.yaml
cd ..

mkdir -p input
cd input
# Genome data
wget -c -O NA12878_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz
wget -c -O NA12878_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_2.fastq.gz
wget -c -O NA24385_60x.bam ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.60x.1.bam
cd ..

mkdir -p work
