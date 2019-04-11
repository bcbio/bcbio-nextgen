#!/bin/bash
#
# Cancer-like mixture of two Genome in a Bottle samples (NA12878 and NA24385)
# ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/use_cases/mixtures/UMCUTRECHT_NA12878_NA24385_mixture_10052016/

set -eu -o pipefail

mkdir -p config
cd config
wget -c https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/cancer-giab-na12878-na24385.yaml
cd ..

base_url=ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/use_cases/mixtures/UMCUTRECHT_NA12878_NA24385_mixture_10052016
mkdir -p input
cd input

# Genome data
wget -c $base_url/README-NA12878_NA24385_mixture.txt
wget -c $base_url/24385-12878-30-200_R1_001.fastq.gz
wget -c $base_url/24385-12878-30-200_R2_001.fastq.gz
wget -c $base_url/24385-200_AH5G7WCCXX_S4_L004_R1_001.fastq.gz
wget -c $base_url/24385-200_AH5G7WCCXX_S4_L004_R2_001.fastq.gz

# Truth sets
wget -c $base_url/na12878-na24385-somatic-truth-regions.bed
wget -c $base_url/na12878-na24385-somatic-truth.vcf.gz
wget -c $base_url/na12878-na24385-somatic-truth.vcf.gz.tbi

cd ..
mkdir -p work
