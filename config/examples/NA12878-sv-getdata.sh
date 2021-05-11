#!/bin/bash
# 
# Prepare a run directory for structural variant evaluations for NA12878 against
# Genome in a Bottle truth sets of deletions and insertions.
# 
# https://bcbio-nextgen.readthedocs.org/en/latest/contents/testing.html#example-pipelines
#
set -eu -o pipefail

mkdir -p config
cd config
wget -c https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/NA12878-sv.yaml
cd ..

mkdir -p input
cd input
# Genome data
wget -c -O NA12878_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz
wget -c -O NA12878_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_2.fastq.gz
# Validation data
wget -O - ftp://ftp.ncbi.nih.gov/giab/ftp/technical/svclassify_Manuscript/Supplementary_Information/Personalis_1000_Genomes_deduplicated_deletions.bed | grep -v ^Chr > giab-svclassify-deletions-2015-05-22.bed
wget -O -  ftp://ftp.ncbi.nih.gov/giab/ftp/technical/svclassify_Manuscript/Supplementary_Information/Spiral_Genetics_insertions.bed | grep -v ^Chr > giab-svclassify-insertions-2015-05-22.bed
cd ..

mkdir -p work
