#!/bin/bash
# Germline small and structural variant analysis of NA24385 against Genome in a Bottle references
# Uses 50x HiSeq x10 dataset from 10x genomics:
# https://support.10xgenomics.com/de-novo-assembly/datasets
# http://biorxiv.org/content/early/2016/08/19/070425
set -eu -o pipefail


CORES=4

mkdir -p input
cd input
# Input prep method from raw files
#
# aws s3 cp --source-region us-west-2 s3://10x.files/samples/assembly/msNA24385/msNA24385_fastqs.tar .
# wget -O - https://s3-us-west-2.amazonaws.com/10x.files/samples/assembly/msNA24385/msNA24385_fastqs.tar | tar -xvp
# # combine, de-interleave and bgzip the fastq files
# # https://gist.github.com/nathanhaigh/3521724
# OUT1=NA24385_1.fastq.gz
# OUT2=NA24385_2.fastq.gz
# zcat 25534_FlowCell0/read-RA_si-*.fastq.gz | paste - - - - - - - - | tee >(cut -f 1-4 | tr "\t" "\n" | pbgzip -c -n $CORES > $OUT1) | cut -f 5-8 | tr "\t" "\n" | pbgzip -c -n $CORES > $OUT2
# zcat 25534_FlowCell0/read-I1_si-*.fastq.gz | pbgzip -c -n $CORES > NA24385_index.fastq.gz

# Generated files moved to
# s3://biodata/giab/na24385/
wget -c --no-check-certificate https://s3.amazonaws.com/biodata/giab/na24385/NA24385_1.fastq.gz
wget -c --no-check-certificate https://s3.amazonaws.com/biodata/giab/na24385/NA24385_2.fastq.gz
cd ..

mkdir -p config
cd config
wget -c https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/NA24385-10x.yaml
cd ..

mkdir -p work
