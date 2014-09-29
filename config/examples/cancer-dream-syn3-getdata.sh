#!/bin/bash
# Retrieve data for doing cancer calling evaluation using synthetic dataset 3
# from the ICGC-TCGA DREAM challenge:
# https://www.synapse.org/#!Synapse:syn312572/wiki/62018
wget --no-check-certificate https://cghub.ucsc.edu/software/downloads/GeneTorrent/3.8.5a/GeneTorrent-download-3.8.5a-94-CentOS5.8.x86_64.tar.gz
tar -xzvpf GeneTorrent-download-3.8.5a-94-CentOS5.8.x86_64.tar.gz
cghub/bin/gtdownload -v -c http://dream.annailabs.com/dream_public.pem -d https://dream.annailabs.com/cghub/data/analysis/download/b19d76a0-a487-4c50-8f9c-3b4d5e53239d
ln -s b19d76a0-a487-4c50-8f9c-3b4d5e53239d/*.bam
ln -s b19d76a0-a487-4c50-8f9c-3b4d5e53239d/*.bam.bai
cghub/bin/gtdownload -v -c http://dream.annailabs.com/dream_public.pem -d https://dream.annailabs.com/cghub/data/analysis/download/8fe6fc33-2daf-4393-929f-7c3493d04bef
ln -s 8fe6fc33-2daf-4393-929f-7c3493d04bef/*.bam
ln -s 8fe6fc33-2daf-4393-929f-7c3493d04bef/*.bam.bai
wget https://s3.amazonaws.com/bcbio_nextgen/dream/synthetic_challenge_set3_tumor_20pctmasked_truth.tar.gz
tar -xzvpf synthetic_challenge_set3_tumor_20pctmasked_truth.tar.gz
wget https://s3.amazonaws.com/bcbio_nextgen/dream/refseq-merged.bed.gz
gunzip refseq-merged.bed
