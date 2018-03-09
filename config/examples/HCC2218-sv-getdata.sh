#!/bin/bash
# 
# Prepare a run directory for somatic SV validations with the HCC2218 breast
# cancer cell line.
set -eu -o pipefail

mkdir -p config
cd config
wget -c https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/HCC2218-sv.yaml
cd ../

mkdir -p input
cd input
wget -c --no-check-certificate https://s3.amazonaws.com/bcbio/HCC2218/HCC2218Truth-clean-prep.vcf.gz
wget -c --no-check-certificate https://s3.amazonaws.com/bcbio/HCC2218/HCC2218Truth-clean-prep.vcf.gz.tbi
wget -c --no-check-certificate https://s3.amazonaws.com/bcbio/HCC2218/NexteraRapidCapture_Exome_TargetedRegions_v1.2Used.bed.gz
wget -c --no-check-certificate https://s3.amazonaws.com/bcbio/HCC2218/HCC2218BL_S1.bam
wget -c --no-check-certificate https://s3.amazonaws.com/bcbio/HCC2218/HCC2218BL_S1.bam.bai
wget -c --no-check-certificate https://s3.amazonaws.com/bcbio/HCC2218/HCC2218C_S1.bam
wget -c --no-check-certificate https://s3.amazonaws.com/bcbio/HCC2218/HCC2218C_S1.bam.bai
cd ../

mkdir -p work
