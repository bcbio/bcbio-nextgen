#!/bin/bash
set -eo

wget -O NA12878-NGv3-LAB1360-A_1.fastq.gz https://s3.amazonaws.com/bcbio_nextgen/NA12878-NGv3-LAB1360-A_1.fastq.gz
wget -O NA12878-NGv3-LAB1360-A_2.fastq.gz https://s3.amazonaws.com/bcbio_nextgen/NA12878-NGv3-LAB1360-A_2.fastq.gz
wget -O NGv3.bed.gz https://s3.amazonaws.com/bcbio_nextgen/NGv3.bed.gz
wget -O GiaB_NIST_RTG_v0_2.vcf.gz ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/variant_calls/GIAB_integration/NIST_RTG_PlatGen_merged_highconfidence_v0.2_Allannotate.vcf.gz
tabix -f -p vcf GiaB_NIST_RTG_v0_2.vcf.gz
wget -O GiaB_NIST_RTG_v0_2_regions.bed.gz ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/variant_calls/GIAB_integration/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.19_2mindatasets_5minYesNoRatio_AddRTGPlatGenConf_filtNISTclustergt9_RemNISTfilt_RemPartComp_RemRep_RemPartComp_v0.2.bed.gz
gunzip *.bed.gz
