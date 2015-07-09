#!/bin/bash
set -eu -o pipefail

wget -O NA12878-NGv3-LAB1360-A_1.fastq.gz https://s3.amazonaws.com/bcbio_nextgen/NA12878-NGv3-LAB1360-A_1.fastq.gz
wget -O NA12878-NGv3-LAB1360-A_2.fastq.gz https://s3.amazonaws.com/bcbio_nextgen/NA12878-NGv3-LAB1360-A_2.fastq.gz
wget -O NGv3.bed.gz https://s3.amazonaws.com/bcbio_nextgen/NGv3.bed.gz
# GiaB calls
wget -O GiaB_v2_19.vcf.gz ftp://ftp.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv2.19/NISTIntegratedCalls_14datasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.19_2mindatasets_5minYesNoRatio_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf.gz
tabix -f -p vcf GiaB_v2_19.vcf.gz
wget -O GiaB_v2_19_regions.bed.gz ftp://ftp.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv2.19/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.19_2mindatasets_5minYesNoRatio.bed.gz
gunzip *.bed.gz
