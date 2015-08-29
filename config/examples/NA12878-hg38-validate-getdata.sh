#!/bin/bash
#
# Retrieve data for hg38/GRCh38 validation
# https://bcbio-nextgen.readthedocs.org/en/latest/contents/testing.html#example-pipelines

set -eu -o pipefail

mkdir -p config
cd config
wget -c https://raw.githubusercontent.com/chapmanb/bcbio-nextgen/master/config/examples/NA12878-hg38-validate.yaml
cd ..

mkdir -p input
cd input
# Genome data
wget -c -O NA12878_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz
wget -c -O NA12878_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_2.fastq.gz

# ## Truth sets

# Genome in a Bottle GRCh37
wget -c -O GiaB_v2_19.vcf.gz ftp://ftp.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv2.19/NISTIntegratedCalls_14datasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.19_2mindatasets_5minYesNoRatio_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf.gz
tabix -f -p vcf GiaB_v2_19.vcf.gz
wget -c -O GiaB_v2_19_regions.bed.gz ftp://ftp.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv2.19/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.19_2mindatasets_5minYesNoRatio.bed.gz
gunzip *.bed.gz

# Remapped Genome in a Bottle for hg38
wget -c https://biodata.s3.amazonaws.com/giab_hg38_remap/GiaB_v2_19-38_crossmap.vcf.gz
wget -c https://biodata.s3.amazonaws.com/giab_hg38_remap/GiaB_v2_19-38_crossmap.vcf.gz.tbi
wget -c https://biodata.s3.amazonaws.com/giab_hg38_remap/GiaB_v2_19-38_crossmap-regions.bed
wget -c https://biodata.s3.amazonaws.com/giab_hg38_remap/GiaB_v2_19-38_remap.vcf.gz
wget -c https://biodata.s3.amazonaws.com/giab_hg38_remap/GiaB_v2_19-38_remap.vcf.gz.tbi
wget -c https://biodata.s3.amazonaws.com/giab_hg38_remap/GiaB_v2_19-38_remap-regions.bed

# Platiumn genomes hg19 and hg38
wget -c -O platgene_NA12878-hg19-8_0_1.vcf.gz ftp://platgene_ro:@ussd-ftp.illumina.com/hg19/8.0.1/NA12878/NA12878.vcf.gz
wget -c -O platgene_NA12878-hg19-8_0_1.vcf.gz.tbi ftp://platgene_ro:@ussd-ftp.illumina.com/hg19/8.0.1/NA12878/NA12878.vcf.gz.tbi
wget -c -O platgene_NA12878-hg19-8_0_1-regions.bed.gz ftp://platgene_ro:@ussd-ftp.illumina.com/hg19/8.0.1/NA12878/ConfidentRegions.bed.gz
wget -c -O platgene_NA12878-hg38-2_0_1.vcf.gz ftp://platgene_ro:@ussd-ftp.illumina.com/hg38/2.0.1/NA12878/NA12878.vcf.gz
wget -c -O platgene_NA12878-hg38-2_0_1.vcf.gz.tbi ftp://platgene_ro:@ussd-ftp.illumina.com/hg38/2.0.1/NA12878/NA12878.vcf.gz.tbi
wget -c -O platgene_NA12878-hg38-2_0_1-regions.bed.gz ftp://platgene_ro:@ussd-ftp.illumina.com/hg38/2.0.1/NA12878/ConfidentRegions.bed.gz 
gunzip *.bed.gz
cd ..

mkdir -p work
