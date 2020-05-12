#!/bin/bash
# validate a vcf file against genome in a bottle calls for NA12878
# $1 - variants.vcf.gz
# #2 - giab file, i.e. bcbio/genomes/Hsapiens/GRCh37/validation/giab-NA12878/truth_small_variants.vcf.gz
# $3 - regions.bed
# $4 - rtg sdf reference, i.e. bcbio/genomes/Hsapiens/GRCh37/rtg/GRCh37.sdf

# rtg manual
# https://github.com/RealTimeGenomics/rtg-tools/blob/master/installer/resources/tools/RTGOperationsManual.pdf

# how bed file is produced in bcbio
# bedtools intersect -nonamecheck -a NA12878-sort-callable_sample.bed -b GiaB_v2_19_regions.bed > NA12878-sort-callable_sample-NA12878-wrm.bed
# uses PASS variants only
# use --vc-score-field='QUAL' for vardict
# memory demanding 10G

# same script works for hg38, use hg38 truth set, bed and SDF

export RTG_JAVA_OPTS='-Xms750m' && export RTG_MEM=9100m && \
   rtg vcfeval --threads 5 \
   -b $2 \
   --bed-regions $3 \
   -c $1 \
   -t $4 \
   -o rtg --vcf-score-field='GQ'
#   --all-records

for f in {tp-baseline,fp,fn}
do
    echo snp $f `bcftools view --types snps rtg/$f.vcf.gz | grep -vc "^#"` >> $1.stat
    echo indels $f `bcftools view --exclude-types snps rtg/$f.vcf.gz | grep -vc "^#"` >> $1.stat
done
