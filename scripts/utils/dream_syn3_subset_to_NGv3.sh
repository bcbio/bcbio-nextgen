#!/bin/bash
#SBATCH -t 0-08:00:00
#SBATCH -p regal
#SBATCH --cpus-per-task=32
#SBATCH --mem=128000
#
# Uses a full run of the somatic DREAM synthetic 3 challenge to subset paired
# end read pairs to exome only regions. Not generalized for multiple systems but
# is a useful template for doing this: get IDs of all reads in exome regions from
# alignment, convert the whole BAM file into fastq and subset to only these
# read pairs using seqtk. This works around the problem of not having a clean
# way to extract read pairs from a subset of the genome if either pair maps into
# a region.
#
set -e
sambamba view -L ../input/NGv3.bed ../work_brad/align/syn3-tumor/2_2014-08-13_dream-syn3-sort.bam | cut -f 1 | sort | uniq | /n/regal/hsph_bioinfo/bcbio_nextgen/anaconda/bin/py -x 'x + "/1\n" + x +"/2"' > tumor_NGv3_ids.txt
sambamba view -L ../input/NGv3.bed  ../work_brad/align/syn3-normal/1_2014-08-13_dream-syn3-sort.bam | cut -f 1 | sort | uniq | /n/regal/hsph_bioinfo/bcbio_nextgen/anaconda/bin/py -x 'x + "/1\n" + x +"/2"' > normal_NGv3_ids.txt
bamtofastq collate=1 F=synthetic_challenge_set3_tumor_1.fq F2=synthetic_challenge_set3_tumor_2.fq filename=../work_brad/align/syn3-tumor/2_2014-08-13_dream-syn3-sort.bam
bamtofastq collate=1 F=synthetic_challenge_set3_normal_1.fq F2=synthetic_challenge_set3_normal_2.fq filename=../work_brad/align/syn3-normal/1_2014-08-13_dream-syn3-sort.bam
seqtk subseq synthetic_challenge_set3_tumor_1.fq tumor_NGv3_ids.txt > synthetic_challenge_set3_tumor_NGv3_1.fq
seqtk subseq synthetic_challenge_set3_tumor_2.fq tumor_NGv3_ids.txt > synthetic_challenge_set3_tumor_NGv3_2.fq
bgzip synthetic_challenge_set3_tumor_NGv3_1.fq
bgzip synthetic_challenge_set3_tumor_NGv3_2.fq
seqtk subseq synthetic_challenge_set3_normal_1.fq normal_NGv3_ids.txt > synthetic_challenge_set3_normal_NGv3_1.fq
seqtk subseq synthetic_challenge_set3_normal_2.fq normal_NGv3_ids.txt > synthetic_challenge_set3_normal_NGv3_2.fq
bgzip synthetic_challenge_set3_normal_NGv3_1.fq
bgzip synthetic_challenge_set3_normal_NGv3_2.fq
