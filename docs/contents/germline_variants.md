# Small germline variants

## Workflow1: validate hg38 calls

This workflow validates variant calls using WES data for NA12878 sample.

### 1. Create project structure
```shell
mkdir validate_giab
cd validate_giab
mkdir config input final work
```

### 2. Prepare input data
```shell
cd input
wget -c ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget -c ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
wget -c ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L002_R1_001.fastq.gz
wget -c ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L002_R2_001.fastq.gz
wget -c ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7086_CGTACTAG_L001_R1_001.fastq.gz
wget -c ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7086_CGTACTAG_L001_R2_001.fastq.gz
wget -c ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7086_CGTACTAG_L002_R1_001.fastq.gz
wget -c ftp://ftp-trace.ncbi.nih.gov/ReferenceSamples/giab/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7086_CGTACTAG_L002_R2_001.fastq.gz
# hg19
wget -c ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/nexterarapidcapture_expandedexome_targetedregions.bed.gz
gunzip nexterarapidcapture_expandedexome_targetedregions.bed.gz
cat *R1* > NA12878_1.fq.gz
cat *R2* > NA12878_2.fq.gz
rm Garvan_NA12878*
```

### 3. Convert capture regions file to hg38 coordinates with [UCSC liftover](https://genome.ucsc.edu/cgi-bin/hgLiftOver).
```
Successfully converted 200993 records: View Conversions
Conversion failed on 78 records.
```
Save file to `validate_giab/input/nexterarapidcapture_expandedexome_targetedregions.hg38.bed`

### 4. Set analysis parameters in `config/NA12878.yaml`:
```yaml
details:
  - files:
    - /full/path/validate_giab/input/NA12878_1.fq.gz
    - /full/path/validate_giab/input/NA12878_2.fq.gz
    description: NA12878
    metadata:
      sex: female
    analysis: variant2
    genome_build: hg38
    algorithm:
      aligner: bwa
      variantcaller: gatk-haplotype
      validate: giab-NA12878/truth_small_variants.vcf.gz
      validate_regions: giab-NA12878/truth_regions.bed
      variant_regions: /full/path/validate_giab/input/nexterarapidcapture_expandedexome_targetedregions.hg38.bed
resources:
  default:
    cores: 7
    jvm_opts:
    - -Xms750m
    - -Xmx7000m
    memory: 8G
upload:
  dir: ../final
```

### 5. Run the project (8cores/64G RAM)
```bash
cd validate_giab/work
bcbio_nextgen.py ../config/NA12878.yaml -n 8
```
Running time ~ 2.4h

### 6. Review results:
- `final/[date]_project/grading-summary-NA12878.csv`:
  ```
  sample,caller,vtype,metric,value
  NA12878,gatk-haplotype,SNPs,tp,36710
  NA12878,gatk-haplotype,Indels,tp,3588
  NA12878,gatk-haplotype,SNPs,fp,285
  NA12878,gatk-haplotype,Indels,fp,981
  NA12878,gatk-haplotype,SNPs,fn,323
  NA12878,gatk-haplotype,Indels,fn,611
  ```
- Here FDR_SNPS = 285 / (285 + 36710) = 0.77%, so precision is 99.23%
- FNR_SNPS = 323 / (323 + 36710) = 0.87%, so sensitivity is 99.13%
- more improvement is possible with higher coverage depth and custom filters,
but overall 99.23% precision / 99.13% sensitivity is a good result for out of the box pipeline.
