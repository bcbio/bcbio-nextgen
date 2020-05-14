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

## Workflow2: Basic germline calling

The best approach to build a bcbio configuration for germline calling is to use
the automated sample configuration with one of the default templates:

* [FreeBayes template](https://github.com/bcbio/bcbio-nextgen/blob/master/config/templates/freebayes-variant.yaml)
--Call variants using FreeBayes with a minimal preparation pipeline. This is a freely available unrestricted pipeline 
fully included in the bcbio installation.
* [GATK HaplotypeCaller template](https://github.com/bcbio/bcbio-nextgen/blob/master/config/templates/gatk-variant.yaml)
--Run GATK best practices, including Base Quality Score Recalibration, realignment and HaplotypeCaller variant calling. This requires a license from Broad for commercial use. You need to manually install GATK along with bcbio using downloads from the GATK Broad site or Appistry.

## Workflow3: Population calling

When calling multiple samples, we recommend calling together to provide improved
sensitivity and a fully squared off final callset. To associate samples together
in a population add a `metadata` `batch` to the sample configuration:

```yaml
- description: Sample1
  metadata:
    batch: Batch1
- description: Sample2
  metadata:
    batch: Batch1
```
Batching samples results in output VCFs and GEMINI databases containing
all merged sample calls.
bcbio has two methods to call samples together:

* Batch or pooled calling -- This calls all samples simultaneously by feeding
them to the variant caller. This works for smaller batch sizes (< 100 samples)
as memory requirements become limiting in larger pools. This is the default approach
taken when you specify a `variantcaller` in the variant calling configuration.

* Joint calling -- This calls samples independently, then combines them together
into a single callset by integrating the individual calls.
This scales to larger population sizes by avoiding the computational bottlenecks
of pooled calling. We recommend joint calling with HaplotypeCaller
but also support joint calling with FreeBayes using a custom implementation.
Specifying a `jointcaller` along with the appropriate `variantcaller` in the variant calling configuration enables this

```yaml
- description: Sample1
  algorithm:
    variantcaller: gatk-haplotype
    jointcaller: gatk-haplotype-joint
  metadata:
    batch: Batch1
- description: Sample2
  algorithm:
    variantcaller: gatk-haplotype
    jointcaller: gatk-haplotype-joint
  metadata:
    batch: Batch1
```

## Workflow4: Whole genome trio (50x) - hg38

This input configuration runs whole genome bwa alignment and GATK variant calling.
It uses a father/mother/child trio from the
[CEPH NA12878 family](https://blog.goldenhelix.com/wp-content/uploads/2013/03/Utah-Pedigree-1463-with-NA12878.png):
NA12891, NA12892, NA12878. Illumina's [Platinum genomes project](https://www.illumina.com/platinumgenomes.html) has
50X whole genome sequencing of the three members.
The analysis compares results against a reference NA12878 callset from NIST's
[Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) initiative.

To run the analysis do:
```shell
mkdir NA12878-trio-eval
cd NA12878-trio-eval
mkdir config input work
cd config
wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/NA12878-trio-wgs-validate.yaml
cd ../input
wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/NA12878-trio-wgs-validate-getdata.sh
bash NA12878-trio-wgs-validate-getdata.sh
cd ../work
bcbio_nextgen.py ../config/NA12878-trio-wgs-validate.yaml -n 16
```
This is a large whole genome analysis and meant to test both pipeline scaling
and validation across the entire genome. It can take multiple days to run depending on available cores.
It requires 300GB for the input files and 1.3TB for the work directory.
Smaller examples below exercise the pipeline with less disk and computational requirements.
