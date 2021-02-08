# PureCN analysis of tumor-only samples

PureCN by [Markus Riester](https://github.com/lima1/) is a reliable tool to estimate
purity and ploidy in cancer samples sequences with high coverage gene panels and exomes.
It provides also segmentation, gene-level copy numbers, LOH, and classification
of variants.

The important aspect of running PureCN is a requirement to have process matched
normals (as opposed to matched T/N).


## Normal DB (panel of normals)
This step could be performed one time for every panel/cohort, and then a PON
could be re-used to analyze many samples.

### 1. Create purecn_pon_template.yaml:
```yaml
details:
  - analysis: variant2
    genome_build: hg38
    algorithm:
      aligner: false
      svcaller: purecn
      variant_regions: /path/to/panel.bed
      variantcaller: mutect2
```

*If you are using bam files as input (aligner: false) make sure that bai indices are located in
the same dir as input bam files. Bcbio generates indices on the fly, but PureCN follows a symlink
to the original bam file and if there is no index, it crashes.*

### 2. Create a sample sheet pon.csv:
You need a minimum of 3 samples for a PON.
See the discussion about the number of sample in PureCN documentation.
The optimal number could be 30 samples (the more the better).
PureCN chooses which samples to include in the PON when provided
many samples. `batch=pon_build` switches PON mode
(which includes creating SNV PON with mutect2 matching PureCN requirements
and PureCN Normal db).

```
samplename,description,batch,phenotype
sample1__N_FFPE,sample1__N_FFPE,pon_build,normal
sample2__N_FFPE,sample2__N_FFPE,pon_build,normal
sample3__N_FFPE,sample3__N_FFPE,pon_build,normal
```

### 3. Create bcbio project structure and a bcbio config pon.yaml:
```bash
$ ls
# place sample sheet and template in the current dir
purecn_pon_template.yaml pon.csv
# create bcbio project structure
$ mkdir -p pon/input pon/config
# copy or symlink bam or fastq files to pon/input
$ ls pon/input
sample1__N_FFPE.bam sample2__N_FFPE.bam sample3__N_FFPE.bam
# create pon.yaml
$ bcbio_nextgen.py -w template purecn_pon_template.yaml pon.csv pon/input/*.bam
# run bcbio (or create a batch script for you job submission system)
$ cd pon/work
$ bcbio_nextgen.py ../config/pon.yaml -n 20
```

### 4. Collect and save Normal DB and SNV PON
Upon successfull bcbio run you may find resulting files in  `pon/final/[date]_pon`
```bash
$ ls -1 pon/final/*_pon
mapping_bias_hg38.rds
normalDB_hg38.rds
# SNV panel of normals
pon_build-mutect2-annotated.vcf.gz
pon_build-mutect2-annotated.vcf.gz.tbi
```

Also, in the folders of individual samples you may find purecn coverage files:
```bash
$ ls -1 pon/final/sample1__N_FFPE/purecn
```

## Tumor-only analysis

### 1. Create project structure
```bash
mkdir -p ton/config ton/input
# copy/link PON files:
cp panel.bed pon_build-mutect2-annotated.vcf.gz normalDB_hg38.rds mapping_bias_hg38.rds ton/config
# copy, move or symlink input files for tumor samples
mv sample1_T_FFPE.bam sample2_T_FFPE ton/input/
```

### 2. Create purecn_ton_template.yaml:
```bash
wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/templates/purecn_ton.yaml
```

```yaml
details:
  - analysis: variant2
    genome_build: hg38
    algorithm:
      aligner: false
      svcaller: purecn
      variant_regions: /path/ton/config/panel.bed
      variantcaller: mutect2
      background:
        variant: /path/ton/config/pon_build-mutect2-annotated.vcf.gz
        cnv_reference:
          purecn_normaldb: /path/ton/config/normalDB_hg38.rds
          purecn_mapping_bias: /path/ton/config/mapping_bias_hg38.rds
    metadata:
        phenotype: tumor
```

### 3. Create a sample sheet ton.csv:
```
samplename,description,batch,phenotype
sample1_T_FFPE,sample1_T_FFPE,sample1_T_FFPE-batch,tumor
sample2_T_FFPE,sample2_T_FFPE,sample2_T_FFPE-batch,tumor
sample3_T_FFPE,sample3_T_FFPE,sample3_T_FFPE-batch,tumor
sample4_T_FFPE,sample4_T_FFPE,sample4_T_FFPE-batch,tumor
sample5_T_FFPE,sample5_T_FFPE,sample5_T_FFPE-batch,tumor
sample6_T_FFPE,sample6_T_FFPE,sample6_T_FFPE-batch,tumor
sample7_T_FFPE,sample7_T_FFPE,sample7_T_FFPE-batch,tumor
sample8_T_FFPE,sample8_T_FFPE,sample8_T_FFPE-batch,tumor
```

### 4. Create ton.yaml:
```bash
$ ls -1
ton
purecn_ton_template.yaml
ton.csv
$ bcbio_nextgen.py -w template purecn_ton_template.yaml ton/input/*.bam
```

### 5. Run bcbio
```bash
$ cd ton/work
$ bcbio_nextgen.py ../config/ton.yaml -n 20
```

### 6. Collect PureCN results
Upon successfull bcbio run PureCN results are saved in the folders of individual
samples
```
$ ls -1 ton/final/sample1_T_FFPE/purecn
sample1_T_FFPE_amplification_pvalues.csv
sample1_T_FFPE_chromosomes.pdf
sample1_T_FFPE.csv
sample1_T_FFPE_dnacopy.seg
sample1_T_FFPE_genes.csv
sample1_T_FFPE_local_optima.pdf
sample1_T_FFPE.log
sample1_T_FFPE_loh.csv
sample1_T_FFPE.pdf
sample1_T_FFPE.rds
sample1_T_FFPE-ready_coverage_loess.png
sample1_T_FFPE-ready_coverage_loess_qc.txt
sample1_T_FFPE-ready_coverage_loess.txt.gz
sample1_T_FFPE-ready_coverage.txt.gz
sample1_T_FFPE_segmentation.pdf
sample1_T_FFPE_variants.csv
```

## T/N analysis

It is possible to run this analysis with T/N pairs (with PON and normal DB),
but the main recommended mode is tumor-only.

For T/N analysis use purecn_ton_template.yaml (with Normal db and SNV PON)
and create a sample sheet similar to:
```
samplename,description,batch,phenotype
sample1_N_FFPE,sample1_N_FFPE,sample1_T_FFPE-batch,normal
sample1_T_FFPE,sample1_T_FFPE,sample1_T_FFPE-batch,tumor
```

PureCN results will be copied to the folder of the tumor sample in the final dir.

## Troubleshooting

We found useful to adjust memory parameters when running large cohorts with ipython (149 samples): [see issue 3230](https://github.com/bcbio/bcbio-nextgen/issues/3230).

To detect deletions in cfDNA samples we found useful to set PureCN parameters in the yaml:
```yaml
resources:
  purecn:
    options: ["--funsegmentation", "PSCBS", "--minpurity", "0.1", "--minaf", "0.01", "--error", "0.0005"]
```

Capturing more SNP markers is also useful for PureCN analysis, the input bed file for the panel or exome capture kit
is usually 100 bp padded on both sides of the probe (it is not padded by bcbio). In addition to that the mutect2
step uses 50bp interval_padding [option](https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/variation/mutect2.py#L126).
The interval padding could be adjusted with:

```yaml
resources:
  mutect2:
    options: ["interval_padding", "100"]
```

## References

- [PureCN publication](https://scfbm.biomedcentral.com/articles/10.1186/s13029-016-0060-z)
- [PureCN github](https://github.com/lima1/PureCN)
- [PureCN in Bioconductor](https://bioconductor.org/packages/release/bioc/html/PureCN.html)
- [PureCN validation](https://ascopubs.org/doi/full/10.1200/CCI.19.00130)
