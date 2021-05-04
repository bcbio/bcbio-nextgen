# Methylation

Whole genome bisulfite sequencing is supported using
the [bismark2](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) pipeline.
It can be turned on by setting `analysis` to `wgbs-seq`.

## Example run

### 1. (skip if installed) Install bismark and bowtie2 references
```bash
bcbio_nextgen.py upgrade \
-u skip \
--genomes hg38 \
--aligners bowtie2 \
--aligners bismark
```

### 2. Create bcbio project structure and download input fastq files
This example run is based on data from [Encode](https://www.encodeproject.org/experiments/ENCSR890UQO/)

```bash
mkdir wgbs_example
cd wgbs_example
mkdir config input final work
cd input
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR423/008/SRR4235788/SRR4235788_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR423/008/SRR4235788/SRR4235788_2.fastq.gz
```

### 3. Create wgbs_example/config/bcbio.yaml
```yaml
details:
- algorithm:
    aligner: bismark
  analysis: wgbs-seq
  description: NA12878BS
  genome_build: hg38
  files:
  - /path/to/wgbs_example/input/SRR4235788_1.fastq.gz
  - /path/to/wgbs_example/input/SRR4235788_2.fastq.gz
resources:
  trim_galore:
    options: ["--clip_r1 4", "--clip_r2 4", "--three_prime_clip_r1 4", "--three_prime_clip_r2 4"]
  bismark:
    bismark_threads: 4
    bowtie_threads: 2
upload:
  dir: ../final
```

### 4. Create and launch bcbio start script wgbs_example/work/bcbio.sh (O2 slurm example)
```bash
#!/bin/bash

# https://slurm.schedmd.com/sbatch.html
# https://wiki.rc.hms.harvard.edu/display/O2

#SBATCH --partition=priority        # Partition (queue) priority
#SBATCH --time=5-00:00              # Runtime in D-HH:MM format, 10:00:00 for hours
#SBATCH --job-name=wgbs             # Job name
#SBATCH -c 32		              	# cores
#SBATCH --mem=100G                  # Memory
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE            # Type of email notification (BEGIN, END, FAIL, ALL)

date
bcbio_nextgen.py ../config/bcbio.yaml -n 32
date
```

## Parameters

### Supported Kits
```eval_rst
+--------+------------+---------+---------+---------+---------+
|Name    |directional?|TrimR1_5'|TrimR1_3'|TrimR2_5'|TrimR2_3'|
+========+============+=========+=========+=========+=========+
|accelngs|yes         | 10 nt   | 10 nt   | 19 nt   | 5       |
+--------+------------+---------+---------+---------+---------+
|nebemseq|yes         | 5 nt    | 5       | 11      | 5       |
+--------+------------+---------+---------+---------+---------+
|truseq  |no          | 8 nt    | 8       | 8       | 8       |
+--------+------------+---------+---------+---------+---------+
```

In the `algorithm` section of the yaml:
- `aligner`: `bismark`
- `kit`: `accelngs`, `nebemseq`, `truseq`; setting a kit automatically applies the proper trimming options.

In the `resources` section of the yaml:
- trim_galore trimming options:
```yaml
resources:
  trim_galore:
    options: ["--clip_r1 8", "--clip_r2 8", "--three_prime_clip_r1 8", "--three_prime_clip_r2 8"]
```
- bismark `--non-directional` option (default is directional mode and you don't have to specify it) and threading options (use with caution, see benchmarking notes below; by default bcbio is trying to calculate the optimal number of bismark threads with this [function](https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/ngsalign/bismark.py#L123); if you alter bismark_threads = X, request 5-7X cores for bcbio; some samples could be processed much faster, but some could fail, reduce threads then). By defaut, we are using `--local --maxins 1000` alignment options (after benchmarking using nebemseq kit).
```yaml
resources:
  bismark:
    options: ["--non_directional"]
    bismark_threads: 4
    bowtie_threads: 2
```

- deduplicate_bismark options:
```yaml
resources:
  deduplicate_bismark:
    options: ["--barcode"]
```

The following configs for the `truseq` kit are equivalent:
```yaml
details:
- analysis: wgbs-seq
  genome_build: hg38
  algorithm:
    aligner: bismark
    kit: truseq
```

```yaml
details:
- analysis: wgbs-seq
  genome_build: hg38
  algorithm:
    aligner: bismark
resources:
  trim_galore:
    options: ["--clip_r1 8", "--clip_r2 8", "--three_prime_clip_r1 8", "--three_prime_clip_r2 8"]
  bismark:
    options: ["--non_directional"]
```

## Output

### Project directory
- multiqc - QC report, including information from Bismark for all samples (alignment rates, deduplication, M-Bias) 

### Sample directory
- sample-bam_report.txt - bismark alignment report
- sample-deduplication_report.txt
- sample-ready.bam - **sorted** bam  (even though we are using the unsorted bam throughout the pipeline).
- bismark - Bismark output
- bismark/sample.html - Bismark processing report

## Steps

* = a step is repeated for every sample

bcbio.yaml
```yaml
details:
- algorithm:
    aligner: bismark
  analysis: wgbs-seq
  description: rep1
  files:
  - /path/to/ENCSR481JIW_rep1_R1.fastq.gz
  - /path/to/ENCSR481JIW_rep1_R2.fastq.gz
  genome_build: hg38
- algorithm:
    aligner: bismark
  analysis: wgbs-seq
  description: rep2
  files:
  - /path/to/ENCSR481JIW_rep2_R1.fastq.gz
  - /path/to/ENCSR481JIW_rep2_R2.fastq.gz
  genome_build: hg38
upload:
  dir: ../final
```

1. [wgbsseqpipeline function](https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/pipeline/main.py#L413)


2. [Read trimmming](https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/wgbsseq/trimming.py#L15) *
```bash
trim_galore  \
--cores 4 \
--length 30 \
--quality 30 \
--fastqc \
--paired \
-o /path/to/work/bcbiotx \
/path/to/ENCSR481JIW_rep1_R1.fastq.gz \ 
/path/to/ENCSR481JIW_rep1_R2.fastq.gz
```

3. [Bismark align](https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/ngsalign/bismark.py#L19) *
```bash
bismark \
--bowtie2 \
--temp_dir /path/to/work/bcbiotx/ \
--gzip \
--parallel 5 \
-o /path/to/work/bcbiotx \
--unmapped \
/path/to/genomes/Hsapiens/hg38/bismark/ \
-1 /path/to/work/trimmed/rep1/ENCSR481JIW_rep1_R1_val_1.fq.gz \
-2 /path/to/work/trimmed/rep1/ENCSR481JIW_rep1_R2_val_2.fq.gz 
```

4. [deduplicate bismark](https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/wgbsseq/deduplication.py#L10) *
```bash
deduplicate_bismark \
--output_dir /path/to/work/dedup/rep1 \
/path/to/work/align/rep1/rep1.bam
```

5. [bismark_calling](https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/wgbsseq/cpg_caller.py#L58) *
   [bismark_methylation_extractor](https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/wgbsseq/cpg_caller.py#L26) *
```bash
bismark_methylation_extractor \
--no_overlap \
--comprehensive \
--cytosine_report \
--genome_folder /path/to/genomes/Hsapiens/hg38/bismark/ \
--merge_non_CpG \
--multicore 1 \
--buffer_size 5G \
--bedGraph \
--gzip /path/to/work/dedup/rep1/rep1.deduplicated.bam
```

6. [bismark2report](https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/wgbsseq/cpg_caller.py#L44) *
```bash
bismark2report \
--alignment_report /path/to/work/align/rep1/rep1_bismark/ENCSR481JIW_rep1_R1_val_1_bismark_bt2_PE_report.txt \
-o /path/to/work/cpg/rep1/bcbiotx/tmpypivttaa/rep1.html \
--mbias_report /path/to/work/cpg/rep1/rep1.deduplicated.M-bias.txt
```

7. [generate QC metrics](https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/pipeline/qcsummary.py#L39);
[samtools sort](https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/pipeline/qcsummary.py#L67) *
```bash
samtools sort -@ 16 -m 3276M -O BAM  \
-T /path/to/work/bcbiotx/tmpj82rrnlc/rep1.sorted-sort \
-o /path/to/work/bcbiotx/tmpj82rrnlc/rep1.sorted.bam \
/path/to/work/align/rep1/rep1.bam
```

8. samtools index *
```bash
samtools \
index -@ 16 \
/path/to/work/align/rep1/rep1.sorted.bam \
/path/to/work/bcbiotx/tmpr02on2ol/rep1.sorted.bam.bai
```

9. samtools stats *
```bash
samtools stats -@ 16 \
/path/to/work/align/rep1/rep1.sorted.bam > \
/path/to/work/bcbiotx/tmp5e6gerdb/rep1.txt
```

10. downsample for fastqc *
```bash
samtools view -O BAM -@ 16 \
-o /path/to/work/bcbiotx/tmp3z8btzek/rep1.sorted-downsample.bam \
-s 42.735 \
/path/to/work/align/rep1/rep1.sorted.bam
```

11. fastqc *
```bash
fastqc \
-d /path/to/work/qc/rep1/bcbiotx \
-t 16 \
--extract \
-o /path/to/work/qc/rep1/bcbiotx \
-f bam /path/to/work/qc/rep1/rep1.sorted-downsample.bam
```

12. samtools idxstats *
```bash
samtools idxstats /path/to/work/align/rep1/rep1.sorted.bam > \
/path/to/work/bcbiotx/rep1-idxstats.txt
```

13. [multiqc summary](https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/qc/multiqc.py#L39)
```bash
multiqc -c /path/to/work/qc/multiqc/multiqc_config.yaml \
-f -l \
/path/to/work/qc/multiqc/list_files.txt \
-o /path/to/work/bcbiotx/
```

14. [upload sample files to final](https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/upload/__init__.py#L29);
[also see here](https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/upload/__init__.py#L128)

15. [upload project file to final](https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/upload/__init__.py#L21);
[also see here](https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/upload/__init__.py#L777)


## Benchmarking
There is an extensive discussion on Bismark and trim_galore performance, [Bismark github](https://github.com/FelixKrueger/Bismark/issues/96).
We ran a test with NA12878 nebemseq data, 125 mln reads (72.5mln read pairs).
We tested performance of bismark/bcbio using `--parallel` (bismark workers) and `-p` (bowtie threads) bismark settings.
We measured the performance only of the alignment step using bcbio-nextgen-commands log timecodes. 16/2/100G RAM was an optimal parameters set, with other having 5X-10X longer runtimes. When running a cohort of samples ~50% passed with 16/2/100G, some processed broken bam files, re-running with 8/2/100G or 4/2/100G solved the issue, see more info [here](https://github.com/FelixKrueger/Bismark/issues/360). For Lambda Phage genome we re-used trimming step results and 4/2/30G settings.

bcbio.yaml:
```
details:
- algorithm:
    aligner: bismark
    kit: nebemseq
  analysis: wgbs-seq
  description: NA12878_1
  files:
  - /path/to/NA12878_1.fq.gz
  - /path/to/NA12878_2.fq.gz
  genome_build: hg38
  metadata:
    batch: NA12878_1-batch
    phenotype: tumor
fc_name: bcbio
resources:
  bismark:
    bismark_threads: 4
    bowtie_threads: 2
upload:
  dir: ../final
```

Tests:

```eval_rst
+--------+-----------------+----------------+----------+----------+-----+
| test_N | bismark_threads | bowtie_threads | time     | bcbio -n | RAM |
+========+=================+================+==========+==========+=====+
| 01     |1                |2               |14h 42min | 16       | 50G |
+--------+-----------------+----------------+----------+----------+-----+
| 02     |1                |4               |>3.5 days | 16       | 50G |
+--------+-----------------+----------------+----------+----------+-----+
| 03     |1                |8               |1d 21h    | 16       | 50G |
+--------+-----------------+----------------+----------+----------+-----+
| 04     |1                |16              |1d 15h    | 32       | 50G |
+--------+-----------------+----------------+----------+----------+-----+
| 05     |1                |32              |>3day 14h | 64       | 100G|
+--------+-----------------+----------------+----------+----------+-----+
| 06     |2                |2               |2days 4h  | 16       | 50G |
+--------+-----------------+----------------+----------+----------+-----+
| **07** |4                |2               |13h       | 16       | 50G |
+--------+-----------------+----------------+----------+----------+-----+
| 08     |8                |2               |3h 34 min | 16       | 50G |
+--------+-----------------+----------------+----------+----------+-----+
| **09** |16               |2               |2h 42 min | 32       | 100G|
+--------+-----------------+----------------+----------+----------+-----+
| 10     |32               |2               |7h 40 min | 64       | 250G|
+--------+-----------------+----------------+----------+----------+-----+
| 11     |2                |4               |1d 23h    | 16       | 50G |
+--------+-----------------+----------------+----------+----------+-----+
| 12     |2                |8               |11h 30 min| 16       | 50G |
+--------+-----------------+----------------+----------+----------+-----+
| 13     |2                |16              |11h 45 min| 32       | 50G |
+--------+-----------------+----------------+----------+----------+-----+
| 14     |2                |32              |>3days 14h| 64       | 100G|
+--------+-----------------+----------------+----------+----------+-----+
| 15     |4                |2               |4h 4 min  | 16       | 50G |
+--------+-----------------+----------------+----------+----------+-----+
| 16     |4                |4               |15 h      | 16       | 50G |
+--------+-----------------+----------------+----------+----------+-----+
| 17     |4                |8               |6 h       | 32       | 50G |
+--------+-----------------+----------------+----------+----------+-----+
| 18     |4                |16              |15h       | 64       | 100G|
+--------+-----------------+----------------+----------+----------+-----+
| 19     |4                |32              |1d 2h     | 128      | 200G|
+--------+-----------------+----------------+----------+----------+-----+
| 20     |24               |2               |8h        | 48       | 192G|
+--------+-----------------+----------------+----------+----------+-----+
```

ipython parallelization is not implemented for `wgbs-seq`, use 1 node / multicore jobs.

## References
- [Bock.2012.Analysing and interpreting DNA methylation data](https://www.nature.com/articles/nrg3273)
- [accelngs](https://swiftbiosci.com/accel-ngs-methyl-seq-dna-library-kit/)
- [nebemseq](https://www.neb.com/products/e7120-nebnext-enzymatic-methyl-seq-kit)
- [truseq](https://www.illumina.com/products/by-type/sequencing-kits/library-prep-kits/truseq-methyl-capture-epic.html)
