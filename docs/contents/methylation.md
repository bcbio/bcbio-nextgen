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
```
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
    bismark_threads: 16
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
#SBATCH -c 16		              	    # cores
#SBATCH --mem=100G                  # Memory
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=NONE             # Type of email notification (BEGIN, END, FAIL, ALL)

date
bcbio_nextgen.py ../config/bcbio.yaml -n 16
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
- `kit`: `accelngs`, `nebemseq`, `truseq`

In the `resources` section of the yaml:
- trim_galore trimming options:
```yaml
resources:
  trim_galore:
    options: ["--clip_r1 8", "--clip_r2 8", "--three_prime_clip_r1 8", "--three_prime_clip_r2 8"]
```
- bismark `--non-directional` option (default is directional mode and you don't have to specify it) and threading options:
```yaml
resources:
  bismark:
    options: ["--non_directional"]
    bismark_threads: 16
    bowtie_threads: 2
```

- deduplicate_bismark options:
```yaml
resources:
  deduplicate_bismark:
    options: ["--barcode"]
```

The following configs for the `truseq` kit are equivalent:
```
details:
- analysis: wgbs-seq
  genome_build: hg38
  algorithm:
    aligner: bismark
    kit: truseq
```

```
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
    bismark_threads: 16
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
| 07     |4                |2               |13h       | 16       | 50G |
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
