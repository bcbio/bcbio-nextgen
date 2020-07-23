# Methylation

Whole genome bisulfite sequencing is supported using
the [bismark2](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) pipeline.
It can be turned on by setting `analysis` to `wgbs-seq`.

## Supported Kits
```eval_rst
+--------+------------+---------+---------+---------+---------+
|Name    |directional?|TrimR1_5'|TrimR1_3'|TrimR2_5'|TrimR2_3'|
+========+============+=========+=========+=========+=========+
|accelngs|yes         | 0 nt    | 19 nt   | 19 nt   | 0       |
+--------+------------+---------+---------+---------+---------+
|nebemseq|yes         | 0 nt    | 0       | 0       | 0       |
+--------+------------+---------+---------+---------+---------+
|truseq  |no          | 8 nt    | 8       | 8       | 8       |
+--------+------------+---------+---------+---------+---------+
```

Vendor links:
- [accelngs](https://swiftbiosci.com/accel-ngs-methyl-seq-dna-library-kit/)
- [nebemseq](https://www.neb.com/products/e7120-nebnext-enzymatic-methyl-seq-kit)
- [truseq](https://www.illumina.com/products/by-type/sequencing-kits/library-prep-kits/truseq-methyl-capture-epic.html)


It is possible to specify the `trim_galore` trimming parameters and `bismark` directionality parameter (default: directional) in the `resources` section of the bcbio config file. The following configs for the `truseq` kit will give the same results:
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

## Parameters
- `aligner`: `bismark`
- `kit`: `accelngs`, `nebemseq`, `truseq`

## Benchmarking
There is an extensive discussion on Bismark and trim_galore performance, [Bismark github](https://github.com/FelixKrueger/Bismark/issues/96).
We ran a test with NA12878 nebemseq data, 125 mln reads (72.5mln read pairs).
We tested performance of bismark/bcbio using `--parallel` (bismark workers) and `-p` (bowtie threads) bismark settings.
We measured the performance only of the alignment step using bcbio-nextgen-commands log timecodes. 16/2/100G RAM was an optimal parameters set, with other having 5X-10X longer runtimes. When running a cohort of samples ~50% passed with 16/2/100G, some processed broken bam files, re-running with 8/2/100G or 4/2/100G solved the issue. For Lambda Phage genome we re-used trimming step results and 4/2/30G settings.

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
  - /path/to/NA12878_2fq.gz
  genome_build: hg38
  metadata:
    batch: NA12878_1-batch
    phenotype: tumor
fc_name: bcbio
resources:
  bismark:
    bismark_threads: 1
    bowtie_threads: 8
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
| 06     |2                |2               |2days 4h  |16        | 50G |
+--------+-----------------+----------------+----------+----------+-----+
| 07     |4                |2               |13h       |16        | 50G |
+--------+-----------------+----------------+----------+----------+-----+
| 08     |8                |2               |3h 34 min |16        | 50G |
+--------+-----------------+----------------+----------+----------+-----+
| 09     |16               |2               |2h 42 min |32        | 100G|
+--------+-----------------+----------------+----------+----------+-----+
| 10     |32               |2               |NA        |32        | 100G|
+--------+-----------------+----------------+----------+----------+-----+
| 11     |2                |4               |1d 23h    |16        | 50G |
+--------+-----------------+----------------+----------+----------+-----+
| 12     |2                |8               |11h 30 min|16        | 50G |
+--------+-----------------+----------------+----------+----------+-----+
| 13     |2                |16              |11h 45 min|32        | 50G |
+--------+-----------------+----------------+----------+----------+-----+
| 14     |2                |32              |>3days 14h|64        | 100G|
+--------+-----------------+----------------+----------+----------+-----+
| 15     |4                |2               |4h 4 min  |16        | 50G |
+--------+-----------------+----------------+----------+----------+-----+
| 16     |4                |4               |15 h      |16        | 50G |
+--------+-----------------+----------------+----------+----------+-----+
| 17     |4                |8               |6 h       |32        | 50G |
+--------+-----------------+----------------+----------+----------+-----+
| 18     |4                |16              |15h       |64        | 100G|
+--------+-----------------+----------------+----------+----------+-----+
| 19     |4                |32              |1d 2h     |128       | 200G|
+--------+-----------------+----------------+----------+----------+-----+
```

## References
- [Bock.2012.Analysing and interpreting DNA methylation data](https://www.nature.com/articles/nrg3273)
