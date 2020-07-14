# Methylation

Whole genome bisulfite sequencing is supported using
the [bismark2](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) pipeline.
It can be turned on by setting `analysis` to `wgbs-seq`.
Right now we support only a few kits and only for paired fastq's. If you want to run a different setup,
please [let us know](https://github.com/bcbio/bcbio-nextgen/issues) and we can add support for it.
Consider this a beta feature.

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
