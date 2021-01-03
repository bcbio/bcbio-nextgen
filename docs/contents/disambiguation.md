# Disambiguation

* `disambiguate` For mixed or explant samples, provide a list of `genome_build`
identifiers to check and remove from alignment. Currently supports cleaning a single organism. For example, with `genome_build: hg19`
and `disambiguate: [mm10]`, it will align to hg19 and mm10, run disambiguation and discard
reads confidently aligned to mm10 and not hg19. Affects fusion detection when `star` is chosen as the aligner. 
Aligner must be set to a non false value for this to run.

Example config:
```
details:
- algorithm:
    aligner: bwa
    background: /path/to/project/config/1000g_pon.hg38.vcf.gz
    disambiguate: mm10
    mark_duplicates: true
    platform: illumina
    quality_format: standard
    realign: false
    recalibrate: false
    remove_lcr: true
    tools_on:
    - noalt_calling
    variantcaller:
    - vardict
    - mutect2
    vcfanno: somatic
  analysis: variant2
  description: SAMPLE
  files:
  - /path/to/project/input/sample_1.fq.gz
  - /path/to/project/input/sample_2.fq.gz
  genome_build: hg38
  metadata:
    batch: bSAMPLE
    phenotype: tumor
upload:
  dir: ../final
```

The resulting `final/project/multiqc/multiqc_report.html` will contain hg38 and mm10 columns with numbers of reads aligned to hg38 and mm10.
Also `final/sample` contains sample-ready.bam - reads aligned to hg38 and sample-disambiguate-mm10.bam - reads aligned to mm10.

