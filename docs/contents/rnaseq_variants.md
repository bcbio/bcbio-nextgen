# Variant calling using bulk RNA-seq data

## Workflow

This workflow demonstrates how to call variants with GATK3.8 using RNA-seq data of RNA-Seq of GM12878:
[SRR307898](https://www.ncbi.nlm.nih.gov/sra/?term=SRR307898). (SRR307897 is of bad quality - don't use it).
The dataset is a bit old (Illumina GAII) but it was used in [Piskol2013](https://www.ncbi.nlm.nih.gov/pubmed/24075185) article,
which is a reliable work on RNA-seq variant calling validation.
GATK3.8 requires an additional installation step:
https://bcbio-nextgen.readthedocs.io/en/latest/contents/installation.html#gatk-and-mutect-mutect2

### 1. Project structure
```
mkdir NA12878_RNA-seq_validation
cd NA12878_RNA-seq_validation
mkdir config input final work
```
### 2. Download input data (~6G total)
```
cd input
wget -c -O NA12878_1.fq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR307/SRR307898/SRR307898_1.fastq.gz
wget -c -O NA12878_2.fq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR307/SRR307898/SRR307898_2.fastq.gz
```

### 3. Config file `input/NA12878.yaml`
```
details:
- algorithm:
    aligner: star
    strandedness: unstranded
    variantcaller: gatk-haplotype
    tools_off:
    - gatk4
  analysis: RNA-seq
  description: NA12878_SRR307898
  files:
  - /path/NA12878/input/NA12878_SRR307898_1.fq.gz
  - /path/NA12878/input/NA12878_SRR307898_2.fq.gz
  genome_build: GRCh37
  metadata:
    batch: NA12878
fc_name: NA12878
resources:
  default:
    cores: 2
    jvm_opts:
    - -Xms750m
    - -Xmx7000m
    memory: 15G
upload:
  dir: ../final
```

### 4. Run bcbio
```
cd work
bcbio_nextgen.py ../config/NA12878.yaml -n 4
```

## Validation (Grch37)
Filter variants passed filters and filter out potential RNA editing events
```
bcftools view -f PASS -e "INFO/DP<10 | INFO/possible_rnaedit==1" NA12878-gatk-haplotype-annotated.vcf.gz |  bgzip -c > NA12878.pass.vcf.gz
tabix NA12878.pass.vcf.gz
wget https://raw.githubusercontent.com/naumenko-sa/cre/master/data/intersect.bed
wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/NA12878.validate.sh
NA12878.validate.sh \
NA12878.pass.vcf.gz \
/path/bcbio/genomes/Hsapiens/GRCh37/validation/giab-NA12878/truth_small_variants.vcf.gz \
intersect.bed \
/path/bcbio/genomes/Hsapiens/GRCh37/rtg/GRCh37.sdf
```

**date**|**type**|**bcbio**|**gatk**|**TP**|**FP**|**FN**|**FDR**|**FNR**|**Target**|**Total\_called**
:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:
2016-12-08|SNP|0.9.9 or 1.0.0|3.6|7,892|266|30,851|3%|80%|38,743|8,158
2016-12-08|INDEL|0.9.9 or 1.0.0|3.6|256|117|2,788|31%|92%|3,044|373
2018-06-06|SNP|1.0.9a0|4.0.1.2|6,817|2,575|31,926|27%|82%|38,743|9,392
2018-06-06|INDEL|1.0.9a0|4.0.1.2|228|24,448|2,816|99%|93%|3,044|24,676
2018-06-14|SNP|1.0.9a0|3.8|7,936|315|30,807|4%|80%|38,743|8,251
2018-06-14|INDEL|1.0.9a0|3.8|306|280|2,738|48%|90%|3,044|586
2019-03-01|SNP|1.1.3|4.1.0.0|6,380|842|32,363|12%|84%|38,743|7,222
2019-03-01|INDEL|1.1.3|4.1.0.0|374|3,131|2,670|89%|88%|3,044|3,505
2019-03-01|SNP|1.1.3|3.8|6,244|558|32,499|8%|84%|38,743|6,802
2019-03-01|INDEL|1.1.3|3.8|192|1,951|2,852|91%|94%|3,044|2,143

## References
- [Piskol2013](https://www.ncbi.nlm.nih.gov/pubmed/24075185)
- [GATK best practices for RNA-seq data variant calling](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-)
- [Validation from Brian Haas using better data](https://github.com/NCIP/ctat-mutations/wiki/Performance-Assessment)
