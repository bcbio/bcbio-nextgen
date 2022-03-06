# Variant calling using bulk RNA-seq data

**2020-05-12: works for validation samples, fails for some random samples due to pybedtools [issue](https://github.com/bcbio/bcbio-nextgen/issues/3078).**
To avoid HaplotypeCallerSpark issue, update to the latest development version.

## Workflow
This workflow demonstrates how to call variants with GATK3.8 using bulk RNA-seq data of GM12878 cell line (blood).

At least 3 RNA-seq datasets exist for NA12878:
- `SRR307897` is of bad quality - don't use it;
- [SRR307898](https://www.ncbi.nlm.nih.gov/sra/?term=SRR307898) is a bit old (Illumina GAII) but it was used in [Piskol2013](https://www.ncbi.nlm.nih.gov/pubmed/24075185) article, which is reliable work on RNA-seq variant calling validation;
- [SRR5665260](https://www.ncbi.nlm.nih.gov/sra/?term=SRR5665260), NextSeq-500.

GATK3.8 requires additional installation [step](https://bcbio-nextgen.readthedocs.io/en/latest/contents/installation.html#gatk-and-mutect-mutect2).

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
or (if ebi is not online)
```
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-1/SRR307898/SRR307898.3
fastq-dump --gzip --split-files SRR307898.3
```

### 3. Config file `input/NA12878.yaml`
Note the use of `batch` even for a single sample.
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
  genome_build: hg38
  metadata:
    batch: NA12878
fc_name: NA12878
resources:
  default:
    cores: 4
    jvm_opts:
    - -Xms750m
    - -Xmx7000m
    memory: 15G
upload:
  dir: ../final
```

### 4. Run bcbio - assuming 60G/4cores machine or cluster node.
```
cd work
bcbio_nextgen.py ../config/NA12878.yaml -n 4
```
## Parameters
- `variantcaller`: gatk-haplotype or vardict. You can use just one variant caller for RNA-seq data in a bcbio project. If you want calls from two callers, run a separate project or edit variantcaller parameter and re-run.
- `tools_off: [gatk4]`: when set, runs gatk3.8, which gives better precision.
- `batch` is required:
   ```
   metadata:
     batch: NA12878
   ```  

## Validation

### How to validate calls from bcbio

Use high quality variants (PASS filters, depth>=10 reads), remove potential RNA editing events.
```
bcftools view -f PASS -e "INFO/DP<10 | INFO/possible_rnaedit==1" NA12878-gatk-haplotype-annotated.vcf.gz |  bgzip -c > NA12878.pass.vcf.gz
tabix NA12878.pass.vcf.gz
wget https://raw.githubusercontent.com/naumenko-sa/cre/master/data/intersect.hg38.bed
wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/NA12878.validate.sh
NA12878.validate.sh \
NA12878.pass.vcf.gz \
/path/bcbio/genomes/Hsapiens/hg38/validation/giab-NA12878/truth_small_variants.vcf.gz \
intersect.hg38.bed \
/path/bcbio/genomes/Hsapiens/hg38/rtg/hg38.sdf
```

Results will be in NA12878.pass.vcf.gz.stat:
```
snp tp-baseline 5720
indels tp-baseline 184
snp fp 1205
indels fp 2141
snp fn 29323
indels fn 2455
```

### Validation results, 2016-2019, GRCh37, SRR307898
```eval_rst
+----------+-----+-----+-------+-----+------+------+---+---+------+------------+
|date      |type |bcbio|gatk   |TP   |FP    |FN    |FDR|FNR|Target|Total called|
+==========+=====+=====+=======+=====+======+======+===+===+======+============+
|2016-12-08|SNP  |1.0.0|3.6    |7,892|266   |30,851|3% |80%|38,743|8,158       |
+----------+-----+-----+-------+-----+------+------+---+---+------+------------+
|2016-12-08|INDEL|1.0.0|3.6    |256  |117   |2,788 |31%|92%|3,044 |373         |
+----------+-----+-----+-------+-----+------+------+---+---+------+------------+
|2018-06-06|SNP  |1.0.9|4.0.1.2|6,817|2,575 |31,926|27%|82%|38,743|9,392       |
+----------+-----+-----+-------+-----+------+------+---+---+------+------------+
|2018-06-06|INDEL|1.0.9|4.0.1.2|228  |24,448|2,816 |99%|93%|3,044 |24,676      |
+----------+-----+-----+-------+-----+------+------+---+---+------+------------+
|2018-06-14|SNP  |1.0.9|3.8    |7,936|315   |30,807|4% |80%|38,743|8,251       |
+----------+-----+-----+-------+-----+------+------+---+---+------+------------+
|2018-06-14|INDEL|1.0.9|3.8    |306  |280   |2,738 |48%|90%|3,044 |586         |
+----------+-----+-----+-------+-----+------+------+---+---+------+------------+
|2019-03-01|SNP  |1.1.3|4.1.0.0|6,380|842   |32,363|12%|84%|38,743|7,222       |
+----------+-----+-----+-------+-----+------+------+---+---+------+------------+
|2019-03-01|INDEL|1.1.3|4.1.0.0|374  |3,131 |2,670 |89%|88%|3,044 |3,505       |
+----------+-----+-----+-------+-----+------+------+---+---+------+------------+
|2019-03-01|SNP  |1.1.3|3.8    |6,244|558   |32,499|8% |84%|38,743|6,802       |
+----------+-----+-----+-------+-----+------+------+---+---+------+------------+
|2019-03-01|INDEL|1.1.3|3.8    |192  |1,951 |2,852 |91%|94%|3,044 |2,143       |
+----------+-----+-----+-------+-----+------+------+---+---+------+------------+
```

### Validation results, 2020, hg38
```eval_rst
+----------+---------+-----+-----+------------------+-----+-----+------+---+---+------+------------+
|date      |data     |type |bcbio|caller            |TP   |FP   |FN    |FDR|FNR|Target|Total called|
+==========+=========+=====+=====+==================+=====+=====+======+===+===+======+============+
|2020-05-12|SRR307898|SNP  |1.2.3|gatk,4.1.6.0      |5,849|1,570|29,194|21%|83%|35,043|7,419       |
+----------+---------+-----+-----+------------------+-----+-----+------+---+---+------+------------+
|2020-05-12|SRR307898|INDEL|1.2.3|gatk,4.1.6.0      |181  |2907 |2458  |94%|93%|2,639 |3,088       |
+----------+---------+-----+-----+------------------+-----+-----+------+---+---+------+------------+
|2020-05-12|SRR307898|SNP  |1.2.3|vardict-java,1.7.0|6,333|916  |28,710|13%|82%|35,043|7,249       |
+----------+---------+-----+-----+------------------+-----+-----+------+---+---+------+------------+
|2020-05-12|SRR307898|INDEL|1.2.3|vardict-java,1.7.0|166  |580  |2,473 |78%|94%|2,639 |746         |
+----------+---------+-----+-----+------------------+-----+-----+------+---+---+------+------------+
|2020-05-13|SRR307898|SNP  |1.2.3|gatk3.8-1.0       |5720 |1205 |29,323|17%|83%|35043 |6925        |
+----------+---------+-----+-----+------------------+-----+-----+------+---+---+------+------------+
|2020-05-13|SRR307898|INDEL|1.2.3|gatk3.8-1.0       |184  |2,141|2,455 |92%|93%|2,639 |2325        |
+----------+---------+-----+-----+------------------+-----+-----+------+---+---+------+------------+
```

## Conclusions
- It is not surprising to see high False Negative rate FN = FN/(FN+TP) in RNA-seq variant calling, i.e. that we are not calling 83% of variants. We validate against intersect.hg38.bed, which is based on exome capture regions and is representing all protein coding genes. Only a minor fraction of genes are expressed in blood, only these genes have read coverage that makes variant calling possible.
- Indel calling is not reliable with RNA-seq data.
- Current recommendation is to use gatk3.8 > vardict > gatk4 for better SNP calling precision in bcbio.
- Annotation of RNA-editing sites and removal of variants around splice junctions improved precision in bcbio, but still
more work is needed to achieve better precision at the level of Piskol2013 (<1% FDR).

## TODO
- update gatk3.8 validation
- validate with SRR5665260
- integrate validation into bcbio
- improve post GATK4 filters for better precision
- paired DNA/RNA-seq variant calling
- somatic variant calling with RNA-seq, see [discusion](https://github.com/bcbio/bcbio-nextgen/issues/3023)
- track pybedtools [release](https://github.com/bcbio/bcbio-nextgen/issues/3078).

## References
- [Piskol2013](https://www.ncbi.nlm.nih.gov/pubmed/24075185)
- [GATK best practices for RNA-seq data variant calling](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-)
- [Variant calling and validation from Brian Haas using CTAT and SRR5665260](https://github.com/NCIP/ctat-mutations/wiki/Performance-Assessment)
- [RTG vcfeval](https://cdn.rawgit.com/RealTimeGenomics/rtg-tools/master/installer/resources/tools/RTGOperationsManual/rtg_command_reference.html#vcfeval)
