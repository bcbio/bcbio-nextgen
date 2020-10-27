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
This is a large whole genome analysis and meant to test both pipeline scaling and validation across the entire genome. It can take multiple days to run depending on available cores. It requires 300GB for the input files, 1.3TB for the work directory, and 48GB of memory (with 16 cores). Smaller examples below exercise the pipeline with less disk and computational requirements.

## Workflow5: Whole genome (10x)

An input configuration for running whole gnome variant calling with bwa
and GATK, using Illumina's [Platinum genomes project](https://www.illumina.com/platinumgenomes.html)
[NA12878-illumina.yaml](https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/NA12878-illumina.yaml).
See this [blog post on whole genome scaling](https://bcb.io/2013/05/22/scaling-variant-detection-pipelines-for-whole-genome-sequencing-analysis/)
for expected run times and more information about the pipeline. To run the analysis:

* Create an input directory structure like:
    ```shell
    ├── config
    │   └── NA12878-illumina.yaml
    ├── input
    └── work
    ```
* Retrieve inputs and comparison calls:
    ```shell
    cd input
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR091/ERR091571/ERR091571_1.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR091/ERR091571/ERR091571_2.fastq.gz
    ```
* Retrieve configuration input file:
    ```shell
    cd config
    wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/NA12878-illumina.yaml
    ```
* Run analysis on 16 core machine:
    ```shell
    cd work
    bcbio_nextgen.py ../config/NA12878-illumina.yaml -n 16
    ```
* Examine summary of concordance and discordance to comparison calls from the `grading-summary.csv` file in the work directory.

## Parameters
* `variantcaller` Variant calling algorithm. Can be a list of multiple options or
false to skip [false, freebayes, gatk-haplotype, haplotyper, platypus, mutect, mutect2, scalpel,
tnhaplotyper, tnscope, vardict, varscan, samtools, gatk]
  * Paired (typically somatic, tumor-normal) variant calling is currently supported by vardict,
  freebayes, mutect2, mutect (see disclaimer below), scalpel (indels only), tnhaplotyper (Sentieon),
  tnscope (Sentieon) and varscan. See somatic variant calling documentation for details on pairing tumor and normal samples.
  * You can generate both somatic and germline calls for paired tumor-normal samples using different sets of callers.
  The pipeline documentation on calling `Somatic with germline variants` details how to do this.
  * mutect, a SNP-only caller, can be combined with indels from scalpel or sid.
  Mutect operates in both tumor-normal and tumor-only modes. In tumor-only mode the indels from scalpel will reflect all indels in the sample, as there is currently no way of separating the germline from somatic indels in tumor-only mode.
* `indelcaller` For the MuTect SNP only variant caller it is possible to add calls from
an indelcaller such as scalpel, pindel and somatic indel detector (for Appistry MuTect users only).
Currently an experimental option that adds these indel calls to MuTect's SNP-only output.
Only one caller supported. Omit to ignore. [scalpel, pindel, sid, false]
* `jointcaller` Joint calling algorithm, combining variants called with the specified `variantcaller`.
Can be a list of multiple options but needs to match with appropriate `variantcaller`.
Joint calling is only needed for larger input sample sizes (>100 samples),
otherwise use standard pooled `population calling`:
  * `gatk-haplotype-joint` [GATK incremental joint discovery](https://gatkforums.broadinstitute.org/gatk/discussion/3896/the-gatk-reference-model-pipeline-for-incremental-joint-discovery-in-full-detail)
  with HaplotypeCaller. Takes individual gVCFs called by `gatk-haploype` and perform combined genotyping.
  * `freebayes-joint` Combine freebayes calls using [bcbio.variation.recall](https://github.com/chapmanb/bcbio.variation.recall)
  with recalling at all positions found in each individual sample. Requires `freebayes` variant calling.
  * `platypus-joint` Combine platypus calls using bcbio.variation.recall
  with squaring off at all positions found in each individual sample. Requires `platypus` variant calling.
  * `samtools-joint` Combine samtools calls using bcbio.variation.recall
  with squaring off at all positions found in each individual sample. Requires `samtools` variant calling.
* `joint_group_size` Specify the maximum number of gVCF samples to feed into joint calling.
Currently applies to GATK HaplotypeCaller joint calling and defaults to the GATK recommendation of 200.
Larger numbers of samples will first get combined prior to genotyping.
* `ploidy` Ploidy of called reads. Defaults to 2 (diploid). You can also tweak
specialty ploidy like mitochondrial calling by setting ploidy as a dictionary.
The defaults are:
```yaml
ploidy:
  default: 2
  mitochondrial: 1
  female: 2
  male: 1
```
* `background` Provide pre-calculated files to use as backgrounds for different processes.
Organized as a dictionary with individual keys for different components of the pipeline. You can enter as many or few as needed:
  * `variant` A VCF file with variants to use as a background reference during variant calling.
  For tumor/normal paired calling use this to supply a panel of normal individuals.
  * `cnv_reference` Background reference file for copy number calling. This can be either
  a single file for one CNV method or a dictionary for multiple methods.
  Supports [CNVkit cnn inputs](https://cnvkit.readthedocs.io/en/stable/fileformats.html#copy-number-reference-profile-cnn),
  [GATK4 HDF5 panel of normals](https://software.broadinstitute.org/gatk/documentation/article?id=11682)
  and [seq2c](https://github.com/AstraZeneca-NGS/Seq2C) combined mapping plus coverage files:
  ```yaml
  background:
    cnv_reference:
      cnvkit: /path/to/background.cnn
      gatk-cnv: /path/to/background_pon.hdf5
      seq2c: /path/to/background.tsv
  ```

### Variant annotation
* `effects` Method used to calculate expected variant effects;
defaults to [snpEff](http://snpeff.sourceforge.net/).
[Ensembl variant effect predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html)
is also available when downloaded using `Customizing data installation`. [snpeff, vep, false]
* `effects_transcripts` Define the transcripts to use for effect prediction annotation.
Options `all`: Standard Ensembl transcript list (the default); `canonical`: Report single canonical transcripts
(`-canon` in snpEff, `-pick` in VEP); `canonical_cancer` Canonical transcripts with hand curated changes for more common cancer transcripts (effects snpEff only).
* `vcfanno` Configuration files for [vcfanno](https://github.com/brentp/vcfanno),
allowing the application of additional annotations to variant calls. By default, bcbio will try and apply:
  * `gemini` -- External population level annotations from [GEMINI](https://gemini.readthedocs.io).
  This is only run for human samples with gemini data installed.
  * `somatic` -- Somatic annotations from COSMIC, ClinVar and friends.
  COSMIC need a custom installation within bcbio. Only added for tumor or tumor/normal somatic calling.
  * `rnaedit` -- RNA editing sites for RNA-seq variant calling runs.
    bcbio installs pre-prepared configuration files in `genomes/build/config/vcfanno`
    or you can specify the full path to a `/path/your/anns.conf` and optionally
    an equivalently named `/path/your/anns.lua` file. This value can be a list for multiple inputs.

## Output
See description in somatic_variants.

## Validation

```eval_rst
+----------+-------------------+-----+-----+------------+------+---+---+-----+-----+-------+------------+
|date      |data               |type |bcbio|caller      |TP    |FP |FN |FDR  |FNR  |Target |Total called|
+==========+===================+=====+=====+============+======+===+===+=====+=====+=======+============+
|2020-05-14|Garvan_NA12878(WES)|SNP  |1.2.3|gatk,4.1.6.0|36,710|285|323|0.77%|0.87%|37,033 |36,995      |
+----------+-------------------+-----+-----+------------+------+---+---+-----+-----+-------+------------+
|2020-05-14|Garvan_NA12878(WES)|INDEL|1.2.3|gatk,4.1.6.0|3,588 |981|611|21%  |15%  |4,199  |4,569       |
+----------+-------------------+-----+-----+------------+------+---+---+-----+-----+-------+------------+
```

bcbio pre-installs standard truth sets for performing validation, and also allows
use of custom local files for assessing reliability of your runs:
* `validate` A VCF file of expected variant calls to perform validation and
grading of small variants (SNPs and indels) from the pipeline. This provides
a mechanism to ensure consistency of calls against a known set of variants,
supporting comparisons to genotyping array data or reference materials.
* `validate_regions` A BED file of regions to evaluate small variant calls in.
This defines specific regions covered by the `validate` VCF file.
* `svvalidate` -- Dictionary of call types and pointer to BED file of known regions.
For example: `DEL: known_deletions.bed` does deletion based validation of outputs
against the BED file.

Each option can be either the path to a local file, or a partial path to a file
in the pre-installed truth sets. For instance, to validate an NA12878 run against
the [Genome in a Bottle](https://github.com/genome-in-a-bottle) truth set:
```yaml
validate: giab-NA12878/truth_small_variants.vcf.gz
validate_regions: giab-NA12878/truth_regions.bed
svvalidate:
  DEL: giab-NA12878/truth_DEL.bed
```

follow the same naming schemes for small variants, regions and different
structural variant types. bcbio has the following validation materials
for germline validations:
* `giab-NA12878` -- [Genome in a Bottle](https://github.com/genome-in-a-bottle)
for NA12878, a Caucasian sample. Truth sets: small_variants, regions, DEL; Builds: GRCh37, hg19, hg38
* `giab-NA24385` -- [Genome in a Bottle](https://github.com/genome-in-a-bottle)
for NA24385, an Ashkenazic Jewish sample. Truth sets: small_variants, regions; Builds: GRCh37, hg19, hg38
* `giab-NA24631` -- [Genome in a Bottle](https://github.com/genome-in-a-bottle)
for NA24631, a Chinese sample. Truth sets: small_variants, regions; Builds: GRCh37, hg19, hg38
* `giab-NA12878-crossmap` -- [Genome in a Bottle](https://github.com/genome-in-a-bottle)
for NA12878 converted to hg38 with CrossMap. Truth sets: small_variants, regions, DEL; Builds: hg38
* `giab-NA12878-remap` -- [Genome in a Bottle](https://github.com/genome-in-a-bottle)
for NA12878 converted to hg38 with Remap. Truth sets: small_variants, regions, DEL; Builds: hg38
* `platinum-genome-NA12878` -- [Illumina Platinum Genome](https://www.illumina.com/platinumgenomes/)
for NA12878. Truth sets: small_variants, regions; Builds: hg19, hg38

For more information on the hg38 truth set preparation see the work
on [validation on build 38 and conversion of human build 37 truth sets to build 38](https://bcb.io/2015/09/17/hg38-validation/).
The [installation recipes](https://github.com/chapmanb/cloudbiolinux/tree/master/ggd-recipes)
contain provenance details about the origins of the installed files.

## References
- [NIST genome in a bottle](https://www.nist.gov/programs-projects/genome-bottle)
- [GIAB data index](https://github.com/genome-in-a-bottle/giab_data_indexes)
- [Illumina Platinum Genomes](https://www.illumina.com/platinumgenomes.html)
- An introduction to the [variant evaluation framework](https://bcb.io/2013/05/06/framework-for-evaluating-variant-detection-methods-comparison-of-aligners-and-callers/).
This includes a comparison of the [bwa mem](http://bio-bwa.sourceforge.net/)
and [novoalign](http://www.novocraft.com) aligners.
We also compared the [FreeBayes](https://github.com/ekg/freebayes),
[GATK HaplotypeCaller](https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php)
and [GATK UnifiedGenotyper](https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_genotyper_UnifiedGenotyper.php) variant callers.
- An in-depth evaluation of
[FreeBayes and BAM post-alignment processing](https://bcb.io/2013/10/21/updated-comparison-of-variant-detection-methods-ensemble-freebayes-and-minimal-bam-preparation-pipelines/).
We found that FreeBayes quality was equal to GATK HaplotypeCaller. Additionally,
a lightweight post-alignment preparation method using only de-duplication was
equivalent to GATK's recommended Base Quality Score Recalibration (BQSR) and
realignment around indels, when using good quality input datasets and callers that do local realignment.
- Additional work to [improve variant filtering](https://bcb.io/2014/05/12/wgs-trio-variant-evaluation/),
providing methods to remove low complexity regions (LCRs) that can bias indel results.
We also tuned [GATK's Variant Quality Score Recalibrator](https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php)
(VQSR) and compared it with cutoff-based soft filtering. VQSR requires
a large number of variants and we use it in bcbio with GATK HaplotypeCaller
when your `algorithm parameters` contain high depth samples (`coverage_depth` is not low)
and you are calling on the whole genome (`coverage_interval` is genome) or have
more than 50 regional or exome samples called concurrently.
- An [evaluation of joint calling](https://bcb.io/2014/10/07/joint-calling/)
with GATK HaplotypeCaller, FreeBayes, Platypus and samtools.
This validates the joint calling implementation, allowing scaling
of large population germline experiments. It also demonstrates improved
performance of new callers: samtools 1.0 and Platypus.
- Support for [build 38 of the human genome](https://bcb.io/2015/09/17/hg38-validation/),
improving precision of detection thanks to the improved genome representation.
- bcbio automates post-variant calling annotation to make the outputs easier to
feed directly into your biological analysis.
We annotate variant effects using [snpEff](http://snpeff.sourceforge.net/) or
[Variant Effect Predictor](https://www.ensembl.org/info/docs/tools/vep/index.html) (VEP),
and prepare a [GEMINI database](https://gemini.readthedocs.io/en/latest/)
that associates variants with multiple external annotations in a SQL-based query interface.
GEMINI databases have the most associated external information for human samples (GRCh37/hg19 and hg38)
but are available for any organism with the database populated using the VCF INFO column and predicted effects.
