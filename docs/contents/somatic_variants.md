# Somatic (cancer) variants

bcbio supports tumor only and tumor-normal small (SNV and indels),
copy number (CNV) variant calling using multiple tools.
We recommend to start with `vardict` and `mutect2`.  
bcbio also supports a majority voting ensemble approach
to combine calls from multiple callers.

## Workflow1 - T/N
This [example](https://github.com/bcbio/bcbio-nextgen/blob/master/config/examples/cancer-dream-syn3.yaml)
runs a tumor-normal calling with mutect2 and vardict.
Supply a consistent batch for tumor/normal pairs and mark them with the phenotype:

```yaml
- description: sample_tumor
  algorithm:
    variantcaller: [vardict, mutect2]
  metadata:
    batch: batch1
    phenotype: tumor
- description: sample_normal
  algorithm:
    variantcaller: [vardict, mutect2]
  metadata:
    batch: batch1
    phenotype: normal
```

This example calls and validates variants using multiple approaches in a paired tumor/normal
cancer sample from the [ICGC-TCGA DREAM challenge](https://www.synapse.org/#!Synapse:syn312572/wiki/58893).
It uses [synthetic dataset 3](https://www.synapse.org/#!Synapse:syn312572/wiki/62018) which has multiple subclones,
enabling detection of lower frequency variants.

The configuration and data file has downloads for exome only and whole genome analyses.
It enables exome by default, but you can use the larger whole genome evaluation
by uncommenting the relevant parts of the configuration and retrieval script.

### 1. Get the data:
```shell
mkdir -p cancer-dream-syn3/config cancer-dream-syn3/input cancer-dream-syn3/work
cd cancer-dream-syn3/config
wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/cancer-dream-syn3.yaml
cd ../input
wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/cancer-dream-syn3-getdata.sh
bash cancer-dream-syn3-getdata.sh
```

### 2. Review parameters in the yaml file:
```yaml
# Cancer tumor/normal calling evaluation using synthetic dataset 3
# from the ICGC-TCGA DREAM challenge:
# https://www.synapse.org/#!Synapse:syn312572/wiki/62018
---
details:
- algorithm:
    aligner: bwa
    mark_duplicates: true
    remove_lcr: true
    variantcaller: [mutect2, vardict]
    variant_regions: ../input/NGv3.bed
    # svcaller: [cnvkit, lumpy, delly]
    # coverage_interval: amplicon
  analysis: variant2
  description: syn3-normal
  #files: ../input/synthetic.challenge.set3.normal.bam
  files:
  - ../input/synthetic_challenge_set3_normal_NGv3_1.fq.gz
  - ../input/synthetic_challenge_set3_normal_NGv3_2.fq.gz
  genome_build: GRCh37
  metadata:
    batch: syn3
    phenotype: normal
 - algorithm:
    aligner: bwa
    mark_duplicates: true
    remove_lcr: true
    variantcaller: [mutect2, vardict]
    variant_regions: ../input/NGv3.bed
    validate: ../input/synthetic_challenge_set3_tumor_20pctmasked_truth.vcf.gz
    validate_regions: ../input/synthetic_challenge_set3_tumor_20pctmasked_truth_regions.bed
    # svcaller: [cnvkit, lumpy, delly]
    # coverage_interval: amplicon
    #   svvalidate:
    #     DEL: ../input/synthetic_challenge_set3_tumor_20pctmasked_truth_sv_DEL.bed
    #     DUP: ../input/synthetic_challenge_set3_tumor_20pctmasked_truth_sv_DUP.bed
    #     INS: ../input/synthetic_challenge_set3_tumor_20pctmasked_truth_sv_INS.bed
    #     INV: ../input/synthetic_challenge_set3_tumor_20pctmasked_truth_sv_INV.bed
  analysis: variant2
  description: syn3-tumor
  #files: ../input/synthetic.challenge.set3.tumor.bam
  files:
  - ../input/synthetic_challenge_set3_tumor_NGv3_1.fq.gz
  - ../input/synthetic_challenge_set3_tumor_NGv3_2.fq.gz
  genome_build: GRCh37
  metadata:
    batch: syn3
    phenotype: tumor
  fc_date: '2014-08-13'
  fc_name: dream-syn3
upload:
  dir: ../final
```

### 3. Run bcbio project
```shell
cd ../work
bcbio_nextgen.py ../config/cancer-dream-syn3.yaml -n 8
```

Cancer calling handles both tumor-normal paired calls and tumor-only calling. To specify a tumor-only sample, provide a single sample labeled with `phenotype: tumor`. Otherwise the configuration and setup is the same as with paired analyses. For tumor-only samples, bcbio will try to remove likely germline variants present in the public databases like 1000 genomes and ExAC, and not in COSMIC. This runs as long as you have a local GEMINI data installation (`--datatarget gemini`) and marks likely germline variants with a `LowPriority` filter. [This post](http://bcb.io/2015/03/05/cancerval/) has more details on the approach and validation.

The standard variant outputs (`sample-caller.vcf.gz`) for tumor calling emphasize somatic differences, those likely variants unique to the cancer. If you have a tumor-only sample and GEMINI data installed, it will also output `sample-caller-germline.vcf.gz`, which tries to identify germline background mutations based on presence in public databases.

We're actively working on improving calling to better account for the heterogeneity and structural variability that define cancer genomes.

## Workflow2: Somatic and germline variants

For tumor/normal somatic samples, bcbio can call both somatic (tumor-specific)
and germline (pre-existing) variants.

To option somatic and germline calls for your tumor/normal inputs, specify
which callers to use for each step in the  configuration:

```yaml
description: your-normal
variantcaller:
   somatic: vardict
   germline: gatk-haplotype
```
bcbio does a single alignment for the normal sample, then splits at the variant
calling steps using this normal sample to do germline calling. In this example, the output files are:

* `your-tumor/your-tumor-vardict.vcf.gz` -- Somatic calls from the tumor samples using the normal as background to subtract existing calls.
* `your-normal/your-normal-freebayes.vcf.gz` -- Germline calls on the normal sample.

Germline calling supports multiple callers, and other configuration options like ensemble and structural variant calling inherit from the remainder configuration. For example, to use 3 callers for somatic and germline calling, create ensemble calls for both and include germline and somatic events from two structural variant callers:
```yaml
variantcaller:
   somatic: [vardict, strelka2, mutect2]
   germline: [freebayes, gatk-haplotype, strelka2]
ensemble:
   numpass: 2
svcaller: [manta, cnvkit]
```
In addition to the somatic and germline outputs attached to the tumor and normal sample outputs as described above, you'll get:
* `your-tumor/your-tumor-manta.vcf.gz` -- Somatic structural variant calls for each specified `svcaller`. These will have genotypes for both the tumor and normal samples, with somatic calls labeled as PASS variants.
* `your-normal/your-normal-manta.vcf.gz` -- Germline structural variant calls for each specified `svcaller`. We expect these to be noisier than the somatic calls due to the lack of a reference sample to help remove technical noise.

## Workflow3: Somatic tumor only CNVs

Copy number variation (CNVs) detection in tumor only samples requires accurately representing the non-somatic capture and sequencing background in the absence of a matched sample. Capture or sequencing specific coverage differences can trigger false positives or negatives. Without a matched normal to remove these artifacts, you can use a process matched set of unrelated samples to build a Panel of Normals (PoN) for the background correction.

To create these, collect all the samples you plan to use for the panel of normals and run through an `automated sample configuration` as a single batch with the background samples set as control and any nonbackground as the non-background. An example sample CSV:
```
samplename,description,svclass,batch
background1.bam,background1,control,pon_build
background2_R1.fq.gz;background2_R2.fq.gz,background2,control,pon_build
testsample.bam,testsample,pon_build
```
and template YAML:
```yaml
details:
  - analysis: variant2
    genome_build: hg38
    algorithm:
      svcaller: [gatk-cnv, seq2c]
      variant_regions: your_regions.bed
```
After running, collect the panel of normal files from each calling method:
* gatk-cnv: _work/structural/testsample/bins/background1-pon_build-pon.hdf5_
* seq2c: This doesn't have a default panel of normals file format so we create a bcbio specific one as a concatenation of the read mapping file (_final/date_project/seq2c-read_mapping.txt_) and coverage file (_final/date_project/seq2c-coverage.tsv_) outputs for the background samples. When fed to future bcbio runs, it will correctly extract and re-use this file as background.
* CNVkit: _final/testsample/testsample-cnvkit-background.cnn_

CNVkit and gatk-cnv cannot be run together, because they require different, incompatible normalization schemes.

Once you have the panel of normals, use them as background in any tumor only project with the same sequencing and capture process in your `variant calling` configuration:
```yaml
svcaller: [gatk-cnv, seq2c]
variant_regions: your_regions.bed
background:
  cnv_reference:
    cnvkit: ../pon/your_regions-cnvkit-pon.cnn
    gatk-cnv: ../pon/your_regions-gatk-cnv-pon.hdf5
    seq2c: ../pon/your_region-seq2c-pon.txt
```

## Workflow4: Cancer-like mixture with Genome in a Bottle samples

This example simulates somatic cancer calling using a mixture of two Genome in a Bottle samples, NA12878 as the "tumor" mixed with NA24385 as the background. The [Hartwig Medical Foundation](https://www.hartwigmedicalfoundation.nl/en/) and [Utrecht Medical Center](https://www.umcutrecht.nl/en/Research/Strategic-themes/Cancer) generated this "tumor/normal" pair by physical mixing of samples prior to sequencing. The GiaB FTP directory has [more details on the design and truth sets](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/use_cases/mixtures/UMCUTRECHT_NA12878_NA24385_mixture_10052016/README-NA12878_NA24385_mixture.txt). The sample has variants at 15% and 30%, providing the ability to look at lower frequency mutations.

To get the data:
```shell
wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/examples/cancer-giab-na12878-na24385-getdata.sh
bash cancer-giab-na12878-na24385-getdata.sh
```
Then run the analysis with:
```shell
cd work
bcbio_nextgen.py ../config/cancer-giab-na12878-na24385.yaml -n 16
```

## Parameters
- `variantcaller`: should be consistent for T/N pairs
* `min_allele_fraction` Minimum allele fraction to detect variants in heterogeneous tumor samples, set as the float or integer __percentage__ to resolve (i.e. 10 = alleles in 10% of the sample). Defaults to 10. Specify this in the tumor sample of a tumor/normal pair. It is percentage, not ratio, it is divided /100.0 when calling vardict!
* `use_lowfreq_filter: false`. When set, forces vardict to report variants with low allelec frequency, useful to call variants in panels with high coverage (>1000x). The default (option is not set to false in the config) is to use low frequency filter, i.e. variants could be underreported (variant VAF is above min_allele_fraction but rejected by the filter).
* also see parameters in `germline_varinants` user story.

### Ensemble variant calling
A simple majority rule ensemble classifier builds a final callset based on the intersection of calls.
It selects variants represented in at least a specified number of callers:
```yaml
variantcaller: [mutect2, varscan, freebayes, vardict]
ensemble:
  numpass: 2
  use_filtered: false
```

This example selects variants present in 2 out of the 4 callers and does not use filtered calls (the default behavior). Because of the difficulties of producing a unified FORMAT/genotype field across callers, the ensemble outputs contains a mix of outputs from the different callers. It picks a representative sample in the order of specified caller, so in the example above would have a MuTect2 call if present, otherwise a VarScan call if present, otherwise a FreeBayes call. This may require custom normalization scripts during post-processing when using these calls. [bcbio.variation.recall](https://github.com/chapmanb/bcbio.variation.recall) implements this approach, which handles speed and file sorting limitations in the [bcbio.variation](https://github.com/chapmanb/bcbio.variation) approach.

This older approach uses the [bcbio.variation](https://github.com/chapmanb/bcbio.variation) toolkit to perform the consolidation. An example configuration in the
`algorithm` section is:
```yaml
variantcaller: [gatk, freebayes, samtools, gatk-haplotype, varscan]
ensemble:
  format-filters: [DP < 4]
  classifier-params:
    type: svm
  classifiers:
    balance: [AD, FS, Entropy]
    calling: [ReadPosEndDist, PL, PLratio, Entropy, NBQ]
  trusted-pct: 0.65
```
The `ensemble` set of parameters configure how to combine calls from the
multiple methods:
* `format-filters` A set of filters to apply to variants before combining. The example removes all calls with a depth of less than 4.
* `classifier-params` Parameters to configure the machine learning approaches used to consolidate calls. The example defines an SVM classifier.
* `classifiers` Groups of classifiers to use for training and evaluating during machine learning. The example defines two set of criteria for distinguishing reads with allele balance issues and those with low calling support.
* `trusted-pct` Define threshold of variants to include in final callset. In the example, variants called by more than 65% of the approaches (4 or more callers) pass without being requiring SVM filtering.

### UMIs
Unique molecular identifiers (UMIs) are short random barcodes used to tag
individual molecules and avoid amplification biased. Both single cell RNA-seq
and variant calling support UMIs.
For variant calling, [fgbio](https://github.com/fulcrumgenomics/fgbio) collapses
sequencing duplicates for each UMI into a single consensus read prior to running
re-alignment and variant calling. This requires `mark_duplicates: true` (the default) since it uses position based duplicates and UMI tags for collapsing duplicate reads into consensus sequences.

To help with preparing fastq files with UMIs bcbio provides a script `bcbio_fastq_umi_prep.py`. This handles two kinds of UMI barcodes:
* Separate UMIs: it converts reads output by an Illumina as 3 files (read 1, read 2, and UMIs).
* Duplex barcodes with tags incorporated at the 5' end of read 1 and read 2

In both cases, these get converted into paired reads with UMIs in the fastq names, allowing specification of `umi_type: fastq_name` in your bcbio YAML configuration. The script runs on a single set of files or autopairs an entire directory of fastq files. To convert a directory with separate UMI files:
```shell
bcbio_fastq_umi_prep.py autopair -c <cores_to_use> <list> <of> <fastq> <files>
```
To convert duplex barcodes present on the ends of read 1 and read 2:
```shell
bcbio_fastq_umi_prep.py autopair -c <cores_to_use> --tag1 5 --tag2 5 <list> <of> <fastq> <files>
```
If you want to prepare your FASTQ files for use with the `umi_type: fastq_name` option and they don't follow the two barcode schemes listed you can pre-transform the FASTQ files yourself by putting `:UMI_yourumisequence` in the read name. Here is an example:
```
@A00574:89:HCLMTDRXX:1:2101:1425:1016:UMI_CGAACGTGTACACG 2:N:0:CTGAAGCT+GGCTCTGA
GTGATATAATTTATTTTCTTAAAATAGCCATGCTGGCTGGAGCCACAGCAGTTTACTCCCAGTTCATTACTCAGCTAACAGACGAAAACCAGT
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
```
For more complex UMI schemes please open up an issue and we can help you write a transformation to prepare your files. We use <https://github.com/vals/umis> to do these more complex transformations.

Configuration options for UMIs:
* `umi_type` The UMI/cellular barcode scheme used for your data. For variant analysis with UMI based consensus calling, supports either `fastq_name` with UMIs in read names or the path to a fastq file with UMIs for each aligned read.
* `correct_umis: [path/to/whitelist_umi.txt]`. For a restricted set of UMIs specify a text file (one UMI per line). UMIs will be corrected with <http://fulcrumgenomics.github.io/fgbio/tools/latest/CorrectUmis.html>

You can adjust the [fgbio default options](https://github.com/bcbio/bcbio-nextgen/blob/8a76c9e546cb79621707082fd763bd643e0e9652/bcbio/ngsalign/postalign.py#L208) by adjusting `resources`. The most common change is modifying the minimum number of reads as input to consensus sequences.
This default to 1 to avoid losing reads, 3 is recommended for high depth panels:
```yaml
resources:
  fgbio:
    options: [--min-reads, 3]
```

## Output

### Project directory:
* `grading-summary.csv` -- Grading details comparing each sample to a reference set of calls. This will only have information when providing a validation callset.
* `BATCH-caller.vcf` -- Variants called for a population/batch of samples by a particular caller.
* `BATCH-caller.db` -- A [GEMINI database](https://github.com/arq5x/gemini) associating variant calls with a wide variety of third party annotations. This provides a queryable framework for assessing variant quality statistics.

### Sample directories:
* `SAMPLE-caller.vcf` -- Variants calls for an individual sample.
* `SAMPLE-gdc-viral-completeness.txt` -- Optional viral contamination estimates. File is of the format **depth, 1x, 5x, 25x**. **depth** is the number of reads aligning to the virus. **1x, 5x, 25x** are percentage of the viral sequence covered by reads of 1x, 5x, 25x depth. Real viral contamination will have broad coverage across the entire genome, so high numbers for these values, depending on sequencing depth. High depth and low viral sequence coverage means a likely false positive.

## Validation
bcbio pre-installs standard truth sets for performing validation, and also allows use of custom local files for assessing reliability of your runs:
* `validate` A VCF file of expected variant calls to perform validation and grading of small variants (SNPs and indels) from the pipeline. This provides a mechanism to ensure consistency of calls against a known set of variants, supporting comparisons to
genotyping array data or reference materials.
* `validate_regions` A BED file of regions to evaluate small variant calls in. This defines specific regions covered by the `validate` VCF file.
* `svvalidate` -- Dictionary of call types and pointer to BED file of known regions. For example: `DEL: known_deletions.bed` does deletion based validation of outputs against the BED file.

Each option can be either the path to a local file, or a partial path to a file in the pre-installed truth sets. For instance, to validate an NA12878 run against the [Genome in a Bottle](https://github.com/genome-in-a-bottle) truth set:
```yaml
validate: giab-NA12878/truth_small_variants.vcf.gz
validate_regions: giab-NA12878/truth_regions.bed
svvalidate:
  DEL: giab-NA12878/truth_DEL.bed
```
For cancer validations:
* `giab-NA12878-NA24385-somatic` -- A [sequenced NA12878/NA24385 mixture](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/use_cases/mixtures/UMCUTRECHT_NA12878_NA24385_mixture_10052016/) providing a somatic-like truth set for detecting low frequency events. Build: Truth sets: small_variants, regions. Builds: GRCh37, hg38
* `dream-syn3` -- Synthetic dataset 3 from the [ICGC-TCGA DREAM mutation calling challenge](https://www.synapse.org/#!Synapse:syn312572/wiki/62018). Truth sets: small_variants, regions, DEL, DUP, INV, INS. Builds: GRCh37.
* `dream-syn4` -- Synthetic dataset 4 from the [ICGC-TCGA DREAM mutation calling challenge](https://www.synapse.org/#!Synapse:syn312572/wiki/62018). Truth sets: small_variants, regions, DEL, DUP, INV. Builds: GRCh37.
* `dream-syn3-crossmap` -- Synthetic dataset 3 from the [ICGC-TCGA DREAM mutation calling challenge](https://www.synapse.org/#!Synapse:syn312572/wiki/62018) converted to human build 38 coordinates with CrossMap. Truth sets: small_variants, regions, DEL, DUP, INV, INS. Builds: hg38.
* `dream-syn4-crossmap` -- Synthetic dataset 4 from the [ICGC-TCGA DREAM mutation calling challenge](https://www.synapse.org/#!Synapse:syn312572/wiki/62018) converted to human build 38 coordinates with CrossMap. Truth sets: small_variants, regions, DEL, DUP, INV. Builds: hg38.


A [full evaluation of cancer calling](http://bcb.io/2015/03/05/cancerval/) validates callers against [synthetic dataset 3 from the ICGC-TCGA DREAM challenge](https://www.synapse.org/#!Synapse:syn312572/wiki/62018).

## TODO
