## Pipelines

### Germline variant calling

bcbio implements configurable SNP, indel and structural variant calling for germline populations. We include whole genome and exome evaluations against reference calls from the [Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) consortium and [Illumina Platinum Genomes](https://www.illumina.com/platinumgenomes.html) project, enabling continuous assessment of new alignment and variant calling algorithms. We regularly report on these comparisons and continue to improve approaches as the community makes new tools available. Here is some of the research that contributes to the current implementation:

* An introduction to the [variant evaluation framework](https://bcb.io/2013/05/06/framework-for-evaluating-variant-detection-methods-comparison-of-aligners-and-callers/). This includes a comparison of the [bwa mem](http://bio-bwa.sourceforge.net/) and [novoalign](http://www.novocraft.com) aligners. We also compared the [FreeBayes](https://github.com/ekg/freebayes), [GATK HaplotypeCaller](https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php) and [GATK UnifiedGenotyper](https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_genotyper_UnifiedGenotyper.php) variant callers.
* An in-depth evaluation of [FreeBayes and BAM post-alignment processing](https://bcb.io/2013/10/21/updated-comparison-of-variant-detection-methods-ensemble-freebayes-and-minimal-bam-preparation-pipelines/). We found that FreeBayes quality was equal to GATK HaplotypeCaller. Additionally, a lightweight post-alignment preparation method using only de-duplication was equivalent to GATK's recommended Base Quality Score Recalibration (BQSR) and realignment around indels, when using good quality input datasets and callers that do local realignment.
* Additional work to [improve variant filtering](https://bcb.io/2014/05/12/wgs-trio-variant-evaluation/), providing methods to remove low complexity regions (LCRs) that can bias indel results. We also tuned [GATK's Variant Quality Score Recalibrator](https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php) (VQSR) and compared it with cutoff-based soft filtering. VQSR requires a large number of variants and we use it in bcbio with GATK HaplotypeCaller when your [algorithm parameters](contents/configuration:algorithm%20parameters) contain high depth samples (`coverage_depth` is not low) and you are calling on the whole genome (`coverage_interval` is genome) or have more than 50 regional or exome samples called concurrently.
* An [evaluation of joint calling](https://bcb.io/2014/10/07/joint-calling/) with GATK HaplotypeCaller, FreeBayes, Platypus and samtools. This validates the joint calling implementation, allowing scaling of large population germline experiments. It also demonstrates improved performance of new callers: samtools 1.0 and Platypus.
* Support for [build 38 of the human genome](https://bcb.io/2015/09/17/hg38-validation/), improving precision of detection thanks to the improved genome representation.

bcbio automates post-variant calling annotation to make the outputs easier to feed directly into your biological analysis. We annotate variant effects using [snpEff](http://snpeff.sourceforge.net/) or [Variant Effect Predictor](https://www.ensembl.org/info/docs/tools/vep/index.html) (VEP), and prepare a [GEMINI database](https://gemini.readthedocs.io/en/latest/) that associates variants with multiple external annotations in a SQL-based query interface. GEMINI databases have the most associated external information for human samples (GRCh37/hg19 and hg38) but are available for any organism with the database populated using the VCF INFO column and predicted effects.

#### Basic germline calling

The best approach to build a bcbio [configuration](contents/configuration:configuration) for germline calling is to use the [automated sample configuration](contents/configuration:automated%20sample%20configuration) with one of the default templates:

* [FreeBayes template](https://github.com/bcbio/bcbio-nextgen/blob/master/config/templates/freebayes-variant.yaml) --Call variants using FreeBayes with a minimal preparation pipeline. This is a freely available unrestricted pipeline fully included in the bcbio installation.
* [GATK HaplotypeCaller template](https://github.com/bcbio/bcbio-nextgen/blob/master/config/templates/gatk-variant.yaml) --Run GATK best practices, including Base Quality Score Recalibration, realignment and HaplotypeCaller variant calling. This requires a license from Broad for commercial use. You need to manually install GATK along with bcbio using downloads from the GATK Broad site or Appistry (see [extra software](contents/installation:extra%20software)).

You may also want to enable [structural variant calling](#structural-variant-calling)
for detection of larger events, which work with either caller. Another good source of inspiration are the configuration files from the [example pipelines](contents/intro:example%20pipelines), which may help identify other configuration variables of interest. A more complex setup with multiple callers and resolution of ensemble calls is generally only useful with a small population where you are especially concerned about
sensitivity. Single caller detection with FreeBayes or GATK HaplotypeCaller provide good resolution of events.

#### Population calling

When calling multiple samples, we recommend calling together to provide improved sensitivity and a fully squared off final callset. To associate samples together in a population add a `metadata` `batch` to the [sample configuration](contents/configuration:samples%20from%20GEO%20or%20SRA):
```yaml
- description: Sample1
  metadata:
    batch: Batch1
- description: Sample2
  metadata:
    batch: Batch1
```
Batching samples results in output VCFs and GEMINI databases containing all merged sample calls. bcbio has two methods to call samples together:

* Batch or pooled calling -- This calls all samples simultaneously by feeding them to the variant caller. This works for smaller batch sizes (< 100 samples) as memory requirements become limiting in larger pools. This is the default approach taken when you specify a `variantcaller` in the [variant calling](contents/configuration:variant%20calling) configuration.

* Joint calling -- This calls samples independently, then combines them together into a single callset by integrating the individual calls. This scales to larger population sizes by avoiding the computational bottlenecks of pooled calling. We recommend joint calling with HaplotypeCaller but also support joint calling with FreeBayes using a custom implementation. Specifying a `jointcaller` along with the appropriate `variantcaller` in the [variant calling](contents/configuration:variant%20calling) configuration enables this:
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

### Cancer variant calling

bcbio supports somatic cancer calling with tumor and optionally matched normal pairs using multiple SNP, indel and structural variant callers. A [full evaluation of cancer calling](http://bcb.io/2015/03/05/cancerval/) validates callers against [synthetic dataset 3 from the ICGC-TCGA DREAM challenge](https://www.synapse.org/#!Synapse:syn312572/wiki/62018). bcbio uses a majority voting ensemble approach to combining calls from multiple SNP and indel callers, and also flattens structural variant calls into a combined representation.

The [example configuration](https://github.com/bcbio/bcbio-nextgen/blob/master/config/examples/cancer-dream-syn3.yaml) for the [whole genome trio (50x) - hg38](contents/intro:whole%20genome%20trio%20(50x)%20-%20hg38) validation is a good starting point for setting up a tumor/normal run on your own dataset. The configuration works similarly to population based calling. Supply a consistent batch for tumor/normal pairs and mark them with the phenotype:
```yaml
- description: your-tumor
  algorithm:
    variantcaller: [vardict, strelka2, mutect2]
  metadata:
    batch: batch1
    phenotype: tumor
- description: your-normal
  algorithm:
    variantcaller: [vardict, strelka2, mutect2]
  metadata:
    batch: batch1
    phenotype: normal
```
Other [somatic variant calling](contents/configuration:somatic%20variant%20calling) configuration
options allow tweaking of the processing parameters. For pairs you want to analyze together, specify a consistent set of `variantcaller` options for both samples.

Cancer calling handles both tumor-normal paired calls and tumor-only calling. To specify a tumor-only sample, provide a single sample labeled with `phenotype: tumor`. Otherwise the configuration and setup is the same as with paired analyses. For tumor-only samples, bcbio will try to remove likely germline variants present in the public databases like 1000 genomes and ExAC, and not in COSMIC. This runs as long as you have a local GEMINI data installation (`--datatarget gemini`) and marks likely germline variants with a `LowPriority` filter. [This post](http://bcb.io/2015/03/05/cancerval/) has more details on the approach and validation.

The standard variant outputs (`sample-caller.vcf.gz`) for tumor calling emphasize somatic differences, those likely variants unique to the cancer. If you have a tumor-only sample and GEMINI data installed, it will also output `sample-caller-germline.vcf.gz`, which tries to identify germline background mutations based on presence in public databases. If you have tumor/normal data and would like to also call likely germline mutations see the [Somatic with germline variants](#somatic-with-germline-variants) section on specifying a germline caller.

We're actively working on improving calling to better account for the heterogeneity and structural variability that define cancer genomes.

### Somatic with germline variants

For tumor/normal somatic samples, bcbio can call both somatic (tumor-specific) and germline (pre-existing) variants. The typical outputs of [cancer variant calling](#cancer-variant-calling) are likely somatic variants acquired by the cancer, but pre-existing germline risk variants are often also diagnostic.

For tumor-only cases we suggest running standard [cancer variant calling](#cancer-variant-calling). Tumor-only inputs mix somatic and germline variants, making it difficult to separate events. For small variants (SNPs and indels) bcbio will attempt to distinguish somatic and germline mutations using the presence of variants in population databases.

To option somatic and germline calls for your tumor/normal inputs, specify which callers to use for each step in the [variant calling](contents/configuration:variant%20calling) configuration:
```yaml
description: your-normal
variantcaller:
   somatic: vardict
   germline: freebayes
```
bcbio does a single alignment for the normal sample, then splits at the variant calling steps using this normal sample to do germline calling. In this example, the output files are:
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

### Somatic tumor only CNVs

Copy number variation (CNVs) detection in tumor only samples requires accurately representing the non-somatic capture and sequencing background in the absence of a matched sample. Capture or sequencing specific coverage differences can trigger false positives or negatives. Without a matched normal to remove these artifacts, you can use a process matched set of unrelated samples to build a Panel of Normals (PoN) for the background correction.

To create these, collect all the samples you plan to use for the panel of normals and run through an [automated sample configuration](contents/configuration:automated%20sample%20configuration) as a single batch with the background samples set as control and any nonbackground as the non-background. An example sample CSV:
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

Once you have the panel of normals, use them as background in any tumor only project with the same sequencing and capture process in your [variant calling](contents/configuration:variant%20calling) configuration:
```yaml
svcaller: [gatk-cnv, seq2c]
variant_regions: your_regions.bed
background:
  cnv_reference:
    cnvkit: ../pon/your_regions-cnvkit-pon.cnn
    gatk-cnv: ../pon/your_regions-gatk-cnv-pon.hdf5
    seq2c: ../pon/your_region-seq2c-pon.txt
```

### Structural variant calling

bcbio can detect larger structural variants like deletions, insertions, inversions and copy number changes for both germline population and cancer variant calling, based on validation against existing truth sets:
* [Validation of germline structural variant detection](https://bcb.io/2014/08/12/validated-whole-genome-structural-variation-detection-using-multiple-callers/) using multiple calling methods to validate against deletions in NA12878. This implements a pipeline that works in tandem with SNP and indel calling to detect larger structural variations like deletions, duplications, inversions and copy number variants (CNVs).
* [Validation of tumor/normal calling](https://bcb.io/2015/03/05/cancerval/) using the synthetic DREAM validation set. This includes validation of additional callers against duplications, insertions and inversions.

To enable structural variant calling, specify `svcaller` options in the algorithm section of your configuration:
```yaml
- description: Sample
  algorithm:
    svcaller: [lumpy, manta, cnvkit]
```
The best supported callers are [Lumpy](https://github.com/arq5x/lumpy-sv) and
[Manta](https://github.com/Illumina/manta), for paired end and split read calling, [CNVkit](https://cnvkit.readthedocs.io/en/latest/) for read-depth based CNV calling, and
[WHAM](https://github.com/jewmanchue/wham) for association testing. We also support [DELLY](https://github.com/tobiasrausch/delly), another excellent paired end and split read caller, although it is slow on large whole genome datasets.

### fast RNA-seq

This mode of `bcbio-nextgen` quantitates transcript expression using [Salmon](https://salmon.readthedocs.io/en/latest/) and does nothing else. It is an order of magnitude faster or more than running the full RNA-seq analysis. The cost of the increased speed is that you will have much less information about your samples at the end of the run, which can make troubleshooting trickier. Invoke with `analysis: fastrna-seq`.

### single-cell RNA-seq

bcbio-nextgen supports universal molecular identifiers (UMI) based single-cell RNA-seq analyses. If your single-cell prep does not use universal molecular identifiers (UMI), you can most likely just run the standard RNA-seq pipeline and use the results from that. The UMI are used to discard reads which are possibly PCR duplicates and is very helpful for removing some of the PCR duplicate noise that can dominate single-cell experiments.

Unlike the standard RNA-seq pipeline, the single-cell pipeline expects the FASTQ input files to not be separated by cellular barcode, so each file is a mix of cells identified by a cellular barcode (CB), and unique reads from a transcript are identified with a UMI. bcbio-nextgen inspects each read, identifies the cellular barcode and UMI and puts them in the read name. Then the reads are aligned to the transcriptome with [RapMap](https://github.com/COMBINE-lab/RapMap) and the number of reads aligning to each transcript is counted for each cellular barcode. The output is a table of counts with transcripts as the rows and columns as the cellular barcodes for each input FASTQ file.

Optionally the reads can be quantitated with `kallisto` to output transcript compatibility counts rather than counts per gene ([TCC paper](https://doi.org/10.1186/s13059-016-0970-8)).

To extract the UMI and cellular barcodes from the read, bcbio-nextgen needs to know where the UMI and the cellular barcode are expected to be in the read. Currently there is support for two schemes, the inDrop system from the Harvard single-cell core facility and CEL-seq. If bcbio-nextgen does not support your UMI and barcoding scheme, please open up an issue and we will help implement support for it.

Most of the heavy lifting for this part of bcbio-nextgen is implemented in the [umis](https://github.com/vals/umis) repository.

### smallRNA-seq

bcbio-nextgen also implements a configurable best-practices pipeline for smallRNA-seq quality controls, adapter trimming, miRNA/isomiR quantification and other small RNA detection.

* Adapter trimming:
  * [atropos](https://atropos.readthedocs.io/en/latest/guide.html)
  * [dnapi](https://github.com/jnktsj/DNApi) for adapter de-novo detection
* Sequence alignment:
  * [STAR](https://code.google.com/archive/p/rna-star) for genome annotation
  * bowtie, _bowtie2_ and [hisat2](https://daehwankimlab.github.io/hisat2/) for genome annotation as an option
* Specific small RNAs quantification (miRNA/tRNAs...):
  * [seqbuster](https://github.com/lpantano/seqbuster) for miRNA annotation
  * [MINTmap](https://github.com/TJU-CMC-Org/MINTmap) for tRNA fragments annotation
  * [miRge2](https://github.com/mhalushka/miRge) for alternative small RNA quantification. To setup this tool, you need to install manually miRge2.0, and download the library data for your species. Read how to install and download the [data](https://github.com/mhalushka/miRge#download-libraries). If you have `human` folder at `/mnt/data/human` the option to pass to resources will be `/mnt/data`. Then setup `resources`:
    ```yaml
    resources:
        mirge:
            options: ["-lib $PATH_TO_PARENT_SPECIES_LIB"]
    ```
* Quality control: [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* Other small RNAs quantification:
  * [seqcluster](https://github.com/lpantano/seqcluster)
  * [mirDeep2](https://www.mdc-berlin.de/content/mirdeep2-documentation) for miRNA prediction

The pipeline generates a _RMD_ template file inside `report` folder that can be rendered with knitr. An example of the report is [here](https://github.com/lpantano/mypubs/blob/master/srnaseq/mirqc/ready_report.md). Count table (`counts_mirna.tst`) from mirbase miRNAs will be inside `mirbase` or final project folder. Input files for [isomiRs](https://github.com/lpantano/isomiRs) package for isomiRs analysis will be inside each sample in `mirbase` folder. If mirdeep2 can run, count table (`counts_mirna_novel.tsv`) for novel miRNAs will be inside `mirdeep2` or final project folder. tdrmapper results will be inside each sample inside `tdrmapper` or final project folder.

### ChIP/ATAC-seq

The bcbio-nextgen implementation of ChIP-seq aligns, removes multimapping reads, calls peaks with a paired input file using MACS2 and outputs a set of greylist regions for filtering possible false peaks in regions of high depth in the input file.

* Adapter trimming: [atropos](https://atropos.readthedocs.io/en/latest/guide.html)
* Sequence alignment: [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [bwa mem](http://bio-bwa.sourceforge.net/)
* Peak calling: [macs2](https://github.com/taoliu/MACS)
* Greylisting: [chipseq-greylist](https://github.com/roryk/chipseq-greylist)
* Quality control: [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

### Methylation

Whole genome bisulfite sequencing is supported using the [bismark2](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) pipeline. It can be turned on by setting `analysis` to `wgbs-seq`. Right now we only support the TruSeq Methyl Capture EPIC kit -- if you want to run a different setup, please [let us know](https://github.com/bcbio/bcbio-nextgen/issues) and we can add support for it. Consider this a beta feature.

### Standard

This pipeline implements `alignment` and `qc` tools. Furthermore, it will run [qsignature](http://sourceforge.net/p/adamajava/wiki/qSignature/) to detect possible duplicated samples, or mislabeling. It uses SNPs signature to create a distance matrix that helps easily to create groups. The project yaml file will show the number of total samples analyzed, the number of very similar samples, and samples that could be duplicated.

#### Configuration

We will assume that you installed bcbio-nextgen with the automated installer, and so your default [bcbio_system.yaml](http://github.com/bcbio/bcbio-nextgen/blob/master/config/bcbio_system.yaml) file is configured correctly with all of the tools pointing to the right places. If that is the case, to run bcbio-nextgen on a set of samples you just need to set up a YAML file that describes your samples and what you would like to do to them.

Let's say that you have a single paired-end control lane, prepared with the [Illumina TruSeq Kit](https://www.illumina.com/products/by-type/sequencing-kits/library-prep-kits/truseq-rna-v2.html) from a human. Here is what a well-formed sample YAML file for that
RNA-seq experiment would look like:
```yaml
fc_date: '070113'
fc_name: control_experiment
upload:
  dir: final
details:
  - files: [/full/path/to/control_1.fastq, /full/path/to/control_2.fastq]
    description: 'Control_rep1'
    genome_build: GRCh37
    analysis: RNA-seq
    algorithm:
      aligner: tophat2
      quality_format: Standard
      trim_reads: read_through
      adapters: [truseq, polya]
      strandedness: unstranded
```
`fc_date` and `fc_name` will be combined to form a prefix to name intermediate files, and can be set to whatever you like. `upload` is explained pretty well in the [configuration
documentation](contents/configuration:upload) and the above will direct bcbio-nextgen to put the output files from the pipeine into the `final` directory. Under `details` is a list of sections each describing a sample to process. You can set many
[parameters](configuration) under each section but most of the time just setting a few like the above is all that is necessary. `analysis` tells bcbio-nextgen to run the best-practice RNA-seq pipeline on this sample.

In the above, since there are two files, `control_1.fastq` and `control_2.fastq` will be automatically run as paired-end data. If you have single end data you can just supply one file and it will run as single-end. The `description` field will be used to eventually rename the files, so make it very evocative since you will be looking at it a lot later. `genome_build` is self-explanatory.

Sometimes you need a little bit more flexibility than the standard pipeline, and the `algorithm` section has many options to fine-tune the behavior of the algorithm. `quality_format` tells bcbio-nextgen what quality format your FASTQ inputs are using, if your samples were sequenced any time past 2009 or so, you probably want to set it to
`Standard`. Adapter read-through is a problem in RNA-seq libraries, so we want to trim off possible adapter sequences on the ends of reads, so `trim_reads` is set to `read_through`, which will also trim off poor quality ends. Since your library is a RNA-seq library prepared with the TruSeq kit, the set of adapters to trim off are the TruSeq adapters and possible polyA tails, so `adapters` is set to both of those.
`strandedness` can be set if your library was prepared in a strand-specific manner and can be set to firststrand, secondstrand or unstranded (the default).

#### Multiple samples

Lets say you have a set of mouse samples to analyze and each sample is a single lane of single-end RNA-seq reads prepared using the NextEra kit. There are two case and two control samples. Here is a sample configuration file for that analysis:
```yaml
fc_date: '070113'
fc_name: mouse_analysis
upload:
  dir: final
details:
  - files: [/full/path/to/control_rep1.fastq]
    description: 'Control_rep1'
    genome_build: mm10
    analysis: RNA-seq
    algorithm:
      aligner: tophat2
      quality_format: Standard
      trim_reads: read_through
      adapters: [nextera, polya]
  - files: [/full/path/to/control_rep2.fastq]
    description: 'Control_rep2'
    genome_build: mm10
    analysis: RNA-seq
    algorithm:
      aligner: tophat2
      quality_format: Standard
      trim_reads: read_through
      adapters: [nextera, polya]
  - files: [/full/path/to/case_rep1.fastq]
    description: 'Case_rep1'
    genome_build: mm10
    analysis: RNA-seq
    algorithm:
      aligner: tophat2
      quality_format: Standard
      trim_reads: read_through
      adapters: [nextera, polya]
  - files: [/full/path/to/case_rep2.fastq]
    description: 'Case_rep2'
    genome_build: mm10
    analysis: RNA-seq
    algorithm:
      aligner: tophat2
      quality_format: Standard
      trim_reads: read_through
      adapters: [nextera, polya]
```
More samples are added just by adding more entries under the details section. This is tedious and error prone to do by hand, so there is an automated [template](contents/configuration:automated%20sample%20configuration) system for common experiments. You could set up the previous experiment by making a mouse version of the [illumina-rnaseq](https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/templates/illumina-rnaseq.yaml) template file and saving it to a local file such as `illumina-mouse-rnaseq.yaml`. Then you can set up the sample file using the templating system:
```shell
bcbio_nextgen.py -w template illumina-mouse-rnaseq.yaml mouse_analysis
/full/path/to/control_rep1.fastq /full/path/to/control_rep2.fastq
/full/path/to/case_rep1.fastq /full/path/to/case_rep2.fastq
```
If you had paired-end samples instead of single-end samples, you can still use the template system as long as the forward and reverse read filenames are the same, barring a _1 and _2. For example: control\_1.fastq and control\_2.fastq will be detected as paired and combined in the YAML file output by the templating system.
