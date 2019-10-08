## 1.1.6 (in progress)

- GATK ApplyBQSRSpark: avoid StreamClosed issue with GATK 4.1+
- RNA-seq: fixes for cufflinks preparation due to python3 transition.
- RNA-seq: output count tables from tximport for genes and transcripts. These
are in `bcbioRNASeq/results/date/genes/counts` and 
`bcbioRNASeq/results/data/transcripts/counts`.
- qualimap (RNA-seq): disable stranded mode for qualimap, as it gives incorrect
results with the hisat2 aligner and for RNA-seq just setting it to unstranded
- Add `quantify_genome_alignments` option to use genome alignments to quantify
  with Salmon.
- Add `--validateMappings` flag to Salmon read quantification mode.
- VEP cache is not installing anymore from bcbio run
- Add support for Salmon SA method when STAR alignments are not available 
  (for hg38).
- Add support for the new read model for filtering in Mutect2. This is
  experimental, and a little flaky, so it can optionally be turned on via:
  `tools_on: mutect2_readmodel`. Thanks to @lbeltrame for implementing this
  feature and doing a ton of work debugging.
- Swap pandas `from_csv` call to `read_csv`.
- Make STAR respect the `transcriptome_gtf` option.
- Prefix regular expression with r. Thanks to @smoe for finding all of these.
- Add informative logging messages at beginning of bcbio run. Includes the version
  and the configuration files being used.
- Swap samtools mpileup to use bcftools mpileup as samtools mpileup is being 
  deprecated (https://github.com/samtools/samtools/releases/tag/1.9).
- Ensure locale is set to one supporting UTF-8 bcbio-wide. This may need to get
  reverted if it introduces issues.
- Added hg38 support for STAR. We did this by taking hg38 and removing the alts,
  decoys and HLA sequences.
- Added support for the arriba fusion caller.
- Added back missing programs from the version provenance file. Fixed formatting
  problems introduced by switch to python3.
- Added initial support for whole genome bisulfite sequencing using bismark. Thanks to
  @hackdna for implementing this and @jnhutchinson for drafting the initial
  pipeline. This is a work in progress in collaboration with @gcampanella, who
  has a similar implementation with some extra features that we will be merging
  in soon.
- qualimap for RNA-seq runs on the downsampled BAM files by default. Set 
  `tools_on: [qualimap_full]` to run on the full BAM files.
- Add STAR junction files to the files captured at the end of a run.
  
## 1.1.5 (12 April 2019)

- Fixes for Python3 incompatibilities on distributed IPython runs.
- Numerous smaller Python3 incompatibilities with strings/unicode and types.
  Thanks to the community for reporting these.
- GATK HaplotypeCaller: correctly apply skipping of marked duplicates only
  for amplicon runs. Thanks to Ben Liesfeld.
- Fix format detection for bzip2 fastq inputs.
- Support latest GATK4 MuTect2 (4.1.1.0) with changes to ploidy and reference
  parameters.
- Support changes to GATK4 for VQSR --resource specification in 4.1.1.0. Thanks
  to Timothee Cezard.
- Support latest bedtools (2.28.0) which expects SAM heads for bgzipped BED
  inputs.

## 1.1.4 (3 April 2019)

- Move to Python 3.6. A python2 environment in the install runs non python3
  compatible programs. The codebase is still compatible with python 2.7 but
  will only get run and tested on python 3 for future releases.
- RNA-seq: fix for race condition when creating the pizzly cache
- RNA-seq: Add Salmon to multiqc report.
- RNA-seq single-cell/DGE: Properly strip transcript versions from GENCODE GTFs.
- RNA-seq: Faster and more flexible rRNA biotype lookup.
- Move to R3.5.1, including updates to all CRAN and Bioconductor packages.
- tumor-only germline prioritization: provide more useful germline filtering
  based on prioritization INFO tag (EPR) rather than filter field.
- Install: do not require fabric for tool and data installs, making full codebase
  compatible with python 3.
- variant: Filter out variants with missing ALT alleles output by GATK4.
- GATK: enable specification of spark specific parameters with `gatk-spark`
  resources.
- RNA-seq single-cell/DGE: added `demultiplexed` option. If set to True, treat the
  data as if it has already been demultiplexed into cells/wells.
- Multiple orders of magnitude faster templating with thousands of input files.

## 1.1.3 (29 January 2019)

- CNV: support background inputs for CNVkit, GATK4 CNV and seq2c. Allows
  pre-computed panel of normals for tumor-only or single sample CNV calling.
- variant: avoid race condition on processing input BED files for variant
  calling when no pre-specific variant_regions available.
- structural variation upload: avoid uploading multiple batched calls into
  sample directories. For lumpy will now have a single output per batch in a
  sample folder.
- install: respect pre-specified bioconda and conda-forge in condarc
  configuration. Allows use of custom package mirrors.
- seq2c: move specialized pre-call calculation upstream to coverage estimation.
  Allows use of seq2c in CWL runs.
- MultiQC upload: fix bug where results from parallel run not moved to final
  directory.
- GATK4 CNV: fix for standardize VCF output, correcting number of columns.
- RNA-seq variation: fix for over-filtering variants near splice junctions with
  STAR.
- Structural variant gene annotation: simplify and handle issues with
  multidirectional comparisons. Handle issues with out of order start/end from CNVkit.
- Catch and report unicode characters in templating or YAML descriptions.

## 1.1.2 (12 December 2018)

- VarDict low frequency somatic filters: generalize strand and mismatch based
  filter based on cross-validation to avoid over filtering on high depth panels.
- strelka2 joint calling: switch to improved gvcfgenotyper approach for calling
  from gVCFs.
- Heterogeneity: initial support for PureCN and TitanCNA heterogeneity analysis
  including reporting on LOH in HLA for human samples. Work in progress validations:
  https://github.com/bcbio/bcbio_validations/tree/master/TCGA-heterogeneity
- CNV: initial support for GATK4 CNV calling as alternative to CNVkit for
  tumor normal analyses
- VarDict RNA-seq variant calling: avoid structural variants with recent vardict-java.
- RNA-seq variation: filter RNA-seq variants close to splice junctions,
  supporting STAR and hisat2.
- RNA-seq variation: add snpEff effects to output variant calls. Thanks to Manasa Surakala.
- RNA-seq: gzip/bgzip FASTQ files in `work/fastq` instead of the original directory.
- use biobambam2 BAM to FASTQ conversion instead of Picard in all cases.
- Trimming: add built-in support for adapters from the SMARTer Universal Low Input RNA Kit
(truseq2) and the Illumina NEXTera DNA prep kit from NEB (nextera2).
- ChIP/ATAC-seq: allow skipping duplicate marking.
- joint calling: ensure correct upload to final directory when no annotations present
- Logging: fix logging in parallel runs with new joblib loky backend. Thanks to
  Ben Liesfeld and Roland Ewald.

## 1.1.1 (6 November 2018)

- single-cell RNA-seq: add built-in support for 10x_v2.
- Fix UMI support for small RNA. Compatible with Qiagen UMI small RNA protocol.
- Ignore .Renviron when running Rscript to head-off PATH conflicts.
- Support SRR ids to download samples with bcbio_prepare_samples script.
- tumor-only prioritization: do not apply LowPriority filter by default, instead
  annotate with external databases. Use `tumoronly_germline_filter` to re-enable
  previous behavior.
- UMIs: apply default filtering based on de-duplicated read depth. Uses
`--min-reads 2` with raw de-duplicated coverage of 800 or more  or `--min-reads 1`
  otherwise. Allows error correction with UMIs for higher depth samples.
- gemini: databases no longer created by default. Use `tools_on: [gemini]` or
  `tools_on: [gemini_orig]` to create a database. We now use a reduced database
  for build 37 to match build 38 and make this forward compatible with CWL.
- vcfanno: run gemini and somatic annotations by default, producing annotated
  VCFs with external information.
- alignment preparation: support a list of split files from multiple sequencing
  lanes, merging into a single fastq
- variant: support octopus variant caller for germline and somatic samples.
- peddy: fix bug where not all files uploaded on first pipeline run
- peddy: For somatic analyses use separate germline calls for tumor/normal, if
  available, or extracted germline calls from supported callers, instead of
  somatic variants.
- GATK: support ploidy specification during joint calling.
- GATK BQSR: bin qualities into static groups (10, 20, 30) to match GATK4
  recommendations. Thanks to Severine Catreux.
- GATK: support 4.0.10.0 which does not use UCSC 2bit references for Spark tools
- variant calling: support bcftools 1.9 which is more strict about duplicated
  key names in INFO and FORMAT.
- seq2c: Upload global calls, coverage and read_mapping files to project directory.
- RNA-seq variant calling: Apply annotations after joint calling for GATK to
  avoid import errors with GenomicsDB. Thanks to Komal Rathi.
- CWL: add `--cwl` target to bcbio_nextgen.py upgrade to add and maintain bcbio-vm.
- CWL: use standard null instead of string "null" for representing None values.
- CWL: support for heterogeneity and structural variant callers that make
  use of variant inputs.
- CWL: support ensemble calling for combining multiple variant callers.
- ensemble: remove no-ALT ref calls that contribute to incorrect ensemble outputs
- RNA-seq: output a matrix of un-deduped UMI counts when doing single-cell/DGE
  for quality control purposes. This is called `tagcounts-dupes.mtx` in the
  final directory.
- single-cell RNA-seq: allow pre-transformed FASTQ files as input to DGE/single-cell
  pipeline.
- single-cell RNA-seq: only create one index per specified genome instead of per
  sample
- fgbio: back compatibility for older quality setting `--min-consensus-base-quality`
- RNA-seq: fix for `fusion_caller` getting interpreted as a path, leading to
  memoization/upload issues.
- RNA-seq: memoize rRNA quality calculations, speeding up reruns.
- RNA-seq: prefix `description` with an X if it starts with a number, for R
  compatibility.
  Thanks to Avinash Reddy and Dan Stetson at AstraZeneca.
- single-cell RNA-seq: respect `--positional` flag with the new tag counting. Thanks to
  Babak Alaei at AstraZeneca.
- RNA-seq: turn on `--seqBias` flag by default for Salmon as early-version overfitting
  issues have been fixed.
- RNA-seq: report insert size from Salmon fragment distribution, not samtools stats.
- RNA-seq: when processing explant samples, produce a combined tx2gene.csv file from
  all organisms processed.

## 1.1.0 (11 July 2018)

- Germline calls: rename outputs to `samplename-germline` to provide easier
  to understand outputs in final directory.
- Add bcbioRNASeq object creation and automatic quality report generation
  with `tools_on: [bcbiornaseq]`
- CWL: Support germline/somatic calling for tumor samples.
- CNVkit: improve whole genome runs. Better speed in normalize_sv_coverage
  through parallelization and avoiding logging. Avoid memory errors in segmentation.
- UMI: upload prepared UMI bam file (pre-consensus) to final output directory
- Add support for bbmap as an aligner
- RNA-seq variant calling: parallelize GATK HaplotypeCaller over regions to
  avoid memory and timeout issues.
- Support joint calling with GATK using pre-prepared gVCF inputs.
- RNA-seq variant calling: allow annotation of output variants with vcfanno
- Support hg38 builds with peddy QC
- QC: support VerifyBamID2 for contamination detection
- CWL: adjust defaults for align_split_size and nomap_split_targets to match
  different parallelization and overhead for these runs
- CWL: support for Cromwell runner
- custom genomes: Unzip GTF file prior to installation.
- Avoid making variant_regions required during processing (by filling with
  coverage) to differentiate targeted and non analyses downstream.
- Avoid attempts to download pre-installed S3 genomes, providing better
  errors with missing genome installs.
- Trimming: add explicit `polyg` option for removing 3' G stretches in NovaSeq
  and NextSeq data. Now defaults to no polyG trimming unless turned on.
- Chip-seq: Add RiP calculation for chip-seq data.
- DeepVariant and Strelka2 support for customized targeted/genome calling models
  per region to handle heterogeneous inputs.
- STAR: enable passing custom options for alignment.
- Add `tools_off: [coverage_qc]` option to skip calculating coverage stats (samtools-stats and picard).
- Adding BAM file for each sample in small-RNAseq pipeline, samtools
  and qualimap qc metrics to multiqc report.
- Allow arbitrary genomes for ChIP-seq. Thanks to @evchambers for pointing out the issue.

## 1.0.9 (10 April 2018)

- Use smoove for lumpy variant calling and genotyping, replacing custom lumpyexpress
  implementation: [validation](https://github.com/bcbio/bcbio_validations/tree/master/NA24385_sv#smoove-validation)
- Generalize exclusion of regions during variant calling with new
  `exclude_regions` target. Includes previously available LCR and high depth
  regions, in addition to removal of polyX and alternative contigs.
- Normalize allele frequency calculation and filtering for Strelka2 and MuTect2.
  Thanks to Vlad Saveliev.
- CNVkit: enable specification of pre-built reference background cnn with
  `background: cnv_reference`.
- CNVkit: handle projects with mixed CNVkit and non-CNVkit usage. Thanks to Luca
  Beltrame.
- Improved Atropos trimming: better use of multicore parallelization in variant
  and RNA-seq pipelines.
- Add support for polyG and polyX trimming to variant calling for NovaSeq 3' end
  cleanup and generally avoiding low complexity reads.
- Structural variant: use SURVIVOR for validation comparisons.
- RNA-seq variant calling: use multiple cores for VarDict.
- Support miRge2.0 for alternative small RNA annotation. Users should
  install the tool manually until compatible with bioconda.
- Add bamCoverage to chip-seq pipeline to calculate bigwig files.
- GATK4: Correctly use GATK4 GatherVcfs when tools_off: [gatk4] specified for
  variant calling. Thanks to Luca Beltrame.
- variant: Default to `mark_duplicates: false` if alignment turned off
  (`aligner: false`).
- variant: Fix race condition when preparing BED files for coverage and
  sv_regions. Thanks to Tristan Lubinski.
- Fix `noalt_calling` to correctly avoid parallelizing on non-standard
  chromosomes without a variant regions file.
- Fix broken `kraken` command. Thanks to @choosehappy.

## 1.0.8 (5 February 2018)

- GATK4 is the new default GATK release used in bcbio when running HaplotypeCaller or
  Base Quality Score Recalibration. Use `tools_off: [gatk4]` to use older GATK
  3.x versions.
- GATK4: Support 4.0 release with changed command line parameters. Re-enable
  multicore calling for CWL runs.
- GATK4: remove older GATK3 based gatk-framework in favor of equivalent GATK4
  commands.
- install: move to using recent IPython parallel to avoid dependency issues
- install: fix resolution issues due to conda 4.4.x (old ipython-cluster-helper,
  missing libquadmath.so with numpy due to libgcc update)
- RNA-seq variant calling: improve joint calling and parallelization with move
  to use GATK4 HaplotypeCaller.
- QC: improve read counting speed by moving to hts-nim-tools, replacing custom
  samtools view counting
- Add `noalt_calling` to avoid variant calling on non standard chromosomes.
  Thanks to Vlad Saveliev and Oliver Hofmann.
- variant alignment: improve core allocations for non-split alignments to avoid
  memory issues on 4Gb/core machines with whole genome samples.
- Add Total number of reads and adapter found to metrics in small RNA-seq pipeline.
- Add mirtop to the tools used in small RNA-seq pipeline for miRNA annotation.
- delly: Support 0.7.8 release which calls all variant types together.
- gVCF: fix basic filtering for GATK and sentieon when running without joint
  calling. Thanks to Tom Morris.

## 1.0.7 (6 January 2018)

- Automatically include bcbio anaconda PATH when running tools. Also allow
  custom BCBIOPATH specification to help with modules integration. Thanks to
  Gabriel Berriz.
- vcfanno: only correct VCF headers to use Number=1 when decomposition takes
  place. Avoids incorrect headers for non-decomposed inputs.
- ensemble: normalize and decompose variants prior to incorporating into
  ensemble calls, handling MNPs called differently across callers. Thanks to
  Vlad Saveliev.
- Avoid bgzipping and grabix indexing fastq inputs when not doing alignment
  splitting to save processing.
- Initial support for minimap2 aligner in variant calling workflows. Still needs
  validation and benchmarking in comparison to bwa.
- Standardize dbSNP annotation to use vcfanno for all variant callers. Remove
  GATK custom annotations for non-GATK callers, which are not present in GATK4.
- CNVkit: drop low coverage contaminating regions in tumor calls. Thanks to
  Eric Talevich.
- Expand `remove_extracontigs` for `bam_clean` to more consistently handle
  compatible pre-aligned BAMs with different extra contigs in reference genome.
- Fix problem collapsing samples for QC when using RNA-seq variant calling with
  gatk-haplotype. Thanks to Neill Gibson.
- Integrate ericscript RNA-seq fusion caller. Thanks to Tetiana Khotiainsteva
  and Vang Le.
- Remove read backed phasing (`phasing: gatk`) for GATK runs in favor of
  HaplotypeCaller internal phasing.
- disambiguation: ensure BAM index present for non-split alignments
- Use only end of reads to detect 3' adapters in small RNA-seq pipeline.
- Fix BCBIO_JAVA_HOME to correctly pass custom Java to GATK and Picard runs.
- ChIP-seq: add generation of greylist regions defined as regions of  high
  depth in the input file on a per sample pair basis.
- RNA-seq: STAR now outputs a MAPQ of 60 for uniquely mapped reads instead of
  255.
- RNA-seq: Ensure BAM files fed into Cufflinks have 255 as the uniquely mapped
  MAPQ instead of 60 as output by hisat2/STAR/etc.
- RNA-seq: omit duplicate files from stringtie assembly merging. Thanks to
  @mmoisse for the bug report.
- Add support for peddy (https://github.com/brentp/peddy) for PED file
  correspondence/ancestry checking.
- ChIP-seq: pass through encode filtered BAMs to upload directory.
- seq2c: pass through mapping_reads.txt file to directory.

## 1.0.6 (5 November 2017)

- Use mosdepth for callability calculations, replacing goleft depth. Centralize
  coverage and QC depth calculations around single mosdepth runs.
- Improve representation of germline and somatic calls in MultiQC report and
  output directory, avoiding confusing "-germline" extension. Thanks
  to Vlad Saveliev.
- Structural variants: return combined tumor/normal calls instead of single
  sample tumor for somatic calls in delly, lumpy, manta, and WHAM.
- VarDict: remove `-v 50` as required option for deep targeted panels (>5000x
  average coverage). Recommend adding if needed by a `var2vcf` resource options.
- Templating: avoid automatically setting flowcell date to maintain consistency
  between runs.
- Add `fusion_caller` as an optional algorithm field to turn on/off fusion
  callers. Currently supports oncofuse and pizzly.
- RNA-seq: better appropriate kmer size estimation for reads < 60 bp for
  Salmon/Rapmap/Sailfish index creation.
- RNA-seq variant calling: require gatk-haplotype instead of gatk as the caller.
- RNA-seq variant calling: support GATK4
- UMIs: move fgbio consensus calling to use filtering, adds `--max-reads` for
  high depth regions and swaps `--min-consensus-base-quality` for `--min-base-quality`
- Correctly re-bgzip fastq inputs even if not using `align_split_size`.
- Fix bug when running with `lumpy_usecnv` that resulted in skipping CNVkit.
- GATK gVCF joint calling: avoid running through bcftools for header fixes,
  using Picard instead. Avoids integer/double conversion incompatibilities.
- CWL: run variantcalling with multiple cores, reducing total jobs and enabling
  mulicore supporting callers.
- CWL: support structural variant calling as part of variant pipelines.
- Add pizzly (http://www.biorxiv.org/content/early/2017/07/20/166322)
  as a fusion caller when fusion mode is enabled.
- VEP: output an effect call per allele for multiallelic positions.
- Define separators for paired fastq files during bcbio_prepare_samples.py
- RNA-seq single-cell/DGE: add `transcriptome_gtf` as an option which will
  collapse single-cell/DGE counts down to the gene level. This is recommended
  for single-cell and DGE experiments.
- ChIP-seq: preliminary support for bwa for ChIP-seq alignment. Compared to bowtie2
  on a test dataset this results in a superset of the bowtie2 peaks, with 95% of the
  common peaks within 50 bases of each other. It calls about 50% more peaks
  though using the bwa alignments, use with care.

## 1.0.5 (11 Sept 2017)

- Add optional downsampling whole genome BAM files to a high maximum coverage
  (200 times the average coverage) to avoid slow runtimes in collapsed repeats
  and poly-ATGC regions. Downsampling happens in parallel with post alignment
  sorting. Currently turned off by default pending runtime improvements.
  Configure using `maxcov_downsample`.
- Separate post alignment recalibration and realignment. Recalibration now
  occurs multicore to support GATK4 implementation. We generally recommend
  skipping realignment.
- Provide multicore read trimming and streaming bgzip fastq output with atropos,
  replacing cutadapt as the default trimmer.
- hg38 runs do not run bwakit's bwa-postproc.js cleanup scripts unless HLA
  calling needed. Avoids slowdowns using this postprocessing script when running
  bwa with multiple cores.
- Tumor-only prioritization uses vcfanno output instead of GEMINI,
  allowing use without needing to build a full GEMINI database.
- Use samtools multicore indexing, replacing sambamba multicore index.
- Replace components of pipeline using single core sambamba view -c with
  parallel samtools equivalents.
- Replace sambamba depth coverage calculations with mosdepth to improve
  speed and parallelization.
- Multicore base quality score recalibration with GATK4 and Sentieon.
- GATK4: add support for gVCF based joint calling.
- GATK4: fix option usage for gVCF creation with HaplotypeCaller
- Allow overriding Java used in bcbio with `BCBIO_JAVA_HOME`
- Do not split individual sample VCFs during pooled batch calling. This
  previously happened only for small batches with less than 5 samples, now we
  avoid it entirely and let users do downstream sample extraction.
- Update OptiType HLA calling to use multicore CBC solver, also avoiding GLPK issues.
- Additional approach to retrieving cluster IP addresses for IPython and
  logging, using the fully qualified domain name.
- Add `archive: [cram-lossless]` to do CRAM archiving of outputs without quality
  score compression. Thanks to Alison Meynert.
- Add `tools_off: [lumpy-genotype]` option to skip Lumpy genotyping.
- CWL/WDL: use single file tarballs for complex collections of files like
  aligner, RTG and snpEff indices.
- GC bias correction is now the default for Salmon read-based quantification.
  See https://github.com/salmonteam/SalmonBlogResponse/blob/master/SalmonBlogResponse.md  for the reasoning behind this change.
- Add kallisto support for non single-cell RNA-seq experiments.
- Salmon can now be run alongside other RNA-seq quantifiers.
- Cufflinks and Stringtie can be run alongside each other as RNA-seq
  quantifiers.
- Check BED input files for coordinates off the ends of contigs.

## 1.0.4 (9 July 2017)

- Initial support for GATK4 variant calling with HaplotypeCaller and MuTect2.
  Requires `tools_on: [gatk4]` https://github.com/bcbio/bcbio_validations/tree/master/gatk4
- Enable adapter trimming for variant calling pipeline.
- Provide `trim_ends` command to quickly do defined end trimming as part of
  variant calling fastq preparation.
- Support duplex UMIs, present as embedded barcodes on read 1 and read 2.
- Sort region based analyses like variant calling by interval size. Ensures
  longest intervals run first avoiding delay at end of sample processing.
- Ensure FreeBayes dbSNP and GATK annotations passed into final file. Thanks
  to @semal.
- Use new Ensembl vep (variant effect predictor) with updated annotations.
  Thanks to Matthias De Smet.
- Accept files from HTTP/FTP as input
- CWL: use json input files for passing inputs instead of flattened command
  line arguments. Improves compatibility with multiple runners.
- Allow subsetting a pre-aligned BAM to only standard chromosomes, removing non
  chr1-22,X,Y for human. This allows runs of pre-aligned data with different
  extra chromosomes than the bcbio reference builds. Thanks to Oliver Hofmann.
- Improved support for pre-aligned BAMs by using contigs in BAM file for
  coverage calculations.
- Avoid grabix race conditions with multiple identical input files. Thanks to
  Andrey Tovchigrechko.
- Remove usage of lxml for qsignature and qualimap to avoid icu library errors.
- CNVkit: merge adjacent calls with identical copy numbers
- Add support for triple-barcoded cellular barcodes.
- Add support for Illumina's SureCell single-cell RNA-seq.

## 1.0.3 (7 May 2017)

- Allow installs to pull a specific git hash or tag revision of bcbio codebase.
- Fix FreeBayes somatic and multi-sample calling order to be consistent between
  chromosome region runs. Thanks to Ho Danliang.
- Fix structural variant output upload for complex batching cases. Correctly
  handle shared normals and other multi-batch by naming outputs using batches.
  Thanks to Sven-Eric Schelhorn.
- Move to samtools/bcftools/htslib 1.4. Provides parallel bgzip, removing need
  for pbgzip and improved concatenation speed for region split VCF files.
- Improve Lumpy prioritization speeds by adjusting location of breakend
  genotyping.
- UMI consensus: reduce runtimes to ~2/3 of previous avoiding unnecessary
  compression and file IO.
- UMI consensus: pass along metrics about consensus read generation as BAM tags
  in final file (cD = depth, cE = error rate)
- Support DNApi for de novo adapter detection in small RNA pipeline
- Several updates to the VarScan support: honor options specified in the
  resource config section; honor min_allele_frac option and set --strand-filter
  flag in the single-sample case; general cleanups. Thanks to Christian Brueffer.
- Update validation plots to support matplotlib 2.0.
- Enable mixed list/string inputs to germline calling. Thanks to Luca Beltrame.
- Fix qsignature outfile parsing. Thanks to Oliver Hofmann.
- Allow structural variant validations with VCF truth sets. Enables more
  flexible comparisons without size and event binning.
- Provide seq2c VCF output and enable validation of calls.
- Allow specification of seq2c options through resources. Thanks to Sally Luke
  and Marisa Cunha.
- Avoid using R_LIBS settings for R runs to limit incompatibilities with
  externally installed R packages.
- Provide absolute paths for relative paths to files in algorithm list inputs.
  Thanks to Matthias De Smet.
- Switch to Salmon from Sailfish as default alignment-free RNA-seq
  quantification algorithm.
- Add `sailfish` as a valid option for `expression_caller`.
- Fix chimeric alignment output option for STAR.
- Remove deprecated tidy counts for Sailfish/Salmon.
- Allow more possible empty/skip inputs in `variantcaller` and `svcaller`: None,
  null and empty lists
- Move DEXSeq to be an opt-in expression caller by default.
- Speed up combination of counts/RPKM/FPKM/TPM of samples into a single table by
  10x.

## 1.0.2 (7 March 2017)

- Fix FreeBayes paired somatic calling by generalizing support for finding
  non-ordered tumor/normal placement in VCF.
- Re-add checks for pre-bgzipped fastq inputs to alignment preparation thanks to
  a fix for grabix to handle Illumina bgzip outputs.
- Provide DNA damage annotation for low frequency sequencing errors in somatic
  samples. Use `tools_on: [damage_filter]`
- Add viral detection for variant calling DNA-seq cancer samples. Uses
  virus sequences from TCGA GDC distribution and provides simple counts of
  unmapped reads against viral sequences in MultiQC report.
- Improve lumpy structural variant runs from pre-aligned BAM files, using
  extract_sv_reads to avoid need to resort input files. Thanks to Neill Gibson.
- Move VCF files from SV prioritization to final upload directory. Thanks to
  Miika Ahdesmaki.
- Provide whole genome coverage plots with goleft indexcov. Thanks to Brent Pedersen.
- Speed up post-alignment callability calculations by using default parameters
  to goleft depth. Thanks to Brent Pedersen.
- Allow custom vcfanno configuration files for variant annotation and
  GEMINI database creation, using `vcfanno` configuration parameter. Optionally
  allows use of `vcfanno` without GEMINI database creation.
- Always use specified cores for analysis re-runs in local multicore mode.
  Avoids confusing core behavior with checkpoints on re-starts of analysis in
  a previous work directory.
- Upload of pipeline results to iRODS. Thanks to Matthias De Smet.
- Add duplicate removal to post-FreeBayes processing. Thanks to Neill Gibson.
- Support latest svtyper (0.1.1) for lumpy to provide speed improvements. Will
  default to 0.1.1 at next release.
  (use `bcbio_conda install -c bioconda svtyper=0.1.1` to test in development)
- Work towards supporting a Python 3 compatible bcbio codebase. Thanks to
  Michael Crusoe.
- Reduce VarDict maximum BED region sizes for better memory usage. Thanks to
  Nikolai Karulin, Oliver Hofmann, Miika Ahdesmaki and Zhongwu Lai.
- Require `tools_on: [lumpy_usecnv]` to pre-run CNVkit as input to Lumpy, allowing
  Lumpy and CNVkit to run in parallel otherwise.
- Avoid issues with CNVkit bin size estimates for normal associated with
  multiple tumors. Thanks to Ho Danliang.
- Fix double uploading of fast RNA-seq quantification.
- Output single-cell RNA-seq counts in annotated MatrixMarket format.

## 1.0.1 (17 January 2017)

- Fix bug in 1.0.0 release with parallel calculations on whole genome samples.
  The release version only parallelizes by chromosome instead of callable
  regions, resulting in less parallelism. Thanks to Sven-Eric Schelhorn and
  Neill Gibson.
- Generalize use of working directories to support runs on S3 mounted
  filesystems. Ensures all work takes place inside transactional directories.
  Thanks to Tetiana Khotiainsteva and Sven-Eric Schelhorn.
- Provide separate germline calling for somatic tumor/normal pairs. Supplements
  somatic calls with standard germline calls on normal samples, including
  ensemble and SV calling.
- Support creating GEMINI databases with new generic mechanism using vcfanno/vcf2db.
  This allows creation of GEMINI output for any organism. Adds support for hg38
  with annotations from dbSNP, Clinvar, ExAC and ESP.
- Support FreeBayes 1.1.0 for improved memory usage and 3-4x speedup.
  Will default to 1.1.0 at next release. Validation work:
  https://github.com/bcbio/bcbio.github.io/blob/master/_posts/2016-11-21-giab-hg38-freebayes.md
- Rework quality control for speed and output directory consistency. Avoid
  re-duplicating calculations and put all output in qc directory to make re-runs
  easier. Thanks to Vlad Saveliev.
- Fixes for Seq2C concurrency problems when preparing BED files. Thanks to Vlad
  Saveliev.
- Update WHAM structural variant caller to support the latest release.
- Update delly structural variant caller to support the latest release.
- Improve dbSNP annotation speeds for adding rs IDs to VarDict output.
  Thanks to Ben Liesfeld.
- Support for VEP 87 with additional plugins and generalization of fields.
  Thanks to Matthias De Smet.
- Deprecate `clinical_reporting` parameter and introduce new
  `effects_transcripts` parameter than enables more control over variant effects
  prediction. Enable HGVS by default for human projects and separates from
  transcript selection.
- For lumpy runs that use samblaster, use samtools sort instead of sambamba
  sort. Avoids segfault issues with samblaster. Thanks to Oliver Hofmann.
- Pre-install capture region BED files and enable short hand specification in
  sample configuration.
- Use vt normalize as part of GEMINI decomposition to clean up complex
  multiallelic variants. Thanks to Sergey Naumenko.
- Testing suite cleanup. Move to py.test and separate integration and unit
  tests. Thanks to Tetiana Khotiainsteva.
- Fix issue with cutadapt hanging on gzipped input. Thanks to Stephen Turner.
- Updated cutadapt to use single-pass trimming for paired-end files, improving
  performance and hitting the disk less.
- Added support for cellular barcode error correction with single-cell RNA-seq
  via the `cellular_barcode_correction` parameter. This corrects edit distances
  up to the set value, defaults to 1.
- Add support for sample-based demultiplexing of single-cell RNA-seq runs.
- Move single-cell RNA-seq results to the upload directory.
- Make positional UMI default to off for single-cell RNA-seq.
- Add support for the Klein lab v3 version of the inDrop protocol.

## 1.0.0 (20 November 2016)

- Default to no calling if `variantcaller` not specified, instead of old GATK
  UnifiedGenotyper default.
- Use samtools depth instead of bedtools genomecov for depth calculations, and
  calculate high depth regions during initial depth calculations.
  Improves speed by more than 6x. Thanks to Brent Pedersen.
- Adjust de-duplication strategy to use bamsormadup from biobambam2 for most
  cases and samblaster when split and discordant reads needed for SV calling
  with lumpy.
- Fix handling of fresh installs with GATK 3.6 only included. Correctly handles
  versioning from bioconda and lack of specifically defined jar directory.
- Unset JAVA_HOME when running gatk-framework and GATK > 3.6, forcing
  use of bcbio installed Java 8. Thanks to Brad Wubbenhorst.
- Fix bug when running realignment without recalibration in GATK 3.6. Thanks to Pär Larsson.
- Get from GEO server, GSM FASTQ samples using bcbio_prepare_samples.py script
- Add seqcluster stats to QC folder
- Allow manual specification of total memory and core usage for machines in
  `resources`. Thanks to Juan Caballero.
- Allow PED based gender specifications (1=male, 2=female). Thanks to Brent
  Pedersen.
- Annotate validation variants with genome context from GA4GH and other sources
  for interpreting true/false positives/negatives.
- Limit GATK cores used for GenotypeGVCFs to avoid excessive memory usage.
- VQSR: allow forcing GATK to try VQSR with tools_on. Generate VQSR plots.
  Thanks to Zhengqiu Cai.
- Support ATAC-seq for chipseq pipeline.
- Remove duplicates after alignment for chipseq pipeline.
- Support for bzip2 input files during variant calling. Thanks to Paulo Silva.
- Allow non-positional UMI Rapmap quantified single-cell RNA-seq.
- Re-enable save_diskspace option to reduce disk usage during alignment
  preparation and split alignments.
- Offload fixing the unmapped Tophat file to tophat-recondition. Thanks to Christian Brueffer.
- Add support for vcfanno (https://github.com/brentp/vcfanno) to annotate VCFs with other
  VCFs/BED files.
- Mark possible RNA-edits for GRCh37/hg19 using RADAR coordinates. Thanks to
  Sergey Naumenko for the suggestion.
- Add `local_controller` option to run the controller alongside the main bcbio
  process. Thanks to Brent Pederson and Sven-Eric Schelhorn.

## 0.9.9 (18 July 2016)

- Change defaults for recalibration and realignment to False. These have been
  the recommended settings (http://bit.ly/bcbio-minimal) and no realignment now
  matches Broad recommendations.
- Use conda installed Java instead of requiring external installation
  for most tools.
- Support GATK 3.6 with Java 8 installed as part of anaconda. Older GATK
  versions for calling and recalibration/realignment require external Java 7.
- Re-organization of variants stats using bcftools and
  cleaning gemini queries to get individual samples metrics.
- Quality control back end revamped to support better parallelization
  and pluggability of new QC metrics.
- Support CNV calling with Seq2C for exome, targeted or amplicon experiments.
  Thanks to Vlad Saveliev.
- Add `fixrg` target to `bam_clean` to accept BAM inputs with correct
  sorting and reads but that need an updated read group.
- More robust file transactions across network filesystems, avoiding failures
  from partially transferred files. Thanks to Sven-Eric Schelhorn.
- Improved checking of BAM files during merge steps. Thanks to Sven-Eric Schelhorn.
- Add SAMPLE and PEDIGREE tags to tumor/normal VCF outputs to enable
  easier post-analysis parsing of results.
- Add single point for annotation following variant calling to improve
  pluggability of new annotation types.
- Add support for running germline and somatic calling with Sentieon
  callers (https://peerj.com/preprints/1672/). Requires license from
  Sentieon.
- Fix fusion calling using Tophat2. Thanks to @csardas for raising the issue.
- Add support for kallisto quantification of single-cell RNA-seq data.
- Add `transcriptome_fasta` option to single-cell RNA-seq. This allows
  the user to provide a transcriptome FASTA file to quantitate against rather
  than use the `bcbio` provided annotation.
- Fix naming of vardict RNA-seq variant calls. Thanks to @csardas.

## 0.9.8 (20 May 2016)

- Correctly install all datatargets on new installation. Previously we'd
  skipped installing default additional data unless specified.
- Use yamllint to find wrong syntaxes in the YAML file that are ignored
  by pyyaml package and can affect the analysis.
- Improve choosing split regions for batch analysis to use the unionized
  intersection of non-callable regions. This enables better use of batches
  with different callable regions. Thanks to Neill Gibson.
- Fix HLA typing issues and handle HLA typing on split alignments.
  Thanks to Miika Ahdesmaki.
- Set `align_split_size` automatically based on input file sizes, trying to
  provide reasonable splits and avoid too many splits for large files.
- Fix high depth identification for whole genome runs, correctly calculating
  it when also inferring coverage estimations. Thanks to Neill Gibson.
- Do not remove duplicates for GATK variant calling when mark_duplicates
  is False or running amplicon sequencing.
- Fix installation of mutect jar via toolplus when mutect not previously
  present in configuration.
- Platypus: revert filtering back to defaults after additional cross-validation:
  http://i.imgur.com/szSo5M6.png
- Enable gVCF output with tools_on: [gvcf] for users who need gVCF output
  for downstream analyses.
- Avoid downscaling memory when recalibrating/realigning with GATK, since we
  should not longer need to work around Java issues. Thanks to Luca Beltrame.
- Do not use samblaster on genomes with greater than 32768 contigs, the
  samblaster maximum. Thanks to morten (@mattingsdal).
- Move to samtools for output CRAM support, using bamUtils for 8-bin compression
  of read quality scores.
- Remove `merge_bamprep` option and always merge realigned BAM files if run.
- Correctly clean up additional problem characters in sample descriptions that
  can confuse shell commands.

## 0.9.7 (29 March 2016)

- Use MultiQC (github.com/ewels/MultiQC) as main package to process all
  QC metrics.
- New install procedure for data: `--datatarget` allows installation of sub-sets
  of supplemental data for smaller installs for small RNA only analysis. Also
  provides a consistent framework for installing larger data types.
- VEP data no longer installed by default. Requires `--datatarget vep`
- During install, `--toolplus` only used for third party tools like GATK and
  MuTect and not data installation, which moved to `--datatarget`
- Provide `data_versions.csv` in output folder that has versions of reference
  data used in the analysis.
- Use sample description for BAM read group IDs, instead of lane index. This
  allows remixing of samples after processing without potential collisions. Thanks
  to Neill Gibson.
- Use sample description for file names instead of lane/flowcall information.
  Makes re-runs more stable when using template and files easier to interpret.
  Back compatible with re-runs of old work directories.
- Finalize support for MuTect2 with validation against the DREAM synthetic 4
  dataset (http://imgur.com/CLqJlNF). Thanks to Alessandro (@apastore).
- Do not bgzip inputs when they are already gzipped and do not require
  parallelization or format conversion. Thanks to Miika Ahdesmaki.
- Use new snpEff annotations (ANN) instead of older approach (EFF). The
  new annotations are more interoperable and supported by GEMINI.
- Lazy import of matplotlib libraries to avoid slow startup times.
- Only apply ploidyfix to all female batches to remove Y chromosome. Avoids
  confusion with file produced in other cases without any changes.
- Improvement to bcbio CWL integration: support parallel alignment and variant
  calling.
- Support for Salmon and RapMap added.
- FastRNA-seq pipeline implemented that does nothing but run Salmon with no QC.
- Singlecell RNA-seq pipeline implemented that uses https://github.com/vals/umis
  to handle the UMI and cellular barcode, aligns with RapMap and quantitates
  by counting, scaling ambiguous reads by the number of transcripts they could have
  come from.
- Migrate bowtie and bowtie2 to handle split input alignments, bgzipped inputs,
  and produce sorted, de-duplicated BAM files. This allows use in additional
  standard pipelines. Thanks to Luca Beltrame.
- Switch final upload directories for salmon and sailfish results to be of the
  form samplename/salmon instead of samplename/salmon/samplename.

## 0.9.6 (12 February 2016)

- Installation uses conda packages from bioconda for Python dependencies and
  third party tools.
- Add macs2 to chipseq pipeline.
- Add germline output files for somatic calling pipelines. The standard variant
  calls identify somatic mutations different from a normal, while the
  germline has pre-existing mutations which might contribute to cancer
  development.
- Use parallel bgzip for preparation of input fastq files for parallelization
  and alignment. Thanks to Guillermo Carrasco.
- Avoid extacting individual sample calls from pooled variant call runs for
  samples with more than 5 individuals in a batch. Avoids slow extraction run
  times. Thanks to Neill Gibson.
- Add explicit check for BED file mismatches with reference genome.
- During validation, report truth counts relative to initial truth set
  representation and pick best metric for plotting ROC scores.
- Remove `--sudo` flag from installer. bcbio requires install into a directory
  structure with user permissions.
- Add ability to tweak fastq preparation for alignment splitting so we can
  explore alternative approaches to bgzip and grabix index.
- Re-enable `stringtie` as an expression caller.
- Allow `stringtie` as a transcript assembler.
- Replace the `assemble_transcriptome` option with `transcript_assembler`, which
  accepts a list of assemblers to run. The output of all the assemblers is
  merged at the end with Cuffmerge.
- Move Picard to use conda installed `picard` single executable instead of
  custom installed java directory of jars.
- Add library type option to Cufflinks assembly. Thanks to Konstantin (@dezzan).
- Tag variants decomposed with vcfallelicprimitives. Thanks to Neill Gibson.
- Fix Platypus problem where we weren't correctly specifying BED regions since
  latest update skips over files not ending with".txt" or ".bed".

## 0.9.5 (12 December 2015)

- Add miRDeep2 to small RNA-seq analysis and quantify the novel miRNAs for
  all samples.
- Enable calling of HLA alleles with human build 38 (hg38). Turn on with the
  `hlacaller` option.
- Structural variant prioritization with BED files of known biologically
  important regions. Extracts SV calls in these regions and produces a tab
  delimited high level summary. Use the `svprioritize` option to enable.
- Add tRNA count and figures by tdrmapper for srna-seq pipeline.
- Avoid running callability checks on smaller chromosomes less than 1 million
  basepairs. Saves computation and disk IO on alt and support regions we don't
  split on.
- Enable nested batch specifications, allowing samples in partially overlapping
  batches.
- Speed improvements for Lumpy genotyping. Move to latest svtyper and avoid
  genotyping breakends.
- Allow use of VEP annotations on non-human analyses.
- Filter VarDict calls with poor mapping quality support (-Q 10) which
  trigger low frequency false positives.
- Remove ENCODE blacklist regions when calling with VarDict and FreeBayes on
  whole genomes. Avoids long run times due to collapsed repeats near centromeres.
- Update VarScan to 2.4.0 and rework support to allow piping between mpileup
  and VarScan to avoid filesystem IO.
- Annotate ensemble calls with information about supporting callers. Thanks to
  Pär Larsson and Son Pham.
- Move eXpress to expression_caller instead of being run by default.
- rRNA calculation uses the count file instead of using counts from GATK.
- Merge STAR fusion calls back into the BAM file. Thanks to Miika Ahdesmaki.
- Added preliminary support for the hisat2 aligner.
- Swapped STAR indexing to use on the fly splice junction indexing.
- Slightly inceased default DEXseq memory requirements in bcbio_system.yaml.
- Add support for RNA-seq for hg38 and hg38-noalt
- Make Sailfish the default for non-count based expression estimation.
  Produces isoform-level (combined.isoform.sf.tpm) and gene-level
  (combined.gene.sf.tpm) TPM expression estimation.
- Move Cufflinks to be off by default for expression estimation (turn on via
  expression_callers if needed).
- Add STAR fusion gene parameters suggested by @felixschlesinger.
- Add disambiguation to Sailfish by creating a master FASTA file of all
  transcripts from all organisms, quantitating each and separating out the
  organism-specific transcripts after.
- Add VarDict support for RNA-seq variant calling. Thanks to Miika Ahdesmaki and
  Sven-Eric Schelhorn.

## 0.9.4 (14 October 2015)

- Ensure genome data sort order is identical to BED files when annotating
  structural variant calls. Thanks To Miika Ahdesmaki.
- Improve low frequency calling for VarDict using vaidation against DREAM
  synthetic dataset 4.
- Install truth sets for germline and cancer calling automatically as part of
  bcbio and make it easy to include them in the configuration files for
  validation.
- Avoid need to set LD_LIBRARY_PATH and PERL5LIB on installations.
- Update Scalpel to latest version (0.5.1) and improve sensitivity for low
  frequency indels: http://imgur.com/a/7Dzd3
- Drop `coverage_depth_max` for downsampling, which no longer works in GATK 3.4.
  The option wasn't supported by other callers so was more confusing than useful.
- Fix missing BAM index when running with `align: false`. Thanks to Stephan
  Pabinger and Severine Catreux.
- Annotate structural variant files with snpEff. Initial steps towards
  summarized structural variant reporting.
- Add ability to specify platform unit (PU) and library (LB) in BAM header.
  Thanks to Brad Wubbenhorst.
- Update gatk-framework to 3.4-46 to avoid errors dealing with new gVCF output.
- Set java.io.tmpdir to avoid filling up global temporary space with snpEff.
  Thanks to Oliver Hofmann.
- Speed up transcriptome-only processing. Thanks to Sven-Eric Schelhorn.
- Add bamtools output to RNA-seq quality metrics. Thanks to Sven-Eric Schelhorn.
- Expand input quality format detection to detect full range of possible Sanger values.

## 0.9.3 (27 September 2015)

- Fix bug when using tumors with multiple normals and no CNV calling. Additional
  tumor sample would get lost due to lack of early (CNV-based) calling. Thanks
  to Miika Ahdesmaki.
- Include R and Rscript in the installation with conda packages and use for
  installing and running R-based tools. Avoids issues with alternative R
  versions and need for a separate installation.
- Fix bug when using CNVkit on disambiguated inputs. Thanks to Miika Ahdesmaki.
- Re-work structural variant infrastructure to provide plug-in parallel ensemble calling,
  removing the previous overlap-based ensemble calls. Currently supports MetaSV for
  ensemble calls. Also re-works validation to not rely on ensemble-overlap calls.
- Default to using Real Time Genomics vcfeval (https://github.com/RealTimeGenomics/rtg-tools)
  for validation instead of bcbio.variation. Improves speed and resolution of
  closely spaced variants. The old funtionality is still available with
  `validate_method: bcbio.variation`.
- Correctly apply BQSR when using recalibration with PrintReads by using GATK
  full instead of the open source GATK framework which silently ignores BQSR
  option. Thanks to Severine Catreux.
- Require larger blocks (250bp, moved from 100bp) to find regions for splitting analysis
  to avoid too tight splitting around small homozygous deletions.
- Adjust mapping quality (MQ) filter for GATK SNP hard filters to improve sensitivity
  http://imgur.com/a/oHRVB
- Ensure memory specification passed to sambamba and samtools sort during
  disambiguation and RNA-seq. Thanks to Sven-Eric Schelhorn.
- Fix compatbility with bedtools groupby in v2.25.0, which needs short
  parameters instead of long parameter names.
- Allow turning off variant quality score recalibration with `tools_off: [vqsr]`
- Generalize group size for batching gVCFs prior to joint calling with
  `joint_group_size`. Thanks to Severine Catreux.
- Support GEMINI 0.17.0, which does not have a --no-bcolz option since that is
the default.
- Remove test_run parameter since it was poorly supported and not used much.
- Fix issue with featureCounts sorting not working in parallel by pre-sorting
and filtering the BAM file.
- Unified stock coverage and experimental coverage reporting.
- Deprecated `report` and `coverage_experimental` as algorithm keys.

## 0.9.2 (1 September 2015)

- Support IPython 4.0 with ipyparallel
- Fix bug in writing BAM and VCF indexes to final directory. Correctly add
  indexes as bam.bai and vcf.gz.tbi.
- Fix bug in queryname sorting on split files for feeding into diambiguation.
  Ensure proper sorting with explicity sambamba sort. Thanks to Sven-Eric
  Schelhorn.
- Ensure extra FreeBayes alleles get removed prior to vcfallelicprimatives,
  avoiding leaving incorrect genotype allele fields. Thanks to Michael
  Schroeder.
- Split CNVkit processing into individual components, enabling better
  parallelization and control over parameters.
- Genotype Lumpy structural variant calls with SVtyper.
- Initial support for small RNA pipeline. Thanks to Lorena Pantano.
- Support for MetaSV to prepare combined structural variant calls.
- Add smallRNA-seq pipeline
- Test automatic report for variants calling and standard pipeline.
- Allow Cufflinks to be turned off via tools_off.

## 0.9.1 (6 August 2015)

- Fix novoalign to work with parallel split alignments. Thanks to Tyler Funnell.
- Move lumpy-sv to latest version which uses lumpyexpress instead of speedseq.
- Support the manta SV caller from Illumina. Validations: http://imgur.com/a/Gajsg
- Remove high depth regions from structural variant calling exclusion file
  to avoid false positives with lumpy. Thanks to Miika Ahdesmaki.
- Move some structural variant calling, like CNV detection, prior to variant
  calling. Allows use of CNV calls as inputs for variant detection tools.
- Generalize support for interaction with blob storage and graphing to support
  alternative cloud providers. Initial support for interacting with Azure.
  Thanks to Alexandru Coman.
- Remove VarDict call lines where reference and alternative allele are
  identical.
- Fix assignment issues during prioritization with new GEMINI and sqlite.
- Support updated versions of sambamba, which provide headers for window depth
  commands.

## 0.9.0 (20 June 2015)

- GATK 3.4: support HaplotypeCaller by avoiding setting downsampling (-dcov)
  option by default.
- Single sample structural variant calling: corectly handle multiple variant
  callers. Thanks to Sven-Eric Schelhorn.
- Make VarDictJava the default caller when `vardict` specified. `vardict-perl`
  is now required to specifically use the Perl version.
- VarDict and VarDictJava: limit regions to 1Mb with overlaps to avoid memory
  errors. Ignore regions without BED reads which can lead to large genomic
  sections and memory errors.
- VarDict and VarDictJava: annotate outputs with dbSNP.
- Add `tools_on` configuration with `svplots` option. This turns off structural
  variant plotting by default, which can be time consuming compared to calling.
- Add a `--only-metadata` argument to template preparation that will only
  include BAM or fastq files in sample YAML if they are present in the metadata
  CSV file.
- samblaster: support -M flag in 0.1.22 release
- Fix VEP/GEMINI incompatibility where empty fields are included in VCF output.
- VarDict: restrict maximum region size within a BED file to 2Mb to avoid high
  memory usage and failures for longer regions.
- Include snpEff effects summary file in output directory when used for effects
  prediction.

## 0.8.9 (10 May 2015)

- Upgrade variant effect predictor (VEP) to the latest Ensembl version (79) with
  support for hg38. The latest VEP has better support for multiple versions
  but incompatible database naming. This requires an update of tools and data in
  a two step process. First `bcbio_nextgen.py upgrade -u stable --tools`
  (or `-u development`) then `bcbio_nextgen.py upgrade --data`.
- Improve de-duplication for split alignments. Do not sort/merge during splits,
  and instead perform a global merge sort and de-duplication of the final set of
  reads.
- Initial support for new human genome build (hg38/GRCh38) including alternative
  alleles. Usage is in place but still requires validation and additional testing.
- Remove alternative alleles from downstream variant calling after using in alignment
  to avoid issues with chromosome names like `HLA*`.
- Enable installation of external conda-managed tools. Adds in builds for
  heterogeneity analysis.
- Clean up preparation process for multi-allelic inputs to GEMINI to avoid
  needing to split/merge. Thanks to Sven-Eric Schelhorn.

## 0.8.8 (29 April 2015)

- Automatically calculate `coverage_interval` based on coverage calculations,
  avoiding need to set this directly in input configuration.
- Update vt decompose to handle additional multi-allelic adjustments including
  all format attributes, providing full support for new GEMINI changes. Thanks
  to Brent Pedersen and Adrian Tan.
- Add `default` configuration target to `bcbio_system.yaml` reducing the need
  to set program specific arguments for everything.
- Ensure `resources` specified in input YAML get passed to global system
  configuration for making parallelization decisions. Thanks to Miika Ahdesmaki.
- Run upload process on distributed machines, allowing upload to S3 on AWS to take
  advantage of machines with multiple cores. Thanks to Lorena Pantano.
- Re-write interactions with external object stores like S3 to be more general
  and incorporate multiple regions and future support for non-S3 storage.
- Scale local jobs by total memory usage when memory constrains resource usage
  jinstead of cores. Thanks to Sven-Eric Schelhorn and Lorena Pantano.
- Disambiguation: improve parallelization by disambiguating on split alignment
  parts prior to merging. Thanks to Sven-Eric Schelhorn.
- Disambiguation: ensure ambiguous and other organism reads are sorted, merged
  and passed to final upload directory. Thanks to Sven-Eric Schelhorn.
- Fix problem with sambamba name sorting not being compatible with samtools.
  Thanks to Sven-Eric Schelhorn.
- FreeBayes: update to latest version (0.9.21-7) with validation
  (http://imgur.com/a/ancGz).
- Allow bz2 files in bcbio_prepare_sample.py script.
- Ensure GEMINI statistics run for project summary file. Thanks to Luca
  Beltrame.
- Better error checking for booleans in input configuration. Thanks to Daryl
  Waggott.
- Implement qualimap for RNAseq QC metrics, but not active yet.
- collect statistics graphing capabilities moved from bcbio-nextgen-vm, enabling
  plotting of resource usage during runs. Thanks to John Morrissey and Lorena
  Pantano.

## 0.8.7 (12 March 2015)

- Run snpEff 4.1 in back-compatibility mode to work with GEMINI database
  loading. Fixes snpEff 4.1/GEMINI effects loading.
- Add PED file to GEMINI database load, containing family, gender and phenotype
  information from bcbio metadata. Thanks to Luca Beltrame and Roy Ronen.
- Enable specification of input PED files into template creation, extracting
  family, gender and phenotype information. Any sample rows from PED files get
  used when creating the GEMINI database.
- Fix preparation of multi-allelic inputs to GEMINI by implementing custom merge
  of bi-allelic and split multi-allelic. Previous implementation using GATK
  CombineVariants re-merged some split multi-allelic, losing effects annotations.
- Skip contig order naming checking with bedtools 2.23.0+ to avoid potential
  issues with complex naming schemes.
- Installation and upgrade: Set pip SSL certificates to point at installed conda
  SSL package if present. Avoids SSL errors when pip can't find system
  certificates. Thanks to Andrew Oler.
- Enable support for PBSPro schedulers through ipython-cluster-helper.

## 0.8.6 (23 February 2015)

- Calculate high depth regions with more than 20x median coverage as targets for
  filtering in structural variants. Attempts to detect and avoid spurious calls
  in repetitive regions.
- Support snpEff 4.1, including re-download of snpEff databases on demand if out
  of sync with older versions.
- Split multi-allelic variants into bi-allelic calls prior to loading into
  GEMINI, since it only handles bi-allelic inputs. Thanks to Pär Larsson.
- Pass ploidy to GATK HaplotypeCaller, supporting multiple ploidies and correct
  calling of X/Y/MT chromosomes. Requires GATK 3.3.
- Remove extra 'none' sample when calling tumor-only samples using
  MuTect. Harmonizes headers with other tumor-only callers and enables
  tumor-only ensemble calling. Thanks to Miika Ahdesmaki.
- Perform variant prioritization as part of tumor-only calling, using population
  based frequencies like 1000 genomes and ExAC and presence in known disease
  causing databases like COSMIC and Clinvar.
- Switch to samtools sort from sambamba sort during alignment streaming. Saves
  steps in processing and conversions on single sample no deduplication inputs.
- On AWS, download inputs for S3 instead of streaming into fastq preparation to
  avoid issues with converting BAM to fasta. Thanks to Roy Ronen.
- Provide better defaults for mincores that packs together multiple single IPython
  processes on a single cluster request -- use core specification from input
  configuration. Thanks to Miika Ahdesmaki.

## 0.8.5 (11 January 2015)

- No longer keep INFO fields with `vcfallelicprimitves` in FreeBayes,
  Platypus and Scalpel calling to prevent introduction of problematic
  fields for multi-allelic MNPs.
- Fix batching problem when using `coverage` and multiple shared batches
  like a global normal in cancer calling. Thanks to Luca Beltrame.
- Use `mincores` specification to ipython-cluster-helper to combine single core
  jobs into a single submission job for better memory shared on resource
  constrained systems.
- Move disambiguation split work inside parallel framework so download and
  preparation occurs on worker nodes or inside Docker containers. Enables on
  demand download of disambiguation genomes.
- Ensure population databases created when some inputs do not have variant calls.
- Switch to seaborn as matplotlib wrapper, from prettplotlib.
- Fixes for ensemble structural variant calling on single samples.
- Fixes for mixing joint and pooled calling in a single configuration file.
- Support for qSNP for tumor-normal calling.
- Add eXpress to RNA-seq pipeline.
- Add transcriptome-only mapping with STAR, bowtie2 or bwa.
- Change logging time stamps to be UTC and set explicitly as ISO 8601 compliant
  output. Improves benchmarking analysis and comparability across runs.
- Add support for RNA-seq variant calling with HaplotypeCaller
- Fix parallelization of DEXSeq.

## 0.8.4 (29 November 2014)

- Improvements in VarDict calling on somatic samples.
- Fix compatibility issue with bedtools 2.22.0 when calculating genome coverage.
- Fix joint calling upload to avoid redundant inclusion of full VCF file in
  individual sample directories.
- Fixes for inclusion of GATK jars inside Docker contains when running
  distributed jobs.
- Enable generation of STAR indexes on demand to handle running STAR on AWS
  instances.
- Re-organize code to prepare samples and reference genomes so it runs inside
  distributed processing components. This isolates process to Docker containers
  on AWS and also enables complex operations like preparing reference genomes on
  demand.

## 0.8.3 (19 November 2014)

- Improve tumor/normal calling with FreeBayes, MuTect, VarDict and VarScan by
  validating against DREAM synthetic 3 data.
- Validate ensemble based calling for somatic analysis using multiple callers.
- Improve ability to run on Amazon AWS, including up to date interaction with
  files originally stored in S3 and transfer to S3 on completion with encryption.
- Avoid race conditions during `bedprep` work on samples with shared input BED
  files. These are now processed sequentially on a single machine to avoid
  conflicts. Thanks to Justin Johnson.
- Add data checks and improved flexibility when specifying
  joint callers. Thanks to Luca Beltrame.
- Default to a reduced number of split regions (`nomap_split_targets` defaults
  to 200 instead of 2000) to avoid controller memory issues with large sample
  sizes.
- Avoid re-calculating depth metrics when running post variant calling
  annotation with GATK to provide accurate metrics on high depth samples.
  Thanks to Miika Ahdesmaki.
- Consistently keep annotations and genotype information for split MNPs from
  vcfallelicprimitives. Thanks to Pär Larsson.
- Enable VQSR for large batches of exome samples (50 or more together) to
  coincide with joint calling availability for large populations.
- Support retrieval of GATK and MuTect jars from S3 to enable integration
  with bcbio inside Docker.
- Bump pybedtools version to avoid potential open file handle issues. Thanks to
  Ryan Dale.
- Move to bgzipped and indexes human_ancestor.fa for LOFTEE to support access
  with new samtools that no longer uses razip.

## 0.8.2 (September 17, 2014)

- Fix bug in creating shared regions for analysis when using a single sample in
  multiple batches: for instance, when using a single normal sample for multiple
  tumors. Thanks to Miika Ahdesmaki.
- Unify approach to creating temporary directories. Allows specification of a
  global temporary directory in `resources: tmp:` used for all
  transactions. This enables full use of local temporary space during
  processing, with results transferred to the shared filesystem on completion.
- Fix issues with concatenating files that fail to work with GATK's
  CatVariants. Fall back to bcftools concat which correctly handles problem
  headers and overlapping segments.
- Enable flexible specification of `indelcaller` for `variantcaller` targets
  that do not have integrated indel methods. Thanks to Miika Ahdesmaki.
- Move to samtools 1.0 release. Update samtools variant calling to support new
  multiallelic approach.
- Improve Platypus integration: correctly pass multiple BAM files, make use of
  assembler, split MNPs, and correctly restrict to variant regions.
- Be more aggressive with system memory usage to try and make better use of
  available resources. The hope is to take advantage of Java memory fixes that
  previously forced us to be conservative.

## 0.8.1 (August 29, 2014)
- Support joint recalling with GATK HapolotypeCaller, FreeBayes and Platypus. The
  `jointcaller` configuration variable enables calling concurrently in large
  populations by independently calling on samples them combining into a final
  combined callset with no-call/reference calls at any position called
  independently.
- Add qsignature tool to standard and variant analyses, which helps identify
  sample swaps. Add `mixup_check` configuration variant to enable.
- Fix issue with merging GATK produced VCF files with vcfcat by swapping to
  GATK's CatVariants. Thanks to Matt De Both.
- Initial support for ensemble calling on cancer tumor/normal calling. Now
  available for initial validation work. Thanks to Miika Ahdesmaki.
- Enable structural variant analyses on shared batches (two tumors with same
  normal). Thanks to Miika Ahdesmaki.
- Avoid Java out of memory errors for large numbers of running processes by
  avoiding Parallel GC collction. Thanks to Justin Johnson and Miika Ahdesmaki.
- Enable streaming S3 input to RNA-seq and variant processing. BAM and fastq
  inputs can stream directly into alignment and trimming steps.
- Speed improvements for re-running samples with large numbers of samples or
  regions.
- Improved cluster cleanup by providing better error handling and removal of
  controllers and engines in additional failure cases.
- Support variant calling for organisms without dbSNP files. Thanks to Mark Rose.
- Support the SNAP aligner, which provides improved speed on systems with
  larger amount of memory (64Gb for human genome alignment).
- Support the Platypus haplotype based variant caller for germline samples with
  both batched and joint calling.
- Fix GATK version detection when `_JAVA_OPTIONS` specified. Thanks to Miika
  Ahdesmaki.
- Use msgpack for ipython serialization to reduce message sizes and IPython
  controller memory instead of homemade json/zlib approach.

## 0.8.0 (July 28, 2014)

- Change defaults for installation: do not use sudo default and require
  `--sudo` flag for installing system packages. No longer includes default
  genomes or aligners to enable more minimal installations. Users install
  genomes by specifically enumerating them on the command line.
- Add support for Ensembl variant effects predictor (VEP). Enables annotation
  of variants with dbNSFP and LOFTEE. Thanks to Daniel MacArthur for VEP
  suggestion.
- Support CADD annotations through new GEMINI database creation support.
- Rework parallelization during variant calling to enable additional multicore
  parallelization for effects prediction with VEP and backfilling/squaring off
  with bcbio-variation-recall.
- Rework calculation of callable regions to use bedtools/pybedtools thanks to
  groupby tricks from Aaron Quinlan. Improves speed and memory usage for
  coverage calculations. Use local temporary directories for
  pybedtools to avoid filling global temporary space.
- Improve parallel region generation to avoid large numbers of segments on
  organisms with many chromosomes.
- Initial support for tumor normal calling with VarDict. Thanks to
  Miika Ahdesmaki and Zhongwu Lai.
- Provide optional support for compressing messages on large IPython jobs to
  reduce memory usage. Enable by adding `compress_msg` to `alogrithm` section of
  `bcbio_system.yaml`. There will be additional testing in future releases
  before making the default, and this may be replaced by new methods like
  transit (https://github.com/cognitect/transit-python).
- Add de-duplication support back for pre-aligned input files. Thanks to
  Severine Catreux.
- Generalize SGE support to handle additional system setups. Thanks to Karl Gutwin.
- Add reference guided transcriptome assembly with Cufflinks along with functions
  to classify novel transcripts as protein coding or not as well as generally clean
  the Cufflinks assembly of low quality transcripts.
- Developer: provide datadict.py with encapsulation functions for looking up and
  setting items in the data dictionary.
- Unit tests fixed. Unit test data moved to external repository:
  https://github.com/roryk/bcbio-nextgen-test-data
- Add exon-level counting with DEXseq.
- Bugfix: Fix for Tophat setting the PI flag as inner-distance-size and not insert size.
- Added kraken support for contamination detection (@lpatano):
  http://ccb.jhu.edu/software/kraken/
- Isoform-level FPKM combined output file generated (@klrl262)
- Use shared conda repository for tricky to install Python packages:
  https://github.com/chapmanb/bcbio-conda
- Added initial chanjo integration for coverage calculation (@kern3020):
  https://github.com/robinandeer/chanjo
- Initial support for automated evaluation of structural variant calling.
- Bugfix: set library-type properly for Cufflinks runs.
- Added `genome_setup.py` a script to prepare your own genome and rnaseq files.

## 0.7.9 (May 19, 2014)

- Redo Illumina sequencer integration to be up to date with current
  code base. Uses external bcl2fastq demultiplexing and new bcbio integrated
  analysis server. Provide documentation on setting up automated infrastructure.
- Perform de-duplication of BAM files as part of streaming alignment process
  using samblaster or biobambam's bammarkduplicates. Removes need for secondary
  split of files and BAM preparation unless recalibration and realignment
  needed. Enables pre-processing of input files for structural variant detection.
- Rework batched regional analysis in variant calling to remove custom cases and
  simplify structure. Filtering now happens explicitly on the combined batch
  file. This is functionally equivalent to previous filters but now the workflow
  is clearer. Avoids special cases for tumor/normal inputs.
- Perform regional splitting of samples grouped by batch instead of globally,
  enabling multiple organisms and experiments within a single input sample YAML.
- Add temporary directory usage to enable use of local high speed scratch disk
  on setups with large enough global temporary storage.
- Update FreeBayes to latest version and provide improved filtering for high
  depth artifacts.
- Update VQSR support for GATK to be up to date with latest best
  practices. Re-organize GATK and filtering to be more modular to help with
  transition to GATK 3.x gVCF approaches.
- Support CRAM files as input to pipeline, including retrieval of reads from
  defined sequence regions.
- Support export of alignment data as CRAM instead of BAM for space storage
  and long term archiving.
- Provide configuration option, `remove_lcr`, to filter out variants in low
  complexity regions.
- Improve Galaxy upload for LIMS supports: enable upload of FastQC as PDF
  reports with wkhtmltopdf installed. Provide tabular summaries of mapped reads.
- Improve checks for pre-aligned BAMs: ensure correct sample names and
  provide more context on errors around mismatching reference genomes.
- GATK HaplotypeCaller: ensure genotype depth annotation with DepthPerSampleHC
  annotation. Enable GATK 3.1 hardware specific optimizations.
- Use bgzipped VCFs for dbSNP, Cosmic and other resources to save disk
  space. Upgrade to Cosmic v68.
- Avoid VCF concatenation errors when first input file is empty. Thanks to
  Jiantao Shi.
- Added preliminary support for oncofuse for calling gene fusion events. Thanks
  to @tanglingfung.

## 0.7.8 (March 21, 2014)

- Add a check for mis-specified FASTQ format in the sample YAML file. Thanks
  to Alla Bushoy.
- Updated RNA-seq integration tests to have more specific tags (singleend, Tophat,
  STAR, explant).
- Fix contig ordering after Tophat alignment which was preventing GATK-based
  tools from running.
- Allow calculation of RPKM on more deeply sampled genes by setting
  `--max-bundle-frags` to 2,000,000. Thanks to Miika Ahdesmaki.
- Provide cleaner installation process for non-distributable tools like
  GATK. The `--tooplus` argument now handles jars from the GATK site or Appistry
  and correctly updates manifest version information.
- Use bgzipped/tabix indexed variant files throughout pipeline instead of raw
  uncompressed VCFs. Reduces space requirements and enables parallelization on
  non-shared filesystems or temporary space by avoiding transferring
  uncompressed outputs.
- Reduce memory usage during post-alignment BAM preparation steps (PrintReads
  downsampling, deduplication and realignment prep) to avoid reaching memory cap
  on limited systems like SLURM. Do not include for IndelRealigner which needs
  memory in high depth regions.
- Provide explicit targets for coverage depth (`coverage_depth_max` and
  `coverage_depth_min`) instead of `coverage_depth` enumeration. Provide
  downsampling of reads to max depth during post-alignment preparation to avoid
  repetitive centromere regions with high depth.
- Ensure read group information correctly supplied with bwa aln. Thanks to Miika
  Ahdesmaki.
- Fix bug in retrieval of snpEff databases on install. Thanks to Matan Hofree.
- Fix bug in normal BAM preparation for tumor/normal variant calling. Thanks to
  Miika Ahdesmaki.
- General removal of GATK for variant manipulation functionality to help focus
  on support for upcoming GATK 3.0. Use bcftools for splitting of variants into
  SNPs and indels instead of GATK. Use vcflib's vcfintersection to combine SNPs
  and indels instead of GATK. Use bcftools for sample selection from
  multi-sample VCFs. Use pysam for calculation of sample coverage.
- Use GATK 3.0 MIT licensed framework for remaining BAM and variant manipulation
  code (PrintReads, CombineVariants) to provide one consistent up to date set of
  functionality for GATK variant manipulation.
- Normalize input variant_regions BED files to avoid overlapping
  segments. Avoids out of order errors with FreeBayes caller which will call in
  each region without flattening the input BED.

## 0.7.7 (February 27, 2014)

- For cancer tumor/normal calling, attach final call information of both to
  the tumor sample. This provides a single downstream file for processing and
  analysis.
- Enable batch specification in metadata to be a list, allowing a single normal
  BAM file to serve as a control for multiple tumor files.
- Re-organization of parallel framework code to enable alternative approaches.
  Document plugging in new parallel frameworks. Does not expose changes to users
  but makes the code cleaner for developers.
- Default to 1Gb/core memory usage when not specified in any programs. Do not
  use default baseline if supplied in input file. Thanks to James Porter.
- Integrate plotting of variant evaluation results using prettyplotlib.
- Add `globals` option to configuration to avoid needing to specify the same
  shared file multiple times in a samples configuration.
- Remove deprecated Celery distributed messaging, replaced in favor of IPython.
- Remove algorithm/custom_algorithm from bcbio_system.yaml, preferring to set
  these directly in the sample YAML files.
- Remove outdated and unused custom B-run trimming.
- Remove ability to guess fastq files from directories with no specification in
  sample YAML. Prefer using generalized template functionality with explicit
  specification of files in sample YAML file.
- Remove deprecated multiplex support, which is outdated and not
  maintained. Prefer approaches in external tools upstream of bcbio-nextgen.
- Add `--tag` argument which labels job names on a cluster to help distinguish
  when multiple bcbio jobs run concurrently. Thanks to Jason Corneveaux.
- Connect min_read_length parameter with read_through trimming in
  RNA-seq. Thanks to James Porter.
- Map `variant` calling specification to `variant2` since original approach
  no longer supported.
- Fix issues with trying to upload directories to Galaxy. Thanks to Jim Peden.
- Made inner distance calculation for Tophat more accurate.
- Added gffutils GFF database to the RNA-seq indices.
- Add gene name annotation from the GFF file instead of from mygene.

## 0.7.6 (January 15, 2014)

- Expand template functionality to provide additional ability to add metadata
  to samples with input CSV. Includes customization of algorithm section and
  better matching of samples using input file names. Improve ability to
  distinguish fastq pairs.
- Generalize snpEff database preparation to use individual databases located
  with each genome. Enables better multi-organism support.
- Enable tumor/normal paired called with FreeBayes. Contributed by Luca Beltrame.
- Provide additional parallelization of bgzip preparation, performing grabix indexing
  in parallel for paired ends.
- Fix downsampling with GATK-lite 2.3.9 releases by moving to sambamba based downsampling.
  Thanks to Przemek Lyszkiewicz.
- Handle Illumina format input files for bwa-mem alignment, and cleanly convert
  these when preparing bgzipped inputs for parallel alignment. Thanks to Miika
  Ahdesmaki.
- Provide better algorithm for distinguishing bwa-mem and bwa-aln usage. Now
  does random sampling of first 2 million reads instead of taking the first set
  of reads which may be non-generalizable. Also lowers requirement to use
  bwa-mem to 75% of reads being smaller than 70bp. Thanks to Paul Tang.
- Enable specification of a GATK key file in the bcbio_system resources
  `keyfile` parameter. Disables callbacks to GATK tracking. Thanks to Severine
  Catreux for keyfile to debug with.
- Correctly handle preparation of pre-aligned BAM files when sorting and
  coordinate specification needed. Thanks to Severine Catreux.
- Fix incorrect quality flag being passed to Tophat. Thanks to Miika Ahdesmaki.
- Fix Tophat not respecting the existing --transcriptome-index. Thanks to Miika Ahdesmaki.
- Keep original gzipped fastq files. Thanks again to Miika Ahdesmaki.
- Fixed incompatibility with complexity calculation and IPython.
- Added strand-specific RNA-seq support via the strandedness option.
- Added Cufflinks support.
- Set stranded flag properly in htseq-count. Thanks to Miika Ahdesmaki.
- Fix to ensure Tophat receives a minimum of 8 gb of memory, regardless of number of cores.
- Remove `hybrid_bait` and `hybrid_target` which were no longer used with new
  lightweight QC framework. Prefer better coverage framework moving forward.
- Added extra summary information to the project-summary.yaml file so downstream tools can
  locate what genome resources were used.
- Added ``test_run`` option to the sample configuration file. Set it to True to run a small
subset of your data through the pipeline to make sure everything is working okay.
- Fusion support added by setting ``fusion_mode: True`` in the algorithim section.
Not officially documented for now until we can come up with best practices for it.
- STAR support re-enabled.
- Fixed issue with the complexity calculation throwing an exception when there
were not enough reads.
- Add disambiguation stats to final project-summary.yaml file. Thanks to Miika Ahdesmaki.
- Remove `Estimated Library Size` and `Complexity` from RNA-seq QC
summary information as they were confusing and unnecessarily alarming,
respectively. Thanks to Miika Ahdesmaki and Sara Dempster.
- Several memory allocation errors resulting in jobs getting killed in
cluster environments for overusing their memory limit fixed.
- Added JVM options by default to Picard to allocate enough memory
for large BAM->FastQ conversion.

## 0.7.5 (November 29, 2013)

- Update overall project metrics summary to move to a flexible YAML format that
  handles multiple analysis types. Re-include target, duplication and variant
  metrics.
- Support disambiguation of mixed samples for RNA-seq pipelines. Handles alignment
  to two genomes, running disambiguation and continuation of disambiguated samples
  through the pipeline. Contributed by Miika Ahdesmaki and AstraZenenca.
- Handle specification of sex in metadata and correctly call X,Y and
  mitochondrial chromosomes.
- Fix issues with open file handles for large population runs. Ensure ZeroMQ contexts
  are closed and enable extension of ulimit soft file and user process limits within
  user available hard limits.
- Avoid calling in regions with excessively deep coverage. Reduces variant calling
  bottlenecks in repetitive regions with 25,000 or more reads.
- Improve `bcbio_nextgen.py upgrade` function to be more consistent on handling of
  code, tools and data. Now each require an implicit specification, while other
  options are remembered. Thanks to Jakub Nowacki.
- Generalize retrieval of RNA-seq resources (GTF files, transcriptome indexes) to use
  genome-resources.yaml. Updates all genome resources files. Contributed by James Porter.
- Use sambamba for indexing, which allows multicore indexing to speed up index
  creation on large BAM processing. Falls back to samtools index if not available.
- Remove custom Picard metrics runs and pdf generation. Eliminates dependencies on
  pdflatex and R for QC metrics.
- Improve memory handling by providing fallbacks during common memory intensive steps.
  Better handle memory on SLURM by explicitly allowing system memory in addition
  to that required for processing.
- Update fastqc runs to use a BAM files downsampled to 10 million reads to avoid
  excessive run times. Part of general speed up of QC step.
- Add Qualimap to generate plots and metrics for BAM alignments. Off by default
  due to speed issues.
- Improve handling of GATK version detection, including support for Appistry versions.
- Allow interruption of read_through trimming with Ctrl-C.
- Improve test suite: use system configuration instead of requiring test specific setup.
  Install and use a local version of nose using the installer provided Python.
- Fix for crash with single-end reads in read_through trimming.
- Added a library complexity calculation for RNA-seq libraries as a QC metric
- Added sorting via sambamba. Internally bcbio-nextgen now inspects the headers
  of SAM/BAM files to find their sorting status, so make sure tools set it correctly.

## 0.7.4 (October 20, 2013)

- Framework for indexing input reads using parallel bgzip and grabix, to handle
  distributed alignment. Enables further distribution of alignment step beyond
  multicore nodes.
- Rework of ensemble calling approach to generalize to population level ensemble
  calls. Provide improved defaults for handle 3 caller consolidation.
- Support for Mouse (mm10) variant calling and RNA-seq.
- For recent versions of Gemini (0.6.3+) do not load filtered variants into
  database, only including passed variants.
- Improve specification of resource parameters, using multiple `-r` flags
  instead of single semi-colon separated input. Allow specification of pename
  resource parameter for selecting correct SGE environment when not
  automatically found.
- Support biobambam's bammarkduplicates2 for duplicate removal.
- Clean up logging handling code to be more resilient to interrupt messages.
- Speed improvements for selecting unanalyzed and unmapped reads to address
  bottlenecks during BAM prep phase.
- Bug fix for algorithm options incorrectly expanded to paths on re-runs. Thanks
  to Brent Pedersen for report.
- Fix for Tophat 2.0.9 support: remove reads with empty read names.
- Save installation and upgrade details to enable cleaner upgrades without
  needing to respecify genomes, tool directory and other options from
  installation.

## 0.7.3 (September 22, 2013)

- Move specification of supporting genome files for variation (dbSNP, training
  files) and RNA-seq (transcript GTF files) analyses into an organism specific
  resources file. Improves ability to support additional organisms and genome
  builds.
- Provide paired tumor/normal variant calling with VarScan. Thanks to Luca Beltrame.
- Require bash shell and use of pipefail for piped commands. Ensures rapid
  detection of failures during piped steps like alignment.
- Use samtools cat for post-BAM merging to avoid issues with bamtools
  requirement for open file handles.
- Add installation/upgrade options to enable commercially restricted and data
  intensive third party tools.
- Support for GATK 2.7
- Fixes for TopHat 2.0.9 support: remove extra non-mate match paired end reads
  from alignment output.
- Pull `description` sample names from BAM files if not present in input
  configuration file. Thanks to Paul Tang for suggestion.
- Bug fixes for non-paired RNA-seq analysis.
- Add custom filtration of FreeBayes samples using bcbio.variation.
- Default to phred33 format for Tophat alignment if none specified.

## 0.7.2 (August 30, 2013)

- Report memory usage for processes to cluster schedulers and use predicted
  memory usage to schedule cores per machine. Gets core and memory information
  for machines and uses to ensure submitted jobs can schedule with available
  resources.
- Provide error checking of input YAML configuration at run start. Avoids
  accidental typos or incorrect settings that won't error out until later in the
  process.
- Drop requirement for fc_name and fc_date in input YAML file. Individual sample
  names are instead used and required to be unique within a processing run.
- Remove original `variant` pipeline, replacing with the all around better
  `variant2` analysis method. Plan for the next version is to automatically
  redirect to `variant2`.
- Improve parallelization of BAM preparation and gemini database creation by
  moving to multicore versions.
- Move variant annotation to work on called sub-regions, to avoid bottlenecks
  when annotating a full whole genome VCF.
- Remove sequencer-specific integration functionality which is poorly maintained
  and better done with third party tools: demultiplexing and statistics from
  Illumina directories.
- Bug fix to re-enable template generation functionality.
- Improve BAM merging on large files using samtools for output sort.
- Uploading results works with the RNA-seq pipeline.
- Rework internals to provide a consistent dictionary of sample attributes up
  front, avoiding lane/sample dichotomy which provided confusing internal code.
- Drop calling htseq-count from the command line in favor of an internal
  implementation.

## 0.7.1 (August 12, 2013)

- Remove requirement for bcbio_system.yaml passed in on command line, defaulting
  to default file prepared by installer unless specified.
- Bug fixes for new approach to parsing *.loc files: handle Galaxy *.loc files
  with mixed tabs and spaces correctly and fall back to previous approaches
  when aligner specific *.loc files are missing.
- Bug fixes for preparing merged BAM files using bamtools: correctly sort after
  merging and avoid duplication of reads in noanalysis files.
- Bug fix for concatenating files when first file in empty.
- Recover from ZeroMQ logging errors, avoiding loss of logging output.

## 0.7.0 (July 30, 2013)

- RNA-seq pipeline updated: deprecate Tophat 1 in favor of Tophat 2. Perform
  automatic adapter trimming of common adapter sequences. STAR aligner support.
  RNA-SeQC support for RNA-seq specific quality control. Transcript quantitation
  with htseq-count.
- Updated installation and upgrade procedures, to make it easier to build an
  initial analysis pipeline and upgrade bcbio-nextgen and third-parts tools and
  data in place.
- Add support for MuTect tumor/normal variant caller, contributed by Luca
  Beltrame.
- Generalize variant calling to support alternative callers like cancer-specific
  calling: provide additional associated files to variant calls and pass along
  sample specific metadata. Document implementation of new variant callers.
- Improve algorithms around post-variant calling preparation. Avoid unnecessary
  tries for VQSR on low coverage whole genome reads, and concatenate VCF files to
  avoid locking penalties.
- Fix logging and memory usage for multicore jobs run within ipython clusters.
- Improve logging for IPython cluster issues, including moving IPython logs
  inside project logging directory for better access.
- Options for improved cluster resiliency: minimize number of clusters started
  during processing with more extensive reuse, flexible timeouts for waiting on
  cluster start up, and expose options to allow job retries. Thanks to Zhengqiu
  Cai for suggestions and testing.

## 0.6.5 (June 07, 2013)

- Improve logging: Detailed debugging logs collect all process standard out and
  error and command lines across distributed systems.
- Piping improvements: provide fully piped analysis with GATK recalibration and
  gkno realignment. Handle smaller reads with novoalign piped analysis.
- Improve collapsing analysis regions into evenly sized blocks to better handle
  large numbers of samples analyzed together.
- Provide template functionality to ease generation of input sample.yaml files
  from lists of BAM of fastq files. Thanks to Brent Pedersen and Paul Tang.
- Updated program support: Improved novoalign support based on evaluation with
  reference genomes. Support GATK 2.5-2. Support VarScan 2.3.5.
- Fix naming of read group information (ID and SM) to be more robust. Identifies
  issues with duplicated read groups up front to avoid downstream errors during
  variant calling. Thanks to Zhengqiu Cai.
- Improve quality control metrics: Cleanup into custom qc directory and ensure
  correct selection of duplicate and other metrics for split post-alignment
  prep, even without merging.
- Fix IPython parallel usage for larger clusters, providing improved resiliency
  for long running jobs.
- Clean up handling of missing programs and input files with better error
  messages. From Brent Pedersen.

## 0.6.4 (May 06, 2013)

- Integrate fully with bcbio.variation to provide automated validation of
  variant calls against reference materials.
- Provide full list of all third party software versions used in analysis.
- Create GEMINI database as part of output process, allowing immediate queries
  of variants with associated population and annotation data.
- Collapse analysis regions into evenly sized blocks separated by non-callable
  regions. Provides better parallelism.
- Documentation and examples for NA12878 exome and whole genome pipelines.
