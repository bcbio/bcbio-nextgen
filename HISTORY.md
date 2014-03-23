## 0.7.9 (in progress)

- GATK HaplotypeCaller: ensure genotype depth annotation with DepthPerSampleHC
  annotation. Enable GATK 3.1 hardware specific optimizations.
- Use bgzipped VCFs for dbSNP, Cosmic and other resources to save disk
  space. Upgrade to Cosmic v68.

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
