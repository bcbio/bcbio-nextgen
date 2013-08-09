## 0.7.1 (in progress)

- Bug fixes for new approach to parsing *.loc files: handle Galaxy *.loc files
  with mixed tabs and spaces correctly and fall back to previous approaches
  when aligner specific *.loc files are missing.
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
