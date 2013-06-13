## 0.7.0 (in development)
- RNA-seq pipeline updated.

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
