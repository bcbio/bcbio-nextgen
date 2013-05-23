## 0.6.5 (in progress)

- Piping improvements: provide fully piped analysis with GATK recalibration and
  gkno realignment. Handle smaller reads with novoalign piped analysis.
- Updated program support: Improved novoalign support based on evaluation with
  reference genomes. Support GATK 2.5-2. Support VarScan 2.3.5.
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
