Current development goals and outstanding tasks for bcbio-nextgen development.
These are roughly ordered by current priority and we welcome contributors.

- Improve summary information. Leverage existing tools like Qualimap for
  alignment metrics and remove custom scripts to call Picard and produce output
  PDF. Remove dependency on pdflatex. ([#140][140]).

[140]: https://github.com/chapmanb/bcbio-nextgen/issues/140

- Integrated structural variant analysis, including CNV prediction. Current
  targets are [lumpy][lumpy] and [delly][delly].

[delly]: http://www.embl.de/~rausch/delly.html
[lumpy]: https://github.com/arq5x/lumpy-sv

- Improved support for cancer tumor/normal paired callers. Suggested callers
  include SomaticSniper ([#66][66], [#109][109]), LoFreq and others. A
  comprehensive discussion is at [#112][112]. Requires improved framework for
  evaluating callers and approaches for handling Ensemble calling with multiple
  inputs ([#67][67]).

[66]: https://github.com/chapmanb/bcbio-nextgen/issues/66
[67]: https://github.com/chapmanb/bcbio-nextgen/issues/67
[109]: https://github.com/chapmanb/bcbio-nextgen/issues/109
[112]: https://github.com/chapmanb/bcbio-nextgen/issues/112

- Explore options for accumulating and displaying summary information from
  multiple runs. Prioritize options which allow accumulation across multiple
  analysis machines and already handle query and visualization.

- Work on GATK-free pipeline with full piping of intermediates during
  preparation. Evaluate against GATK best-practice using NA12878 reference data.

- Performance improvements and testing on Amazon EC2. Make use of high speed
  local ephemeral storage for temporary space.

- Implement and evaluate GATK's ReducedReads for use in large scale variant
  calling projects ([#90][90]).

[90]: https://github.com/chapmanb/bcbio-nextgen/issues/90

- Correctly handle haplotype chromosomes (chrM, chrY) and setup calling based on
  sample sex ([#33][33]).

[33]: https://github.com/chapmanb/bcbio-nextgen/issues/33

