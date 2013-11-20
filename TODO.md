Current development goals and outstanding tasks for bcbio-nextgen development.
These are roughly ordered by current priority and we welcome contributors.

- Integrated structural variant analysis, including CNV prediction. Current
  targets are [lumpy][lumpy] and [cn.mops][cn.mops].

[cn.mops]: http://www.bioconductor.org/packages/release/bioc/html/cn.mops.html
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

- Once initial structural variation analysis and evaluation is in place,
  incorporate and evaluate additional CNV and structural variant callers. Some
  current targets are the [VarScan2 CNV caller][vs2] and [Control-FREEC][cfc].

[cfc]: http://bioinfo-out.curie.fr/projects/freec/
[vs2]: http://varscan.sourceforge.net/copy-number-calling.html

- Explore options for accumulating and displaying summary information from
  multiple runs. Prioritize options which allow accumulation across multiple
  analysis machines and already handle query and visualization.

- Performance improvements and testing on Amazon EC2. Make use of high speed
  local ephemeral storage for temporary space.

- Implement and evaluate GATK's ReducedReads for use in large scale variant
  calling projects ([#90][90]).

[90]: https://github.com/chapmanb/bcbio-nextgen/issues/90

