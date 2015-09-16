Current development goals and outstanding tasks for bcbio-nextgen development.
These are roughly ordered by current priority and we welcome contributors.

- Improved deployment experience using [docker][docker] containers to provide a
  fully isolated bcbio-nextgen installation. Requires re-working of installation
  process to be a two step process: download docker + add external biological
  data.  Also requires adjustment of the pipeline and distributed processing to
  involve starting and using code isolated inside docker container. Work in
  progress is at [bcbio-nextgen-vm][bcbio-nextgen-vm].

[docker]: http://www.docker.io/
[bcbio-nextgen-vm]: https://github.com/chapmanb/bcbio-nextgen-vm

- Enable processing on Amazon EC2 with use of spot instances and no shared
  filesystem. Store file intermediates in S3 object storage instead of
  globally shared filesystem and make use of high speed local ephemeral storage.

- Integrated structural variant analysis, including CNV prediction. Current
  targets are [lumpy][lumpy], [delly][delly] and [cn.mops][cn.mops].

[cn.mops]: http://www.bioconductor.org/packages/release/bioc/html/cn.mops.html
[lumpy]: https://github.com/arq5x/lumpy-sv
[delly]: https://github.com/tobiasrausch/delly

- Improved support for cancer tumor/normal paired callers. Suggested callers
  include SomaticSniper ([#66][66], [#109][109]), LoFreq and others. A
  comprehensive discussion is at [#112][112]. FreeBayes supports tumor/normal
  calling: see [this mailing list discussion][fb-somatic] for the suggested
  parameters.  Requires improved framework for evaluating callers and approaches
  for handling Ensemble calling with multiple inputs ([#67][67]).

[66]: https://github.com/chapmanb/bcbio-nextgen/issues/66
[67]: https://github.com/chapmanb/bcbio-nextgen/issues/67
[109]: https://github.com/chapmanb/bcbio-nextgen/issues/109
[112]: https://github.com/chapmanb/bcbio-nextgen/issues/112
[fb-somatic]: https://groups.google.com/d/msg/freebayes/beLYRuHMkQE/RwFMniDmBYoJ

- Improve analysis of coverage, especially in targeted sequencing
  experiments. Plan to integrate with [chanjo]. See [#249][249] for more
  discussion.

[chanjo]: https://github.com/robinandeer/chanjo
[249]: https://github.com/chapmanb/bcbio-nextgen/issues/249

- Support [gVCF and incremental join discovery][gatk3-ijd] approach for calling
  variants. Switches batch approaches to calling independently, then combining
  in a final step. Also integrate [bcbio.variation.recall] for performing in
  non-GATK, non-gVCF scenarios.

[gatk3-ijd]: http://gatkforums.broadinstitute.org/discussion/3896/the-gatk-reference-model-pipeline-for-incremental-joint-discovery-in-full-detail

- Explore options for accumulating and displaying summary information from
  multiple runs. Prioritize options which allow accumulation across multiple
  analysis machines and already handle query and visualization.

- Once initial structural variation analysis and evaluation is in place,
  incorporate and evaluate additional CNV and structural variant callers. Some
  current targets are the [VarScan2 CNV caller][vs2] and [Control-FREEC][cfc].

[cfc]: http://bioinfo-out.curie.fr/projects/freec/
[vs2]: http://varscan.sourceforge.net/copy-number-calling.html

- Document and expand [Ensemble calling][ensemble] functionality with work on
  speed ups and parallelization. Integrate development work on
  [bcbio.variation.recall] using recalling with local realignment.

[ensemble]: http://bcb.io/2013/10/21/updated-comparison-of-variant-detection-methods-ensemble-freebayes-and-minimal-bam-preparation-pipelines/
[bcbio.variation.recall]: https://github.com/chapmanb/bcbio.variation.recall

- Add in methylation analysis approaches. See
  [[#618][https://github.com/chapmanb/bcbio-nextgen/issues/618]]
  for discussion.

- Handle split inputs across multiple sequencing lanes, handling merging of
  multiple fastq/BAM inputs and correctly maintaining lane information in BAM
  read group headers.

- Test to see if [less strict quality trimming][quality] results in better RNA-seq DE results.

- Evaluate RNA-seq fusion analysis callers and implement support for one if we can find one with
  reliable results ([#210][210]).

[210]: https://github.com/chapmanb/bcbio-nextgen/issues/210
[quality]: http://biorxiv.org/content/early/2013/12/23/000422
