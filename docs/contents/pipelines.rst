Pipelines
---------

Germline variant calling
~~~~~~~~~~~~~~~~~~~~~~~~

bcbio implements configurable SNP, indel and structural variant calling for
germline populations. We include whole genome and exome evaluations against
reference calls from the `Genome in a Bottle`_ consortium and `Illumina Platinum
Genomes <http://www.illumina.com/platinumgenomes/>`_ project, enabling continuous
assessment of new alignment and variant calling algorithms. We regularly report
on these comparisons and continue to improve approaches as the community makes
new tools available. Here is some of the research that contributes to the
current implementation:

- An introduction to the `variant evaluation framework`_. This includes a
  comparison of the `bwa mem`_ and `novoalign`_ aligners. We also compared the
  `FreeBayes`_, `GATK HaplotypeCaller`_ and `GATK UnifiedGenotyper`_ variant
  callers.

- An in-depth evaluation of `FreeBayes and BAM post-alignment processing`_. We
  found that FreeBayes quality was equal to GATK HaplotypeCaller. Additionally,
  a lightweight post-alignment preparation method using only de-duplication was
  equivalent to GATK's recommended Base Quality Score Recalibration (BQSR) and
  realignment around indels, when using good quality input datasets and callers
  that do local realignment.

- Additional work to `improve variant filtering`_, providing methods to
  remove low complexity regions (LCRs) that can bias indel results. We also
  tuned `GATK's Variant Quality Score Recalibrator`_ (VQSR) and compared it with
  hard filtering. VQSR requires a large number of variants and we use
  it in bcbio with GATK HaplotypeCaller when your :ref:`algorithm-config`
  contains high depth samples (``coverage_depth`` is not low) and you are
  calling on the whole genome (``coverage_interval`` is genome) or have more
  than 50 regional or exome samples called concurrently.

- An `evaluation of joint calling`_ with GATK HaplotypeCaller, FreeBayes,
  Platypus and samtools. This validates the joint calling implementation,
  allowing scaling of large population germline experiments. It also
  demonstrates improved performance of new callers: samtools 1.0 and Platypus.

- Support for `build 38 of the human genome
  <http://bcb.io/2015/09/17/hg38-validation/>`_, improving precision of
  detection thanks to the improved genome representation.

bcbio automates post-variant calling annotation to make
the outputs easier to feed directly into your biological analysis. We annotate
variant effects using `snpEff`_ or `Variant Effect Predictor`_ (VEP), and
prepare a `GEMINI database`_ that associates variants with multiple
external annotations in a SQL-based query interface.

.. _Genome in a Bottle: http://www.genomeinabottle.org/
.. _variant evaluation framework: https://bcb.io/2013/05/06/framework-for-evaluating-variant-detection-methods-comparison-of-aligners-and-callers/
.. _FreeBayes and BAM post-alignment processing: https://bcb.io/2013/10/21/updated-comparison-of-variant-detection-methods-ensemble-freebayes-and-minimal-bam-preparation-pipelines/
.. _improve variant filtering: http://bcb.io/2014/05/12/wgs-trio-variant-evaluation/
.. _FreeBayes: https://github.com/ekg/freebayes
.. _GATK UnifiedGenotyper: http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html
.. _GATK HaplotypeCaller: http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html
.. _samtools mpileup: http://samtools.sourceforge.net/mpileup.shtml
.. _GATK's Variant Quality Score Recalibrator: http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html
.. _bwa mem: http://bio-bwa.sourceforge.net/
.. _novoalign: http://www.novocraft.com
.. _snpEff: http://snpeff.sourceforge.net/
.. _GEMINI database: http://gemini.readthedocs.org/en/latest/
.. _Variant Effect Predictor: http://www.ensembl.org/info/docs/tools/vep/index.html
.. _evaluation of joint calling: http://bcb.io/2014/10/07/joint-calling/

Basic germline calling
======================

The best approach to build a bcbio :ref:`docs-config` for germline calling is to use
the :ref:`automated-sample-config` with one of the default templates:

- `FreeBayes template
  <https://github.com/chapmanb/bcbio-nextgen/blob/master/config/templates/freebayes-variant.yaml>`_ --
  Call variants using FreeBayes with a minimal preparation pipeline. This is a
  freely available unrestricted pipeline fully included in the bcbio installation.

- `GATK HaplotypeCaller template
  <https://github.com/chapmanb/bcbio-nextgen/blob/master/config/templates/gatk-variant.yaml>`_ --
  Run GATK best practices, including Base Quality Score Recalibration,
  realignment and HaplotypeCaller variant calling. This requires a license from
  Broad for commercial use. You need to manually install GATK along with bcbio
  using downloads from the GATK Broad site or Appistry (see :ref:`extra-install`).

You may also want to enable :ref:`svs-pipeline` for detection of larger events,
which work with either caller. Another good source of inspiration are the
configuration files from the :ref:`example-pipelines`, which may help identify
other configuration variables of interest. A more complex setup with multiple
callers and resolution of ensemble calls is generally only useful with a small
population where you are especially concerned about sensitivity. Single
caller detection with FreeBayes or GATK HaplotypeCaller provide good resolution
of events.

Population calling
==================

When calling multiple samples, we recommend calling together to provide improved
sensitivity and a fully squared off final callset. To associate samples together
in a population add a ``metadata`` ``batch`` to the :ref:`sample-configuration`::

    - description: Sample1
      metadata:
        batch: Batch1
    - description: Sample2
      metadata:
        batch: Batch1

Batching samples results in output VCFs and GEMINI databases containing
all merged sample calls. bcbio has two methods to call samples together:

- Batch or pooled calling -- This calls all samples simultaneously by feeding
  them to the variant caller. This works for smaller batch sizes (< 50 samples)
  as memory requirements become limiting in larger pools. This is the default
  approach taken when you specify a ``variantcaller`` in the
  :ref:`variant-config` configuration.

- Joint calling -- This calls samples independently, then combines them together
  into a single callset by integrating the individual calls. This scales to
  larger population sizes by avoiding the computational bottlenecks of pooled
  calling. We recommend joint calling with HaplotypeCaller if you have a
  license for GATK usage, but also support joint calling with FreeBayes using a
  custom implementation. Specifying a ``jointcaller`` along with the appropriate
  ``variantcaller`` in the :ref:`variant-config` configuration enables this::

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

Cancer variant calling
~~~~~~~~~~~~~~~~~~~~~~
bcbio supports somatic cancer calling with tumor and optionally matched normal pairs using
multiple SNP, indel and structural variant callers. A `full evaluation of cancer calling`_
validates callers against `synthetic dataset 3 from the ICGC-TCGA DREAM challenge`_.
bcbio uses a majority voting ensemble approach to combining calls from
multiple SNP and indel callers, and also flattens structural variant calls into a
combined representation.

The `example configuration <https://github.com/chapmanb/bcbio-nextgen/blob/master/config/examples/cancer-dream-syn3.yaml>`_
for the :ref:`example-cancer` validation is a good starting point for setting up
a tumor/normal run on your own dataset. The configuration works similarly to
population based calling. Supply a consistent batch for tumor/normal pairs and
mark them with the phenotype::

    - description: your-tumor
      metadata:
        batch: batch1
        phenotype: tumor
    - description: your-normal
      metadata:
        batch: batch1
        phenotype: normal

Other :ref:`config-cancer` configuration options allow tweaking of the
processing parameters.

We're actively working on improving calling to better account for the
heterogeneity and structural variability that define cancer genomes.

.. _full evaluation of cancer calling: http://bcb.io/2015/03/05/cancerval/
.. _synthetic dataset 3 from the ICGC-TCGA DREAM challenge: https://www.synapse.org/#!Synapse:syn312572/wiki/62018

.. _svs-pipeline:

Structural variant calling
~~~~~~~~~~~~~~~~~~~~~~~~~~
bcbio can detect larger structural variants like deletions, insertions, inversions
and copy number changes for both germline population and cancer variant calling,
based on validation against existing truth sets:

- `Validation of germline structural variant detection`_ using multiple calling methods
  to validate against deletions in NA12878. This implements a pipeline that
  works in tandem with SNP and indel calling to detect larger structural
  variations like deletions, duplications, inversions and copy number variants
  (CNVs).

- `Validation of tumor/normal calling <http://bcb.io/2015/03/05/cancerval/>`_
  using the synthetic DREAM validation set. This includes validation of
  additional callers against duplications, insertions and inversions.

To enable structural variant calling, specify ``svcaller`` options in the
algorithm section of your configuration::

    - description: Sample
      algorithm:
        svcaller: [lumpy, manta, cnvkit]

The best supported callers are `Lumpy <https://github.com/arq5x/lumpy-sv>`_ and
`Manta <https://github.com/Illumina/manta>`_, for paired end and split read
calling, `CNVkit <http://cnvkit.readthedocs.org/en/latest/>`_ for read-depth
based CNV calling, and `WHAM <https://github.com/jewmanchue/wham>`_ for
association testing. We also support `DELLY
<https://github.com/tobiasrausch/delly>`_, another excellent paired end and
split read calling, although it is slow on large whole genome datasets.

In addition to results from individual callers, bcbio can create a summarized
ensemble callset using `MetaSV <https://github.com/bioinform/metasv>`_. We're
actively working on improved structural variant reporting to highlight potential
variants of interest.

.. _Validation of germline structural variant detection: http://bcb.io/2014/08/12/validated-whole-genome-structural-variation-detection-using-multiple-callers/

RNA-seq
~~~~~~~

bcbio-nextgen also implements a configurable best-practice pipeline for RNA-seq
quality control, adapter trimming, alignment and post-alignment quantitation

- Adapter trimming:
  - `cutadapt`_

- Sequence alignment:
  - `tophat2`_
  - `STAR`_
  - `hisat2`_

- Quality control:
  - `qualimap`_
  - `FastQC`_

- Quantitation:
  - `featureCounts`_
  - `DEXSeq`_
  - `Sailfish`_

After a run you will have in the ``upload`` directory a directory for each
sample which contains a BAM file of the aligned and unaligned reads, a
``Sailfish`` directory with the output of Sailfish, including TPM values,
and a ``qc`` directory with plots from FastQC and qualimap.

In addition to directories for each sample, in the ``upload`` directory there is
a project directory which contains a YAML file describing some summary statistics
about each sample and some provenance data. In that directory is also a
``combined.counts`` file which can be used as a starting point for performing
differential expression calling using any count-based method such as EdgeR,
DESeq2 or voom+limma, etc.

smallRNA-seq
~~~~~~~~~~~~

bcbio-nextgen also implements a configurable best-practices pipeline for smallRNA-seq
quality controls, adapter trimming, miRNA/isomiR quantification and other small RNA
detection.

- Adapter trimming:
  - `cutadapt`_

- Sequence alignment:
  - `STAR`_ for genome annotation
  - bowtie, `bowtie2` and  hisat2 for genome annotation as an option
  - `seqbuster <https://github.com/lpantano/seqbuster>`_ for miRNA annotation

- Quality control:
  - `FastQC`_

- Other small RNAs:
  - `seqcluster <https://github.com/lpantano/seqcluster>`_

The pipeline generates a RMD template file that can be rendered with knitr.
An example of the report can be seen `here <https://github.com/lpantano/mypubs/blob/master/srnaseq/mirqc/ready_report.md>`_.

ChIP-seq
~~~~~~~~
bcbio-nextgen implements the first steps of a ChIP-seq analysis up to aligning with
bowtie2. It doesn't do anything other than get the samples into a state
where a peak caller like MACS2 can be used.

- Adapter trimming:
  - `cutadapt`_

- Sequence alignment:
  - `bowtie2`_

- Quality control:
  - `FastQC`_

Standard
~~~~~~~~

This pipeline implements ``alignment`` and ``qc`` tools. Furthermore, it will
run `qsignature`_ to detect possible duplicated samples, or miss-labeling. It
uses SNPs signature to create a distance matrix that helps easily to create
groups. The project yaml file will show number of total samples analyzed, number
of very similar samples, and samples that could be duplicated.

.. _qsignature: http://sourceforge.net/p/adamajava/wiki/qSignature/

Configuration
=============
We will assume that you installed bcbio-nextgen with the automated installer,
and so your default `bcbio_system.yaml`_ file is configured correctly with all
of the tools pointing to the right places. If that is the case, to run
bcbio-nextgen on a set of samples you just need to set up a YAML file that
describes your samples and what you would like to do to them. Let's say that you
have a single paired-end control lane, prepared with the Illumina `TruSeq`_ Kit
from a human. Here is what a well-formed sample YAML file for that RNA-seq
experiment would look like::

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

``fc_date`` and ``fc_name`` will be combined to form a prefix to name
intermediate files, you can set them to whatever you like.  ``upload`` is
explained pretty well in the `configuration documentation`_ and the above will
direct bcbio-nextgen to put the output files from the pipeine into the ``final``
directory.  Under ``details`` is a list of sections each describing a sample to
process.  You can set many `parameters`_ under each section but most of
the time just setting a few like the above is all that is necessary.
``analysis`` tells bcbio-nextgen to run the best-practice RNA-seq pipeline on
this sample.

In the above, since there are two files, ``control_1.fastq`` and
``control_2.fastq`` will be automatically run as paired-end data. If you have
single end data you can just supply one file and it will run as single-end. The
``description`` field will be used to eventually rename the files, so make it
very evocative since you will be looking at it a lot later. ``genome_build`` is
self-explanatory.

Sometimes you need a little bit more flexibility than the standard pipeline, and
the ``algorithm`` section has many options to fine-tune the behavior of the
algorithm. ``quality_format`` tells bcbio-nextgen what quality format your FASTQ
inputs are using, if your samples were sequenced any time past 2009 or so, you
probably want to set it to ``Standard``. Adapter read-through is a problem in
RNA-seq libraries, so we want to trim off possible adapter sequences on the ends
of reads, so ``trim_reads`` is set to ``read_through``, which will also trim off
poor quality ends. Since your library is a RNA-seq library prepared with the
TruSeq kit, the set of adapters to trim off are the TruSeq adapters and possible
polyA tails, so ``adapters`` is set to the both of those. ``strandedness``
can be set if your library was prepared in a strand-specific manner and can
be set to firststrand, secondstrand or unstranded (the default).

Multiple samples
================
Lets say you have a set of mouse samples to analyze and each sample is a single
lane of single-end RNA-seq reads prepared using the NextEra kit.  There are
two case and two control samples. Here is a
sample configuration file for that analysis::

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

More samples are added just by adding more entries under the details section.
This is tedious and error prone to do by hand, so there is an automated
`template_` system for common experiments. You could set up the previous
experiment by making a mouse version of the `illumina-rnaseq`_ template
file and saving it to a local file such as ``illumina-mouse-rnaseq.yaml``. Then
you can set up the sample file using the templating system::

    bcbio_nextgen.py -w template illumina-mouse-rnaseq.yaml mouse_analysis
    /full/path/to/control_rep1.fastq /full/path/to/control_rep2.fastq
    /full/path/to/case_rep1.fastq /full/path/to/case_rep2.fastq


If you had paired-end samples instead of single-end samples, you can still use
the template system as long as the forward and reverse read filenames are
the same, barring a _1 and _2. For example: control_1.fastq and control_2.fastq
will be detected as paired and combined in the YAML file output by the
templating system.


.. _bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
.. _tophat2: http://tophat.cbcb.umd.edu/
.. _STAR: http://code.google.com/p/rna-star/
.. _cutadapt: http://cutadapt.readthedocs.org/en/latest/guide.html
.. _qualimap: http://qualimap.bioinfo.cipf.es
.. _FastQC: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _HTSeq: http://www-huber.embl.de/users/anders/HTSeq/doc/index.html
.. _TruSeq: http://www.illumina.com/products/truseq_rna_sample_prep_kit_v2.ilmn
.. _bcbio_system.yaml: http://github.com/chapmanb/bcbio-nextgen/blob/master/config/bcbio_system.yaml
.. _configuration documentation: http://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html#upload
.. _parameters: http://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html
.. _template: http://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html#automated-sample-configuration
.. _illumina-rnaseq: http://raw.github.com/chapmanb/bcbio-nextgen/master/config/templates/illumina-rnaseq.yaml
.. _eXpress: http://bio.math.berkeley.edu/eXpress/overview.html
.. _featureCounts: http://bioinf.wehi.edu.au/featureCounts/
.. _DEXSeq: https://bioconductor.org/packages/release/bioc/html/DEXSeq.html
.. _Sailfish: https://github.com/kingsfordgroup/sailfish
.. _hisat2: https://ccb.jhu.edu/software/hisat2/index.shtml
