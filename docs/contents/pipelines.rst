Pipelines
---------

Germline variant calling
~~~~~~~~~~~~~~~~~~~~~~~~

bcbio implements configurable SNP, indel and structural variant calling. We
include whole genome and exome evaluations against reference calls from
the `Genome in a Bottle`_ consortium, enabling continuous assessment of new
alignment and variant calling algorithms. We regularly report on these
comparisons and continue to improve approaches as the community makes new
tools available. Here is some of the research that contributes to the
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
  calling on the whole genome (``coverage_interval`` is genome).

- `Validation of structural variant detection`_ using multiple calling
  methods. This implements a pipeline that works in tandem with SNP and indel
  calling to detect larger structural variations like deletions, duplications,
  inversions and copy number variants (CNVs).

- An `evaluation of joint calling`_ with GATK HaplotypeCaller, FreeBayes,
  Platypus and samtools. This validates the joint calling implementation,
  allowing scaling of large population germline experiments. It also
  demonstrates improved performance of new callers: samtools 1.0 and Platypus.

bcbio automates post-variant calling annotation to make
the outputs easier to feed directly into your biological analysis. We annotate
variant effects using `snpEff`_ or `Variant Effect Predictor`_ (VEP), and
prepare a `GEMINI database`_ that associates variants with multiple
external annotations in a SQL-based query interface.

.. _Genome in a Bottle: http://www.genomeinabottle.org/
.. _variant evaluation framework: https://bcbio.wordpress.com/2013/05/06/framework-for-evaluating-variant-detection-methods-comparison-of-aligners-and-callers/
.. _FreeBayes and BAM post-alignment processing: https://bcbio.wordpress.com/2013/10/21/updated-comparison-of-variant-detection-methods-ensemble-freebayes-and-minimal-bam-preparation-pipelines/
.. _improve variant filtering: http://bcbio.wordpress.com/2014/05/12/wgs-trio-variant-evaluation/
.. _Validation of structural variant detection: http://bcbio.wordpress.com/2014/08/12/validated-whole-genome-structural-variation-detection-using-multiple-callers/

.. _GATK UnifiedGenotyper: http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html
.. _GATK HaplotypeCaller: http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html
.. _samtools mpileup: http://samtools.sourceforge.net/mpileup.shtml
.. _GATK's Variant Quality Score Recalibrator: http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html
.. _bwa mem: http://bio-bwa.sourceforge.net/
.. _novoalign: http://www.novocraft.com
.. _snpEff: http://snpeff.sourceforge.net/
.. _GEMINI database: http://gemini.readthedocs.org/en/latest/
.. _Variant Effect Predictor: http://www.ensembl.org/info/docs/tools/vep/index.html
.. _evaluation of joint calling: http://bcbio.wordpress.com/2014/10/07/joint-calling/

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

Another good source of inspiration are the configuration files from the
:ref:`example-pipelines`, which may help identify other configuration variables
of interest.

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
  calling. Specifying a ``jointcaller`` along with the appropriate
  ``variantcaller`` in the :ref:`variant-config` configuration enables this::

    - description: Sample1
      algorithm:
        variantcaller: freebayes
        jointcaller: freebayes-joint
      metadata:
        batch: Batch1
    - description: Sample2
      algorithm:
        variantcaller: freebayes
        jointcaller: freebayes-joint
      metadata:
        batch: Batch1

Cancer variant calling
~~~~~~~~~~~~~~~~~~~~~~

We support cancer calling with tumor and optionally matched normal pairs using
multiple callers:

   - `MuTect`_ (version 1.1.5 and above)
   - `VarScan`_
   - `FreeBayes`_
   - `VarDict`_

We're actively working on evaluating and tuning these callers using `synthetic
dataset 3 from the ICGC-TCGA DREAM challenge`_ and will be reporting on those
results and suggesting recommended approaches.

.. _MuTect: http://www.broadinstitute.org/cancer/cga/mutect
.. _VarScan: http://varscan.sourceforge.net
.. _VarDict: https://github.com/AstraZeneca-NGS/VarDict
.. _FreeBayes: https://github.com/ekg/freebayes
.. _synthetic dataset 3 from the ICGC-TCGA DREAM challenge: https://www.synapse.org/#!Synapse:syn312572/wiki/62018

RNA-seq
~~~~~~~

bcbio-nextgen also implements a configurable best-practice pipeline for RNA-seq
quality control, adapter trimming, alignment and post-alignment quantitation

- Adapter trimming:
  - `AlienTrimmer`_

- Sequence alignment:
  - `tophat2`_
  - `STAR`_

- Quality control:
  - `RNA-SeQC`_
  - `FastQC`_

- Quantitation:
  - `HTSeq`_

After a run you will have in the ``upload`` directory a directory for each
sample which contains a BAM file of the aligned and unaligned reads, a
``cufflinks`` directory with the output of Cufflinks, including FPKM values,
and a ``qc`` directory with plots from FastQC and RNA-SeQC. It is useful to look
at the fastqc report an the RNA-SeQC report for each of your samples to ensure
nothing looks abnormal.

In addition to directories for each sample, in the ``upload`` directory there is
a project directory which contains a YAML file describing some summary statistics
about each sample and some provenance data. In that directory is also a
``combined.counts`` file which can be used as a starting point for performing
differential expression calling using any count-based method such as EdgeR,
DESeq2 or voom+limma, etc.

Standard
~~~~~~~~

This pipeline implements `alignment` and `qc` tools. Furthermore, it will run `qsignature`_ to detect possible duplicated samples, or miss-labeling. It uses SNPs signature to create a distance matrix that helps easily to create groups. The project yaml file will show number of total samples analyzed, number of very similar samples, and samples that could be duplicated.

.. _qsignature: http://sourceforge.net/p/adamajava/wiki/qSignature/

Configuration
=============
We will assume that you installed bcbio-nextgen with the automated installer,
and so your default `bcbio_system.yaml`_ file is configured correctly with all
of the tools pointing to the right places. If that is the case to run
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
process.  You can set a great many `parameters`_ under each section but most of
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
the template system as long as you the forward and reverse read filenames are
the same, barring a _1 and _2. For example: control_1.fastq and control_2.fastq
will be detected as paired and combined in the YAML file output by the
templating system.


.. _bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
.. _tophat2: http://tophat.cbcb.umd.edu/
.. _STAR: http://code.google.com/p/rna-star/
.. _cutadapt: http://code.google.com/p/cutadapt/
.. _RNA-SeQC: https://www.broadinstitute.org/cancer/cga/rna-seqc
.. _FastQC: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _HTSeq: http://www-huber.embl.de/users/anders/HTSeq/doc/index.html
.. _TruSeq: http://www.illumina.com/products/truseq_rna_sample_prep_kit_v2.ilmn
.. _bcbio_system.yaml: http://github.com/chapmanb/bcbio-nextgen/blob/master/config/bcbio_system.yaml
.. _configuration documentation: http://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html#upload
.. _parameters: http://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html
.. _template: http://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html#automated-sample-configuration
.. _illumina-rnaseq: http://raw.github.com/chapmanb/bcbio-nextgen/master/config/templates/illumina-rnaseq.yaml
.. _AlienTrimmer: http://www.ncbi.nlm.nih.gov/pubmed/23912058
