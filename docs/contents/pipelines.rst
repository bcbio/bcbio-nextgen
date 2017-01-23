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
  cutoff-based soft filtering. VQSR requires a large number of variants and we use
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
external annotations in a SQL-based query interface. GEMINI databases have the
most associated external information for human samples (GRCh37/hg19 and hg38)
but are available for any organism with the database populated using the VCF
INFO column and predicted effects.

.. _Genome in a Bottle: http://www.genomeinabottle.org/
.. _variant evaluation framework: https://bcb.io/2013/05/06/framework-for-evaluating-variant-detection-methods-comparison-of-aligners-and-callers/
.. _FreeBayes and BAM post-alignment processing: https://bcb.io/2013/10/21/updated-comparison-of-variant-detection-methods-ensemble-freebayes-and-minimal-bam-preparation-pipelines/
.. _improve variant filtering: http://bcb.io/2014/05/12/wgs-trio-variant-evaluation/
.. _FreeBayes: https://github.com/ekg/freebayes
.. _GATK UnifiedGenotyper: https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_genotyper_UnifiedGenotyper.php
.. _GATK HaplotypeCaller: https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
.. _samtools mpileup: http://samtools.sourceforge.net/mpileup.shtml
.. _GATK's Variant Quality Score Recalibrator: https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php
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
  using downloads from the GATK Broad site or Appistry (see :ref:`toolplus-install`).

You may also want to enable :ref:`svs-pipeline` for detection of larger events,
which work with either caller. Another good source of inspiration are the
configuration files from the :ref:`example-pipelines`, which may help identify
other configuration variables of interest. A more complex setup with multiple
callers and resolution of ensemble calls is generally only useful with a small
population where you are especially concerned about sensitivity. Single
caller detection with FreeBayes or GATK HaplotypeCaller provide good resolution
of events.

.. _population-calling:

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
  them to the variant caller. This works for smaller batch sizes (< 100 samples)
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

.. _cancer-calling:

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

Cancer calling handles both tumor-normal paired calls and tumor-only calling.
For tumor-only samples, bcbio will try to remove likely germline variants
present in the public databases like 1000 genomes and ExAC, and not in COSMIC.
This runs as long as you have a local GEMINI installation and marks likely
germline variants with a ``LowPriority`` filter. `This post has more details
<http://bcb.io/2015/03/05/cancerval/>`_ on the approach and validation.

The standard variant outputs (``sample-caller.vcf.gz``) for tumor calling
emphasize somatic differences, those likely variants unique to the cancer. In
addition to this file, we also produce a ``sample-caller-germline.vcf.gz`` file
containing likely germline mutations. These are useful for identifying
pre-existing genomic changes that can contribute to cancer development, or in
paired cases like pre and post treatment where you may want to identify
maintained mutations after treatment.

We're actively working on improving calling to better account for the
heterogeneity and structural variability that define cancer genomes.

.. _full evaluation of cancer calling: http://bcb.io/2015/03/05/cancerval/
.. _synthetic dataset 3 from the ICGC-TCGA DREAM challenge: https://www.synapse.org/#!Synapse:syn312572/wiki/62018

.. _somatic-w-germline-variants:

Somatic with germline variants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
For tumor/normal somatic samples, bcbio can call both somatic (tumor-specific)
and germline (pre-existing) variants. The typical outputs of
:ref:`cancer-calling` are likely somatic variants acquired by the cancer, but
pre-existing germline risk variants are often also diagnostic.

bcbio enables calling both somatic and germline variants within a single
pipeline run. Since the algorithms for calling both are different, you specify
which callers to use for each step in the :ref:`variant-config` configuration::

    description: your-normal
    variantcaller:
       somatic: vardict
       germline: freebayes

bcbio does a single alignment for the normal sample, then splits at the variant
calling steps. In this example, you'd get FreeBayes germline calls labeled as
``your-normal-germline`` and VarDict somatic calls for the tumor sample linked
to this normal.

Germline calling supports multiple callers, and other configuration options like
ensemble and structural variant calling inherit from the remainder configuration. For
example, to use 3 callers for somatic and germline calling, create ensemble calls
for both and include germline and somatic events from two structural variant
callers::

    variantcaller:
       somatic: [vardict, varscan, mutect2]
       germline: [freebayes, gatk-haplotype, platypue]
    ensemble:
       numpass: 2
    svcaller: [manta, cnvkit]

Tumor-only inputs mix somatic and germline variants, making it difficult to
separate events. For tumor-only cases we suggest running standard
:ref:`cancer-calling`. bcbio will attempt to distinguish somatic and germline
mutations using the presence of variants in population databases.

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
bcbio can also be use to analyze RNA-seq data. It includes steps for quality
control, adapter trimming, alignment, variant calling, transcriptome
reconstruction and post-alignment quantitation at the level of the gene
and isoform.

We recommend using the STAR aligner for all genomes where there are no alt
alleles. For genomes such as hg38 that have alt alleles, hisat2 should be used
as it handles the alts correctly and STAR does not yet. Use Tophat2 only
if you do not have enough RAM available to run STAR (about 30 GB).

Our current recommendation is to run adapter trimming only if using the Tophat2
aligner. Adapter trimming is very slow, and aligners that soft clip the ends of
reads such as STAR and hisat2, or algorithms using pseudoalignments like
Sailfish handle contaminant sequences at the ends properly. This makes trimming
unnecessary. Tophat2 does not perform soft clipping so if using Tophat2,
trimming must still be done.

Sailfish, which is an extremely fast alignment-free method of quantitation, is
run for all experiments. Sailfish can accurately quantitate the expression of
genes, even ones which are hard to quantitate with other methods (see `this
paper <http://www.genomebiology.com/2015/16/1/177>`_ for example). It also
quantitates at the transcript level which can help gene-level analyses (see
`this paper <http://f1000research.com/articles/4-1521/v1>`_ for example).
We recommend using the Sailfish quantitation rather than the counts from
featureCounts to perform downstream quantification.

Although we do not recommend using the featureCount based counts, the alignments
are still useful because they give you many more quality metrics than the
pseudoalignments from Sailfish.

After a bcbio RNA-seq run there will be in the ``upload`` directory a directory
for each sample which contains a BAM file of the aligned and unaligned reads, a
``Sailfish`` directory with the output of Sailfish, including TPM values, and a
``qc`` directory with plots from FastQC and qualimap.

In addition to directories for each sample, in the ``upload`` directory there is
a project directory which contains a YAML file describing some summary
statistics for each sample and some provenance data about the bcbio run. In that
directory is also a ``combined.counts`` file with the featureCounts derived
counts per cell.

fast RNA-seq
~~~~~~~~~~~~
This mode of ``bcbio-nextgen`` quantitates transcript expression using `Salmon
<http://salmon.readthedocs.org/en/latest/>`_ and does nothing else. It is an
order of magnitude faster or more than running the full RNA-seq analysis. The
cost of the increased speed is that you will have much less information about
your samples at the end of the run, which can make troubleshooting trickier.
Invoke with ``analysis: fastrna-seq``.

single-cell RNA-seq
~~~~~~~~~~~~~~~~~~~
bcbio-nextgen supports universal molecular identifiers (UMI) based single-cell
RNA-seq analyses. If your single-cell prep does not use universal molecular
identifiers (UMI), you can most likely just run the standard RNA-seq pipeline
and use the results from that. The UMI are used to discard reads which
are possibly PCR duplicates and is very helpful for removing some of the
PCR duplicate noise that can dominate single-cell experiments.

Unlike the standard RNA-seq pipeline, the single-cell pipeline expects the FASTQ
input files to not be separated by cellular barcode, so each file is a mix of
cells identified by a cellular barcode (CB), and unique reads from a transcript
are identified with a UMI. bcbio-nextgen inspects each read, identifies the
cellular barcode and UMI and puts them in the read name. Then the reads are
aligned to the transcriptome with `RapMap <https://github.com/COMBINE-lab/RapMap>`_
and the number of reads aligning to each transcript is counted for each cellular
barcode. The output is a table of counts with transcripts as the rows and
columns as the cellular barcodes for each input FASTQ file.

Optionally the reads can be quantitated with ``kallisto`` to output transcript
compatibility counts rather than counts per gene
(`TCC paper <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0970-8>`_)
``kallisto`` is free for academic use, but if you are a commerical entity,
you need a `license <https://pachterlab.github.io/kallisto/download>`_ from
UC Berkeley.

To extract the UMI and cellular barcodes from the read, bcbio-nextgen
needs to know where the UMI and the cellular barcode are expected to be
in the read. Currently there is support for two schemes, the inDrop system from
the Harvard single-cell core facility and CEL-seq. If bcbio-nextgen does not
support your UMI and barcoding scheme, please open up an issue and we will
help implement support for it.

Most of the heavy lifting for this part of `bcbio-nextgen` is implemented in
the `umis <https://github.com/vals/umis>`_ repository.

smallRNA-seq
~~~~~~~~~~~~

bcbio-nextgen also implements a configurable best-practices pipeline for smallRNA-seq
quality controls, adapter trimming, miRNA/isomiR quantification and other small RNA
detection.

- Adapter trimming:
  - `cutadapt`_

- Sequence alignment:
  - `STAR`_ for genome annotation
  - bowtie, `bowtie2` and  `hisat2`_ for genome annotation as an option

- Known small RNAs quantification:
  - `seqbuster <https://github.com/lpantano/seqbuster>`_ for miRNA annotation
  - `tdrmapper`_ for tRNA fragments annotation

- Quality control:
  - `FastQC`_

- Other small RNAs quantification:
  - `seqcluster <https://github.com/lpantano/seqcluster>`_
  - `mirDeep2`_ for miRNA prediction

The pipeline generates a _RMD_ template file inside ``report`` folder
that can be rendered with knitr. An example of the report is at `here <https://github.com/lpantano/mypubs/blob/master/srnaseq/mirqc/ready_report.md>`_.
Count table (``counts_mirna.tst``) from mirbase miRNAs will be
inside ``mirbase`` or final project folder.
Input files for `isomiRs`_ package for isomiRs analysis will be
inside each sample in ``mirbase`` folder..
If mirdeep2 can run, count table (``counts_mirna_novel.tsv``)
for novel miRNAs will be inside
``mirdeep2`` or final project folder.
tdrmapper results will be inside each sample
inside ``tdrmapper`` or final project folder.

.. _tdrmapper: https://github.com/sararselitsky/tDRmapper
.. _miRDeep2: https://www.mdc-berlin.de/8551903/en/
.. _isomiRs: https://github.com/lpantano/isomiRs

ChIP-seq
~~~~~~~~
bcbio-nextgen implements the first steps of a ChIP-seq analysis up to aligning with
bowtie2. It does alignment and peak calling with MACS2.

- Adapter trimming:
  - `cutadapt`_

- Sequence alignment:
  - `bowtie2`_

- Peak calling:
  - `macs2`_

- Quality control:
  - `FastQC`_

.. _macs2: https://github.com/taoliu/MACS

Standard
~~~~~~~~

This pipeline implements ``alignment`` and ``qc`` tools. Furthermore, it will
run `qsignature`_ to detect possible duplicated samples, or miss-labeling. It
uses SNPs signature to create a distance matrix that helps easily to create
groups. The project yaml file will show the number of total samples analyzed,
the number of very similar samples, and samples that could be duplicated.

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
