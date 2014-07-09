Pipelines
---------

Variant calling
~~~~~~~~~~~~~~~

bcbio-nextgen implements configurable best-practice pipelines for SNP
and small indel calling:

-  Sequence alignment:

   - `bwa mem`_
   - `novoalign`_
   - `bowtie2`_
   - `mosaik`_

-  Base Quality Recalibration
-  Realignment around indels
-  Variant calling:

   -  `GATK Unified Genotyper`_ (supports both GATK-lite in GATK 2.3
      and commercial restricted version in GATK 2.4+)
   -  `GATK Haplotype caller`_ (part of the commercially restricted GATK 2.4+)
   -  `FreeBayes`_
   -  `samtools mpileup`_
   -  `cortex\_var`_
   -  `VarScan`_

-  Paired tumor / normal variant calling:

   - `MuTect`_ (version 1.1.5 and above)
   - `VarScan`_

-  Quality filtering, using either hard filtering or
   `GATK's Variant Quality Score Recalibrator`_ (VQSR). VQSR
   requires a large number of variants. Practically this means high
   depth whole genome variant calling experiments. bcbio-nextgen
   attempts VQSR with the following :ref:`algorithm-config`

   - ``variantcaller`` is gatk or gatk-haplotype
   - ``coverage_depth`` is not low
   - ``coverage_interval`` is genome

-  Annotation of variant effects, using `snpEff`_
-  Variant exploration and prioritization, using `GEMINI`_

It follows approaches from:

- `GATK best practice`_ guidelines for variant calling
- Marth Lab's `gkno pipelines`_


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
	genome_build: GRCm38
	analysis: RNA-seq
	algorithm:
             aligner: tophat2
	     quality_format: Standard
	     trim_reads: read_through
	     adapters: [nextera, polya]
      - files: [/full/path/to/control_rep2.fastq]
	description: 'Control_rep2'
	genome_build: GRCm38
	analysis: RNA-seq
	algorithm:
             aligner: tophat2
	     quality_format: Standard
	     trim_reads: read_through
	     adapters: [nextera, polya]
      - files: [/full/path/to/case_rep1.fastq]
	description: 'Case_rep1'
	genome_build: GRCm38
	analysis: RNA-seq
	algorithm:
             aligner: tophat2
	     quality_format: Standard
	     trim_reads: read_through
	     adapters: [nextera, polya]
      - files: [/full/path/to/case_rep2.fastq]
	description: 'Case_rep2'
	genome_build: GRCm38
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


.. _GATK best practice: http://gatkforums.broadinstitute.org/discussion/1186/best-practice-variant-detection-with-the-gatk-v4-for-release-2-0
.. _GATK Unified Genotyper: http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html
.. _GATK Haplotype caller: http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html
.. _FreeBayes: https://github.com/ekg/freebayes
.. _samtools mpileup: http://samtools.sourceforge.net/mpileup.shtml
.. _cortex\_var: http://cortexassembler.sourceforge.net/index_cortex_var.html
.. _GATK's Variant Quality Score Recalibrator: http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html
.. _snpEff: http://snpeff.sourceforge.net/
.. _bwa mem: http://bio-bwa.sourceforge.net/
.. _bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
.. _mosaik: https://github.com/wanpinglee/MOSAIK
.. _novoalign: http://www.novocraft.com
.. _gkno pipelines: http://gkno.me/pipelines.html
.. _GEMINI: http://gemini.readthedocs.org/en/latest/
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
.. _VarScan: http://varscan.sourceforge.net
.. _MuTect: http://www.broadinstitute.org/cancer/cga/mutect
.. _AlienTrimmer: http://www.ncbi.nlm.nih.gov/pubmed/23912058
