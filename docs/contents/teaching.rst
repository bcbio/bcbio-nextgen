.. _teaching:

Teaching 
--------

Single cell RNA-seq analysis
~~~~~~~~~~~~~~~~~~~~~~~~
`Setting up bcbio single cell RNA-seq analysis <https://github.com/hbc/tutorials/blob/master/scRNAseq/scRNAseq_analysis_tutorial/lessons/01_bcbio_run.md>`_
tutorial outlines the steps needed to run bcbio in that use case.

Cancer tumor-normal variant calling
~~~~~~~~~~~~~~~~~~~~~~~~
This is a teaching orientated example of using bcbio from the Cold Spring Harbor
Laboratory's `Advanced Sequencing Technology and Applications course
<http://meetings.cshl.edu/courses.aspx?course=C-SEQTEC&year=15>`_. This uses
cancer tumor normal data from the `ICGC-TCGA DREAM synthetic 3 challenge
<https://www.synapse.org/#!Synapse:syn312572/wiki/58893>`_, subset to exomes on
chromosome 6 to reduce runtimes. It demonstrates:

- Running a cancer tumor/normal workflow through bcbio.
- Analysis with human genome build 38.
- SNP and indel detection, with 3 variant callers and an ensemble method.
- Structural variant calling, with 2 callers.
- Prioritization of structural variants for cancer associated genes in
  `CIViC <https://civic.genome.wustl.edu/#/home>`_.
- HLA typing.
- Validation of both small and structural variants against truth sets.

Loading pre-run analysis
========================

To save downloading the genome data and running the analysis, we have a
pre-prepared AMI with the data and analysis run. Use the `AWS Console
<https://console.aws.amazon.com/ec2>`_ to launch the pre-built AMI -- search
Community AMIs for ami-5e84fe34. Any small instance type is fine for exploring
the configuration, run directory and output files. Make sure you associate a
public IP and a security group that allows remote ssh.

Once launched, ssh into the remote machine with ``ssh -i your-keypair
ubuntu@public.ip.address`` to explore the inputs and outputs.
The default PATH contains bcbio and third party programs in ``/usr/local/bin``,
with the biological data installed in ``/usr/local/share/bcbio``. The run is in
a ``~/run/cancer-syn3-chr6``.

Input configuration file
========================

To run bcbio, you prepare a small configuration file describing your analysis.
You can `prepare it manually or use an automated configuration method <https://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html>`_.
The example has a pre-written configuration file with tumor/normal data located
in the ``config`` directory and this section walks through the settings.

You define the type of analysis (variant calling) along with the input files and
genome build::

    analysis: variant2
    files: [../input/cancer-syn3-chr6-tumor-1.fq.gz, ../input/cancer-syn3-chr6-tumor-2.fq.gz]
    genome_build: hg38

Sample description and assignment as a tumor sample, called together with a
matched normal::

    description: syn3-tumor
    metadata:
      batch: syn3
      phenotype: tumor
      sex: female

Next it defines parameters for running the analysis. First we pick our aligner
(bwa mem)::

    algorithm:
      aligner: bwa

Post-alignment, we mark duplicates but do not perform recalibration and realignment::

      mark_duplicates: true
      recalibrate: false
      realign: false

We call variants in exome regions on chromosome 6 using a BED file input, call
variants as low as 2% in the tumor sample, and use 3 variant callers along with
an ensemble method that combines results for any found in 2 out of 3::

      variant_regions: ../input/NGv3-chr6-hg38.bed
      min_allele_fraction: 2
      variantcaller: [vardict, freebayes, varscan]
      ensemble:
        numpass: 2

For structural variant calling, we use two callers and prioritize variants to
those found in the CIViC database::

      svcaller: [lumpy, manta]
      svprioritize: cancer/civic

Call HLA types with OptiType::

      hlacaller: optitype

Finally, we validate both the small variants and structural variants. These use
pre-installed validation sets that come with bcbio. We limit validation regions
to avoid low complexity regions, which cause bias in `validating indels
<http://bcb.io/2014/05/12/wgs-trio-variant-evaluation/>`_::

      exclude_regions: [lcr]
      validate: dream-syn3-crossmap/truth_small_variants.vcf.gz
      validate_regions: dream-syn3-crossmap/truth_regions.bed
      svvalidate:
        DEL: dream-syn3-crossmap/truth_DEL.bed
        DUP: dream-syn3-crossmap/truth_DUP.bed
        INV: dream-syn3-crossmap/truth_INV.bed

Output files
============

Output files are in ``~/run/cancer-syn3-chr6/final``, extracted from the full
work directory in ``~/run/cancer-syn3-chr6/work``.

The directories with sample information are in ``syn3-tumor/``. Aligned BAMs
include a ``-ready.bam`` file with all of the original reads (including split
and discordants) and separate files with only the split (``-sr.bam``) and
discordant (``-disc.bam``) reads::

    syn3-tumor-ready.bam
    syn3-tumor-ready.bam.bai
    syn3-tumor-sr.bam
    syn3-tumor-sr.bam.bai
    syn3-tumor-disc.bam
    syn3-tumor-disc.bam.bai

SNP and indel calls for 3 callers, plus combined ensemble calls::

    syn3-tumor-ensemble.vcf.gz
    syn3-tumor-ensemble.vcf.gz.tbi
    syn3-tumor-freebayes.vcf.gz
    syn3-tumor-freebayes.vcf.gz.tbi
    syn3-tumor-varscan.vcf.gz
    syn3-tumor-varscan.vcf.gz.tbi
    syn3-tumor-vardict.vcf.gz
    syn3-tumor-vardict.vcf.gz.tbi

Structural variant calls for 2 callers, plus a simplified list of structural
variants in cancer genes of interest::

    syn3-tumor-sv-prioritize.tsv
    syn3-tumor-lumpy.vcf.gz
    syn3-tumor-lumpy.vcf.gz.tbi
    syn3-tumor-manta.vcf.gz
    syn3-tumor-manta.vcf.gz.tbi

HLA typing results::

    syn3-tumor-hla-optitype.csv

Validation results from comparisons against truth set, including plots::

    syn3-tumor-sv-validate.csv
    syn3-tumor-sv-validate-DEL.png
    syn3-tumor-sv-validate-df.csv
    syn3-tumor-sv-validate-DUP.png
    syn3-tumor-sv-validate-INV.png
    syn3-tumor-validate.png

The top level directory for the project, ``2015-11-18_syn3-cshl/`` has files
relevant to the entire run. There is a consolidated quality control report::

    multiqc/multiqc_report.html

Povenance information, with log files of all commands run and program versions used::

    bcbio-nextgen.log
    bcbio-nextgen-commands.log
    programs.txt
    data_versions.csv

A top level summary of metrics for alignment, variant calling and coverage that
is useful downstream::

    project-summary.yaml

Preparing and Running
=====================
The steps to prepare an AMI from a bare machine and run the analysis. These are
pre-done on the teaching AMI to save time:

1. Use the `AWS Console <https://console.aws.amazon.com/ec2>`_ to launch
   a Ubuntu Server 14.04 (ami-d05e75b8). Start an m4.4xlarge instance with a
   100Gb SSD. Make sure you associate a public IP and can ssh in externally.

2. SSH to your instance::

     ssh -i ~/.ec2/your-key.pem ubuntu@public-ip

3. Install bcbio with hg38 data::

     sudo apt-get update
     sudo apt-get install -y build-essential zlib1g-dev wget curl python-setuptools git \
                             openjdk-7-jdk openjdk-7-jre ruby libncurses5-dev libcurl4-openssl-dev \
                             libbz2-dev unzip pigz bsdmainutils
     wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
     python bcbio_nextgen_install.py /usr/local/share/bcbio --tooldir /usr/local \
            --genomes hg38 --aligners bwa --sudo --isolate -u development

4. Install the analysis data::

     mkdir -p run
     cd run
     wget https://raw.githubusercontent.com/bcbio/bcbio-nextgen/master/config/teaching/cancer-syn3-chr6-prep.sh
     bash cancer-syn3-chr6-prep.sh

5. Run the analysis::

     cd cancer-syn3-chr6/work
     bcbio_nextgen.py ../config/cancer-syn3-chr6.yaml -n 16
