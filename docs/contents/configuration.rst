.. _docs-config:

Configuration
-------------

Two configuration files, in easy to write `YAML format`_, specify
details about your system and samples to run:

- ``bcbio_system.yaml`` High level information about the system, including
  locations of installed programs like GATK and cores and memory usage (see
  :ref:`tuning-cores`). These apply across multiple runs. The automated
  installer creates a ready to go system configuration file that can be manually
  edited to match the system. Find the file in the galaxy sub-directory within
  your installation data location (ie.
  ``/usr/local/share/bcbio-nextgen/galaxy``). By default, the pipeline uses the
  standard pre-created configuration file but multiple system configurations can
  be independently maintained and passed as the first argument to
  ``bcbio_nextgen.py`` commands.

- ``bcbio_sample.yaml`` Details about a set of samples to process,
  including input files and analysis options. You configure these for
  each set of samples to process.

Commented `system`_ and `sample`_ example files are available in the
``config`` directory. The :ref:`example-pipelines` section contains
additional examples of ready to run sample files.

.. _automated-sample-config:

Automated sample configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bcbio-nextgen provides a utility to create configuration files for
multiple sample inputs using a base template. Start with one of
the `best-practice templates`_, or define your own, then apply to
multiple samples using the template workflow command::

    bcbio_nextgen.py -w template freebayes-variant project1.csv sample1.bam sample2_1.fq sample2_2.fq

- ``freebayes-variant`` is the name of the standard ``freebayes-variant.yaml``
  input, which the script fetches from GitHub. This argument can also
  be a path to a locally customized YAML configuration. In both cases,
  the script replicates the single sample template configuration to
  all input samples.

- ``project1.csv`` is a comma separated value file containing sample
  metadata, descriptions and algorithm tweaks::

        samplename,description,batch,phenotype,sex,variant_regions
        sample1,ERR256785,batch1,normal,female,/path/to/regions.bed
        sample2,ERR256786,batch1,tumor,,/path/to/regions.bed

  The first column links the metadata to a specific input file. The
  template command tries to identify the ``samplename`` from read group
  information in a BAM file, or uses the base filename if no read group
  information is present. For BAM files, this would be the filename without the
  extension and path (``/path/to/yourfile.bam => yourfile``). For fastq
  files, the template functionality will identify pairs using standard
  conventions (``_1`` and ``_2``, including Illumina extensions like ``_R1``),
  so use the base filename without these (``/path/to/yourfile_R1.fastq => yourfile``).

    The remaining columns can contain:

   - ``description`` Changes the sample description, originally
     supplied by the file name or BAM read group, to this value. You can also
     set the ``lane``, although this is less often done as the default
     sequential numbering works here. See the documentation for
     :ref:`sample-configuration` on how these map to BAM read groups.

   - Algorithm parameters specific for this sample. If the column name matches
     an available :ref:`algorithm-config`, then this value substitutes
     into the sample ``algorithm``, replacing the defaults from the template.
     You can also change other information in the BAM read group through the
     ``algorithm`` parameters. See :ref:`alignment-config` configuration
     documentation for details on how these map to read group information.

   -  :ref:`sample-configuration` metadata key/value pairs. Any columns not
      falling into the above cases will go into the metadata section. A ``ped``
      specification will allow bcbio to read family, gender and phenotype
      information from a PED input file.

  Individual column items can contain booleans (true or false), integers, or
  lists (separated by semi-colons). These get converted into the expected time
  in the output YAML file. For instance, to specify a sample that should go into
  multiple batches::

       samplename,description,phenotype,batch
       normal.bam,two_normal,normal,Batch1;Batch2

  The name of the metadata file, minus the ``.csv`` extension, is a
  short name identifying the current project. The script creates a
  ``project1`` directory containing the sample configuration in
  ``project1/config/project1.yaml``.

- The remaining arguments are input BAM or fastq files. The script
  pairs fastq files (identified by ``_1`` and ``_2``) and extracts
  sample names from input BAMs, populating the ``files`` and
  ``description`` field in the final configuration file. Specify the
  full path to sample files on your current machine.

To make it easier to define your own project specific template, an
optional first step is to download and edit a local template. First
retrieve a standard template::

    bcbio_nextgen -w template freebayes-variant project1

This pulls the current GATK best practice variant calling template
into your project directory in
``project1/config/project1-template.yaml``. Manually edit this file to
define your options, then run the full template creation for your
samples, pointing to this custom configuration file::


    bcbio_nextgen -w template project1/config/project1-template.yaml project1.csv folder/*

If your sample folder contains additional BAM or fastq files you do not wish to
include in the sample YAML configuration, you can restrict the output to only
include samples in the metadata CSV with ``--only-metadata``. The output will
print warnings about samples not present in the metadata file, then leave these
out of the final output YAML::

    bcbio_nextgen -w template --only-metadata project1/config/project1-template.yaml project1.csv folder/*


.. _best-practice templates: https://github.com/chapmanb/bcbio-nextgen/tree/master/config/templates

.. _multi-files-sample-configuration:

Multiple files per sample
~~~~~~~~~~~~~~~~~~~~~~~~~

In case you have multiple FASTQ or BAM files for each sample you can use ``bcbio_prepare_samples.py``.
The main parameters are:

- ``--out``: the folder where the merged files will be
- ``--csv``: the CSV file that is exactly the same than described previously, but having as many duplicate lines for each samples as files to be merged::


        samplename,description,batch,phenotype,sex,variant_regions
        file1.fastq,sample1,batch1,normal,female,/path/to/regions.bed
        file2.fastq,sample1,batch1,normal,female,/path/to/regions.bed
        file1.fastq,sample2,batch1,tumor,,/path/to/regions.bed

An example of usage is::


    bcbio_prepare_samples.py --out merged --csv project1.csv

The script will create the ``sample1.fastq,sample2.fastq`` in the ``merged`` folder, and a new CSV file
in the same folder than the input CSV :``project1-merged.csv``. Later, it can be used for bcbio::


    bcbio_nextgen -w template project1/config/project1-template.yaml project1-merged.csv merged/*fastq

The new CSV file will look like::

        samplename,description,batch,phenotype,sex,variant_regions
        sample1.fastq,sample1,batch1,normal,female,/path/to/regions.bed
        sample2.fastq,sample2,batch1,tumor,,/path/to/regions.bed

It supports parallelization the same way ``bcbio_nextgen.py`` does::


    python $BCBIO_PATH/scripts/utils/bcbio_prepare_samples.py --out merged --csv project1.csv -t ipython -q queue_name -s lsf -n 1

See more examples at `parallelize pipeline`_.

.. _parallelize pipeline: https://bcbio-nextgen.readthedocs.org/en/latest/contents/parallel.html

In case of paired reads, the CSV file should contain all files::

        samplename,description,batch,phenotype,sex,variant_regions
        file1_R1.fastq,sample1,batch1,normal,female,/path/to/regions.bed
        file2_R1.fastq,sample1,batch1,normal,female,/path/to/regions.bed
        file1_R2.fastq,sample1,batch1,normal,femela,/path/to/regions.bed
        file2_R2.fastq,sample1,batch1,normal,female,/path/to/regions.bed

The script will try to guess the paired files the same way than ``bcbio_nextgen.py -w template`` does. It would detect paired files if the difference among two files is only
``_R1/_R2`` or ``-1/-2`` or ``_1/_2`` or ``.1/.2``

The output CSV will look like and is compatible with bcbio::

        samplename,description,batch,phenotype,sex,variant_regions
        sample1,sample1,batch1,normal,female,/path/to/regions.bed


.. _sample-configuration:

Sample information
~~~~~~~~~~~~~~~~~~

The sample configuration file defines ``details`` of each sample to process::

    details:
      - analysis: variant2
        lane: 1
        description: Example1
        files: [in_pair_1.fq, in_pair_2.fq]
        genome_build: hg19
        algorithm:
          platform: illumina
        metadata:
          batch: Batch1
          sex: female
          platform_unit: flowcell-barcode.lane
          library: library_type


- ``analysis`` Analysis method to use [variant2, RNA-seq, smallRNA-seq]

- ``lane`` A unique number within the project. Corresponds to the
  ``ID`` parameter in the BAM read group.

- ``description`` Unique name for this sample, corresponding to the
  ``SM`` parameter in the BAM read group. Required.

- ``files`` A list of files to process. This currently supports either a single
  end or two paired end fastq files, or a single BAM file. It does not yet
  handle merging BAM files or more complicated inputs.

- ``genome_build`` Genome build to align to, which references a genome
  keyword in Galaxy to find location build files.

- ``algorithm`` Parameters to configure algorithm inputs. Options
  described in more detail below:

  - ``platform`` Sequencing platform used. Corresponds to the ``PL``
    parameter in BAM read groups. Optional, defaults to ``illumina``.

- ``metadata`` Additional descriptive metadata about the sample:

   - ``batch`` defines a group that the sample falls in. We perform
     multi-sample variant calling on all samples with the same batch
     name. This can also be a list, allowing specification of a single normal
     sample to pair with multiple tumor samples in paired cancer variant
     calling (``batch: [MatchWithTumor1, MatchWithTumor2]``).

   - ``sex`` specifies the sample gender used to correctly prepare X/Y
     chromosomes.

   -  ``phenotype`` stratifies cancer samples into ``tumor`` and ``normal`` or
      case/controls into ``affected`` and ``unaffected``.

   - ``ped`` provides a `PED phenotype file
     <http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped>`_
     containing sample phenotype and family information. Template creation uses
     this to extract ``sex`` and ``phenotype`` information. GEMINI database
     creation uses the PED file.

   - ``platform_unit`` -- Unique identifier for sample. Optional, defaults to
     ``lane`` if not specified.

   - ``library`` -- Name of library preparation used. Optional, empty if not
     present.

   - ``validate_batch`` -- Specify a batch name to group samples together for
     preparing validation plots. This is useful if you want to process samples
     in specific batches, but include multiple batches into the same
     validation plot.

.. _upload-configuration:

Upload
~~~~~~

The ``upload`` section of the sample configuration file describes where to put
the final output files of the pipeline. At its simplest, you can configure
bcbio-nextgen to upload results to a local directory, for example a folder
shared amongst collaborators or a Dropbox account. You can also configure
it to upload results automatically to a Galaxy instance or to
`Amazon S3`_. Here is the simplest configuration, uploading to a local
directory::

     upload:
       dir: /local/filesystem/directory

General parameters, always required:

- ``method`` Upload method to employ. Defaults to local filesystem.
  [filesystem, galaxy, s3]
- ``dir`` Local filesystem directory to copy to.

Galaxy parameters:

- ``galaxy_url`` URL of the Galaxy instance to upload to. Upload
  assumes you are able to access a shared directory also present on
  the Galaxy machine.
- ``galaxy_api_key`` User API key to access Galaxy: see the
  `Galaxy API`_ documentation.
- ``galaxy_library`` Name of the Galaxy Data Library to upload to. You
  can specify this globally for a project in ``upload`` or for
  individual samples in the sample details section.
- ``galaxy_role`` Specific Galaxy access roles to assign to the
  uploaded datasets. This is optional and will default to the access
  of the parent data library if not supplied. You can specify this
  globally for a project in ``upload`` or for individual samples in
  the sample details section. The `Galaxy Admin`_ documentation
  has more details about roles.

Here is an example configuration for uploading to a Galaxy instance. This
assumes you have a shared mounted filesystem that your Galaxy instance can
also access::

      upload:
        method: galaxy
        dir: /path/to/shared/galaxy/filesystem/folder
        galaxy_url: http://url-to-galaxy-instance
        galaxy_api_key: YOURAPIKEY
        galaxy_library: data_library_to_upload_to

Your Galaxy universe_wsgi.ini configuration needs to have
``allow_library_path_paste = True`` set to enable uploads.

S3 parameters:

- ``bucket`` AWS bucket to direct output.
- ``folder`` A folder path within the AWS bucket to prefix the output.
- ``region`` AWS region name to use. Defaults to us-east-1
- ``reduced_redundancy`` Flag to determine if we should store S3 data
  with reduced redundancy: cheaper but less reliable [false, true]

For S3 access credentials, set the standard environmental variables,
``AWS_ACCESS_KEY_ID``, ``AWS_SECRET_ACCESS_KEY``, and ``AWS_DEFAULT_REGION``
or use `IAM access roles <http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/iam-roles-for-amazon-ec2.html>`_
with an instance profile on EC2 to give your instances permission to create
temporary S3 access.

Globals
~~~~~~~
You can define files used multiple times in the ``algorithm`` section of your
configuration in a top level ``globals`` dictionary. This saves copying and
pasting across the configuration and makes it easier to manually adjust the
configuration if inputs change::

  globals:
    my_custom_locations: /path/to/file.bed
  details:
    - description: sample1
      algorithm:
        variant_regions: my_custom_locations
    - description: sample2
      algorithm:
        variant_regions: my_custom_locations

.. _algorithm-config:

Algorithm parameters
~~~~~~~~~~~~~~~~~~~~

The YAML configuration file provides a number of hooks to customize
analysis in the sample configuration file. Place these under the
``algorithm`` keyword.

.. _alignment-config:

Alignment
=========

- ``platform`` Sequencing platform used. Corresponds to the ``PL``
  parameter in BAM read groups. Default 'Illumina'.
-  ``aligner`` Aligner to use: [bwa, bowtie, bowtie2, novoalign, snap, star,
   false] To use pre-aligned BAM files as inputs to the pipeline, set to
   ``false``. Using pre-aligned inputs requires proper assignment of BAM read
   groups and sorting. The ``bam_clean`` argument can often resolve issues with
   problematic input BAMs.
-  ``bam_clean`` Clean an input BAM when skipping alignment step. This
   handles adding read groups, sorting to a reference genome and
   filtering problem records that cause problems with GATK. Set to
   ``picard`` to do Picard/GATK based cleaning. To fix misencoded input BAMs
   with non-standard scores, set ``quality_format`` to ``illumina``.
-  ``bam_sort`` Allow sorting of input BAMs when skipping alignment
   step (``aligner`` set to false). Options are coordinate or
   queryname. For additional processing through standard pipelines
   requires coordinate sorted inputs. The default is to not do
   additional sorting and assume pre-sorted BAMs.
- ``disambiguate`` For mixed or explant samples, provide a list of
  ``genome_build``  identifiers to check and remove from alignment. Currently
  supports cleaning a single organism. For example, with ``genome_build: hg19``
  and ``disambiguate: [mm10]``, it will align to hg19 and mm10, run
  disambiguation and continue with reads confidently aligned to hg19. Affects
  fusion detection when ``star`` is chosen as the aligner. Aligner must be
  set to a non false value for this to run.
- ``trim_reads`` Can be set to trim low quality ends or to also trim off,
  in conjunction with the ``adapters`` field a set of adapter sequences or
  poly-A tails that could appear on the ends of reads. Only used in RNA-seq
  pipelines, not variant calling. [False, read_through]. Default to False,
  recommended to leave as False unless running Tophat2.
- ``min_read_length`` Minimum read length to maintain when
  ``read_through`` trimming set in ``trim_reads``. Defaults to 20.
-  ``adapters`` If trimming adapter read through, trim a set of stock
   adapter sequences. Allows specification of multiple items in a list,
   for example [truseq, polya] will trim both TruSeq adapter sequences
   and polyA tails. Valid items are [truseq, illumina, nextera, polya]
-  ``custom_trim`` A list of sequences to trim from the end of reads,
   for example: [AAAATTTT, GGGGCCCC]
- ``align_split_size``: Increase parallelization of alignment. As of 0.9.8,
  bcbio will try to determine a useful parameter and you don't need to set this.
  If you manually set it, bcbio will respect for you specification. Set to false
  to avoid splitting entirely. If setting, this defines the number of records to
  feed into each independent parallel step (for example, 5000000 = 5 million
  reads per chunk). It converts the original inputs into bgzip grabix indexed
  FASTQ files, and then retrieves chunks for parallel alignment. Following
  alignment, it combines all chunks back into the final merged alignment file.
  This allows parallelization at the cost of additional work of preparing inputs
  and combining split outputs. The tradeoff makes sense when you have large
  files and lots of distributed compute. When you have fewer large multicore
  machines this parameter may not help speed up processing.
-  ``quality_format`` Quality format of fastq or BAM inputs [standard, illumina]
-  ``strandedness`` For RNA-seq libraries, if your library is strand
   specific, set the appropriate flag form [unstranded, firststrand, secondstrand].
   Defaults to unstranded. For dUTP marked libraries, firststrand is correct; for
   Scriptseq prepared libraries, secondstrand is correct.

Coverage information
====================
- ``coverage_interval`` Regions covered by sequencing. bcbio calculates this
  automatically from alignment coverage information, so you only need to
  specify it in the input configuration if you have specific needs or bcbio
  does not determine coverage correctly. ``genome`` specifies full genome
  sequencing, ``regional`` identifies partial-genome pull down sequencing like
  exome analyses, and ``amplicon`` is partial-genome sequencing from
  PCR amplicon sequencing. This influences GATK options for filtering: we use
  Variant Quality Score Recalibration when set to ``genome``, otherwise we
  apply hard filters. Also affects copy number calling with CNVkit, structural
  variant calling and deep panel calling in cancer samples, where we tune
  regional/amplicon analyses to maximize sensitivity.
  [genome, regional, amplicon]
-  ``coverage_depth_min`` Minimum depth of coverage. Regions with less reads
   will not get called. Defaults to 4. Setting lower than 4 will trigger
   low-depth calling options for GATK.
- ``coverage`` A BED file of regions to check for coverage. Coverage
  and completeness are calculated over these regions and a Rmarkdown
  report is generated in the `report` directory.


Experimental information
========================

-  ``ploidy`` Ploidy of called reads. Defaults to 2 (diploid).

.. _variant-config:

Variant calling
===============

-  ``variantcaller`` Variant calling algorithm. Can be a list of
   multiple options or false to skip [gatk, freebayes, gatk-haplotype, platypus,
   mutect, mutect2, scalpel, vardict, varscan, samtools, false]

    - Paired (typically somatic, tumor-normal) variant calling is currently
      supported by vardict, freebayes, mutect2, mutect (see disclaimer below),
      scalpel (indels only) and varscan. See ``phenotype`` below for how to pair tumor
      and normal samples.
    - Selecting mutect (SNP caller) can also be combined by indels from scalpel or sid and
      combine the output. Mutect operates in both tumor-normal and tumor-only modes.
      In tumor-only mode the indels from scalpel will reflect all indels in the sample,
      as there is currently no way of separating the germline from somatic indels in
      tumor-only mode.
- ``indelcaller`` For the MuTect SNP only variant caller it is possible to add
   calls from an indelcaller such as scalpel, pindel and somatic indel detector
   (for Appistry MuTect users only). Currently an experimental option that adds
   these indel calls to MuTect's SNP-only output. Only one caller supported.
   Omit to ignore. [scalpel, pindel, sid, false]
-  ``jointcaller`` Joint calling algorithm, combining variants called with the
   specified ``variantcaller``. Can be a list of multiple options but needs to
   match with appropriate ``variantcaller``

     - ``gatk-haplotype-joint`` `GATK incremental joint discovery
       <http://www.broadinstitute.org/gatk/guide/article?id=3893>`_ with
       HaplotypeCaller. Takes individual gVCFs called by ``gatk-haploype`` and
       perform combined genotyping.
     - ``freebayes-joint`` Combine freebayes calls using
       `bcbio.variation.recall`_ with recalling at
       all positions found in each individual sample. Requires ``freebayes``
       variant calling.
     - ``platypus-joint`` Combine platypus calls using bcbio.variation.recall
       with squaring off at all positions found in each individual
       sample. Requires ``platypus`` variant calling.
     - ``samtools-joint`` Combine platypus calls using bcbio.variation.recall
       with squaring off at all positions found in each individual
       sample. Requires ``samtools`` variant calling.
-  ``variant_regions`` BED file of regions to call variants in.
-  ``mark_duplicates`` Identify and remove variants [true, false]
   If true, will perform streaming duplicate marking with `samblaster`_ for
   paired reads and `biobambam's bammarkduplicates`_ for single end reads.
-  ``recalibrate`` Perform base quality score recalibration on the
   aligned BAM file. Defaults to gatk. [false, gatk]
-  ``realign`` Perform realignment around indels on the aligned BAM
   file. Defaults to no realignment since realigning callers like FreeBayes and
   GATK HaplotypeCaller handle this as part of the calling process. [false, gatk]
- ``effects`` Method used to calculate expected variant effects. Defaults to
  `snpEff`_ and `Ensembl variant effect predictor (VEP)`_ is also available
  with support for `dbNSFP`_ annotation, when downloaded using
  :ref:`datatarget-install`. [snpeff, vep, false]
-  ``phasing`` Do post-call haplotype phasing of variants. Defaults to
   no phasing [false, gatk]
-  ``remove_lcr`` Remove variants in low complexity regions (LCRs)
   for human variant calling. `Heng Li's variant artifacts paper`_ provides
   these regions, which cover ~2% of the genome but contribute to a large
   fraction of problematic calls due to the difficulty of resolving variants
   in repetitive regions. Removal can help facilitate comparisons between
   methods and reduce false positives if you don't need calls in LCRs for your
   biological analysis. [false, true]
- ``joint_group_size`` Specify the maximum number of gVCF samples to feed into
  joint calling. Currently applies to GATK HaplotypeCaller joint calling and
  defaults to the GATK recommendation of 200. Larger numbers of samples will
  first get combined prior to genotyping.
- ``clinical_reporting`` Tune output for clinical reporting.
  Modifies snpEff parameters to use HGVS notational on canonical
  transcripts [false, true].
- ``background`` Provide a VCF file with variants to use as a background
  reference during variant calling. For tumor/normal paired calling use this to
  supply a panel of normal individuals.

.. _snpEff: http://snpeff.sourceforge.net/
.. _Ensembl variant effect predictor (VEP): http://www.ensembl.org/info/docs/tools/vep/index.html
.. _dbNSFP: https://sites.google.com/site/jpopgen/dbNSFP
.. _samblaster: https://github.com/GregoryFaust/samblaster
.. _biobambam's bammarkduplicates: https://github.com/gt1/biobambam
.. _Heng Li's variant artifacts paper: http://arxiv.org/abs/1404.0929

Structural variant calling
==========================

- ``svcaller`` -- List of structural variant callers to use. [lumpy, manta,
  cnvkit]. LUMPY, Manta and DELLY require paired end reads.
- ``svprioritize`` --  Produce a tab separated summary file of structural
  variants in regions of interest. This complements the full VCF files of
  structural variant calls to highlight changes in known genes. This can be
  either the path to a BED file (with ``chrom start end gene_name``) or the name
  of one of the pre-installed prioritization files:

     - ``cancer/civic`` (hg19, GRCh37, hg38) -- Known cancer associated genes from
       `CIViC <https://civic.genome.wustl.edu>`_.
     - ``cancer/az300`` (hg19, GRCh37, hg38) -- 300 cancer associated genes
       contributed by `AstraZeneca oncology <https://www.astrazeneca.com/our-focus-areas/oncology.html>`_.
- ``sv_regions`` -- A specification of regions to target during structural
  variant calling. By default, bcbio uses regions specified in
  ``variant_regions`` but this allows custom specification for structural
  variant calling. This can be a pointer to a bed file or special inputs:
  ``exons`` for only exon regions, ``transcripts`` for transcript regions (the
  min start and max end of exons) or ``transcriptsXXXX`` for transcripts plus a
  window of XXXX size around it. The size can be an integer (``transcripts1000``)
  or exponential (``transcripts1e5``). This applies to CNVkit and heterogeneity
  analysis.
- ``fusion_mode`` Enable fusion detection in RNA-seq when using STAR (recommended)
  or Tophat (not recommended) as the aligner. OncoFuse is used to summarise the fusions
  but currently only supports ``hg19`` and ``GRCh37``. For explant samples
  ``disambiguate`` enables disambiguation of ``STAR`` output [false, true].

HLA typing
==========
- ``hlacaller`` -- Perform identification of highly polymorphic HLAs with human
  build 38 (hg38). The recommended options is ``optitype``, using the `OptiType
  <https://github.com/FRED-2/OptiType>`_ caller. Also supports using the `bwa
  HLA typing implementation
  <https://github.com/lh3/bwa/blob/master/README-alt.md#hla-typing>`_ with ``bwakit``

Validation
===========

bcbio pre-installs standard truth sets for performing validation,
and also allows use of custom local files for assessing reliability of your
runs:

-  ``validate`` A VCF file of expected variant calls to perform
   validation and grading of small variants (SNPs and indels) from the pipeline.
   This provides a mechanism to ensure consistency of calls against
   a known set of variants, supporting comparisons to genotyping
   array data or reference materials.
- ``validate_regions`` A BED file of regions to evaluate small variant calls in. This
  defines specific regions covered by the ``validate`` VCF  file.
- ``svvalidate`` -- Dictionary of call types and pointer to BED file of known
  regions. For example: ``DEL: known_deletions.bed`` does deletion based
  validation of outputs against the BED file.

Each option can be either the path to a local file, or a partial path to a file
in the pre-installed truth sets. For instance, to validate an NA12878 run
against the `Genome in a Bottle <https://github.com/genome-in-a-bottle>`_ truth set::

    validate: giab-NA12878/truth_small_variants.vcf.gz
    validate_regions: giab-NA12878/truth_regions.bed
    svvalidate:
      DEL: giab-NA12878/truth_DEL.bed

follow the same naming schemes for small variants, regions and
different structural variant types. bcbio has the following validation materials
for germline validations:

- ``giab-NA12878`` --  `Genome in a Bottle
  <https://github.com/genome-in-a-bottle>`_ for NA12878. Truth sets: small_variants,
  regions, DEL; Builds: GRCh37, hg19.
- ``giab-NA12878-crossmap`` --  `Genome in a Bottle
  <https://github.com/genome-in-a-bottle>`_ for NA12878 converted to hg38 with CrossMap. Truth sets: small_variants,
  regions, DEL; Builds: hg38.
- ``giab-NA12878-remap`` --  `Genome in a Bottle
  <https://github.com/genome-in-a-bottle>`_ for NA12878 converted to hg38 with Remap. Truth sets: small_variants,
  regions, DEL; Builds: hg38.
- ``platinum-genome-NA12878`` -- `Illumina Platinum Genome
  <http://www.illumina.com/platinumgenomes/>`_ for NA12878. Truth sets:
  small_variants, regions; Builds: hg19, hg38.

and for cancer validations:

- ``dream-syn3`` -- Synthetic dataset 3 from the `ICGC-TCGA DREAM mutation
  calling challenge <https://www.synapse.org/#!Synapse:syn312572/wiki/62018>`_.
  Truth sets: small_variants, regions, DEL, DUP, INV, INS. Builds: GRCh37.
- ``dream-syn4`` -- Synthetic dataset 4 from the `ICGC-TCGA DREAM mutation
  calling challenge <https://www.synapse.org/#!Synapse:syn312572/wiki/62018>`_.
  Truth sets: small_variants, regions, DEL, DUP, INV. Builds: GRCh37.
- ``dream-syn3-crossmap`` -- Synthetic dataset 3 from the `ICGC-TCGA DREAM mutation
  calling challenge <https://www.synapse.org/#!Synapse:syn312572/wiki/62018>`_
  converted to human build 38 coordinates with CrossMap.
  Truth sets: small_variants, regions, DEL, DUP, INV, INS. Builds: hg38.
- ``dream-syn4-crossmap`` -- Synthetic dataset 4 from the `ICGC-TCGA DREAM mutation
  calling challenge <https://www.synapse.org/#!Synapse:syn312572/wiki/62018>`_
  converted to human build 38 coordinates with CrossMap.
  Truth sets: small_variants, regions, DEL, DUP, INV. Builds: hg38.

For more information on the hg38 truth set preparation see the work on `validation on build
38 and converstion of human build 37 truth sets to build 38
<http://bcb.io/2015/09/17/hg38-validation/>`_. The `installation recipes
<https://github.com/chapmanb/cloudbiolinux/tree/master/ggd-recipes>`_ contain
provenance details about the origins of the installed files.

.. _config-cancer:

Cancer variant calling
======================

- ``min_allele_fraction`` Minimum allele fraction to detect variants in
  heterogeneous tumor samples, set as the float or integer percentage to
  resolve (i.e. 10 = alleles in 10% of the sample). Defaults to 10. Specify this
  in the tumor sample of a tumor/normal pair.

RNA sequencing
======================

- ``transcript_assembler`` If set, will assemble novel genes and transcripts and
  merge the results into the known annotation. Can have multiple values set in a
  list. Supports ['cufflinks', 'sailfish'].
- ``transcriptome_align`` If set to True, will also align reads to just the
  transcriptome, for use with EBSeq and others.
- ``expression_caller`` A list of optional expression callers to turn on.
  Supports ['cufflinks', 'express', 'stringtie']. Sailish and count based
  expression estimation are run by default.
-  ``variantcaller`` Variant calling algorithm to call variants on RNA-seq data.
  Supports [gatk] or [vardict].

Single-cell RNA sequencing
==========================

- ``umi_type`` The UMI/cellular barcode scheme used for your data. Supports
  [harvard-indrop, harvard-indrop-v2, cel-seq].
- ``minimum_barcode_depth`` Cellular barcodes with less than this many reads
  assigned to them are discarded (default 100,000).
- ``cellular_barcodes`` An optional list of one or two files which have the
  valid cellular barcodes. Provide one file if there is only one barcode and
  two files if the barcodes are split. If no file is provided, all cellular
  barcodes passing the ``minimum_barcode_depth`` filter are kept.

smallRNA sequencing
===================

- ``adapter`` The 3' end adapter that needs to be remove.
- ``species`` 3 letters code to indicate the species in mirbase classification (i.e. hsa for human).
- ``aligner`` Currently STAR is the only one tested although bowtie can be used as well.

ChIP sequencing
===============

- ``peakcaller`` bcbio only accepts ``macs2``
- ``aligner`` Currently ``bowtie2`` is the only one tested
The ``phenotype`` and ``batch`` tags need to be set under ``metadata`` in the config YAML file. The ``phenotype`` tag will specify the chip (``phenotype: chip``) and input samples (``phenotype: input``). The ``batch`` tag will specify the input-chip pairs of samples for example, ``batch: pair1``. Same input can be used for different chip samples giving a list of distinct values: ``batch: [sample1, sample2]``.

You can pass different parameters for ``macs2`` adding to :ref:`config-resources`::


        resources:
          macs2:
            options: ["--broad"]

Quality control
===============

- ``mixup_check`` Detect potential sample mixups. Currently supports
  `qSignature <https://sourceforge.net/p/adamajava/wiki/qSignature/>`_.
  ``qsignature_full`` runs a larger analysis while ``qsignature`` runs a smaller
  subset on chromosome 22.  [False, qsignature, qsignature_full]
- ``kraken`` Turn on kraken algorithm to detect possible contamination. You can add `kraken: minikraken` and it will use a minimal database to detect possible `contaminants`_. As well, you can point to a `custom database`_ directory and kraken will use it. You will find the results in the `qc` directory. This tool only run during `rnaseq` pipeline.

.. _contaminants: https://ccb.jhu.edu/software/kraken/
.. _custom database: https://github.com/DerrickWood/kraken

Post-processing
===============

- ``archive`` Specify targets for long term archival. ``cram`` does 8-bin
  compression of BAM files into `CRAM format`_.
  Default: [] -- no archiving.
- ``tools_off`` Specify third party tools to skip as part of analysis
  pipeline. Enables turning off specific components of pipelines if not
  needed. ``gemini`` provides a `GEMINI database`_ of variants for downstream
  query during variant calling pipelines. ``vardict_somatic_filter`` disables
  running a post calling filter for VarDict to remove variants found in normal
  samples. Without ``vardict_somatic_filter`` in paired analyses no soft
  filtering of germline variants is performed but all high quality variants pass.
  ``bwa-mem`` forces use of original ``bwa aln`` alignment. Without this,
  we use bwa mem with 70bp or longer reads. ``fastqc`` turns off quality
  control FastQC usage. ``pbgzip`` turns off use of parallel bgzip
  during preparation of alignment inputs.
  ``vqsr`` turns off variant quality score recalibration for all samples.
  Default: [] -- all tools on.
- ``tools_on`` Specify functionality to enable that is off by default.
  ``svplots`` adds additional coverage and summary plots for CNVkit and detected
  ensemble variants. ``qualimap`` runs `Qualimap
  <http://qualimap.bioinfo.cipf.es/>`_ (qualimap uses downsampled files and
  numbers here are an estimation of 1e7 reads.). ``qualimap_full`` uses the full
  bam files but it may be slow. ``bwa-mem`` forces use of bwa mem even for
  samples with less than 70bp reads.  ``bnd-genotype`` enables genotyping
  of breakends in Lumpy calls, which improves accuracy but can be slow.

.. _CRAM format: http://www.ebi.ac.uk/ena/about/cram_toolkit
.. _GEMINI database: https://github.com/arq5x/gemini

parallelization
===============

- ``nomap_split_size`` Unmapped base pair regions required to split
  analysis into blocks. Creates islands of mapped reads surrounded by
  unmapped (or N) regions, allowing each mapped region to run in
  parallel. (default: 250)

- ``nomap_split_targets`` Number of target intervals to attempt to
  split processing into. This picks unmapped regions evenly spaced
  across the genome to process concurrently. Limiting targets prevents
  a large number of small targets. (default: 200)

Ensemble variant calling
========================

In addition to single method variant calling, we support calling with
multiple calling methods and consolidating into a final Ensemble
callset.

The recommended method to do this uses a simple majority rule ensemble
classifier that builds a final callset based on the intersection of calls. It
selects variants represented in at least a specified number of callers::

    variantcaller: [mutect2, varscan, freebayes, vardict]
    ensemble:
      numpass: 2
      use_filtered: false

This example selects variants present in 2 out of the 4 callers and does not use
filtered calls (the default behavior).
`bcbio.variation.recall`_ implements this approach, which handles speed and file
sorting limitations in the `bcbio.variation`_ approach.

This older approach uses the `bcbio.variation`_
toolkit to perform the consolidation. An example configuration in the
``algorithm`` section is::

    variantcaller: [gatk, freebayes, samtools, gatk-haplotype, varscan]
    ensemble:
      format-filters: [DP < 4]
      classifier-params:
        type: svm
      classifiers:
        balance: [AD, FS, Entropy]
        calling: [ReadPosEndDist, PL, PLratio, Entropy, NBQ]
      trusted-pct: 0.65

The ``ensemble`` set of parameters configure how to combine calls from
the multiple methods:

- ``format-filters`` A set of filters to apply to variants before
  combining. The example removes all calls with a depth of less than
  4.
- ``classifier-params`` Parameters to configure the machine learning
  approaches used to consolidate calls. The example defines an SVM
  classifier.
- ``classifiers`` Groups of classifiers to use for training and
  evaluating during machine learning. The example defines two set of
  criteria for distinguishing reads with allele balance issues and
  those with low calling support.
- ``trusted-pct`` Define threshold of variants to include in final
  callset. In the example, variants called by more than 65% of the
  approaches (4 or more callers) pass without being requiring SVM
  filtering.

.. _config-resources:

Resources
~~~~~~~~~

The ``resources`` section allows customization of locations of programs
and memory and compute resources to devote to them::

    resources:
      bwa:
        cores: 12
        cmd: /an/alternative/path/to/bwa
      samtools:
        cores: 16
        memory: 2G
      gatk:
        jvm_opts: ["-Xms2g", "-Xmx4g"]
        dir: /usr/share/java/gatk

- ``cmd`` Location of an executable. By default, we assume executables
  are on the path.
- ``dir`` For software not distributed as a single executable, like
  files of Java jars, the location of the base directory.
- ``cores`` Cores to use for multi-proccessor enabled software. This is how
  many cores will be allocated per job. For example if you are running
  10 samples and passed -n 40 to bcbio-nextgen and the step you are running
  has cores: 8 set, a maximum of five samples will run in parallel, each using
  8 cores.
- ``jvm_opts`` Specific memory usage options for Java software. For
  memory usage on programs like GATK, specify the maximum usage per
  core. On multicore machines, that's machine-memory divided by cores.
  This avoids memory errors when running multiple jobs simultaneously,
  while the framework will adjust memory up when running multicore
  jobs.
- ``memory`` Specify the memory per core used by a process. For programs
  where memory control is available, like ``samtools sort``,
  this limits memory usage. For other programs this is an estimate of
  usage, used by :ref:`memory-management` to avoid over-scheduling
  memory. Always specify this as the memory usage for a single core,
  and the pipeline handles scaling this when a process uses multiple
  cores.
- ``keyfile`` Specify the location of a program specific key file, obtained from
  the third party software tool. Include the path to a GATK supplied key file
  to disable the `GATK phone home`_ feature.

For GATK you can individually control memory for variant calling (which uses the
``gatk`` memory target) and for framework usage like merging and variant file
preparation (which can optionally use the the ``gatk-framework`` target). If
you only set ``gatk``, that specification gets used for framework calls as well.

Temporary directory
===================

You also use the resource section to specify system specific parameters like
global temporary directories::

    resources:
      tmp:
        dir: /scratch

This is useful on cluster systems with large attached local storage, where you
can avoid some shared filesystem IO by writing temporary files to the local
disk. When setting this keep in mind that the global temporary disk must have
enough space to handle intermediates. The space differs between steps but
generally you'd need to have 2 times the largest input file per sample and
account for samples running simultaneously on multiple core machines.

To handle clusters that specify local scratch space with an environmental
variable, bcbio will resolve environmental variables like::


    resources:
      tmp:
        dir: $YOUR_SCRATCH_LOCATION

.. _sample-resources:

Sample or run specific resources
================================

To override any of the global resource settings in a sample specific manner, you
write a resource section within your sample YAML configuration. For example, to
create a sample specific temporary directory and pass a command line option to
novoalign, write a sample resource specification like::

    - description: Example
      analysis: variant2
      resources:
        novoalign:
          options: [-o, FullNW]
        tmp:
          dir: tmp/sampletmpdir

To adjust resources for an entire run, you can add this resources specification
at the top level of your sample YAML::

     details:
       - description: Example
     resources:
       default:
         cores: 16

.. _bcbio.variation: https://github.com/chapmanb/bcbio.variation
.. _bcbio.variation.recall: https://github.com/chapmanb/bcbio.variation.recall
.. _CloudBioLinux: https://github.com/chapmanb/cloudbiolinux
.. _YAML format: https://en.wikipedia.org/wiki/YAML#Examples
.. _GATK: http://www.broadinstitute.org/gatk/
.. _system: https://github.com/chapmanb/bcbio-nextgen/blob/master/config/bcbio_system.yaml
.. _sample: https://github.com/chapmanb/bcbio-nextgen/blob/master/config/bcbio_sample.yaml
.. _Galaxy API: http://wiki.galaxyproject.org/Learn/API
.. _Amazon S3: http://aws.amazon.com/s3/
.. _Galaxy Admin: http://wiki.galaxyproject.org/Admin/DataLibraries/LibrarySecurity
.. _GATK phone home: http://gatkforums.broadinstitute.org/discussion/1250/what-is-phone-home-and-how-does-it-affect-me

Genome configuration files
~~~~~~~~~~~~~~~~~~~~~~~~~~
Each genome build has an associated ``buildname-resources.yaml``
configuration file which contains organism specific naming and
resource files. bcbio-nextgen expects a resource file present next to
the genome FASTA file. `Example genome configuration files`_ are available, and
automatically installed for natively supported genomes. Create these
by hand to support additional organisms or builds.

The major sections of the file are:

- ``aliases`` -- Names for third-party programs used as part of the
  analysis, since naming expectations can differ between software
  programs.

- ``variation`` -- Supporting data files for variant analysis. For human
  analyses, the dbSNP and training files are from the `GATK resource bundle`_.
  These are inputs into the training models for
  recalibration. The automated `CloudBioLinux`_ data scripts will
  download and install these in the variation subdirectory relative to
  the genome files.

- ``rnaseq`` -- Supporting data files for RNA-seq analysis. The
  automated installer and updater handles retrieval and installation
  of these resources for supported genome builds.

- ``srnaseq`` -- Supporting data files for smallRNA-seq analysis. Same as in
  rnaseq, the automated installer and updater handle this for supported genome
  builds.

By default, we place the ``buildname-resources.yaml`` files next to
the genome FASTA files in the reference directory. For custom setups,
you specify an alternative directory in the ref:`config-resources`
section of your ``bcbio_system.yaml`` file::

  resources:
    genome:
      dir: /path/to/resources/files

.. _Example genome configuration files: https://github.com/chapmanb/bcbio-nextgen/tree/master/config/genomes
.. _GATK resource bundle: http://www.broadinstitute.org/gatk/guide/article.php?id=1213

Reference genome files
~~~~~~~~~~~~~~~~~~~~~~

The pipeline requires access to reference genomes, including the raw
FASTA sequence and pre-built indexes for aligners. The automated installer
will install reference files and indexes for commonly used genomes (see the
:ref:`upgrade-install` documentation for command line options). For human,
GRCh37 and hg19, we use the 1000 genome references provided in the
`GATK resource bundle`_.

You can use pre-existing data and reference indexes by pointing bcbio-nextgen at
these resources. We use the `Galaxy .loc files`_ approach to describing the
location of the sequence and index data, as described in
:ref:`data-requirements`. This does not require a Galaxy installation since the
installer sets up a minimal set of ``.loc`` files. It finds these by identifying
the root ``galaxy`` directory, in which it expects a ``tool-data`` sub-directory
with the ``.loc`` files. It can do this in two ways:

- Using the directory of your ``bcbio-system.yaml``. This is the
  default mechanism setup by the automated installer and requires no additional
  work.

- From the path specified by the ``galaxy_config`` option in your
  ``bcbio-system.yaml``. If you'd like to move your system YAML file,
  add the full path to your ``galaxy`` directory here. This is useful if you
  have a pre-existing Galaxy installation with reference data.

To manually make genomes available to bcbio-nextgen, edit the individual
``.loc`` files with locations to your reference and index genomes. You need to
edit ``sam_fa_indices.loc`` to point at the FASTA files and then any genome
indexes corresponding to aligners you'd like to use (for example:
``bwa_index.loc`` for bwa and ``bowtie2_indices.loc`` for bowtie2). The database
key names used (like ``GRCh37`` and ``mm10``) should match those used in the
``genome_build`` of your sample input configuration file.

.. _Galaxy .loc files: http://wiki.galaxyproject.org/Admin/NGS%20Local%20Setup

Adding custom genomes
~~~~~~~~~~~~~~~~~~~~~~
``bcbio_setup_genome.py`` will help you to install a custom genome and apply all changes needed
to the configuration files. It needs the genome in FASTA format, and the annotation file
in GTF or GFF3 format. It can create index for all aligners used by bcbio. Moreover, it will create
the folder `rnaseq` to allow you run the RNAseq pipeline without further configuration.

::

    bcbio_setup_genome.py -f genome.fa -g annotation.gtf -i bowtie2 star seq -n Celegans -b WBcel135

If you want to add smallRNA-seq data files, you will need to add the 3 letters code of mirbase
for your genome (i.e hsa for human) and the GTF file for the annotation of smallRNA data.
Here you can use the same file than the transcriptome if no other available.

::

    bcbio_setup_genome.py -f genome.fa -g annotation.gtf -i bowtie2 star seq -n Celegans -b WBcel135 --species cel --srna_gtf another_annotation.gtf

To use that genome just need to configure your YAML files as::

    genome_build: WBcel135
