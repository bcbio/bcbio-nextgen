.. _docs-config:

Configuration
-------------

Two configuration files, in easy to write `YAML format`_, specify
details about your system and samples to run:

- ``bcbio_system.yaml`` High level information about the system,
  including locations of installed programs like Picard and GATK.
  These apply across multiple runs. The automated installer creates
  a ready to go system configuration file that can be manually
  edited to match the system. Find the file in the galaxy sub-directory
  within your installation data location
  (ie. ``/usr/local/share/bcbio-nextgen/galaxy``). By default, the
  pipeline uses the standard pre-created configuration file but
  multiple system configurations can be independently maintained
  and passed as the first argument to ``bcbio_nextgen.py`` commands.

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

        samplename,description,batch,phenotype,sex,coverage_interval
        sample1,ERR256785,batch1,normal,female,genome
        sample2,ERR256786,batch1,tumor,,genome

  The first column links the metadata to a specific input file. The
  template command tries to identify the ``samplename`` from read group
  information in a BAM file, or uses the base filename if no read group
  information is present.  The remaining columns can contain:

   - ``description`` Changes the sample description, originally
     supplied by the file name or BAM read group, to this value.

   - Algorithm parameters specific for this sample. If the column name matches
     an available :ref:`algorithm-config`, then this value substitutes
     into the sample ``algorithm``, replacing the defaults from the template.

   -  :ref:`sample-configuration` metadata key/value pairs. Any columns not
      falling into the above cases will go into the metadata section.

  Individual column items can contain booleans (true or false), integers, or
  lists (separated by semi-colons). These get converted into the expected time
  in the output YAML file.

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
samples, pointing to this custom configuration file.

.. _best-practice templates: https://github.com/chapmanb/bcbio-nextgen/tree/master/config/templates
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
          platform:illumina
        metadata:
          batch: Batch1
          sex: female

- ``analysis`` Analysis method to use [variant2, RNA-seq]
- ``lane`` A unique number within the project. Corresponds to the
  ``ID`` parameter in the BAM read group. Required.
- ``description`` Unique name for this sample, corresponding to the
  ``SM`` parameter in the BAM read group.
- ``files`` A list of files to process. This currently supports either a single
  end or two paired end fastq files, or a single BAM file. It does not yet
  handle merging BAM files or more complicated inputs.
- ``genome_build`` Genome build to align to, which references a genome
  keyword in Galaxy to find location build files.

- ``algorithm`` Parameters to configure algorithm inputs. Options
  described in more detail below.
- ``metadata`` Additional descriptive metadata about the sample:

    - ``batch`` defines a group that the sample falls in. We perform
       multi-sample variant calling on all samples with the same batch
       name. This can also be a list, allowing specification of a single normal
       sample to pair with multiple tumor samples in paired cancer variant
       calling (``batch: [MatchWithTumor1, MatchWithTumor2]``).

    - ``sex`` specifies the sample sex used to correctly prepare X/Y
      chromosomes.

Setting up a test run
~~~~~~~~~~~~~~~~~~~~~
The if you set the ``test_run`` option to ``True`` at the top of your sample
configuration file like this::

  test_run: True

bcbio-nextgen will downsample your input files to 500,000 entries before
running the pipeline.

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

- ``bucket`` AWS bucket to upload to
- ``access_key_id`` AWS access key ID from Amazon credentials page
- ``secret_access_key`` AWS secret key ID from Amazon credentials page
- ``reduced_redundancy`` Flag to determine if we should store S3 data
  with reduced redundancy: cheaper but less reliable [false, true]

.. _algorithm-config:

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

Algorithm parameters
~~~~~~~~~~~~~~~~~~~~

The YAML configuration file provides a number of hooks to customize
analysis in the sample configuration file. Place these under the
``algorithm`` keyword.

Alignment
=========

- ``platform`` Sequencing platform used. Corresponds to the ``PL``
  parameter in BAM read groups. Default 'Illumina'.
-  ``aligner`` Aligner to use: [bwa, bowtie, bowtie2, mosaik, novoalign, star,
   false]
-  ``bam_clean`` Clean an input BAM when skipping alignment step. This
   handles adding read groups, sorting to a reference genome and
   filtering problem records that cause problems with GATK. Set to
   ``picard`` to do Picard/GATK based cleaning.
-  ``bam_sort`` Allow sorting of input BAMs when skipping alignment
   step (``aligner`` set to false). Options are coordinate or
   queryname. For additional processing through standard pipelines
   requires coordinate sorted inputs. The default is to not do
   additional sorting and assume pre-sorted BAMs.
- ``disambiguate`` For mixed or explant samples, provide a list of
  ``genome_build``  identifiers to check and remove from alignment. Currently
  supports cleaning a single organism. For example, with ``genome_build: hg19``
  and ``disambiguate: [mm10]``, it will align to hg19 and mm10, run
  disambiguation and continue with reads confidently aligned to hg19.
-  ``trim_reads`` Can be set to trim low quality ends or to also trim off,
    in conjunction with the ``adapters`` field a set of adapter sequences or
    poly-A tails that could appear on the ends of reads. Only used in RNA-seq
    pipelines, not variant calling. [False, read_through]
- ``min_read_length`` Minimum read length to maintain when
  ``read_through`` trimming set in ``trim_reads``. Defaults to 20.
-  ``adapters`` If trimming adapter read through, trim a set of stock
   adapter sequences. Allows specification of multiple items in a list,
   for example [truseq, polya] will trim both TruSeq adapter sequences
   and polyA tails. Valid items are [truseq, illumina, nextera, polya]
-  ``custom_trim`` A list of sequences to trim from the end of reads,
   for example: [AAAATTTT, GGGGCCCC]
-  ``align_split_size``: Split FASTQ files into specified number of
   records per file. Allows parallelization at the cost of increased
   temporary disk space usage.
-  ``quality_bin``: Perform binning of quality scores with CRAM to
   reduce file sizes. Uses the Illumina 8-bin approach. Supply a list
   of times to perform binning: [prealignment, postrecal]
-  ``quality_format`` Quality format of fastq inputs [illumina,
   standard]
-  ``merge_bamprep`` Merge regional BAM prepped files into a final
   prepared BAM. false avoids the time consuming merge when you only
   want variant calls [true, false]
-  ``coverage_bigwig`` Generate a bigwig file of coverage, for loading
   into the UCSC genome browser [true, false]
-  ``strandedness`` For RNA-seq libraries, if your library is strand
   specific, set the appropriate flag form [unstranded, firststrand, secondstrand].
   Defaults to unstranded. For dUTP marked libraries, firststrand is correct; for
   Scriptseq prepared libraries, secondstrand is correct.

Experimental information
========================

-  ``coverage_interval`` Regions covered by sequencing. Influences GATK
   options for filtering. GATK will use Variant Quality Score Recalibration
   when set to 'genome', otherwise we apply hard filters. [exome, genome, regional]
- ``coverage_depth_max`` Maximum depth of coverage. We downsample coverage
   regions with more than this value to approximately the specified
   coverage. Actual coverage depth per position will be higher since we
   downsample reads based on shared start positions, although some callers like
   GATK can also downsample to exactly this coverage per position. We avoid
   calling entirely in super high depth regions with more than 7 times coverage
   for this parameter. This controls memory usage in highly repetitive regions
   like centromeres. Defaults to 10000. Set to 0 to perform no downsampling.
-  ``coverage_depth_min`` Minimum depth of coverage. Regions will less reads
   will not get called. Defaults to 4. Setting lower than 4 will trigger
   low-depth calling options for GATK.
-  ``ploidy`` Ploidy of called reads. Defaults to 2 (diploid).

Variant calling
===============

-  ``variantcaller`` Variant calling algorithm. Can be a list of
   multiple options [gatk, freebayes, varscan, samtools,
   gatk-haplotype, cortex, mutect]
-  ``variant_regions`` BED file of regions to call variants in.
-  ``mark_duplicates`` Identify and remove variants [picard,
   biobambam, samtools, false]
-  ``recalibrate`` Perform base quality score recalibration on the
   aligned BAM file. [gatk, false]
-  ``realign`` Perform realignment around indels on the aligned BAM
   file. [gatk, gkno, false]
-  ``phasing`` Do post-call haplotype phasing of variants. Defaults to
   no phasing [false, gatk]
-  ``validate`` A VCF file of expected variant calls to perform
    validation and grading of output variants from the pipeline.
    This provides a mechanism to ensure consistency of calls against
    a known set of variants, supporting comparisons to genotyping
    array data or reference materials.
- ``validate_regions`` A BED file of regions to evaluate in. This
  defines specific regions covered by the ``validate`` VCF  file.
- ``validate_genome_build``: Genome build of the validation file, if
  different than the samples genome build. Helps manage hg19/GRCh37
  chromosome naming differences.
- ``clinical_reporting`` Tune output for clinical reporting.
  Modifies snpEff parameters to use HGVS notational on canonical
  transcripts [false, true].
- ``background`` Provide a VCF file with variants to use as a background
  reference during variant calling. For tumor/normal paired calling use this to
  supply a panel of normal individuals.

Parallelization
===============

- ``nomap_split_size`` Unmapped base pair regions required to split
  analysis into blocks. Creates islands of mapped reads surrounded by
  unmapped (or N) regions, allowing each mapped region to run in
  parallel. (default: 100)

- ``nomap_split_targets`` Number of target intervals to attempt to
  split processing into. This picks unmapped regions evenly spaced
  across the genome to process concurrently. Limiting targets prevents
  a large number of small targets. (default: 2000)

Ensemble variant calling
========================

In addition to single method variant calling, we support calling with
multiple calling methods and consolidating into a final Ensemble
callset. This requires the `bcbio.variation`_ toolkit to perform the
consolidation. An example configuration in the ``algorithm`` section is::

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

.. _bcbio.variation: https://github.com/chapmanb/bcbio.variation
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

- `aliases` -- Names for third-party programs used as part of the
  analysis, since naming expectations can differ between software
  programs.

- `variation` -- Supporting data files for variant analysis. For human
  analyses, the dbSNP and training files are from the `GATK resource bundle`_.
  These are inputs into the training models for
  recalibration. The automated `CloudBioLinux`_ data scripts will
  download and install these in the variation subdirectory relative to
  the genome files.

- `rnaseq` -- Supporting data files for RNA-seq analysis. The
  automated installer and updater handles retrieval and installation
  of these resources for supported genome builds.

By default, we place the ``buildname-resources.yaml`` files next to
the genome FASTA files in the reference directory. For custom setups,
you specify an alternative directory in the ref:`config-resources`
section of your ``bcbio_system.yaml`` file::

  resources:
    genome:
      dir: /path/to/resources/files

.. _Example genome configuration files: https://github.com/chapmanb/bcbio-nextgen/tree/master/config/genomes
.. _GATK resource bundle: http://www.broadinstitute.org/gsa/wiki/index.php/GATK_resource_bundle

Reference genome files
~~~~~~~~~~~~~~~~~~~~~~

The pipeline requires access to reference genomes, including the raw
FASTA sequence and pre-built indexes for aligners. For human genomes, the
automated installer provides hg19 and GRCh37 1000 genomes references as
provided in the `GATK resource bundle`_.
:ref:`data-requirements` section describes the expected layout of
`Galaxy .loc files`_ pointing to the actual sequence and index
files.

The pipeline identifies the root ``galaxy`` directory, in which it
expects a ``tool-data`` sub-directory with the ``.loc`` files, in two
ways:

- Using the directory of your ``bcbio-system.yaml``. This is the
  default mechanism setup by the automated installer.

- From the path specified by the ``galaxy_config`` option in your
  ``bcbio-system.yaml``. If you'd like to move your system YAML file,
  add the full path to your ``galaxy`` directory here.

.. _Galaxy .loc files: http://wiki.galaxyproject.org/Admin/NGS%20Local%20Setup
