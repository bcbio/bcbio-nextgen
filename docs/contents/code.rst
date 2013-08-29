Code
----
This section provides useful concepts for getting started digging into
the code and contributing new functionality. We welcome contributors
and hope these notes help make it easier to get started.

Development goals
=================

During development we seek to maximize functionality and usefulness,
while avoiding complexity. Since these goals are sometimes in
conflict, it's useful to understand the design approaches:

- Support high level configurability but avoid exposing all program
  options. Since pipelines support a wide variety of tools, each with
  a large number of options, we try to define configuration variables
  at high level based on biological intent and then translate these
  into best-practice options for each tool. The goal is to avoid
  having an overwhelming number of input configuration options.

- Provide best-practice pipelines that make recommended decisions for
  processing. Coupled with goal of minimizing configuration
  parameters, this requires trust and discussion around algorithm
  choices. An example is bwa alignment, which uses ``bwa aln`` for
  reads shorter than 75bp and ``bwa mem`` for longer reads, based on
  recommendations from Heng Li. Our general goal is to encourage
  discussion and development of best-practices to make it easy to do
  the right thing.

- Support extensive debugging output. In complex distributed systems,
  programs fail in unexpected ways even during production runs. We try
  to maximize logging to help identify and diagnose these type of
  unexpected problems.

- Avoid making mistakes. This results in being conservative about
  decisions like deleting file intermediates. Coupled with extensive
  logging, we trade off disk usage for making it maximally
  easy to restart and debug problems. If you'd like to delete work or
  log directories automatically, we recommend doing this as part of
  your batch scripts wrapping bcbio-nextgen.

- Strive for a clean, readable code base. We strive to make the code a
  secondary source of information after hand written docs.
  Practically, this means maximizing information content in source
  files while using in-line documentation to clarify as needed.

- Focus on a functional coding style with minimal use of global
  mutable objects. This approach works well with distributed code and
  isolates debugging to individual functions rather than globally
  mutable state.

Overview
========

The most useful modules inside ``bcbio``, ordered by likely interest:

- ``pipeline`` -- Top level functionality that drives the analysis
  pipeline. ``main.py`` contains top level definitions of pipelines
  like variant calling and RNAseq, and is the best place to start
  understanding the overall organization of the code.
- ``ngsalign`` -- Integration with aligners for high-throughput
  sequencing data. We support individual aligners with their own
  separate modules.
- ``variation`` -- Tools for variant calling. Individual variant
  calling and processing approaches each have their own submodules.
- ``rnaseq`` -- Run RNA-seq pipelines, currently supporting TopHat/Cufflinks.
- ``provenance`` -- Track third party software versions, command lines
  and program flow. Handle writing of debugging details.
- ``distributed`` -- Handle distribution of programs across multiple
  cores, or across multiple machines using IPython.
- ``workflow`` -- Provide high level tools to run customized analyses.
  They tie into specialized analyses or visual front ends to make
  running bcbio-nextgen easier for specific common tasks.
- ``broad`` -- Code to handle calling Broad tools like GATK and
  Picard, as well as other Java-based programs.

Adding tools
============

Aligner
~~~~~~~
Write new aligners within their own submodule inside the ``ngsalign``
directory. `bwa.py`_ is a good example to follow along with. There are
three functions to implement, based on which type of alignment you'd
like to allow:

- ``align_bam`` -- Performs alignment given an input BAM file.
  Expected to return a sorted BAM output file.

- ``align_pipe`` -- Performs alignment given FASTQ inputs (gzipped or
  not). Expected to implemented an approach with unix-pipe that
  minimizes intermediates and disk IO. Expected to return a sorted BAM
  output file.

- ``align`` -- Performs alignment given FASTQ inputs, returning a text
  based SAM file.

``align_bam`` and ``align_pipe`` are most commonly used now, which
``align`` provides functionality for older aligners that do not easily
support a piped approach. See the :ref:`names-codedetails` section for more
details on arguments.

Other required implementation details include:

- ``galaxy_loc_file`` -- Provides the name of the `Galaxy loc file`_
  used to identify locations of indexes for this aligner. The
  automated installer sets up these loc files automatically.

- ``remap_index_fn`` -- A function that remaps an index from the
  Galaxy location file into the exact one for this aligner. This is
  useful for tools which aren't supported by a Galaxy .loc file but
  you can locate them relative to another index.

.. _bwa.py: https://github.com/chapmanb/bcbio-nextgen/blob/master/bcbio/ngsalign/bwa.py
.. _Galaxy loc file: http://wiki.galaxyproject.org/Admin/Data%20Integration

Once implemented, plug the aligner into the pipeline by defining it as
a ``_tool`` in `bcbio/pipeline/alignment.py`_. You can then use it as
normal by specifying the name of the aligner in the `aligner` section
of your configuration input.

.. _bcbio/pipeline/alignment.py: https://github.com/chapmanb/bcbio-nextgen/blob/master/bcbio/pipeline/alignment.py

Variant caller
~~~~~~~~~~~~~~

New variant calling approaches live within their own module inside
``bcbio/variation``. The `freebayes.py`_ implementation is a good
example to follow for providing your own variant caller. Implement a
function to run variant calling on multiple BAMs in an input region
that takes the following inputs:

- ``align_bams`` -- A list of BAM files to call simultaneously.
- ``items`` -- List of ``data`` dictionaries associated with each of the
  samples in ``align_bams``. Enables customization of variant calling
  based on sample configuration inputs. See documentation on the
  :ref:`data-codedetails` dictionary for all of the information
  contained inside each ``data`` item. Having multiple
  configurations allows customization of sample specific variant calls
  using parameters supplied to :ref:`sample-configuration`.
- ``ref_file`` -- Fasta reference genome file.
- ``assoc_files`` -- Useful associated files for variant calling. This
  includes the DbSNP VCF file. It's a named tuple mapping to files
  specified in the configuration. `bcbio/pipeline/shared.py`_ has the
  available inputs.
- ``region`` -- A tuple of (chromosome, start, end) specifying the
  region to call in.
- ``out_file``-- The output file to write to. This should contain calls
  for all input samples in the supplied region.

Once implemented, add the variant caller into the pipeline by updating
``caller_fns`` in the ``variantcall_sample`` function in
`bcbio/variation/genotype.py`_. You can use it by specifying it in the
``variantcaller`` parameter of your sample configuration.

.. _freebayes.py: https://github.com/chapmanb/bcbio-nextgen/blob/master/bcbio/variation/freebayes.py
.. _bcbio/variation/genotype.py: https://github.com/chapmanb/bcbio-nextgen/blob/master/bcbio/variation/genotype.py#L548
.. _bcbio/pipeline/shared.py: https://github.com/chapmanb/bcbio-nextgen/blob/master/bcbio/pipeline/shared.py#L176

Standard function arguments
===========================

.. _names-codedetails:

names
~~~~~
This dictionary provides lane and other `BAM run group`_ naming
information used to correctly build BAM files. We use the ``rg``
attribute as the ID within a BAM file::

    {'lane': '7_100326_FC6107FAAXX',
     'pl': 'illumina',
     'pu': '7_100326_FC6107FAAXX',
     'rg': '7',
     'sample': 'Test1'}

.. _BAM run group: http://samtools.sourceforge.net/SAM1.pdf

.. _data-codedetails:

data
~~~~

The `data` dictionary is a large dictionary representing processing,
configuration and files associated with a sample. The standard
work flow is to pass this dictionary between functions, updating with
associated files from the additional processing. Populating this
dictionary only with standard types allows serialization to JSON for
distributed processing.

The dictionary is dynamic throughout the workflow depending on the
step, but some of the most useful key/values available throughout are:

- ``config`` -- Input configuration variables about how to process in
  the ``algorithm`` section and locations of programs in the ``resources``
  section.
- ``dirs`` -- Useful directories for building output files or retrieving
  inputs.
- ``metadata`` -- Top level metadata associated with a sample, specified
  in the initial configuration.

It also contains information the genome build, sample name and
reference genome file throughout. Here's an example of these inputs::

    {'config': {'algorithm': {'aligner': 'bwa',
                              'callable_regions': 'analysis_blocks.bed',
                              'coverage_depth': 'low',
                              'coverage_interval': 'regional',
                              'dbsnp': 'variation/dbsnp_132.vcf',
                              'mark_duplicates': 'samtools',
                              'max_errors': 2,
                              'nomap_split_size': 50,
                              'nomap_split_targets': 20,
                              'num_cores': 1,
                              'platform': 'illumina',
                              'quality_format': 'Standard',
                              'realign': 'gkno',
                              'recalibrate': 'gatk',
                              'save_diskspace': True,
                              'train_1000g_omni': 'variation/1000G_omni2.5.vcf',
                              'train_hapmap': 'variation/hapmap_3.3.vcf',
                              'train_indels': 'variation/Mills_Devine_2hit.indels.vcf',
                              'upload_fastq': False,
                              'validate': '../reference_material/7_100326_FC6107FAAXX-grade.vcf',
                              'variant_regions': '../data/automated/variant_regions-bam.bed',
                              'variantcaller': 'freebayes'},
                'resources': {'bcbio_variation': {'dir': '/usr/share/java/bcbio_variation'},
                              'bowtie': {'cores': None},
                              'bwa': {'cores': 4},
                              'cortex': {'dir': '~/install/CORTEX_release_v1.0.5.14'},
                              'cram': {'dir': '/usr/share/java/cram'},
                              'gatk': {'cores': 2,
                                       'dir': '/usr/share/java/gatk',
                                       'jvm_opts': ['-Xms750m', '-Xmx2000m'],
                                       'version': '2.4-9-g532efad'},
                              'gemini': {'cores': 4},
                              'novoalign': {'cores': 4,
                                            'memory': '4G',
                                            'options': ['-o', 'FullNW']},
                              'picard': {'cores': 1,
                                         'dir': '/usr/share/java/picard'},
                              'snpEff': {'dir': '/usr/share/java/snpeff',
                                         'jvm_opts': ['-Xms750m', '-Xmx3g']},
                              'stampy': {'dir': '~/install/stampy-1.0.18'},
                              'tophat': {'cores': None},
                              'ucsc_bigwig': {'memory': '36g'},
                              'varscan': {'dir': '/usr/share/java/varscan'},
                              'vcftools': {'dir': '~/install/vcftools_0.1.9'}}},
     'dirs': {'fastq': 'input fastq directory',
              'galaxy': 'directory with galaxy loc and other files',
              'work': 'base work directory'},
     'metadata': {'batch': 'TestBatch1'},
     'genome_build': 'hg19',
     'name': ('', 'Test1'),
     'sam_ref': '/path/to/hg19.fa'}

Processing also injects other useful key/value pairs. Here's an example of
additional information supplied during a variant calling workflow::

    {'prep_recal': 'Test1/7_100326_FC6107FAAXX-sort.grp',
     'summary': {'metrics': [('Reference organism', 'hg19', ''),
                             ('Total', '39,172', '76bp paired'),
                             ('Aligned', '39,161', '(100.0\\%)'),
                             ('Pairs aligned', '39,150', '(99.9\\%)'),
                             ('Pair duplicates', '0', '(0.0\\%)'),
                             ('Insert size', '152.2', '+/- 31.4')],
                 'pdf': '7_100326_FC6107FAAXX-sort-prep-summary.pdf',
                 'project': 'project-summary.csv'},
     'validate': {'concordant': 'Test1-ref-eval-concordance.vcf',
                  'discordant': 'Test1-eval-ref-discordance-annotate.vcf',
                  'grading': 'validate-grading.yaml',
                  'summary': 'validate-summary.csv'},
     'variants': [{'population': {'db': 'gemini/TestBatch1-freebayes.db',
                                  'vcf': None},
                   'validate': None,
                   'variantcaller': 'freebayes',
                   'vrn_file': '7_100326_FC6107FAAXX-sort-variants-gatkann-filter-effects.vcf'}],
     'vrn_file': '7_100326_FC6107FAAXX-sort-variants-gatkann-filter-effects.vcf',
     'work_bam': '7_100326_FC6107FAAXX-sort-prep.bam'}
