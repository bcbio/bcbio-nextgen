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

Development infrastructure
==========================

bcbio-nextgen uses GitHub for code development, and we welcome
pull requests. GitHub makes it easy to establish custom forks of the
code and contribute those back. The Biopython documentation has great
information on `using git and GitHub`_ for a community developed
project.

The automated bcbio-nextgen installer creates an isolated Python
environment using `Anaconda`_. This will be a subdirectory of your
installation root, like ``/usr/local/share/bcbio_nextgen/anaconda``.
You can manually use the ``bin/python`` within this subdirectory, or
setup a Python virtual environment management system like
`virtualenv-burrito`_ or `Python Env Wrapper`_ to switch in and out of
this environment.

You generally will want to make changes to your local copy of the
bcbio-nextgen code and then install these into the code directory
using ``/path/to/anaconda/bin/python setup.py install``. One tricky
part that I don't yet know how to work around is that pip and standard
``setup.py install`` have different ideas about how to write Python
eggs. ``setup.py install`` will create an isolated python egg
directory like ``bcbio_nextgen-0.7.5a-py2.7.egg``, while pip creates
an egg pointing to a top level ``bcbio`` directory. Where this gets
tricky is that the top level ``bcbio`` directory takes precedence. The
best way to work around this problem is to manually remove the current
pip installed bcbio-nextgen code (``rm -rf /path/to/anaconda/lib/python2.7/site-packages/bcbio*``)
before managing it manually with ``python setup.py install``. We'd
welcome tips about ways to force consistent installation across
methods.

.. _using git and GitHub: http://biopython.org/wiki/GitUsage
.. _Anaconda: http://docs.continuum.io/anaconda/index.html
.. _virtualenv-burrito: https://github.com/brainsik/virtualenv-burrito
.. _Python Env Wrapper: https://github.com/berdario/invewrapper

Adding tools
============

Aligner
~~~~~~~
Write new aligners within their own submodule inside the ``ngsalign``
directory. `bwa.py`_ is a good example to follow along with. There are
two functions to implement, based on which type of alignment you'd
like to allow:

- ``align_bam`` -- Performs alignment given an input BAM file.
  Expected to return a sorted BAM output file.

- ``align`` -- Performs alignment given FASTQ inputs (gzipped or not). This is
  generally expected to implement an approach with unix-pipe that minimizes
  intermediates and disk IO, returning a sorted BAM output file. For
  back-compatibility this can also return a text based SAM file.

See the :ref:`names-codedetails` section for more details on arguments.

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
- ``genome_resources`` -- Naming aliases and associated files
  associated with the current genome build. Retrieved from organism
  specific configuration files (``buildname-resources.yaml``) this
  specifies the location of supplemental organism specific files like
  support files for variation and RNA-seq analysis.

It also contains information the genome build, sample name and
reference genome file throughout. Here's an example of these inputs::

    {'config': {'algorithm': {'aligner': 'bwa',
                              'callable_regions': 'analysis_blocks.bed',
                              'coverage_depth': 'low',
                              'coverage_interval': 'regional',
                              'mark_duplicates': 'samtools',
                              'nomap_split_size': 50,
                              'nomap_split_targets': 20,
                              'num_cores': 1,
                              'platform': 'illumina',
                              'quality_format': 'Standard',
                              'realign': 'gkno',
                              'recalibrate': 'gatk',
                              'save_diskspace': True,
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
    'genome_resources': {'aliases': {'ensembl': 'human',
                                      'human': True,
                                      'snpeff': 'hg19'},
                          'rnaseq': {'transcripts': '/path/to/rnaseq/ref-transcripts.gtf',
                                     'transcripts_mask': '/path/to/rnaseq/ref-transcripts-mask.gtf'},
                          'variation': {'dbsnp': '/path/to/variation/dbsnp_132.vcf',
                                        'train_1000g_omni': '/path/to/variation/1000G_omni2.5.vcf',
                                        'train_hapmap': '/path/to/hg19/variation/hapmap_3.3.vcf',
                                        'train_indels': '/path/to/variation/Mills_Devine_2hit.indels.vcf'},
                          'version': 1},
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
                 'project': 'project-summary.yaml'},
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

Parallelization framework
=========================

bcbio-nextgen supports parallel runs on local machines using multiple cores and
distributed on a cluster using IPython using a general framework.

The first parallelization step starts up a set of resources for processing. On a
cluster this spawns a IPython parallel controller and set of engines for
processing. The `prun (parallel run)`_ ``start`` function is the entry point to
spawning the cluster and the main argument is a ``parallel`` dictionary which
contains arguments to the engine processing command. Here is an example input
from an IPython parallel run::

    {'cores': 12,
     'type': 'ipython'
     'progs': ['aligner', 'gatk'],
     'ensure_mem': {'star': 30, 'tophat': 8, 'tophat2': 8},
     'module': 'bcbio.distributed',
     'queue': 'batch',
     'scheduler': 'torque',
     'resources': [],
     'retries': 0,
     'tag': '',
     'timeout': 15}

The ``cores`` and ``type`` arguments must be present, identifying the total
cores to use and type of processing, respectively. Following that are arguments
to help identify the resources to use. ``progs`` specifies the programs used,
here the aligner, which bcbio looks up from the input sample file, and
gatk. ``ensure_mem`` is an optional argument that specifies minimum memory
requirements to programs if used in the workflow. The remaining
arguments are all specific to IPython to help it spin up engines on the
appropriate computing cluster.

A shared component of all processing runs is the identification of used programs
from the ``progs`` argument. The run creation process looks up required memory
and CPU resources for each program from the :ref:`config-resources` section of
your ``bcbio_system.yaml`` file. It combines these resources into required
memory and cores using the logic described in the :ref:`memory-management`
section of the parallel documentation. Passing these requirements to the cluster
creation process ensures the available machines match program requirements.

bcbio-nextgen's `pipeline.main`_ code contains examples of starting and using
set of available processing engines. This example starts up machines that use
samtools, gatk and cufflinks then runs an RNA-seq expression analysis::

    with prun.start(_wprogs(parallel, ["samtools", "gatk", "cufflinks"]),
                    samples, config, dirs, "rnaseqcount") as run_parallel:
        samples = rnaseq.estimate_expression(samples, run_parallel)

The pipelines often reuse a single set of machines for multiple distributed
functions to avoid the overhead of starting up and tearing down machines and
clusters.

The ``run_parallel`` function returned from the ``prun.start`` function enables
running on jobs in the parallel on the created machines. The `ipython wrapper`_
code contains examples of implementing this. It is a simple function that takes
two arguments, the name of the function to run and a set of multiple arguments
to pass to that function::

    def run(fn_name, items):

The ``items`` arguments need to be strings, lists and dictionaries to allow
serialization to JSON format. The internals of the run function take care of
running all of the code in parallel and returning the results back to the caller
function.

In this setup, the main processing code is fully independent from the parallel
method used so running on a single multicore machine or in parallel on a cluster
return identical results and require no changes to the logical code defining the
pipeline.

During re-runs, we avoid the expense of spinning up processing clusters for
completed tasks using simple checkpoint files in the ``checkpoints_parallel``
directory. The ``prun.start`` wrapper writes these on completion of processing
for a group of tasks with the same parallel architecture, and on subsequent runs
will go through these on the local machine instead of parallelizing. The
processing code supports these quick re-runs by checking for and avoiding
re-running of tasks when it finds output files.

Plugging new parallelization approaches into this framework involves writing
interface code that handles the two steps. First, create a cluster of ready to
run machines given the ``parallel`` function with expected core and memory
utilization:

- ``num_jobs`` -- Total number of machines to start.
- ``cores_per_job`` -- Number of cores available on each machine.
- ``mem`` -- Expected memory needed for each machine. Divide by ``cores_per_job`` to
  get the memory usage per core on a machine.

Second, implement a ``run_parallel`` function that handles using these resources
to distribute jobs and return results. The `multicore wrapper`_ and
`ipython wrapper`_ are useful starting points for understanding the current
implementations.

.. _prun (parallel run): https://github.com/chapmanb/bcbio-nextgen/blob/master/bcbio/distributed/prun.py
.. _pipeline.main: https://github.com/chapmanb/bcbio-nextgen/blob/master/bcbio/pipeline/main.py
.. _ipython wrapper: https://github.com/chapmanb/bcbio-nextgen/blob/master/bcbio/distributed/ipython.py
.. _multicore wrapper: https://github.com/chapmanb/bcbio-nextgen/blob/master/bcbio/distributed/multi.py
