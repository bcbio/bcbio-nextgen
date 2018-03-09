Common Workflow Language (CWL)
------------------------------

bcbio supports running with `Common Workflow Language (CWL)
<https://github.com/common-workflow-language/common-workflow-language>`_
compatible parallelization software. bcbio generates a CWL workflow from a
`sample YAML description file
<https://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html>`_.
Any tool that supports CWL input can run this workflow. CWL-based tools do the
work of managing files and workflows, and bcbio performs the biological analysis
using either a Docker container or a local installation.

This is a work in progress and not yet a complete production implementation. The
documentation orients anyone interested in helping with development.

Current status
~~~~~~~~~~~~~~

bcbio currently supports creation of CWL for alignment, small variant
calls (SNPs and indels), coverage assessment, HLA typing and quality
control. It generates a `CWL v1.0.2 <http://www.commonwl.org/v1.0/>`_ compatible
workflow. The actual biological code execution during runs works with
either the bcbio docker container
(`bcbio/bcbio <https://hub.docker.com/r/bcbio/bcbio/>`_) or a local
installation of bcbio.

The implementation includes bcbio's approaches to splitting and batching
analyses. At the top level workflow, we parallelize by samples. Using
sub-workflows, we split fastq inputs into sections for parallel alignment over
multiple machines following by merging. We also use sub-workflows, along with
CWL records, to batch multiple samples and run in parallel. This enables pooled
and tumor/normal cancer calling with parallelization by chromosome regions based
on coverage calculations.

.. figure:: http://i.imgur.com/iyU8VIZ.png
   :width: 600
   :height: 700
   :align: center
   :alt: Variant calling overview

bcbio supports these CWL-compatible tools:

- `toil <https://github.com/BD2KGenomics/toil>`_ -- parallel local and
  distributed cluster runs on schedulers like SLURM, SGE and PBSPro.

- `rabix bunny <https://github.com/rabix/bunny>`_ -- multicore local runs.

- `Arvados <https://arvados.org/>`_ -- fully parallel distributed analyses. We
  include an example below of running on the `public Curoverse
  <https://cloud.curoverse.com/>`_ instance running on
  `Microsoft Azure <https://azure.microsoft.com>`_.

- `Seven Bridges <https://www.sevenbridges.com/>`_ -- parallel distributed
  analyses on the Seven Bridges platform and `Cancer Genomics Cloud
  <http://www.cancergenomicscloud.org/>`_.

- `cwltool <https://github.com/common-workflow-language/cwltool>`_ -- a single
  core analysis engine, primarily used for testing.

We plan to continue to expand CWL support to include more components of bcbio,
and also need to evaluate the workflow on larger, real life analyses. This
includes supporting additional CWL runners. We're working on supporting
`DNAnexus <https://www.dnanexus.com/>`_, evaluating `Galaxy/Planemo
<https://github.com/galaxyproject/planemo>`_ for integration with the Galaxy
community, and generating inputs for `Broad's Cromwell WDL runner
<http://gatkforums.broadinstitute.org/wdl/discussion/8454/feedback-on-initial-version-of-bcbio-wdl-converted-from-cwl>`_.

Getting started
~~~~~~~~~~~~~~~

`bcbio-vm <https://github.com/bcbio/bcbio-nextgen-vm>`_ installs all
dependencies required to generate CWL and run bcbio, along with supported CWL
runners. To install using `Miniconda <http://conda.pydata.org/miniconda.html>`_
and `bioconda packages <https://bioconda.github.io/>`_::

    wget http://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
    bash Miniconda2-latest-Linux-x86_64.sh -b -p ~/install/bcbio-vm/anaconda
    ~/install/bcbio-vm/anaconda/bin/conda install --yes -c conda-forge -c bioconda bcbio-nextgen-vm
    ln -s ~/install/bcbio-vm/anaconda/bin/bcbio_vm.py /usr/local/bin/bcbio_vm.py
    ln -s ~/install/bcbio-vm/anaconda/bin/conda /usr/local/bin/bcbiovm_conda

If you have `Docker <https://www.docker.com/>`_ present on your system this is
all you need to get started running examples. If you instead prefer to use a
local installation, `install bcbio
<https://bcbio-nextgen.readthedocs.io/en/latest/contents/installation.html#automated>`_
and make it available in your path. To only run the tests, you don't need a full
data installation so can install with ``--nodata``.

To make it easy to get started, we have a pre-built CWL description that
uses test data. This will run in under 5 minutes on a local machine and
doesn't require a bcbio installation if you have Docker available on
your machine:

1. Download and unpack the `test repository <https://github.com/bcbio/test_bcbio_cwl>`_::

     wget -O test_bcbio_cwl.tar.gz https://github.com/bcbio/test_bcbio_cwl/archive/master.tar.gz
     tar -xzvpf test_bcbio_cwl.tar.gz
     cd test_bcbio_cwl-master/somatic

2. Run the analysis using either Toil or Rabix bunny. If you have Docker
   available on your machine, the runner will download the correct `bcbio
   container <https://github.com/bcbio/bcbio_docker>`_ and you don't need to
   install anything else to get started. If you have an old version of the
   container you want to update to the latest with ``docker pull
   quay.io/bcbio/bcbio-vc``. There are shell scripts that provide the command
   lines for running::

     bash run_toil.sh
     bash run_bunny.sh

   Or you can run directly using the ``bcbio_vm.py`` wrappers::

     bcbio_vm.py cwlrun toil somatic-workflow
     bcbio_vm.py cwlrun bunny somatic-workflow

   Thes wrappers automatically handle temporary directories, permissions,
   logging and re-starts. If running without Docker, use a `local installation of
   bcbio
   <https://bcbio-nextgen.readthedocs.org/en/latest/contents/installation.html>`_
   add ``--no-container`` to the commands in the shell scripts.

Generating CWL for local or cluster runs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first step in running your analysis project in bcbio is to generate CWL. The
inputs to this are:

- A `standard bcbio sample configuration file
  <https://bcbio-nextgen.readthedocs.io/en/latest/contents/configuration.html>`_
  defining the samples. This can either be a full prepared YAML file or a
  `template file and CSV with sample data <http://bcbio-nextgen.readthedocs.io/en/latest/contents/configuration.html#automated-sample-configuration>`_.

- A ``bcbio_system.yaml`` file defining the system environment for running the
  program. This includes the resource specification with `cores and memory per
  core for your machines
  <http://bcbio-nextgen.readthedocs.io/en/latest/contents/configuration.html#resources>`_.
  You generally want to set this to match the parameters of a single machine
  either for a local run or on a cluster. It also includes paths to the
  reference biodata and optionally input files if you want to avoid specifying
  full paths in your inputs. Here is an example for a 16 core machine with 3.5Gb
  of memory per core::

      local:
        ref: /path/to/bcbio/genomes/Hsapiens
        inputs:
          - /path/to/input/files
      resources:
        default:
          cores: 16
          memory: 3500M
          jvm_opts: [-Xms1g, -Xmx3500m]

Generate CWL with::

    bcbio_vm.py template --systemconfig bcbio_system.yaml template.yaml samples.csv
    bcbio_vm.py cwl --systemconfig bcbio_system.yaml samples/config/samples.yaml

producing a ``sample-workflow`` output directory with the CWL. You can run this
with any CWL compatible runner. The ``bcbio_vm.py cwlrun`` wrappers described
above make this easier for local runs with Toil or Bunny.

Running bcbio CWL on Toil
~~~~~~~~~~~~~~~~~~~~~~~~~

The `Toil pipeline management system <https://github.com/BD2KGenomics/toil>`_
runs CWL workflows in parallel on a local machine, on a cluster or at AWS.
Toil comes pre-installed with bcbio-vm.

To run a bcbio CWL workflow locally with Toil using Docker::

    bcbio_vm.py cwlrun toil sample-workflow

If you want to run from a locally installed bcbio add ``--no-container`` to the
commandline.

To run distributed on a Slurm cluster::

    bcbio_vm.py cwlrun toil sample-workflow -- --batchSystem slurm

Running bcbio CWL on Arvados
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We're actively testing bcbio generated CWL workflows on
`Arvados <https://arvados.org/>`_. These instructions detail how to run
on the `Arvdos public instance <https://cloud.curoverse.com/>`_.
`Arvados cwl-runner <https://github.com/curoverse/arvados>`_ comes
pre-installed with
`bcbio-vm <https://github.com/bcbio/bcbio-nextgen-vm#installation>`_.

Retrieve API keys from the `Arvados public
instance <https://cloud.curoverse.com/>`_. Login, then go to `'User
Icon -> Personal Token' <https://cloud.curoverse.com/current_token>`_.
Copy and paste the commands given there into your shell. You'll
specifically need to set ``ARVADOS_API_HOST`` and ``ARVADOS_API_TOKEN``.

To run an analysis:

1. Create a new project from the web interface (Projects -> Add a new
   project). Note the project ID from the URL of the project (an
   identifier like ``qr1hi-j7d0g-7t73h4hrau3l063``).

2. Upload reference data to Arvados Keep. Note the genome collection
   UUID::

     arv-put --name hg19-testdata --project-uuid $PROJECT_ID testdata/genomes

3. Upload input data to Arvados Keep. Note the collection UUID::

     arv-put --name input-testdata --project-uuid $PROJECT_ID testdata/100326_FC6107FAAXX testdata/automated testdata/reference_material

4. Create an Arvados section in a ``bcbio_system.yaml`` file specifying
   locations to look for reference and input data. ``input`` can be one or more
   collections containing files or associated files in the original sample YAML::

     arvados:
       reference: qr1hi-4zz18-kuz1izsj3wkfisq
       input: [qr1hi-j7d0g-h691y6104tlg8b4]
     resources:
       default: {cores: 4, memory: 2G, jvm_opts: [-Xms750m, -Xmx2500m]}

5. Generate the CWL to run your samples. If you're using multiple input
   files with a `CSV metadata file and template <https://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html#automated-sample-configuration>`_
   start with creation of a configuration file::

     bcbio_vm.py template --systemconfig bcbio_system_arvados.yaml testcwl_template.yaml testcwl.csv

   To generate the CWL from the system and sample configuration files::

     bcbio_vm.py cwl --systemconfig bcbio_system_arvados.yaml testcwl/config/testcwl.yaml

6. Import bcbio Docker image to your Arvados project::

     docker pull quay.io/bcbio/bcbio-vc
     arv-keepdocker --project-$PROJECT_ID -- quay.io/bcbio/bcbio-vc latest

7. Run the CWL on the Arvados public cloud using the Arvados cwl-runner::

     bcbio_vm.py cwlrun arvados arvados_testcwl-workflow -- --project-uuid qr1hi-your-projectuuid

Running bcbio CWL on DNAnexus
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bcbio runs on the `DNAnexus platform <https://www.dnanexus.com/>`_ by converting
bcbio generated CWL into DNAnexus workflows and apps using
`dx-cwl <https://github.com/dnanexus/dx-cwl>`_. This describes the process
using the
'Create and Run bcbio workflow applet <https://platform.dnanexus.com/projects/F541fX00f5v9vKJjJ34gvgbv/data/applets>`_
in the public `bcbio_resources
<https://platform.dnanexus.com/projects/F541fX00f5v9vKJjJ34gvgbv/data/>`_
project, Secondarily, we also show how to install and prepare things locally for
additional control and debugging.

0. Set some useful environmental variables:

   - ``$PNAME`` -- The name of the project you're analyzing. For convenience
     here we keep this the same for your local files and remote DNAnexus
     project, although that does not have to be true.
   - ``$DX_AUTH_TOKEN`` -- The DNAnexus authorization token for access, used for
     the ``dx`` command line tool and bcbio scripts.
   - ``$DX_PROJECT_ID`` -- The DNAnexus GUID identifier for your project
      (similar to ``project-F8Q7fJj0XFJJ3XbBPQYXP4B9``). You can get this from
      ``dx env`` after creating/selecting a project in steps 1 and 2.

1. Create an analysis project::

     dx new project $PNAME

2. Upload sample data to the project::

     dx select $PNAME
     dx upload -p --path /data/input *.bam

3. Create bcbio system file with projects, locations of files and
   desired core and memory usage for jobs. bcbio uses the core and memory
   specifications to ::

     dnanexus:
       project: PNAME
       ref:
         project: bcbio_resources
         folder: /reference_genomes
       inputs:
         - /data/input
         - /data/input/regions
     resources:
       default: {cores: 8, memory: 3500M, jvm_opts: [-Xms1g, -Xmx3500m]}

4. Create bcbio sample YAML file referencing samples to run. The files can be
   relative to the ``inputs`` directory specified above; bcbio will search
   recursively for files, so you don't need to specify full paths if your file
   names are unique. Start with a template and sample specification::

       samplename,description,batch,phenotype
       file1.bam,sample1,b1,tumor
       file2.bam,sample2,b1,normal
       file3.bam,sample3,b2,tumor
       file4.bam,sample4,b2,normal

5. Pick a template file that describes the `bcbio configuration
   <http://bcbio-nextgen.readthedocs.io/en/latest/contents/configuration.html>`_
   variables. You can define parameters either globally (in the template) file
   or by sample (in the csv) using the `standard bcbio templating
   <http://bcbio-nextgen.readthedocs.io/en/latest/contents/configuration.html#automated-sample-configuration>`_.
   An example template for GATK4 germline variant calling is::

      details:
       - algorithm:
           aligner: bwa
           variantcaller: gatk-haplotype
         analysis: variant2
         genome_build: hg38

6. Supply the three inputs (``bcbio_system.yaml``, ``project.csv`` and
   ``template.yaml``) to the `Create and run bcbio workflow` applet
   <https://platform.dnanexus.com/projects/F541fX00f5v9vKJjJ34gvgbv/data/applets>`_.
   You can do this using the web interface or via the command line with a small
   script like::

      dx select $DX_PROJECT_ID
      dx mkdir -p $PNAME
      for F in your-template.yaml $PNAME.csv bcbio_system-dnanexus.yaml
      do
              dx rm -a /$PNAME/$F || true
              dx upload --path /$PNAME/ $F

      done
      dx ls $PNAME
      dx rm -a -r /$PNAME/dx-cwl-run || true
      dx run bcbio_resources:/applets/bcbio-run-workflow -iyaml_template=/$PNAME/qc-template.yaml -isample_spec=/$PNAME/$PNAME.csv -isystem_configuration=/$PNAME/bcbio_system-dnanexus.yaml -ioutput_folder=/$PNAME/dx-cwl-run

The applet will lookup all files, prepare a bcbio CWL workflow, convert into a
DNAnexus workflow, and submit to the platform. The workflow runs as a standard
DNAnexus workflow and you can monitor through the command line (with ``dx find
executions --root job-YOURJOBID`` and ``dx watch``) or the web interface
(``Monitor`` tab).

If you prefer not to use the DNAnexus app you can run locally by installing
`bcbio-vm <https://github.com/bcbio/bcbio-nextgen-vm#installation>`_ on your
local machine:


1. Follow the :ref:`automated-sample-config` workflow to generate a full configuration::

       bcbio_vm.py template --systemconfig bcbio_system-dnanexus.yaml your-template.yaml $PNAME.csv

2. Generate a CWL description of the workflow from the full generated configuration::

       bcbio_vm.py cwl --systemconfig bcbio_system-dnanexus.yaml $PNAME/config/$PNAME.yaml

3. Determine project information and login credentials. You'll want to note the
   ``Auth token used`` and ``Current workspace`` project ID::

       dx env

4. Compile the CWL workflow into a DNAnexus workflow::

       dx-cwl compile-workflow $PNAME-workflow/main-$PNAME.cwl --project PROJECT_ID --token $DX_AUTH_TOKEN

5. Upload sample information from generated CWL and run workflow::

       dx mkdir -p $DX_PROJECT_ID:/$PNAME-workflow
       dx upload -p --path $DX_PROJECT_ID:/$PNAME-workflow $PNAME-workflow/main-$PNAME-samples.json
       dx-cwl run-workflow /dx-cwl-run/main-$PNAME/main-$PNAME \
              /$PNAME-workflow/main-$PNAME-samples.json \
              --project PROJECT_ID --token $DX_AUTH_TOKEN

Development notes
~~~~~~~~~~~~~~~~~

bcbio generates a common workflow language description. Internally,
bcbio represents the files and information related to processing as `a
comprehensive
dictionary <https://bcbio-nextgen.readthedocs.org/en/latest/contents/code.html#data>`_.
This world object describes the state of a run and associated files, and
new processing steps update or add information to it. The world object
is roughly equivalent to CWL's JSON-based input object, but CWL enforces
additional annotations to identify files and models new inputs/outputs
at each step. The work in bcbio is to move from our laissez-faire
approach to the more structured CWL model.

The generated CWL workflow is in ``run_info-cwl-workflow``:

-  ``main-*.cwl`` -- the top level CWL file describing the workflow
   steps
-  ``main*-samples.json`` -- the flattened bcbio world structure
   represented as CWL inputs
-  ``wf-*.cwl`` -- CWL sub-workflows, describing sample level parallel
   processing of a section of the workflow, with potential internal
   parallelization.
-  ``steps/*.cwl`` -- CWL descriptions of sections of code run inside
   bcbio. Each of these are potential parallelization points and make up
   the nodes in the workflow.

To help with defining the outputs at each step, there is a
``WorldWatcher`` object that can output changed files and world
dictionary objects between steps in the pipeline when running a bcbio in
the standard way. The `variant
pipeline <https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/pipeline/main.py>`_
has examples using it. This is useful when preparing the CWL definitions
of inputs and outputs for new steps in the `bcbio CWL step
definitions <https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/cwl/workflow.py>`_.

ToDo
~~~~

-  Support the full variant calling workflow with additional steps like
   ensemble calling, structural variation, heterogeneity detection and
   disambiguation.

-  Port RNA-seq and small RNA workflows to CWL.

-  Replace the custom python code in the `bcbio step
   definitions <https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/cwl/defs.py>`_
   with a higher level DSL in YAML we can parse and translate to CWL.
