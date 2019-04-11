.. _docs-cwl:

Common Workflow Language (CWL)
------------------------------

bcbio runs with `Common Workflow Language (CWL)
<https://github.com/common-workflow-language/common-workflow-language>`_
compatible parallelization software. bcbio generates a CWL workflow from a
`standard bcbio sample YAML description file
<https://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html>`_
and any tool that supports CWL input can run the workflow. CWL-based tools do the
work of managing files and workflows, and bcbio performs the biological analysis
using either a Docker container or a local installation.

Current status
~~~~~~~~~~~~~~

bcbio creates CWL for alignment, small variant calls (SNPs and indels), coverage
assessment, HLA typing, quality control and structural variant calling. It
generates a `CWL v1.0.2 <http://www.commonwl.org/v1.0/>`_ compatible workflow.
The actual biological code execution during runs works with either a `bcbio
docker container <https://github.com/bcbio/bcbio_docker>`_ or a `local
installation of bcbio
<https://bcbio-nextgen.readthedocs.io/en/latest/contents/installation.html>`_.

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

- `Cromwell <http://cromwell.readthedocs.io>`_ -- multicore local runs and
  distributed runs on HPC systems with shared filesystems and schedulers like
  SLURM, SGE and PBSPro.

- `Arvados <https://arvados.org/>`_ -- a hosted platform that runs on top of
  parallel cloud environments. We include an example below of running on the
  `public Curoverse <https://cloud.curoverse.com/>`_ instance running on
  `Microsoft Azure <https://azure.microsoft.com>`_.

- `DNANexus <https://www.dnanexus.com/>`_ -- a hosted platform running
  distributed jobs on cloud environments, working with both AWS and Azure.

- `Seven Bridges <https://www.sevenbridges.com/>`_ -- parallel distributed
  analyses on the Seven Bridges platform and `Cancer Genomics Cloud
  <http://www.cancergenomicscloud.org/>`_.

- `Toil <https://github.com/BD2KGenomics/toil>`_ -- parallel local and
  distributed cluster runs on schedulers like SLURM, SGE and PBSPro.

- `rabix bunny <https://github.com/rabix/bunny>`_ -- multicore local runs.

- `cwltool <https://github.com/common-workflow-language/cwltool>`_ -- a single
  core analysis engine, primarily used for testing.

We plan to continue to expand CWL support to include more components of bcbio,
and also need to evaluate the workflow on larger, real life analyses. This
includes supporting additional CWL runners. We're working on evaluating
`Galaxy/Planemo <https://github.com/galaxyproject/planemo>`_ for integration
with the Galaxy community.

.. _docs-cwl-installation:

Installation
~~~~~~~~~~~~

`bcbio-vm <https://github.com/bcbio/bcbio-nextgen-vm>`_ installs all
dependencies required to generate CWL and run bcbio, along with supported CWL
runners. There are two install choices, depending on your usage of bcbio:
running CWL with a existing local bcbio install, or running with containers.

Install bcbio-vm with a local bcbio
===================================

To run bcbio without using containers, first `install bcbio
<https://bcbio-nextgen.readthedocs.io/en/latest/contents/installation.html#automated>`_
and make it available in your path. You'll need both the bcbio code and tools.
To only run the tests and bcbio validations, you don't need a full data
installation so can install with ``--nodata``.

To then install bcbio-vm, add the ``--cwl`` flag to the install::

    bcbio_nextgen.py upgrade --cwl

Adding this to any future upgrades will also update the bcbio-vm wrapper code
and tools.

When you begin running your own analysis and need the data available,
pre-prepare your bcbio data directory with ``bcbio_nextgen.py upgrade --data
--cwl``.

Install bcbio-vm with containers
================================

If you don't have an existing local bcbio installation and want to run with CWL
using the tools and data embedded in containers, you can do a stand along
install of just bcbio-vm. To install using `Miniconda
<http://conda.pydata.org/miniconda.html>`_ and `bioconda packages
<https://bioconda.github.io/>`_ on Linux::

    export TARGETDIR=~/install/bcbio-vm/anaconda
    export BINDIR=/usr/local/bin
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p $TARGETDIR
    $TARGETDIR/bin/conda install --yes -c conda-forge -c bioconda python=3 bcbio-nextgen
    $TARGETDIR/bin/conda install --yes -c conda-forge -c bioconda python=3 bcbio-nextgen-vm
    mkdir -p $BINDIR
    ln -s $TARGETDIR/bin/bcbio_vm.py $BINDIR/bcbio_vm.py
    ln -s $TARGETDIR/bin/conda $BINDIR/bcbiovm_conda
    ln -s $TARGETDIR/bin/python $BINDIR/bcbiovm_python

In the above commands, the `bcbio-vm` install goes in ``$TARGETDIR``.
The example is in your home directory but set it anywhere you have space.
Also, as an alternative to symbolic linking to a ``$BINDIR``, you can
add the install bin directory to your PATH::

    export PATH=$TARGETDIR/bin:$PATH

This install includes bcbio-nextgen libraries, used in generating CWL and
orchestrating runs, but is not a full bcbio installation. It requires
`Docker <https://www.docker.com/>`_ present on your
system this is all you need to get started running examples, since the CWL
runners will pull in Docker containers with the bcbio tools.

Getting started
~~~~~~~~~~~~~~~

To make it easy to get started, we have pre-built CWL descriptions that
use test data. These run in under 5 minutes on a local machine and
don't require a bcbio installation if you have Docker available on
your machine:

1. Download and unpack the `test repository <https://github.com/bcbio/test_bcbio_cwl>`_::

     wget -O test_bcbio_cwl.tar.gz https://github.com/bcbio/test_bcbio_cwl/archive/master.tar.gz
     tar -xzvpf test_bcbio_cwl.tar.gz
     cd test_bcbio_cwl-master/somatic

2. Run the analysis using either Cromwell, Rabix bunny or Toil. If you have Docker
   available on your machine, the runner will download the correct `bcbio
   container <https://github.com/bcbio/bcbio_docker>`_ and you don't need to
   install anything else to get started. If you have an old version of the
   container you want to update to the latest with ``docker pull
   quay.io/bcbio/bcbio-vc``. There are shell scripts that provide the command
   lines for running::

     bash run_cromwell.sh
     bash run_bunny.sh
     bash run_toil.sh

   Or you can run directly using the ``bcbio_vm.py`` wrappers::

     bcbio_vm.py cwlrun cromwell somatic-workflow
     bcbio_vm.py cwlrun toil somatic-workflow
     bcbio_vm.py cwlrun bunny somatic-workflow

   Thes wrappers automatically handle temporary directories, permissions,
   logging and re-starts. If running without Docker, use a `local installation of
   bcbio
   <https://bcbio-nextgen.readthedocs.org/en/latest/contents/installation.html>`_
   add ``--no-container`` to the commands in the shell scripts.

.. _docs-cwl-generate:

Generating CWL for input to a tool
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first step in running your analysis project in bcbio is to generate CWL. If
you're already familiar with bcbio, the `process of preparing information about
your sample inputs and analysis <bcbio-nextgen.readthedocs.io/en/latest/contents/configuration.html>`_
are almost identical:

- A `standard bcbio sample configuration file
  <https://bcbio-nextgen.readthedocs.io/en/latest/contents/configuration.html>`_
  defining the samples. This can either be a full prepared YAML file or a
  `template file and CSV with sample data <http://bcbio-nextgen.readthedocs.io/en/latest/contents/configuration.html#automated-sample-configuration>`_.

- A ``bcbio_system.yaml`` file defining the system environment for running the
  program. This includes the resource specification with `cores and memory per
  core for your machines
  <http://bcbio-nextgen.readthedocs.io/en/latest/contents/configuration.html#resources>`_.
  For choosing cores and memory per cores, you generally want to set this to
  match the parameters of a single machine either for a local run or on a
  cluster.

  In addition to `resources
  <http://bcbio-nextgen.readthedocs.io/en/latest/contents/configuration.html#resources>`_
  specifications, the bcbio system file now also includes paths to the
  reference biodata and optionally input file directories if you want to avoid
  specifying full paths to your inputs in the ``bcbio_vm.py template`` command.
  bcbio will recursively look up file locations within those ``inputs``, and
  this has the advantage of working identically for non-local file locations.
  Here is an example for a 16 core machine with 3.5Gb of memory per core::

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

    bcbio_vm.py template --systemconfig bcbio_system.yaml template.yaml samples.csv [optional list of fastq or BAM inputs]
    bcbio_vm.py cwl --systemconfig bcbio_system.yaml samples/config/samples.yaml

producing a ``sample-workflow`` output directory with the CWL.


On a first CWL generation run with a new genome, this process will run for a
longer time as it needs to make your reference compatible with CWL. This
includes creating single tar.gz files from some reference directories so they
can get passed to CWL steps where they'll get unpacked. This process only
happens a single time and keeps unpacked versions so your reference setup is
compatible with both old bcbio IPython and new CWL runs.

You can now run this with any CWL compatible runner and the ``bcbio_vm.py
cwlrun`` wrappers standardize running across multiple tools in different
environments.

Running with Cromwell (local, HPC, cloud)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The `Cromwell <http://cromwell.readthedocs.io/>`_ workflow management system runs
bcbio either locally on a single machine or distributed on a cluster using a
scheduler like SLURM, SGE or PBSPro.

To run a bcbio CWL workflow locally using Docker::

    bcbio_vm.py cwlrun cromwell sample-workflow

If you want to run from a locally installed bcbio add ``--no-container`` to the
commandline.

To run distributed on a SLURM cluster::

    bcbio_vm.py cwlrun cromwell sample-workflow --no-container -q your_queue -s slurm -r timelimit=0-12:00

Tweak scheduler parameters using the
`same options as the older bcbio IPython approach <http://bcbio-nextgen.readthedocs.io/en/latest/contents/parallel.html#ipython-parallel>`_.

To control the resources used Cromwell, set `--joblimit` to the allowed jobs
allocated concurrently. This isn't total cores used, but rather the number of jobs
either locally or remotely scheduled concurrently. Since CWL steps are
heterogeneous and use only cores necessary for that job, the total cores used
will max out at joblimit times maximum cores for an individual process. Setting
this helps avoid over-committing jobs to a shared scheduler during highly
parallel processes like variant calling.

Cromwell can also run directly on cloud resources: :ref:`docs-cloud-gcp`.

Running with Toil (local, HPC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The `Toil pipeline management system <https://github.com/BD2KGenomics/toil>`_
runs CWL workflows in parallel on a local machine, on a cluster or at AWS.

To run a bcbio CWL workflow locally with Toil using Docker::

    bcbio_vm.py cwlrun toil sample-workflow

If you want to run from a locally installed bcbio add ``--no-container`` to the
commandline.

To run distributed on a Slurm cluster::

    bcbio_vm.py cwlrun toil sample-workflow -- --batchSystem slurm

Running on Arvados (hosted cloud)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bcbio generated CWL workflows run on `Arvados <https://arvados.org/>`_ and these
instructions detail how to run on the `Arvdos public instance
<https://cloud.curoverse.com/>`_. `Arvados cwl-runner
<https://github.com/curoverse/arvados>`_ comes pre-installed with `bcbio-vm
<https://github.com/bcbio/bcbio-nextgen-vm#installation>`_.
We have a publicly accessible project, called `bcbio_resources
<https://workbench.qr1hi.arvadosapi.com/projects/qr1hi-j7d0g-8g1u4lh8mwev36n>`_
that contains the latest Docker images, test data and genome references you can
use for runs.

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
   UUID. You can also use the existing genomes pre-installed in the
   ``bcbio_resources`` project if using the public Arvados playground::

     arv-put --name testdata_genomes --project-uuid $PROJECT_ID testdata/genomes/hg19

3. Upload input data to Arvados Keep. Note the collection UUID::

     arv-put --name testdata_inputs --project-uuid $PROJECT_ID testdata/100326_FC6107FAAXX testdata/automated testdata/reference_material

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

6. In most cases, Arvados should directly pick up the Docker images you need
   from the public bcbio_resources project in your instance. If you need to
   manually add to your project, you can copy latest bcbio Docker image into
   your project from bcbio_resources using `arv-copy
   <https://doc.arvados.org/user/topics/arv-copy.html>`_. You'll need to find
   the UUID of ``quay.io/bcbio/bcbio-vc`` and ``arvados/jobs``::

     arv-copy $JOBS_ID --project-uuid $PROJECT_ID --src qr1hi --dst qr1hi
     arv-copy $BCBIO_VC_ID --project-uuid $PROJECT_ID --src qr1hi --dst qr1hi

   or import local Docker images to your Arvados project::

     docker pull arvados/jobs:1.0.20180216164101
     arv-keepdocker --project $PROJECT_ID -- arvados/jobs 1.0.20180216164101
     docker pull quay.io/bcbio/bcbio-vc
     arv-keepdocker --project $PROJECT_ID -- quay.io/bcbio/bcbio-vc latest

7. Run the CWL on the Arvados public cloud using the Arvados cwl-runner::

     bcbio_vm.py cwlrun arvados arvados_testcwl-workflow -- --project-uuid $PROJECT_ID

Running on DNAnexus (hosted cloud)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bcbio runs on the `DNAnexus platform <https://www.dnanexus.com/>`_ by converting
bcbio generated CWL into DNAnexus workflows and apps using
`dx-cwl <https://github.com/dnanexus/dx-cwl>`_. This describes the process
using the bcbio workflow app (bcbio-run-workflow) and
`bcbio workflow applet (bcbio_resources:/applets/bcbio-run-workflow) <https://platform.dnanexus.com/projects/F541fX00f5v9vKJjJ34gvgbv/data/applets>`_
in the public `bcbio_resources
<https://platform.dnanexus.com/projects/F541fX00f5v9vKJjJ34gvgbv/data/>`_
project, Both are `regularly updated and maintained on the DNAnexus
platform <https://github.com/bcbio/bcbio-dnanexus-wrapper>`_. Secondarily, we
also show how to install and create workflows locally for
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

3. Create a bcbio system YAML file with projects, locations of files and
   desired core and memory usage for jobs. bcbio uses the core and memory
   specifications to determine machine instance types to use::

     dnanexus:
       project: PNAME
       ref:
         project: bcbio_resources
         folder: /reference_genomes
       inputs:
         - /data/input
         - /data/input/regions
     resources:
       default: {cores: 8, memory: 3000M, jvm_opts: [-Xms1g, -Xmx3000m]}

4. Create a bcbio sample CSV file referencing samples to run. The files can be
   relative to the ``inputs`` directory specified above; bcbio will search
   recursively for files, so you don't need to specify full paths if your file
   names are unique. Start with a sample specification::

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
   ``template.yaml``) to the either the bcbio-run-workflow app or applet. This
   example uses a specific version of the bcbio app for full reproducibility;
   any future re-runs will always use the exact same versioned tools and
   workflows. You can do this using the web interface or via the command line with a small
   script like::

      TEMPLATE=germline
      APP_VERSION=0.0.2
      FOLDER=/bcbio/$PNAME
      dx select "$PROJECT"
      dx mkdir -p $FOLDER
      for F in $TEMPLATE-template.yaml $PNAME.csv bcbio_system-dnanexus.yaml
      do
              dx rm -a /$FOLDER/$F || true
              dx upload --path /$FOLDER/ $F
      done
      dx ls $FOLDER
      dx rm -a -r /$FOLDER/dx-cwl-run || true
      dx run bcbio-run-workflow/$APP_VERSION -iyaml_template=/$FOLDER/$TEMPLATE-template.yaml -isample_spec=/$FOLDER/$PNAME.csv -isystem_configuration=/$FOLDER/bcbio_system-dnanexus.yaml -ioutput_folder=/$FOLDER/dx-cwl-run

   Alternatively if you want the latest bcbio code, change the final command to
   use the applet. Everything else in the script is identical::

       dx run bcbio_resources:/applets/bcbio-run-workflow -iyaml_template=/$FOLDER/$TEMPLATE-template.yaml -isample_spec=/$FOLDER/$PNAME.csv -isystem_configuration=/$FOLDER/bcbio_system-dnanexus.yaml -ioutput_folder=/$FOLDER/dx-cwl-run

The app will lookup all files, prepare a bcbio CWL workflow, convert into a
DNAnexus workflow, and submit to the platform. The workflow runs as a standard
DNAnexus workflow and you can monitor through the command line (with ``dx find
executions --root job-YOURJOBID`` and ``dx watch``) or the web interface
(``Monitor`` tab).

If you prefer not to use the DNAnexus app, you can also submit jobs locally by
installing `bcbio-vm <https://github.com/bcbio/bcbio-nextgen-vm#installation>`_
on your local machine. This can also be useful to test generation of CWL and
manually ensure identification of all your samples and associated files on the
DNAnexus platform.

1. Follow the :ref:`automated-sample-config` workflow to generate a full
   configuration, and generate a CWL description of the workflow::

       TEMPLATE=germline
       rm -rf $PNAME $PNAME-workflow
       bcbio_vm.py template --systemconfig bcbio_system-dnanexus.yaml $TEMPLATE-template.yaml $PNAME.csv
       bcbio_vm.py cwl --systemconfig bcbio_system-dnanexus.yaml $PNAME/config/$PNAME.yaml

2. Determine project information and login credentials. You'll want to note the
   ``Auth token used`` and ``Current workspace`` project ID::

       dx env

3. Compile the CWL workflow into a DNAnexus workflow::

       dx-cwl compile-workflow $PNAME-workflow/main-$PNAME.cwl \
          --project PROJECT_ID --token $DX_AUTH_TOKEN \
          --rootdir $FOLDER/dx-cwl-run

4. Upload sample information from generated CWL and run workflow::

       FOLDER=/bcbio/$PNAME
       dx mkdir -p $DX_PROJECT_ID:$FOLDER/$PNAME-workflow
       dx upload -p --path $DX_PROJECT_ID:$FOLDER/$PNAME-workflow $PNAME-workflow/main-$PNAME-samples.json
       dx-cwl run-workflow $FOLDER/dx-cwl-run/main-$PNAME/main-$PNAME \
              $FOLDER/$PNAME-workflow/main-$PNAME-samples.json \
              --project PROJECT_ID --token $DX_AUTH_TOKEN \
              --rootdir $FOLDER/dx-cwl-run

Running on Seven Bridges (hosted cloud)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bcbio runs on the `Seven Bridges <https://www.sevenbridges.com/>`_
including the main platform and specialized data sources like the
`Cancer Genomics Cloud <https://www.cancergenomicscloud.org/>`_ and
`Cavatica <http://www.cavatica.org/>`_. Seven Bridges uses generated CWL
directly and bcbio has utilities to query your remote data on the platform and
prepare CWL for direct submission.

1. Since Seven Bridges is available on multiple platforms and data access
   points, we authenticate with a configuration file in
   ``$HOME/.sevenbridges/credentials`` with potentially `multiple profiles defining
   API access URLs and authentication keys
   <https://sevenbridges-python.readthedocs.io/en/latest/quickstart/#initialize-the-library-using-a-configuration-file>`_.
   We reference the `specified credentials
   <https://docs.sevenbridges.com/docs/store-credentials-to-access-seven-bridges-client-applications-and-libraries#section-unified-configuration-file>`_
   when setting up a ``bcbio_system-sbg.yaml`` file to ensure correct authentication.

2. Upload your inputs and bcbio reference data using the `Seven Bridges command
   line uploader
   <https://docs.sevenbridges.com/docs/upload-via-the-command-line>`_. We plan
   to host standard bcbio reference data in a public project so you should only
   need to upload your project specific data::

       sbg-uploader.sh -p chapmanb/bcbio-test --folder inputs --preserve-folder fastq_files regions

3. Create ``bcbio_system-sbg.yaml`` file defining locations of inputs::

       sbgenomics:
         profile: default
         project: chapmanb/bcbio-test
         inputs:
           - /testdata/100326_FC6107FAAXX
           - /testdata/automated
           - /testdata/genomes
           - /testdata/reference_material
       resources:
         default:
           cores: 2
           memory: 3G
           jvm_opts: [-Xms750m, -Xmx3000m]

4. Follow the :ref:`automated-sample-config` workflow to generate a full
   configuration, and generate a CWL description of the workflow::

       PNAME=somatic
       bcbio_vm.py template --systemconfig=bcbio_system-sbg.yaml ${PNAME}_template.yaml $PNAME.csv
       bcbio_vm.py cwl --systemconfig=bcbio_system-sbg.yaml $PNAME/config/$PNAME.yaml

5. Run the job on the Seven Bridges platform::

       PNAME=somatic
       SBG_PROJECT=bcbio-test
       bcbio_vm.py cwlrun sbg ${PNAME}-workflow -- --project ${SBG_PROJECT}

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
   ensemble calling, heterogeneity detection and disambiguation.

-  Port RNA-seq and small RNA workflows to CWL.
