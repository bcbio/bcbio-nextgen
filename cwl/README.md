## bcbio with the Common Workflow Language (CWL)

This directory contains details on work to run bcbio using tools that support
the [Common Workflow Language (CWL)][0]. bcbio generates a CWL workflow from a
[sample YAML description file](https://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html).
Any tool that supports CWL input can run this workflow. CWL-based tools do the
work of managing files and workflows, and bcbio performs the biological
analysis using a Docker container or local install.

This is a work in progress and not yet a complete production implementation. The
documentation orients anyone interested in helping with development.

### Current status

bcbio currently supports creation of CWL for alignment, small variant calls
(SNPs and indels), coverage assessment, HLA typing and quality control. It
generates a [CWL draft 3](http://common-workflow-language.github.io/draft-3/)
compatible workflow. The actual biological code execution during runs works with
either the bcbio docker container
([bcbio/bcbio](https://hub.docker.com/r/bcbio/bcbio/)) or a local installation.

The implementation includes bcbio's approaches to parallelization and batching.
At the top level workflow, we parallelize by samples. Using sub-workflows, we
split samples into sections for processing, which enables parallel alignment
over multiple machines following by merging. We also use sub-workflows, along
with CWL records, to batch multiple samples and run in parallel. This enables
pooled and tumor/normal cancer calling with parallelization by chromosome
regions based on coverage calculations.

<img src="http://i.imgur.com/iyU8VIZ.png" width="600" height="700" />

We plan to continue to expand CWL support to include more components of bcbio,
and also need to evaluate the workflow on larger, real life analyses. This includes:

- Running on CWL runners like [Arvados](https://arvados.org/),
  [Galaxy/Planemo](https://github.com/galaxyproject/planemo) and
  [toil](https://github.com/BD2KGenomics/toil). Our current evaluations use
  [cwltool][1].

- Enable checkpointing of analysis so we can re-start a failed workflow at the
  previous end point.

- Test parallelization on single machines and over multiple machines with
  distributed CWL runners.

### Running

To prepare your system for running you first need to install cwltool
from the [common workflow language reference implementation][1]. The easiest way
is to use conda. You can install [Miniconda](http://conda.pydata.org/miniconda.html)
if you don't have conda present, and cwltool is also included with
[bcbio-vm](https://github.com/chapmanb/bcbio-nextgen-vm) installs:

```
conda install -c bioconda cwltool
```
cwltool uses javascript for data manipulation and requires either a local
installation of [nodejs](https://nodejs.org) or having
[Docker](https://www.docker.com/) installed and running `docker pull node:slim`.

To make it easy to get started, we have a pre-built CWL description that uses
test data. This will run in under 5 minutes on a local machine and doesn't require
a bcbio installation if you have Docker available on your machine:

1. Download and unpack the test data:
   ```
   wget https://s3.amazonaws.com/bcbio/cwl/test_bcbio_cwl.tar.gz
   tar -xzvpf test_bcbio_cwl.tar.gz
   cd test_bcbio_cwl
   ```

2. Run the analysis using `cwltool`. If you have Docker, cwltool will download the
   `bcbio/bcbio` container and you don't need to install anything else to get
   started. You can use the `run_cwl.sh` script or run directly from the command
   line:
   ```
   cwltool --verbose run_info-cwl-workflow/main-run_info-cwl.cwl run_info-cwl-workflow/main-run_info-cwl-samples.json
   ```
   If you don't have Docker, you can also use a
   [local installation of bcbio][3]. You don't need to install genome data since
   the tests use small local data. Then run with:
   ```
   cwltool --verbose --preserve-environment PATH HOME --no-container \
     run_info-cwl-workflow/main-run_info-cwl.cwl \
     run_info-cwl-workflow/main-run_info-cwl-samples.json
   ```

### Generating CWL from bcbio

To generate CWL for a different bcbio sample YAML file, you need to install
[bcbio-vm](https://github.com/chapmanb/bcbio-nextgen-vm), which
generates the CWL command line wrapper. You can install with `conda install -c
bioconda bcbio-nextgen-vm`.

As an example, to generate the test data show above, clone the
[bcbio GitHub repository locally][2] to get the test suite and run a minimal CWL
workflow generated automatically by bcbio from the inputs:
```
git clone https://github.com/chapmanb/bcbio-nextgen.git
cd bcbio-nextgen/tests
./run_tests.sh cwl_local
./run_tests.sh cwl_docker
```
This will create a CWL workflow inside `test_automated_output` which you can run
again manually with either a local bcbio installation:

To generate CWL directly from a sample input and the test bcbio system file:
```
bcbio_vm.py cwl ../data/automated/run_info-cwl.yaml --systemconfig bcbio_system.yaml
```

### Running bcbio CWL on Arvados

We're actively testing bcbio generated CWL workflows on
[Arvados](https://arvados.org/). These instructions detail how to run this on
the [Arvdos public instance](https://cloud.curoverse.com/).

[Arvados cwl-runner](https://github.com/curoverse/arvados) comes pre-installed
with bcbio-vm. If you don't have it you can install using conda:
```
conda install -c bioconda arvados-cwl-runner
```

Retrieve API keys from the
[Arvados public instance](https://cloud.curoverse.com/). Login, then go to ['User
Icon -> Personal Token'](https://cloud.curoverse.com/current_token). Copy and
paste the commands given there into your shell. You'll specifically need to set
`ARVADOS_API_HOST` and `ARVADOS_API_TOKEN`.

To run an analysis:

1. Create a new project from the web interface (Projects -> Add a new project).
   Note the project ID from the URL of the project (an identifier like
   `qr1hi-j7d0g-7t73h4hrau3l063`).

2. Upload reference data to Aravdos Keep. Note the genome collection
portable data hash:
```
arv-put --portable-data-hash --name hg19-testdata --project-uuid qr1hi-j7d0g-7t73h4hrau3l063 testdata/genomes
```

3. Upload input data to Arvados Keep. Note the collection portable data hash:
```
arv-put --portable-data-hash --name input-testdata --project-uuid qr1hi-j7d0g-7t73h4hrau3l063 testdata/100326_FC6107FAAXX testdata/automated testdata/reference_material
```

4. Create an Arvados section in a `bcbio_system.yaml` file specifying locations
to look. The token and host sections are optional and can instead be in
environmental variables:
```
arvados:
  token: your_token
  host: qr1hi.arvadosapi.com
  project: qr1hi-j7d0g-7t73h4hrau3l063
  reference: a84e575534ef1aa756edf1bfb4cad8ae+1927
  input: a1d976bc7bcba2b523713fa67695d715+464
resources:
  default:
    cores: 4
    memory: 1G
  bwa:
    cores: 4
    memory: 2G
  gatk:
    jvm_opts: [-Xms750m, -Xmx2500m]
```

5. Generate the CWL to run your samples. If you're using multiple input files with a
   [CSV metadata file and template](https://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html#automated-sample-configuration
   then start with creation of a configuration file:
```
bcbio_vm.py template --systemconfig bcbio_system_arvados.yaml testcwl_template.yaml testcwl.csv
```
To generate the CWL from the system and sample configuration files:
```
bcbio_vm.py cwl --systemconfig bcbio_system_arvados.yaml testcwl/config/testcwl.yaml

```

6. Run the CWL on the Arvados public cloud using the Arvados cwl-runner:
```
arvados-cwl-runner --project-uuid qr1hi-your-projectuuid --enable-reuse testcwl-workflow/main-testcwl.cwl testcwl-workflow/main-testcwl-samples.json
```

### Development notes

bcbio generates a common workflow language description. Internally, bcbio
represents the files and information related to processing as
[a comprehensive dictionary](https://bcbio-nextgen.readthedocs.org/en/latest/contents/code.html#data).
This world object describes the state of a run and associated files, and new
processing steps update or add information to it. The world object is roughly
equivalent to CWL's JSON-based input object, but CWL enforces additional
annotations to identify files and models new inputs/outputs at each step. The
work in bcbio is to move from our laissez-faire approach to the more structured
CWL model.

The generated CWL workflow is in `run_info-cwl-workflow`:

- `main-*.cwl` -- the top level CWL file describing the workflow steps
- `main*-samples.json` -- the flattened bcbio world structure represented as
  CWL inputs
- `wf-*.cwl` -- CWL sub-workflows, describing sample level parallel processing
  of a section of the workflow, with potential internal parallelization.
- `steps/*.cwl` -- CWL descriptions of sections of code run inside bcbio. Each
  of these are potential parallelization points and make up the nodes in the
  workflow.

To help with defining the outputs at each step, there is a `WorldWatcher` object
that can output changed files and world dictionary objects between steps in the
pipeline when running a bcbio in the standard way. The
[variant pipeline](https://github.com/chapmanb/bcbio-nextgen/blob/master/bcbio/pipeline/main.py)
has examples using it. This is useful when preparing the CWL
definitions of inputs and outputs for new steps in the [bcbio CWL step definitions](https://github.com/chapmanb/bcbio-nextgen/blob/master/bcbio/cwl/workflow.py).

### ToDo

These are additional items of emphasis beyond the development work listed under
`Current Status`:

- Support the full variant calling workflow with additional steps like
  ensemble calling, structural variation, heterogeneity detection and
  disambiguation.

- Port RNA-seq and small RNA workflows to CWL.

- Determine when we should skip steps based on configuration to avoid writing
  them to the CWL file. For instance, right now we include HLA typing even if
  it's not defined and have an extra do-nothing step in the CWL output. We
  should have a clean way to skip writing this step if not needed based on the
  configuration.

- Replace the custom python code in the
  [bcbio step definitions](https://github.com/chapmanb/bcbio-nextgen/blob/master/bcbio/cwl/workflow.py)
  with a higher level DSL in YAML we can parse and translate to CWL.

[0]: https://github.com/common-workflow-language/common-workflow-language
[1]: https://github.com/common-workflow-language/cwltool
[2]: https://github.com/chapmanb/bcbio-nextgen
[3]: https://bcbio-nextgen.readthedocs.org/en/latest/contents/installation.html
