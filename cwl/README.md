## bcbio with the Common Workflow Language (CWL)

This directory contains exploratory work to run bcbio using the
[Common Workflow Language (CWL)][0]. In the initial implementation, bcbio
generates a CWL workflow and starts a variant calling processing run. This is
still a work in progress and not close to a complete implementation, and this
documentation is to orient anyone interested in helping with development.

### Running

To prepare your system for running this you first need to install cwltool
from the [common workflow language reference implementation][1]:
```
git clone https://github.com/common-workflow-language/cwltool.git
cd cwltool
python setup.py install
```
Then you need a working installation of
[bcbio-vm](https://github.com/chapmanb/bcbio-nextgen-vm), which provides the CWL
command line wrapper. The actual code runs can be from either bcbio_vm.py, with
the Docker container running tools, or a [local installation of bcbio][3].
You don't need to install genome data since the tests use small local data.

Next clone the [bcbio GitHub repository locally][2] to get the test suite:
```
git clone https://github.com/chapmanb/bcbio-nextgen.git
cd bcbio-nextgen/tests
./run_tests.sh cwl_local
```
This will create a CWL workflow inside `test_automated_output` which you can run
with either a local bcbio installation:
```
cwltool --verbose --preserve-environment PATH HOME --no-container \
  run_info-bam-workflow/run_info-bam-main.cwl \
  run_info-bam-workflow/run_info-bam-main-samples.json
```
or with bcbio inside a Docker container:

```
cwltool --verbose run_info-bam-workflow/run_info-bam-main.cwl run_info-bam-workflow/run_info-bam-main-samples.json
```

### Development notes

bcbio does the work of auto-generating a common workflow language description.
Internally, bcbio represents the files and information related to processing as
[a dictionary](https://bcbio-nextgen.readthedocs.org/en/latest/contents/code.html#data).
This world object describes the state of a run and associated files, and new
processing steps update or add information to it. The world object is roughly
equivalent to CWL's JSON-based input object, but CWL enforces additional
annotations to identify files and models new inputs/outputs at each step. The
work in bcbio is to move from our laissez-faire approach to the more structured
CWL model.

The generated `run_info-bam-workflow/run_info-bam-main-samples.json` file
contains the flattened bcbio structure represented as CWL.

Current items to do:

- Handle parallelization of multiple samples by using arrays of values and a
  [dotproduct scatter](http://common-workflow-language.github.io/draft-2/#workflowstep).
  Example
  [JSON input](https://github.com/common-workflow-language/common-workflow-language/blob/master/draft-2/draft-2/scatter-job2.json)
  and
  [CWL description](https://github.com/common-workflow-language/common-workflow-language/blob/master/draft-2/draft-2/scatter-wf4.cwl).

- Remove requirement and assumptions for a shared filesystem/work directory in
  bcbio (`dirs` object) or model with the local CWL temporary directory.
  Specifically do this for the first step `prep_align_inputs` so we can get a
  single step fully running.

- Create new docker container and test with the new CWL structure/approach.

- Explicitly define inputs/outputs at each bcbio step to avoid needing to pass
  around and return the entire world object. One approach would be to write
  something that checks the state of the world object and shared filesystem
  after each step in the standard test run, and writes an annotated file of
  changes. We could turn this into an curated list of changes associated with
  each step that bcbio uses to prepare the CWL.

[0]: https://github.com/common-workflow-language/common-workflow-language
[1]: https://github.com/common-workflow-language/cwltool
[2]: https://github.com/chapmanb/bcbio-nextgen
[3]: https://bcbio-nextgen.readthedocs.org/en/latest/contents/installation.html
