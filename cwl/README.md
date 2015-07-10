## bcbio with the Common Workflow Language (CWL)

This directory contains exploratory work to run bcbio using the
[Common Workflow Language (CWL)][0]. The initial implementation parses a bcbio
system YAML file and run info YAML file to produce a CWL workflow and run that
handles multiple samples.

To prepare your system for running this you first need to install cwl-runner
from the [common workflow language reference implementation][1].

Then you either need a working [local installation of bcbio][3] or you can run
off a minimal installation of bcbio python code and test data without the third
party dependencies and reference data:

- Clone the [bcbio GitHub repository locally][2].
- Install the Python code so `bcbio_nextgen.py` is available in your PATH.
  Either use `python setup.py install` or Anaconda:
  `conda install -c bcbio bcbio-nextgen`
- Run at least one test to download the testing reference genomes. The test will
  not work unless you have all of bcbio installed, but for the initial
  CWL testing you shouldn't need third party tools so this is okay:
  `cd tests && ./run_tests.sh devel`

Then run locally with:
```
cwl-runner --verbose --no-container bcbio2cwl.cwl testinput-args.json
```
The runner should report success and you'll see a `world.yaml` output file that
contains the parsed and prepped reference data along with the input sample
information.`

Docker support is a work in progress and not expected to work right now. To test
and develop, run:
```
cwl-runner --verbose bcbio2cwl.cwl testinput-args.json
```
[0]: https://github.com/common-workflow-language/common-workflow-language
[1]: https://github.com/common-workflow-language/common-workflow-language/tree/master/reference
[2]: https://github.com/chapmanb/bcbio-nextgen
[3]: https://bcbio-nextgen.readthedocs.org/en/latest/contents/installation.html
