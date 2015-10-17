## bcbio with the Common Workflow Language (CWL)

This directory contains exploratory work to run bcbio using the
[Common Workflow Language (CWL)][0]. The initial implementation parses a bcbio
system YAML file and run info YAML file to produce a CWL workflow and run that
handles multiple samples.

To prepare your system for running this you first need to install cwltool
from the [common workflow language reference implementation][1]:
```
git clone https://github.com/common-workflow-language/cwltool.git
cd cwltool
python setup.py install
```
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

Then run locally using the tests:
```
cd tests && ./run_tests.sh cwl
```
or manually with:
```
cwltool --verbose --preserve-env HOME PATH --no-container bcbio2cwl.cwl testinput-args.json
```
or
```
cwltool --verbose --no-container bcbio2cwl.cwl \
  --run_config ../tests/data/automated/run_info-bam.yaml \
  --system_config testinput-bcbio_system.yaml
```
The runner should report success and you'll see a `world.yaml` output file that
contains the parsed and prepped reference data along with the input sample
information.`

Docker support is a work in progress and not expected to work right now. To test
and develop, run:
```
cwltool --verbose bcbio2cwl.cwl testinput-args.json
```
[0]: https://github.com/common-workflow-language/common-workflow-language
[1]: https://github.com/common-workflow-language/cwltool
[2]: https://github.com/chapmanb/bcbio-nextgen
[3]: https://bcbio-nextgen.readthedocs.org/en/latest/contents/installation.html
