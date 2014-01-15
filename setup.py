#!/usr/bin/env python
"""Setup file and install script for NextGen sequencing analysis scripts.
"""
import os
from setuptools import setup, find_packages

version = "0.7.6"

def write_version_py():
    version_py = os.path.join(os.path.dirname(__file__), 'bcbio', 'pipeline',
                              'version.py')
    try:
        import subprocess
        p = subprocess.Popen(["git", "rev-parse", "--short", "HEAD"],
                             stdout=subprocess.PIPE)
        githash = p.stdout.read().strip()
    except:
        githash = ""
    with open(version_py, "w") as out_handle:
        out_handle.write("\n".join(['__version__ = "%s"' % version,
                                    '__git_revision__ = "%s"' % githash]))

with open("requirements.txt", "r") as f:
    install_requires = [x.strip() for x in f.readlines() if not x.startswith("bcbio-nextgen")]

write_version_py()
setup(name="bcbio-nextgen",
      version=version,
      author="Brad Chapman and bcbio-nextgen contributors",
      author_email="chapmanb@50mail.com",
      description="Best-practice pipelines for fully automated high throughput sequencing analysis",
      license="MIT",
      url="https://github.com/chapmanb/bcbio-nextgen",
      namespace_packages=["bcbio"],
      packages=find_packages(),
      zip_safe=False,
      scripts=['scripts/bcbio_nextgen.py',
               'scripts/bam_to_wiggle.py',
               'scripts/barcode_sort_trim.py',
               'scripts/illumina_finished_msg.py',
               'scripts/nextgen_analysis_server.py',
               'scripts/solexa_qseq_to_fastq.py'],
      install_requires=install_requires)
