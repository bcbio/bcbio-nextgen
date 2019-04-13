#!/usr/bin/env python
"""Setup file and install script for NextGen sequencing analysis scripts.
"""
import os
from setuptools import setup, find_packages

version = "1.1.6a"

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

install_requires = [] # install dependencies via conda
zip_safe = False
scripts = ['scripts/bcbio_nextgen.py', 'scripts/bcbio_setup_genome.py', 'scripts/bcbio_prepare_samples.py',
           'scripts/bcbio_fastq_umi_prep.py', 'scripts/cwltool2wdl.py']

write_version_py()
setup(name="bcbio-nextgen",
      version=version,
      author="bcbio community",
      author_email="biovalidation@googlegroups.com",
      description="Best-practice pipelines for fully automated high throughput sequencing analysis",
      long_description=(open('README.rst').read()),
      license="MIT",
      url="https://github.com/bcbio/bcbio-nextgen",
      packages=find_packages(exclude=["tests"]),
      zip_safe=zip_safe,
      scripts=scripts,
      install_requires=install_requires,
      include_package_data=True)
