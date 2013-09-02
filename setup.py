#!/usr/bin/env python
"""Setup file and install script for NextGen sequencing analysis scripts.
"""
import os
from setuptools import setup, find_packages

version_py = os.path.join(os.path.dirname(__file__), 'bcbio', 'pipeline',
                          'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"', '')

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
                 'scripts/solexa_qseq_to_fastq.py',
                 ],
      install_requires=[
          "bioblend >= 0.3.3",
          "biopython >= 1.61",
          "boto >= 2.9.6",
          "cutadapt >= 1.2.1",
          "Cython >= 0.19",
          "fabric >= 1.7.0",
          "ipython >= 1.0",
          "ipython-cluster-helper >= 0.1.12",
          "Logbook >= 0.4.1",
          "Mako >= 0.7.3",
          "psutil >= 0.6.1",
          "pybedtools >= 0.6.2",
          "py_descriptive_statistics >= 0.2",
          "pysam >= 0.6",
          "PyYAML >= 3.10",
          "pyzmq >= 2.1.11",
          "joblib >= 0.7.0d",
          "sh >= 1.07",
          "HTSeq >= 0.5.3p5"
          #"paramiko >= 1.9.0",
          #"celery >= 2.2.7,<3.0.0",
          #"rpy2 >= 2.0.7"
      ])
