#!/usr/bin/env python
"""Setup file and install script for NextGen sequencing analysis scripts.
"""
from setuptools import setup, find_packages

setup(name = "bcbio-nextgen",
      version = "0.5a",
      author = "Brad Chapman",
      author_email = "chapmanb@50mail.com",
      description = "Automated, distributed next-gen sequencing pipeline; includes Galaxy interaction",
      license = "MIT",
      url = "https://github.com/chapmanb/bcbb/tree/master/nextgen",
      namespace_packages = ["bcbio"],
      packages = find_packages(),
      scripts = ['scripts/bcbio_nextgen.py',
                 'scripts/bam_to_wiggle.py',
                 'scripts/barcode_sort_trim.py',
                 'scripts/illumina_finished_msg.py',
                 'scripts/nextgen_analysis_server.py',
                 'scripts/solexa_qseq_to_fastq.py',
                 'scripts/upload_to_galaxy.py',
                 ],
      install_requires = [
          "numpy >= 1.6.2",
          "biopython >= 1.60",
          "Mako >= 0.3.6",
          "PyYAML >= 3.09",
          "sh >= 1.07",
          "Logbook >= 0.3",
          "Cython >= 0.17.3",
          "pysam >= 0.7",
          "fabric >= 1.5",
          "pyzmq >= 2.2.0.1",
          "ipython >= 0.13.1",
          "paramiko >= 1.9.0",
          "celery >= 2.2.7",
          #"rpy2 >= 2.0.7"
      ])
