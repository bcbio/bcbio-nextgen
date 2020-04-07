#!/usr/bin/env python

"""Setup file and install script for NextGen sequencing analysis scripts"""

import os
import subprocess

import setuptools

VERSION = '1.2.3'

# add bcbio version number and git commit hash of the current revision to version.py
try:
    git_run = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], stdout=subprocess.PIPE)
    git_run.check_returncode()
except subprocess.SubprocessError:
    commit_hash = ''
else:
    commit_hash = git_run.stdout.strip().decode()

here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, 'bcbio', 'pipeline', 'version.py'), 'w') as version_file:
    version_file.writelines([f'__version__ = "{VERSION}"\n',
                             f'__git_revision__ = "{commit_hash}"\n'])

# dependencies are installed via Conda from
# https://github.com/chapmanb/cloudbiolinux/blob/master/contrib/flavor/ngs_pipeline_minimal/packages-conda.yaml
setuptools.setup(version=VERSION)
