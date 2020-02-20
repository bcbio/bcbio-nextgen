#!/usr/bin/env python

"""Setup file and install script for NextGen sequencing analysis scripts"""

import os
import subprocess

import setuptools

__version__ = '1.2.0'

with open('README.rst', 'r') as readme_file:
    long_description = readme_file.read()

scripts = ['scripts/bcbio_nextgen.py',
           'scripts/bcbio_setup_genome.py',
           'scripts/bcbio_prepare_samples.py',
           'scripts/bcbio_fastq_umi_prep.py',
           'scripts/cwltool2wdl.py']


def write_version_py():
    version_py = os.path.join(os.path.dirname(__file__), 'bcbio', 'pipeline', 'version.py')
    try:
        git = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], stdout=subprocess.PIPE)
        git.check_returncode()
    except subprocess.SubprocessError:
        commit_hash = ''
    else:
        commit_hash = git.stdout.strip().decode()
    with open(version_py, 'w') as version_file:
        version_file.writelines([f'__version__ = "{__version__}"\n',
                                 f'__git_revision__ = "{commit_hash}"\n'])


write_version_py()

setuptools.setup(
    name='bcbio-nextgen',
    version=__version__,
    author='bcbio community',
    author_email='biovalidation@googlegroups.com',
    description='Best-practice pipelines for fully automated high throughput sequencing analysis',
    long_description=long_description,
    license='MIT',
    url='https://github.com/bcbio/bcbio-nextgen',
    packages=setuptools.find_packages(exclude=['tests']),
    zip_safe=False,
    scripts=scripts,
    install_requires=[],  # install dependencies via conda
    include_package_data=True
)
