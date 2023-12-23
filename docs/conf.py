# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os

from recommonmark.transform import AutoStructify

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'bcbio-nextgen'
copyright = '2021, bcbio-nextgen contributors'
author = 'bcbio-nextgen contributors'

# The full version, including alpha/beta/rc tags
version = release = '1.2.9'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.viewcode',
    'sphinx.ext.autosectionlabel',
    'myst_parser',
]

# Prefix document path to section labels, otherwise autogenerated labels would look like 'heading'
# rather than 'path/to/file:heading'
autosectionlabel_prefix_document = True

# source_suffix = {
#    '.md': 'markdown',
#}

# Add any paths that contain templates here, relative to this directory.
templates_path = ['templates']

# The master toctree document (required by Read the Docs).
master_doc = 'index'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build']


# -- Options for HTML output -------------------------------------------------

# to use Read The Docs theme when building docs locally
if not os.getenv('READTHEDOCS'):
    try:
        import sphinx_rtd_theme
    except ModuleNotFoundError:
        pass
    else:
        html_theme = 'sphinx_rtd_theme'
        html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]


# enable automatic table of contents on the index page
def setup(app):
    app.add_config_value('recommonmark_config', {
        'auto_toc_tree_section': 'Contents',
        'auto_toc_maxdepth': 2,
    }, True)
    app.add_transform(AutoStructify)
