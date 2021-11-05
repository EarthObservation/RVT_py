# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import shutil

# Needed for removing credits
import sphinx.ext.autodoc

sys.path.insert(0, os.path.abspath('..'))  # Source code dir relative to this file

# -- Project information -----------------------------------------------------

project = 'Relief Visualization Toolbox in Python'
copyright = '2010-2021, ZRC SAZU and University of Ljubljana'
author = 'Žiga Kokalj, Krištof Oštir, Klemen Zakšek and Žiga Maroh'

# -- General configuration ---------------------------------------------------

master_doc = 'index'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',  # Core library for html generation from docstrings
    'sphinx.ext.autosummary',  # Create neat summary tables
    'sphinx.ext.viewcode',  # Add links to highlighted source code
    'sphinx.ext.napoleon',  # Use NumPy docstrings
    'sphinxcontrib.bibtex', # BibTeX support
    'nbsphinx'  # Jupyter Notebook support
]

autodoc_member_order = 'bysource'  # Content is in the same order as in module
autosummary_generate = True  # Turn on sphinx.ext.autosummary
autoclass_content = "both"  # Add __init__ doc (ie. params) to class summaries
html_show_sourcelink = False  # Remove 'view source code' from top of page (for html, not python)
autodoc_inherit_docstrings = True  # If no class summary, inherit base class summary

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# Ignore import errors
# nbsphinx_allow_errors = True

# BibTeX files
bibtex_bibfiles = ['RVT.bib']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'sphinx_rtd_theme'
html_theme = 'furo'

# Logo
html_logo = './figures/RVT_head.png'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = []

# Remove module docstring credits
def setup(app):
    app.connect('autodoc-process-docstring', sphinx.ext.autodoc.between('Credits:', what=['module'], exclude=True))

# Autommatiyally remove GDAL from imports
autodoc_mock_imports = ['gdal', 'osgeo']

# Copy examples to docs
examples_folder = './examples'
shutil.rmtree(examples_folder, ignore_errors=True)

try:
    shutil.copytree('../examples', examples_folder)
except:
    raise Exception('Error: Cannot copy examples to Docs')
