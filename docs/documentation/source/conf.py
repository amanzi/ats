# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'ATS'
copyright = '202X, jointly shared by contributor institutions'
author = 'Ethan Coon'

# The short X.Y version
version = '1.4'
# The full version, including alpha/beta/rc tags
release = '1.4'


# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
#    'sphinx.ext.autodoc',
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.todo',
    'sphinx.ext.intersphinx',
#    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
#    'sphinx.ext.githubpages',
#    'sphinx.ext.napoleon',
    'sphinxcontrib.jquery',
    'nbsphinx',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
language = 'en'

exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'sphinx_book_theme'
html_theme = 'pydata_sphinx_theme'

html_logo = "_static/images/logo_full.png"
html_title = "Amanzi-ATS"

html_sidebars = {
    "**" : ["sidebar-nav-bs.html",
            "page-toc.html",
            ]
}

html_theme_options = {
  "secondary_sidebar_items": []
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


html_css_files = [
    'https://cdn.datatables.net/2.0.8/css/dataTables.dataTables.css',
    'styles/custom_theme.css',
]

html_js_files = [
    'https://cdn.datatables.net/2.0.8/js/dataTables.js',
    'main.js',
]



# -- Options for intersphinx extension ---------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {'amanzi': ('https://amanzi.github.io/amanzi', None),
                       'ats_demos': ('https://amanzi.github.io/ats_demos', None),
                      }


