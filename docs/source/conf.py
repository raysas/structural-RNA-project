# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))  # allows autodoc to find your lib


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'rna_score'
copyright = '2026, Yazid Hoblos, Joelle Assy, Denys Buryi, Raul Duran De Alba, Rayane Adam'
author = 'Yazid Hoblos, Joelle Assy, Denys Buryi, Raul Duran De Alba, Rayane Adam'
release = '0.1.0'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',          # generate docs from docstrings
    'sphinx.ext.viewcode',         # adds links to source code
    'sphinx.ext.napoleon',         # support Google/NumPy style docstrings
    'sphinx.ext.mathjax',          # render LaTeX math
    'sphinx_autodoc_typehints',    # show type hints nicely
    'myst_parser',                 # enable markdown
    'nbsphinx',                    # enable notebooks
    'sphinx_design',               # cards, buttons, dropdowns
]

templates_path = ['_templates']
exclude_patterns = ['_build', '**.ipynb_checkpoints']

# -- Options for HTML output -------------------------------------------------
html_theme = "pydata_sphinx_theme"

# Optional: configure the theme
html_theme_options = {
    "logo": {
        "image_light": "logo-light.png",
        "image_dark": "logo-dark.png",
        "text_logo_height": "72px",  # larger brand text for fallback
    },
    "navigation_depth": 3,
    "show_prev_next": True,
    "use_edit_page_button": False,
    "secondary_sidebar_items": {"*": ["page-toc"]},
}

html_context = {
    "github_user": "raysas",        # replace
    "github_repo": "structural-RNA-project",             # replace
    "github_version": "main",                     # or branch name
    "doc_path": "docs/source",                   # path to your docs folder
}



html_static_path = ['_static']  # keep this for CSS or images
html_css_files = ["custom.css"]

# -- MyST settings (Markdown) -----------------------------------------------

myst_enable_extensions = [
    "colon_fence",       # supports ::: blocks
    "deflist",           # definition lists
    "dollarmath",        # enable $...$ and $$...$$ math
]