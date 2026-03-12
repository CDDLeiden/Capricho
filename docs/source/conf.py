import os
import sys

sys.path.insert(0, os.path.abspath("../../src"))

# Project information
project = "CAPRICHO"
copyright = "2024, David Araripe"
author = "David Araripe"

# The full version, including alpha/beta/rc tags
release = "0.1.0"

# Extensions
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx_autodoc_typehints",
    "myst_parser",
]

# MyST settings
myst_enable_extensions = [
    "colon_fence",
    "fieldlist",
]

# Templates path
templates_path = ["_templates"]

# List of patterns to ignore when looking for source files
exclude_patterns = []

# Theme
html_theme = "furo"
html_title = " "
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_theme_options = {
    "light_logo": "logo-light.svg",
    "dark_logo": "logo-dark.svg",
}

# Intersphinx mapping
intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False

# Source settings
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
