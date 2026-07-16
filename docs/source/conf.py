"""Configuration file for the Sphinx documentation builder.
For the full list of built-in configuration values, see the documentation:
https://www.sphinx-doc.org/en/master/usage/configuration.html
"""

import os
import sys

import versioningit
from packaging.version import Version

sys.path.insert(0, os.path.abspath("../../src/lr_reduction"))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "lr_reduction"
copyright = "2024, ORNL"  # noqa A001
author = "ORNL"
version = versioningit.get_version("../../")
# The full version (major.minor.patch) without pre-/post-release metadata
release = Version(version).base_version

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.coverage",
    # myst_parser enables markdown support
    "myst_parser",
]

autodoc_mock_imports = [
    "mantid",
    "mantid.api",
    "mantid.kernel",
    "mantid.utils",
    "mantid.utils.logging",
    "mantid.simpleapi",
    "mantid.geometry",
    "mantidqt.widgets",
    "mantidqt.widgets.algorithmprogress",
    "qtpy",
    "qtpy.uic",
    "qtpy.QtWidgets",
    "mantidqt",
    "mantid.plots",
    "mantid.plots.plotfunctions",
    "mantid.plots.datafunctions",
    "mantid.plots.utility",
    "requests",
]

master_doc = "index"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "sphinx": ("https://www.sphinx-doc.org/en/master/", None),
    "mantid": ("https://docs.mantidproject.org/nightly/", None),
    "requests": ("https://requests.readthedocs.io/en/latest/", None),
}
intersphinx_disabled_domains = ["std"]

# Suppress warnings for references that can't be resolved due to mocked imports
suppress_warnings = [
    "ref.class",  # Suppress class reference warnings for mocked modules like mantid
]

# Nitpicky mode (-n) reports unresolvable references as warnings. TypeVars used
# in Generic[...] base classes are rendered as py:obj cross-references by autodoc
# but are not documented objects themselves, so ignore them explicitly here.
nitpick_ignore = [
    ("py:obj", "lr_reduction.processing.interfaces.DataT"),
    ("py:obj", "lr_reduction.processing.interfaces.ConfigT"),
    ("py:obj", "lr_reduction.processing.interfaces.OutT"),
]

templates_path = ["_templates"]
exclude_patterns = ["_build"]
pygments_style = "sphinx"

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# MyST parser configuration
myst_enable_extensions = [
    "deflist",
    "dollarmath",
    "smartquotes",
    "html_admonition",
    "html_image",
]

# Allow heading anchors
myst_heading_anchors = 3

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"  # pylint: disable=C0103

# collapse_navigation=False keeps nested toctree entries (e.g. the "Workflow"
# sub-group) expanded in the sidebar at all times, instead of only expanding
# them when browsing a page within that section.
html_theme_options = {
    "style_nav_header_background": "#472375",
    "collapse_navigation": False,
}

epub_show_urls = "footnote"  # pylint: disable=C0103
