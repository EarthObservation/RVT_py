"""
Sphinx configuration for the Relief Visualization Toolbox documentation.

The version is read from the installed ``rvt`` package, so the docs always
describe the code they were built against.
"""

from datetime import datetime

import rvt

# -- Project information -----------------------------------------------------

project = "Relief Visualization Toolbox"
author = "ZRC SAZU and University of Ljubljana (UL FGG)"
copyright = f"2010-{datetime.now().year}, {author}"
release = rvt.__version__
version = release

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
root_doc = "index"

rst_epilog = f"""
.. |year| replace:: {datetime.now().year}
"""

# -- Autodoc -----------------------------------------------------------------

# Content is documented in the order it appears in the source, and the frozen
# configs' constructor parameters (their reST field lists) are part of the
# class documentation.
autodoc_member_order = "bysource"
autoclass_content = "both"
autodoc_inherit_docstrings = True
autodoc_typehints = "description"
autodoc_typehints_description_target = "documented_params"

# -- Intersphinx -------------------------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "rasterio": ("https://rasterio.readthedocs.io/en/stable/", None),
}

# -- HTML output -------------------------------------------------------------

html_theme = "furo"
html_title = f"{project} {release}"
html_logo = "figures/RVT_head.png"
html_static_path = []
html_show_sourcelink = False
