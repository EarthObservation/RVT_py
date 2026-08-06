.. _install:

Installation
============

RVT can be installed as a Python package for use in Python scripts and Jupyter Notebooks. Separate RVT distributions are available as custom raster functions for ArcGIS Pro and as a plugin for QGIS.

Python package requirements
---------------------------

The RVT Python package requires Python 3.11 or newer and these runtime libraries:

* NumPy
* SciPy
* GDAL
* Matplotlib

GDAL includes native libraries. Conda is therefore the recommended installation method for most users. A PyPI installation requires a compatible GDAL native library and matching Python bindings.

The notebook and repository examples use additional optional packages. They are not required for the installed ``rvt`` package.

----

**CONTENTS**

.. toctree::

   install_python
   install_arcgis
   install_qgis
