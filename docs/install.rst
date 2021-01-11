.. _install:

Installation
============

conda
-----

You can install all required libraries and rvt-py with Anaconda environment from Anaconda cloud (`conda rvt_py <https://anaconda.org/zmigyyy/rvt_py>`_). To do that open Anaconda Prompt (activate conda environment) and run:

`conda install -c zmigyyy rvt_py`

pypi
----

Another option is to install required libraries and rvt-py with Python Package Index, pypi (`pypi rvt-py <https://pypi.org/project/rvt-py>`_). To do that open command prompt (terminal) and run:

`pip install rvt-py`

This might not work, because pypi usually has problems installing gdal. To solve that first try to install gdal then run above command.

requirements
------------

We suggest using an Anaconda environment (easiest gdal installation) and Python 3.6 or higher. Required libraries with tested versions (could also work with other versions):

*   numpy 1.19.2
*   scipy 1.5.2
*   gdal 3.0.2


You can also clone the repository (`github rvt_py <https://github.com/EarthObservation/RVT_py>`_).

Library `rvt-py` can be used in Python scripts, Jupyter Notebooks and in ArcGIS Pro.


ArcGIS Raster functions
-----------------------

ArcGIS Raster functions implementation is described in :ref:`arcgis`.

QGIS Plugin
-----------

QGIS RVT plugin installation and usage is described in :ref:`qgis`.
