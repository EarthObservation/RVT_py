Installation
============

The RVT Python package can be installed using Conda or PyPI, and can be used in
Python scripts, Jupyter Notebooks and ArcGIS Pro. You can also clone the
repository (`GitHub rvt_py <https://github.com/EarthObservation/RVT_py>`_).

Requirements
------------

Version 3 requires Python 3.12 and the following libraries (rasterio replaces
the ``gdal`` dependency of earlier versions):

*   numpy
*   scipy
*   rasterio
*   matplotlib
*   affine

Conda
-----

The ``rvt`` package is `available from the Anaconda Cloud repository <https://anaconda.org/rvtpy/rvt_py>`_. Using Conda to install the ``rvt`` package will include all required libraries.

First `install Anaconda and Conda <https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html>`_.

Then open Anaconda Prompt (Windows) or Terminal (MacOS) and run:

``conda install -c rvtpy rvt_py``

PyPI
----

Another option is to install the ``rvt-py`` package and required libraries `using the Python Package Index <https://pypi.org/project/rvt-py>`_ (PyPI).

Open Command Prompt (Windows) or Terminal (MacOS) and run:

``pip install rvt-py``

Development install
-------------------

The project uses `uv <https://docs.astral.sh/uv/>`_:

.. code-block:: bash

   uv sync
   uv run pytest -n auto
   uv run ruff check .
   uv run ruff format .
   uv run pyrefly check

The test suite verifies the implementation against golden-master fixtures in
``tests/golden/``, generated with the version 2 code;
``tests/golden/manifest.json`` records the exact parameters and provenance of
every fixture.

GIS applications
----------------

Besides the Python library, RVT is distributed as a set of custom raster
functions for ArcGIS and as a QGIS plugin. Both are separate projects with
their own release cycle and are built on the version 2 code base; they are
documented in their own repositories.

*   `ArcGIS Raster Functions repository <https://github.com/EarthObservation/rvt-arcgis-pro>`_ — download it with ``Code → Download ZIP``, unzip and *rename* the folder to ``rvt-arcgis-pro``, then copy it to ``<ArcGIS Pro install path>/Resources/Raster/Functions/Custom``.
*   `RVT QGIS plugin <https://plugins.qgis.org/plugins/rvt-qgis/>`_ — in QGIS select ``Plugins → Manage and Install Plugins → All``, search for ``RVT`` or ``Relief Visualization Toolbox`` and select ``Install Plugin``.
