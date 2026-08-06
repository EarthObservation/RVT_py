.. _install_python:

Python installation
===================

RVT requires Python 3.11 or newer. Because GDAL includes native libraries, Conda is the recommended installation method for most users.

Conda
-----

The current RVT Conda package is published on the `rvtpy channel <https://anaconda.org/rvtpy/rvt_py>`_. Its dependencies should be installed from conda-forge.

Create a new environment:

.. code-block:: console

   conda create --name rvt --override-channels --channel conda-forge --channel rvtpy python=3.11 rvt_py
   conda activate rvt

PyPI
----

The RVT distribution name on the `Python Package Index <https://pypi.org/project/rvt-py>`_ is ``rvt-py``. The import package name is ``rvt``.

A standard pip installation works when a compatible GDAL native library and its Python bindings are already available:

.. code-block:: console

   python -m pip install rvt-py

If GDAL is not already installed, use the Conda method above. Alternatively, let Conda supply the runtime dependencies and install RVT itself from PyPI:

.. code-block:: console

   conda create --name rvt-pip --channel conda-forge python=3.11 gdal numpy scipy matplotlib pip
   conda activate rvt-pip
   python -m pip install --no-deps rvt-py

The ``--no-deps`` option prevents pip from replacing the packages supplied by Conda.

Optional example dependencies
-----------------------------

The installed ``rvt`` package does not require the notebook and example dependencies. In a pip-managed environment, install them only when needed:

.. code-block:: console

   python -m pip install "rvt-py[examples]"

Installation check
------------------

Verify that RVT and GDAL import from the active environment:

.. code-block:: console

   python -c "from osgeo import gdal; import rvt; import rvt.vis; import rvt.blend"
