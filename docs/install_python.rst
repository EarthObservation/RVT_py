.. _install_python:

Python installation
===================

Conda
-----

The ``rvt`` package is `available from the Anaconda Cloud repository <https://anaconda.org/zmigyyy/rvt_py>`_. Using Conda to install the ``rvt`` package will include all required libraries.

First `install Anaconda and Conda <https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html>`_.

Then open Anaconda Prompt (Windows) or Terminal (MacOS) and run:

``conda install -c zmigyyy rvt_py``

PyPI
----

Another option is to install the ``rvt-py`` package and required libraries `using the Python Package Index <https://pypi.org/project/rvt-py>`_ (PyPI).

PyPI usually has problems installing ``gdal``, so `install gdal first <https://pypi.org/project/GDAL/>`_.

Then open Command Prompt (Windows) or Terminal (MacOS) and run:

``pip install rvt-py``