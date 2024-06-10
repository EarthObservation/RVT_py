.. _examples:

Example notebooks
=================

This section contains example code for using the ``rvt.vis`` and ``rvt.default`` modules in two Jupyter Notebooks. 

Both notebooks are available in the `examples folder in the repository <https://github.com/EarthObservation/RVT_py/tree/master/examples>`_ .


Tiled processing
++++
We have prepared two Jupyter Notebooks for processing of large datasets (to large to fit in RAM).
The notebooks can also be used for processing large number of GeoTIFF files with one click.

The underlying Python code splits large datasets into tiles and uses multiprocessing to speed up the computation.
The code can be executed through an `GUI <https://github.com/EarthObservation/RVT_py/tree/master/examples/tiled_processing/rvt_tiled_GUI.ipynb>`_  based on `ipywidgets`.
If you want to integrate the processing into your Python Workflow, see the `rvt_tiled.ipynb <https://github.com/EarthObservation/RVT_py/tree/master/examples/tiled_processing/rvt_tiled.ipynb>`_ Jupyter Notebook.


Both notebooks for **tiled processing** are available in the `tiled_processing <https://github.com/EarthObservation/RVT_py/tree/master/examples/tiled_processing>`_ folder.

----

**CONTENTS**

.. toctree::

    examples/rvt_vis_example.ipynb
    examples/rvt_default_example.ipynb
    examples/tiled_processing/rvt_tiled.ipynb
    examples/tiled_processing/rvt_tiled_GUI.ipynb
