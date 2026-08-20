.. _examples:

Example notebooks
=================

This section contains four Jupyter notebooks. 
The rvt_vis_example.ipynb and rvt_default_example.ipynb contain code for using the basic functions of the ``rvt.vis`` and ``rvt.default`` modules and are available in the `examples folder in the repository <https://github.com/EarthObservation/RVT_py/tree/master/examples>`_ .

The notebooks for **tiled processing** can be used for easy visualisation of large datasets and for computing many visualisations on a list of GeoTIFF files with one click. They are available in the `tiled_processing <https://github.com/EarthObservation/RVT_py/tree/master/examples/tiled_processing>`_ folder.

Tiled processing
++++++++++++++++
The Python code splits large datasets into tiles and uses multiprocessing to speed up the computation. 

The code can be executed through a `GUI <https://github.com/EarthObservation/RVT_py/tree/master/examples/tiled_processing/rvt_tiled_GUI.ipynb>`_ based on `ipywidgets`.

If you want to integrate the processing into your Python Workflow, see the `rvt_tiled.ipynb <https://github.com/EarthObservation/RVT_py/tree/master/examples/tiled_processing/rvt_tiled.ipynb>`_.


----

**CONTENTS**

.. toctree::

    examples/rvt_vis_example.ipynb
    examples/rvt_default_example.ipynb
    examples/tiled_processing/rvt_tiled.ipynb
    examples/tiled_processing/rvt_tiled_GUI.ipynb
