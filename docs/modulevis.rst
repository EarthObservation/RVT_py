.. _module_vis:

Vis module
==========

The module ``rvt.vis`` contains the ``rvt`` visualization functions. 

Every function takes a DEM (as 2D numpy array) with parameters and outputs a visualization (as 2D numpy array).

For example , to calculate a hillshade, first import the module, read the DEM (:ref:`Reading and saving raster`), call the ``rvt.vis.hillshade()`` function with its parameters, and save the visualization (:ref:`Reading and saving raster`). See example below for calculating a hillshade with sun azimuth 315° and sun elevation 35°:

.. code-block:: python

    import rvt.vis

    hillshade_arr = rvt.vis.hillshade(dem=dem_arr, sun_azimuth=315, sun_elevation=35, resolution_x=dem_x_resolution, resolution_y=dem_y_resolution, no_data=dem_no_data)

----

Find out more about visualization functions and their parameters in :ref:`rvt.vis`.
