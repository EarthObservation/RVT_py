.. _creating_vis:

Creating a visualization
========================

Visualizations with rvt.vis
---------------------------

The module ``rvt.vis`` contains the ``rvt`` visualization functions. 

Every function takes a DEM (as 2D numpy array) with parameters and outputs a visualization (as 2D numpy array).

For example , to calculate a hillshade, first import the module, read the DEM (:ref:`Reading and saving raster`), call the ``rvt.vis.hillshade()`` function with its parameters, and save the visualization (:ref:`Reading and saving raster`). 

Example (calculating a hillshade with sun azimuth 315° and sun elevation 35°):

.. code-block:: python

    import rvt.vis

    hillshade_arr = rvt.vis.hillshade(dem=dem_arr, sun_azimuth=315, sun_elevation=35, resolution_x=dem_x_resolution, resolution_y=dem_y_resolution, no_data=dem_no_data)

.. note:: Find out more about visualization functions and their parameters in :ref:`rvt.vis`.

----

Visualizations with rvt.default (beginner)
------------------------------------------

For beginner python users we suggest using ``rvt.default`` instead of ``rvt.vis`` to calculate and store visualizations.

As well as containing functions to read and save rasters, the ``rvt.default`` also contains the class ``DefaultValues()`` where we can store our visualization functions parameters. We can call the methods of this class for saving and computing visualizations with those parameters (these methods use ``rvt.vis`` for computing visualizations).

For example, to calculate and get or save a hillshade with the default module, first import the module and create a ``DefaultValues()`` class instance. Then we can change the default parameters for a hillshade (they are attributes of ``DefaultValues()``, their name starts with ``hs_``). After that, call the method to get the hillshade array or to save hillshade to GeoTIFF. 

Example:

.. code-block:: python

    import rvt.default

    # create DefaultValues() instance
    default = rvt.default.DefaultValues()
    # change hillshade parameters default values to our needs
    default.hs_sun_el = 45
    default.hs_sun_azi = 300
    # call default.get_hillshade() method which uses set parameters and returns hillshade numpy array
    hillshade_arr = default.get_hillshade(dem_arr=dem_arr, resolution_x=dem_x_resolution, resolution_y=dem_y_resolution, no_data=dem_no_data)
    # if we don't need hillshade array and we just want to save hillshade we can directly call default.save_hillshade() method
    # this method also uses set hillshade parameters and saves visualization as GeoTIFF in dem_path directory
    default.save_hillshade(dem_path=dem_path, save_float=True, save_8bit=True)  # if we want also 8bit version of result we set save_8bit=True

Parameters of a ``DefaultValues()`` instance can be saved to a ``JSON`` configuration file which can be edited. You can then load this file back and overwrite the attribute values (visualization functions parameters).

Example:

.. code-block:: python

    import rvt.default

    default = rvt.default.DefaultValues()
    config_json_path = r"C:/rvt_default_values.json"  # change path to where you would like to save config file
    # save set attributes values to JSON configuration file
    default.save_default_to_file(file_path=config_json_path)
    # overwrite DefaultValues() instance (default) attributes values from config file
    default.read_default_from_file(file_path=config_json_path)
    
The ``DefaultValues()`` class also contains the methods: ``get_slope()``, ``save_slope()``, ``get_multi_hillshade()``, ``save_multi_hillshade()``, ``get_slrm()``,
``save_slrm()``, ``get_sky_view_factor()``, ``save_sky_view_factor()``, ``get_neg_opns()``, ``save_neg_opns()``, ``get_local_dominance()``, ``save_local_dominance()``,
``get_sky_illumination()``, ``save_sky_illumination()``.

Find out more about the methods and attributes of the ``DefaultValues()`` class in :ref:`rvt.default`.

