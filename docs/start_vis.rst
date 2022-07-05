.. _start_vis:

Creating a visualization
========================

Visualizations can be created with both the ``rvt.vis`` and ``rvt.default`` modules.

----

Visualizations with rvt.vis
---------------------------

The module ``rvt.vis`` contains the ``rvt`` visualization functions. 

Every function takes a DEM (as a 2D numpy array) with parameters, and outputs a visualization (as a 2D numpy array).

**Example**

Let's say we need to calculate a hillshade with sun azimuth 315° and sun elevation 35°:

.. code-block:: python

    # import the module
    import rvt.vis
    
    # read the DEM 
    # Follow the steps in `Reading raster data`
    
    # call the rvt.vis.hillshade() function with its parameters
    hillshade_arr = rvt.vis.hillshade(
        dem=dem_arr, 
        sun_azimuth=315, 
        sun_elevation=35, 
        resolution_x=dem_x_resolution, 
        resolution_y=dem_y_resolution, 
        no_data=dem_no_data
        )
    
    # save the visualization 
    # Follow the steps in `Saving raster data`

.. seealso:: Find out more about visualization functions and their parameters in :ref:`rvt.vis`.

----

Visualizations with rvt.default (beginner)
------------------------------------------

For beginner Python users we suggest using ``rvt.default`` instead of ``rvt.vis`` to calculate and store visualizations.

As well as containing functions to read and save rasters, ``rvt.default`` also contains the class ``DefaultValues()`` where we can store our visualization functions parameters. We can call the methods of this class for saving and computing visualizations with those parameters (these methods use ``rvt.vis`` for computing visualizations).

**Example**

To calculate, get or save a hillshade using ``rvt.default``:

.. code-block:: python

    #  import the module 
    import rvt.default

    # create a DefaultValues() class instance
    default = rvt.default.DefaultValues()
    
    # change hillshade parameters default values to our needs 
    # (they are attributes of DefaultValues(), their name starts with hs_)
    default.hs_sun_el = 45
    default.hs_sun_azi = 300
    
    # call the method default.get_hillshade() which uses the set parameters and returns the hillshade numpy array
    hillshade_arr = default.get_hillshade(
        dem_arr=dem_arr, 
        resolution_x=dem_x_resolution, 
        resolution_y=dem_y_resolution, 
        no_data=dem_no_data
        )
    
    # if we don't need the hillshade array and we just want to save the
    # hillshade, we can directly call the default.save_hillshade() method
    # this method also uses the set hillshade parameters and saves the 
    # visualization as a GeoTIFF in the dem_path directory
    # to save the 8bit version of the result, set save_8bit=True
    default.save_hillshade(
        dem_path=dem_path, 
        save_float=True, 
        save_8bit=True
        )  
   
Configuring visualization parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Parameters of a ``DefaultValues()`` instance can be saved to a ``JSON`` configuration file which can be edited. You can then load this file back and overwrite the attribute values (or visualization functions parameters).

**Example**

.. code-block:: python

    # import the module
    import rvt.default

    # create a DefaultValues() class instance
    default = rvt.default.DefaultValues()
    
    # change this path to where you would like to save the config file
    config_json_path = r"C:/rvt_default_values.json"
    
    # save set attributes values to a JSON configuration file
    default.save_default_to_file(file_path=config_json_path)
    
    # overwrite the DefaultValues() instance (default) attributes values from the config file
    default.read_default_from_file(file_path=config_json_path)
   
DefaultValues() class methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
The ``DefaultValues()`` class also contains the methods: ``get_slope()``, ``save_slope()``, ``get_multi_hillshade()``, ``save_multi_hillshade()``, ``get_slrm()``, ``save_slrm()``, ``get_sky_view_factor()``, ``save_sky_view_factor()``, ``get_neg_opns()``, ``save_neg_opns()``, ``get_local_dominance()``, ``save_local_dominance()``, ``get_sky_illumination()``, ``save_sky_illumination()``.

.. seealso:: Find out more about the methods and attributes of the ``DefaultValues()`` class in :ref:`rvt.default`.
