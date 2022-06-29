Default values
==============

For beginner python users we suggest using ``rvt.default`` instead of ``rvt.vis`` to calculate and store visualizations.

The ``rvt.default`` module contains functions to read and save rasters. This module was initially developed for GUI backend use. The default module also contains the class ``DefaultValues()`` where we can store our visualization functions parameters. We can then call the methods of this class for saving and computing visualizations with those parameters (these methods use ``rvt.vis`` for computing visualizations).

For example, to calculate and get or save a hillshade with the default module, first import the module and create a ``DefaultValues()`` class instance. Then we can change the default parameters for a hillshade (they are attributes of ``DefaultValues()``, their name starts with ``hs_``). After that, call the method to get the hillshade array or to save hillshade to GeoTIFF. See example below:

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

Parameters of a ``DefaultValues()`` instance can be saved to a ``JSON`` configuration file which can be edited. You can then load this file back and overwrite the attribute values (visualization functions parameters). See example below:

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

----

Find out more about the methods and attributes of the ``DefaultValues()`` class in :ref:`rvt.default`.
