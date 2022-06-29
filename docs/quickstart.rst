.. _quickstart:

Getting Started
===============

**CONTENTS**

.. toctree::
    quickstart
    
----

Introduction
------------
This section gives a basic explanation of how to use ``rvt``. 

More detailed explanations of how to use ``rvt`` can be found in :ref:`Examples` and also in the following files in the GitHub repository:

.. code-block:: python

    test_vis.py
    test_blend.py
    test_default.py
    test_custom_color_scheme.py

If you need some data to get started, a `small sample dataset <https://github.com/EarthObservation/RVT_py/blob/master/test_data/TM1_564_146.tif>`_ (2.43 MB) and a `larger sample dataset <https://rebrand.ly/rvt_demo>`_ (152 MB) are available for download. If you have cloned the repository, save the data in the ``test_data`` directory.

----

.. _Reading and saving raster:

Reading and saving raster data
------------------------------

For reading raster data (DEMs) from files (GeoTIFFs) to a numpy array we suggest using the ``rvt.default`` module (which uses ``gdal``).
You can also use ``rasterio``, ``gdal`` or any other module for reading and saving geo rasters.

To read rasters with ``rvt.default``, first import the module. Then call the function ``rvt.default.get_raster_arr()`` to get a dictionary with keys array (contains numpy array of raster), resolution (contains tuple(x resolution, y resolution)) and no_data (contains value of no_data). See example below:

.. code-block:: python

    import rvt.default

    dem_path = r"C:/data/dem.tif"  # change path to your GeoTIFF
    dem_dict = rvt.default.get_raster_arr(dem_path)  # returns dictionary: {"array": array, "resolution": (x_res, y_res), "no_data": no_data}
    dem_arr = dem_dict["array"]  # numpy array
    dem_resolution_tuple = dem_dict["resolution"]  # resolution tuple (x direction resolution, y direction resolution)
    dem_x_resolution = dem_resolution_tuple[0]  # first element of resolution tuple is x direction resolution
    dem_y_resolution = dem_resolution_tuple[1]  # second element of resolution tuple is y direction resolution
    dem_no_data = dem_dict["no_data"]  # returns value of no_data stored in DEM


To save a raster with ``rvt.default``, import the module then call the function ``rvt.default.save_raster()``. You have to define the function parameters: ``src_raster_path``: source raster path (dem_path) to copy metadata, ``out_raster_path``: path to new file (visualization tif), ``out_raster_arr``: vizualization numpy array, ``no_data``: value of no_data (visualizations return no data as np.nan).

For example, the example below would compute a hillshade (with rvt.vis) in hillshade_arr from a DEM stored in dem_path location, storing this hillshade visualization to hillshade_path:

.. code-block:: python

    import rvt.default
    import numpy as np

    rvt.default.save_raster(src_raster_path=dem_path, out_raster_path=hillshade_path, out_raster_arr=hillshade_arr, no_data=np.nan)

----

.. _module_vis:

Vis module
----------

The module ``rvt.vis`` contains visualization functions. Every function takes a DEM (as 2D numpy array) with parameters and outputs a visualization (as 2D numpy array).

For example , to calculate a hillshade, first import the module, read the DEM (:ref:`Reading and saving raster`), call the ``rvt.vis.hillshade()`` function with its parameters, and save the visualization (:ref:`Reading and saving raster`). See example below for calculating a hillshade with sun azimuth 315° and sun elevation 35°:

.. code-block:: python

    import rvt.vis

    hillshade_arr = rvt.vis.hillshade(dem=dem_arr, sun_azimuth=315, sun_elevation=35, resolution_x=dem_x_resolution, resolution_y=dem_y_resolution, no_data=dem_no_data)

Find out more about visualization functions and their parameters in :ref:`rvt.vis`.

----

.. _module_default:

Default module
--------------

For beginner python users we suggest using ``rvt.default`` instead of ``rvt.vis`` to calculate and store visualizations.

As mentioned before, ``rvt.default`` module contains functions to read and save rasters. This module was initially developed for GUI backend use. The default module also contains the class ``DefaultValues()`` where we can store our visualization functions parameters. We can then call the methods of this class for saving and computing visualizations with those parameters (these methods use ``rvt.vis`` for computing visualizations).

For example, to calculate and get or save a hillshade with the default module, we have to import the module and create a ``DefaultValues()`` class instance. Than we can change the default parameters for a hillshade (they are attributes of ``DefaultValues()``, their name starts with ``hs_``). After that we call the method to get the hillshade array or to save hillshade to GeoTIFF. See example below:

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

The ``DefaultValues()`` class also contains the methods: ``get_slope()``, ``save_slope()``, ``get_multi_hillshade()``, ``save_multi_hillshade()``, ``get_slrm()``,
``save_slrm()``, ``get_sky_view_factor()``, ``save_sky_view_factor()``, ``get_neg_opns()``, ``save_neg_opns()``, ``get_local_dominance()``, ``save_local_dominance()``,
``get_sky_illumination()``, ``save_sky_illumination()``.

Find out more about the methods and attributes of the ``DefaultValues()`` class in :ref:`rvt.default`.

Parameters of a ``DefaultValues()`` instance can be saved to a ``JSON`` configuration file which can be edited. You can then load this file back and overwrite the attribute values (visualization functions parameters). See example below:

.. code-block:: python

    import rvt.default

    default = rvt.default.DefaultValues()
    config_json_path = r"C:/rvt_default_values.json"  # change path to where you would like to save config file
    # save set attributes values to JSON configuration file
    default.save_default_to_file(file_path=config_json_path)
    # overwrite DefaultValues() instance (default) attributes values from config file
    default.read_default_from_file(file_path=config_json_path)

----

.. _module_blend:

Blend module
------------

You can blend manually or automatically. When blending manually you have to define each layer (visualization) in python. Manually blending allows you to use visualizations that are not part of ``rvt``. Automatically blending automatically computes visualizations (they need to be a part of ``rvt``) and blends them together from configuration ``JSON`` file (can be edited).

The main class of the ``rvt.blend`` module for blending is ``BlenderCombination`` which has list attribute ``layers`` where instances of class ``BlenderLayer`` are stored.

In ``BlenderLayer`` instances in ``layers`` we store a specific visualization and its parameters for blending. The ``BlenderCombination`` class has the method ``render_all_images()``, which blends together all ``BlenderLayer`` instances (visualizations) in the ``BlenderCombination.layers`` list and outputs the blended image.

Find out more about blending in :ref:`rvt.blend`.

Manual blending
^^^^^^^^^^^^^^^

When blending, you have to import the ``rvt.blend`` module and create a ``BlenderCombination`` instance. 

For adding a layer (visualization) with parameters to a combination, you can call ``BlenderCombination.create_layer()``. This creates a ``BlenderLayer`` instance and adds it to ``BlenderCombination.layers``.

For example, let's say you have already calculated the simple local relief model (slrm_arr), slope (slope_arr) and hillshade (hillshade_arr), and now you want to blend all the calculated visualizations together:

.. code-block:: python

    import rvt.blend

    # create combination class which will hold layers (visualizations)
    combination_manual = rvt.blend.BlenderCombination()

    # 1st layer: Add slrm layer with 2% perc cuttoff on both sides, multiply blend mode and 25% opacity
    combination_manual.create_layer(vis_method="Simple local relief model", normalization="perc", minimum=2, maximum=2,
                              blend_mode="multiply", opacity=25, image=slrm_arr)
    # 2nd layer: Add slope layer with value stretch from 0 to 51, luminosity blend mode and 50% opacity
    combination_manual.create_layer(vis_method="Slope gradient", normalization="value", minimum=0, maximum=51,
                              blend_mode="luminosity", opacity=50, image=slope_arr)
    # 3rd layer: Add hillshade layer with value stretch from 0 to 1, normal blend mode and 100% opacity
    combination_manual.create_layer(vis_method="Hillshade", normalization="value", minimum=0, maximum=1,
                              blend_mode="normal", opacity=100, image=hillshade_arr)

    # if we wish to save blended image in file we have to add dem_path to combination (for metadata, geodata)
    combination_manual.add_dem_path(dem_path=input_dem_path)

    # blend them all together, you can save blend to GeoTIFF if save_render_path presented (and dem_path is added) else it only returns array
    render_arr = combination_manual.render_all_images(save_render_path=output_blend_path)

You can also let the ``BlenderCombination`` class automatically compute the visualization or give the path to a visualization. 

If you don't provide parameter image, and vis_method parameter is correct (existing rvt visualization function), blender automatically calculates the visualization. 

If you provide parameter image_path and not image (if you provide both image will be used), blender will read visualization from image_path. 

If you don't input image and image_path parameter, you have to add an ``rvt.default.DefaultValues`` instance as a parameter to ``BlenderCombination.render_all_images()``. Blender then takes parameters set in this class when calculating specific visualization.
You also have to add dem array and its resolution. See example below which uses all three methods:

.. code-block:: python

    import rvt.blend
    import rvt.default

    # create combination class which will hold layers (visualizations)
    combination_manual = rvt.blend.BlenderCombination()

    # we will let blender to calculate slrm so we need to create rvt.default.DefaultValues() and change parameters of
    # slrm, we will later add default to combination_manual.render_all_images() method
    default = rvt.default.DefaultValues()
    default.slrm_rad_cell = 15

    # 1st layer: Add slrm layer with 2% perc cuttoff on both sides, multiply blend mode and 25% opacity
    # slrm is calculated automatically, because we didn't provide image and image_path parameters
    combination_manual.create_layer(vis_method="Simple local relief model", normalization="perc", minimum=2, maximum=2,
                              blend_mode="multiply", opacity=25)
    # 2nd layer: Add slope layer with value stretch from 0 to 51, luminosity blend mode and 50% opacity
    # we provide image_path to slope, so slope is read from file
    combination_manual.create_layer(vis_method="Slope gradient", normalization="value", minimum=0, maximum=51,
                              blend_mode="luminosity", opacity=50, image_path=slope_path)
    # 3rd layer: Add hillshade layer with value stretch from 0 to 1, normal blend mode and 100% opacity
    # we provide image
    combination_manual.create_layer(vis_method="Hillshade", normalization="value", minimum=0, maximum=1,
                              blend_mode="normal", opacity=100, image=hillshade_arr)

    # we have to add dem array and resolution so that slrm can be computed
    combination_manual.add_dem_arr(dem_arr=input_dem_arr, dem_resolution=resolution)

    # blend them all together and add default where are defined slrm parameters
    render_arr = combination_manual.render_all_images(default=default)


You can always add as many layers as you want.

Automatic blending
^^^^^^^^^^^^^^^^^^

Automatic blending is blending from a configuration ``JSON`` file. You can create an example file and change it to suit your needs.
To blend from a file we create the ``BlenderCombination()`` class, call the method ``read_from_file()`` and then ``render_all_images()``. In the ``render_all_images()`` method we can save (to dem_path directory) a specific visualization if we set the parameter ``save_visualization`` to True.

.. code-block:: python

    import rvt.blend

    combination_auto = rvt.blend.BlenderCombination()
    # to create JSON blender combination configuration file example you can change
    blender_combination_path = r"settings\blender_file_example.txt"  # change path to where you wish to save
    rvt.blend.create_blender_file_example(file_path=blender_combination_path)

    # set parameters of visualizations you will be using
    default = rvt.default.DefaultValues()
    # for example default.hs_sun_el=40

    # read json combination file from JSON
    combination_auto.read_from_file(file_path=blender_combination_path)

    layers_auto.add_dem_path(input_dem_path)  # needed when save_visualizations is True and save_rander_path is not None

    layers_auto.render_all_images(default=default, save_visualizations=True, save_render_path=output_blend_path,
                              save_float=True, save_8bit=True)  # if you also wish to save 8bit version
