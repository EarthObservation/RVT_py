.. _usage:

Usage
=====

Reading and saving raster
----------

For reading raster (DEM) from files (GeoTIFF) to numpy array we suggest you to use our ``rvt.default`` module (uses gdal).
You can also use ``rasterio``, ``gdal`` or any other module for reading and saving geo rasters.

To read raster with ``rvt.default`` you first have to import module.
Than call function ``rvt.default.get_raster_arr()`` to get dictionary with keys array (contains numpy array of raster),
resolution (contains tuple(x resolution, y resolution)) and no_data (contains value of no_data). See example below:

.. code-block:: python

    import rvt.default

    dem_path = r"C:/data/dem.tif"  # change path to your GeoTIFF
    dem_dict = rvt.default.get_raster_arr(dem_path)  # returns dictionary: {"array": array, "resolution": (x_res, y_res), "no_data": no_data}
    dem_arr = dem_dict["array"]  # numpy array
    dem_resolution_tuple = dem_dict["resolution"]  # resolution tuple (x direction resolution, y direction resolution)
    dem_x_resolution = dem_resolution_tuple[0]  # first element of resolution tuple is x direction resolution
    dem_y_resolution = dem_resolution_tuple[1]  # second element of resolution tuple is y direction resolution
    dem_no_data = dem_dict["no_data"]  # returns value of no_data stored in DEM


To save raster with ``rvt.default`` you also have to have imported module. Then you call function ``rvt.default.save_raster()``.
You have to define function parameters: src_raster_path: source raster path (dem_path) to copy metadata, out_raster_path: path to new file (visualization tif), out_raster_arr: vizualization numpy array, no_data: value of no_data (visualizations return no data as np.nan).
For example if you compute hillshade (with rvt.vis) in hillshade_arr from dem stored in dem_path location. And you wish to store this hillshade visualiztion to hillshade_path. Saving hillshade would look like:

.. code-block:: python

    import rvt.default
    import numpy as np

    rvt.default.save_raster(src_raster_path=dem_path, out_raster_path=hillshade_path, out_raster_arr=hillshade_arr, no_data=np.nan)


.. _module_vis:

Module vis
----------

Module ``rvt.vis`` contains visualization functions. Every function takes dem (as 2D numpy array) with parameters and outputs visualization (as 2D numpy array).
For example to calculate hillshade you first have to import modules, read DEM (:ref:`Reading and saving raster`) and than call ``rvt.vis.hillshade()`` function with its parameters, then you can save visualization (:ref:`Reading and saving raster`).
Example of calculating Hillshade with sun azimuth 315° and sun elevation 35°:

.. code-block:: python

    import rvt.vis

    hillshade_arr = rvt.vis.hillshade(dem=dem_arr, sun_azimuth=315, sun_elevation=35, resolution_x=dem_x_resolution, resolution_y=dem_y_resolution, no_data=dem_no_data)


To find more about visualization functions and to learn about their parameters look into :ref:`rvt.vis`.

.. _module_default:

Module default
--------------

For beginner python users we suggest using ``rvt.default`` instead of ``rvt.vis`` to calculate and store visualizations.
As mentioned before ``rvt.default`` module contains functions to read and save rasters. This module was initially developed for GUI backend use.
Default module also contains class ``DefaultValues()`` where we can store our visualization functions parameters.
We can than call methods of this class for saving and computing visualizations with that parameters (methods use ``rvt.vis`` for computing visualizations).


For example to calculate and get or save hillshade with default module, we have to import module and create ``DefaultValues()`` class instance. Than we can change default parameters for hillshade (they are attributes of ``DefaultValues()``, their name starts with ``hs_``).
After that we call method to get hillshade array or to save hillshade to GeoTIFF. Example:

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


Class ``DefaultValues()`` also contains methods: ``get_slope()``, ``save_slope()``, ``get_multi_hillshade()``, ``save_multi_hillshade()``, ``get_slrm()``,
``save_slrm()``, ``get_sky_view_factor()``, ``save_sky_view_factor()``, ``get_neg_opns()``, ``save_neg_opns()``, ``get_local_dominance()``, ``save_local_dominance()``,
``get_sky_illumination()``, ``save_sky_illumination()``. Additional info (about methods and attributes of ``DefaultValues()`` class) is in :ref:`rvt.default`.


Parameters of ``DefaultValues()`` instance can be saved to JSON configuration file which can be edited. Then you can load this file back and overwrite attributes (visualization functions parameters) values.
Example how to do that:

.. code-block:: python

    import rvt.default

    default = rvt.default.DefaultValues()
    config_json_path = r"C:/rvt_default_values.json"  # change path to where you would like to save config file
    # save set attributes values to JSON configuration file
    default.save_default_to_file(file_path=config_json_path)
    # overwrite DefaultValues() instance (default) attributes values from config file
    default.read_default_from_file(file_path=config_json_path)


.. _module_blend:

Module blend
------------

TODO

Additional info is in :ref:`rvt.blend`.

Manual blending
^^^^^^^^^^^^^^^

.. code-block:: python

    layers_manual = rvt.blend.BlenderCombination()  # create class which will hold layers
    # you have two options to add layer:
    # option 1, create with method
    layers_manual.create_layer(vis_method="Sky-View Factor", normalization="value", minimum=0.7, maximum=1,
                              blend_mode="multiply", opacity=25,
                              image=svf_arr)  # automatically creates BlenderLayer() and appends it to BlenderCombination()
    # option 2, create class BlenderLayer instance and then add with method
    layer1 = rvt.blend.BlenderLayer(vis_method="Sky-View Factor", normalization="value", minimum=0.7, maximum=1,
                                    blend_mode="multiply", opacity=25,
                                    image=svf_arr)
    layers_manual.add_layer(layer1)

You can add as many layers as you need. When adding/creating layers you can define image or image_path parameter or none of them. If you define ``image_path`` (you have to save image first) and not ``image`` then blending will work faster because it will not hold all images (from all layers) in memory. It will read them simultaneously. If both ``image`` and ``image_path`` are None (not defined) then when calling method ``render_all_images()`` visualizations will be calculated automatically when needed (``vis_method`` parameter has to be correct).

.. code-block:: python

    # you can input calculated image (preferred method for non rvt visualizations)
    layers_manual.create_layer(vis_method="Sky-View Factor", normalization="value", minimum=0.7, maximum=1,
                              blend_mode="multiply", opacity=25,
                              image=svf_arr)
    # or you can input image_path
    layers_manual.create_layer(vis_method="Sky-View Factor", normalization="value", minimum=0.7, maximum=1,
                              blend_mode="multiply", opacity=25,
                              image_path=svf_path)
    # or you don't define them (None), vis_method has to be correct (rvt, suggested method)
    layers_manual.create_layer(vis_method="Sky-View Factor", normalization="value", minimum=0.7, maximum=1,
                              blend_mode="multiply", opacity=25)

After you added all the layers you would like to blend. You call method ``render_all_images()`` to create blended image. If both ``image`` and ``image_path`` are None, you can define parameters for specific visualisation function with parameter ``default``. If you call method ``add_dem_path()`` (needed for profile) and define method parameter ``save_render_path``, result will be saved in that path, else it will only return result raster array.

.. code-block:: python

    layers_manual.add_dem_path(dem_path=input_dem_path)  # needed when you wish to save render (save_render_path defined in render_all_images())
    render_arr = layers_manual.render_all_images(save_render_path=save_render_path)  # to save rendered array in save_render_path
    render_arr = layers_manual.render_all_images()  # to only get result render array (render_arr)

Automatic blending
^^^^^^^^^^^^^^^^^^

Automatic blending depends on ``rvt.default``, so you have to import ``rvt.default``.

.. code-block:: python

    import rvt.blend
    import rvt.default

Automatic blending is filling ``rvt.blender.BlenderCombination`` from file. To create example file where we can later change parameters we call function ``create_blender_file_example()``.

.. code-block:: python

    blender_file = rvt.blend.create_blender_file_example(file_path=r"settings\blender_file_example.txt")

To blend from file we also need visualization function parameters values which we define in   class ``rvt.default.DefaultValues()`` (see :ref:`module_default`).

.. code-block:: python

    default = rvt.default.DefaultValues()

To blend from file we create ``BlenderCombination()`` class, call method ``read_from_file()`` and then ``render_all_images()``. In ``render_all_images()`` method we can save (to dem_path directory) specific visualization if we set parameter ``save_visualization`` to True.

.. code-block:: python

    layers_auto = rvt.blend.BlenderCombination()
    layers_auto.read_from_file(file_path=blender_file)   # we can make our own blender_file (change example)
    layers_auto.add_dem_path(input_dem_path) # needed when save_visualizations is True, and we wish to save render (save_render_path is set)
    layers_auto.add_dem_arr(dem_arr=input_dem_arr, dem_resolution=x_res)  # needed when save_visualizations is False
    render_arr = layers_auto.render_all_images(save_visualizations=False, save_render_path=output_blend_path)

Sample dataset
--------------

A sample dataset for trying RVT python is available here in the package ``/test_data/TM1_564_146.tif``. Additional files are available here:

`RVT Demo Data <https://rebrand.ly/rvt_demo>`_

Download it, save it in ``test_data`` directory and try the visualisations.

Examples on how to use are in the following files:

.. code-block:: python

    test_vis.py
    test_blend.py
    test_default.py
    test_custom_color_scheme.py
