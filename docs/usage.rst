.. _usage:

Usage
=====

.. _module_vis:

Module vis
----------

For reading raster data (DEM) from files (GeoTIFF) we suggest to use ``rasterio``.

.. code-block:: python

    import rasterio as rio
    
    input_dem_path = "test_data/TM1_564_146.tif"  # path to your DEM
    input_dem_dataset = rio.open(input_dem_path)  # open with rasterio
    t = input_dem_dataset.transform
    x_res = t[0]  # to get x resolution
    y_res = -t[4] # to get y resolution
    input_dem_arr = input_dem_dataset.read()[0]  # to get 2D numpy array of DEM
    input_dem_dataset.close()


Importing ``rvt.vis`` module.

.. code-block:: python

    import rvt.vis


Additional info is in :ref:`rvt.vis`.

.. _module_default:

Module default
--------------

Default module contains class ``DefaultValues()`` where we can store our visualization functions parameters, this class also has methods for saving and computing visualization functions with that parameters. To import it we use:

.. code-block:: python

    import rvt.default


We have to create class instance first.

.. code-block:: python

    default = rvt.default.DefaultValues()

When we create instance parameter values are already populated with default values. We can store them in file, change them in file and read them back from file. We can also change them in program.

Additional info is in :ref:`rvt.default`.

.. _module_blend:

Module blend
------------

Blend is module for blending different visualizations together. To import it we use:

.. code-block:: python

    import rvt.blend

We could use manual blending (we compute visualisations) or we could use automatic blending from file (automatically computed visualisations with ``rvt.vis`` and stored in location where DEM is).

Manual blending depends on module ``rvt.default`` (see :ref:`module_default`) and ``rvt.vis``. To start blending we first need to create instance of class ``BlenderLayers()`` which contains list of layers (``BlenderLayer()``). Single layer is defined in ``BlenderLayer()`` class.

Additional info is in :ref:`rvt.blend`.

Manual blending
^^^^^^^^^^^^^^^

.. code-block:: python

    layers_manual = rvt.blend.BlenderLayers()  # create class which will hold layers
    # you have two options to add layer:
    # option 1, create with method
    layers_manual.create_layer(vis_method="Sky-View Factor", normalization="value", minimum=0.7, maximum=1,
                              blend_mode="multiply", opacity=25,
                              image=svf_arr)  # automatically creates BlenderLayer() and appends it to BlenderLayers()
    # option 2, create class BlenderLayer instance and then add with method
    layer1 = rvt.blend.BlenderLayer(vis_method="Sky-View Factor", normalization="value", minimum=0.7, maximum=1,
                                    blend_mode="multiply", opacity=25,
                                    image=svf_arr)
    layers_manual.add_layer(layer1)

You can add as many layers as you need. When adding/creating layers you can define image or image_path parameter. If you define ``image_path`` (you have to save image first) and not image then blending will work faster because it will not hold all images (from all layers) in memory. It will read them simultaneously.

.. code-block:: python

    layers_manual.create_layer(vis_method="Sky-View Factor", normalization="value", minimum=0.7, maximum=1,
                              blend_mode="multiply", opacity=25,
                              image_path=svf_path)  # image_path instead of image

After you added all the layers you would like to blend. You call method ``render_all()``. If you define method parameters ``input_dem_path`` (needed for profile info) and ``output_blend_path``. Result will be saved in ``output_blend_path`` else it will only return result raster array.

.. code-block:: python

    render_arr = layers_manual.render_all_images(input_dem_path, output_blend_path)  # to save rendered array in output_blend_path
    render_arr = layers_manual.render_all_images()  # to only get result render array (render_arr)

Automatic blending
^^^^^^^^^^^^^^^^^^

Automatic blending depends on ``rvt.default``, so you have to import ``rvt.default``.

.. code-block:: python

    import rvt.blend
    import rvt.default

Automatic blending is blending from file. All the visualisations we need for blending will be computed and stored in the directory where input raster is. To create example file where we can later change parameters we call function ``create_blender_file_example``.

.. code-block:: python

    blender_file = rvt.blend.create_blender_file_example(file_path=r"settings\blender_file_example.txt")

To blend from file we also need function parameters values which we define in   class ``rvt.default.DefaultValues()`` (see :ref:`module_default`).

.. code-block:: python

    default = rvt.default.DefaultValues()

To blend from file we create ``BlenderLayers()`` class, call method ``build_blender_layers_from_file`` and then ``render_all_images``.

.. code-block:: python

    layers_auto = rvt.blend.BlenderLayers()
    layers_auto.build_blender_layers_from_file(dem_path=input_dem_path, file_path=blender_file, default=default)  # we can make our own blender_file (change example), we can change default
    render_arr = layers_auto.render_all_images(dem_path=input_dem_path, save_render_path=output_blend_path)  # if we dont input dem_path and save_rander_path result won't be stored as raster (tif)

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
