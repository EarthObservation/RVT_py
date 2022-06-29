.. _module_blend:

Blending
========

You can blend manually or automatically.  

Manually blending allows you to use visualizations that are not part of ``rvt``. When blending manually you have to define each layer (visualization) in python.

Automatically blending automatically computes ``rvt`` visualizations and blends them together according to the configuration ``JSON`` file, which can be edited.

The main class of the ``rvt.blend`` module for blending is ``BlenderCombination``, which has the list attribute ``layers`` where instances of class ``BlenderLayer`` are stored. In ``BlenderLayer`` instances in ``layers`` we store a specific visualization and its parameters for blending.

The ``BlenderCombination`` class has the method ``render_all_images()``, which blends together all ``BlenderLayer`` instances (visualizations) in the ``BlenderCombination.layers`` list and outputs the blended image.

Manual blending
---------------

Import the ``rvt.blend`` module and create a ``BlenderCombination`` instance. 

To add a layer (visualization) with parameters to a combination, you can call ``BlenderCombination.create_layer()``. This creates a ``BlenderLayer`` instance and adds it to ``BlenderCombination.layers``.

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

If you don't provide the image parameter, and the vis_method parameter is correct (existing rvt visualization function), blender automatically calculates the visualization. 

If you provide the image_path parameter and not the image parameter (if you provide both image will be used), blender will read the visualization from image_path. 

If you don't input the image and image_path parameters, you have to add an ``rvt.default.DefaultValues`` instance as a parameter to ``BlenderCombination.render_all_images()``. Blender then takes the parameters set in this class when calculating specific visualizations.
You also have to add dem array and its resolution. 

See example below which uses all three methods:

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

You can add as many layers as you want.

Automatic blending
------------------

Automatic blending is blending from a configuration ``JSON`` file. You can create an example file and change it to suit your needs.

To blend from a file, create the ``BlenderCombination()`` class, call the method ``read_from_file()`` and then ``render_all_images()``. In the ``render_all_images()`` method we can save a specific visualization (to dem_path directory) if we set the parameter ``save_visualization`` to True.

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
                              
----

Find out more about blending in :ref:`rvt.blend`.
