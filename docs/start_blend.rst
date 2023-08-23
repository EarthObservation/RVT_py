.. _start_blend:

Blending visualizations
=======================

You can blend manually or automatically.

**Manual blending** allows you to use visualizations that are not part of ``rvt``. When blending manually you have to define each layer (visualization) in Python.

**Automatic blending** automatically computes ``rvt`` visualizations and blends them together according to a configuration ``JSON`` file, which can be edited.

The main class of the ``rvt.blend`` module for blending is ``BlenderCombination``, which has the list attribute ``layers`` where instances of class ``BlenderLayer`` are stored. In ``BlenderLayer`` instances in ``layers`` we store a specific visualization and its parameters for blending.

The ``BlenderCombination`` class has the method ``render_all_images()``, which blends together all ``BlenderLayer`` instances (visualizations) in the ``BlenderCombination.layers`` list and outputs the blended image.

You can blend as many layers as you want.

----

Manual blending
---------------

**Example**

Let's say we have already calculated the simple local relief model (slrm_arr), slope (slope_arr) and hillshade (hillshade_arr), and now need to blend all the calculated visualizations together:

.. code-block:: python
    
    # import the module
    import rvt.blend

    # create the BlenderCombination() class instance which will hold the layers (visualizations)
    combination_manual = rvt.blend.BlenderCombination()

    # call BlenderCombination.create_layer() to add a layer
    # this creates a BlenderLayer instance and adds it to BlenderCombination.layers

    # 1st layer
    # add slrm layer with 2% perc cuttoff on both sides, multiply blend mode and 25% opacity
    combination_manual.create_layer(
        vis_method="Simple local relief model", 
        normalization="perc", 
        minimum=2, 
        maximum=2,
        blend_mode="multiply", 
        opacity=25, 
        image=slrm_arr
        )
                              
    # 2nd layer
    # add slope layer with value stretch from 0 to 51, luminosity blend mode and 50% opacity
    combination_manual.create_layer(
        vis_method="Slope gradient", 
        normalization="value", 
        minimum=0, 
        maximum=51,
        blend_mode="luminosity", 
        opacity=50, 
        image=slope_arr
        )
                              
    # 3rd layer
    # add hillshade layer with value stretch from 0 to 1, normal blend mode and 100% opacity
    combination_manual.create_layer(
        vis_method="Hillshade", 
        normalization="value", 
        minimum=0, 
        maximum=1,
        blend_mode="normal", 
        opacity=100, 
        image=hillshade_arr
        )

    # if we want to save the blended image to a file, we need to add dem_path to the 
    # combination (for metadata, geodata)
    combination_manual.add_dem_path(dem_path=input_dem_path)

    # blend them all together
    # you can save the blend to GeoTIFF if save_render_path presented 
    # (and dem_path is added), otherwise it only returns array
    render_arr = combination_manual.render_all_images(save_render_path=output_blend_path)

**Example**

You can also let the ``BlenderCombination`` class automatically compute the visualization or give the path to a visualization. 

If you **don't** provide the **image** parameter, and the vis_method parameter is correct (an existing ``rvt.vis`` function), blender automatically calculates the visualization. 

If you **don't** provide the **image** parameter, but **do** provide the **image_path** parameter (if you provide both image will be used), blender will read the visualization from image_path.

If you **don't** provide the **image and image_path** parameters, you have to add an ``rvt.default.DefaultValues`` instance as a parameter to ``BlenderCombination.render_all_images()``. Blender then takes the parameters set in this class when calculating specific visualizations.
You also have to add dem array and its resolution. 

The example below uses all three methods:

.. code-block:: python

    # import all required modules
    import rvt.blend
    import rvt.default

    # create the BlenderCombination() class instance which will hold the layers (visualizations)
    combination_manual = rvt.blend.BlenderCombination()

    # we will let blender compute the slrm visualization. so, we need to create 
    # rvt.default.DefaultValues() and change the parameters for slrm. we will later 
    # add default to the combination_manual.render_all_images() method
    default = rvt.default.DefaultValues()
    default.slrm_rad_cell = 15

    # 1st layer
    # add slrm layer with 2% perc cuttoff on both sides, multiply blend mode and 25% opacity
    # image and image_path parameters both not provided, so slrm is calculated automatically
    combination_manual.create_layer(
        vis_method="Simple local relief model",
        normalization="perc", 
        minimum=2, 
        maximum=2,
        blend_mode="multiply", 
        opacity=25
        )
                              
    # 2nd layer
    # add slope layer with value stretch from 0 to 51, luminosity blend mode and 50% opacity
    # image_path parameter provided to slope, so slope is read from file
    combination_manual.create_layer(
        vis_method="Slope gradient", 
        normalization="value", 
        minimum=0, 
        maximum=51,
        blend_mode="luminosity", 
        opacity=50, 
        image_path=slope_path
        )
                              
    # 3rd layer
    # add hillshade layer with value stretch from 0 to 1, normal blend mode and 100% opacity
    # image parameter provided
    combination_manual.create_layer(
        vis_method="Hillshade", 
        normalization="value", 
        minimum=0, 
        maximum=1,
        blend_mode="normal",
        opacity=100,
        image=hillshade_arr
        )

    # we have to add dem array and resolution so that slrm can be computed
    combination_manual.add_dem_arr(dem_arr=input_dem_arr, dem_resolution=resolution)

    # blend them all together and add default where slrm parameters are defined
    render_arr = combination_manual.render_all_images(default=default)

----

Automatic blending
------------------

Automatic blending is blending from a configuration ``JSON`` file. You can create a ``JSON`` file and change it to suit your needs.

**Example**

.. code-block:: python

    # import the module
    import rvt.blend

    # create the BlenderCombination() class
    combination_auto = rvt.blend.BlenderCombination()
    
    # to create the JSON blender combination configuration file example, change the 
    # path to where you wish to save the file
    blender_combination_path = r"settings\blender_file_example.txt"
    rvt.blend.create_blender_file_example(file_path=blender_combination_path)

    # set the parameters of the visualizations you will be using
    default = rvt.default.DefaultValues()
    # for example default.hs_sun_el=40

    # read the JSON combination configuration file
    combination_auto.read_from_json_file(file_path=blender_combination_path)

    # needed when save_visualizations is True and save_render_path is not None
    layers_auto.add_dem_path(input_dem_path)

    # call the method render_all_images() and its parameters
    # we can save a specific visualization (to dem_path directory) if we set the 
    # parameter ``save_visualization`` to True
    layers_auto.render_all_images(
        default=default, 
        save_visualizations=True, 
        save_render_path=output_blend_path,
        save_float=True, 
        save_8bit=True # set save_8bit=True if you also wish to save an 8bit version
        )  

.. seealso:: Find out more about blending in :ref:`rvt.blend`.
