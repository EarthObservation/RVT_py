# rvt.blend quick TEST

import rvt.default
import rvt.visualizations
import rvt.blender
import numpy as np

# test blend combination archeological (VAT), general

#####
# manual blending, custom raster numpy arrays

# if you create_layer and don't input image or image_path then if vis_method is correct it automatically
# calculates visualization in render_all_images

layers_manual = rvt.blend.BlenderCombination()
input_dem_path = r"test_data\TM1_564_146.tif"
layers_manual.add_dem_path(dem_path=input_dem_path)
output_blend_path = r"test_data\TM1_564_146_test_blend_manual.tif"
dict_arr_res = rvt.default.get_raster_arr(input_dem_path)
input_dem_arr = dict_arr_res["array"]
x_res = dict_arr_res["resolution"][0]
y_res = dict_arr_res["resolution"][1]

# layers;vis_method;norm;min;max;blending_mode;opacity
# 1;svf;value;0.7;1.0;multiply;25
svf_dict = rvt.vis.sky_view_factor(dem=input_dem_arr, resolution=x_res, compute_svf=True, compute_opns=True)
svf_arr = svf_dict["svf"]
layers_manual.create_layer(vis_method="Sky-View Factor", normalization="value", minimum=0.7, maximum=1,
                           blend_mode="multiply", opacity=25,
                           image=svf_arr)  # you could also input image_path
# 2;opns_pos;value;68;93;overlay;50
opns_arr = svf_dict["opns"]
layers_manual.create_layer(vis_method="Openness - Positive", normalization="value", minimum=68, maximum=93,
                           blend_mode="overlay",
                           opacity=50, image=opns_arr)
# 3;slope;value;0;50;luminosity;50
slope_dict = rvt.vis.slope_aspect(dem=input_dem_arr, resolution_x=x_res, resolution_y=y_res, output_units="degree",
                                  ve_factor=1)
slope_arr = slope_dict["slope"]
layers_manual.create_layer(vis_method="Slope gradient", normalization="value", minimum=0, maximum=50,
                           blend_mode="luminosity",
                           opacity=50, image=slope_arr)
# # 4;hillshade;value;0;1;normal;100
hillshade_arr = rvt.vis.hillshade(dem=input_dem_arr, resolution_x=x_res, resolution_y=y_res)
layers_manual.create_layer(vis_method="Hillshade", normalization="value", minimum=0, maximum=1, blend_mode="normal",
                           opacity=100, image=hillshade_arr)
# 5;None
layers_manual.create_layer(vis_method=None)
# you can save to GeoTif if save_render_path presented else it only returns array
render_arr = layers_manual.render_all_images(save_render_path=output_blend_path)

# you can save layers combination to .json file, be aware image and image_path won't be saved
# this is a problem when vis_method is non rvt visualization(is not correct)!
layers_manual.save_to_json_file(r"settings\blender_custom_layers.json")

#####

#####
# automatic blending, blending from blender_file with values from default.DefaultValues class
# when save_visualizations=False, blending save every needed visualization in GeoTif in dem_path directory
input_dem_path = r"test_data\TM1_564_146.tif"
# Example file (for file_path) in dir settings: blender_file_example.txt
blender_file = r"settings\blender_file_example.json"
output_blend_path = r"test_data\TM1_564_146_test_blend_automatic.tif"
layers_auto = rvt.blend.BlenderCombination()
default = rvt.default.DefaultValues()
default.read_default_from_file(r"settings\default_settings.json")
layers_auto.read_from_json_file(file_path=blender_file)  # build BlenderCombination from file
# when building_blender from file single BlenderLayer image and image_path are None
layers_auto.add_dem_path(input_dem_path)  # needed when save_visualizations is True and save_rander_path is not None
# render_all_images reads images simultaneously if layer (BlenderLayer) image is None and image_path is None it
# calculates them
layers_auto.render_all_images(default=default, save_visualizations=True, save_render_path=output_blend_path,
                              save_float=True, save_8bit=True)
#####

#####
# automatic blending, blending from blender_file with values from default.DefaultValues class
# when save_visualizations=False, blending doesn't save every visualization, it calculates it when needed
layers_auto.add_dem_arr(dem_arr=input_dem_arr, dem_resolution=x_res)  # needed when save_visualizations is False
rendered_arr = layers_auto.render_all_images(save_visualizations=False)
#####
