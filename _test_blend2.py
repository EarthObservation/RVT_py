# simple blending test

import rvt.default
import rvt.visualizations
import rvt.blender

layers = rvt.blend.BlenderCombination()
input_dem_path = r"test_data\TM1_564_146.tif"
layers.add_dem_path(dem_path=input_dem_path)
output_blend_path = r"test_data\test_blend.tif"
dict_arr_res = rvt.default.get_raster_arr(input_dem_path)
input_dem_arr = dict_arr_res["array"]
x_res = dict_arr_res["resolution"][0]
y_res = dict_arr_res["resolution"][1]

layers.add_dem_arr(input_dem_arr, x_res)

# layers;vis_method;norm;min;max;blending_mode;opacity
# 1;multi-hillshade;value;0;1;normal;50
layers.create_layer(vis_method="Multiple directions hillshade", normalization="value", minimum=0, maximum=1,
                    blend_mode="normal", opacity=50)

# 2;svf;value;0.7;1.0;multiply;100
layers.create_layer(vis_method="Sky-View Factor", normalization="value", minimum=0.7, maximum=1,
                    blend_mode="normal", opacity=100)  # you could also input image_path

layers.render_all_images(save_render_path=output_blend_path)
