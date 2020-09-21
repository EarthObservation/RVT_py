# rvt.blend quick TEST


import rasterio as rio
import rvt.vis
import rvt.blend
import numpy as np

# read dem
input_dem_path = r"D:\RVT_py\test\TM1_564_146.tif"
output_blend_path = r"D:\RVT_py\test\TM1_564_146_test_blend.tif"
input_dem_dataset = rio.open(input_dem_path)
t = input_dem_dataset.transform
x_res = t[0]
y_res = -t[4]
input_dem_arr = input_dem_dataset.read()[0]

# test blend combination archeological (VAT), general

layers = rvt.blend.BlenderLayers()

# layers;vis_method;norm;min;max;blending_mode;opacity
# 1;svf;value;0.7;1.0;multiply;25
svf_dict = rvt.vis.sky_view_factor(dem=input_dem_arr, resolution=x_res, compute_svf=True, compute_opns=True)
svf_arr = svf_dict["svf"]
layers.create_layer(vis_method="svf", normalization="value", minimum=0.7, maximum=1, blend_mode="multiply", opacity=25,
                    image=np.array(svf_arr))

# 2;opns_pos;value;68;93;overlay;50
opns_arr = svf_dict["opns"]
layers.create_layer(vis_method="opns_pos", normalization="value", minimum=68, maximum=93, blend_mode="overlay",
                    opacity=50, image=np.array(opns_arr))

# 3;slope;value;0;50;luminosity;50
slope_dict = rvt.vis.slope_aspect(dem=input_dem_arr, resolution_x=x_res, resolution_y=y_res, output_units="degree",
                                  ve_factor=1)
slope_arr = slope_dict["slope"]
layers.create_layer(vis_method="slope", normalization="value", minimum=0, maximum=50, blend_mode="normal",
                    opacity=50, image=np.array(slope_arr))

# 4;hillshade;value;0;1;normal;100
hillshade_arr = rvt.vis.hillshade(dem=input_dem_arr, resolution_x=x_res, resolution_y=y_res)
layers.create_layer(vis_method="hillshade", normalization="value", minimum=0, maximum=1, blend_mode="normal",
                    opacity=100, image=np.array(hillshade_arr))

# 5;None
layers.create_layer(vis_method=None)


layers.normalize_images_on_layers()
rendered_imgs = layers.render_all_images()
rendered_imgs = rendered_imgs.astype('float64')
profile = input_dem_dataset.profile
profile.update(dtype='float64')
output_blend_arr_dataset = rio.open(output_blend_path, "w", **profile)
output_blend_arr_dataset.write(np.array([rendered_imgs]))
