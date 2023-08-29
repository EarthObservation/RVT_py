import rvt.default
import rvt.blender_functions
import rvt.blender
import rvt.visualizations
import numpy as np


# Test color relief image map (CRIM)
dem_path = r"D:\RVT_py\test\CRIM_TEST\Skolj_dem_05m.tif"
dem_dict = rvt.default.get_raster_arr(dem_path)
dem_arr = dem_dict["array"]
dem_res = dem_dict["resolution"][0]

output_srrim_path = "{}_{}.tif".format(dem_path.rstrip(".tif"), "crim_reds")
rendered_image = rvt.blend.color_relief_image_map(dem=dem_arr, resolution=dem_res, colormap="Reds_r",
                                                  min_colormap_cut=0.5, max_colormap_cut=1)
rvt.default.save_raster(src_raster_path=dem_path, out_raster_arr=rendered_image, out_raster_path=output_srrim_path,
                        e_type=6)

output_srrim_path = "{}_{}.tif".format(dem_path.rstrip(".tif"), "crim_blues")
rendered_image = rvt.blend.color_relief_image_map(dem=dem_arr, resolution=dem_res, colormap="Blues_r",
                                                  min_colormap_cut=0.5, max_colormap_cut=1)
rvt.default.save_raster(src_raster_path=dem_path, out_raster_arr=rendered_image, out_raster_path=output_srrim_path,
                        e_type=6)

output_srrim_path = "{}_{}.tif".format(dem_path.rstrip(".tif"), "crim_purples")
rendered_image = rvt.blend.color_relief_image_map(dem=dem_arr, resolution=dem_res, colormap="Purples_r",
                                                  min_colormap_cut=0.5, max_colormap_cut=1)
rvt.default.save_raster(src_raster_path=dem_path, out_raster_arr=rendered_image, out_raster_path=output_srrim_path,
                        e_type=6)

output_srrim_path = "{}_{}.tif".format(dem_path.rstrip(".tif"), "crim_oranges")
rendered_image = rvt.blend.color_relief_image_map(dem=dem_arr, resolution=dem_res, colormap="Oranges_r",
                                                  min_colormap_cut=0.5, max_colormap_cut=1)
rvt.default.save_raster(src_raster_path=dem_path, out_raster_arr=rendered_image, out_raster_path=output_srrim_path,
                        e_type=6)

output_srrim_path = "{}_{}.tif".format(dem_path.rstrip(".tif"), "crim_YlOrBr")
rendered_image = rvt.blend.color_relief_image_map(dem=dem_arr, resolution=dem_res, colormap="YlOrBr_r")
rvt.default.save_raster(src_raster_path=dem_path, out_raster_arr=rendered_image, out_raster_path=output_srrim_path,
                        e_type=6)

output_srrim_path = "{}_{}.tif".format(dem_path.rstrip(".tif"), "crim_coolwarm")
rendered_image = rvt.blend.color_relief_image_map(dem=dem_arr, resolution=dem_res, colormap="coolwarm_r",
                                                  min_colormap_cut=0.2, max_colormap_cut=0.8)
rvt.default.save_raster(src_raster_path=dem_path, out_raster_arr=rendered_image, out_raster_path=output_srrim_path,
                        e_type=6)

output_srrim_path = "{}_{}.tif".format(dem_path.rstrip(".tif"), "crim_hot")
rendered_image = rvt.blend.color_relief_image_map(dem=dem_arr, resolution=dem_res, colormap="hot_r",
                                                  min_colormap_cut=0.4, max_colormap_cut=0.6)
rvt.default.save_raster(src_raster_path=dem_path, out_raster_arr=rendered_image, out_raster_path=output_srrim_path,
                        e_type=6)

output_srrim_path = "{}_{}.tif".format(dem_path.rstrip(".tif"), "crim_RdYlBu")
rendered_image = rvt.blend.color_relief_image_map(dem=dem_arr, resolution=dem_res, colormap="RdYlBu",
                                                  min_colormap_cut=0.1, max_colormap_cut=0.9)
rvt.default.save_raster(src_raster_path=dem_path, out_raster_arr=rendered_image, out_raster_path=output_srrim_path,
                        e_type=6)
