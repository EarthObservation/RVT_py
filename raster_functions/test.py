# RVT_vis_fn TEST

import rasterio as rio
import RVT_vis_fn
import numpy as np


def test_slope_aspect(input_DEM_path, resolution, ve_factor, output_units, output_path):
    input_DEM_dataset = rio.open(input_DEM_path)
    input_DEM_arr = input_DEM_dataset.read()[0]

    slope_arr, aspect_arr = RVT_vis_fn.slope_aspect(input_DEM_arr=input_DEM_arr, resolution_x=resolution,
                                                    resolution_y=resolution, ve_factor=ve_factor,
                                                    is_padding_applied=False, output_units=output_units)
    output_slope_dataset = rio.open(output_path, "w", **input_DEM_dataset.meta)
    output_slope_dataset.write(np.array([slope_arr]))


def test_analytical_hillshading(input_DEM_path, output_path):
    input_DEM_dataset = rio.open(input_DEM_path)
    t = input_DEM_dataset.transform
    x_res = t[0]
    y_res = -t[4]
    input_DEM_arr = input_DEM_dataset.read()[0]

    hillshading_arr = RVT_vis_fn.analytical_hillshading(input_DEM_arr=input_DEM_arr, resolution_x=x_res,
                                                        resolution_y=y_res, sun_azimuth=315, sun_elevation=35,
                                                        is_padding_applied=False, bytscl=False)
    output_hillshading_dataset = rio.open(output_path, "w", **input_DEM_dataset.meta)
    output_hillshading_dataset.write(np.array([hillshading_arr]))

def test_multiple_directions_hillshading(input_DEM_path, output_path):
    input_DEM_dataset = rio.open(input_DEM_path)
    t = input_DEM_dataset.transform
    x_res = t[0]
    y_res = -t[4]
    input_DEM_arr = input_DEM_dataset.read()[0]
    multi_hillshading_arr = RVT_vis_fn.multiple_directions_hillshading(input_DEM_arr=input_DEM_arr, resolution_x=x_res,
                                                        resolution_y=y_res, nr_directions=16, sun_elevation=35,
                                                        is_padding_applied=False)

    profile = input_DEM_dataset.profile
    profile.update(count=16)
    profile.update(dtype='uint8')
    output_multi_hillshading_dataset = rio.open(output_path, "w", **profile)
    output_multi_hillshading_dataset.write(multi_hillshading_arr)


test_slope_aspect(input_DEM_path=r"D:\RVT_py\test\TM1_564_146.tif", resolution=1, ve_factor=1, output_units="degree",
                  output_path=r"D:\RVT_py\test\TM1_564_146_test_slope.tif")
test_analytical_hillshading(input_DEM_path=r"D:\RVT_py\test\TM1_564_146.tif",
                            output_path=r"D:\RVT_py\test\TM1_564_146_test_hillsahade.tif")
test_multiple_directions_hillshading(input_DEM_path=r"D:\RVT_py\test\TM1_564_146.tif",
                            output_path=r"D:\RVT_py\test\TM1_564_146_test_multi_hillsahade.tif")
