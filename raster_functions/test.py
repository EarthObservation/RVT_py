# RVT_vis_fn TEST

import rasterio as rio
import RVT_vis_fn

def test_slope_aspect(input_DEM_path, resolution, ve_factor, output_units, output_path):
    input_DEM_dataset = rio.open(input_DEM_path)
    input_DEM_arr = input_DEM_dataset.read()

    slope_arr, aspect_arr = RVT_vis_fn.slope_aspect(input_DEM_arr, resolution, ve_factor, output_units)

    output_slope_dataset = rio.open(output_path, "w", **input_DEM_dataset.meta)
    output_slope_dataset.profile["dtype"] = "uint32"
    output_slope_dataset.write(slope_arr)

def test_analytical_hillshading(input_DEM_path, resolution, output_path, sun_azimuth=315, sun_elevation=35):
    input_DEM_dataset = rio.open(input_DEM_path)
    input_DEM_arr = input_DEM_dataset.read()

    hillshading_arr = RVT_vis_fn.analytical_hillshading(input_DEM_arr, resolution, sun_azimuth, sun_elevation)

    output_hillshading_dataset = rio.open(output_path, "w", **input_DEM_dataset.meta)
    output_hillshading_dataset.profile["dtype"] = "uint32"
    output_hillshading_dataset.write(hillshading_arr)

#test_slope_aspect(input_DEM_path=r"D:\RVT_py\test\TM1_564_146.tif", resolution=1 , ve_factor=1, output_units="degree",
#                  output_path=r"D:\RVT_py\test\TM1_564_146_test_slope.tif")
#test_analytical_hillshading(input_DEM_path=r"D:\RVT_py\test\TM1_564_146.tif", resolution=1,
#                            output_path=r"D:\RVT_py\test\TM1_564_146_test_hillsahade.tif")