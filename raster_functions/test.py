# RVT_vis_fn TEST
import time
import rasterio as rio
import RVT_vis_fn
import numpy as np
import scipy.ndimage


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
                                                                       resolution_y=y_res, nr_directions=16,
                                                                       sun_elevation=35,
                                                                       is_padding_applied=False)

    profile = input_DEM_dataset.profile
    profile.update(count=16)
    profile.update(dtype='uint8')
    output_multi_hillshading_dataset = rio.open(output_path, "w", **profile)
    output_multi_hillshading_dataset.write(multi_hillshading_arr)


def test_SLRM(input_DEM_path, output_path):
    input_DEM_dataset = rio.open(input_DEM_path)
    input_DEM_arr = input_DEM_dataset.read()[0]

    slrm_arr = RVT_vis_fn.SLRM(input_DEM_arr, bytscl=False)
    profile = input_DEM_dataset.profile
    # profile.update(dtype='uint8')
    output_slrm_arr_dataset = rio.open(output_path, "w", **profile)
    output_slrm_arr_dataset.write(np.array([slrm_arr]))


def test_sky_view_factor(input_DEM_path, output_svf_path, output_asvf_path, output_opns_path):
    input_DEM_dataset = rio.open(input_DEM_path)
    t = input_DEM_dataset.transform
    x_res = t[0]
    y_res = -t[4]
    input_DEM_arr = input_DEM_dataset.read()[0]

    svf_arr, asvf_arr, opns_arr = RVT_vis_fn.sky_view_factor(input_DEM_arr=input_DEM_arr, resolution=x_res,
                                                             bytscl=False)

    profile = input_DEM_dataset.profile
    profile.update(dtype='float64')
    output_svf = rio.open(output_svf_path, "w", **profile)
    output_svf.write(np.array([svf_arr]))

    output_asvf = rio.open(output_asvf_path, "w", **profile)
    output_asvf.write(np.array([asvf_arr]))

    output_opns = rio.open(output_opns_path, "w", **profile)
    output_opns.write(np.array([opns_arr]))


def test_sky_illumination(input_DEM_path, output_path):
    input_DEM_dataset = rio.open(input_DEM_path)
    t = input_DEM_dataset.transform
    x_res = t[0]
    y_res = -t[4]
    input_DEM_arr = input_DEM_dataset.read()[0]
    start_time = time.time()
    skyillumination_arr = RVT_vis_fn.sky_illumination(input_DEM_arr=input_DEM_arr, resolution=x_res)
    end_time = time.time()
    print(end_time - start_time)
    profile = input_DEM_dataset.profile
    profile.update(dtype='float64')
    output_skyilum = rio.open(output_path, "w", **profile)
    output_skyilum.write(np.array([skyillumination_arr]))


def test_local_dominance(input_DEM_path, output_path):
    input_DEM_dataset = rio.open(input_DEM_path)
    input_DEM_arr = input_DEM_dataset.read()[0]
    localdominance_arr = RVT_vis_fn.local_dominance(input_DEM_arr)

    profile = input_DEM_dataset.profile
    #profile.update(dtype='float64')
    output_localdom = rio.open(output_path, "w", **profile)
    output_localdom.write(np.array([localdominance_arr]))


# test_slope_aspect(input_DEM_path=r"D:\RVT_py\test\TM1_564_146.tif", resolution=1, ve_factor=1, output_units="degree",
#                   output_path=r"D:\RVT_py\test\TM1_564_146_test_slope.tif")
# test_analytical_hillshading(input_DEM_path=r"D:\RVT_py\test\TM1_564_146.tif",
#                             output_path=r"D:\RVT_py\test\TM1_564_146_test_hillsahade.tif")
# test_multiple_directions_hillshading(input_DEM_path=r"D:\RVT_py\test\TM1_564_146.tif",
#                             output_path=r"D:\RVT_py\test\TM1_564_146_test_multi_hillsahade.tif")
#
# test_SLRM(input_DEM_path=r"D:\RVT_py\test\TM1_564_146.tif",
#                   output_path=r"D:\RVT_py\test\TM1_564_146_test_SLRM.tif")

# test_sky_view_factor(input_DEM_path=r"D:\RVT_py\test\TM1_564_146.tif",
#                      output_svf_path=r"D:\RVT_py\test\TM1_564_146_test_SVF.tif",
#                      output_asvf_path=r"D:\RVT_py\test\TM1_564_146_test_ASVF.tif",
#                      output_opns_path=r"D:\RVT_py\test\TM1_564_146_test_OPNS.tif")

test_sky_illumination(input_DEM_path=r"D:\RVT_py\test\TM1_564_146.tif",
                      output_path=r"D:\RVT_py\test\TM1_564_146_test_skyillum.tif")

# test_local_dominance(input_DEM_path=r"D:\RVT_py\test\TM1_564_146.tif",
#                       output_path=r"D:\RVT_py\test\TM1_564_146_test_localdom.tif")