# rvt_vis quick TEST

import time
import rasterio as rio
import rvt_vis
import numpy as np


def test_slope_aspect(input_dem_path, output_path):
    input_dem_dataset = rio.open(input_dem_path)
    t = input_dem_dataset.transform
    x_res = t[0]
    y_res = -t[4]
    input_dem_arr = input_dem_dataset.read()[0]

    dict_slp_asp = rvt_vis.slope_aspect(dem=input_dem_arr, resolution_x=x_res,
                                        resolution_y=y_res, ve_factor=1,
                                        is_padding_applied=False, output_units="degree")
    slope_arr = dict_slp_asp["slope"]
    slope_arr = slope_arr.astype('float64')
    profile = input_dem_dataset.profile
    profile.update(dtype='float64')
    output_slope_dataset = rio.open(output_path, "w", **profile)
    output_slope_dataset.write(np.array([slope_arr]))


def test_hillshade(input_dem_path, output_path):
    input_dem_dataset = rio.open(input_dem_path)
    t = input_dem_dataset.transform
    x_res = t[0]
    y_res = -t[4]
    input_dem_arr = input_dem_dataset.read()[0]

    hillshading_arr = rvt_vis.hillshade(dem=input_dem_arr, resolution_x=x_res,
                                        resolution_y=y_res, sun_azimuth=315, sun_elevation=35,
                                        is_padding_applied=False)
    hillshading_arr = hillshading_arr.astype('float64')
    profile = input_dem_dataset.profile
    profile.update(dtype='float64')
    output_hillshading_dataset = rio.open(output_path, "w", **profile)
    output_hillshading_dataset.write(np.array([hillshading_arr]))


def test_multi_hillshade(input_dem_path, output_path):
    input_dem_dataset = rio.open(input_dem_path)
    t = input_dem_dataset.transform
    x_res = t[0]
    y_res = -t[4]
    input_dem_arr = input_dem_dataset.read()[0]
    multi_hillshading_arr = rvt_vis.multi_hillshade(dem=input_dem_arr, resolution_x=x_res,
                                                    resolution_y=y_res, nr_directions=16,
                                                    sun_elevation=35,
                                                    is_padding_applied=False)
    multi_hillshading_arr = multi_hillshading_arr.astype('float64')
    profile = input_dem_dataset.profile
    profile.update(dtype='float64')
    profile.update(count=16)
    output_multi_hillshading_dataset = rio.open(output_path, "w", **profile)
    output_multi_hillshading_dataset.write(multi_hillshading_arr)


def test_slrm(input_dem_path, output_path):
    input_dem_dataset = rio.open(input_dem_path)
    input_dem_arr = input_dem_dataset.read()[0]

    slrm_arr = rvt_vis.slrm(input_dem_arr)
    slrm_arr = slrm_arr.astype('float64')
    profile = input_dem_dataset.profile
    profile.update(dtype='float64')
    output_slrm_arr_dataset = rio.open(output_path, "w", **profile)
    output_slrm_arr_dataset.write(np.array([slrm_arr]))


def test_sky_view_factor(input_dem_path, output_svf_path, output_asvf_path, output_opns_path):
    input_dem_dataset = rio.open(input_dem_path)
    t = input_dem_dataset.transform
    x_res = t[0]
    y_res = -t[4]
    input_dem_arr = input_dem_dataset.read()[0]

    dict_svf_asvf_opns = rvt_vis.sky_view_factor(dem=input_dem_arr, resolution=x_res, compute_svf=True,
                                                 compute_asvf=True, compute_opns=True)
    svf_arr = dict_svf_asvf_opns["svf"]
    svf_arr = svf_arr.astype('float64')

    profile = input_dem_dataset.profile
    profile.update(dtype='float64')
    output_svf = rio.open(output_svf_path, "w", **profile)
    output_svf.write(np.array([svf_arr]))

    asvf_arr = dict_svf_asvf_opns["asvf"]
    asvf_arr = asvf_arr.astype('float64')

    output_asvf = rio.open(output_asvf_path, "w", **profile)
    output_asvf.write(np.array([asvf_arr]))

    opns_arr = dict_svf_asvf_opns["opns"]
    opns_arr = opns_arr.astype('float64')

    output_opns = rio.open(output_opns_path, "w", **profile)
    output_opns.write(np.array([opns_arr]))


def test_sky_illumination(input_dem_path, output_path):
    input_dem_dataset = rio.open(input_dem_path)
    t = input_dem_dataset.transform
    x_res = t[0]
    y_res = -t[4]
    input_dem_arr = input_dem_dataset.read()[0]
    skyillumination_arr = rvt_vis.sky_illumination(dem=input_dem_arr, resolution=x_res)
    skyillumination_arr = skyillumination_arr.astype('float64')
    profile = input_dem_dataset.profile
    profile.update(dtype='float64')
    output_skyilum = rio.open(output_path, "w", **profile)
    output_skyilum.write(np.array([skyillumination_arr]))


def test_local_dominance(input_dem_path, output_path):
    input_dem_dataset = rio.open(input_dem_path)
    input_dem_arr = input_dem_dataset.read()[0]
    localdominance_arr = rvt_vis.local_dominance(input_dem_arr)
    localdominance_arr = localdominance_arr.astype('float64')
    profile = input_dem_dataset.profile
    profile.update(dtype='float64')
    output_localdom = rio.open(output_path, "w", **profile)
    output_localdom.write(np.array([localdominance_arr]))

###

# test_slope_aspect(input_dem_path=r"D:\RVT_py\test\TM1_564_146.tif",
#                   output_path=r"D:\RVT_py\test\TM1_564_146_test_slope.tif")

# test_hillshade(input_dem_path=r"D:\RVT_py\test\TM1_564_146.tif",
#                output_path=r"D:\RVT_py\test\TM1_564_146_test_hillsahade.tif")

# test_multi_hillshade(input_dem_path=r"D:\RVT_py\test\TM1_564_146.tif",
#                      output_path=r"D:\RVT_py\test\TM1_564_146_test_multi_hillsahade.tif")

# test_slrm(input_dem_path=r"D:\RVT_py\test\TM1_564_146.tif",
#           output_path=r"D:\RVT_py\test\TM1_564_146_test_SLRM.tif")

# test_sky_view_factor(input_dem_path=r"D:\RVT_py\test\TM1_564_146.tif",
#                      output_svf_path=r"D:\RVT_py\test\TM1_564_146_test_SVF.tif",
#                      output_asvf_path=r"D:\RVT_py\test\TM1_564_146_test_ASVF.tif",
#                      output_opns_path=r"D:\RVT_py\test\TM1_564_146_test_OPNS.tif")

# test_sky_illumination(input_dem_path=r"D:\RVT_py\test\TM1_564_146.tif",
#                       output_path=r"D:\RVT_py\test\TM1_564_146_test_skyillum.tif")

# test_local_dominance(input_dem_path=r"D:\RVT_py\test\TM1_564_146.tif",
#                      output_path=r"D:\RVT_py\test\TM1_564_146_test_localdom.tif")

###

# test_analytical_hillshading(input_dem_path=r"D:\RVT_py\test\manhattan\Manhattan_DSM_1m_clip.tif",
#                             output_path=r"D:\RVT_py\test\manhattan\Manhattan_DSM_1m_clip_test_hillshade.tif")
