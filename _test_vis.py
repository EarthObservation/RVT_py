# rvt.vis quick TEST
# test files: https://www.dropbox.com/sh/p7ia8fk6mywa8y3/AABWuw4wFUvULU7SeNXyWhjka?dl=0

import time
import rasterio as rio
import rvt.visualizations
import numpy as np


def test_slope_aspect(input_dem_path, output_path, ve_factor=1, output_units="degree"):
    input_dem_dataset = rio.open(input_dem_path)
    t = input_dem_dataset.transform
    x_res = t[0]
    y_res = -t[4]
    input_dem_arr = input_dem_dataset.read()[0]

    dict_slp_asp = rvt.vis.slope_aspect(dem=input_dem_arr, resolution_x=x_res,
                                        resolution_y=y_res, ve_factor=ve_factor,
                                        output_units=output_units)
    slope_arr = dict_slp_asp["slope"]
    slope_arr = slope_arr.astype('float32')
    profile = input_dem_dataset.profile
    profile.update(dtype='float32')
    output_slope_dataset = rio.open(output_path, "w", **profile)
    output_slope_dataset.write(np.array([slope_arr]))
    input_dem_dataset.close()
    output_slope_dataset.close()


def test_hillshade(input_dem_path, output_path, sun_azimuth=315, sun_elevation=35):
    input_dem_dataset = rio.open(input_dem_path)
    t = input_dem_dataset.transform
    x_res = t[0]
    y_res = -t[4]
    input_dem_arr = input_dem_dataset.read()[0]

    hillshade_arr = rvt.vis.hillshade(dem=input_dem_arr, resolution_x=x_res,
                                      resolution_y=y_res, sun_azimuth=sun_azimuth, sun_elevation=sun_elevation)
    hillshade_arr = hillshade_arr.astype('float32')
    profile = input_dem_dataset.profile
    profile.update(dtype='float32')
    output_hillshade_dataset = rio.open(output_path, "w", **profile)
    output_hillshade_dataset.write(np.array([hillshade_arr]))
    input_dem_dataset.close()
    output_hillshade_dataset.close()


def test_multi_hillshade(input_dem_path, output_path, nr_directions=16, sun_elevation=35):
    input_dem_dataset = rio.open(input_dem_path)
    t = input_dem_dataset.transform
    x_res = t[0]
    y_res = -t[4]
    input_dem_arr = input_dem_dataset.read()[0]
    multi_hillshade_arr = rvt.vis.multi_hillshade(dem=input_dem_arr, resolution_x=x_res,
                                                  resolution_y=y_res, nr_directions=nr_directions,
                                                  sun_elevation=sun_elevation)
    multi_hillshade_arr = multi_hillshade_arr.astype('float32')
    profile = input_dem_dataset.profile
    profile.update(dtype='float32', count=16)
    output_multi_hillshade_dataset = rio.open(output_path, "w", **profile)
    output_multi_hillshade_dataset.write(multi_hillshade_arr)
    input_dem_dataset.close()
    output_multi_hillshade_dataset.close()


def test_slrm(input_dem_path, output_path):
    input_dem_dataset = rio.open(input_dem_path)
    input_dem_arr = input_dem_dataset.read()[0]

    slrm_arr = rvt.vis.slrm(input_dem_arr)
    slrm_arr = slrm_arr.astype('float32')
    profile = input_dem_dataset.profile
    profile.update(dtype='float32')
    output_slrm_dataset = rio.open(output_path, "w", **profile)
    output_slrm_dataset.write(np.array([slrm_arr]))
    input_dem_dataset.close()
    output_slrm_dataset.close()


def test_sky_view_factor(
        input_dem_path,
        output_svf_path=None,
        output_asvf_path=None,
        output_opns_path=None
):
    input_dem_dataset = rio.open(input_dem_path)
    t = input_dem_dataset.transform
    x_res = t[0]
    y_res = -t[4]
    input_dem_arr = input_dem_dataset.read()[0]
    no_data = input_dem_dataset.nodata

    if output_svf_path:
        compute_svf = True
    else:
        compute_svf = False

    if output_asvf_path:
        compute_asvf = True
    else:
        compute_asvf = False

    if output_opns_path:
        compute_opns = True
    else:
        compute_opns = False

    dict_svf_asvf_opns = rvt.vis.sky_view_factor(
        dem=input_dem_arr,
        resolution=x_res,
        compute_svf=compute_svf,
        compute_asvf=compute_asvf,
        compute_opns=compute_opns,
        no_data=no_data
    )

    svf_arr = dict_svf_asvf_opns["svf"]
    svf_arr = svf_arr.astype('float32')

    profile = input_dem_dataset.profile
    profile.update(dtype='float32')
    output_svf = rio.open(output_svf_path, "w", **profile)
    output_svf.write(np.array([svf_arr]))
    output_svf.close()

    asvf_arr = dict_svf_asvf_opns["asvf"]
    asvf_arr = asvf_arr.astype('float32')

    output_asvf = rio.open(output_asvf_path, "w", **profile)
    output_asvf.write(np.array([asvf_arr]))
    output_asvf.close()

    opns_arr = dict_svf_asvf_opns["opns"]
    opns_arr = opns_arr.astype('float32')

    output_opns = rio.open(output_opns_path, "w", **profile)
    output_opns.write(np.array([opns_arr]))
    output_opns.close()
    input_dem_dataset.close()


def test_sky_illumination(input_dem_path, output_path):
    input_dem_dataset = rio.open(input_dem_path)
    t = input_dem_dataset.transform
    x_res = t[0]
    y_res = -t[4]
    input_dem_arr = input_dem_dataset.read()[0]
    skyillumination_arr = rvt.vis.sky_illumination(dem=input_dem_arr, resolution=x_res, compute_shadow=True)
    skyillumination_arr = skyillumination_arr.astype('float32')
    profile = input_dem_dataset.profile
    profile.update(dtype='float32')
    output_skyilum = rio.open(output_path, "w", **profile)
    output_skyilum.write(np.array([skyillumination_arr]))
    output_skyilum.close()
    input_dem_dataset.close()


def test_local_dominance(input_dem_path, output_path):
    input_dem_dataset = rio.open(input_dem_path)
    input_dem_arr = input_dem_dataset.read()[0]
    localdominance_arr = rvt.vis.local_dominance(input_dem_arr)
    localdominance_arr = localdominance_arr.astype('float32')
    profile = input_dem_dataset.profile
    profile.update(dtype='float32')
    output_localdom = rio.open(output_path, "w", **profile)
    output_localdom.write(np.array([localdominance_arr]))
    output_localdom.close()
    input_dem_dataset.close()


###

# test_slope_aspect(input_dem_path=r"test_data\TM1_564_146.tif",
#                   output_path=r"test_data\TM1_564_146_test_slope.tif")

# test_hillshade(input_dem_path=r"test_data\TM1_564_146.tif",
#                output_path=r"test_data\TM1_564_146_test_hillsahade.tif")

# test_multi_hillshade(input_dem_path=r"test_data\TM1_564_146.tif",
#                      output_path=r"test_data\TM1_564_146_test_multi_hillsahade.tif")

# test_slrm(input_dem_path=r"test_data\TM1_564_146.tif",
#           output_path=r"test_data\TM1_564_146_test_SLRM.tif")

# test_sky_view_factor(input_dem_path=r"test_data\TM1_564_146.tif",
#                      output_svf_path=r"test_data\TM1_564_146_test_SVF.tif",
#                      output_asvf_path=r"test_data\TM1_564_146_test_ASVF.tif",
#                      output_opns_path=r"test_data\TM1_564_146_test_OPNS.tif")

# test_sky_illumination(input_dem_path=r"test_data\TM1_564_146.tif",
#                       output_path=r"test_data\TM1_564_146_test_skyillum.tif")

# test_local_dominance(input_dem_path=r"test_data\TM1_564_146.tif",
#                      output_path=r"test_data\TM1_564_146_test_localdom.tif")
