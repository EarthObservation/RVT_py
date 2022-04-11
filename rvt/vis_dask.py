"""
Relief Visualization Toolbox – Wrapper functions for Visualization Functions

Module for handling dask version of computing the visualizations.

Credits:
    Žiga Kokalj (ziga.kokalj@zrc-sazu.si)
    Krištof Oštir (kristof.ostir@fgg.uni-lj.si)
    Klemen Zakšek
    Peter Pehani
    Klemen Čotar
    Maja Somrak
    Žiga Maroh

Copyright:
    2010-2020 Research Centre of the Slovenian Academy of Sciences and Arts
    2016-2020 University of Ljubljana, Faculty of Civil and Geodetic Engineering
"""

from logging import raiseExceptions
import warnings
import rvt.vis
import numpy as np
from functools import partial
import dask
import dask.array as da
from typing import Union, Dict, List, Any, Tuple, Optional
from nptyping import NDArray


def _dask_byte_scale_wrapper(data: NDArray[np.float32],
                             c_min: Union[int, float],
                             c_max: Union[int, float],
                             high: int,
                             low: int,
                             no_data: Union[int, float, None]) -> NDArray[np.float32]:
    scaled_image = rvt.vis.byte_scale(data = data, c_min = c_min, c_max = c_max, high = high, 
                              low = low, no_data = no_data )
    return scaled_image

def dask_byte_scale(data: da.Array, 
                    c_min = None,
                    c_max = None,
                    high = 255,
                    low = 0,
                    no_data = None) -> da.Array:
    """Maps byte_scale function over dask.array (no overlapping computations).

    :param data: The input dask array (temp 2D).
    :param c_min: Integer or float - bias scaling of small values.
    :param c_max: Integer or float - bias scaling of large values.
    :param high: The integer to scale the max value to.
    :param low: The integer to scale the min value to.
    :param no_data: The value that represents no data.
    :return: A dask array of byte scaled image `data`.

    TODO: Multiple bands."""
    _func = partial(_dask_byte_scale_wrapper, 
                    c_min = c_min, c_max = c_max,
                    high = high, low = low,
                    no_data = no_data)
    out_scaled_image = da.map_blocks(_func, 
                                     data, 
                                     dtype = np.float32)
    return out_scaled_image


def _slope_aspect_wrapper(np_chunk: NDArray[np.float32], resolution_x: Union[int, float],
                            resolution_y: Union[int, float], ve_factor: Union[int, float], output_units: str,
                            no_data : Union[int, None]) -> NDArray[np.float32]: 
    """Wrapper function for vis.dask_slope_aspect. Calculates `slope` and `aspect` for each dask array chunk (np.array). 
    Returns np.array dim (2, y, x).""" 
    result_dict = rvt.vis.slope_aspect(dem = np_chunk[0,:,:].squeeze(), resolution_x = resolution_x, 
                                        resolution_y= resolution_y, ve_factor = ve_factor, 
                                        output_units= output_units, no_data= no_data)
    slope_out = result_dict["slope"]
    aspect_out = result_dict["aspect"]
    # return slope_out
    output_for_dask_slo_asp = np.stack((slope_out, aspect_out))
    return output_for_dask_slo_asp

def dask_slope_aspect(input_dem: da.Array, 
                    resolution_x,
                    resolution_y,
                    ve_factor,
                    output_units,
                    no_data = None) -> da.Array:
    """Maps slope_aspect function over dask.array (with overlap of `depth`).

    :param input_dem: The input dask array.
    :param resolution_x:  Integer or float of DEM resolution in X direction.
    :param resolution_y: Integer or float of DEM resolution in Y direction.
    :param ve_factor: Integer or float of vertical exaggeration factor.
    :param output_units: Output units: percent, degree, radian.
    :param no_data: The value that represents no data.
    :return: A 3D dask array of calculated `slope` and `aspect` of the `input_dem` raster.
    FIXME: Most outter edges (set in vis.slope_aspect) not ok for both outputs."""

    input_dem = input_dem.astype(np.float32)
    data_volume = da.stack([input_dem, input_dem])[[0,0]]  ##magic numbers to avoid rechunking, fix when map_overlap supports input and output of different shapes
    _func = partial(_slope_aspect_wrapper,
                    resolution_x = resolution_x,
                    resolution_y = resolution_y,
                    ve_factor = ve_factor,
                    output_units = output_units,
                    no_data = no_data)
    depth = { 0:0, 1: 1, 2: 1}
    boundary = {0: 0, 1: 'periodic', 2: 'periodic'}         # (outter edges of) slope_out not ok
    # boundary = {0: np.nan, 1: np.nan}                     # (outter edges of) aspect_out not ok
    out_slp = data_volume.map_overlap(_func,
                            depth=depth,
                            boundary=boundary,
                            # new_axis = 0,                 #works only with map_blocks
                            meta=np.array((), dtype=np.float32))
    return out_slp


def _hillshade_wrapper(np_chunk: NDArray[np.float32],
                    resolution_x: Union[int,float],
                    resolution_y: Union[int,float],
                    sun_azimuth: Union[int, float],
                    sun_elevation: Union[int, float],
                    ve_factor : Union[int, float],
                    no_data: Union[int, None]) -> NDArray[np.float32]: 
    """Wrapper function for vis.dask_hillshade. Calculates `hillshade` for each dask array chunk (np.array). 
    Returns np.array dim (x, y).""" 
    result_out = rvt.vis.hillshade(dem = np_chunk, resolution_x = resolution_x, resolution_y = resolution_y,
                                  sun_azimuth = sun_azimuth, sun_elevation = sun_elevation,
                                  ve_factor=ve_factor, no_data = no_data)
    output_for_dask_hs = result_out
    return output_for_dask_hs 
    
def dask_hillshade(input_dem: da.Array,
                    resolution_x,
                    resolution_y ,
                    sun_azimuth,
                    sun_elevation,
                    ve_factor,
                    no_data = None) -> da.Array:
    """Maps hillshade function over dask.array (with overlap of `depth`).

    :param input_dem: The input dask array.
    :param resolution_x:  Integer or float of DEM resolution in X direction.
    :param resolution_y: Integer or float of DEM resolution in Y direction.
    :param sun_azimuth: Solar azimuth angle in degrees.
    :param sun_elevation: Solar vertical angle in degrees.
    :param ve_factor: Integer or float of vertical exaggeration factor.
    :param no_data: The value that represents no data.
    :return: A 2D dask array of calculated `hillshade` of the `input_dem` raster.
    TODO: Can  params `slope` or `aspect` be calculated input?!"""

    input_dem = input_dem.astype(np.float32)
    data_volume = input_dem 
    _func = partial(_hillshade_wrapper, 
                    resolution_x = resolution_x, 
                    resolution_y = resolution_y,
                    sun_azimuth = sun_azimuth,
                    sun_elevation = sun_elevation,
                    ve_factor = ve_factor,
                    no_data = no_data)
    depth =  {0: 1, 1: 1}
    boundary = {0: 'periodic', 1: 'periodic'}    
    out_hs = data_volume.map_overlap(_func,
                                    depth=depth,
                                    boundary=boundary,
                                    meta=np.array((), dtype=np.float32))
    return out_hs


def _sky_view_factor_wrapper(np_chunk: NDArray[np.float32], resolution: Union[int, float],  
                            compute_svf:bool, compute_opns: bool, compute_asvf: bool,
                            svf_n_dir : int, svf_r_max:Union[int,float],
                            svf_noise : Union[int, float], asvf_dir:int, 
                            asvf_level: Union[int, float], ve_factor:Union[int,float],
                            no_data: Union[int, None]) -> NDArray[np.float32]: 
    """Wrapper function for vis.sky_viw_factor. Calculates `svf`, `opns` or `asvf` for each dask array chunk (np.array). 
    Returns np.array dim (x, y). Assumes calculation of only one visualization at a time.""" 
    result_dict = rvt.vis.sky_view_factor(dem = np_chunk, resolution = resolution, compute_svf = compute_svf,
                                                compute_opns = compute_opns, compute_asvf= compute_asvf,
                                                svf_n_dir=svf_n_dir, svf_r_max=svf_r_max,
                                                svf_noise=svf_noise, asvf_dir=asvf_dir,
                                                asvf_level=asvf_level, ve_factor=ve_factor,
                                                no_data = no_data)
    svf_out, = result_dict.values()            
    return svf_out   

def dask_sky_view_factor(input_dem: da.Array,
                        resolution, compute_svf,
                        compute_opns, compute_asvf,
                        svf_n_dir, svf_r_max,
                        svf_noise, asvf_dir,
                        asvf_level, ve_factor,
                        no_data)-> da.Array: 
    """Maps sky_view_factor function over dask.array (with overlap of `depth`).

    :param input_dem: The input dask array.
    :param resolution:  Integer or float of DEM resolution (pixel resolution?).
    :param compute_svf: True if computing `svf` visualization, else False.
    :param compute_opns: True if computing `opns` visualization, else False.
    :param compute_asvf: True if computing `asvf` visualization, else False.
    :param svf_n_dir: Number of directions of calculating `svf` visualization.
    :param svf_r_max: Maximal search radius in pixels.
    :param svf_noise: The level of noise removal.
    :param asvf_dir: Dirction of anisotropy.
    :param asvf_level: Level of anisotropy.
    :param ve_factor: Integer or float of vertical exaggeration factor.
    :param no_data: The value that represents no data.
    :return: A 2D dask array of calculated `svf`, `opns` or `asvf` of the `input_dem` raster."""

    input_dem = input_dem.astype(np.float32)
    data_volume = input_dem 
    _func = partial(_sky_view_factor_wrapper,
                    resolution = resolution, compute_svf = compute_svf,
                    compute_opns = compute_opns, compute_asvf = compute_asvf,
                    svf_n_dir=svf_n_dir, svf_r_max=svf_r_max,
                    svf_noise=svf_noise, asvf_dir=asvf_dir,
                    asvf_level=asvf_level, ve_factor=ve_factor,
                    no_data = no_data) 
    radius_max = svf_r_max
    depth = { 0: radius_max, 1: radius_max} 
    boundary = { 0: "reflect", 1: "reflect"}
    out_svf = data_volume.map_overlap(_func,
                                        depth = depth,
                                        boundary = boundary,
                                        meta=np.array((), dtype = np.float32))
    return out_svf