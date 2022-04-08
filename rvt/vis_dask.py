from logging import raiseExceptions
import warnings
import rvt.vis
import numpy as np
from functools import partial
import dask
import dask.array as da
from typing import Union, Dict, List, Any, Tuple, Optional
from nptyping import NDArray

def _slope_aspect_wrapper(np_chunk: NDArray[np.float32], resolution_x: Union[int, float],
                            resolution_y: Union[int, float], ve_factor: Union[int, float], output_units: str,
                            no_data : Union[int, None]) -> NDArray[np.float32]: 
    result_dict = rvt.vis.slope_aspect(dem = np_chunk[0,:,:].squeeze(), resolution_x = resolution_x, 
                                        resolution_y= resolution_y, ve_factor = ve_factor, 
                                        output_units= output_units, no_data= no_data)
    slope_out = result_dict["slope"]
    aspect_out = result_dict["aspect"]
    # return slope_out
    output_for_dask_slo_asp = np.stack((slope_out, aspect_out))
    return output_for_dask_slo_asp

def dask_slope_aspect( input_dem: da.Array, 
                    resolution_x,
                    resolution_y,
                    ve_factor,
                    output_units,
                    no_data = None) -> da.Array:
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
                            svf_n_dir : Union[int, float], svf_r_max:Union[int,float],
                            svf_noise : Union[int, float], asvf_dir:Union[int,float], 
                            asvf_level: Union[int, float], ve_factor:Union[int,float],
                            no_data: Union[int, None]) -> NDArray[np.float32]:  
    result_dict = rvt.vis.sky_view_factor(dem = np_chunk, resolution = resolution, compute_svf = compute_svf,
                                                compute_opns = compute_opns, compute_asvf= compute_asvf,
                                                svf_n_dir=svf_n_dir, svf_r_max=svf_r_max,
                                                svf_noise=svf_noise, asvf_dir=asvf_dir,
                                                asvf_level=asvf_level, ve_factor=ve_factor,
                                                no_data = no_data)
    svf_out, = result_dict.values()             #assume calculation of ONLY ONE vis at a time and return
    return svf_out   

def dask_sky_view_factor(input_dem: da.Array,
                        resolution, compute_svf,
                        compute_opns, compute_asvf,
                        svf_n_dir, svf_r_max,
                        svf_noise, asvf_dir,
                        asvf_level, ve_factor,
                        no_data)-> da.Array: 
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