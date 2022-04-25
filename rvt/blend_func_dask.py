"""
Relief Visualization Toolbox – Visualization Functions - Dask

Contains core functions for blending. Wrapped for Dask.

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

# TODO: more testing, find and fix bugs if they exists, wrap additional functions from rvt.blend_func

# python libraries
import rvt.blend_func
import matplotlib as mpl
import matplotlib.cm
import matplotlib.colors
import numpy as np
import warnings
import dask.array as da
from functools import partial
from typing import Union, Dict, List, Any, Tuple, Optional
from nptyping import NDArray


def _gray_scale_to_color_ramp_wrapper(np_chunk: NDArray[np.float32],
                                      colormap: str,
                                      min_colormap_cut: float,
                                      max_colormap_cut: float, 
                                      alpha: bool,
                                      output_8bit: bool) -> NDArray[np.float32]:
    norm_image = rvt.blend_func.gray_scale_to_color_ramp(gray_scale = np_chunk, colormap = colormap, min_colormap_cut = min_colormap_cut, 
                                max_colormap_cut = max_colormap_cut, alpha = alpha, output_8bit = output_8bit)
    return norm_image
    
def dask_gray_scale_to_color_ramp(gray_scale:da.Array,
                                colormap,
                                min_colormap_cut,
                                max_colormap_cut, 
                                alpha,
                                output_8bit) -> da.Array:
    """Turns normalized gray scale np.array to rgba.
    :returns da.array rgba_out: (3D: red 0-255, green 0-255, blue 0-255). 
                If alpha False: (4D: red 0-255, green 0-255, blue 0-255, alpha 0-255) """
    gray_scale = gray_scale.astype(np.float32)
    # data_volume = gray_scale
    _func = partial(_gray_scale_to_color_ramp_wrapper,
                                colormap = colormap,
                                min_colormap_cut = min_colormap_cut,
                                max_colormap_cut = max_colormap_cut, 
                                alpha = alpha,
                                output_8bit = output_8bit) 
    out_normalize = da.map_blocks(_func, 
                                  gray_scale,
                                  new_axis = 0,
                                  dtype = np.float32)
    return out_normalize


def _normalize_image_wrapper(np_chunk: NDArray[np.float32],
                                visualization: str,
                                min_norm : Union[int, float],
                                max_norm : Union[int, float],
                                normalization: str) -> NDArray[np.float32]:
    norm_image = rvt.blend_func.normalize_image(visualization = visualization, image = np_chunk, 
                                min_norm = min_norm, max_norm = max_norm, normalization = normalization)
    return norm_image

def dask_normalize_image(image:da.Array,
                        visualization,
                        min_norm,
                        max_norm,
                        normalization) -> da.Array:
    image = image.astype(np.float32)
    # data_volume = image
    _func = partial(_normalize_image_wrapper,
                    visualization = visualization,
                    min_norm = min_norm, max_norm = max_norm,
                    normalization = normalization) 
    out_normalize = da.map_blocks(_func, 
                                  image,
                                  dtype = np.float32)
    return out_normalize


def _blend_images_wrapper(np_chunk_act: NDArray[np.float32],
                          np_chunk_bkg: NDArray[np.float32],
                          blend_mode: str,
                          min_c: Union[int, float, None],
                          max_c: Union[int, float, None]) -> NDArray[np.float32]: #numpy 2d CHUNK out from vis.function
    top = rvt.blend_func.blend_images(active = np_chunk_act, background = np_chunk_bkg,
                                      blend_mode = blend_mode, min_c = min_c, max_c = max_c)
    return top

def dask_blend_images(active:da.Array,
                      background:da.Array,
                      blend_mode,
                      min_c = None,
                      max_c = None) -> da.Array:
    active = active.astype(np.float32)
    background = background.astype(np.float32)
    _func = partial(_blend_images_wrapper,
                    blend_mode = blend_mode,
                    min_c = min_c,
                    max_c = max_c)
    out_top = da.map_blocks(_func, 
                            active, 
                            background,
                            dtype = np.float32)
    return out_top


def _render_images_wrapper(np_chunk_act: NDArray[np.float32],
                           np_chunk_bkg: NDArray[np.float32],
                           opacity = Union[int,float]) -> NDArray[np.float32]:
    rendered_image = rvt.blend_func.render_images(active = np_chunk_act, background = np_chunk_bkg , opacity = opacity)
    return rendered_image

def dask_render_images(active:da.Array,
                       background:da.Array, 
                       opacity) -> da.Array:
    active = active.astype(np.float32)
    background = background.astype(np.float32)
    _func = partial(_render_images_wrapper,
                    opacity = opacity)
    out_rendered_image = da.map_blocks(_func, 
                                       active, 
                                       background,
                                       dtype = np.float32)
    return out_rendered_image


def _normalize_lin_wrapper(np_chunk: NDArray[np.float32], 
                            minimum: Union[int, float],
                            maximum: Union[int, float]) -> NDArray[np.float32]:
    norm_lin_image = rvt.blend_func.normalize_lin(image = np_chunk, minimum = minimum, maximum = maximum)
    return norm_lin_image
    
def dask_normalize_lin(image: da.Array, 
                       minimum, 
                       maximum) -> da.Array:
    image.astype(np.float32)
    _func = partial(_normalize_lin_wrapper,
                    minimum = minimum, 
                    maximum = maximum)
    out_norm_lin_image = da.map_blocks(_func, 
                                       image, 
                                       dtype = np.float32)
    return out_norm_lin_image


def _normalize_perc_wrapper(np_chunk: NDArray[np.float32], 
                            minimum: Union[int, float],
                            maximum: Union[int, float]) -> NDArray[np.float32]:
    norm_perc_image = rvt.blend_func.normalize_perc(image = np_chunk, minimum = minimum, maximum = maximum)
    return norm_perc_image

def dask_normalize_perc(image: da.Array, 
                       minimum, 
                       maximum) -> da.Array:
    image.astype(np.float32)
    _func = partial(_normalize_perc_wrapper,
                    minimum = minimum, 
                    maximum = maximum)
    out_norm_perc_image = da.map_blocks(_func, 
                                       image, 
                                       dtype = np.float32)
    return out_norm_perc_image
