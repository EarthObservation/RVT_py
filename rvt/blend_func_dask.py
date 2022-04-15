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

def _normalize_image_wrapper(np_chunk: NDArray[np.float32],
                                visualization: str,
                                min_norm : Union[int, float],
                                max_norm : Union[int, float],
                                normalization: str) -> NDArray[np.float32]:
    norm_image = rvt.blend_func.normalize_image(visualization = visualization, image = np_chunk, 
                                min_norm = min_norm, max_norm = max_norm, normalization = normalization)
    return norm_image

def dask_normalize_image(image:da.Array,
                        visualization =  str,
                        min_norm = Union[int, float],
                        max_norm = Union[int, float],
                        normalization = str) -> da.Array:
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


def _blend_images_wrapper(active: NDArray[np.float32],
                          background: NDArray[np.float32],
                          blend_mode: str,
                          min_c: Union[int, float, None],
                          max_c: Union[int, float, None]) -> NDArray[np.float32]: #numpy 2d CHUNK out from vis.function
    top = rvt.blend_func.blend_images(blend_mode = blend_mode, active = active, background = background, 
                                      min_c = min_c, max_c = max_c)
    return top

def dask_blend_images(image:da.Array,
                      image2:da.Array,
                      blend_mode: str,
                      min_c = Union[int, float, None],
                      max_c = Union[int, float, None]) -> da.Array:
    image = image.astype(np.float32)
    # data_volume = image 
    _func = partial(_blend_images_wrapper,
                    blend_mode = blend_mode,
                    min_c = min_c,
                    max_c = max_c)
    out_top = da.map_blocks(_func, 
                            image, 
                            image2,
                            dtype = np.float32)
    return out_top


def _render_images_wrapper(active: NDArray[np.float32],
                           background: NDArray[np.float32],
                           opacity = Union[int,float]) -> NDArray[np.float32]:
    rendered_image = rvt.blend_func.render_images(active = active, background = background , opacity = opacity)
    return rendered_image

def dask_render_images(image:da.Array,
                       image2:da.Array, 
                       opacity: Union[int,float]) -> da.Array:
    image = image.astype(np.float32)
    # data_volume = image 
    _func = partial(_render_images_wrapper,
                    opacity = opacity)
    out_rendered_image = da.map_blocks(_func, 
                                       image, 
                                       image2,
                                       dtype = np.float32)
    return out_rendered_image