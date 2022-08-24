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

def get_da_min_max(image:da.Array)-> Tuple[Union[int, float]]:
    return da.nanmin(image), da.nanmax(image)

def get_da_distribution(image:da.Array, minimum:Union[int, float], maximum:Union[int, float]) -> NDArray[np.float32]:
    image_valid = image[~da.isnan(image)].compute_chunk_sizes()
    distribution = da.percentile(a=image_valid, q=np.array([minimum, 100 - maximum])) 
    #min_lin, max_lin = distribution.compute()
    return distribution

def dask_lin_cutoff_calc_from_perc(image:da.Array, minimum:Union[float, int], maximum:Union[float, int]) -> Dict[str, Union[int, float]]:
    if minimum < 0 or maximum < 0 or minimum > 100 or maximum > 100:
        raise Exception("rvt.blend_funct.lin_cutoff_calc_from_perc: minimum, maximum are percent and have to be in "
                        "range 0-100!")
    if minimum + maximum > 100:
        raise Exception("rvt.blend_funct.lin_cutoff_calc_from_perc: if minimum + maximum > 100% then there are no"
                        " values left! You can't cutoff whole image!")
    min_lin, max_lin = get_da_distribution(image, minimum, maximum).compute() #ok?
    if min_lin == max_lin:
        min_lin, max_lin = get_da_min_max(image)
    return min_lin, max_lin
    # return {"min_lin": min_lin, "max_lin": max_lin}

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


def _blend_images_wrapper(np_chunk_act: NDArray[np.float32],
                          np_chunk_bkg: NDArray[np.float32],
                          blend_mode: str,
                          min_c: Union[int, float, None],
                          max_c: Union[int, float, None]) -> NDArray[np.float32]:
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
                            #drop_axis = 0, new_axis (TODO: Blend mode Normal changes dimension, fix based on condition)
                            meta=np.array((), dtype=np.float32))
    return out_top


def _render_images_wrapper(np_chunk_act: NDArray[np.float32],
                           np_chunk_bkg: NDArray[np.float32],
                           opacity: Union[int,float],
                           bkg_abs_max: Union[int,float],
                           bkg_abs_min: Union[int,float],
                           act_abs_max: Union[int,float],
                           act_abs_min: Union[int,float]) -> NDArray[np.float32]:
    rendered_image = rvt.blend_func.render_images(active = np_chunk_act, background = np_chunk_bkg , opacity = opacity, bkg_abs_max= bkg_abs_max, bkg_abs_min=bkg_abs_min, act_abs_max= act_abs_max, act_abs_min=act_abs_min)
    return rendered_image

def dask_render_images(active:da.Array,
                       background:da.Array, 
                       opacity) -> da.Array:
    active = active.astype(np.float32)
    background = background.astype(np.float32)
    bkg_abs_min, bkg_abs_max = get_da_min_max(background)  ##TODO: If multiple layers, this has to be calculated per layer
    act_abs_min, act_abs_max = get_da_min_max(active)
    _func = partial(_render_images_wrapper,
                    opacity = opacity)
    out_rendered_image = da.map_blocks(_func, 
                                       active, 
                                       background,
                                       bkg_abs_max = bkg_abs_max,
                                       bkg_abs_min = bkg_abs_min,
                                       act_abs_max = act_abs_max,
                                       act_abs_min = act_abs_min,
                                       meta=np.array((), dtype=np.float32))
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
    norm_perc_image = rvt.blend_func.normalize_lin(image = np_chunk, minimum = minimum, maximum = maximum)
    return norm_perc_image

def dask_normalize_perc(image: da.Array, 
                       minimum, 
                       maximum) -> da.Array:
    image.astype(np.float32)
    min_lin, max_lin = dask_lin_cutoff_calc_from_perc(image = image, minimum = minimum, maximum = maximum)
    _func = partial(_normalize_perc_wrapper,
                    minimum = min_lin, 
                    maximum = max_lin)
    out_norm_perc_image = da.map_blocks(_func, 
                                       image, 
                                       dtype = np.float32)
    return out_norm_perc_image


def dask_advanced_normalization(image:da.Array,
                        min_norm,
                        max_norm,
                        normalization) -> Union[None, da.Array]:
    image = image.astype(np.float32)
    if normalization is None:
        out_normalize = image
    else:
        if normalization.lower() == 'value':
            out_normalize = dask_normalize_lin(image = image,
                    minimum= min_norm, maximum = max_norm) 
        elif normalization.lower() == 'percent':
            #min_lin, max_lin = dask_lin_cutoff_calc_from_perc(image = image, minimum = min_norm, maximum = max_norm)
            out_normalize = dask_normalize_perc(image = image, minimum = min_norm, maximum = max_norm) 
    return out_normalize


def dask_normalize_image(image:da.Array, 
                        visualization, 
                        min_norm,
                        max_norm,
                        normalization):
    if visualization is None:  #workaround: If return None -> NoneType has no .compute()
        def _func_none():
            return None
        out_none= da.map_blocks(_func_none, 
                                meta=np.array((), dtype=np.float32))
        return out_none
    else:
        norm_image =  dask_advanced_normalization(image=image,
                            min_norm=min_norm,
                            max_norm=max_norm,
                            normalization=normalization)
        abs_min, abs_max = get_da_min_max(norm_image)
        if abs_max > 1: 
            if visualization.lower() == "multiple directions hillshade" or visualization == "mhs":
                out_norm_image = dask_scale_0_to_1(image=norm_image, abs_max=abs_max, abs_min=abs_min)
            else:
                out_norm_image = dask_scale_0_to_1(image=norm_image, abs_max=abs_max, abs_min=abs_min)
                warnings.warn("rvt.blend.normalize_images_on_layers: unexpected values! max > 1")
            if abs_min < 0:
                warnings.warn("rvt.blend.normalize_images_on_layers: unexpected values! min < 0")

        if visualization.lower() == "slope gradient" or visualization.lower() == "openness - negative" or \
                visualization == "slp" or visualization == "neg_opns":
            def _func_invert(np_chunk):
                return 1 - np_chunk
            out_norm_image = da.map_blocks(_func_invert,
                                            norm_image,
                                            meta=np.array((), dtype=np.float32))
        else:
            out_norm_image=norm_image
    return out_norm_image


def _scale_0_to_1_wrapper(np_chunk: NDArray[np.float32], 
                          abs_max: Union[int, float],
                          abs_min: Union[int, float]) -> NDArray[np.float32]:
    scaled_image = rvt.blend_func.scale_0_to_1(numeric_value= np_chunk, abs_max=abs_max, abs_min=abs_min)
    return scaled_image

def dask_scale_0_to_1(image:da.Array,
                      abs_max,
                      abs_min) -> da.Array:
    image = image.astype(np.float32)
    _func = partial(_scale_0_to_1_wrapper, 
                    abs_max = abs_max, 
                    abs_min=abs_min)
    out_scaled = da.map_blocks(_func, 
                               image,
                               meta=np.array((), dtype=np.float32))
    return out_scaled
