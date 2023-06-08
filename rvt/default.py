"""
Relief Visualization Toolbox – Visualization Functions

Contains all default values for visualisation functions, which can be changed.
Allows computing from rvt.visualization with using defined default values and saving
output rasters with default names (dependent on default values).

Credits:
    Žiga Kokalj (ziga.kokalj@zrc-sazu.si)
    Krištof Oštir (kristof.ostir@fgg.uni-lj.si)
    Klemen Zakšek
    Peter Pehani
    Klemen Čotar
    Maja Somrak
    Žiga Maroh
    Nejc Čož

Copyright:
    2010-2022 Research Centre of the Slovenian Academy of Sciences and Arts
    2016-2022 University of Ljubljana, Faculty of Civil and Geodetic Engineering
"""

import warnings
from enum import Enum
from pathlib import Path
from typing import Optional, Tuple

import rvt.vis
import rvt.blend_func
import rvt.tile
import os
from osgeo import gdal
import numpy as np
import json
import datetime
import time


class RVTVisualization(Enum):
    SLOPE = "slp"
    HILLSHADE = "hs"
    SHADOW = "shd"
    MULTI_HILLSHADE = "mhs"
    SIMPLE_LOCAL_RELIEF_MODEL = "slrm"
    SKY_VIEW_FACTOR = "svf"
    ANISOTROPIC_SKY_VIEW_FACTOR = "asvf"
    POSITIVE_OPENNESS = "pos_opns"
    NEGATIVE_OPENNESS = "neg_opns"
    SKY_ILLUMINATION = "sim"
    LOCAL_DOMINANCE = "ld"
    MULTI_SCALE_RELIEF_MODEL = "msrm"
    MULTI_SCALE_TOPOGRAPHIC_POSITION = "mstp"


class DefaultValues:
    """
    Class which define layer for blending. BlenderLayer is basic element in BlenderCombination.layers list.

    Attributes
    ----------
    overwrite : bool
        When saving visualisation functions and file already exists, if 0 it doesn't compute it, if 1 it overwrites it.
    ve_factor : float
        For all visualization functions. Vertical exaggeration.
    slp_compute : bool
        If compute Slope. Parameter for GUIs.
    slp_output_units : str
        Slope. Output units [radian, degree, percent].
    hs_compute : bool
        If compute Hillshade. Parameter for GUIs.
    hs_sun_azi : int
        Hillshade. Solar azimuth angle (clockwise from North) in degrees.
    hs_sun_el : int
        Hillshade. Solar vertical angle (above the horizon) in degrees.
    hs_shadow : bool
        Hillshade. If 1 (Ture) computes binary shadow raster, if 0 (False) it doesn't.
    mhs_compute : bool
        If compute Multi directional hillshade. Parameter for GUIs.
    mhs_nr_dir : int
        Multi directional hillshade. Number of solar azimuth angles (clockwise from North).
    mhs_sun_el : int
        Multi directional hillshade. Solar vertical angle (above the horizon) in degrees.
    slrm_compute : bool
        If compute Simple local relief model. Parameter for GUIs.
    slrm_rad_cell : int
        Simple local relief model. Radius for trend assessment in pixels.
    svf_compute : bool
        If compute Sky-View Factor. Parameter for GUIs.
    svf_n_dir : int
        Sky-View Factor (Anisotropic Sky-View Factor, Openness). Number of directions.
    svf_r_max : int
        Sky-View Factor (Anisotropic Sky-View Factor, Openness). Maximal search radius in pixels.
    svf_noise : int
        Sky-View Factor (Anisotropic Sky-View Factor, Openness). The level of noise remove [0-don't remove, 1-low, 2-med
        , 3-high].
    asvf_compute : bool
        If compute Anisotropic Sky-View Factor. Parameter for GUIs.
    asvf_dir : int
        Anisotropic Sky-View Factor. Direction of anisotropy in degrees.
    asvf_level : int
        Anisotropic Sky-View Factor. Level of anisotropy [1-low, 2-high].
    pos_opns_compute : bool
        If compute Positive Openness. Parameter for GUIs.
    neg_opns_compute : bool
        If compute Negative Openness. Parameter for GUIs.
    sim_compute : bool
        If compute Sky illumination. Parameter for GUIs.
    sim_sky_mod : str
        Sky illumination. Sky model [overcast, uniform].
    sim_compute_shadow : bool
        Sky illumination. If 1 it computes shadows, if 0 it doesn't.
    sim_nr_dir : int
        Sky illumination. Number of directions.
    sim_shadow_dist : int
        Sky illumination. Max shadow modeling distance in pixels.
    sim_shadow_az : int
        Sky illumination. Shadow azimuth in degrees.
    sim_shadow_el : int
        Sky illumination. Shadow elevation in degrees.
    ld_compute : bool
        If compute Local dominance. Parameter for GUIs.
    ld_min_rad : int
        Local dominance. Minimum radial distance (in pixels) at which the algorithm starts with visualization
        computation.
    ld_max_rad : int
        Local dominance. Maximum radial distance (in pixels) at which the algorithm ends with visualization computation.
    ld_rad_inc : int
        Local dominance. Radial distance steps in pixels.
    ld_anglr_res : int
        Local dominance. Angular step for determination of number of angular directions.
    ld_observer_h : float
        Local dominance. Height at which we observe the terrain.
    msrm_compute : bool
        If compute Multi-scale relief model. Parameter for GUIs.
    msrm_feature_min : float
        Minimum size of the feature you want to detect in meters.
    msrm_feature_max : float
        Maximum size of the feature you want to detect in meters.
    msrm_scaling_factor : int
        Scaling factor.
    mstp_compute : bool
        If compute Multi-scale topographic position (MSTP).
    mstp_local_scale : tuple(min_radius, max_radius, step)
        Local scale minimum radius, maximum radius and step in pixels to calculate maximum mean deviation from
        elevation.
        All have to be integers!
    mstp_meso_scale : tuple(min_radius, max_radius, step)
        Meso scale minimum radius, maximum radius and step in pixels to calculate maximum mean deviation from elevation.
        All have to be integers!
    mstp_broad_scale : tuple(min_radius, max_radius, step)
        Broad scale minimum radius, maximum radius and step in pixels to calculate maximum mean deviation from
        elevation.
        All have to be integers!
    mstp_lightness : float
        Lightness of image.
    slp_save_float : bool
        Slope. If 1 (True) it saves float, if 0 (False) it doesn't.
    hs_save_float : bool
        Hillshade. If 1 (True) it saves float, if 0 (False) it doesn't.
    mhs_save_float : bool
        Multi hillshade. If 1 (True) it saves float, if 0 (False) it doesn't.
    slrm_save_float : bool
        Simplified local relief model. If 1 (True) it saves float, if 0 (False) it doesn't.
    svf_save_float : bool
        Sky-view factor (asvf, pos_opns). If 1 (True) it saves float, if 0 (False) it doesn't.
    neg_opns_save_float : bool
        Negative openness. If 1 (True) it saves float, if 0 (False) it doesn't.
    sim_save_float : bool
        Sky illumination. If 1 (True) it saves float, if 0 (False) it doesn't.
    ld_save_float : bool
        Local dominance. If 1 (True) it saves float, if 0 (False) it doesn't.
    msrm_save_float : bool
        Multi-scale relief model. If 1 (True) it saves float, if 0 (False) it doesn't.
    slp_save_8bit : bool
        Slope. If 1 (True) it saves 8bit, if 0 (False) it doesn't.
    hs_save_8bit : bool
        Hillshade. If 1 (True) it saves 8bit, if 0 (False) it doesn't.
    mhs_save_8bit : bool
        Multi hillshade. If 1 (True) it saves 8bit, if 0 (False) it doesn't.
    slrm_save_8bit : bool
        Simplified local relief model. If 1 (True) it saves 8bit, if 0 (False) it doesn't.
    svf_save_8bit : bool
        Sky-view factor (asvf, pos_opns). If 1 (True) it saves 8bit, if 0 (False) it doesn't.
    neg_opns_save_8bit : bool
        Negative openness. If 1 (True) it saves 8bit, if 0 (False) it doesn't.
    sim_save_8bit : bool
        Sky illumination. If 1 (True) it saves 8bit, if 0 (False) it doesn't.
    ld_save_8bit : bool
        local dominance. If 1 (True) it saves 8bit, if 0 (False) it doesn't.
    msrm_save_8bit : bool
        Multi-scale relief model. If 1 (True) it saves 8bit, if 0 (False) it doesn't.
    slp_bytscl : tuple(mode, min, max)
        Slope, bytescale (0-255) for 8bit raster. Mode can be 'value' (linear stretch) or 'percent'
        (histogram equalization, cut-off). Values min and max define stretch/cut-off borders.
    hs_bytscl : tuple(mode, min, max)
        Hillshade, linear stretch, bytescale (0-255) for 8bit raster. Mode can be 'value' or 'percent'
        (cut-off units). Values min and max define stretch borders (in mode units).
    mhs_bytscl : tuple(mode, min, max)
        Multi directional hillshade, linear stretch, bytescale (0-255) for 8bit raster. Mode can be 'value' or 'percent'
        (cut-off units). Values min and max define stretch borders (in mode units).
    slrm_bytscl : tuple(mode, min, max)
        Simplified local relief model, linear stretch, bytescale (0-255) for 8bit raster. Mode can be 'value' or
        'percent' (cut-off units). Values min and max define stretch borders (in mode units).
    svf_bytscl : tuple(mode, min, max)
        Sky-view factor, linear stretch, bytescale (0-255) for 8bit raster. Mode can be 'value' or 'percent'
        (cut-off units). Values min and max define stretch borders (in mode units).
    asvf_bytscl : tuple(mode, min, max)
        Anisotropic Sky-view factor, linear stretch, bytescale (0-255) for 8bit raster. Mode can be 'value' or 'percent'
        (cut-off units). Values min and max define stretch borders (in mode units).
    pos_opns_bytscl : tuple(mode, min, max)
        Positive Openness, linear stretch, bytescale (0-255) for 8bit raster. Mode can be 'value' or 'percent'
        (cut-off units). Values min and max define stretch borders (in mode units).
    neg_opns_bytscl : tuple(mode, min, max)
        Negative Openness, linear stretch, bytescale (0-255) for 8bit raster. Mode can be 'value' or 'percent'
        (cut-off units). Values min and max define stretch borders (in mode units).
    sim_bytscl : tuple(mode, min, max)
        Sky illumination, linear stretch, bytescale (0-255) for 8bit raster. Mode can be 'value' or 'percent'
        (cut-off units). Values min and max define stretch borders (in mode units).
    ld_bytscl : tuple(mode, min, max)
        Local dominance, linear stretch, bytescale (0-255) for 8bit raster. Mode can be 'value' or 'percent'
        (cut-off units). Values min and max define stretch borders (in mode units).
    msrm_bytscl : tuple(mode, min, max)
        Multi-scale relief model, linear stretch, bytescale (0-255) for 8bit raster. Mode can be 'value' or 'percent'
        (cut-off units). Values min and max define stretch borders (in mode units).
    tile_size_limit : int
        If array size bigger than tile_size_limit it uses saving tile by tile (rvt.tile module).
    tile_size : tuple(x_size, y_size)
        Size of single tile when saving tile by tile.
    """

    def __init__(self):
        self.overwrite = 0  # (0=False, 1=True)
        self.ve_factor = 1
        # slope gradient
        self.slp_compute = 0
        self.slp_output_units = "degree"
        # hillshade
        self.hs_compute = 1
        self.hs_sun_azi = 315
        self.hs_sun_el = 35
        self.hs_shadow = 0
        # multi hillshade
        self.mhs_compute = 0
        self.mhs_nr_dir = 16
        self.mhs_sun_el = 35
        # simple local relief model
        self.slrm_compute = 0
        self.slrm_rad_cell = 20
        # sky view factor
        self.svf_compute = 0
        self.svf_n_dir = 16
        self.svf_r_max = 10
        self.svf_noise = 0
        # anisotropic sky-view factor
        self.asvf_compute = 0
        self.asvf_dir = 315
        self.asvf_level = 1
        # positive openness
        self.pos_opns_compute = 0
        # negative openness
        self.neg_opns_compute = 0
        # sky_illum
        self.sim_compute = 0
        self.sim_sky_mod = "overcast"
        self.sim_compute_shadow = 0
        self.sim_shadow_dist = 100
        self.sim_nr_dir = 32
        self.sim_shadow_az = 315
        self.sim_shadow_el = 35
        # local dominance
        self.ld_compute = 0
        self.ld_min_rad = 10
        self.ld_max_rad = 20
        self.ld_rad_inc = 1
        self.ld_anglr_res = 15
        self.ld_observer_h = 1.7
        # multi-scale relief model
        self.msrm_compute = 0
        self.msrm_feature_min = 0
        self.msrm_feature_max = 20
        self.msrm_scaling_factor = 2
        # multi-scale topographic position
        self.mstp_compute = 0
        self.mstp_local_scale = (3, 21, 2)
        self.mstp_meso_scale = (23, 203, 18)
        self.mstp_broad_scale = (223, 2023, 180)
        self.mstp_lightness = 1.2
        # save float
        self.slp_save_float = 1
        self.hs_save_float = 1
        self.mhs_save_float = 1
        self.slrm_save_float = 1
        self.svf_save_float = 1
        self.neg_opns_save_float = 1
        self.sim_save_float = 1
        self.ld_save_float = 1
        self.msrm_save_float = 1
        self.mstp_save_float = 1
        # save 8bit
        self.slp_save_8bit = 0
        self.hs_save_8bit = 0
        self.mhs_save_8bit = 0
        self.slrm_save_8bit = 0
        self.svf_save_8bit = 0
        self.neg_opns_save_8bit = 0
        self.sim_save_8bit = 0
        self.ld_save_8bit = 0
        self.msrm_save_8bit = 0
        self.mstp_save_8bit = 0
        # 8-bit bytescale parameters
        self.slp_bytscl = ("value", 0.00, 51.00)
        self.hs_bytscl = ("value", 0.00, 1.00)
        self.mhs_bytscl = ("value", 0.00, 1.00)
        self.slrm_bytscl = ("value", -2.00, 2.00)
        self.svf_bytscl = ("value", 0.6375, 1.00)
        self.asvf_bytscl = ("value", 0.70, 0.90)
        self.pos_opns_bytscl = ("value", 60.00, 95.00)
        self.neg_opns_bytscl = ("value", 60.00, 95.00)
        self.sim_bytscl = ("percent", 0.25, 0.00)
        self.ld_bytscl = ("value", 0.50, 1.80)
        self.msrm_bytscl = ("value", -2.50, 2.50)
        self.mstp_bytscl = ("value", 0.00, 1.00)
        # tile
        self.tile_size_limit = (
            10000 * 10000
        )  # if arr size > tile_size limit, it uses tile module
        self.tile_size = (
            4000,
            4000,
        )  # size of single tile when using tile module (x_size, y_size)

    def save_default_to_file(self, file_path=None):
        """Saves default attributes into .json file."""
        data = {
            "default_settings": {
                "overwrite": {
                    "value": self.overwrite,
                    "description": "When saving visualisation functions and file already exists, if 0 "
                    "it doesn't compute it, if 1 it overwrites it.",
                },
                "ve_factor": {
                    "value": self.ve_factor,
                    "description": "Vertical exaggeration.",
                },
                "Hillshade": {
                    "hs_compute": {
                        "value": self.hs_compute,
                        "description": "If compute Hillshade. Parameter for GUIs.",
                    },
                    "hs_sun_azi": {
                        "value": self.hs_sun_azi,
                        "description": "Solar azimuth angle (clockwise from North) in "
                        "degrees.",
                    },
                    "hs_sun_el": {
                        "value": self.hs_sun_el,
                        "description": "Solar vertical angle (above the horizon) in "
                        "degrees.",
                    },
                    "hs_save_float": {
                        "value": self.hs_save_float,
                        "description": "If 1 it saves float raster, if 0 it doesn't.",
                    },
                    "hs_save_8bit": {
                        "value": self.hs_save_8bit,
                        "description": "If 1 it saves 8bit raster, if 0 it doesn't.",
                    },
                    "hs_shadow": {
                        "value": self.hs_shadow,
                        "description": "If 1 it saves shadow binary raster, if 0 it doesn't.",
                    },
                    "hs_bytscl": {
                        "mode": self.hs_bytscl[0],
                        "min": self.hs_bytscl[1],
                        "max": self.hs_bytscl[2],
                        "description": "Linear stretch and byte scale (0-255) for 8bit raster. "
                        "Mode can be 'value' or 'percent' (cut-off units). "
                        "Values min and max define stretch borders (in mode units).",
                    },
                },
                "Multiple directions hillshade": {
                    "mhs_compute": {
                        "value": self.mhs_compute,
                        "description": "If compute Multiple directions hillshade."
                        " Parameter for GUIs.",
                    },
                    "mhs_nr_dir": {
                        "value": self.mhs_nr_dir,
                        "description": "Number of solar azimuth angles (clockwise "
                        "from North).",
                    },
                    "mhs_sun_el": {
                        "value": self.mhs_sun_el,
                        "description": "Solar vertical angle (above the horizon) in "
                        "degrees.",
                    },
                    "mhs_save_float": {
                        "value": self.mhs_save_float,
                        "description": "If 1 it saves float raster, if 0 it doesn't.",
                    },
                    "mhs_save_8bit": {
                        "value": self.mhs_save_8bit,
                        "description": "If 1 it saves 8bit raster, if 0 it doesn't.",
                    },
                    "mhs_bytscl": {
                        "mode": self.mhs_bytscl[0],
                        "min": self.mhs_bytscl[1],
                        "max": self.mhs_bytscl[2],
                        "description": "Linear stretch and byte scale (0-255) for 8bit raster. "
                        "Mode can be 'value' or 'percent' (cut-off units). "
                        "Values min and max define stretch borders (in mode units).",
                    },
                },
                "Slope gradient": {
                    "slp_compute": {
                        "value": self.slp_compute,
                        "description": "If compute Slope. Parameter for GUIs.",
                    },
                    "slp_output_units": {
                        "value": self.slp_output_units,
                        "description": "Slope output units [radian, degree, "
                        "percent].",
                    },
                    "slp_save_float": {
                        "value": self.slp_save_float,
                        "description": "If 1 it saves float raster, if 0 it doesn't.",
                    },
                    "slp_save_8bit": {
                        "value": self.slp_save_8bit,
                        "description": "If 1 it saves 8bit raster, if 0 it doesn't.",
                    },
                    "slp_bytscl": {
                        "mode": self.slp_bytscl[0],
                        "min": self.slp_bytscl[1],
                        "max": self.slp_bytscl[2],
                        "description": "Linear stretch and byte scale (0-255) for 8bit raster. "
                        "Mode can be 'value' or 'percent' (cut-off units). "
                        "Values min and max define stretch borders (in mode units).",
                    },
                },
                "Simple local relief model": {
                    "slrm_compute": {
                        "value": self.slrm_compute,
                        "description": "If compute Simple local relief model. "
                        "Parameter for GUIs.",
                    },
                    "slrm_rad_cell": {
                        "value": self.slrm_rad_cell,
                        "description": "Radius for trend assessment in pixels.",
                    },
                    "slrm_save_float": {
                        "value": self.slrm_save_float,
                        "description": "If 1 it saves float raster, if 0 it doesn't.",
                    },
                    "slrm_save_8bit": {
                        "value": self.slrm_save_8bit,
                        "description": "If 1 it saves 8bit raster, if 0 it doesn't.",
                    },
                    "slrm_bytscl": {
                        "mode": self.slrm_bytscl[0],
                        "min": self.slrm_bytscl[1],
                        "max": self.slrm_bytscl[2],
                        "description": "Linear stretch and byte scale (0-255) for 8bit raster. "
                        "Mode can be 'value' or 'percent' (cut-off units). "
                        "Values min and max define stretch borders (in mode units).",
                    },
                },
                "Multi-scale relief model": {
                    "msrm_compute": {
                        "value": self.msrm_compute,
                        "description": "If compute Multi-scale relief model. "
                        "Parameter for GUIs.",
                    },
                    "msrm_feature_min": {
                        "value": self.msrm_feature_min,
                        "description": "Minimum size of the feature you want to detect in meters.",
                    },
                    "msrm_feature_max": {
                        "value": self.msrm_feature_max,
                        "description": "Maximum size of the feature you want to detect in meters.",
                    },
                    "msrm_scaling_factor": {
                        "value": self.msrm_scaling_factor,
                        "description": "Scaling factor.",
                    },
                    "msrm_save_float": {
                        "value": self.msrm_save_float,
                        "description": "If 1 it saves float raster, if 0 it doesn't.",
                    },
                    "msrm_save_8bit": {
                        "value": self.msrm_save_8bit,
                        "description": "If 1 it saves 8bit raster, if 0 it doesn't.",
                    },
                    "msrm_bytscl": {
                        "mode": self.msrm_bytscl[0],
                        "min": self.msrm_bytscl[1],
                        "max": self.msrm_bytscl[2],
                        "description": "Linear stretch and byte scale (0-255) for 8bit raster. "
                        "Mode can be 'value' or 'percent' (cut-off units). "
                        "Values min and max define stretch borders (in mode units).",
                    },
                },
                "Sky-View Factor": {
                    "svf_compute": {
                        "value": self.svf_compute,
                        "description": "If compute Sky-View Factor."
                        " Parameter for GUIs.",
                    },
                    "svf_n_dir": {
                        "value": self.svf_n_dir,
                        "description": "Number of directions.",
                    },
                    "svf_r_max": {
                        "value": self.svf_r_max,
                        "description": "Maximal search " "radious in pixels.",
                    },
                    "svf_noise": {
                        "value": self.svf_noise,
                        "description": "The level of noise remove [0-don't remove, "
                        "1-low, 2-med, 3-high].",
                    },
                    "svf_save_float": {
                        "value": self.svf_save_float,
                        "description": "If 1 it saves float raster, if 0 it doesn't.",
                    },
                    "svf_save_8bit": {
                        "value": self.svf_save_8bit,
                        "description": "If 1 it saves 8bit raster, if 0 it doesn't.",
                    },
                    "svf_bytscl": {
                        "mode": self.svf_bytscl[0],
                        "min": self.svf_bytscl[1],
                        "max": self.svf_bytscl[2],
                        "description": "Linear stretch and byte scale (0-255) for 8bit raster. "
                        "Mode can be 'value' or 'percent' (cut-off units). "
                        "Values min and max define stretch borders (in mode units).",
                    },
                },
                "Anisotropic Sky-View Factor": {
                    "asvf_compute": {
                        "value": self.asvf_compute,
                        "description": "If compute Anisotropic Sky-View Factor."
                        " Parameter for GUIs.",
                    },
                    "asvf_dir": {
                        "value": self.asvf_dir,
                        "description": "Direction of anisotropy in degrees.",
                    },
                    "asvf_level": {
                        "value": self.asvf_level,
                        "description": "Level of anisotropy [1-low, 2-high].",
                    },
                    "asvf_bytscl": {
                        "mode": self.asvf_bytscl[0],
                        "min": self.asvf_bytscl[1],
                        "max": self.asvf_bytscl[2],
                        "description": "Linear stretch and byte scale (0-255) for 8bit raster. "
                        "Mode can be 'value' or 'percent' (cut-off units). "
                        "Values min and max define stretch borders (in mode units).",
                    },
                },
                "Openness - Positive": {
                    "pos_opns_compute": {
                        "value": self.pos_opns_compute,
                        "description": "If compute Openness - Positive. "
                        "Parameter for GUIs.",
                    },
                    "pos_opns_bytscl": {
                        "mode": self.pos_opns_bytscl[0],
                        "min": self.pos_opns_bytscl[1],
                        "max": self.pos_opns_bytscl[2],
                        "description": "Linear stretch and byte scale (0-255) for 8bit raster. "
                        "Mode can be 'value' or 'percent' (cut-off units). "
                        "Values min and max define stretch borders (in mode units).",
                    },
                },
                "Openness - Negative": {
                    "neg_opns_compute": {
                        "value": self.neg_opns_compute,
                        "description": "If compute Openness - Negative. "
                        "Parameter for GUIs.",
                    },
                    "neg_opns_save_float": {
                        "value": self.neg_opns_save_float,
                        "description": "If 1 it saves float raster, if 0 it doesn't.",
                    },
                    "neg_opns_save_8bit": {
                        "value": self.neg_opns_save_8bit,
                        "description": "If 1 it saves 8bit raster, if 0 it doesn't.",
                    },
                    "neg_opns_bytscl": {
                        "mode": self.neg_opns_bytscl[0],
                        "min": self.neg_opns_bytscl[1],
                        "max": self.neg_opns_bytscl[2],
                        "description": "Linear stretch and byte scale (0-255) for 8bit raster. "
                        "Mode can be 'value' or 'percent' (cut-off units). "
                        "Values min and max define stretch borders (in mode units).",
                    },
                },
                "Sky illumination": {
                    "sim_compute": {
                        "value": self.sim_compute,
                        "description": "If compute Sky illumination. Parameter for "
                        "GUIs.",
                    },
                    "sim_sky_mod": {
                        "value": self.sim_sky_mod,
                        "description": "Sky model [overcast, uniform].",
                    },
                    "sim_compute_shadow": {
                        "value": self.sim_compute_shadow,
                        "description": "If 1 it computes shadows, if 0 it doesn't.",
                    },
                    "sim_shadow_dist": {
                        "value": self.sim_shadow_dist,
                        "description": "Max shadow modeling distance in pixels.",
                    },
                    "sim_nr_dir": {
                        "value": self.sim_nr_dir,
                        "description": "Number of directions to search for horizon",
                    },
                    "sim_shadow_az": {
                        "value": self.sim_shadow_az,
                        "description": "Shadow " "azimuth in " "degrees.",
                    },
                    "sim_shadow_el": {
                        "value": self.sim_shadow_el,
                        "description": "Shadow " "elevation in " "degrees.",
                    },
                    "sim_save_float": {
                        "value": self.sim_save_float,
                        "description": "If 1 it saves float raster, if 0 it doesn't.",
                    },
                    "sim_save_8bit": {
                        "value": self.sim_save_8bit,
                        "description": "If 1 it saves 8bit raster, if 0 it doesn't.",
                    },
                    "sim_bytscl": {
                        "mode": self.sim_bytscl[0],
                        "min": self.sim_bytscl[1],
                        "max": self.sim_bytscl[2],
                        "description": "Linear stretch and byte scale (0-255) for 8bit raster. "
                        "Mode can be 'value' or 'percent' (cut-off units). "
                        "Values min and max define stretch borders (in mode units).",
                    },
                },
                "Local dominance": {
                    "ld_compute": {
                        "value": self.ld_compute,
                        "description": "If compute Local dominance. Parameter for "
                        "GUIs.",
                    },
                    "ld_min_rad": {
                        "value": self.ld_min_rad,
                        "description": "Minimum radial distance (in pixels) at which "
                        "the algorithm starts with visualization "
                        "computation.",
                    },
                    "ld_max_rad": {
                        "value": self.ld_max_rad,
                        "description": "Maximum radial distance (in pixels) at which "
                        "the algorithm ends with visualization "
                        "computation.",
                    },
                    "ld_rad_inc": {
                        "value": self.ld_rad_inc,
                        "description": "Radial distance " "steps in pixels.",
                    },
                    "ld_anglr_res": {
                        "value": self.ld_anglr_res,
                        "description": "Angular step for determination of number of "
                        "angular directions.",
                    },
                    "ld_observer_h": {
                        "value": self.ld_observer_h,
                        "description": "Height at which we observe the terrain.",
                    },
                    "ld_save_float": {
                        "value": self.ld_save_float,
                        "description": "If 1 it saves float raster, if 0 it doesn't.",
                    },
                    "ld_save_8bit": {
                        "value": self.ld_save_8bit,
                        "description": "If 1 it saves 8bit raster, if 0 it doesn't.",
                    },
                    "ld_bytscl": {
                        "mode": self.ld_bytscl[0],
                        "min": self.ld_bytscl[1],
                        "max": self.ld_bytscl[2],
                        "description": "Linear stretch and byte scale (0-255) for 8bit raster. "
                        "Mode can be 'value' or 'percent' (cut-off units). "
                        "Values min and max define stretch borders (in mode units).",
                    },
                },
                "Multi-scale topographic position": {
                    "mstp_compute": {
                        "value": self.mstp_compute,
                        "description": "If compute Multi-scale topographic position. "
                        "Parameter for GUIs.",
                    },
                    "mstp_local_scale": {
                        "min": self.mstp_local_scale[0],
                        "max": self.mstp_local_scale[1],
                        "step": self.mstp_local_scale[2],
                        "description": "Local scale minimum radius, maximum radius and step in pixels to"
                        " calculate maximum mean deviation from elevation."
                        " All have to be integers!",
                    },
                    "mstp_meso_scale": {
                        "min": self.mstp_meso_scale[0],
                        "max": self.mstp_meso_scale[1],
                        "step": self.mstp_meso_scale[2],
                        "description": "Meso scale minimum radius, maximum radius and step in pixels to"
                        " calculate maximum mean deviation from elevation."
                        " All have to be integers!",
                    },
                    "mstp_broad_scale": {
                        "min": self.mstp_broad_scale[0],
                        "max": self.mstp_broad_scale[1],
                        "step": self.mstp_broad_scale[2],
                        "description": "Broad scale minimum radius, maximum radius and step in pixels to"
                        " calculate maximum mean deviation from elevation."
                        " All have to be integers!",
                    },
                    "mstp_lightness": {
                        "value": self.mstp_lightness,
                        "description": "Lightness factor to adjust MSTP visibility.",
                    },
                    "mstp_save_float": {
                        "value": self.mstp_save_float,
                        "description": "If 1 it saves float raster, if 0 it doesn't.",
                    },
                    "mstp_save_8bit": {
                        "value": self.mstp_save_8bit,
                        "description": "If 1 it saves 8bit raster, if 0 it doesn't.",
                    },
                    "mstp_bytscl": {
                        "mode": self.mstp_bytscl[0],
                        "min": self.mstp_bytscl[1],
                        "max": self.mstp_bytscl[2],
                        "description": "Linear stretch and byte scale (0-255) for 8bit raster. "
                        "Mode can be 'value' or 'percent' (cut-off units). "
                        "Values min and max define stretch borders (in mode units).",
                    },
                },
            }
        }
        if file_path is None:
            file_path = r"settings\default_settings.json"
            if os.path.isfile(file_path):
                pass
            else:
                if not os.path.exists(os.path.dirname(file_path)):
                    os.makedirs(os.path.dirname(file_path))
        dat = open(file_path, "w")
        dat.write(json.dumps(data, indent=4))
        dat.close()

    def read_default_from_file(self, file_path):
        """Reads default attributes from file."""
        # Example file in dir settings: default_settings.json
        extension = os.path.splitext(file_path)[1]
        if extension == ".txt":
            dat = open(file_path, "r")
            remove_noise = 0
            for line in dat:
                line = line.strip()
                if line == "":
                    continue
                if line[0] == "#":
                    continue
                line_list = line.split("#")  # remove comments
                line_list = line_list[0].split("=")
                if len(line_list) == 2:
                    parameter_name = line_list[0].strip().lower()
                    parameter_value = line_list[1].strip()
                    if parameter_name == "overwrite":
                        self.overwrite = int(parameter_value)
                    elif parameter_name == "exaggeration_factor":
                        self.ve_factor = float(parameter_value)
                    elif parameter_name == "slope_gradient":
                        self.slp_compute = int(parameter_value)
                    elif parameter_name == "hillshading":
                        self.hs_compute = int(parameter_value)
                    elif parameter_name == "sun_azimuth":
                        self.hs_sun_azi = int(parameter_value)
                    elif parameter_name == "sun_elevation":
                        self.hs_sun_el = int(parameter_value)
                    elif parameter_name == "multiple_hillshading":
                        self.mhs_compute = int(parameter_value)
                    elif parameter_name == "hillshade_directions":
                        self.mhs_nr_dir = int(parameter_value)
                    elif parameter_name == "sun_elevation":
                        self.mhs_sun_el = int(parameter_value)
                    elif parameter_name == "simple_local_relief":
                        self.slrm_compute = int(parameter_value)
                    elif parameter_name == "trend_radius":
                        self.slrm_rad_cell = int(parameter_value)
                    elif parameter_name == "sky_view_factor":
                        self.svf_compute = int(parameter_value)
                    elif parameter_name == "svf_directions":
                        self.svf_n_dir = int(parameter_value)
                    elif parameter_name == "search_radius":
                        self.svf_r_max = int(parameter_value)
                    elif parameter_name == "remove_noise":
                        if int(parameter_value) == 0:
                            remove_noise = 0
                            self.svf_noise = 0
                        else:
                            remove_noise = 1
                    elif parameter_name == "noise_removal":
                        if remove_noise != 0:
                            if parameter_value == "low":
                                self.svf_noise = 1
                            elif parameter_value == "medium":
                                self.svf_noise = 2
                            elif parameter_value == "high":
                                self.svf_noise = 3
                    elif parameter_name == "anisotropic_svf":
                        self.asvf_compute = int(parameter_value)
                    elif parameter_name == "anisotropy_direction":
                        self.asvf_dir = int(parameter_value)
                    elif parameter_name == "anisotropy_level":
                        if parameter_value == "low":
                            self.asvf_level = 1
                        elif parameter_value == "high":
                            self.asvf_level = 2
                    elif parameter_name == "positive_openness":
                        self.pos_opns_compute = int(parameter_value)
                    elif parameter_name == "negative_openness":
                        self.neg_opns_compute = int(parameter_value)
                    elif parameter_name == "sky_illumination":
                        self.sim_compute = int(parameter_value)
                    elif parameter_name == "sky_model":
                        self.sim_sky_mod = parameter_value
                    elif parameter_name == "max_shadow_dist":
                        self.sim_shadow_dist = int(parameter_value)
                    elif parameter_name == "local_dominance":
                        self.ld_compute = int(parameter_value)
                    elif parameter_name == "min_radius":
                        self.ld_min_rad = int(parameter_value)
                    elif parameter_name == "max_radius":
                        self.ld_max_rad = int(parameter_value)
                    # else:
                    #     warnings.warn("rvt.default.read_default_from_file: Line '{}' not used.".format(line))
                else:
                    warnings.warn(
                        "rvt.default.read_default_from_file: Wrong line '{}'".format(
                            line
                        )
                    )
                    continue
            dat.close()
        elif extension == ".json":
            dat = open(file_path, "r")
            data = json.load(dat)
            default_data = data["default_settings"]
            self.overwrite = int(default_data["overwrite"]["value"])
            self.ve_factor = float(default_data["ve_factor"]["value"])
            # Slope gradient
            self.slp_compute = int(
                default_data["Slope gradient"]["slp_compute"]["value"]
            )
            self.slp_output_units = str(
                default_data["Slope gradient"]["slp_output_units"]["value"]
            )
            self.slp_save_float = int(
                default_data["Slope gradient"]["slp_save_float"]["value"]
            )
            self.slp_save_8bit = int(
                default_data["Slope gradient"]["slp_save_8bit"]["value"]
            )
            self.slp_bytscl = (
                str(default_data["Slope gradient"]["slp_bytscl"]["mode"]),
                float(default_data["Slope gradient"]["slp_bytscl"]["min"]),
                float(default_data["Slope gradient"]["slp_bytscl"]["max"]),
            )
            # Hillshade
            self.hs_compute = int(default_data["Hillshade"]["hs_compute"]["value"])
            self.hs_sun_azi = int(default_data["Hillshade"]["hs_sun_azi"]["value"])
            self.hs_sun_el = int(default_data["Hillshade"]["hs_sun_el"]["value"])
            self.hs_save_float = int(
                default_data["Hillshade"]["hs_save_float"]["value"]
            )
            self.hs_save_8bit = int(default_data["Hillshade"]["hs_save_8bit"]["value"])
            self.hs_shadow = int(default_data["Hillshade"]["hs_shadow"]["value"])
            self.hs_bytscl = (
                str(default_data["Hillshade"]["hs_bytscl"]["mode"]),
                float(default_data["Hillshade"]["hs_bytscl"]["min"]),
                float(default_data["Hillshade"]["hs_bytscl"]["max"]),
            )
            # Multiple directions hillshade
            self.mhs_compute = int(
                default_data["Multiple directions hillshade"]["mhs_compute"]["value"]
            )
            self.mhs_nr_dir = int(
                default_data["Multiple directions hillshade"]["mhs_nr_dir"]["value"]
            )
            self.mhs_sun_el = int(
                default_data["Multiple directions hillshade"]["mhs_sun_el"]["value"]
            )
            self.mhs_save_float = int(
                default_data["Multiple directions hillshade"]["mhs_save_float"]["value"]
            )
            self.mhs_save_8bit = int(
                default_data["Multiple directions hillshade"]["mhs_save_8bit"]["value"]
            )
            self.mhs_bytscl = (
                str(
                    default_data["Multiple directions hillshade"]["mhs_bytscl"]["mode"]
                ),
                float(
                    default_data["Multiple directions hillshade"]["mhs_bytscl"]["min"]
                ),
                float(
                    default_data["Multiple directions hillshade"]["mhs_bytscl"]["max"]
                ),
            )
            # Simple local relief model
            self.slrm_compute = int(
                default_data["Simple local relief model"]["slrm_compute"]["value"]
            )
            self.slrm_rad_cell = int(
                default_data["Simple local relief model"]["slrm_rad_cell"]["value"]
            )
            self.slrm_save_float = int(
                default_data["Simple local relief model"]["slrm_save_float"]["value"]
            )
            self.slrm_save_8bit = int(
                default_data["Simple local relief model"]["slrm_save_8bit"]["value"]
            )
            self.slrm_bytscl = (
                str(default_data["Simple local relief model"]["slrm_bytscl"]["mode"]),
                float(default_data["Simple local relief model"]["slrm_bytscl"]["min"]),
                float(default_data["Simple local relief model"]["slrm_bytscl"]["max"]),
            )
            # Mulit-scale relief model
            self.msrm_compute = int(
                default_data["Multi-scale relief model"]["msrm_compute"]["value"]
            )
            self.msrm_feature_min = float(
                default_data["Multi-scale relief model"]["msrm_feature_min"]["value"]
            )
            self.msrm_feature_max = float(
                default_data["Multi-scale relief model"]["msrm_feature_max"]["value"]
            )
            self.msrm_scaling_factor = int(
                default_data["Multi-scale relief model"]["msrm_scaling_factor"]["value"]
            )
            self.msrm_save_float = int(
                default_data["Multi-scale relief model"]["msrm_save_float"]["value"]
            )
            self.msrm_save_8bit = int(
                default_data["Multi-scale relief model"]["msrm_save_8bit"]["value"]
            )
            self.msrm_bytscl = (
                str(default_data["Multi-scale relief model"]["msrm_bytscl"]["mode"]),
                float(default_data["Multi-scale relief model"]["msrm_bytscl"]["min"]),
                float(default_data["Multi-scale relief model"]["msrm_bytscl"]["max"]),
            )
            # Sky-View Factor
            self.svf_compute = int(
                default_data["Sky-View Factor"]["svf_compute"]["value"]
            )
            self.svf_n_dir = int(default_data["Sky-View Factor"]["svf_n_dir"]["value"])
            self.svf_r_max = int(default_data["Sky-View Factor"]["svf_r_max"]["value"])
            self.svf_noise = int(default_data["Sky-View Factor"]["svf_noise"]["value"])
            self.svf_save_float = int(
                default_data["Sky-View Factor"]["svf_save_float"]["value"]
            )
            self.svf_save_8bit = int(
                default_data["Sky-View Factor"]["svf_save_8bit"]["value"]
            )
            self.svf_bytscl = (
                str(default_data["Sky-View Factor"]["svf_bytscl"]["mode"]),
                float(default_data["Sky-View Factor"]["svf_bytscl"]["min"]),
                float(default_data["Sky-View Factor"]["svf_bytscl"]["max"]),
            )
            # Anisotropic Sky-View Factor
            self.asvf_compute = int(
                default_data["Anisotropic Sky-View Factor"]["asvf_compute"]["value"]
            )
            self.asvf_dir = int(
                default_data["Anisotropic Sky-View Factor"]["asvf_dir"]["value"]
            )
            self.asvf_level = int(
                default_data["Anisotropic Sky-View Factor"]["asvf_level"]["value"]
            )
            self.asvf_bytscl = (
                str(default_data["Anisotropic Sky-View Factor"]["asvf_bytscl"]["mode"]),
                float(
                    default_data["Anisotropic Sky-View Factor"]["asvf_bytscl"]["min"]
                ),
                float(
                    default_data["Anisotropic Sky-View Factor"]["asvf_bytscl"]["max"]
                ),
            )
            # Openness - Positive
            self.pos_opns_compute = int(
                default_data["Openness - Positive"]["pos_opns_compute"]["value"]
            )
            self.pos_opns_bytscl = (
                str(default_data["Openness - Positive"]["pos_opns_bytscl"]["mode"]),
                float(default_data["Openness - Positive"]["pos_opns_bytscl"]["min"]),
                float(default_data["Openness - Positive"]["pos_opns_bytscl"]["max"]),
            )
            # Openness - Negative
            self.neg_opns_compute = int(
                default_data["Openness - Negative"]["neg_opns_compute"]["value"]
            )
            self.neg_opns_save_float = int(
                default_data["Openness - Negative"]["neg_opns_save_float"]["value"]
            )
            self.neg_opns_save_8bit = int(
                default_data["Openness - Negative"]["neg_opns_save_8bit"]["value"]
            )
            self.neg_opns_bytscl = (
                str(default_data["Openness - Negative"]["neg_opns_bytscl"]["mode"]),
                float(default_data["Openness - Negative"]["neg_opns_bytscl"]["min"]),
                float(default_data["Openness - Negative"]["neg_opns_bytscl"]["max"]),
            )
            # Sky illumination
            self.sim_compute = int(
                default_data["Sky illumination"]["sim_compute"]["value"]
            )
            self.sim_sky_mod = str(
                default_data["Sky illumination"]["sim_sky_mod"]["value"]
            )
            self.sim_compute_shadow = int(
                default_data["Sky illumination"]["sim_compute_shadow"]["value"]
            )
            self.sim_nr_dir = int(
                default_data["Sky illumination"]["sim_nr_dir"]["value"]
            )
            self.sim_shadow_dist = int(
                default_data["Sky illumination"]["sim_shadow_dist"]["value"]
            )
            self.sim_shadow_az = int(
                default_data["Sky illumination"]["sim_shadow_az"]["value"]
            )
            self.sim_shadow_el = int(
                default_data["Sky illumination"]["sim_shadow_el"]["value"]
            )
            self.sim_save_float = int(
                default_data["Sky illumination"]["sim_save_float"]["value"]
            )
            self.sim_save_8bit = int(
                default_data["Sky illumination"]["sim_save_8bit"]["value"]
            )
            self.sim_bytscl = (
                str(default_data["Sky illumination"]["sim_bytscl"]["mode"]),
                float(default_data["Sky illumination"]["sim_bytscl"]["min"]),
                float(default_data["Sky illumination"]["sim_bytscl"]["max"]),
            )
            # Local dominance
            self.ld_compute = int(
                default_data["Local dominance"]["ld_compute"]["value"]
            )
            self.ld_min_rad = int(
                default_data["Local dominance"]["ld_min_rad"]["value"]
            )
            self.ld_max_rad = int(
                default_data["Local dominance"]["ld_max_rad"]["value"]
            )
            self.ld_rad_inc = int(
                default_data["Local dominance"]["ld_rad_inc"]["value"]
            )
            self.ld_anglr_res = int(
                default_data["Local dominance"]["ld_anglr_res"]["value"]
            )
            self.ld_observer_h = float(
                default_data["Local dominance"]["ld_observer_h"]["value"]
            )
            self.ld_save_float = int(
                default_data["Local dominance"]["ld_save_float"]["value"]
            )
            self.ld_save_8bit = int(
                default_data["Local dominance"]["ld_save_8bit"]["value"]
            )
            self.ld_bytscl = (
                str(default_data["Local dominance"]["ld_bytscl"]["mode"]),
                float(default_data["Local dominance"]["ld_bytscl"]["min"]),
                float(default_data["Local dominance"]["ld_bytscl"]["max"]),
            )
            # Multi-scale topographic position
            self.mstp_compute = int(
                default_data["Multi-scale topographic position"]["mstp_compute"][
                    "value"
                ]
            )
            self.mstp_local_scale = (
                int(
                    default_data["Multi-scale topographic position"][
                        "mstp_local_scale"
                    ]["min"]
                ),
                int(
                    default_data["Multi-scale topographic position"][
                        "mstp_local_scale"
                    ]["max"]
                ),
                int(
                    default_data["Multi-scale topographic position"][
                        "mstp_local_scale"
                    ]["step"]
                ),
            )
            self.mstp_meso_scale = (
                int(
                    default_data["Multi-scale topographic position"]["mstp_meso_scale"][
                        "min"
                    ]
                ),
                int(
                    default_data["Multi-scale topographic position"]["mstp_meso_scale"][
                        "max"
                    ]
                ),
                int(
                    default_data["Multi-scale topographic position"]["mstp_meso_scale"][
                        "step"
                    ]
                ),
            )
            self.mstp_broad_scale = (
                int(
                    default_data["Multi-scale topographic position"][
                        "mstp_broad_scale"
                    ]["min"]
                ),
                int(
                    default_data["Multi-scale topographic position"][
                        "mstp_broad_scale"
                    ]["max"]
                ),
                int(
                    default_data["Multi-scale topographic position"][
                        "mstp_broad_scale"
                    ]["step"]
                ),
            )
            self.mstp_lightness = float(
                default_data["Multi-scale topographic position"]["mstp_lightness"][
                    "value"
                ]
            )
            self.mstp_save_float = int(
                default_data["Multi-scale topographic position"]["mstp_save_float"][
                    "value"
                ]
            )
            self.mstp_save_8bit = int(
                default_data["Multi-scale topographic position"]["mstp_save_8bit"][
                    "value"
                ]
            )
            self.mstp_bytscl = (
                str(
                    default_data["Multi-scale topographic position"]["mstp_bytscl"][
                        "mode"
                    ]
                ),
                float(
                    default_data["Multi-scale topographic position"]["mstp_bytscl"][
                        "min"
                    ]
                ),
                float(
                    default_data["Multi-scale topographic position"]["mstp_bytscl"][
                        "max"
                    ]
                ),
            )
            dat.close()

    def get_shadow_file_name(self, dem_path):
        """Returns shadow name, with added hillshade parameters (hs_sun_azi == shadow azimuth,
        hs_sun_el == shadow_elevation)."""
        dem_name = os.path.basename(dem_path).split(".")[
            0
        ]  # base name without extension
        return "{}_shadow_A{}_H{}.tif".format(dem_name, self.hs_sun_azi, self.hs_sun_el)

    def get_shadow_path(self, dem_path):
        """Returns path to Shadow. Generates shadow name (uses default attributes and dem name from dem_path) and
        adds dem directory (dem_path) to it."""
        return os.path.normpath(
            os.path.join(os.path.dirname(dem_path), self.get_shadow_file_name(dem_path))
        )

    def get_hillshade_file_name(self, dem_path, bit8=False):
        """Returns Hillshade name, dem name (from dem_path) with added hillshade parameters.
        If bit8 it returns 8bit file name."""
        dem_name = os.path.basename(dem_path).split(".")[
            0
        ]  # base name without extension
        if bit8:
            return "{}_HS_A{}_H{}_8bit.tif".format(
                dem_name, self.hs_sun_azi, self.hs_sun_el
            )
        else:
            return "{}_HS_A{}_H{}.tif".format(dem_name, self.hs_sun_azi, self.hs_sun_el)

    def get_hillshade_path(self, dem_path, bit8=False):
        """Returns path to Hillshade. Generates hillshade name (uses default attributes and dem name from dem_path) and
        adds dem directory (dem_path) to it. If bit8 it returns 8bit file path."""
        return os.path.normpath(
            os.path.join(
                os.path.dirname(dem_path), self.get_hillshade_file_name(dem_path, bit8)
            )
        )

    def get_slope_file_name(self, dem_path, bit8=False):
        """Returns Slope name, dem name (from dem_path) with added slope parameters.
        If bit8 it returns 8bit file name."""
        dem_name = os.path.basename(dem_path).split(".")[
            0
        ]  # base name without extension
        if bit8:
            return "{}_SLOPE_8bit.tif".format(dem_name)
        else:
            return "{}_SLOPE.tif".format(dem_name)

    def get_slope_path(self, dem_path, bit8=False):
        """Returns path to slope. Generates slope name and adds dem directory (dem_path) to it.
        If bit8 it returns 8bit file path."""
        return os.path.normpath(
            os.path.join(
                os.path.dirname(dem_path), self.get_slope_file_name(dem_path, bit8)
            )
        )

    def get_multi_hillshade_file_name(self, dem_path, bit8=False):
        """Returns Multiple directions hillshade name, dem name (from dem_path) with added
        multi hillshade parameters. If bit8 it returns 8bit file name."""
        dem_name = os.path.basename(dem_path).split(".")[
            0
        ]  # base name without extension
        if bit8:
            return "{}_MULTI-HS_D{}_H{}_8bit.tif".format(
                dem_name, self.mhs_nr_dir, self.mhs_sun_el
            )
        else:
            return "{}_MULTI-HS_D{}_H{}.tif".format(
                dem_name, self.mhs_nr_dir, self.mhs_sun_el
            )

    def get_multi_hillshade_path(self, dem_path, bit8=False):
        """Returns path to Multiple directions hillshade. Generates multi hillshade name (uses default attributes and
        dem name from dem_path) and adds dem directory (dem_path) to it. If bit8 it returns 8bit file path."""
        return os.path.normpath(
            os.path.join(
                os.path.dirname(dem_path),
                self.get_multi_hillshade_file_name(dem_path, bit8),
            )
        )

    def get_slrm_file_name(self, dem_path, bit8=False):
        """Returns Simple local relief model name, dem name (from dem_path) with added slrm parameters.
        If bit8 it returns 8bit file name."""
        dem_name = os.path.basename(dem_path).split(".")[
            0
        ]  # base name without extension
        if bit8:
            return "{}_SLRM_R{}_8bit.tif".format(dem_name, self.slrm_rad_cell)
        else:
            return "{}_SLRM_R{}.tif".format(dem_name, self.slrm_rad_cell)

    def get_slrm_path(self, dem_path, bit8=False):
        """Returns path to Simple local relief model. Generates slrm name (uses default attributes and dem name from
        dem_path) and adds dem directory (dem_path) to it. If bit8 it returns 8bit file path."""
        return os.path.normpath(
            os.path.join(
                os.path.dirname(dem_path), self.get_slrm_file_name(dem_path, bit8)
            )
        )

    def get_svf_file_name(self, dem_path, bit8=False):
        """Returns Sky-view factor name, dem name (from dem_path) with added svf parameters.
        If bit8 it returns 8bit file name."""
        dem_name = os.path.basename(dem_path).split(".")[
            0
        ]  # base name without extension
        out_name = "{}_SVF_R{}_D{}".format(dem_name, self.svf_r_max, self.svf_n_dir)
        if self.svf_noise == 1:
            out_name += "_NRlow"
        elif self.svf_noise == 2:
            out_name += "_NRmedium"
        elif self.svf_noise == 3:
            out_name += "_NRhigh"
        if bit8:
            out_name += "_8bit"
        return out_name + ".tif"

    def get_svf_path(self, dem_path, bit8=False):
        """Returns path to Sky-view factor. Generates svf name (uses default attributes and dem name from
        dem_path) and adds dem directory (dem_path) to it. If bit8 it returns 8bit file path."""
        return os.path.normpath(
            os.path.join(
                os.path.dirname(dem_path), self.get_svf_file_name(dem_path, bit8)
            )
        )

    def get_asvf_file_name(self, dem_path, bit8=False):
        """Returns Anisotropic Sky-view factor name, dem name (from dem_path) with added asvf parameters.
        If bit8 it returns 8bit file name."""
        dem_name = os.path.basename(dem_path).split(".")[
            0
        ]  # base name without extension
        out_name = "{}_SVF-A_R{}_D{}_A{}".format(
            dem_name, self.svf_r_max, self.svf_n_dir, self.asvf_dir
        )
        if self.asvf_level == 1:
            out_name += "_ALlow"
        elif self.asvf_level == 2:
            out_name += "_ALhigh"
        if self.svf_noise == 1:
            out_name += "_NRlow"
        elif self.svf_noise == 2:
            out_name += "_NRmedium"
        elif self.svf_noise == 3:
            out_name += "_NRhigh"
        if bit8:
            out_name += "_8bit"
        return out_name + ".tif"

    def get_asvf_path(self, dem_path, bit8=False):
        """Returns path to Anisotropic Sky-view factor. Generates asvf name (uses default attributes and dem name from
        dem_path) and adds dem directory (dem_path) to it. If bit8 it returns 8bit file path."""
        return os.path.normpath(
            os.path.join(
                os.path.dirname(dem_path), self.get_asvf_file_name(dem_path, bit8)
            )
        )

    def get_opns_file_name(self, dem_path, bit8=False):
        """Returns Positive Openness name, dem name (from dem_path) with added pos opns parameters.
        If bit8 it returns 8bit file name."""
        dem_name = os.path.basename(dem_path).split(".")[
            0
        ]  # base name without extension
        out_name = "{}_OPEN-POS_R{}_D{}".format(
            dem_name, self.svf_r_max, self.svf_n_dir
        )
        if self.svf_noise == 1:
            out_name += "_NRlow"
        elif self.svf_noise == 2:
            out_name += "_NRmedium"
        elif self.svf_noise == 3:
            out_name += "_NRhigh"
        if bit8:
            out_name += "_8bit"
        return out_name + ".tif"

    def get_opns_path(self, dem_path, bit8=False):
        """Returns path to Positive Openness. Generates pos opns name (uses default attributes and dem name from
        dem_path) and adds dem directory (dem_path) to it. If bit8 it returns 8bit file path."""
        return os.path.normpath(
            os.path.join(
                os.path.dirname(dem_path), self.get_opns_file_name(dem_path, bit8)
            )
        )

    def get_neg_opns_file_name(self, dem_path, bit8=False):
        """Returns Negative Openness name, dem name (from dem_path) with added neg opns parameters.
        If bit8 it returns 8bit file name."""
        dem_name = os.path.basename(dem_path).split(".")[
            0
        ]  # base name without extension
        out_name = "{}_OPEN-NEG_R{}_D{}".format(
            dem_name, self.svf_r_max, self.svf_n_dir
        )
        if self.svf_noise == 1:
            out_name += "_NRlow"
        elif self.svf_noise == 2:
            out_name += "_NRmedium"
        elif self.svf_noise == 3:
            out_name += "_NRstrong"
        if bit8:
            out_name += "_8bit"
        return out_name + ".tif"

    def get_neg_opns_path(self, dem_path, bit8=False):
        """Returns path to Negative Openness. Generates pos neg name (uses default attributes and dem name from
        dem_path) and adds dem directory (dem_path) to it. If bit8 it returns 8bit file path."""
        return os.path.normpath(
            os.path.join(
                os.path.dirname(dem_path), self.get_neg_opns_file_name(dem_path, bit8)
            )
        )

    def get_sky_illumination_file_name(self, dem_path, bit8=False):
        """Returns Sky illumination name, dem name (from dem_path) with added sim parameters.
        If bit8 it returns 8bit file name."""
        dem_name = os.path.basename(dem_path).split(".")[
            0
        ]  # base name without extension
        if bit8:
            return "{}_SIM_{}_D{}_{}px_8bit.tif".format(
                dem_name, self.sim_sky_mod, self.sim_nr_dir, self.sim_shadow_dist
            )
        else:
            return "{}_SIM_{}_D{}_{}px.tif".format(
                dem_name, self.sim_sky_mod, self.sim_nr_dir, self.sim_shadow_dist
            )

    def get_sky_illumination_path(self, dem_path, bit8=False):
        """Returns path to Sky illumination. Generates sim name (uses default attributes and dem name from
        dem_path) and adds dem directory (dem_path) to it. If bit8 it returns 8bit file path."""
        return os.path.normpath(
            os.path.join(
                os.path.dirname(dem_path),
                self.get_sky_illumination_file_name(dem_path, bit8),
            )
        )

    def get_local_dominance_file_name(self, dem_path, bit8=False):
        """Returns Local dominance name, dem name (from dem_path) with added ld parameters.
        If bit8 it returns 8bit file name."""
        dem_name = os.path.basename(dem_path).split(".")[
            0
        ]  # base name without extension
        if bit8:
            return "{}_LD_R_M{}-{}_DI{}_A{}_OH{}_8bit.tif".format(
                dem_name,
                self.ld_min_rad,
                self.ld_max_rad,
                self.ld_rad_inc,
                self.ld_anglr_res,
                self.ld_observer_h,
            )
        else:
            return "{}_LD_R_M{}-{}_DI{}_A{}_OH{}.tif".format(
                dem_name,
                self.ld_min_rad,
                self.ld_max_rad,
                self.ld_rad_inc,
                self.ld_anglr_res,
                self.ld_observer_h,
            )

    def get_local_dominance_path(self, dem_path, bit8=False):
        """Returns path to Local dominance. Generates ld name (uses default attributes and dem name from
        dem_path) and adds dem directory (dem_path) to it. If bit8 it returns 8bit file path."""
        return os.path.normpath(
            os.path.join(
                os.path.dirname(dem_path),
                self.get_local_dominance_file_name(dem_path, bit8),
            )
        )

    def get_msrm_file_name(self, dem_path, bit8=False):
        """Returns Multi-scale relief model name, dem name (from dem_path) with added msrm parameters.
        If bit8 it returns 8bit file name."""
        dem_name = os.path.basename(dem_path).split(".")[
            0
        ]  # base name without extension
        if bit8:
            return "{}_MSRM_F_M{}-{}_S{}_8bit.tif".format(
                dem_name,
                self.msrm_feature_min,
                self.msrm_feature_max,
                self.msrm_scaling_factor,
            )
        else:
            return "{}_MSRM_F_M{}-{}_S{}.tif".format(
                dem_name,
                self.msrm_feature_min,
                self.msrm_feature_max,
                self.msrm_scaling_factor,
            )

    def get_msrm_path(self, dem_path, bit8=False):
        """Returns path to Multi-scale relief model. Generates msrm name (uses default attributes and dem name from
        dem_path) and adds dem directory (dem_path) to it. If bit8 it returns 8bit file path."""
        return os.path.normpath(
            os.path.join(
                os.path.dirname(dem_path), self.get_msrm_file_name(dem_path, bit8)
            )
        )

    def get_mstp_file_name(self, dem_path, bit8=False):
        """Returns Multi-scale topographic position name, dem name (from dem_path) with added mstp parameters.
        If bit8 it returns 8bit file name."""
        dem_name = os.path.basename(dem_path).split(".")[
            0
        ]  # base name without extension

        mstp_file_name = "{}_MSTP_{}_{}_{}_L{}.tif".format(
            dem_name,
            self.mstp_local_scale[1],
            self.mstp_meso_scale[1],
            self.mstp_broad_scale[1],
            self.mstp_lightness,
        )

        if bit8:
            mstp_file_name = mstp_file_name.replace(".tif", "_8bit.tif")

        return mstp_file_name

    def get_mstp_path(self, dem_path, bit8=False):
        return os.path.normpath(
            os.path.join(
                os.path.dirname(dem_path), self.get_mstp_file_name(dem_path, bit8)
            )
        )

    def get_visualization_file_name(
        self, rvt_visualization: RVTVisualization, dem_path: Path, path_8bit: bool
    ) -> str:
        """ "Return visualization path."""
        if rvt_visualization == rvt.default.RVTVisualization.SLOPE:
            return self.get_slope_file_name(dem_path=dem_path, bit8=path_8bit)
        elif rvt_visualization == rvt.default.RVTVisualization.HILLSHADE:
            return self.get_hillshade_file_name(dem_path=dem_path, bit8=path_8bit)
        elif rvt_visualization == rvt.default.RVTVisualization.SHADOW:
            return self.get_shadow_file_name(dem_path=dem_path)
        elif rvt_visualization == rvt.default.RVTVisualization.MULTI_HILLSHADE:
            return self.get_multi_hillshade_file_name(dem_path=dem_path, bit8=path_8bit)
        elif (
            rvt_visualization == rvt.default.RVTVisualization.SIMPLE_LOCAL_RELIEF_MODEL
        ):
            return self.get_slrm_file_name(dem_path=dem_path, bit8=path_8bit)
        elif rvt_visualization == rvt.default.RVTVisualization.SKY_VIEW_FACTOR:
            return self.get_svf_file_name(dem_path=dem_path, bit8=path_8bit)
        elif (
            rvt_visualization
            == rvt.default.RVTVisualization.ANISOTROPIC_SKY_VIEW_FACTOR
        ):
            return self.get_asvf_file_name(dem_path=dem_path, bit8=path_8bit)
        elif rvt_visualization == rvt.default.RVTVisualization.POSITIVE_OPENNESS:
            return self.get_opns_file_name(dem_path=dem_path, bit8=path_8bit)
        elif rvt_visualization == rvt.default.RVTVisualization.NEGATIVE_OPENNESS:
            return self.get_neg_opns_file_name(dem_path=dem_path, bit8=path_8bit)
        elif rvt_visualization == rvt.default.RVTVisualization.SKY_ILLUMINATION:
            return self.get_sky_illumination_file_name(
                dem_path=dem_path, bit8=path_8bit
            )
        elif rvt_visualization == rvt.default.RVTVisualization.LOCAL_DOMINANCE:
            return self.get_local_dominance_file_name(dem_path=dem_path, bit8=path_8bit)
        elif rvt_visualization == rvt.default.RVTVisualization.MULTI_SCALE_RELIEF_MODEL:
            return self.get_msrm_file_name(dem_path=dem_path, bit8=path_8bit)
        elif (
            rvt_visualization
            == rvt.default.RVTVisualization.MULTI_SCALE_TOPOGRAPHIC_POSITION
        ):
            return self.get_mstp_file_name(dem_path=dem_path, bit8=path_8bit)

    def get_visualization_path(
        self,
        rvt_visualization: RVTVisualization,
        dem_path: Path,
        output_dir_path: Path,
        path_8bit: bool,
    ) -> Path:
        """ "Return visualization path."""
        if rvt_visualization == rvt.default.RVTVisualization.SLOPE:
            return output_dir_path / Path(
                self.get_slope_file_name(dem_path=dem_path, bit8=path_8bit)
            )
        elif rvt_visualization == rvt.default.RVTVisualization.HILLSHADE:
            return output_dir_path / Path(
                self.get_hillshade_file_name(dem_path=dem_path, bit8=path_8bit)
            )
        elif rvt_visualization == rvt.default.RVTVisualization.SHADOW:
            return output_dir_path / Path(self.get_shadow_file_name(dem_path=dem_path))
        elif rvt_visualization == rvt.default.RVTVisualization.MULTI_HILLSHADE:
            return output_dir_path / Path(
                self.get_multi_hillshade_file_name(dem_path=dem_path, bit8=path_8bit)
            )
        elif (
            rvt_visualization == rvt.default.RVTVisualization.SIMPLE_LOCAL_RELIEF_MODEL
        ):
            return output_dir_path / Path(
                self.get_slrm_file_name(dem_path=dem_path, bit8=path_8bit)
            )
        elif rvt_visualization == rvt.default.RVTVisualization.SKY_VIEW_FACTOR:
            return output_dir_path / Path(
                self.get_svf_file_name(dem_path=dem_path, bit8=path_8bit)
            )
        elif (
            rvt_visualization
            == rvt.default.RVTVisualization.ANISOTROPIC_SKY_VIEW_FACTOR
        ):
            return output_dir_path / Path(
                self.get_asvf_file_name(dem_path=dem_path, bit8=path_8bit)
            )
        elif rvt_visualization == rvt.default.RVTVisualization.POSITIVE_OPENNESS:
            return output_dir_path / Path(
                self.get_opns_file_name(dem_path=dem_path, bit8=path_8bit)
            )
        elif rvt_visualization == rvt.default.RVTVisualization.NEGATIVE_OPENNESS:
            return output_dir_path / Path(
                self.get_neg_opns_file_name(dem_path=dem_path, bit8=path_8bit)
            )
        elif rvt_visualization == rvt.default.RVTVisualization.SKY_ILLUMINATION:
            return output_dir_path / Path(
                self.get_sky_illumination_file_name(dem_path=dem_path, bit8=path_8bit)
            )
        elif rvt_visualization == rvt.default.RVTVisualization.LOCAL_DOMINANCE:
            return output_dir_path / Path(
                self.get_local_dominance_file_name(dem_path=dem_path, bit8=path_8bit)
            )
        elif rvt_visualization == rvt.default.RVTVisualization.MULTI_SCALE_RELIEF_MODEL:
            return output_dir_path / Path(
                self.get_msrm_file_name(dem_path=dem_path, bit8=path_8bit)
            )
        elif (
            rvt_visualization
            == rvt.default.RVTVisualization.MULTI_SCALE_TOPOGRAPHIC_POSITION
        ):
            return output_dir_path / Path(
                self.get_mstp_file_name(dem_path=dem_path, bit8=path_8bit)
            )

    def float_to_8bit(
        self,
        float_arr: np.array,
        visualization: RVTVisualization,
        x_res: float = None,
        y_res: float = None,
        no_data: Optional[float] = None,
    ):
        """Converts (byte scale) float visualization to 8bit. Resolution (x_res, y_res) and no_data needed only for
        multiple directions hillshade! Method first normalize then byte scale (0-255)."""
        if visualization == RVTVisualization.HILLSHADE:
            norm_arr = rvt.blend_func.normalize_image(
                visualization="hs",
                image=float_arr,
                min_norm=self.hs_bytscl[1],
                max_norm=self.hs_bytscl[2],
                normalization=self.hs_bytscl[0],
            )
            return rvt.vis.byte_scale(data=norm_arr, no_data=np.nan, c_min=0, c_max=1)
        elif visualization == RVTVisualization.SLOPE:
            norm_arr = rvt.blend_func.normalize_image(
                visualization="slp",
                image=float_arr,
                min_norm=self.slp_bytscl[1],
                max_norm=self.slp_bytscl[2],
                normalization=self.slp_bytscl[0],
            )
            return rvt.vis.byte_scale(data=norm_arr, no_data=np.nan, c_min=0, c_max=1)
        elif visualization == RVTVisualization.SHADOW:
            return float_arr
        elif visualization == RVTVisualization.MULTI_HILLSHADE:
            # Be careful when multihillshade we input dem, because we have to calculate hillshade in 3 directions
            red_band_arr = rvt.vis.hillshade(
                dem=float_arr,
                resolution_x=x_res,
                resolution_y=y_res,
                sun_elevation=self.mhs_sun_el,
                sun_azimuth=315,
                no_data=no_data,
            )
            green_band_arr = rvt.vis.hillshade(
                dem=float_arr,
                resolution_x=x_res,
                resolution_y=y_res,
                sun_elevation=self.mhs_sun_el,
                sun_azimuth=22.5,
                no_data=no_data,
            )
            blue_band_arr = rvt.vis.hillshade(
                dem=float_arr,
                resolution_x=x_res,
                resolution_y=y_res,
                sun_elevation=self.mhs_sun_el,
                sun_azimuth=90,
                no_data=no_data,
            )
            if (
                self.mhs_bytscl[0].lower() == "percent"
                or self.slp_bytscl[0].lower() == "perc"
            ):
                red_band_arr = rvt.blend_func.normalize_perc(
                    image=red_band_arr,
                    minimum=self.mhs_bytscl[1],
                    maximum=self.mhs_bytscl[2],
                )
                red_band_arr = rvt.vis.byte_scale(
                    data=red_band_arr, no_data=np.nan, c_min=0, c_max=1
                )
                green_band_arr = rvt.blend_func.normalize_perc(
                    image=green_band_arr,
                    minimum=self.mhs_bytscl[1],
                    maximum=self.mhs_bytscl[2],
                )
                green_band_arr = rvt.vis.byte_scale(
                    data=green_band_arr, no_data=np.nan, c_min=0, c_max=1
                )
                blue_band_arr = rvt.blend_func.normalize_perc(
                    image=blue_band_arr,
                    minimum=self.mhs_bytscl[1],
                    maximum=self.mhs_bytscl[2],
                )
                blue_band_arr = rvt.vis.byte_scale(
                    data=blue_band_arr, no_data=np.nan, c_min=0, c_max=1
                )
            else:  # self.mhs_bytscl[0] == "value"
                red_band_arr = rvt.blend_func.normalize_lin(
                    image=red_band_arr,
                    minimum=self.mhs_bytscl[1],
                    maximum=self.mhs_bytscl[2],
                )
                red_band_arr = rvt.vis.byte_scale(
                    data=red_band_arr, no_data=np.nan, c_min=0, c_max=1
                )
                green_band_arr = rvt.blend_func.normalize_lin(
                    image=green_band_arr,
                    minimum=self.mhs_bytscl[1],
                    maximum=self.mhs_bytscl[2],
                )
                green_band_arr = rvt.vis.byte_scale(
                    data=green_band_arr, no_data=np.nan, c_min=0, c_max=1
                )
                blue_band_arr = rvt.blend_func.normalize_lin(
                    image=blue_band_arr,
                    minimum=self.mhs_bytscl[1],
                    maximum=self.mhs_bytscl[2],
                )
                blue_band_arr = rvt.vis.byte_scale(
                    data=blue_band_arr, no_data=np.nan, c_min=0, c_max=1
                )
            multi_hillshade_8bit_arr = np.array(
                [red_band_arr, green_band_arr, blue_band_arr]
            )
            return multi_hillshade_8bit_arr
        elif visualization == RVTVisualization.SIMPLE_LOCAL_RELIEF_MODEL:
            norm_arr = rvt.blend_func.normalize_image(
                visualization="slrm",
                image=float_arr,
                min_norm=self.slrm_bytscl[1],
                max_norm=self.slrm_bytscl[2],
                normalization=self.slrm_bytscl[0],
            )
            return rvt.vis.byte_scale(data=norm_arr, no_data=np.nan, c_min=0, c_max=1)
        elif visualization == RVTVisualization.SKY_VIEW_FACTOR:
            norm_arr = rvt.blend_func.normalize_image(
                visualization="svf",
                image=float_arr,
                min_norm=self.svf_bytscl[1],
                max_norm=self.svf_bytscl[2],
                normalization=self.svf_bytscl[0],
            )
            return rvt.vis.byte_scale(data=norm_arr, no_data=np.nan, c_min=0, c_max=1)
        elif visualization == RVTVisualization.ANISOTROPIC_SKY_VIEW_FACTOR:
            norm_arr = rvt.blend_func.normalize_image(
                visualization="asvf",
                image=float_arr,
                min_norm=self.asvf_bytscl[1],
                max_norm=self.asvf_bytscl[2],
                normalization=self.asvf_bytscl[0],
            )
            return rvt.vis.byte_scale(data=norm_arr, no_data=np.nan, c_min=0, c_max=1)
        elif visualization == RVTVisualization.POSITIVE_OPENNESS:
            norm_arr = rvt.blend_func.normalize_image(
                visualization="pos_opns",
                image=float_arr,
                min_norm=self.pos_opns_bytscl[1],
                max_norm=self.pos_opns_bytscl[2],
                normalization=self.pos_opns_bytscl[0],
            )
            return rvt.vis.byte_scale(data=norm_arr, no_data=np.nan, c_min=0, c_max=1)
        elif visualization == RVTVisualization.NEGATIVE_OPENNESS:
            norm_arr = rvt.blend_func.normalize_image(
                visualization="neg_opns",
                image=float_arr,
                min_norm=self.neg_opns_bytscl[1],
                max_norm=self.neg_opns_bytscl[2],
                normalization=self.neg_opns_bytscl[0],
            )
            return rvt.vis.byte_scale(data=norm_arr, no_data=np.nan, c_min=0, c_max=1)
        elif visualization == RVTVisualization.SKY_ILLUMINATION:
            norm_arr = rvt.blend_func.normalize_image(
                visualization="sim",
                image=float_arr,
                min_norm=self.sim_bytscl[1],
                max_norm=self.sim_bytscl[2],
                normalization=self.sim_bytscl[0],
            )
            return rvt.vis.byte_scale(data=norm_arr, no_data=np.nan, c_min=0, c_max=1)
        elif visualization == RVTVisualization.LOCAL_DOMINANCE:
            norm_arr = rvt.blend_func.normalize_image(
                visualization="ld",
                image=float_arr,
                min_norm=self.ld_bytscl[1],
                max_norm=self.ld_bytscl[2],
                normalization=self.ld_bytscl[0],
            )
            return rvt.vis.byte_scale(data=norm_arr, no_data=np.nan, c_min=0, c_max=1)
        elif visualization == RVTVisualization.MULTI_SCALE_RELIEF_MODEL:
            norm_arr = rvt.blend_func.normalize_image(
                visualization="msrm",
                image=float_arr,
                min_norm=self.msrm_bytscl[1],
                max_norm=self.msrm_bytscl[2],
                normalization=self.msrm_bytscl[0],
            )
            return rvt.vis.byte_scale(data=norm_arr, no_data=np.nan, c_min=0, c_max=1)
        elif visualization == RVTVisualization.MULTI_SCALE_TOPOGRAPHIC_POSITION:
            # This might not be necessary, as all mstp data should already be between 0 and 1
            norm_arr = rvt.blend_func.normalize_image(
                visualization="mstp",
                image=float_arr,
                min_norm=self.mstp_bytscl[1],
                max_norm=self.mstp_bytscl[2],
                normalization=self.mstp_bytscl[0],
            )

            return rvt.vis.byte_scale(data=norm_arr, no_data=np.nan, c_min=0, c_max=1)
        else:
            raise Exception(
                "rvt.default.DefaultValues.float_to_8bit: Wrong visualization (visualization) parameter!"
            )

    def get_slope(self, dem_arr, resolution_x, resolution_y, no_data=None):
        slope_arr = rvt.vis.slope_aspect(
            dem=dem_arr,
            resolution_x=resolution_x,
            resolution_y=resolution_y,
            ve_factor=self.ve_factor,
            output_units=self.slp_output_units,
            no_data=no_data,
        )["slope"]
        return slope_arr

    def save_slope(self, dem_path, custom_dir=None, save_float=None, save_8bit=None):
        """Calculates and saves Slope from dem (dem_path) with default parameters. If custom_dir is None it saves
        in dem directory else in custom_dir. If path to file already exists we can overwrite file (overwrite=0) or
        not (overwrite=1). If save_float is True method creates Gtiff with real values,
        if save_8bit is True method creates GTiff with bytescaled values (0-255)."""

        # if save_float is None it takes boolean from default (self)
        if save_float is None:
            save_float = self.slp_save_float
        # if save_8bit is None it takes boolean from default (self)
        if save_8bit is None:
            save_8bit = self.slp_save_8bit

        if not save_float and not save_8bit:
            raise Exception(
                "rvt.default.DefaultValues.save_slope: Both save_float and save_8bit are False,"
                " at least one of them has to be True!"
            )
        if not os.path.isfile(dem_path):
            raise Exception(
                "rvt.default.DefaultValues.save_slope: dem_path doesn't exist!"
            )

        if custom_dir is None:
            slope_path = self.get_slope_path(dem_path)
            slope_8bit_path = self.get_slope_path(dem_path, bit8=True)
        else:
            slope_path = os.path.join(custom_dir, self.get_slope_file_name(dem_path))
            slope_8bit_path = os.path.join(
                custom_dir, self.get_slope_file_name(dem_path, bit8=True)
            )

        # if file already exists and overwrite=0
        if save_float and save_8bit:
            if (
                os.path.isfile(slope_8bit_path)
                and os.path.isfile(slope_path)
                and not self.overwrite
            ):
                return 0
        elif save_float and not save_8bit:
            if os.path.isfile(slope_path) and not self.overwrite:
                return 0
        elif not save_float and not save_8bit:
            if os.path.isfile(slope_8bit_path) and not self.overwrite:
                return 0

        dem_size = get_raster_size(raster_path=dem_path)
        if dem_size[0] * dem_size[1] > self.tile_size_limit:  # tile by tile calculation
            if custom_dir is None:
                custom_dir = Path(dem_path).parent
            rvt.tile.save_rvt_visualization_tile_by_tile(
                rvt_visualization=RVTVisualization.SLOPE,
                rvt_default=self,
                dem_path=Path(dem_path),
                output_dir_path=Path(custom_dir),
                save_float=save_float,
                save_8bit=save_8bit,
            )
            return 1
        else:  # singleprocess
            dict_arr_res = get_raster_arr(raster_path=dem_path)
            dem_arr = dict_arr_res["array"]
            no_data = dict_arr_res["no_data"]
            x_res = dict_arr_res["resolution"][0]
            y_res = dict_arr_res["resolution"][1]
            slope_arr = self.get_slope(
                dem_arr=dem_arr, resolution_x=x_res, resolution_y=y_res, no_data=no_data
            )
            if save_float:
                if (
                    os.path.isfile(slope_path) and not self.overwrite
                ):  # file exists and overwrite=0
                    pass
                else:
                    save_raster(
                        src_raster_path=dem_path,
                        out_raster_path=slope_path,
                        out_raster_arr=slope_arr,
                        no_data=np.nan,
                    )
            if save_8bit:
                if (
                    os.path.isfile(slope_8bit_path) and not self.overwrite
                ):  # file exists and overwrite=0
                    pass
                else:
                    slope_8bit_arr = self.float_to_8bit(
                        float_arr=slope_arr, visualization=RVTVisualization.SLOPE
                    )
                    save_raster(
                        src_raster_path=dem_path,
                        out_raster_path=slope_8bit_path,
                        out_raster_arr=slope_8bit_arr,
                        e_type=1,
                    )
            return 1

    def get_shadow(self, dem_arr, resolution, no_data=None):
        shadow_arr = rvt.vis.shadow_horizon(
            dem=dem_arr,
            resolution=resolution,
            shadow_az=self.hs_sun_azi,
            shadow_el=self.hs_sun_el,
            ve_factor=self.ve_factor,
            no_data=no_data,
        )["shadow"]
        return shadow_arr

    def get_hillshade(self, dem_arr, resolution_x, resolution_y, no_data=None):
        hillshade_arr = rvt.vis.hillshade(
            dem=dem_arr,
            resolution_x=resolution_x,
            resolution_y=resolution_y,
            sun_azimuth=self.hs_sun_azi,
            sun_elevation=self.hs_sun_el,
            ve_factor=self.ve_factor,
            no_data=no_data,
        )
        return hillshade_arr

    def save_hillshade(
        self,
        dem_path,
        custom_dir=None,
        save_float=None,
        save_8bit=None,
        save_shadow=None,
    ):
        """Calculates and saves Hillshade from dem (dem_path) with default parameters. If custom_dir is None it saves
        in dem directory else in custom_dir. If path to file already exists we can overwrite file (overwrite=1)
        or not (overwrite=0). If save_float is True method creates Gtiff with real values,
        if save_8bit is True method creates GTiff with bytescaled values (0-255)."""

        # if save_float is None it takes boolean from default (self)
        if save_float is None:
            save_float = self.hs_save_float
        # if save_8bit is None it takes boolean from default (self)
        if save_8bit is None:
            save_8bit = self.hs_save_8bit
        # if save_shadow is None it takes boolean from default (self)
        if save_shadow is None:
            save_shadow = self.hs_shadow

        if not save_float and not save_8bit:
            raise Exception(
                "rvt.default.DefaultValues.save_hillshade: Both save_float and save_8bit are False,"
                " at least one of them has to be True!"
            )
        if not os.path.isfile(dem_path):
            raise Exception(
                "rvt.default.DefaultValues.save_hillshade: dem_path doesn't exist!"
            )

        if custom_dir is None:
            hillshade_path = self.get_hillshade_path(dem_path)
            hillshade_8bit_path = self.get_hillshade_path(dem_path, bit8=True)
            shadow_path = self.get_shadow_path(dem_path)
        else:
            hillshade_path = os.path.join(
                custom_dir, self.get_hillshade_file_name(dem_path)
            )
            hillshade_8bit_path = os.path.join(
                custom_dir, self.get_hillshade_file_name(dem_path, bit8=True)
            )
            shadow_path = os.path.join(custom_dir, self.get_shadow_path(dem_path))

        # if file already exists and overwrite=0
        if save_float and save_8bit and save_shadow:
            if (
                os.path.isfile(hillshade_8bit_path)
                and os.path.isfile(hillshade_path)
                and os.path.isfile(shadow_path)
                and not self.overwrite
            ):
                return 0
        elif save_float and not save_8bit and not save_shadow:
            if os.path.isfile(hillshade_path) and not self.overwrite:
                return 0
        elif not save_float and not save_8bit and not save_shadow:
            if os.path.isfile(hillshade_8bit_path) and not self.overwrite:
                return 0
        elif save_float and not save_8bit and save_shadow:
            if (
                os.path.isfile(hillshade_path)
                and os.path.isfile(shadow_path)
                and not self.overwrite
            ):
                return 0
        elif not save_float and not save_8bit and save_shadow:
            if (
                os.path.isfile(hillshade_8bit_path)
                and os.path.isfile(shadow_path)
                and not self.overwrite
            ):
                return 0

        dem_size = get_raster_size(raster_path=dem_path)
        if dem_size[0] * dem_size[1] > self.tile_size_limit:  # tile by tile
            if custom_dir is None:
                custom_dir = Path(dem_path).parent
            rvt.tile.save_rvt_visualization_tile_by_tile(
                rvt_visualization=RVTVisualization.HILLSHADE,
                rvt_default=self,
                dem_path=Path(dem_path),
                output_dir_path=Path(custom_dir),
                save_float=save_float,
                save_8bit=save_8bit,
            )
            if save_shadow:
                rvt.tile.save_rvt_visualization_tile_by_tile(
                    rvt_visualization=RVTVisualization.SHADOW,
                    rvt_default=self,
                    dem_path=Path(dem_path),
                    output_dir_path=Path(custom_dir),
                    save_float=True,
                    save_8bit=False,
                )
            return 1
        else:  # singleprocess
            dict_arr_res = get_raster_arr(raster_path=dem_path)
            dem_arr = dict_arr_res["array"]
            no_data = dict_arr_res["no_data"]
            x_res = dict_arr_res["resolution"][0]
            y_res = dict_arr_res["resolution"][1]
            hillshade_arr = self.get_hillshade(
                dem_arr=dem_arr, resolution_x=x_res, resolution_y=y_res, no_data=no_data
            ).astype("float32")
            if save_float:
                if (
                    os.path.isfile(hillshade_path) and not self.overwrite
                ):  # file exists and overwrite=0
                    pass
                else:
                    save_raster(
                        src_raster_path=dem_path,
                        out_raster_path=hillshade_path,
                        out_raster_arr=hillshade_arr,
                        no_data=np.nan,
                    )
            if save_8bit:
                if (
                    os.path.isfile(hillshade_8bit_path) and not self.overwrite
                ):  # file exists and overwrite=0
                    pass
                else:
                    hillshade_8_bit_arr = self.float_to_8bit(
                        float_arr=hillshade_arr,
                        visualization=RVTVisualization.HILLSHADE,
                    )
                    save_raster(
                        src_raster_path=dem_path,
                        out_raster_path=hillshade_8bit_path,
                        out_raster_arr=hillshade_8_bit_arr,
                        e_type=1,
                    )
            if save_shadow:
                shadow_arr = self.get_shadow(dem_arr=dem_arr, resolution=x_res)
                if (
                    os.path.isfile(shadow_path) and not self.overwrite
                ):  # file exists and overwrite=0
                    pass
                else:
                    save_raster(
                        src_raster_path=dem_path,
                        out_raster_path=shadow_path,
                        out_raster_arr=shadow_arr,
                        no_data=np.nan,
                    )
            return 1

    def get_multi_hillshade(self, dem_arr, resolution_x, resolution_y, no_data=None):
        multi_hillshade_arr = rvt.vis.multi_hillshade(
            dem=dem_arr,
            resolution_x=resolution_x,
            resolution_y=resolution_y,
            nr_directions=self.mhs_nr_dir,
            sun_elevation=self.mhs_sun_el,
            ve_factor=self.ve_factor,
            no_data=no_data,
        )
        return multi_hillshade_arr

    def save_multi_hillshade(
        self, dem_path, custom_dir=None, save_float=None, save_8bit=None
    ):
        """Calculates and saves Multidirectional hillshade from dem (dem_path) with default parameters.
        If custom_dir is None it saves in dem directory else in custom_dir. If path to file already exists we can
        overwrite file (overwrite=1) or not (overwrite=0). If save_float is True method creates Gtiff with real values,
        if save_8bit is True method creates GTiff with bytescaled values (0-255)."""

        # if save_float is None it takes boolean from default (self)
        if save_float is None:
            save_float = self.mhs_save_float
        # if save_8bit is None it takes boolean from default (self)
        if save_8bit is None:
            save_8bit = self.mhs_save_8bit

        if not save_float and not save_8bit:
            raise Exception(
                "rvt.default.DefaultValues.save_multi_hillshade: Both save_float and save_8bit are False,"
                " at least one of them has to be True!"
            )
        if not os.path.isfile(dem_path):
            raise Exception(
                "rvt.default.DefaultValues.save_multi_hillshade: dem_path doesn't exist!"
            )

        if custom_dir is None:
            multi_hillshade_path = self.get_multi_hillshade_path(dem_path)
            multi_hillshade_8bit_path = self.get_multi_hillshade_path(
                dem_path, bit8=True
            )
        else:
            multi_hillshade_path = os.path.join(
                custom_dir, self.get_multi_hillshade_file_name(dem_path)
            )
            multi_hillshade_8bit_path = os.path.join(
                custom_dir, self.get_multi_hillshade_file_name(dem_path, bit8=True)
            )

        # if file already exists and overwrite=0
        if save_float and save_8bit:
            if (
                os.path.isfile(multi_hillshade_8bit_path)
                and os.path.isfile(multi_hillshade_path)
                and not self.overwrite
            ):
                return 0
        elif save_float and not save_8bit:
            if os.path.isfile(multi_hillshade_path) and not self.overwrite:
                return 0
        elif not save_float and not save_8bit:
            if os.path.isfile(multi_hillshade_8bit_path) and not self.overwrite:
                return 0

        dem_size = get_raster_size(raster_path=dem_path)
        if dem_size[0] * dem_size[1] > self.tile_size_limit:  # tile by tile
            if custom_dir is None:
                custom_dir = Path(dem_path).parent
            rvt.tile.save_rvt_visualization_tile_by_tile(
                rvt_visualization=RVTVisualization.MULTI_HILLSHADE,
                rvt_default=self,
                dem_path=Path(dem_path),
                output_dir_path=Path(custom_dir),
                save_float=save_float,
                save_8bit=save_8bit,
            )
            return 1
        else:  # singleprocess
            dict_arr_res = get_raster_arr(raster_path=dem_path)
            dem_arr = dict_arr_res["array"]
            no_data = dict_arr_res["no_data"]
            x_res = dict_arr_res["resolution"][0]
            y_res = dict_arr_res["resolution"][1]
            if save_float:
                if (
                    os.path.isfile(multi_hillshade_path) and not self.overwrite
                ):  # file exists and overwrite=0
                    pass
                else:
                    multi_hillshade_arr = self.get_multi_hillshade(
                        dem_arr=dem_arr,
                        resolution_x=x_res,
                        resolution_y=y_res,
                        no_data=no_data,
                    ).astype("float32")
                    save_raster(
                        src_raster_path=dem_path,
                        out_raster_path=multi_hillshade_path,
                        out_raster_arr=multi_hillshade_arr,
                        no_data=np.nan,
                    )
            if save_8bit:
                if (
                    os.path.isfile(multi_hillshade_8bit_path) and not self.overwrite
                ):  # file exists and overwrite=0
                    pass
                else:
                    multi_hillshade_8bit_arr = self.float_to_8bit(
                        float_arr=dem_arr,
                        visualization=RVTVisualization.MULTI_HILLSHADE,
                        x_res=x_res,
                        y_res=y_res,
                        no_data=no_data,
                    )
                    save_raster(
                        src_raster_path=dem_path,
                        out_raster_path=multi_hillshade_8bit_path,
                        out_raster_arr=multi_hillshade_8bit_arr,
                        e_type=1,
                    )
            return 1

    def get_slrm(self, dem_arr, no_data=None):
        slrm_arr = rvt.vis.slrm(
            dem=dem_arr,
            radius_cell=self.slrm_rad_cell,
            ve_factor=self.ve_factor,
            no_data=no_data,
        )
        return slrm_arr

    def save_slrm(self, dem_path, custom_dir=None, save_float=None, save_8bit=None):
        """Calculates and saves Simple local relief model from dem (dem_path) with default parameters.
        If custom_dir is None it saves in dem directory else in custom_dir. If path to file already exists we can
        overwrite file (overwrite=1) or not (overwrite=0). If save_float is True method creates Gtiff with real values,
        if save_8bit is True method creates GTiff with bytescaled values (0-255)."""

        # if save_float is None it takes boolean from default (self)
        if save_float is None:
            save_float = self.slrm_save_float
        # if save_8bit is None it takes boolean from default (self)
        if save_8bit is None:
            save_8bit = self.slrm_save_8bit

        if not save_float and not save_8bit:
            raise Exception(
                "rvt.default.DefaultValues.save_slrm: Both save_float and save_8bit are False,"
                " at least one of them has to be True!"
            )
        if not os.path.isfile(dem_path):
            raise Exception(
                "rvt.default.DefaultValues.save_slrm: dem_path doesn't exist!"
            )

        if custom_dir is None:
            slrm_path = self.get_slrm_path(dem_path)
            slrm_8bit_path = self.get_slrm_path(dem_path, bit8=True)
        else:
            slrm_path = os.path.join(custom_dir, self.get_slrm_file_name(dem_path))
            slrm_8bit_path = os.path.join(
                custom_dir, self.get_slrm_file_name(dem_path, bit8=True)
            )

        # if file already exists and overwrite=0
        if save_float and save_8bit:
            if (
                os.path.isfile(slrm_8bit_path)
                and os.path.isfile(slrm_path)
                and not self.overwrite
            ):
                return 0
        elif save_float and not save_8bit:
            if os.path.isfile(slrm_path) and not self.overwrite:
                return 0
        elif not save_float and not save_8bit:
            if os.path.isfile(slrm_8bit_path) and not self.overwrite:
                return 0

        dem_size = get_raster_size(raster_path=dem_path)
        if dem_size[0] * dem_size[1] > self.tile_size_limit:  # tile by tile
            if custom_dir is None:
                custom_dir = Path(dem_path).parent
            rvt.tile.save_rvt_visualization_tile_by_tile(
                rvt_visualization=RVTVisualization.SIMPLE_LOCAL_RELIEF_MODEL,
                rvt_default=self,
                dem_path=Path(dem_path),
                output_dir_path=Path(custom_dir),
                save_float=save_float,
                save_8bit=save_8bit,
            )
            return 1
        else:  # singleprocess
            dict_arr_res = get_raster_arr(raster_path=dem_path)
            dem_arr = dict_arr_res["array"]
            no_data = dict_arr_res["no_data"]
            slrm_arr = self.get_slrm(dem_arr=dem_arr, no_data=no_data).astype("float32")
            if save_float:
                if (
                    os.path.isfile(slrm_path) and not self.overwrite
                ):  # file exists and overwrite=0
                    pass
                else:
                    save_raster(
                        src_raster_path=dem_path,
                        out_raster_path=slrm_path,
                        out_raster_arr=slrm_arr,
                        no_data=np.nan,
                    )
            if save_8bit:
                if (
                    os.path.isfile(slrm_8bit_path) and not self.overwrite
                ):  # file exists and overwrite=0
                    pass
                else:
                    slrm_8bit_arr = self.float_to_8bit(
                        float_arr=slrm_arr,
                        visualization=RVTVisualization.SIMPLE_LOCAL_RELIEF_MODEL,
                    )
                    save_raster(
                        src_raster_path=dem_path,
                        out_raster_path=slrm_8bit_path,
                        out_raster_arr=slrm_8bit_arr,
                        e_type=1,
                    )
            return 1

    def get_sky_view_factor(
        self,
        dem_arr,
        resolution,
        compute_svf=True,
        compute_asvf=False,
        compute_opns=False,
        no_data=None,
    ):
        dict_svf_asvf_opns = rvt.vis.sky_view_factor(
            dem=dem_arr,
            resolution=resolution,
            compute_svf=compute_svf,
            compute_opns=compute_opns,
            compute_asvf=compute_asvf,
            svf_n_dir=self.svf_n_dir,
            svf_r_max=self.svf_r_max,
            svf_noise=self.svf_noise,
            asvf_dir=self.asvf_dir,
            asvf_level=self.asvf_level,
            ve_factor=self.ve_factor,
            no_data=no_data,
        )
        return dict_svf_asvf_opns

    def save_sky_view_factor(
        self,
        dem_path,
        save_svf=True,
        save_asvf=False,
        save_opns=False,
        custom_dir=None,
        save_float=None,
        save_8bit=None,
    ):
        """Calculates and saves Sky-view factor(save_svf=True), Anisotropic Sky-view factor(save_asvf=True) and
        Positive Openness(save_opns=True) from dem (dem_path) with default parameters.
        If custom_dir is None it saves in dem directory else in custom_dir. If path to file already exists we can
        overwrite file (overwrite=1) or not (overwrite=0). If save_float is True method creates Gtiff with real values,
        if save_8bit is True method creates GTiff with bytescaled values (0-255)."""

        # if save_float is None it takes boolean from default (self)
        if save_float is None:
            save_float = self.svf_save_float
        # if save_8bit is None it takes boolean from default (self)
        if save_8bit is None:
            save_8bit = self.svf_save_8bit

        if not save_float and not save_8bit:
            raise Exception(
                "rvt.default.DefaultValues.save_sky_view_factor: Both save_float and save_8bit are False,"
                " at least one of them has to be True!"
            )
        if not os.path.isfile(dem_path):
            raise Exception(
                "rvt.default.DefaultValues.save_sky_view_factor: dem_path doesn't exist!"
            )

        svf_path = ""
        asvf_path = ""
        opns_path = ""
        svf_8bit_path = ""
        asvf_8bit_path = ""
        opns_8bit_path = ""
        if custom_dir is None:
            if save_svf:
                svf_path = self.get_svf_path(dem_path)
                svf_8bit_path = self.get_svf_path(dem_path, bit8=True)
            if save_asvf:
                asvf_path = self.get_asvf_path(dem_path)
                asvf_8bit_path = self.get_asvf_path(dem_path, bit8=True)
            if save_opns:
                opns_path = self.get_opns_path(dem_path)
                opns_8bit_path = self.get_opns_path(dem_path, bit8=True)
        else:
            if save_svf:
                svf_path = os.path.join(custom_dir, self.get_svf_file_name(dem_path))
                svf_8bit_path = os.path.join(
                    custom_dir, self.get_svf_file_name(dem_path, bit8=True)
                )
            if save_asvf:
                asvf_path = os.path.join(custom_dir, self.get_asvf_file_name(dem_path))
                asvf_8bit_path = os.path.join(
                    custom_dir, self.get_asvf_file_name(dem_path, bit8=True)
                )
            if save_opns:
                opns_path = os.path.join(custom_dir, self.get_opns_file_name(dem_path))
                opns_8bit_path = os.path.join(
                    custom_dir, self.get_opns_file_name(dem_path, bit8=True)
                )

        # if file already exists and overwrite=0
        if save_float and save_8bit:
            if (
                os.path.isfile(svf_path)
                and os.path.isfile(asvf_path)
                and os.path.isfile(opns_path)
                and os.path.isfile(svf_8bit_path)
                and os.path.isfile(asvf_8bit_path)
                and os.path.isfile(opns_8bit_path)
                and not self.overwrite
            ):
                return 0
        elif save_float and not save_8bit:
            if (
                os.path.isfile(svf_path)
                and os.path.isfile(asvf_path)
                and os.path.isfile(opns_path)
                and not self.overwrite
            ):
                return 0
        elif not save_float and save_8bit:
            if (
                os.path.isfile(svf_8bit_path)
                and os.path.isfile(asvf_8bit_path)
                and os.path.isfile(opns_8bit_path)
                and not self.overwrite
            ):
                return 0

        dem_size = get_raster_size(raster_path=dem_path)
        if dem_size[0] * dem_size[1] > self.tile_size_limit:  # tile by tile
            if custom_dir is None:
                custom_dir = Path(dem_path).parent
            if save_svf:
                rvt.tile.save_rvt_visualization_tile_by_tile(
                    rvt_visualization=RVTVisualization.SKY_VIEW_FACTOR,
                    rvt_default=self,
                    dem_path=Path(dem_path),
                    output_dir_path=Path(custom_dir),
                    save_float=save_float,
                    save_8bit=save_8bit,
                )
            if save_asvf:
                rvt.tile.save_rvt_visualization_tile_by_tile(
                    rvt_visualization=RVTVisualization.ANISOTROPIC_SKY_VIEW_FACTOR,
                    rvt_default=self,
                    dem_path=Path(dem_path),
                    output_dir_path=Path(custom_dir),
                    save_float=save_float,
                    save_8bit=save_8bit,
                )
            if save_opns:
                rvt.tile.save_rvt_visualization_tile_by_tile(
                    rvt_visualization=RVTVisualization.POSITIVE_OPENNESS,
                    rvt_default=self,
                    dem_path=Path(dem_path),
                    output_dir_path=Path(custom_dir),
                    save_float=save_float,
                    save_8bit=save_8bit,
                )
            return 1
        else:
            dict_arr_res = get_raster_arr(raster_path=dem_path)
            dem_arr = dict_arr_res["array"]
            no_data = dict_arr_res["no_data"]
            x_res = dict_arr_res["resolution"][0]
            y_res = dict_arr_res["resolution"][1]
            dict_svf_asvf_opns = self.get_sky_view_factor(
                dem_arr=dem_arr,
                resolution=x_res,
                compute_svf=save_svf,
                compute_asvf=save_asvf,
                compute_opns=save_opns,
                no_data=no_data,
            )
            if save_float:
                if save_svf:
                    if (
                        os.path.isfile(svf_path) and not self.overwrite
                    ):  # file exists and overwrite=0
                        pass
                    else:  # svf_path, file doesn't exists or exists and overwrite=1
                        save_raster(
                            src_raster_path=dem_path,
                            out_raster_path=svf_path,
                            out_raster_arr=dict_svf_asvf_opns["svf"].astype("float32"),
                            no_data=np.nan,
                        )
                if save_asvf:
                    if (
                        os.path.isfile(asvf_path) and not self.overwrite
                    ):  # file exists and overwrite=0
                        pass
                    else:  # asvf_path, file doesn't exists or exists and overwrite=1
                        save_raster(
                            src_raster_path=dem_path,
                            out_raster_path=asvf_path,
                            out_raster_arr=dict_svf_asvf_opns["asvf"].astype("float32"),
                            no_data=np.nan,
                        )
                if save_opns:
                    if (
                        os.path.isfile(opns_path) and not self.overwrite
                    ):  # file exists and overwrite=0
                        pass
                    else:  # opns_path, file doesn't exists or exists and overwrite=1
                        save_raster(
                            src_raster_path=dem_path,
                            out_raster_path=opns_path,
                            out_raster_arr=dict_svf_asvf_opns["opns"].astype("float32"),
                            no_data=np.nan,
                        )
            if save_8bit:
                if save_svf:
                    if (
                        os.path.isfile(svf_8bit_path) and not self.overwrite
                    ):  # file exists and overwrite=0
                        pass
                    else:  # svf_8bit_path, file doesn't exists or exists and overwrite=1
                        svf_8bit_arr = self.float_to_8bit(
                            float_arr=dict_svf_asvf_opns["svf"],
                            visualization=RVTVisualization.SKY_VIEW_FACTOR,
                        )
                        save_raster(
                            src_raster_path=dem_path,
                            out_raster_path=svf_8bit_path,
                            out_raster_arr=svf_8bit_arr,
                            e_type=1,
                        )
                if save_asvf:
                    if (
                        os.path.isfile(asvf_8bit_path) and not self.overwrite
                    ):  # file exists and overwrite=0
                        pass
                    else:  # asvf_8bit_path, file doesn't exists or exists and overwrite=1
                        asvf_8bit_arr = self.float_to_8bit(
                            float_arr=dict_svf_asvf_opns["asvf"],
                            visualization=RVTVisualization.ANISOTROPIC_SKY_VIEW_FACTOR,
                        )
                        save_raster(
                            src_raster_path=dem_path,
                            out_raster_path=asvf_8bit_path,
                            out_raster_arr=asvf_8bit_arr,
                            e_type=1,
                        )
                if save_opns:
                    if (
                        os.path.isfile(opns_8bit_path) and not self.overwrite
                    ):  # file exists and overwrite=0
                        pass
                    else:  # opns_8bit_path, file doesn't exists or exists and overwrite=1
                        opns_8bit_arr = self.float_to_8bit(
                            float_arr=dict_svf_asvf_opns["opns"],
                            visualization=RVTVisualization.POSITIVE_OPENNESS,
                        )
                        save_raster(
                            src_raster_path=dem_path,
                            out_raster_path=opns_8bit_path,
                            out_raster_arr=opns_8bit_arr,
                            e_type=1,
                        )
            return 1

    def get_neg_opns(self, dem_arr, resolution, no_data=None):
        dem_arr = -1 * dem_arr
        dict_neg_opns = rvt.vis.sky_view_factor(
            dem=dem_arr,
            resolution=resolution,
            svf_n_dir=self.svf_n_dir,
            svf_r_max=self.svf_r_max,
            svf_noise=self.svf_noise,
            compute_svf=False,
            compute_asvf=False,
            compute_opns=True,
            ve_factor=self.ve_factor,
            no_data=no_data,
        )
        neg_opns_arr = dict_neg_opns["opns"]
        return neg_opns_arr

    def save_neg_opns(self, dem_path, custom_dir=None, save_float=None, save_8bit=None):
        """Calculates and saves Negative Openness from dem (dem_path) with default parameters. If custom_dir is None
        it saves in dem directory else in custom_dir. If path to file already exists we can
        overwrite file (overwrite=1) or not (overwrite=0). If save_float is True method creates Gtiff with real values,
        if save_8bit is True method creates GTiff with bytescaled values (0-255)."""

        # if save_float is None it takes boolean from default (self)
        if save_float is None:
            save_float = self.neg_opns_save_float
        # if save_8bit is None it takes boolean from default (self)
        if save_8bit is None:
            save_8bit = self.neg_opns_save_8bit

        if not save_float and not save_8bit:
            raise Exception(
                "rvt.default.DefaultValues.save_neg_opns: Both save_float and save_8bit are False,"
                " at least one of them has to be True!"
            )
        if not os.path.isfile(dem_path):
            raise Exception(
                "rvt.default.DefaultValues.save_neg_opns: dem_path doesn't exist!"
            )

        if custom_dir is None:
            neg_opns_path = self.get_neg_opns_path(dem_path)
            neg_opns_8bit_path = self.get_neg_opns_path(dem_path, bit8=True)
        else:
            neg_opns_path = os.path.join(
                custom_dir, self.get_neg_opns_file_name(dem_path)
            )
            neg_opns_8bit_path = os.path.join(
                custom_dir, self.get_neg_opns_file_name(dem_path, bit8=True)
            )

        # if file already exists and overwrite=0
        if save_float and save_8bit:
            if (
                os.path.isfile(neg_opns_8bit_path)
                and os.path.isfile(neg_opns_path)
                and not self.overwrite
            ):
                return 0
        elif save_float and not save_8bit:
            if os.path.isfile(neg_opns_path) and not self.overwrite:
                return 0
        elif not save_float and not save_8bit:
            if os.path.isfile(neg_opns_8bit_path) and not self.overwrite:
                return 0

        dem_size = get_raster_size(raster_path=dem_path)
        if dem_size[0] * dem_size[1] > self.tile_size_limit:  # tile by tile
            if custom_dir is None:
                custom_dir = Path(dem_path).parent
            rvt.tile.save_rvt_visualization_tile_by_tile(
                rvt_visualization=RVTVisualization.NEGATIVE_OPENNESS,
                rvt_default=self,
                dem_path=Path(dem_path),
                output_dir_path=Path(custom_dir),
                save_float=save_float,
                save_8bit=save_8bit,
            )
            return 1
        else:  # singleprocess
            dict_arr_res = get_raster_arr(raster_path=dem_path)
            dem_arr = dict_arr_res["array"]
            no_data = dict_arr_res["no_data"]
            x_res = dict_arr_res["resolution"][0]
            y_res = dict_arr_res["resolution"][1]

            neg_opns_arr = self.get_neg_opns(
                dem_arr=dem_arr, resolution=x_res, no_data=no_data
            ).astype("float32")
            if save_float:
                if (
                    os.path.isfile(neg_opns_path) and not self.overwrite
                ):  # file exists and overwrite=0
                    pass
                else:
                    save_raster(
                        src_raster_path=dem_path,
                        out_raster_path=neg_opns_path,
                        out_raster_arr=neg_opns_arr,
                        no_data=np.nan,
                    )
            if save_8bit:
                if (
                    os.path.isfile(neg_opns_8bit_path) and not self.overwrite
                ):  # file exists and overwrite=0
                    pass
                else:
                    neg_opns_8bit_arr = self.float_to_8bit(
                        float_arr=neg_opns_arr,
                        visualization=RVTVisualization.NEGATIVE_OPENNESS,
                    )
                    save_raster(
                        src_raster_path=dem_path,
                        out_raster_path=neg_opns_8bit_path,
                        out_raster_arr=neg_opns_8bit_arr,
                        e_type=1,
                    )
            return 1

    def get_sky_illumination(self, dem_arr, resolution, no_data=None):
        sky_illumination_arr = rvt.vis.sky_illumination(
            dem=dem_arr,
            resolution=resolution,
            sky_model=self.sim_sky_mod,
            compute_shadow=bool(self.sim_compute_shadow),
            max_fine_radius=self.sim_shadow_dist,
            num_directions=self.sim_nr_dir,
            shadow_az=self.sim_shadow_az,
            shadow_el=self.sim_shadow_el,
            ve_factor=self.ve_factor,
            no_data=no_data,
        )
        return sky_illumination_arr

    def save_sky_illumination(
        self, dem_path, custom_dir=None, save_float=None, save_8bit=None
    ):
        """Calculates and saves Sky illumination from dem (dem_path) with default parameters. If custom_dir is None
        it saves in dem directory else in custom_dir. If path to file already exists we can
        overwrite file (overwrite=1) or not (overwrite=0). If save_float is True method creates Gtiff with real values,
        if save_8bit is True method creates GTiff with bytescaled values (0-255)."""

        # if save_float is None it takes boolean from default (self)
        if save_float is None:
            save_float = self.sim_save_float
        # if save_8bit is None it takes boolean from default (self)
        if save_8bit is None:
            save_8bit = self.sim_save_8bit

        if not save_float and not save_8bit:
            raise Exception(
                "rvt.default.DefaultValues.save_sky_illumination: Both save_float and save_8bit are False,"
                " at least one of them has to be True!"
            )
        if not os.path.isfile(dem_path):
            raise Exception(
                "rvt.default.DefaultValues.save_sky_illumination: dem_path doesn't exist!"
            )

        if custom_dir is None:
            sky_illumination_path = self.get_sky_illumination_path(dem_path)
            sky_illumination_8bit_path = self.get_sky_illumination_path(
                dem_path, bit8=True
            )
        else:
            sky_illumination_path = os.path.join(
                custom_dir, self.get_sky_illumination_file_name(dem_path)
            )
            sky_illumination_8bit_path = os.path.join(
                custom_dir, self.get_sky_illumination_file_name(dem_path, bit8=True)
            )

        # if file already exists and overwrite=0
        if save_float and save_8bit:
            if (
                os.path.isfile(sky_illumination_8bit_path)
                and os.path.isfile(sky_illumination_path)
                and not self.overwrite
            ):
                return 0
        elif save_float and not save_8bit:
            if os.path.isfile(sky_illumination_path) and not self.overwrite:
                return 0
        elif not save_float and not save_8bit:
            if os.path.isfile(sky_illumination_8bit_path) and not self.overwrite:
                return 0

        dem_size = get_raster_size(raster_path=dem_path)
        if dem_size[0] * dem_size[1] > self.tile_size_limit:  # tile by tile
            if custom_dir is None:
                custom_dir = Path(dem_path).parent
            rvt.tile.save_rvt_visualization_tile_by_tile(
                rvt_visualization=RVTVisualization.SKY_ILLUMINATION,
                rvt_default=self,
                dem_path=Path(dem_path),
                output_dir_path=Path(custom_dir),
                save_float=save_float,
                save_8bit=save_8bit,
            )
            return 1
        else:  # singleprocess
            dict_arr_res = get_raster_arr(raster_path=dem_path)
            dem_arr = dict_arr_res["array"]
            no_data = dict_arr_res["no_data"]
            x_res = dict_arr_res["resolution"][0]
            y_res = dict_arr_res["resolution"][1]

            sky_illumination_arr = self.get_sky_illumination(
                dem_arr=dem_arr, resolution=x_res, no_data=no_data
            ).astype("float32")
            if save_float:
                if (
                    os.path.isfile(sky_illumination_path) and not self.overwrite
                ):  # file exists and overwrite=0
                    pass
                else:
                    save_raster(
                        src_raster_path=dem_path,
                        out_raster_path=sky_illumination_path,
                        out_raster_arr=sky_illumination_arr,
                        no_data=np.nan,
                    )
            if save_8bit:
                if (
                    os.path.isfile(sky_illumination_8bit_path) and not self.overwrite
                ):  # file exists and overwrite=0
                    pass
                else:
                    sky_illumination_8bit_arr = self.float_to_8bit(
                        float_arr=sky_illumination_arr,
                        visualization=RVTVisualization.SKY_ILLUMINATION,
                    )
                    save_raster(
                        src_raster_path=dem_path,
                        out_raster_path=sky_illumination_8bit_path,
                        out_raster_arr=sky_illumination_8bit_arr,
                        e_type=1,
                    )
            return 1

    def get_local_dominance(self, dem_arr, no_data=None):
        local_dominance_arr = rvt.vis.local_dominance(
            dem=dem_arr,
            min_rad=self.ld_min_rad,
            max_rad=self.ld_max_rad,
            rad_inc=self.ld_rad_inc,
            angular_res=self.ld_anglr_res,
            observer_height=self.ld_observer_h,
            ve_factor=self.ve_factor,
            no_data=no_data,
        )
        return local_dominance_arr

    def save_local_dominance(
        self, dem_path, custom_dir=None, save_float=None, save_8bit=None
    ):
        """Calculates and saves Local dominance from dem (dem_path) with default parameters. If custom_dir is None
        it saves in dem directory else in custom_dir. If path to file already exists we can
        overwrite file (overwrite=1) or not (overwrite=0). If save_float is True method creates Gtiff with real values,
        if save_8bit is True method creates GTiff with bytescaled values (0-255)."""

        # if save_float is None it takes boolean from default (self)
        if save_float is None:
            save_float = self.ld_save_float
        # if save_8bit is None it takes boolean from default (self)
        if save_8bit is None:
            save_8bit = self.ld_save_8bit

        if not save_float and not save_8bit:
            raise Exception(
                "rvt.default.DefaultValues.save_local_dominance: Both save_float and save_8bit are False,"
                " at least one of them has to be True!"
            )
        if not os.path.isfile(dem_path):
            raise Exception(
                "rvt.default.DefaultValues.save_local_dominance: dem_path doesn't exist!"
            )

        if custom_dir is None:
            local_dominance_path = self.get_local_dominance_path(dem_path)
            local_dominance_8bit_path = self.get_local_dominance_path(
                dem_path, bit8=True
            )
        else:
            local_dominance_path = os.path.join(
                custom_dir, self.get_local_dominance_file_name(dem_path)
            )
            local_dominance_8bit_path = os.path.join(
                custom_dir, self.get_local_dominance_file_name(dem_path, bit8=True)
            )

        # if file already exists and overwrite=0
        if save_float and save_8bit:
            if (
                os.path.isfile(local_dominance_8bit_path)
                and os.path.isfile(local_dominance_path)
                and not self.overwrite
            ):
                return 0
        elif save_float and not save_8bit:
            if os.path.isfile(local_dominance_path) and not self.overwrite:
                return 0
        elif not save_float and not save_8bit:
            if os.path.isfile(local_dominance_8bit_path) and not self.overwrite:
                return 0

        dem_size = get_raster_size(raster_path=dem_path)
        if dem_size[0] * dem_size[1] > self.tile_size_limit:  # tile by tile
            if custom_dir is None:
                custom_dir = Path(dem_path).parent
            rvt.tile.save_rvt_visualization_tile_by_tile(
                rvt_visualization=RVTVisualization.LOCAL_DOMINANCE,
                rvt_default=self,
                dem_path=Path(dem_path),
                output_dir_path=Path(custom_dir),
                save_float=save_float,
                save_8bit=save_8bit,
            )
            return 1
        else:  # singleprocess
            dict_arr_res = get_raster_arr(raster_path=dem_path)
            dem_arr = dict_arr_res["array"]
            no_data = dict_arr_res["no_data"]
            local_dominance_arr = self.get_local_dominance(
                dem_arr=dem_arr, no_data=no_data
            ).astype("float32")
            if save_float:
                if (
                    os.path.isfile(local_dominance_path) and not self.overwrite
                ):  # file exists and overwrite=0
                    pass
                else:
                    save_raster(
                        src_raster_path=dem_path,
                        out_raster_path=local_dominance_path,
                        out_raster_arr=local_dominance_arr,
                        no_data=np.nan,
                    )
            if save_8bit:
                if (
                    os.path.isfile(local_dominance_8bit_path) and not self.overwrite
                ):  # file exists and overwrite=0
                    pass
                else:
                    local_dominance_8bit_arr = self.float_to_8bit(
                        float_arr=local_dominance_arr,
                        visualization=RVTVisualization.LOCAL_DOMINANCE,
                    )
                    save_raster(
                        src_raster_path=dem_path,
                        out_raster_path=local_dominance_8bit_path,
                        out_raster_arr=local_dominance_8bit_arr,
                        e_type=1,
                    )
            return 1

    def get_msrm(self, dem_arr, resolution, no_data=None):
        msrm_arr = rvt.vis.msrm(
            dem=dem_arr,
            resolution=resolution,
            feature_min=self.msrm_feature_min,
            feature_max=self.msrm_feature_max,
            scaling_factor=self.msrm_scaling_factor,
            ve_factor=self.ve_factor,
            no_data=no_data,
        )
        return msrm_arr

    def save_msrm(self, dem_path, custom_dir=None, save_float=None, save_8bit=None):
        """Calculates and saves Multi-scale relief model from dem (dem_path) with default parameters.
        If custom_dir is None it saves in dem directory else in custom_dir. If path to file already exists we can
        overwrite file (overwrite=1) or not (overwrite=0). If save_float is True method creates Gtiff with real values,
        if save_8bit is True method creates GTiff with bytescaled values (0-255)."""

        # if save_float is None it takes boolean from default (self)
        if save_float is None:
            save_float = self.msrm_save_float
        # if save_8bit is None it takes boolean from default (self)
        if save_8bit is None:
            save_8bit = self.msrm_save_8bit

        if not save_float and not save_8bit:
            raise Exception(
                "rvt.default.DefaultValues.save_msrm: Both save_float and save_8bit are False,"
                " at least one of them has to be True!"
            )
        if not os.path.isfile(dem_path):
            raise Exception(
                "rvt.default.DefaultValues.save_msrm: dem_path doesn't exist!"
            )

        if custom_dir is None:
            msrm_path = self.get_msrm_path(dem_path)
            msrm_8bit_path = self.get_msrm_path(dem_path, bit8=True)
        else:
            msrm_path = os.path.join(custom_dir, self.get_msrm_file_name(dem_path))
            msrm_8bit_path = os.path.join(
                custom_dir, self.get_msrm_file_name(dem_path, bit8=True)
            )

        # if file already exists and overwrite=0
        if save_float and save_8bit:
            if (
                os.path.isfile(msrm_8bit_path)
                and os.path.isfile(msrm_path)
                and not self.overwrite
            ):
                return 0
        elif save_float and not save_8bit:
            if os.path.isfile(msrm_path) and not self.overwrite:
                return 0
        elif not save_float and not save_8bit:
            if os.path.isfile(msrm_8bit_path) and not self.overwrite:
                return 0

        dem_size = get_raster_size(raster_path=dem_path)
        if dem_size[0] * dem_size[1] > self.tile_size_limit:  # tile by tile
            if custom_dir is None:
                custom_dir = Path(dem_path).parent
            rvt.tile.save_rvt_visualization_tile_by_tile(
                rvt_visualization=RVTVisualization.MULTI_SCALE_RELIEF_MODEL,
                rvt_default=self,
                dem_path=Path(dem_path),
                output_dir_path=Path(custom_dir),
                save_float=save_float,
                save_8bit=save_8bit,
            )
            return 1
        else:  # singleprocess
            dict_arr_res = get_raster_arr(raster_path=dem_path)
            dem_arr = dict_arr_res["array"]
            no_data = dict_arr_res["no_data"]
            x_res = dict_arr_res["resolution"][0]
            y_res = dict_arr_res["resolution"][1]

            msrm_arr = self.get_msrm(
                dem_arr=dem_arr, resolution=x_res, no_data=no_data
            ).astype("float32")
            if save_float:
                if (
                    os.path.isfile(msrm_path) and not self.overwrite
                ):  # file exists and overwrite=0
                    pass
                else:
                    save_raster(
                        src_raster_path=dem_path,
                        out_raster_path=msrm_path,
                        out_raster_arr=msrm_arr,
                        no_data=np.nan,
                    )
            if save_8bit:
                if (
                    os.path.isfile(msrm_8bit_path) and not self.overwrite
                ):  # file exists and overwrite=0
                    pass
                else:
                    msrm_8bit_arr = self.float_to_8bit(
                        float_arr=msrm_arr,
                        visualization=RVTVisualization.MULTI_SCALE_RELIEF_MODEL,
                    )
                    save_raster(
                        src_raster_path=dem_path,
                        out_raster_path=msrm_8bit_path,
                        out_raster_arr=msrm_8bit_arr,
                        e_type=1,
                    )
            return 1

    def get_mstp(self, dem_arr, no_data=None):
        mstp_arr = rvt.vis.mstp(
            dem=dem_arr,
            local_scale=self.mstp_local_scale,
            meso_scale=self.mstp_meso_scale,
            broad_scale=self.mstp_broad_scale,
            lightness=self.mstp_lightness,
            no_data=no_data,
        )
        return mstp_arr

    def save_mstp(self, dem_path, custom_dir=None, save_float=None, save_8bit=None):
        """Calculates and saves Multi-scale topographic position from dem (dem_path) with default parameters.
        If custom_dir is None it saves in dem directory else in custom_dir. If path to file already exists we can
        overwrite file (overwrite=1) or not (overwrite=0)."""

        # if save_float is None it takes boolean from default (self)
        if save_float is None:
            save_float = self.mstp_save_float
        # if save_8bit is None it takes boolean from default (self)
        if save_8bit is None:
            save_8bit = self.mstp_save_8bit

        if not save_float and not save_8bit:
            raise Exception(
                "rvt.default.DefaultValues.save_mstp: Both save_float and save_8bit are False,"
                " at least one of them has to be True!"
            )

        if not os.path.isfile(dem_path):
            raise Exception(
                "rvt.default.DefaultValues.save_mstp: dem_path doesn't exist!"
            )

        if custom_dir is None:
            mstp_path = self.get_mstp_path(dem_path)
            mstp_8bit_path = self.get_mstp_path(dem_path, bit8=True)
        else:
            mstp_path = os.path.join(custom_dir, self.get_mstp_file_name(dem_path))
            mstp_8bit_path = os.path.join(
                custom_dir, self.get_mstp_file_name(dem_path, bit8=True)
            )

        # if file already exists and overwrite=0
        if save_float and save_8bit:
            if (
                os.path.isfile(mstp_8bit_path)
                and os.path.isfile(mstp_path)
                and not self.overwrite
            ):
                return 0
        elif save_float and not save_8bit:
            if os.path.isfile(mstp_path) and not self.overwrite:
                return 0
        elif not save_float and not save_8bit:
            if os.path.isfile(mstp_8bit_path) and not self.overwrite:
                return 0

        dem_size = get_raster_size(raster_path=dem_path)
        if dem_size[0] * dem_size[1] > self.tile_size_limit:  # tile by tile calculation
            if custom_dir is None:
                custom_dir = Path(dem_path).parent
            rvt.tile.save_rvt_visualization_tile_by_tile(
                rvt_visualization=RVTVisualization.MULTI_SCALE_TOPOGRAPHIC_POSITION,
                rvt_default=self,
                dem_path=Path(dem_path),
                output_dir_path=Path(custom_dir),
                save_float=save_float,
                save_8bit=save_8bit,
            )
            return 1
        else:  # singleprocess
            dict_arr_res = get_raster_arr(raster_path=dem_path)
            dem_arr = dict_arr_res["array"]
            no_data = dict_arr_res["no_data"]

            mstp_arr = self.get_mstp(dem_arr=dem_arr, no_data=no_data)

            if save_float:
                if (
                    os.path.isfile(mstp_path) and not self.overwrite
                ):  # file exists and overwrite=0
                    pass
                else:
                    save_raster(
                        src_raster_path=dem_path,
                        out_raster_path=mstp_path,
                        out_raster_arr=mstp_arr,
                        no_data=np.nan,
                        e_type=6,
                    )
            if save_8bit:
                if (
                    os.path.isfile(mstp_8bit_path) and not self.overwrite
                ):  # file exists and overwrite=0
                    pass
                else:
                    mstp_8bit_arr = self.float_to_8bit(
                        float_arr=mstp_arr,
                        visualization=RVTVisualization.MULTI_SCALE_TOPOGRAPHIC_POSITION,
                    )
                    save_raster(
                        src_raster_path=dem_path,
                        out_raster_path=mstp_8bit_path,
                        out_raster_arr=mstp_8bit_arr,
                        no_data=np.nan,
                        e_type=1,
                    )

            return 1

    def save_visualizations(self, dem_path, custom_dir=None):
        """Save all visualizations where self.'visualization'_compute = True also saves float where self.'visualization'
        _save_float = True and 8bit where self.'visualization'_save_8bit = True. In the end method creates log file."""
        start_time = time.time()
        if self.slp_compute:
            self.save_slope(dem_path, custom_dir=custom_dir)
        if self.hs_compute:
            self.save_hillshade(dem_path, custom_dir=custom_dir)
        if self.mhs_compute:
            self.save_multi_hillshade(dem_path, custom_dir=custom_dir)
        if self.slrm_compute:
            self.save_slrm(dem_path, custom_dir=custom_dir)
        if self.svf_compute or self.asvf_compute or self.pos_opns_compute:
            self.save_sky_view_factor(
                dem_path,
                save_svf=bool(self.svf_compute),
                save_asvf=bool(self.asvf_compute),
                save_opns=bool(self.pos_opns_compute),
                custom_dir=custom_dir,
            )
        if self.neg_opns_compute:
            self.save_neg_opns(dem_path, custom_dir=custom_dir)
        if self.sim_compute:
            self.save_sky_illumination(dem_path, custom_dir=custom_dir)
        if self.ld_compute:
            self.save_local_dominance(dem_path, custom_dir=custom_dir)
        if self.msrm_compute:
            self.save_msrm(dem_path, custom_dir=custom_dir)
        if self.mstp_compute:
            self.save_mstp(dem_path, custom_dir=custom_dir)
        end_time = time.time()
        compute_time = end_time - start_time
        self.create_log_file(
            dem_path=dem_path, custom_dir=custom_dir, compute_time=compute_time
        )

    def calculate_visualization(
        self,
        visualization: RVTVisualization,
        dem: np.array,
        resolution_x: float,
        resolution_y: float,
        no_data: Optional[float] = None,
        save_float: bool = True,
        save_8bit: bool = False,
    ) -> Optional[Tuple[np.array, np.array]]:  # tuple[vis_float_arr, vis_8bit_arr]
        vis_arr = None
        vis_float_arr = None
        vis_8bit_arr = None
        if visualization == RVTVisualization.SLOPE:
            vis_arr = self.get_slope(
                dem_arr=dem,
                resolution_x=resolution_x,
                resolution_y=resolution_y,
                no_data=no_data,
            )
        elif visualization == RVTVisualization.SHADOW:
            vis_arr = self.get_shadow(
                dem_arr=dem, resolution=resolution_x, no_data=no_data
            )
        elif visualization == RVTVisualization.HILLSHADE:
            vis_arr = self.get_hillshade(
                dem_arr=dem,
                resolution_x=resolution_x,
                resolution_y=resolution_y,
                no_data=no_data,
            )
        elif visualization == RVTVisualization.MULTI_HILLSHADE:
            vis_arr = self.get_multi_hillshade(
                dem_arr=dem,
                resolution_x=resolution_x,
                resolution_y=resolution_y,
                no_data=no_data,
            )
        elif visualization == RVTVisualization.SIMPLE_LOCAL_RELIEF_MODEL:
            vis_arr = self.get_slrm(dem_arr=dem, no_data=no_data)
        elif visualization == RVTVisualization.SKY_VIEW_FACTOR:
            vis_arr = self.get_sky_view_factor(
                dem_arr=dem,
                resolution=resolution_x,
                compute_svf=True,
                compute_asvf=False,
                compute_opns=False,
                no_data=no_data,
            )["svf"]
        elif visualization == RVTVisualization.ANISOTROPIC_SKY_VIEW_FACTOR:
            vis_arr = self.get_sky_view_factor(
                dem_arr=dem,
                resolution=resolution_x,
                compute_svf=False,
                compute_asvf=True,
                compute_opns=False,
                no_data=no_data,
            )["asvf"]
        elif visualization == RVTVisualization.POSITIVE_OPENNESS:
            vis_arr = self.get_sky_view_factor(
                dem_arr=dem,
                resolution=resolution_x,
                compute_svf=False,
                compute_asvf=False,
                compute_opns=True,
                no_data=no_data,
            )["opns"]
        elif visualization == RVTVisualization.NEGATIVE_OPENNESS:
            vis_arr = self.get_neg_opns(
                dem_arr=dem, resolution=resolution_x, no_data=no_data
            )
        elif visualization == RVTVisualization.SKY_ILLUMINATION:
            vis_arr = self.get_sky_illumination(
                dem_arr=dem, resolution=resolution_x, no_data=no_data
            )
        elif visualization == RVTVisualization.LOCAL_DOMINANCE:
            vis_arr = self.get_local_dominance(dem_arr=dem, no_data=no_data)
        elif visualization == RVTVisualization.MULTI_SCALE_RELIEF_MODEL:
            vis_arr = self.get_msrm(
                dem_arr=dem, resolution=resolution_x, no_data=no_data
            )
        elif visualization == RVTVisualization.MULTI_SCALE_TOPOGRAPHIC_POSITION:
            vis_arr = self.get_mstp(dem_arr=dem, no_data=no_data)
        if save_float:
            vis_float_arr = vis_arr
        if save_8bit:
            if visualization == RVTVisualization.MULTI_HILLSHADE:
                vis_8bit_arr = self.float_to_8bit(
                    float_arr=dem,
                    visualization=visualization,
                    x_res=resolution_x,
                    y_res=resolution_y,
                    no_data=no_data,
                )
            else:
                vis_8bit_arr = self.float_to_8bit(
                    float_arr=vis_arr, visualization=visualization
                )
        return vis_float_arr, vis_8bit_arr

    def create_log_file(self, dem_path, custom_dir=None, compute_time=None):
        """Creates log file in custom_dir, if custom_dir=None it creates it in dem directory (dem_path).
        Be aware, all default parameters have to be right! Parameter compute_time is in seconds."""
        dict_arr_res = get_raster_arr(raster_path=dem_path)
        resolution = dict_arr_res["resolution"]
        arr_shape = np.array(dict_arr_res["array"]).shape
        del dict_arr_res
        nr_bands = 0
        nr_cols = 0
        nr_rows = 0
        if len(arr_shape) == 3:
            nr_bands = arr_shape[0]
            nr_rows = arr_shape[1]
            nr_cols = arr_shape[2]
        elif len(arr_shape) == 2:
            nr_bands = 1
            nr_rows = arr_shape[0]
            nr_cols = arr_shape[1]
        dem_dir = os.path.dirname(dem_path)
        log_dir = dem_dir
        if custom_dir is not None:
            log_dir = custom_dir
        dem_name = os.path.splitext(os.path.basename(dem_path))[0]
        log_file_time = datetime.datetime.now()
        log_file_time_str = log_file_time.strftime("%Y-%m-%d_%H-%M-%S")
        log_name = "{}_vis_log_{}.txt".format(dem_name, log_file_time_str)
        log_path = os.path.join(log_dir, log_name)
        dat = open(log_path, "w")
        dat.write(
            "===============================================================================================\n"
            "Relief Visualization Toolbox (python), visualizations log\n"
            "Copyright:\n"
            "\tResearch Centre of the Slovenian Academy of Sciences and Arts\n"
            "\tUniversity of Ljubljana, Faculty of Civil and Geodetic Engineering\n"
            "===============================================================================================\n"
        )
        dat.write("\n\n\n")

        dat.write(
            "Processing info about visualizations\n"
            "===============================================================================================\n\n"
        )
        dat.write("# Metadata of the input file\n\n")
        dat.write("\tInput filename:\t\t{}\n".format(dem_path))
        dat.write("\tNumber of rows:\t\t{}\n".format(nr_rows))
        dat.write("\tNumber of columns:\t{}\n".format(nr_cols))
        dat.write("\tNumber of bands:\t{}\n".format(nr_bands))
        dat.write("\tResolution (x, y):\t{}, {}\n".format(resolution[0], resolution[1]))
        dat.write("\n")

        dat.write("# Selected visualization parameters\n")
        dat.write("\tOverwrite: {}\n".format(self.overwrite))
        dat.write("\tVertical exaggeration factor: {}\n".format(self.ve_factor))
        if nr_rows * nr_cols > self.tile_size_limit:
            dat.write("\tCalculating tile by tile: {}\n".format("ON"))
            dat.write(
                "\t\tTile block size: {}x{}\n".format(
                    self.tile_size[0], self.tile_size[1]
                )
            )
        else:
            dat.write("\tCalculating tile by tile: {}\n".format("OFF"))

        dat.write("\n")

        dat.write("# The following visualizations have been preformed:\n\n")
        if self.hs_compute:
            dat.write("\tHillshade\n")
            dat.write("\t\ths_sun_el=\t\t{}\n".format(self.hs_sun_el))
            dat.write("\t\ths_sun_azi=\t\t{}\n".format(self.hs_sun_azi))
            if self.hs_save_float:
                dat.write("\t\t>> Output file:\n")
                dat.write(
                    "\t\t\t{}\n".format(
                        os.path.abspath(
                            os.path.join(
                                log_dir, self.get_hillshade_file_name(dem_path)
                            )
                        )
                    )
                )
            if self.hs_save_8bit:
                dat.write("\t\t>> Output 8bit file:\n")
                dat.write(
                    "\t\ths_bytscl=\t\t({}, {}, {})\n".format(
                        self.hs_bytscl[0], self.hs_bytscl[1], self.hs_bytscl[2]
                    )
                )
                dat.write(
                    "\t\t\t{}\n".format(
                        os.path.abspath(
                            os.path.join(
                                log_dir,
                                self.get_hillshade_file_name(dem_path, bit8=True),
                            )
                        )
                    )
                )
            if self.hs_shadow:
                dat.write("\t\t>> Output shadow file:\n")
                dat.write(
                    "\t\t\t{}\n".format(
                        os.path.abspath(
                            os.path.join(log_dir, self.get_shadow_file_name(dem_path))
                        )
                    )
                )
            dat.write("\n")
        if self.mhs_compute:
            dat.write("\tMultiple directions hillshade\n")
            dat.write("\t\tmhs_sun_el=\t\t{}\n".format(self.mhs_sun_el))
            dat.write("\t\tmhs_nr_dir=\t\t{}\n".format(self.mhs_nr_dir))
            if self.mhs_save_float:
                dat.write("\t\t>> Output file:\n")
                dat.write(
                    "\t\t\t{}\n".format(
                        os.path.abspath(
                            os.path.join(
                                log_dir, self.get_multi_hillshade_file_name(dem_path)
                            )
                        )
                    )
                )
            if self.mhs_save_8bit:
                dat.write("\t\t>> Output 8bit file:\n")
                dat.write(
                    "\t\tmhs_bytscl=\t\t({}, {}, {})\n".format(
                        self.mhs_bytscl[0], self.mhs_bytscl[1], self.mhs_bytscl[2]
                    )
                )
                dat.write(
                    "\t\t\t{}\n".format(
                        os.path.abspath(
                            os.path.join(
                                log_dir,
                                self.get_multi_hillshade_file_name(dem_path, bit8=True),
                            )
                        )
                    )
                )
            dat.write("\n")
        if self.slp_compute:
            dat.write("\tSlope gradient\n")
            dat.write("\t\tslp_output_units=\t\t{}\n".format(self.slp_output_units))
            if self.slp_save_float:
                dat.write("\t\t>> Output file:\n")
                dat.write(
                    "\t\t\t{}\n".format(
                        os.path.abspath(
                            os.path.join(log_dir, self.get_slope_file_name(dem_path))
                        )
                    )
                )
            if self.slp_save_8bit:
                dat.write("\t\t>> Output 8bit file:\n")
                dat.write(
                    "\t\tslp_bytscl=\t\t({}, {}, {})\n".format(
                        self.slp_bytscl[0], self.slp_bytscl[1], self.slp_bytscl[2]
                    )
                )
                dat.write(
                    "\t\t\t{}\n".format(
                        os.path.abspath(
                            os.path.join(
                                log_dir, self.get_slope_file_name(dem_path, bit8=True)
                            )
                        )
                    )
                )
            dat.write("\n")
        if self.slrm_compute:
            dat.write("\tSimple local relief model\n")
            dat.write("\t\tslrm_rad_cell=\t\t{}\n".format(self.slrm_rad_cell))
            if self.slrm_save_float:
                dat.write("\t\t>> Output file:\n")
                dat.write(
                    "\t\t\t{}\n".format(
                        os.path.abspath(
                            os.path.join(log_dir, self.get_slrm_file_name(dem_path))
                        )
                    )
                )
            if self.slrm_save_8bit:
                dat.write("\t\t>> Output 8bit file:\n")
                dat.write(
                    "\t\tslrm_bytscl=\t\t({}, {}, {})\n".format(
                        self.slrm_bytscl[0], self.slrm_bytscl[1], self.slrm_bytscl[2]
                    )
                )
                dat.write(
                    "\t\t\t{}\n".format(
                        os.path.abspath(
                            os.path.join(
                                log_dir, self.get_slrm_file_name(dem_path, bit8=True)
                            )
                        )
                    )
                )
            dat.write("\n")
        if self.msrm_compute:
            dat.write("\tMulti-scale relief model\n")
            dat.write("\t\tmsrm_feature_min=\t\t{}\n".format(self.msrm_feature_min))
            dat.write("\t\tmsrm_feature_max=\t\t{}\n".format(self.msrm_feature_max))
            dat.write(
                "\t\tmsrm_scaling_factor=\t\t{}\n".format(self.msrm_scaling_factor)
            )
            if self.msrm_save_float:
                dat.write("\t\t>> Output file:\n")
                dat.write("\t\t>> Output file:\n")
                dat.write(
                    "\t\t\t{}\n".format(
                        os.path.abspath(
                            os.path.join(log_dir, self.get_msrm_file_name(dem_path))
                        )
                    )
                )
            if self.msrm_save_8bit:
                dat.write("\t\t>> Output 8bit file:\n")
                dat.write(
                    "\t\tmsrm_bytscl=\t\t({}, {}, {})\n".format(
                        self.msrm_bytscl[0], self.msrm_bytscl[1], self.msrm_bytscl[2]
                    )
                )
                dat.write(
                    "\t\t\t{}\n".format(
                        os.path.abspath(
                            os.path.join(
                                log_dir, self.get_msrm_file_name(dem_path, bit8=True)
                            )
                        )
                    )
                )
            dat.write("\n")
        if self.svf_compute:
            dat.write("\tSky-View Factor\n")
            dat.write("\t\tsvf_n_dir=\t\t{}\n".format(self.svf_n_dir))
            dat.write("\t\tsvf_noise=\t\t{}\n".format(self.svf_noise))
            dat.write("\t\tsvf_r_max=\t\t{}\n".format(self.svf_r_max))
            if self.svf_save_float:
                dat.write("\t\t>> Output file:\n")
                dat.write(
                    "\t\t\t{}\n".format(
                        os.path.abspath(
                            os.path.join(log_dir, self.get_svf_file_name(dem_path))
                        )
                    )
                )
            if self.svf_save_8bit:
                dat.write("\t\t>> Output 8bit file:\n")
                dat.write(
                    "\t\tsvf_bytscl=\t\t({}, {}, {})\n".format(
                        self.svf_bytscl[0], self.svf_bytscl[1], self.svf_bytscl[2]
                    )
                )
                dat.write(
                    "\t\t\t{}\n".format(
                        os.path.abspath(
                            os.path.join(
                                log_dir, self.get_svf_file_name(dem_path, bit8=True)
                            )
                        )
                    )
                )
            dat.write("\n")
        if self.asvf_compute:
            dat.write("\tAnisotropic Sky-View Factor\n")
            dat.write("\t\tsvf_n_dir=\t\t{}\n".format(self.svf_n_dir))
            dat.write("\t\tsvf_noise=\t\t{}\n".format(self.svf_noise))
            dat.write("\t\tsvf_r_max=\t\t{}\n".format(self.svf_r_max))
            dat.write("\t\tasvf_level=\t\t{}\n".format(self.asvf_level))
            dat.write("\t\tasvf_dir=\t\t{}\n".format(self.asvf_dir))
            if self.svf_save_float:
                dat.write("\t\t>> Output file:\n")
                dat.write(
                    "\t\t\t{}\n".format(
                        os.path.abspath(
                            os.path.join(log_dir, self.get_asvf_file_name(dem_path))
                        )
                    )
                )
            if self.svf_save_8bit:
                dat.write("\t\t>> Output 8bit file:\n")
                dat.write(
                    "\t\tasvf_bytscl=\t\t({}, {}, {})\n".format(
                        self.asvf_bytscl[0], self.asvf_bytscl[1], self.asvf_bytscl[2]
                    )
                )
                dat.write(
                    "\t\t\t{}\n".format(
                        os.path.abspath(
                            os.path.join(
                                log_dir, self.get_asvf_file_name(dem_path, bit8=True)
                            )
                        )
                    )
                )
            dat.write("\n")
        if self.pos_opns_compute:
            dat.write("\tOpenness - Positive\n")
            dat.write("\t\tsvf_n_dir=\t\t{}\n".format(self.svf_n_dir))
            dat.write("\t\tsvf_noise=\t\t{}\n".format(self.svf_noise))
            dat.write("\t\tsvf_r_max=\t\t{}\n".format(self.svf_r_max))
            if self.svf_save_float:
                dat.write("\t\t>> Output file:\n")
                dat.write(
                    "\t\t\t{}\n".format(
                        os.path.abspath(
                            os.path.join(log_dir, self.get_opns_file_name(dem_path))
                        )
                    )
                )
            if self.svf_save_8bit:
                dat.write("\t\t>> Output 8bit file:\n")
                dat.write(
                    "\t\tpos_opns_bytscl=\t\t({}, {}, {})\n".format(
                        self.pos_opns_bytscl[0],
                        self.pos_opns_bytscl[1],
                        self.pos_opns_bytscl[2],
                    )
                )
                dat.write(
                    "\t\t\t{}\n".format(
                        os.path.abspath(
                            os.path.join(
                                log_dir, self.get_opns_file_name(dem_path, bit8=True)
                            )
                        )
                    )
                )
            dat.write("\n")
        if self.neg_opns_compute:
            dat.write("\tOpenness - Negative\n")
            dat.write("\t\tsvf_n_dir=\t\t{}\n".format(self.svf_n_dir))
            dat.write("\t\tsvf_noise=\t\t{}\n".format(self.svf_noise))
            dat.write("\t\tsvf_r_max=\t\t{}\n".format(self.svf_r_max))
            if self.neg_opns_save_float:
                dat.write("\t\t>> Output file:\n")
                dat.write(
                    "\t\t\t{}\n".format(
                        os.path.abspath(
                            os.path.join(log_dir, self.get_neg_opns_file_name(dem_path))
                        )
                    )
                )
            if self.neg_opns_save_8bit:
                dat.write("\t\t>> Output 8bit file:\n")
                dat.write(
                    "\t\tneg_opns_bytscl=\t\t({}, {}, {})\n".format(
                        self.neg_opns_bytscl[0],
                        self.neg_opns_bytscl[1],
                        self.neg_opns_bytscl[2],
                    )
                )
                dat.write(
                    "\t\t\t{}\n".format(
                        os.path.abspath(
                            os.path.join(
                                log_dir,
                                self.get_neg_opns_file_name(dem_path, bit8=True),
                            )
                        )
                    )
                )
            dat.write("\n")
        if self.sim_compute:
            dat.write("\tSky illumination\n")
            dat.write("\t\tsim_sky_mod=\t\t{}\n".format(self.sim_sky_mod))
            dat.write("\t\tsim_compute_shadow=\t\t{}\n".format(self.sim_compute_shadow))
            dat.write("\t\tsim_shadow_az=\t\t{}\n".format(self.sim_shadow_az))
            dat.write("\t\tsim_shadow_el=\t\t{}\n".format(self.sim_shadow_el))
            dat.write("\t\tsim_nr_dir=\t\t{}\n".format(self.sim_nr_dir))
            dat.write("\t\tsim_shadow_dist=\t\t{}\n".format(self.sim_shadow_dist))
            if self.sim_save_float:
                dat.write("\t\t>> Output file:\n")
                dat.write(
                    "\t\t\t{}\n".format(
                        os.path.abspath(
                            os.path.join(
                                log_dir, self.get_sky_illumination_file_name(dem_path)
                            )
                        )
                    )
                )
            if self.sim_save_8bit:
                dat.write("\t\t>> Output 8bit file:\n")
                dat.write(
                    "\t\tsim_bytscl=\t\t({}, {}, {})\n".format(
                        self.sim_bytscl[0], self.sim_bytscl[1], self.sim_bytscl[2]
                    )
                )
                dat.write(
                    "\t\t\t{}\n".format(
                        os.path.abspath(
                            os.path.join(
                                log_dir,
                                self.get_sky_illumination_file_name(
                                    dem_path, bit8=True
                                ),
                            )
                        )
                    )
                )
            dat.write("\n")
        if self.ld_compute:
            dat.write("\tLocal dominance\n")
            dat.write("\t\tld_rad_inc=\t\t{}\n".format(self.ld_rad_inc))
            dat.write("\t\tld_min_rad=\t\t{}\n".format(self.ld_min_rad))
            dat.write("\t\tld_max_rad=\t\t{}\n".format(self.ld_max_rad))
            dat.write("\t\tld_anglr_res=\t\t{}\n".format(self.ld_anglr_res))
            dat.write("\t\tld_observer_h=\t\t{}\n".format(self.ld_observer_h))
            if self.ld_save_float:
                dat.write("\t\t>> Output file:\n")
                dat.write(
                    "\t\t\t{}\n".format(
                        os.path.abspath(
                            os.path.join(
                                log_dir, self.get_local_dominance_file_name(dem_path)
                            )
                        )
                    )
                )
            if self.ld_save_8bit:
                dat.write("\t\t>> Output 8bit file:\n")
                dat.write(
                    "\t\tld_bytscl=\t\t({}, {}, {})\n".format(
                        self.ld_bytscl[0], self.ld_bytscl[1], self.ld_bytscl[2]
                    )
                )
                dat.write(
                    "\t\t\t{}\n".format(
                        os.path.abspath(
                            os.path.join(
                                log_dir,
                                self.get_local_dominance_file_name(dem_path, bit8=True),
                            )
                        )
                    )
                )
            dat.write("\n")
        if self.mstp_compute:
            dat.write("\tMulti-scale topographic position\n")
            dat.write(
                "\t\tmstp_local_scale=\t({}, {}, {})\n".format(
                    self.mstp_local_scale[0],
                    self.mstp_local_scale[1],
                    self.mstp_local_scale[2],
                )
            )
            dat.write(
                "\t\tmstp_meso_scale=\t({}, {}, {})\n".format(
                    self.mstp_meso_scale[0],
                    self.mstp_meso_scale[1],
                    self.mstp_meso_scale[2],
                )
            )
            dat.write(
                "\t\tmstp_broad_scale=\t({}, {}, {})\n".format(
                    self.mstp_broad_scale[0],
                    self.mstp_broad_scale[1],
                    self.mstp_broad_scale[2],
                )
            )
            dat.write("\t\tmstp_lightness=\t\t{}\n".format(self.mstp_lightness))
            dat.write("\t\t>> Output file:\n")
            dat.write(
                "\t\t\t{}\n".format(
                    os.path.abspath(
                        os.path.join(log_dir, self.get_mstp_file_name(dem_path))
                    )
                )
            )
            dat.write("\n")

        if compute_time is not None:
            dat.write("# Computation time: {:.3f}s".format(compute_time))
        dat.close()


def get_raster_arr(raster_path):
    """
    Reads raster from raster_path and returns its array(value) and resolution.

    Parameters
    ----------
    raster_path : str
        Path to raster

    Returns
    -------
    dict_out : dict
        Returns {"array": array, "resolution": (x_res, y_res), "no_data": no_data} : dict("array": np.array,
        "resolution": tuple(float, float), "no_data": float).
        Returns dictionary with keys: array, resolution and no_data. Key resolution is tuple where first element is x
        resolution and second is y resolution. Key no_data represent value of no data.
    """
    data_set = gdal.Open(raster_path)
    gt = data_set.GetGeoTransform()
    x_res = abs(gt[1])
    y_res = abs(-gt[5])
    bands = []
    no_data = data_set.GetRasterBand(
        1
    ).GetNoDataValue()  # we assume that all the bands have same no_data val
    if data_set.RasterCount == 1:  # only one band
        array = np.array(data_set.GetRasterBand(1).ReadAsArray())
        data_set = None
        return {"array": array, "resolution": (x_res, y_res), "no_data": no_data}
    else:  # multiple bands
        for i_band in range(data_set.RasterCount):
            i_band += 1
            band = np.array(data_set.GetRasterBand(i_band).ReadAsArray())
            if band is None:
                continue
            else:
                bands.append(band)
        data_set = None  # close dataset
        return {
            "array": np.array(bands),
            "resolution": (x_res, y_res),
            "no_data": no_data,
        }


def get_raster_size(raster_path, band=1):
    """Opens raster path and returns selected band size.

    Parameters
    ----------
    raster_path : str
        Path to raster.
    band : int
        Selected band number.

    Returns
    -------
    tuple(x_size, y_size)
    """
    data_set = gdal.Open(raster_path)  # Open dem raster
    band = data_set.GetRasterBand(band)
    x_size = band.XSize  # number of columns
    y_size = band.YSize  # number of rows
    del band
    data_set = None  # close data_set
    return x_size, y_size


def save_raster(
    src_raster_path, out_raster_path, out_raster_arr: np.ndarray, no_data=None, e_type=6
):
    """Saves raster array (out_rast_arr) to out_raster_path (GTiff), using src_rast_path information.

    Parameters
    ----------
    src_raster_path : str
        Path to source raster.
    out_raster_path : str
        Path to new file, where to save raster (GTiff).
    out_raster_arr : np.array (2D - one band, 3D - multiple bands)
        Array with raster data.
    no_data : float
        Value that represents no data pixels.
    e_type : GDALDataType
        https://gdal.org/api/raster_c_api.html#_CPPv412GDALDataType, (GDT_Float32 = 6, GDT_UInt8 = 1, ...)
    """
    src_data_set = gdal.Open(src_raster_path)
    gtiff_driver = gdal.GetDriverByName("GTiff")
    if len(out_raster_arr.shape) == 2:  # 2D array, one band
        out_data_set = gtiff_driver.Create(
            out_raster_path,
            xsize=out_raster_arr.shape[1],
            ysize=out_raster_arr.shape[0],
            bands=1,
            eType=e_type,  # eType: 6 = GDT_Float32
            options=["COMPRESS=LZW"],
        )
        out_data_set.SetProjection(src_data_set.GetProjection())
        out_data_set.SetGeoTransform(src_data_set.GetGeoTransform())
        out_data_set.GetRasterBand(1).WriteArray(out_raster_arr)
        if no_data is not None:
            out_data_set.GetRasterBand(1).SetNoDataValue(no_data)

    elif len(out_raster_arr.shape) == 3:  # 3D array, more bands
        out_data_set = gtiff_driver.Create(
            out_raster_path,
            xsize=out_raster_arr.shape[2],
            ysize=out_raster_arr.shape[1],
            bands=out_raster_arr.shape[0],
            eType=e_type,  # eType: 6 = GDT_Float32
            options=["COMPRESS=LZW"],
        )
        out_data_set.SetProjection(src_data_set.GetProjection())
        out_data_set.SetGeoTransform(src_data_set.GetGeoTransform())
        for i_band in range(out_raster_arr.shape[0]):
            out_data_set.GetRasterBand(i_band + 1).WriteArray(
                out_raster_arr[i_band, :, :]
            )
        if no_data is not None:
            out_data_set.GetRasterBand(1).SetNoDataValue(no_data)
    else:
        raise Exception(
            "rvt.default.save_raster: You have to input 2D or 3D numpy array!"
        )
    out_data_set.FlushCache()
    src_data_set = None  # Close source data set
    out_data_set = None  # Close output data set
