"""
Relief Visualization Toolbox – Visualization Functions

Contains all default values for visualisation functions, which can be changed.
Allows computing from rvt.vis with using defined default values and saving
output rasters with default names (dependent on default values).

Credits:
    Žiga Kokalj (ziga.kokalj@zrc-sazu.si)
    Krištof Oštir (kristof.ostir@fgg.uni-lj.si)
    Klemen Zakšek
    Klemen Čotar
    Maja Somrak
    Žiga Maroh

Copyright:
    2010-2020 Research Centre of the Slovenian Academy of Sciences and Arts
    2016-2020 University of Ljubljana, Faculty of Civil and Geodetic Engineering
"""

import warnings
import rvt.vis
import os
import gdal
import numpy as np
import json


class DefaultValues:
    """
    Class which define layer for blending. BlenderLayer is basic element in BlenderLayers.layers list.

    Attributes
    ----------
    overwrite : bool
        When saving visualisation functions and file already exists, if 0 it doesn't compute it, if 1 it overwrites it.
    ve_factor : float
        For all vis functions. Vertical exaggeration.
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
    sim_samp_pnts : int
        Sky illumination. Number of sampling points [250 or 500].
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
    """

    def __init__(self):
        self.overwrite = 0  # (0=False, 1=True)
        self.ve_factor = 1
        # slope gradient
        self.slp_compute = False
        self.slp_output_units = "degree"
        # hillshade
        self.hs_compute = True
        self.hs_sun_azi = 315
        self.hs_sun_el = 35
        # multi hillshade
        self.mhs_compute = False
        self.mhs_nr_dir = 16
        self.mhs_sun_el = 35
        # simple local relief model
        self.slrm_compute = False
        self.slrm_rad_cell = 20
        # sky view factor
        self.svf_compute = False
        self.svf_n_dir = 16
        self.svf_r_max = 10
        self.svf_noise = 0
        # anisotropic sky-view factor
        self.asvf_compute = False
        self.asvf_dir = 315
        self.asvf_level = 1
        # positive openness
        self.pos_opns_compute = False
        # negative openness
        self.neg_opns_compute = False
        # sky_illum
        self.sim_compute = False
        self.sim_sky_mod = "overcast"
        self.sim_samp_pnts = 250
        self.sim_shadow_dist = 100
        self.sim_shadow_az = 315
        self.sim_shadow_el = 35
        # local dominance
        self.ld_compute = False
        self.ld_min_rad = 10
        self.ld_max_rad = 20
        self.ld_rad_inc = 1
        self.ld_anglr_res = 15
        self.ld_observer_h = 1.7

    def save_default_to_file(self, file_path=None):
        """Saves default attributes into .json file."""
        data = {"default_settings": {"overwrite": {
                                    "value": self.overwrite,
                                    "description": "When saving visualisation functions and file already exists, if 0 "
                                    "it doesn't compute it, if 1 it overwrites it."
                                    },
                                    "ve_factor": {
                                    "value": self.ve_factor,
                                    "description": "Vertical exaggeration."},
                                     "Hillshade": {
                                         "hs_compute": {"value": self.hs_compute,
                                                        "description": "If compute Hillshade. Parameter for GUIs."},
                                         "hs_sun_azi": {"value": self.hs_sun_azi,
                                                        "description": "Solar azimuth angle (clockwise from North) in "
                                                                       "degrees."},
                                         "hs_sun_el": {"value": self.hs_sun_el,
                                                       "description": "Solar vertical angle (above the horizon) in "
                                                                      "degrees."}
                                     },
                                     "Multiple directions hillshade": {
                                         "mhs_compute": {"value": self.mhs_compute,
                                                         "description": "If compute Multiple directions hillshade."
                                                                        " Parameter for GUIs."},
                                         "mhs_nr_dir": {"value": self.mhs_nr_dir,
                                                        "description": "Number of solar azimuth angles (clockwise "
                                                                       "from North)."},
                                         "mhs_sun_el": {"value": self.mhs_sun_el,
                                                        "description": "Solar vertical angle (above the horizon) in "
                                                                       "degrees."}
                                     },
                                     "Slope gradient": {
                                         "slp_compute": {"value": self.slp_compute,
                                                         "description": "If compute Slope. Parameter for GUIs."},
                                         "slp_output_units": {"value": self.slp_output_units,
                                                              "description": "Slope output units [radian, degree, "
                                                                             "percent]."}
                                     },
                                     "Simple local relief model": {
                                         "slrm_compute": {"value": self.slrm_compute,
                                                          "description": "If compute Simple local relief model. "
                                                                         "Parameter for GUIs."},
                                         "slrm_rad_cell": {"value": self.slrm_rad_cell,
                                                           "description": "Radius for trend assessment in pixels."}
                                     },
                                     "Sky-View Factor": {
                                         "svf_compute": {"value": self.svf_compute,
                                                         "description": "If compute Sky-View Factor."
                                                                        " Parameter for GUIs."},
                                         "svf_n_dir": {"value": self.svf_n_dir, "description": "Number of directions."},
                                         "svf_r_max": {"value": self.svf_r_max, "description": "Maximal search "
                                                                                               "radious in pixels."},
                                         "svf_noise": {"value": self.svf_noise,
                                                       "description": "The level of noise remove [0-don't remove, "
                                                                      "1-low, 2-med, 3-high]."}
                                     },
                                     "Anisotropic Sky-View Factor": {
                                         "asvf_compute": {"value": self.asvf_compute,
                                                          "description": "If compute Anisotropic Sky-View Factor."
                                                                         " Parameter for GUIs."},
                                         "asvf_dir": {"value": self.asvf_dir,
                                                      "description": "Direction of anisotropy in degrees."},
                                         "asvf_level": {"value": self.asvf_level,
                                                        "description": "Level of anisotropy [1-low, 2-high]."}
                                     },
                                     "Openness - Positive": {
                                        "pos_opns_compute": {"value": self.pos_opns_compute,
                                                             "description": "If compute Openness - Positive. "
                                                                            "Parameter for GUIs."}
                                    },
                                     "Openness - Negative": {
                                        "neg_opns_compute": {"value": self.neg_opns_compute,
                                                             "description": "If compute Openness - Negative. "
                                                                            "Parameter for GUIs."}
                                    },
                                     "Sky illumination": {
                                         "sim_compute": {"value": self.sim_compute,
                                                         "description": "If compute Sky illumination. Parameter for "
                                                                        "GUIs."},
                                         "sim_sky_mod": {"value": self.sim_sky_mod,
                                                         "description": "Sky model [overcast, uniform]."},
                                         "sim_samp_pnts": {"value": self.sim_samp_pnts,
                                                           "description": "Number of sampling points [250 or 500]."},
                                         "sim_shadow_dist": {"value": self.sim_shadow_dist,
                                                             "description": "Max shadow modeling distance in pixels."},
                                         "sim_shadow_az": {"value": self.sim_shadow_az, "description": "Shadow "
                                                                                                       "azimuth in "
                                                                                                       "degrees."},
                                         "sim_shadow_el": {"value": self.sim_shadow_el, "description": "Shadow "
                                                                                                       "elevation in "
                                                                                                       "degrees."}
                                     },
                                     "Local dominance": {
                                         "ld_compute": {"value": self.ld_compute,
                                                        "description": "If compute Local dominance. Parameter for "
                                                                       "GUIs."},
                                         "ld_min_rad": {"value": self.ld_min_rad,
                                                        "description": "Minimum radial distance (in pixels) at which "
                                                                       "the algorithm starts with visualization "
                                                                       "computation."},
                                         "ld_max_rad": {"value": self.ld_max_rad,
                                                        "description": "Maximum radial distance (in pixels) at which "
                                                                       "the algorithm ends with visualization "
                                                                       "computation."},
                                         "ld_rad_inc": {"value": self.ld_rad_inc, "description": "Radial distance "
                                                                                                 "steps in pixels."},
                                         "ld_anglr_res": {"value": self.ld_anglr_res,
                                                          "description": "Angular step for determination of number of "
                                                                         "angular directions."},
                                         "ld_observer_h": {"value": self.ld_observer_h,
                                                           "description": "Height at which we observe the terrain."}
                                     }}}
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
                    elif parameter_name == "number_points":
                        self.sim_samp_pnts = int(parameter_value)
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
                    warnings.warn("rvt.default.read_default_from_file: Wrong line '{}'".format(line))
                    continue
            dat.close()
        elif extension == ".json":
            dat = open(file_path, "r")
            data = json.load(dat)
            default_data = data["default_settings"]
            self.overwrite = int(default_data["overwrite"]["value"])
            self.ve_factor = float(default_data["ve_factor"]["value"])
            # Slope gradient
            self.slp_compute = int(default_data["Slope gradient"]["slp_compute"]["value"])
            self.slp_output_units = str(default_data["Slope gradient"]["slp_output_units"]["value"])
            # Hillshade
            self.hs_compute = int(default_data["Hillshade"]["hs_compute"]["value"])
            self.hs_sun_azi = int(default_data["Hillshade"]["hs_sun_azi"]["value"])
            self.hs_sun_el = int(default_data["Hillshade"]["hs_sun_el"]["value"])
            # Multiple directions hillshade
            self.mhs_compute = int(default_data["Multiple directions hillshade"]["mhs_compute"]["value"])
            self.mhs_nr_dir = int(default_data["Multiple directions hillshade"]["mhs_nr_dir"]["value"])
            self.mhs_sun_el = int(default_data["Multiple directions hillshade"]["mhs_sun_el"]["value"])
            # Simple local relief model
            self.slrm_compute = int(default_data["Simple local relief model"]["slrm_compute"]["value"])
            self.slrm_rad_cell = int(default_data["Simple local relief model"]["slrm_rad_cell"]["value"])
            # Sky-View Factor
            self.svf_compute = int(default_data["Sky-View Factor"]["svf_compute"]["value"])
            self.svf_n_dir = int(default_data["Sky-View Factor"]["svf_n_dir"]["value"])
            self.svf_r_max = int(default_data["Sky-View Factor"]["svf_r_max"]["value"])
            self.svf_noise = int(default_data["Sky-View Factor"]["svf_noise"]["value"])
            # Anisotropic Sky-View Factor
            self.asvf_compute = int(default_data["Anisotropic Sky-View Factor"]["asvf_compute"]["value"])
            self.asvf_dir = int(default_data["Anisotropic Sky-View Factor"]["asvf_dir"]["value"])
            self.asvf_level = int(default_data["Anisotropic Sky-View Factor"]["asvf_level"]["value"])
            # Openness - Positive
            self.pos_opns_compute = int(default_data["Openness - Positive"]["pos_opns_compute"]["value"])
            # Openness - Negative
            self.neg_opns_compute = int(default_data["Openness - Negative"]["neg_opns_compute"]["value"])
            # Sky illumination
            self.sim_compute = int(default_data["Sky illumination"]["sim_compute"]["value"])
            self.sim_sky_mod = str(default_data["Sky illumination"]["sim_sky_mod"]["value"])
            self.sim_samp_pnts = int(default_data["Sky illumination"]["sim_samp_pnts"]["value"])
            self.sim_shadow_dist = int(default_data["Sky illumination"]["sim_shadow_dist"]["value"])
            self.sim_shadow_az = int(default_data["Sky illumination"]["sim_shadow_az"]["value"])
            self.sim_shadow_el = int(default_data["Sky illumination"]["sim_shadow_el"]["value"])
            # Local dominance
            self.ld_compute = int(default_data["Local dominance"]["ld_compute"]["value"])
            self.ld_min_rad = int(default_data["Local dominance"]["ld_min_rad"]["value"])
            self.ld_max_rad = int(default_data["Local dominance"]["ld_max_rad"]["value"])
            self.ld_rad_inc = int(default_data["Local dominance"]["ld_rad_inc"]["value"])
            self.ld_anglr_res = int(default_data["Local dominance"]["ld_anglr_res"]["value"])
            self.ld_observer_h = float(default_data["Local dominance"]["ld_observer_h"]["value"])
            dat.close()

    def get_hillshade_file_name(self, dem_path):
        """Returns Hillshade name, dem name (from dem_path) with added hillshade parameters."""
        dem_name = os.path.basename(dem_path).split(".")[0]  # base name without extension
        return "{}_HS_A{}_H{}.tif".format(dem_name, self.hs_sun_azi, self.hs_sun_el)

    def get_hillshade_path(self, dem_path):
        """Returns path to Hillshade. Generates hillshade name (uses default attributes and dem name from dem_path) and
        adds dem directory (dem_path) to it."""
        return os.path.join(os.path.dirname(dem_path), self.get_hillshade_file_name(dem_path))

    def get_slope_file_name(self, dem_path):
        """Returns Slope name, dem name (from dem_path) with added slope parameters."""
        dem_name = os.path.basename(dem_path).split(".")[0]  # base name without extension
        return "{}_SLOPE.tif".format(dem_name)

    def get_slope_path(self, dem_path):
        """Returns path to slope. Generates slope name and adds dem directory (dem_path) to it."""
        return os.path.join(os.path.dirname(dem_path), self.get_slope_file_name(dem_path))

    def get_multi_hillshade_file_name(self, dem_path):
        """Returns Multiple directions hillshade name, dem name (from dem_path) with added
        multi hillshade parameters."""
        dem_name = os.path.basename(dem_path).split(".")[0]  # base name without extension
        return "{}_MULTI-HS_D{}_H{}.tif".format(dem_name, self.mhs_nr_dir, self.mhs_sun_el)

    def get_multi_hillshade_path(self, dem_path):
        """Returns path to Multiple directions hillshade. Generates multi hillshade name (uses default attributes and
        dem name from dem_path) and adds dem directory (dem_path) to it."""
        return os.path.join(os.path.dirname(dem_path), self.get_multi_hillshade_file_name(dem_path))

    def get_slrm_file_name(self, dem_path):
        """Returns Simple local relief model name, dem name (from dem_path) with added slrm parameters."""
        dem_name = os.path.basename(dem_path).split(".")[0]  # base name without extension
        return "{}_SLRM_R{}.tif".format(dem_name, self.slrm_rad_cell)

    def get_slrm_path(self, dem_path):
        """Returns path to Simple local relief model. Generates slrm name (uses default attributes and dem name from
        dem_path) and adds dem directory (dem_path) to it."""
        return os.path.join(os.path.dirname(dem_path), self.get_slrm_file_name(dem_path))

    def get_svf_file_name(self, dem_path):
        """Returns Sky-view factor name, dem name (from dem_path) with added svf parameters."""
        dem_name = os.path.basename(dem_path).split(".")[0]  # base name without extension
        return "{}_SVF_R{}_D{}.tif".format(dem_name, self.svf_r_max, self.svf_n_dir)

    def get_svf_path(self, dem_path):
        """Returns path to Sky-view factor. Generates svf name (uses default attributes and dem name from
        dem_path) and adds dem directory (dem_path) to it."""
        return os.path.join(os.path.dirname(dem_path), self.get_svf_file_name(dem_path))

    def get_asvf_file_name(self, dem_path):
        """Returns Anisotropic Sky-view factor name, dem name (from dem_path) with added asvf parameters."""
        dem_name = os.path.basename(dem_path).split(".")[0]  # base name without extension
        return "{}_SVF-A_R{}_D{}_A{}.tif".format(dem_name, self.svf_r_max, self.svf_n_dir, self.asvf_dir)

    def get_asvf_path(self, dem_path):
        """Returns path to Anisotropic Sky-view factor. Generates asvf name (uses default attributes and dem name from
        dem_path) and adds dem directory (dem_path) to it."""
        return os.path.join(os.path.dirname(dem_path), self.get_asvf_file_name(dem_path))

    def get_opns_file_name(self, dem_path):
        """Returns Positive Openness name, dem name (from dem_path) with added pos opns parameters."""
        dem_name = os.path.basename(dem_path).split(".")[0]  # base name without extension
        return "{}_OPEN-POS_R{}_D{}.tif".format(dem_name, self.svf_r_max, self.svf_n_dir)

    def get_opns_path(self, dem_path):
        """Returns path to Positive Openness. Generates pos opns name (uses default attributes and dem name from
        dem_path) and adds dem directory (dem_path) to it."""
        return os.path.join(os.path.dirname(dem_path), self.get_opns_file_name(dem_path))

    def get_neg_opns_file_name(self, dem_path):
        """Returns Negative Openness name, dem name (from dem_path) with added neg opns parameters."""
        dem_name = os.path.basename(dem_path).split(".")[0]  # base name without extension
        return "{}_OPEN-NEG_R{}_D{}.tif".format(dem_name, self.svf_r_max, self.svf_n_dir)

    def get_neg_opns_path(self, dem_path):
        """Returns path to Negative Openness. Generates pos neg name (uses default attributes and dem name from
        dem_path) and adds dem directory (dem_path) to it."""
        return os.path.join(os.path.dirname(dem_path), self.get_neg_opns_file_name(dem_path))

    def get_sky_illumination_file_name(self, dem_path):
        """Returns Sky illumination name, dem name (from dem_path) with added sim parameters."""
        dem_name = os.path.basename(dem_path).split(".")[0]  # base name without extension
        return "{}_SIM_{}_{}sp_{}px.tif".format(dem_name, self.sim_sky_mod, self.sim_samp_pnts, self.sim_shadow_dist)

    def get_sky_illumination_path(self, dem_path):
        """Returns path to Sky illumination. Generates sim name (uses default attributes and dem name from
        dem_path) and adds dem directory (dem_path) to it."""
        return os.path.join(os.path.dirname(dem_path), self.get_sky_illumination_file_name(dem_path))

    def get_local_dominance_file_name(self, dem_path):
        """Returns Local dominance name, dem name (from dem_path) with added ld parameters."""
        dem_name = os.path.basename(dem_path).split(".")[0]  # base name without extension
        return "{}_LD_R_M{}-{}_DI{}_A{}_OH{}.tif".format(dem_name, self.ld_min_rad, self.ld_max_rad,
                                                         self.ld_rad_inc, self.ld_anglr_res, self.ld_observer_h)

    def get_local_dominance_path(self, dem_path):
        """Returns path to Local dominance. Generates ld name (uses default attributes and dem name from
        dem_path) and adds dem directory (dem_path) to it."""
        return os.path.join(os.path.dirname(dem_path), self.get_local_dominance_file_name(dem_path))

    # get slope, aspect dict
    def get_slope(self, dem_arr, resolution_x, resolution_y):
        dict_slp_asp = rvt.vis.slope_aspect(dem=dem_arr, resolution_x=resolution_x, resolution_y=resolution_y,
                                            ve_factor=self.ve_factor, output_units=self.slp_output_units)
        return dict_slp_asp

    # save default slope gradient
    def save_slope(self, dem_path, custom_dir=None):
        """Calculates and saves Slope from dem (dem_path) with default parameters. If custom_dir is None it saves
        in dem directory else in custom_dir. If path to file already exists we can overwrite file (overwrite=0) or
        not (overwrite=1)."""
        if custom_dir is None:
            slope_path = self.get_slope_path(dem_path)
        else:
            slope_path = os.path.join(custom_dir, self.get_slope_file_name(dem_path))
        if os.path.isfile(slope_path) and not self.overwrite:  # if file already exists and overwrite=0
            return 0
        if not os.path.isfile(dem_path):
            raise Exception("rvt.default.DefaultValues.save_slope: dem_path doesn't exist!")
        dict_arr_res = get_raster_arr(raster_path=dem_path)
        dem_arr = dict_arr_res["array"]
        x_res = dict_arr_res["resolution"][0]
        y_res = dict_arr_res["resolution"][1]
        dict_slp_asp = self.get_slope(dem_arr=dem_arr, resolution_x=x_res, resolution_y=y_res)
        slope_arr = dict_slp_asp["slope"].astype('float32')
        save_raster(src_raster_path=dem_path, out_raster_path=slope_path, out_raster_arr=slope_arr)
        return 1

    # get hillshade array
    def get_hillshade(self, dem_arr, resolution_x, resolution_y):
        hillshade_arr = rvt.vis.hillshade(dem=dem_arr, resolution_x=resolution_x, resolution_y=resolution_y,
                                          sun_azimuth=self.hs_sun_azi, sun_elevation=self.hs_sun_el,
                                          ve_factor=self.ve_factor)
        return hillshade_arr

    # save default hillshade
    def save_hillshade(self, dem_path, custom_dir=None):
        """Calculates and saves Hillshade from dem (dem_path) with default parameters. If custom_dir is None it saves
        in dem directory else in custom_dir. If path to file already exists we can overwrite file (overwrite=1)
        or not (overwrite=0)."""
        if custom_dir is None:
            hillshade_path = self.get_hillshade_path(dem_path)
        else:
            hillshade_path = os.path.join(custom_dir, self.get_hillshade_file_name(dem_path))
        if os.path.isfile(hillshade_path) and not self.overwrite:  # if file already exists and overwrite=0
            return 0
        if not os.path.isfile(dem_path):
            raise Exception("rvt.default.DefaultValues.save_hillshade: dem_path doesn't exist!")
        dict_arr_res = get_raster_arr(raster_path=dem_path)
        dem_arr = dict_arr_res["array"]
        x_res = dict_arr_res["resolution"][0]
        y_res = dict_arr_res["resolution"][1]
        hillshade_arr = self.get_hillshade(dem_arr=dem_arr, resolution_x=x_res, resolution_y=y_res).astype('float32')
        save_raster(src_raster_path=dem_path, out_raster_path=hillshade_path, out_raster_arr=hillshade_arr)
        return 1

    # get multi hillshade array
    def get_multi_hillshade(self, dem_arr, resolution_x, resolution_y):
        multi_hillshade_arr = rvt.vis.multi_hillshade(dem=dem_arr, resolution_x=resolution_x, resolution_y=resolution_y,
                                                      nr_directions=self.mhs_nr_dir, sun_elevation=self.mhs_sun_el,
                                                      ve_factor=self.ve_factor)
        return multi_hillshade_arr

    def save_multi_hillshade(self, dem_path, custom_dir=None):
        """Calculates and saves Multidirectional hillshade from dem (dem_path) with default parameters.
        If custom_dir is None it saves in dem directory else in custom_dir. If path to file already exists we can
        overwrite file (overwrite=1) or not (overwrite=0)."""
        if custom_dir is None:
            multi_hillshade_path = self.get_multi_hillshade_path(dem_path)
        else:
            multi_hillshade_path = os.path.join(custom_dir, self.get_multi_hillshade_file_name(dem_path))
        if os.path.isfile(multi_hillshade_path) and not self.overwrite:  # if file already exists and overwrite=0
            return 0
        if not os.path.isfile(dem_path):
            raise Exception("rvt.default.DefaultValues.save_multi_hillshade: dem_path doesn't exist!")
        dict_arr_res = get_raster_arr(raster_path=dem_path)
        dem_arr = dict_arr_res["array"]
        x_res = dict_arr_res["resolution"][0]
        y_res = dict_arr_res["resolution"][1]
        multi_hillshade_arr = self.get_multi_hillshade(dem_arr=dem_arr, resolution_x=x_res,
                                                       resolution_y=y_res).astype('float32')
        save_raster(src_raster_path=dem_path, out_raster_path=multi_hillshade_path, out_raster_arr=multi_hillshade_arr)
        return 1

    def get_slrm(self, dem_arr):
        slrm_arr = rvt.vis.slrm(dem=dem_arr, radius_cell=self.slrm_rad_cell, ve_factor=self.ve_factor)
        return slrm_arr

    def save_slrm(self, dem_path, custom_dir=None):
        """Calculates and saves Simple local relief model from dem (dem_path) with default parameters.
        If custom_dir is None it saves in dem directory else in custom_dir. If path to file already exists we can
        overwrite file (overwrite=1) or not (overwrite=0)."""
        if custom_dir is None:
            slrm_path = self.get_slrm_path(dem_path)
        else:
            slrm_path = os.path.join(custom_dir, self.get_slrm_file_name(dem_path))
        if os.path.isfile(slrm_path) and not self.overwrite:  # if file already exists and overwrite=0
            return 0
        if not os.path.isfile(dem_path):
            raise Exception("rvt.default.DefaultValues.save_slrm: dem_path doesn't exist!")
        dict_arr_res = get_raster_arr(raster_path=dem_path)
        dem_arr = dict_arr_res["array"]
        slrm_arr = self.get_slrm(dem_arr=dem_arr).astype('float32')
        save_raster(src_raster_path=dem_path, out_raster_path=slrm_path, out_raster_arr=slrm_arr)
        return 1

    # get svf, asvf, opns dict
    def get_sky_view_factor(self, dem_arr, resolution, compute_svf=True, compute_asvf=False, compute_opns=False):
        dict_svf_asvf_opns = rvt.vis.sky_view_factor(dem=dem_arr, resolution=resolution, compute_svf=compute_svf,
                                                     compute_opns=compute_opns, compute_asvf=compute_asvf,
                                                     svf_n_dir=self.svf_n_dir, svf_r_max=self.svf_r_max,
                                                     svf_noise=self.svf_noise, asvf_dir=self.asvf_dir,
                                                     asvf_level=self.asvf_level, ve_factor=self.ve_factor)
        return dict_svf_asvf_opns

    def save_sky_view_factor(self, dem_path, save_svf=True, save_asvf=False, save_opns=False, custom_dir=None):
        """Calculates and saves Sky-view factor(save_svf=True), Anisotropic Sky-view factor(save_asvf=True) and
        Positive Openness(save_opns=True) from dem (dem_path) with default parameters.
        If custom_dir is None it saves in dem directory else in custom_dir. If path to file already exists we can
        overwrite file (overwrite=1) or not (overwrite=0)."""
        svf_path = ""
        asvf_path = ""
        opns_path = ""
        if custom_dir is None:
            if save_svf:
                svf_path = self.get_svf_path(dem_path)
            if save_asvf:
                asvf_path = self.get_asvf_path(dem_path)
            if save_opns:
                opns_path = self.get_opns_path(dem_path)
        else:
            if save_svf:
                svf_path = os.path.join(custom_dir, self.get_svf_file_name(dem_path))
            if save_asvf:
                asvf_path = os.path.join(custom_dir, self.get_asvf_file_name(dem_path))
            if save_opns:
                opns_path = os.path.join(custom_dir, self.get_opns_file_name(dem_path))
        if os.path.isfile(svf_path) and os.path.isfile(asvf_path) and os.path.isfile(opns_path)\
                and not self.overwrite:  # if file already exists and overwrite=0
            return 0
        if not os.path.isfile(dem_path):
            raise Exception("rvt.default.DefaultValues.save_sky_view_factor: dem_path doesn't exist!")
        dict_arr_res = get_raster_arr(raster_path=dem_path)
        dem_arr = dict_arr_res["array"]
        x_res = dict_arr_res["resolution"][0]
        y_res = dict_arr_res["resolution"][1]
        if x_res != y_res:
            raise Exception("rvt.default.DefaultValues.save_sky_view_factor: dem resolution is not the same in x and y"
                            " directions!")
        dict_svf_asvf_opns = self.get_sky_view_factor(dem_arr=dem_arr, resolution=x_res, compute_svf=save_svf,
                                                      compute_asvf=save_asvf, compute_opns=save_opns)
        if save_svf:
            if os.path.isfile(svf_path) and self.overwrite:  # svf file exists and overwrite=1
                pass
            else:  # svf_path, file doesn't exists or exists and overwrite=1
                save_raster(src_raster_path=dem_path, out_raster_path=svf_path,
                            out_raster_arr=dict_svf_asvf_opns["svf"].astype('float32'))
        if save_asvf:
            if os.path.isfile(svf_path) and self.overwrite:  # svf file exists and overwrite=1
                pass
            else:  # asvf_path, file doesn't exists or exists and overwrite=1
                save_raster(src_raster_path=dem_path, out_raster_path=asvf_path,
                            out_raster_arr=dict_svf_asvf_opns["asvf"].astype('float32'))
        if save_opns:
            if os.path.isfile(svf_path) and self.overwrite:  # svf file exists and overwrite=1
                pass
            else:  # opns_path, file doesn't exists or exists and overwrite=1
                save_raster(src_raster_path=dem_path, out_raster_path=opns_path,
                            out_raster_arr=dict_svf_asvf_opns["opns"].astype('float32'))
        return 1

    def get_neg_opns(self, dem_arr, resolution):
        dem_arr = -1 * dem_arr
        dict_neg_opns = rvt.vis.sky_view_factor(dem=dem_arr, resolution=resolution, svf_n_dir=self.svf_n_dir,
                                                svf_r_max=self.svf_r_max, svf_noise=self.svf_noise,
                                                compute_svf=False, compute_asvf=False, compute_opns=True,
                                                ve_factor=self.ve_factor)
        neg_opns_arr = dict_neg_opns["opns"]
        return neg_opns_arr

    def save_neg_opns(self, dem_path, custom_dir=None):
        """Calculates and saves Negative Openness from dem (dem_path) with default parameters. If custom_dir is None
        it saves in dem directory else in custom_dir. If path to file already exists we can
        overwrite file (overwrite=1) or not (overwrite=0)."""
        if custom_dir is None:
            neg_opns_path = self.get_neg_opns_path(dem_path)
        else:
            neg_opns_path = os.path.join(custom_dir, self.get_neg_opns_file_name(dem_path))
        if os.path.isfile(neg_opns_path) and not self.overwrite:  # if file already exists and overwrite=0
            return 0
        if not os.path.isfile(dem_path):
            raise Exception("rvt.default.DefaultValues.save_neg_opns: dem_path doesn't exist!")
        dict_arr_res = get_raster_arr(raster_path=dem_path)
        dem_arr = dict_arr_res["array"]
        x_res = dict_arr_res["resolution"][0]
        y_res = dict_arr_res["resolution"][1]
        if x_res != y_res:
            raise Exception("rvt.default.DefaultValues.save_neg_opns: dem resolution is not the same in x and y"
                            " directions!")
        neg_opns_arr = self.get_neg_opns(dem_arr=dem_arr, resolution=x_res).astype('float32')
        save_raster(src_raster_path=dem_path, out_raster_path=neg_opns_path, out_raster_arr=neg_opns_arr)
        return 1

    def get_sky_illumination(self, dem_arr, resolution):
        sky_illumination_arr = rvt.vis.sky_illumination(dem=dem_arr, resolution=resolution, sky_model=self.sim_sky_mod,
                                                        sampling_points=self.sim_samp_pnts,
                                                        shadow_dist=self.sim_shadow_dist,
                                                        shadow_az=self.sim_shadow_az, shadow_el=self.sim_shadow_el,
                                                        ve_factor=self.ve_factor)
        return sky_illumination_arr

    def save_sky_illumination(self, dem_path, custom_dir=None):
        """Calculates and saves Sky illumination from dem (dem_path) with default parameters. If custom_dir is None
        it saves in dem directory else in custom_dir. If path to file already exists we can
        overwrite file (overwrite=1) or not (overwrite=0)."""
        if custom_dir is None:
            sky_illumination_path = self.get_sky_illumination_path(dem_path)
        else:
            sky_illumination_path = os.path.join(custom_dir, self.get_sky_illumination_file_name(dem_path))
        if os.path.isfile(sky_illumination_path) and not self.overwrite:  # if file already exists and overwrite=0
            return 0
        if not os.path.isfile(dem_path):
            raise Exception("rvt.default.DefaultValues.save_sky_illumination: dem_path doesn't exist!")
        dict_arr_res = get_raster_arr(raster_path=dem_path)
        dem_arr = dict_arr_res["array"]
        x_res = dict_arr_res["resolution"][0]
        y_res = dict_arr_res["resolution"][1]
        if x_res != y_res:
            raise Exception("rvt.default.DefaultValues.save_sky_illumination: dem resolution is not the same in x and y"
                            " directions!")
        sky_illumination_arr = self.get_sky_illumination(dem_arr=dem_arr, resolution=x_res).astype('float32')
        save_raster(src_raster_path=dem_path, out_raster_path=sky_illumination_path,
                    out_raster_arr=sky_illumination_arr)
        return 1

    def get_local_dominance(self, dem_arr):
        local_dominance_arr = rvt.vis.local_dominance(dem=dem_arr, min_rad=self.ld_min_rad, max_rad=self.ld_max_rad,
                                                      rad_inc=self.ld_rad_inc, angular_res=self.ld_anglr_res,
                                                      observer_height=self.ld_observer_h, ve_factor=self.ve_factor)
        return local_dominance_arr

    def save_local_dominance(self, dem_path, custom_dir=None):
        """Calculates and saves Local dominance from dem (dem_path) with default parameters. If custom_dir is None
        it saves in dem directory else in custom_dir. If path to file already exists we can
        overwrite file (overwrite=1) or not (overwrite=0)."""
        if custom_dir is None:
            local_dominance_path = self.get_local_dominance_path(dem_path)
        else:
            local_dominance_path = os.path.join(custom_dir, self.get_local_dominance_file_name(dem_path))
        if os.path.isfile(local_dominance_path) and not self.overwrite:  # if file already exists and overwrite=0
            return 0
        if not os.path.isfile(dem_path):
            raise Exception("rvt.default.DefaultValues.save_local_dominance: dem_path doesn't exist!")
        dict_arr_res = get_raster_arr(raster_path=dem_path)
        dem_arr = dict_arr_res["array"]
        local_dominance_arr = self.get_local_dominance(dem_arr=dem_arr).astype('float32')
        save_raster(src_raster_path=dem_path, out_raster_path=local_dominance_path, out_raster_arr=local_dominance_arr)
        return 1

    def save_visualizations(self, dem_path, custom_dir=None, sav_slope=True, sav_hillshade=True,
                            sav_mulit_hillshade=True, sav_slrm=True, sav_svf=True, sav_asvf=True, sav_opns=True,
                            sav_neg_opns=True, sav_sky_illumination=True, sav_local_dominance=True):
        if sav_slope:
            self.save_slope(dem_path, custom_dir=custom_dir)
        if sav_hillshade:
            self.save_hillshade(dem_path, custom_dir=custom_dir)
        if sav_mulit_hillshade:
            self.save_multi_hillshade(dem_path, custom_dir=custom_dir)
        if sav_slrm:
            self.save_slrm(dem_path, custom_dir=custom_dir)
        if sav_svf or sav_asvf or sav_opns:
            self.save_sky_view_factor(dem_path, save_svf=sav_svf, save_asvf=sav_asvf, save_opns=sav_opns,
                                      custom_dir=custom_dir)
        if sav_neg_opns:
            self.save_neg_opns(dem_path, custom_dir=custom_dir)
        if sav_sky_illumination:
            self.save_neg_opns(dem_path, custom_dir=custom_dir)
        if sav_local_dominance:
            self.save_local_dominance(dem_path, custom_dir=custom_dir)


def get_raster_arr(raster_path):
    """
    Reads raster from raster_path and returns its array(value) and resolution.

    Parameters
    ----------
    raster_path : str
        Path to raster

    Returns
    {"array": array, "resolution": (x_res, y_res)} : dict("array": np.array, "resolution": tuple(float, float))
        Returns dictionary with keys array and resolution, resolution is tuple where first element is x resolution and
        second is y resolution.
    -------

    """
    data_set = gdal.Open(raster_path)
    gt = data_set.GetGeoTransform()
    x_res = gt[1]
    y_res = -gt[5]
    bands = []
    if data_set.RasterCount == 1:  # only one band
        array = np.array(data_set.GetRasterBand(1).ReadAsArray())
        data_set = None
        return {"array": array, "resolution": (x_res, y_res)}
    else:  # multiple bands
        for i_band in range(data_set.RasterCount):
            i_band += 1
            band = np.array(data_set.GetRasterBand(i_band).ReadAsArray())
            if band is None:
                continue
            else:
                bands.append(band)
        data_set = None  # close dataset
        return {"array": np.array(bands), "resolution": (x_res, y_res)}


def save_raster(src_raster_path, out_raster_path, out_raster_arr, e_type=6):
    """Saves raster array (out_rast_arr) to out_raster_path (GTiff), using src_rast_path information.

    Parameters
    ----------
    src_raster_path : str
        Path to source raster.
    out_raster_path : str
        Path to new file, where to save raster (GTiff).
    out_raster_arr : np.array (2D - one band, 3D - multiple bands)
        Array with raster data.
    e_type : GDALDataType
        https://gdal.org/api/raster_c_api.html#_CPPv412GDALDataType, (GDT_Float32 = 6, GDT_UInt16 = 2, ...)
    """
    src_data_set = gdal.Open(src_raster_path)
    gtiff_driver = gdal.GetDriverByName("GTiff")
    if len(out_raster_arr.shape) == 2:  # 2D array, one band
        out_data_set = gtiff_driver.Create(out_raster_path, xsize=out_raster_arr.shape[1],
                                           ysize=out_raster_arr.shape[0],
                                           bands=1,
                                           eType=e_type)  # eType: 6 = GDT_Float32
        out_data_set.SetProjection(src_data_set.GetProjection())
        out_data_set.SetGeoTransform(src_data_set.GetGeoTransform())
        out_data_set.GetRasterBand(1).WriteArray(out_raster_arr)

    elif len(out_raster_arr.shape) == 3:  # 3D array, more bands
        out_data_set = gtiff_driver.Create(out_raster_path, xsize=out_raster_arr.shape[2],
                                           ysize=out_raster_arr.shape[1],
                                           bands=out_raster_arr.shape[0],
                                           eType=6)  # eType: 6 = GDT_Float32
        out_data_set.SetProjection(src_data_set.GetProjection())
        out_data_set.SetGeoTransform(src_data_set.GetGeoTransform())
        for i_band in range(out_raster_arr.shape[0]):
            out_data_set.GetRasterBand(i_band + 1).WriteArray(out_raster_arr[i_band, :, :])

    out_data_set.FlushCache()
    src_data_set = None  # Close source data set
    out_data_set = None  # Close output data set
