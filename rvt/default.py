"""
Relief Visualization Toolbox Default Values for Visualization Functions

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
2010-2018 Space-SI
2016-2020 University of Ljubljana, Faculty of Civil and Geodetic Engineering
"""
"""
NAME:

PROJECT MANAGER:
    Žiga Kokalj (ziga.kokalj@zrc-sazu.si)

AUTHORS:
    Klemen Zakšek
    Krištof Oštir
    Klemen Čotar
    Maja Somrak
    Žiga Maroh

COPYRIGHT:
    Research Centre of the Slovenian Academy of Sciences and Arts
    Space-SI
    University of Ljubljana, Faculty of Civil and Geodetic Engineering
"""

import warnings
import rvt.vis
import os
import rasterio as rio
import numpy as np
import json


class DefaultValues():
    """
    Class which define layer for blending. BlenderLayer is basic element in BlenderLayers.layers list.

    Attributes
    ----------
    ve_factor : float
        For all vis functions. Vertical exaggeration.
    slp_output_units : str
        Slope. Output units [radian, degree, percent].
    hs_sun_azi : int
        Hillshade. Solar azimuth angle (clockwise from North) in degrees.
    hs_sun_el : int
        Hillshade. Solar vertical angle (above the horizon) in degrees.
    mhs_nr_dir : int
        Multi directional hillshade. Number of solar azimuth angles (clockwise from North).
    mhs_sun_el : int
        Multi directional hillshade. Solar vertical angle (above the horizon) in degrees.
    slrm_rad_cell : int
        Simple local relief model. Radius for trend assessment in pixels.
    svf_n_dir : int
        Sky-View Factor (Anisotropic Sky-View Factor, Openness). Number of directions.
    svf_r_max : int
        Sky-View Factor (Anisotropic Sky-View Factor, Openness). Maximal search radius in pixels.
    svf_noise : int
        Sky-View Factor (Anisotropic Sky-View Factor, Openness). The level of noise remove [0-don't remove, 1-low, 2-med, 3-high].
    asvf_dir : int
        Anisotropic Sky-View Factor. Direction of anisotropy in degrees.
    asvf_level : int
        Anisotropic Sky-View Factor. Level of anisotropy [1-low, 2-high].
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
    ld_min_rad : int
        Local dominance. Minimum radial distance (in pixels) at which the algorithm starts with visualization computation.
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
        self.ve_factor = 1
        # slope gradient
        self.slp_output_units = "degree"
        # hillshade
        self.hs_sun_azi = 315
        self.hs_sun_el = 35
        # multi hillshade
        self.mhs_nr_dir = 16
        self.mhs_sun_el = 35
        # simple local relief model
        self.slrm_rad_cell = 20
        # sky_view_factor
        self.svf_n_dir = 16
        self.svf_r_max = 10
        self.svf_noise = 0
        self.asvf_dir = 315
        self.asvf_level = 1
        # sky_illum
        self.sim_sky_mod = "overcast"
        self.sim_samp_pnts = 250
        self.sim_shadow_dist = 100
        self.sim_shadow_az = 315
        self.sim_shadow_el = 35
        # local dominance
        self.ld_min_rad = 10
        self.ld_max_rad = 20
        self.ld_rad_inc = 1
        self.ld_anglr_res = 15
        self.ld_observer_h = 1.7

    def save_default_to_file(self, file_path=None):
        """Saves default attributes into .json file."""
        data = {"default_settings": {"ve_factor": {
                                    "value" : self.ve_factor,
                                    "description" : "Vertical exaggeration."},
                                     "Hillshade": {
                                         "hs_sun_azi": {"value": self.hs_sun_azi,
                                                        "description": "Solar azimuth angle (clockwise from North) in "
                                                                       "degrees."},
                                         "hs_sun_el": {"value": self.hs_sun_el,
                                                       "description": "Solar vertical angle (above the horizon) in "
                                                                      "degrees."}
                                     },
                                     "Multiple directions hillshade": {
                                         "mhs_nr_dir": {"value": self.mhs_nr_dir,
                                                        "description": "Number of solar azimuth angles (clockwise "
                                                                       "from North)."},
                                         "mhs_sun_el": {"value": self.mhs_sun_el,
                                                        "description": "Solar vertical angle (above the horizon) in "
                                                                       "degrees."}
                                     },
                                     "Slope gradient": {
                                         "slp_output_units": {"value": self.slp_output_units,
                                                              "description": "Slope output units [radian, degree, "
                                                                             "percent]."}
                                     },
                                     "Simple local relief model": {
                                         "slrm_rad_cell": {"value": self.slrm_rad_cell,
                                                           "description": "Radius for trend assessment in pixels."}
                                     },
                                     "Sky-View Factor": {
                                         "svf_n_dir": {"value": self.svf_n_dir, "description": "Number of directions."},
                                         "svf_r_max": {"value": self.svf_r_max, "description": "Maximal search "
                                                                                               "radious in pixels."},
                                         "svf_noise": {"value": self.svf_noise,
                                                       "description": "The level of noise remove [0-don't remove, "
                                                                      "1-low, 2-med, 3-high]."}
                                     },
                                     "Anisotropic Sky-View Factor": {
                                         "asvf_dir": {"value": self.asvf_dir,
                                                      "description": "Direction of anisotropy in degrees."},
                                         "asvf_level": {"value": self.asvf_level,
                                                        "description": "Level of anisotropy [1-low, 2-high]."}
                                     },
                                     "Sky illumination": {
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
                else:
                    warnings.warn("rvt.default.read_default_from_file: Wrong line '{}'".format(line))
                if parameter_name == "exaggeration_factor":
                    self.ve_factor = float(parameter_value)
                elif parameter_name == "sun_azimuth":
                    self.hs_sun_azi = int(parameter_value)
                elif parameter_name == "sun_elevation":
                    self.hs_sun_el = int(parameter_value)
                elif parameter_name == "hillshade_directions":
                    self.mhs_nr_dir = int(parameter_value)
                elif parameter_name == "sun_elevation":
                    self.mhs_sun_el = int(parameter_value)
                elif parameter_name == "trend_radius":
                    self.slrm_rad_cell = int(parameter_value)
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
                elif parameter_name == "anisotropy_direction":
                    self.asvf_dir = int(parameter_value)
                elif parameter_name == "anisotropy_level":
                    if parameter_value == "low":
                        self.asvf_level = 1
                    elif parameter_value == "high":
                        self.asvf_level = 2
                elif parameter_name == "sky_model":
                    self.sim_sky_mod = parameter_value
                elif parameter_name == "number_points":
                    self.sim_samp_pnts = int(parameter_value)
                elif parameter_name == "max_shadow_dist":
                    self.sim_shadow_dist = int(parameter_value)
                elif parameter_name == "min_radius":
                    self.ld_min_rad = int(parameter_value)
                elif parameter_name == "max_radius":
                    self.ld_max_rad = int(parameter_value)
                # else:
                #     warnings.warn("rvt.default.read_default_from_file: Line '{}' not used.".format(line))
            dat.close()
        elif extension == ".json":
            dat = open(file_path, "r")
            data = json.load(dat)
            default_data = data["default_settings"]
            self.ve_factor = float(default_data["ve_factor"]["value"])
            # Slope gradient
            self.slp_output_units = str(default_data["Slope gradient"]["slp_output_units"]["value"])
            # Hillshade
            self.hs_sun_azi = int(default_data["Hillshade"]["hs_sun_azi"]["value"])
            self.hs_sun_el = int(default_data["Hillshade"]["hs_sun_el"]["value"])
            # Multiple directions hillshade
            self.mhs_nr_dir = int(default_data["Multiple directions hillshade"]["mhs_nr_dir"]["value"])
            self.mhs_sun_el = int(default_data["Multiple directions hillshade"]["mhs_sun_el"]["value"])
            # Simple local relief model
            self.slrm_rad_cell = int(default_data["Simple local relief model"]["slrm_rad_cell"]["value"])
            # Sky-View Factor
            self.svf_n_dir = int(default_data["Sky-View Factor"]["svf_n_dir"]["value"])
            self.svf_r_max = int(default_data["Sky-View Factor"]["svf_r_max"]["value"])
            self.svf_noise = int(default_data["Sky-View Factor"]["svf_noise"]["value"])
            # Anisotropic Sky-View Factor
            self.asvf_dir = int(default_data["Anisotropic Sky-View Factor"]["asvf_dir"]["value"])
            self.asvf_level = int(default_data["Anisotropic Sky-View Factor"]["asvf_level"]["value"])
            # Sky illumination
            self.sim_sky_mod = str(default_data["Sky illumination"]["sim_sky_mod"]["value"])
            self.sim_samp_pnts = int(default_data["Sky illumination"]["sim_samp_pnts"]["value"])
            self.sim_shadow_dist = int(default_data["Sky illumination"]["sim_shadow_dist"]["value"])
            self.sim_shadow_az = int(default_data["Sky illumination"]["sim_shadow_az"]["value"])
            self.sim_shadow_el = int(default_data["Sky illumination"]["sim_shadow_el"]["value"])
            # Local dominance
            self.ld_min_rad = int(default_data["Local dominance"]["ld_min_rad"]["value"])
            self.ld_max_rad = int(default_data["Local dominance"]["ld_max_rad"]["value"])
            self.ld_rad_inc = int(default_data["Local dominance"]["ld_rad_inc"]["value"])
            self.ld_anglr_res = int(default_data["Local dominance"]["ld_anglr_res"]["value"])
            self.ld_observer_h = float(default_data["Local dominance"]["ld_observer_h"]["value"])
            dat.close()

    def get_hillshade_path(self, dem_path):
        dem_path_split = dem_path.split(".")  # split file type
        return "{}_HS_A{}_H{}.tif".format(dem_path_split[0], self.hs_sun_azi, self.hs_sun_el)

    def get_slope_path(self, dem_path):
        dem_path_split = dem_path.split(".")  # split file type
        return "{}_SLOPE.tif".format(dem_path_split[0])

    def get_multi_hillshade_path(self, dem_path):
        dem_path_split = dem_path.split(".")  # split file type
        return "{}_MULTI-HS_D{}_H{}.tif".format(dem_path_split[0], self.mhs_nr_dir, self.mhs_sun_el)

    def get_slrm_path(self, dem_path):
        dem_path_split = dem_path.split(".")  # split file type
        return "{}_SLRM_R{}.tif".format(dem_path_split[0], self.slrm_rad_cell)

    def get_svf_path(self, dem_path):
        dem_path_split = dem_path.split(".")  # split file type
        return "{}_SVF_R{}_D{}.tif".format(dem_path_split[0], self.svf_r_max, self.svf_n_dir)

    def get_asvf_path(self, dem_path):
        dem_path_split = dem_path.split(".")  # split file type
        return "{}_SVF-A_R{}_D{}_A{}.tif".format(dem_path_split[0], self.svf_r_max, self.svf_n_dir, self.asvf_dir)

    def get_opns_path(self, dem_path):
        dem_path_split = dem_path.split(".")  # split file type
        return "{}_OPEN-POS_R{}_D{}.tif".format(dem_path_split[0], self.svf_r_max, self.svf_n_dir)

    def get_neg_opns_path(self, dem_path):
        dem_path_split = dem_path.split(".")  # split file type
        return "{}_OPEN-NEG_R{}_D{}.tif".format(dem_path_split[0], self.svf_r_max, self.svf_n_dir)

    def get_sky_illumination_path(self, dem_path):
        dem_path_split = dem_path.split(".")  # split file type
        return "{}_SIM_{}_{}sp_{}px.tif".format(dem_path_split[0], self.sim_sky_mod, self.sim_samp_pnts,
                                                self.sim_shadow_dist)

    def get_local_dominance_path(self, dem_path):
        dem_path_split = dem_path.split(".")  # split file type
        return "{}_LD_R_M{}-{}_DI{}_A{}_OH{}.tif".format(dem_path_split[0], self.ld_min_rad, self.ld_max_rad,
                                                         self.ld_rad_inc, self.ld_anglr_res, self.ld_observer_h)

    # get slope, aspect dict
    def get_slope(self, dem_arr, resolution_x, resolution_y):
        dict_slp_asp = rvt.vis.slope_aspect(dem=dem_arr, resolution_x=resolution_x, resolution_y=resolution_y,
                                            ve_factor=self.ve_factor, output_units=self.slp_output_units)
        return dict_slp_asp

    # save default slope gradient
    def save_slope(self, dem_path):
        slope_path = self.get_slope_path(dem_path)
        if os.path.isfile(slope_path):
            return 0
        if not os.path.isfile(dem_path):
            raise Exception("rvt.default.DefaultValues.save_slope: dem_path doesn't exist!")
        dem_dataset = rio.open(dem_path)
        t = dem_dataset.transform
        x_res = t[0]
        y_res = -t[4]
        dem_arr = dem_dataset.read()[0]
        dict_slp_asp = self.get_slope(dem_arr=dem_arr, resolution_x=x_res, resolution_y=y_res)
        slope_arr = dict_slp_asp["slope"].astype('float32')
        profile = dem_dataset.profile
        profile.update(dtype='float32')
        slope_dataset = rio.open(slope_path, "w", **profile)
        slope_dataset.write(np.array([slope_arr]))
        dem_dataset.close()
        slope_dataset.close()
        return 1

    # get hillshade array
    def get_hillshade(self, dem_arr, resolution_x, resolution_y):
        hillshade_arr = rvt.vis.hillshade(dem=dem_arr, resolution_x=resolution_x, resolution_y=resolution_y,
                                          sun_azimuth=self.hs_sun_azi, sun_elevation=self.hs_sun_el,
                                          ve_factor=self.ve_factor)
        return hillshade_arr

    # save default hillshade
    def save_hillshade(self, dem_path):
        hillshade_path = self.get_hillshade_path(dem_path)
        if os.path.isfile(hillshade_path):
            return 0
        if not os.path.isfile(dem_path):
            raise Exception("rvt.default.DefaultValues.save_hillshade: dem_path doesn't exist!")
        dem_dataset = rio.open(dem_path)
        t = dem_dataset.transform
        x_res = t[0]
        y_res = -t[4]
        dem_arr = dem_dataset.read()[0]
        hillshade_arr = self.get_hillshade(dem_arr=dem_arr, resolution_x=x_res, resolution_y=y_res).astype('float32')
        profile = dem_dataset.profile
        profile.update(dtype='float32')
        hillshade_dataset = rio.open(hillshade_path, "w", **profile)
        hillshade_dataset.write(np.array([hillshade_arr]))
        dem_dataset.close()
        hillshade_dataset.close()
        return 1

    # get multi hillshade array
    def get_multi_hillshade(self, dem_arr, resolution_x, resolution_y):
        multi_hillshade_arr = rvt.vis.multi_hillshade(dem=dem_arr, resolution_x=resolution_x, resolution_y=resolution_y,
                                                      nr_directions=self.mhs_nr_dir, sun_elevation=self.mhs_sun_el,
                                                      ve_factor=self.ve_factor)
        return multi_hillshade_arr

    def save_multi_hillshade(self, dem_path):
        multi_hillshade_path = self.get_multi_hillshade_path(dem_path)
        if os.path.isfile(multi_hillshade_path):
            return 0
        if not os.path.isfile(dem_path):
            raise Exception("rvt.default.DefaultValues.save_multi_hillshade: dem_path doesn't exist!")
        dem_dataset = rio.open(dem_path)
        t = dem_dataset.transform
        x_res = t[0]
        y_res = -t[4]
        dem_arr = dem_dataset.read()[0]
        multi_hillshade_arr = self.get_multi_hillshade(dem_arr=dem_arr, resolution_x=x_res,
                                                       resolution_y=y_res).astype('float32')
        profile = dem_dataset.profile
        profile.update(dtype='float32')
        profile.update(count=16)
        multi_hillshade_dataset = rio.open(multi_hillshade_path, "w", **profile)
        multi_hillshade_dataset.write(multi_hillshade_arr)
        dem_dataset.close()
        multi_hillshade_dataset.close()
        return 1

    def get_slrm(self, dem_arr):
        slrm_arr = rvt.vis.slrm(dem=dem_arr, radius_cell=self.slrm_rad_cell, ve_factor=self.ve_factor)
        return slrm_arr

    def save_slrm(self, dem_path):
        slrm_path = self.get_slrm_path(dem_path)
        if os.path.isfile(slrm_path):
            return 0
        if not os.path.isfile(dem_path):
            raise Exception("rvt.default.DefaultValues.save_slrm: dem_path doesn't exist!")
        dem_dataset = rio.open(dem_path)
        dem_arr = dem_dataset.read()[0]
        slrm_arr = self.get_slrm(dem_arr=dem_arr).astype('float32')
        profile = dem_dataset.profile
        profile.update(dtype='float32')
        slrm_dataset = rio.open(slrm_path, "w", **profile)
        slrm_dataset.write(np.array([slrm_arr]))
        dem_dataset.close()
        slrm_dataset.close()
        return 1

    # get svf, asvf, opns dict
    def get_sky_view_factor(self, dem_arr, resolution, compute_svf=True, compute_asvf=False, compute_opns=False):
        dict_svf_asvf_opns = rvt.vis.sky_view_factor(dem=dem_arr, resolution=resolution, compute_svf=compute_svf,
                                                     compute_opns=compute_opns, compute_asvf=compute_asvf,
                                                     svf_n_dir=self.svf_n_dir, svf_r_max=self.svf_r_max,
                                                     svf_noise=self.svf_noise, asvf_dir=self.asvf_dir,
                                                     asvf_level=self.asvf_level, ve_factor=self.ve_factor)
        return dict_svf_asvf_opns

    def save_sky_view_factor(self, dem_path, save_svf=True, save_asvf=False, save_opns=False):
        if save_svf:
            svf_path = self.get_svf_path(dem_path)
        if save_asvf:
            asvf_path = self.get_asvf_path(dem_path)
        if save_opns:
            opns_path = self.get_opns_path(dem_path)

        if save_svf and os.path.isfile(svf_path):  # svf already exists, don't need to compute
            save_svf = False
        if save_asvf and os.path.isfile(asvf_path):  # asvf already exists
            save_asvf = False
        if save_opns and os.path.isfile(opns_path):  # opns already exists
            save_opns = False

        if not save_svf and not save_asvf and not save_opns:
            return 0
        if not os.path.isfile(dem_path):
            raise Exception("rvt.default.DefaultValues.save_sky_view_factor: dem_path doesn't exist!")
        dem_dataset = rio.open(dem_path)
        t = dem_dataset.transform
        x_res = t[0]
        y_res = -t[4]
        if x_res != y_res:
            raise Exception("rvt.default.DefaultValues.save_sky_view_factor: dem x resolution not equal y resolution")
        dem_arr = dem_dataset.read()[0]
        dict_svf_asvf_opns = self.get_sky_view_factor(dem_arr=dem_arr, resolution=x_res, compute_svf=save_svf,
                                                      compute_asvf=save_asvf, compute_opns=save_opns)
        if save_svf:
            svf_arr = dict_svf_asvf_opns["svf"]
            svf_arr = svf_arr.astype('float32')
            profile = dem_dataset.profile
            profile.update(dtype='float32')
            svf_dataset = rio.open(svf_path, "w", **profile)
            svf_dataset.write(np.array([svf_arr]))
            svf_dataset.close()
        if save_asvf:
            asvf_arr = dict_svf_asvf_opns["asvf"]
            asvf_arr = asvf_arr.astype('float32')
            profile = dem_dataset.profile
            profile.update(dtype='float32')
            asvf_dataset = rio.open(asvf_path, "w", **profile)
            asvf_dataset.write(np.array([asvf_arr]))
            asvf_dataset.close()
        if save_opns:
            opns_arr = dict_svf_asvf_opns["opns"]
            opns_arr = opns_arr.astype('float32')
            profile = dem_dataset.profile
            profile.update(dtype='float32')
            opns_dataset = rio.open(opns_path, "w", **profile)
            opns_dataset.write(np.array([opns_arr]))
            opns_dataset.close()
        dem_dataset.close()
        return 1

    def get_neg_opns(self, dem_arr, resolution):
        dem_arr = -1 * dem_arr
        dict_neg_opns = rvt.vis.sky_view_factor(dem=dem_arr, resolution=resolution, svf_n_dir=self.svf_n_dir,
                                                svf_r_max=self.svf_r_max, svf_noise=self.svf_noise,
                                                compute_svf=False, compute_asvf=False, compute_opns=True,
                                                ve_factor=self.ve_factor)
        neg_opns_arr = dict_neg_opns["opns"]
        return neg_opns_arr

    def save_neg_opns(self, dem_path):  # negative openness is openness with negative dem (-dem)
        neg_opns_path = self.get_neg_opns_path(dem_path)
        if os.path.isfile(neg_opns_path):
            return 0
        if not os.path.isfile(dem_path):
            raise Exception("rvt.default.DefaultValues.save_neg_opns: dem_path doesn't exist!")
        dem_dataset = rio.open(dem_path)
        t = dem_dataset.transform
        x_res = t[0]
        y_res = -t[4]
        if x_res != y_res:
            raise Exception("rvt.default.DefaultValues.save_sky_view_factor: dem x resolution not equal y resolution")
        dem_arr = dem_dataset.read()[0]
        neg_opns_arr = self.get_neg_opns(dem_arr=dem_arr, resolution=x_res).astype('float32')
        profile = dem_dataset.profile
        profile.update(dtype='float32')
        neg_opns_dataset = rio.open(neg_opns_path, "w", **profile)
        neg_opns_dataset.write(np.array([neg_opns_arr]))
        dem_dataset.close()
        neg_opns_dataset.close()
        return 1

    def get_sky_illumination(self, dem_arr, resolution):
        sky_illumination_arr = rvt.vis.sky_illumination(dem=dem_arr, resolution=resolution, sky_model=self.sim_sky_mod,
                                                        sampling_points=self.sim_samp_pnts,
                                                        shadow_dist=self.sim_shadow_dist,
                                                        shadow_az=self.sim_shadow_az, shadow_el=self.sim_shadow_el,
                                                        ve_factor=self.ve_factor)
        return sky_illumination_arr

    def save_sky_illumination(self, dem_path):
        sky_illumination_path = self.get_sky_illumination_path(dem_path)
        if os.path.isfile(sky_illumination_path):
            return 0
        if not os.path.isfile(dem_path):
            raise Exception("rvt.default.DefaultValues.save_sky_view_factor: dem_path doesn't exist!")
        dem_dataset = rio.open(dem_path)
        t = dem_dataset.transform
        x_res = t[0]
        y_res = -t[4]
        if x_res != y_res:
            raise Exception("rvt.default.DefaultValues.save_sky_view_factor: dem x resolution not equal y resolution")
        dem_arr = dem_dataset.read()[0]
        sky_illumination_arr = self.get_sky_illumination(dem_arr=dem_arr, resolution=x_res).astype('float32')
        profile = dem_dataset.profile
        profile.update(dtype='float32')
        sky_illumination_dataset = rio.open(sky_illumination_path, "w", **profile)
        sky_illumination_dataset.write(np.array([sky_illumination_arr]))
        dem_dataset.close()
        sky_illumination_dataset.close()
        return 1

    def get_local_dominance(self, dem_arr):
        local_dominance_arr = rvt.vis.local_dominance(dem=dem_arr, min_rad=self.ld_min_rad, max_rad=self.ld_max_rad,
                                                      rad_inc=self.ld_rad_inc, angular_res=self.ld_anglr_res,
                                                      observer_height=self.ld_observer_h, ve_factor=self.ve_factor)
        return local_dominance_arr

    def save_local_dominance(self, dem_path):
        local_dominance_path = self.get_local_dominance_path(dem_path)
        if os.path.isfile(local_dominance_path):
            return 0
        if not os.path.isfile(dem_path):
            raise Exception("rvt.default.DefaultValues.save_local_dominance: dem_path doesn't exist!")
        dem_dataset = rio.open(dem_path)
        dem_arr = dem_dataset.read()[0]
        local_dominance_arr = self.get_local_dominance(dem_arr=dem_arr).astype('float32')
        profile = dem_dataset.profile
        profile.update(dtype='float32')
        local_dominance_dataset = rio.open(local_dominance_path, "w", **profile)
        local_dominance_dataset.write(np.array([local_dominance_arr]))
        dem_dataset.close()
        local_dominance_dataset.close()
        return 1

    def save_visualizations(self, dem_path, sav_slope=True, sav_hillshade=True, sav_mulit_hillshade=True, sav_slrm=True,
                            sav_svf=True, sav_asvf=True, sav_opns=True, sav_neg_opns=True, sav_sky_illumination=True,
                            sav_local_dominance=True):
        if sav_slope:
            self.save_slope(dem_path)
        if sav_hillshade:
            self.save_hillshade(dem_path)
        if sav_mulit_hillshade:
            self.save_multi_hillshade(dem_path)
        if sav_slrm:
            self.save_slrm(dem_path)
        if sav_svf or sav_asvf or sav_opns:
            self.save_sky_view_factor(dem_path, save_svf=sav_svf, save_asvf=sav_asvf, save_opns=sav_opns)
        if sav_neg_opns:
            self.save_neg_opns(dem_path)
        if sav_sky_illumination:
            self.save_neg_opns(dem_path)
        if sav_local_dominance:
            self.save_local_dominance(dem_path)
