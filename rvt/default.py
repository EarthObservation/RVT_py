"""
NAME:
    rvt_py, rvt.default - rvt default values for visualization functions

DESCRIPTION:
    Contains all default values for visualisation functions, which can be changed.
    Allows computing from rvt.vis with using defined default values
    and saving output rasters with default names (dependent on default values).

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


class DefaultValues():
    def __init__(self):
        # slope gradient
        self.slp_ve_factor = 1
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

    def read_default_from_file(self, file_path):
        # Example file in dir settings: default_settings.txt
        dat = open(file_path, "r")
        for line in dat:
            line = line.strip()
            if line == "":
                continue
            if line[0] == "#":
                continue
            line_list = line.split("#")  # remove comments
            line_list = line_list[0].split("=")
            if len(line_list) == 2:
                parameter_name = line_list[0].strip()
                parameter_value = line_list[1].strip()
            else:
                warnings.warn("RVT read_default_from_file: Wrong line '{}'".format(line))
            if parameter_name == "hs_sun_azi":
                self.hs_sun_azi = int(parameter_value)
            elif parameter_name == "hs_sun_el":
                self.hs_sun_el = int(parameter_value)
            elif parameter_name == "mhs_nr_dir":
                self.mhs_nr_dir = int(parameter_value)
            elif parameter_name == "mhs_sun_el":
                self.mhs_sun_el = int(parameter_value)
            elif parameter_name == "slrm_rad_cell":
                self.slrm_rad_cell = int(parameter_value)
            elif parameter_name == "svf_n_dir":
                self.svf_n_dir = int(parameter_value)
            elif parameter_name == "svf_r_max":
                self.svf_r_max = int(parameter_value)
            elif parameter_name == "svf_noise":
                self.svf_noise = int(parameter_value)
            elif parameter_name == "asvf_dir":
                self.asvf_dir = int(parameter_value)
            elif parameter_name == "asvf_level":
                self.asvf_level = int(parameter_value)
            elif parameter_name == "sim_sky_mod":
                self.sim_sky_mod = parameter_value
            elif parameter_name == "sim_samp_pnts":
                self.sim_samp_pnts = int(parameter_value)
            elif parameter_name == "sim_shadow_dist":
                self.sim_shadow_dist = int(parameter_value)
            elif parameter_name == "sim_shadow_az":
                self.sim_shadow_az = int(parameter_value)
            elif parameter_name == "sim_shadow_el":
                self.sim_shadow_el = int(parameter_value)
            elif parameter_name == "ld_min_rad":
                self.ld_min_rad = int(parameter_value)
            elif parameter_name == "ld_max_rad":
                self.ld_max_rad = int(parameter_value)
            elif parameter_name == "ld_rad_inc":
                self.ld_rad_inc = int(parameter_value)
            elif parameter_name == "ld_anglr_res":
                self.ld_anglr_res = int(parameter_value)
            elif parameter_name == "ld_observer_h":
                self.ld_observer_h = float(parameter_value)
            else:
                warnings.warn("RVT read_default_from_file: Wrong line '{}'".format(line))
        dat.close()

    def get_hillshade_path(self, dem_path):
        dem_path_split = dem_path.split(".")  # split file type
        return "{}_HS_A{}_H{}.tif".format(dem_path_split[0], self.hs_sun_azi, self.hs_sun_el)

    def get_slope_path(self, dem_path):
        dem_path_split = dem_path.split(".")  # split file type
        return "{}_SLOPE.tif".format(dem_path_split[0])

    # save default slope gradient
    def save_slope(self, dem_path):  # if return 0 - exists, if return 1 - computed and saved
        slope_path = self.get_slope_path(dem_path)
        if os.path.isfile(slope_path):
            return 0
        if not os.path.isfile(dem_path):
            raise Exception("RVT DefaultValues.save_slope: dem_path doesn't exist!")
        dem_dataset = rio.open(dem_path)
        t = dem_dataset.transform
        x_res = t[0]
        y_res = -t[4]
        dem_arr = dem_dataset.read()[0]
        dict_slp_asp = rvt.vis.slope_aspect(dem=dem_arr, resolution_x=x_res,
                                            resolution_y=y_res, ve_factor=self.slp_ve_factor,
                                            output_units=self.slp_output_units)
        slope_arr = dict_slp_asp["slope"]
        slope_arr = slope_arr.astype('float32')
        profile = dem_dataset.profile
        profile.update(dtype='float32')
        slope_dataset = rio.open(slope_path, "w", **profile)
        slope_dataset.write(np.array([slope_arr]))
        dem_dataset.close()
        slope_dataset.close()
        return 1

    # save default hillshade
    def save_hillshade(self, dem_path):  # if return 0 - exists, if return 1 - computed and saved
        hillshade_path = self.get_hillshade_path(dem_path)
        if os.path.isfile(hillshade_path):
            return 0
        if not os.path.isfile(dem_path):
            raise Exception("RVT DefaultValues.save_hillshade: dem_path doesn't exist!")
        dem_dataset = rio.open(dem_path)
        t = dem_dataset.transform
        x_res = t[0]
        y_res = -t[4]
        dem_arr = dem_dataset.read()[0]
        hillshade_arr = rvt.vis.hillshade(dem=dem_arr, resolution_x=x_res, resolution_y=y_res,
                                          sun_azimuth=self.hs_sun_azi, sun_elevation=self.hs_sun_el)
        hillshade_arr = hillshade_arr.astype('float32')
        profile = dem_dataset.profile
        profile.update(dtype='float32')
        hillshade_dataset = rio.open(hillshade_path, "w", **profile)
        hillshade_dataset.write(np.array([hillshade_arr]))
        dem_dataset.close()
        hillshade_dataset.close()
        return 1

    # TODO: make save and get_path methods for all vis functions



