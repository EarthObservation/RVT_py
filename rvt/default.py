"""
rvt.default - rvt default values for visualization functions

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
        """Saves default attributes into file."""
        if file_path is None:
            file_path = r"settings\default_settings.txt"
            if os.path.isfile(file_path):
                pass
            else:
                if not os.path.exists(os.path.dirname(file_path)):
                    os.makedirs(os.path.dirname(file_path))
            dat = open(file_path, "w")
        else:
            dat = open(file_path, "w")
        dat.write("# --------------------------------------------------------------------------------------------------"
                  "--------\n# ----------------------------------------------------------------------------------------"
                  "------------------\n#\n#	Relief Visualization Toolbox python, settings file to change default"
                  " values in rvt.default.DefaultValues class\n#	This class stores default values for visualisation"
                  " functions\n#	method: rvt.default.DefaultValues.read_default_from_file(file_path)\n#\n"
                  "# --------------------------------------------------------------------------------------------------"
                  "--------\n# ----------------------------------------------------------------------------------------"
                  "------------------\n#\n# VISUALIZATION METHODS:\n#\n# Hillshade\n# Multiple directions hillshade\n"
                  "# Slope gradient\n# Simple local relief model\n# Sky-View Factor\n# Anisotropic Sky-View Factor\n"
                  "# Openness - Positive\n# Openness - Negative\n# Sky illumination\n# Local dominance\n\n"
                  "# --------------------------------------------------------------------------------------------------"
                  "--------\n\n")
        dat.write("ve_factor = {}  # vertical exaggeration (positive number!)\n\n".format(self.ve_factor))
        dat.write("# Slope gradient\n")
        dat.write("slp_output_units = {}  # slope output units"
                  " [radian, degree, percent]\n\n".format(self.slp_output_units))
        dat.write("# Hillshade\n")
        dat.write("hs_sun_azi = {}  # solar azimuth angle (clockwise from North) in degrees\n".format(self.hs_sun_azi))
        dat.write("hs_sun_el = {}  # solar vertical angle (above the horizon) in degrees\n\n".format(self.hs_sun_el))
        dat.write("# Multiple directions hillshade\n")
        dat.write("mhs_nr_dir = {}  # number of solar azimuth angles (clockwise from North)\n".format(self.mhs_nr_dir))
        dat.write("mhs_sun_el = {}  # solar vertical angle (above the horizon) in degrees\n\n".format(self.mhs_sun_el))
        dat.write("# Simple local relief model\n")
        dat.write("slrm_rad_cell = {}  # radius for trend assessment in pixels\n\n".format(self.slrm_rad_cell))
        dat.write("# Sky-View Factor\n")
        dat.write("svf_n_dir = {}  # number of directions\n".format(self.svf_n_dir))
        dat.write("svf_r_max = {}  # maximal search radius in pixels\n".format(self.svf_r_max))
        dat.write("svf_noise = {}  # the level of noise remove [0-don't remove, 1-low, 2-med, 3-high]"
                  "\n\n".format(self.svf_noise))
        dat.write("# Anisotropic Sky-View Factor\n")
        dat.write("asvf_dir = {}  # direction of anisotropy in degrees\n".format(self.asvf_dir))
        dat.write("asvf_level = {}  # level of anisotropy [1-low, 2-high]\n\n".format(self.asvf_level))
        dat.write("# Sky illumination\n")
        dat.write("sim_sky_mod = {}  # sky model [overcast, uniform]\n".format(self.sim_sky_mod))
        dat.write("sim_samp_pnts = {}  # number of sampling points [250 or 500]\n".format(self.sim_samp_pnts))
        dat.write("sim_shadow_dist = {}  # max shadow modeling distance in pixels\n".format(self.sim_shadow_dist))
        dat.write("sim_shadow_az = {}  # shadow azimuth in degrees\n".format(self.sim_shadow_az))
        dat.write("sim_shadow_el = {}  # shadow elevation in degrees\n\n".format(self.sim_shadow_el))
        dat.write("# Local dominance\n")
        dat.write("ld_min_rad = {}  # minimum radial distance (in pixels) at which the algorithm starts with visualizat"
                  "ion computation\n".format(self.ld_min_rad))
        dat.write("ld_max_rad = {}  # maximum radial distance (in pixels) at which the algorithm ends with visualizati"
                  "on computation\n".format(self.ld_max_rad))
        dat.write("ld_rad_inc = {}  # radial distance steps in pixels\n".format(self.ld_rad_inc))
        dat.write("ld_anglr_res = {}  # angular step for determination of number of angular directions"
                  "\n".format(self.ld_anglr_res))
        dat.write("ld_observer_h = {}  # height at which we observe the terrain\n".format(self.ld_observer_h))

    def read_default_from_file(self, file_path):
        """Reads default attributes from file."""
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
                warnings.warn("rvt.default.read_default_from_file: Wrong line '{}'".format(line))
            if parameter_name == "ve_factor":
                self.ve_factor = float(parameter_value)
            elif parameter_name == "slp_output_units":
                self.slp_output_units = str(parameter_value)
            elif parameter_name == "hs_sun_azi":
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
                warnings.warn("rvt.default.read_default_from_file: Wrong line '{}'".format(line))
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
