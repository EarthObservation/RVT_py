"""
Relief Visualization Toolbox – Visualization Functions

Contains classes and methods for blending.

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

# TODO: more testing, find and fix bugs if they exists

# python libraries
import numpy as np
import warnings
import rvt.default
import rvt.vis
import rvt.vis_dask
import rvt.blend_func_dask
from rvt.blend_func import *
import os
import json
import datetime
import dask.array as da
from typing import Union, Dict, List, Any, Tuple, Optional
from nptyping import NDArray
import dask


def create_blender_file_example(file_path=None):
    """Create blender .json file example (can be changed and read). Example is VAT - Archaeological combination"""
    data = {"combination": {"name": "VAT - Archaeological",
                            "layers":
                                [
                                    {"layer": "1", "visualization_method": "Sky-View Factor", "norm": "Value",
                                     "min": 0.7, "max": 1.0,
                                     "blend_mode": "Multiply", "opacity": 25},
                                    {"layer": "2", "visualization_method": "Openness - Positive",
                                     "norm": "Value", "min": 68, "max": 93,
                                     "blend_mode": "Overlay", "opacity": 50},
                                    {"layer": "3", "visualization_method": "Slope gradient", "norm": "Value",
                                     "min": 0, "max": 50,
                                     "blend_mode": "Luminosity", "opacity": 50},
                                    {"layer": "4", "visualization_method": "Hillshade", "norm": "Value",
                                     "min": 0, "max": 1,
                                     "blend_mode": "Normal", "opacity": 100},
                                    {"layer": "5", "visualization_method": "None"}
                                ]
                            }}
    if file_path is None:
        file_path = r"settings\blender_file_example.json"
        if os.path.isfile(file_path):
            pass
        else:
            if not os.path.exists(os.path.dirname(file_path)):
                os.makedirs(os.path.dirname(file_path))

    dat = open(file_path, "w")
    dat.write(json.dumps(data, indent=4))
    dat.close()


def _check_min_max(image_chunk: NDArray[np.float32]) -> NDArray[np.float32]:
    """Returns (un)scaled image_chunk. To be mapped over dask.array. Comparision is evaluated only when calling compute."""
    if np.nanmin(image_chunk) < 0 or np.nanmax(image_chunk) > 1:
        new_image_chunk = scale_0_to_1(image_chunk)
        return new_image_chunk
    else:
        return image_chunk

def _check_min_max_rendered(image_chunk: NDArray[np.float32]) -> NDArray[np.float32]:
    """Returns image_chunk (and issues a warning). To be mapped over dask.array. Comparision is evaluated only when calling compute."""
    if np.nanmin(image_chunk) < 0 or np.nanmax(image_chunk) > 1:
        warnings.warn("rvt.blend.BlenderCombination.render_all_images: Rendered image scale distorted")
    return image_chunk


class BlenderLayer:
    """
    Class which define layer for blending. BlenderLayer is basic element in BlenderCombination.layers list.

    Attributes
    ----------
    vis : str
        Visualization method.
    normalization : str
        Normalization type, could be "Value" or "Percent".
    min : float
        Normalization minimum.
    max : float
        Normalization maximum.
    opacity : integer
        Image (visualization) opacity.
    colormap : str
        Colormap form matplotlib (https://matplotlib.org/3.3.2/tutorials/colors/colormaps.html).
    min_colormap_cut : float
        What lower part of colormap to cut to select part of colormap.
        Valid values are between 0 and 1, if 0.2 it cuts off (deletes) 20% of lower colors in colormap.
        If None cut is not applied.
    max_colormap_cut : float
        What upper part of colormap to cut to select part of colormap.
        Valid values are between 0 and 1, if 0.8 it cuts off (deletes) 20% of upper colors in colormap.
        If None cut is not applied.
    image_path : str
        Path to DEM. Doesn't matter if image is not None. Leave None if you would like for blender to compute it.
    image : numpy.array (2D)
        Visualization raster. Leave None if you would like for blender to compute it.
    """

    def __init__(self, vis_method=None, normalization="value", minimum=None, maximum=None,
                 blend_mode="normal", opacity=100, colormap=None, min_colormap_cut=None, max_colormap_cut=None,
                 image=None, image_path=None):
        self.vis = vis_method
        self.normalization = normalization
        self.min = minimum
        self.max = maximum
        self.blend_mode = blend_mode
        self.opacity = opacity
        self.colormap = colormap
        self.min_colormap_cut = min_colormap_cut
        self.max_colormap_cut = max_colormap_cut
        self.image_path = image_path
        self.image = image

    def check_data(self):
        """ Check Attributes """
        if self.normalization == "percent":
            self.normalization = "perc"
        if self.vis is None:  # if visualization is None everything is None
            self.normalization = None
            self.min = None
            self.max = None
            self.blend_mode = None
            self.opacity = None
            self.colormap = None  # leave None if you don't want to apply colormap to grayscale
            self.min_colormap_cut = None  # if None it equals to 0
            self.max_colormap_cut = None  # if None it equals to 1
            self.image = None  # leave None if you wish for blender to compute visualization
            self.image_path = None  # leave None if you wish for blender to compute visualization
        else:
            if not (self.normalization.lower() == "value" or self.normalization.lower() == "perc"):
                raise Exception("rvt.blend.BlenderLayer.check_data: normalization value incorrect!")
            if self.normalization.lower() == "perc" and (self.min + self.max) >= 100:
                raise Exception("rvt.blend.BlenderLayer.check_data: when normalization is perc, min + max has "
                                "to smaller then 100%!")
            if self.min > self.max and self.normalization.lower() == "value":
                raise Exception("rvt.blend.BlenderLayer.check_data: min bigger than max!")
            if self.blend_mode.lower() != "normal" and self.blend_mode.lower() != "multiply" and \
                    self.blend_mode.lower() != "overlay" and self.blend_mode.lower() != "luminosity" and \
                    self.blend_mode.lower() != "screen" and self.blend_mode.lower() != "soft_light":
                raise Exception("rvt.blend.BlenderLayer.check_data: blend_mode incorrect!")
            if 0 > self.opacity > 100:
                raise Exception("rvt.blend.BlenderLayer.check_data: opacity incorrect [0-100]!")
            if self.colormap is None and (self.min_colormap_cut is not None or self.max_colormap_cut is not None):
                raise Exception("rvt.blend.BlenderLayer.check_data: colormap not defined but min_colormap_cut or "
                                "max_colormap_cut are!")
            if self.image is None and self.image_path is None:
                if self.vis.lower() != "slope gradient" and self.vis.lower() != "hillshade" and \
                        self.vis.lower() != "shadow" and self.vis.lower() != "multiple directions hillshade" and \
                        self.vis.lower() != "simple local relief model" and self.vis.lower() != "sky-view factor" and \
                        self.vis.lower() != "anisotropic sky-view factor" and \
                        self.vis.lower() != "openness - positive" and self.vis.lower() != "openness - negative" and \
                        self.vis.lower() != "sky illumination" and self.vis.lower() != "local dominance" and \
                        self.vis.lower() != "multi-scale relief model" and \
                        self.vis.lower() != "multi-scale topographic position":
                    raise Exception("rvt.blend.BlenderLayer.check_data: Incorrect vis, if you don't input image or "
                                    "image_path you have to input known visualization method (vis)!")


class BlenderCombination:
    """
    Class for storing layers (rasters, parameters  for blending) and rendering(blending) into blended raster.

    Attributes
    ----------
    dem_arr : np.array (2D)
        Array with DEM data, needed for calculating visualization functions in memory.
    dem_path : str
        Path to DEM, needed for calculating visualization functions and saving them.
    name : str
        Name of BlenderCombination combination.
    layers : [BlenderLayer]
        List of BlenderLayer instances which will be blended together.
    """

    def __init__(self, dem_arr=None, dem_resolution=None, dem_path=None):
        self.dem_arr = dem_arr
        self.dem_resolution = dem_resolution
        self.dem_path = dem_path
        self.name = ""
        self.layers = []

    def add_dem_arr(self, dem_arr, dem_resolution):
        """Add or change dem_arr attribute and its resolution dem_resolution attribute."""
        self.dem_arr = dem_arr
        self.dem_resolution = dem_resolution

    def add_dem_path(self, dem_path):
        """Add or change dem_path attribute."""
        self.dem_path = dem_path

    def create_layer(self, vis_method=None, normalization="value", minimum=None, maximum=None,
                     blend_mode="normal", opacity=100, colormap=None, min_colormap_cut=None, max_colormap_cut=None,
                     image=None, image_path=None):
        """Create BlenderLayer and adds it to layers attribute."""
        if vis_method is not None:
            layer = BlenderLayer(vis_method=vis_method, normalization=normalization, minimum=minimum, maximum=maximum,
                                 blend_mode=blend_mode, opacity=opacity, colormap=colormap,
                                 min_colormap_cut=min_colormap_cut, max_colormap_cut=max_colormap_cut,
                                 image=image, image_path=image_path)
            self.layers.append(layer)

    def add_layer(self, layer: BlenderLayer):
        """Add BlenderLayer instance to layers attribute."""
        if layer.vis is not None:
            self.layers.append(layer)

    def remove_all_layers(self):
        """Empties layers attribute."""
        self.layers = []

    def layers_info(self):
        layer_nr = 1
        layers_info = []
        for layer in self.layers:
            layers_info.append("layer {}, visualization = {}, normalization = {}, min = {}, max = {}, blend_mode = {},"
                               " opacity = {}, colormap = {},"
                               " min_colormap_cut = {}, max_colormap_cut".format(layer_nr, layer.vis,
                                                                                 layer.normalization,
                                                                                 layer.min,
                                                                                 layer.max, layer.blend_mode,
                                                                                 layer.opacity,
                                                                                 layer.colormap,
                                                                                 layer.min_colormap_cut,
                                                                                 layer.max_colormap_cut))
            layer_nr += 1
        return layers_info

    def read_from_file(self, file_path):
        """Reads class attributes from .json file."""
        # Example file (for file_path) in dir settings: blender_file_example.txt
        dat = open(file_path, "r")
        json_data = json.load(dat)
        self.read_from_json(json_data)
        dat.close()

    def read_from_json(self, json_data):
        """Fill class attributes with json data."""
        self.layers = []
        self.name = json_data["combination"]["name"]
        layers_data = json_data["combination"]["layers"]
        for layer in layers_data:
            layer_name = layer["layer"]
            if layer["visualization_method"] is None:
                continue
            if layer["visualization_method"].lower() == "none" or layer["visualization_method"].lower() == "null":
                continue
            else:
                vis_method = str(layer["visualization_method"])
            norm = str(layer["norm"])
            norm_min = float(layer["min"])
            norm_max = float(layer["max"])
            blend_mode = str(layer["blend_mode"])
            opacity = int(layer["opacity"])
            colormap = None
            min_colormap_cut = None
            max_colormap_cut = None
            try:
                colormap = str(layer["colormap"])
                if colormap.lower() == "null" or colormap.lower() == "none":
                    colormap = None
            except:
                colormap = None
            if colormap is not None:
                try:
                    min_colormap_cut = str(layer["min_colormap_cut"])
                    if min_colormap_cut.lower() == "null" or min_colormap_cut.lower() == "none":
                        min_colormap_cut = None
                except:
                    min_colormap_cut = None
                try:
                    max_colormap_cut = str(layer["max_colormap_cut"])
                    if max_colormap_cut.lower() == "null" or max_colormap_cut.lower() == "none":
                        max_colormap_cut = None
                except:
                    max_colormap_cut = None

            self.add_layer(BlenderLayer(vis_method=vis_method, normalization=norm, minimum=norm_min, maximum=norm_max,
                                        blend_mode=blend_mode, opacity=opacity, colormap=colormap,
                                        min_colormap_cut=min_colormap_cut, max_colormap_cut=max_colormap_cut))

    def save_to_file(self, file_path):
        """Save layers (manually) to .json file. Parameters image and image_path in each layer have to be None,
         visualization has to be correct!"""
        json_data = self.to_json()
        dat = open(file_path, "w")
        dat.write(json.dumps(json_data, indent=4))
        dat.close()

    def to_json(self):
        """Outputs class attributes as json."""
        json_data = {"combination": {"name": self.name,
                                     "layers": []
                                     }}
        i_layer = 1
        for layer in self.layers:
            if layer.colormap is None:
                json_data["combination"]["layers"].append({"layer": str(i_layer), "visualization_method": layer.vis,
                                                           "norm": layer.normalization, "min": layer.min,
                                                           "max": layer.max, "blend_mode": layer.blend_mode,
                                                           "opacity": layer.opacity})
            else:
                if layer.min_colormap_cut is None and layer.max_colormap_cut is None:
                    json_data["combination"]["layers"].append({"layer": str(i_layer), "visualization_method": layer.vis,
                                                               "norm": layer.normalization, "min": layer.min,
                                                               "max": layer.max, "blend_mode": layer.blend_mode,
                                                               "opacity": layer.opacity, "colormap": layer.colormap})
                elif layer.min_colormap_cut is not None and layer.max_colormap_cut is None:
                    json_data["combination"]["layers"].append({"layer": str(i_layer), "visualization_method": layer.vis,
                                                               "norm": layer.normalization, "min": layer.min,
                                                               "max": layer.max, "blend_mode": layer.blend_mode,
                                                               "opacity": layer.opacity, "colormap": layer.colormap,
                                                               "min_colormap_cut": layer.min_colormap_cut})
                elif layer.min_colormap_cut is None and layer.max_colormap_cut is not None:
                    json_data["combination"]["layers"].append({"layer": str(i_layer), "visualization_method": layer.vis,
                                                               "norm": layer.normalization, "min": layer.min,
                                                               "max": layer.max, "blend_mode": layer.blend_mode,
                                                               "opacity": layer.opacity, "colormap": layer.colormap,
                                                               "max_colormap_cut": layer.max_colormap_cut})
                elif layer.min_colormap_cut is not None and layer.max_colormap_cut is not None:
                    json_data["combination"]["layers"].append({"layer": str(i_layer), "visualization_method": layer.vis,
                                                               "norm": layer.normalization, "min": layer.min,
                                                               "max": layer.max, "blend_mode": layer.blend_mode,
                                                               "opacity": layer.opacity, "colormap": layer.colormap,
                                                               "min_colormap_cut": layer.min_colormap_cut,
                                                               "max_colormap_cut": layer.max_colormap_cut})
            i_layer += 1
        return json_data

    def check_data(self):
        for layer in self.layers:
            layer.check_data()

    def render_all_images(self, default=None, save_visualizations=False, save_render_path=None, save_float=True,
                          save_8bit=False, no_data=None):
        """Render all layers and returns blended image. If specific layer (BlenderLayer) in layers has image
        (is not None), method uses this image, if image is None and layer has image_path method reads image from
        path. If both image and image_path are None method calculates visualization. If save_visualization is True
        method needs dem_path and saves each visualization (if it doesn't exists) in directory of dem_path,
        else (save_visualization=False) method needs dem_arr, dem_resolution and calculates each visualization
        simultaneously (in memory). Be careful save_visualisation applies only if specific BlenderLayer
        image and image_path are None. Parameter no_data changes all pixels with this values to np.nan,
        if save_visualizations is Ture it is not needed."""

        # check data
        self.check_data()

        if save_render_path is not None and self.dem_path is None:
            raise Exception(
                "rvt.blend.BlenderCombination.render_all_images: If you would like to save rendered image (blender), "
                "you have to define dem_path (BlenderCombination.add_dem_path())!")
        if not save_float and not save_8bit and save_render_path:
            raise Exception(
                "rvt.blend.BlenderCombination.render_all_images: If you would like to save rendered image (blender), "
                "you have to set save_float or save_8bit to True!")
        if save_render_path is not None:
            save_render_directory = os.path.abspath(os.path.dirname(save_render_path))
            save_render_8bit_path = r"{}\{}_8bit.tif".format(save_render_directory,
                                                             os.path.splitext(os.path.basename(save_render_path))[0])
        else:
            save_render_directory = None

        # default is rvt.default.DefaultValues class
        if default is None:  # if not defined, predefined values are used
            default = rvt.default.DefaultValues()

        # rendering across all layers - form last to first layer
        rendered_image = []
        for i_img in range(len(self.layers) - 1, -1, -1):
            visualization = self.layers[i_img].vis
            if visualization is None:  # empty layer, skip
                continue

            min_norm = self.layers[i_img].min
            max_norm = self.layers[i_img].max
            normalization = self.layers[i_img].normalization
            image = self.layers[i_img].image
            image_path = self.layers[i_img].image_path

            if save_visualizations and self.dem_path is None and image_path is None and image is None:
                raise Exception(
                    "rvt.blend.BlenderCombination.render_all_images: If you would like to save visualizations, "
                    "you have to define dem_path (BlenderCombination.add_dem_path())!")
            if not save_visualizations and self.dem_arr is None and self.dem_resolution is None and \
                    image_path is None and image is None:
                raise Exception(
                    "rvt.blend.BlenderCombination.render_all_images: If you would like to compute visualizations, "
                    "you have to define dem_arr and its resolution (BlenderCombination.add_dem_arr())!")

            # normalize images
            # if image is not presented and image_path is
            norm_image = None
            if image is None and image_path is not None:
                norm_image = rvt.blend_func_dask.dask_normalize_image(image = rvt.default.get_raster_dask_arr(image_path)["array"], 
                                                  visualization = visualization, min_norm = min_norm, max_norm = max_norm, normalization = normalization)
                                             
            # if image is presented
            elif image is not None:
                norm_image = rvt.blend_func_dask.dask_normalize_image(image = image, visualization = visualization, 
                                                  min_norm = min_norm, max_norm = max_norm, normalization = normalization)
            # they are both none
            else:
                if self.layers[i_img].vis.lower() == "slope gradient":
                    if save_visualizations:
                        default.save_dask_slope(dem_path=self.dem_path, custom_dir=save_render_directory, save_float=True,
                                           save_8bit=False)
                        image_path = default.get_slope_path(self.dem_path)
                        norm_image = rvt.blend_func_dask.dask_normalize_image(image = rvt.default.get_raster_dask_arr(image_path)["array"],
                                                                         visualization =  visualization, min_norm = min_norm, max_norm = max_norm, normalization = normalization)
                    else:
                        image = default.get_dask_slope(dem_arr=self.dem_arr, resolution_x=self.dem_resolution,
                                                        resolution_y=self.dem_resolution, no_data=no_data)
                        norm_image = rvt.blend_func_dask.dask_normalize_image(image= image, visualization = visualization,
                                                                min_norm = min_norm, max_norm = max_norm, normalization = normalization)
                elif self.layers[i_img].vis.lower() == "hillshade":
                    if save_visualizations:
                        default.save_dask_hillshade(dem_path=self.dem_path, custom_dir=save_render_directory,
                                               save_float=True, save_8bit=False)
                        image_path = default.get_hillshade_path(self.dem_path)
                        norm_image = rvt.blend_func_dask.dask_normalize_image(image = rvt.default.get_raster_dask_arr(image_path)["array"],
                                                                         visualization =  visualization, min_norm = min_norm, max_norm = max_norm, normalization = normalization)
                    else:
                        image = default.get_dask_hillshade(dem_arr=self.dem_arr, resolution_x=self.dem_resolution,
                                                      resolution_y=self.dem_resolution, no_data=no_data)
                        norm_image = rvt.blend_func_dask.dask_normalize_image(image= image, visualization = visualization, min_norm = 
                                                                            min_norm, max_norm = max_norm, normalization = normalization)
                elif self.layers[i_img].vis.lower() == "shadow":
                    raise Exception("Not implemented with dask yet!")
                    if save_visualizations:
                        default.save_hillshade(dem_path=self.dem_path, custom_dir=save_render_directory,
                                               save_float=True, save_8bit=False, save_shadow=True)
                        image_path = default.get_shadow_path(self.dem_path)
                        norm_image = normalize_image(visualization, rvt.default.get_raster_arr(image_path)["array"],
                                                     min_norm, max_norm, normalization)
                    else:
                        image = default.get_shadow(dem_arr=self.dem_arr, resolution=self.dem_resolution,
                                                   no_data=no_data)
                        norm_image = normalize_image(visualization, image, min_norm, max_norm, normalization)

                elif self.layers[i_img].vis.lower() == "multiple directions hillshade":
                    warnings.warn("Multiple directions hillshade may contain bugs.")
                    if save_visualizations:
                        default.save_dask_multi_hillshade(dem_path=self.dem_path, custom_dir=save_render_directory,
                                                     save_float = False, save_8bit=True) 
                                                     ## save_float (vis multihillshade): nr_dir bands, save_8bit (3 times vis hillshade):
                                                     ## 3 bands
                        image_path = default.get_multi_hillshade_path(self.dem_path, bit8=True)
                        ## This image is already normalized and saved in rvt.default. Why is it normalized again here?
                        ## This norm_image (save_visualizations = True) is also not equal to the one below (save_visualizations = False). 
                        ## Normalization in rvt.default is differnt no normalization here...
                        norm_image = rvt.blend_func_dask.dask_normalize_image(image = rvt.default.get_raster_dask_arr(image_path)["array"],
                                                                         visualization =  visualization, min_norm = min_norm, max_norm = max_norm, normalization = normalization)
                    else:
                        red_band_arr = rvt.vis_dask.dask_hillshade(input_dem=self.dem_arr, resolution_x=self.dem_resolution,
                                                         resolution_y=self.dem_resolution,
                                                         sun_elevation=default.mhs_sun_el, sun_azimuth=315,
                                                         no_data=no_data)
                        green_band_arr = rvt.vis_dask.dask_hillshade(input_dem=self.dem_arr, resolution_x=self.dem_resolution,
                                                           resolution_y=self.dem_resolution,
                                                           sun_elevation=default.mhs_sun_el, sun_azimuth=22.5,
                                                           no_data=no_data)
                        blue_band_arr = rvt.vis_dask.dask_hillshade(input_dem=self.dem_arr, resolution_x=self.dem_resolution,
                                                          resolution_y=self.dem_resolution,
                                                          sun_elevation=default.mhs_sun_el, sun_azimuth=90,
                                                          no_data=no_data)
                        image = da.stack((red_band_arr, green_band_arr, blue_band_arr)).rechunk({0: -1}) #check performance of rechunking
                        image, = dask.optimize(image)
                        norm_image = rvt.blend_func_dask.dask_normalize_image(image= image, visualization = visualization, min_norm = 
                                                                            min_norm, max_norm = max_norm, normalization = normalization)
                elif self.layers[i_img].vis.lower() == "simple local relief model":
                    if save_visualizations:
                        default.save_dask_slrm(dem_path=self.dem_path, custom_dir=save_render_directory, save_float=True,
                                          save_8bit=False)
                        image_path = default.get_slrm_path(self.dem_path)
                        norm_image = rvt.blend_func_dask.dask_normalize_image(visualization, rvt.default.get_raster_dask_arr(image_path)["array"], min_norm, max_norm, normalization)
                    else:
                        image = default.get_dask_slrm(dem_arr=self.dem_arr, no_data=no_data)
                        norm_image = rvt.blend_func_dask.dask_normalize_image(visualization, image, min_norm, max_norm, normalization)
                elif self.layers[i_img].vis.lower() == "sky-view factor":
                    if save_visualizations:
                        default.save_dask_sky_view_factor(dem_path=self.dem_path, save_svf=True, save_asvf=False,
                                                     save_opns=False, custom_dir=save_render_directory, save_float=True,
                                                     save_8bit=False)
                        image_path = default.get_svf_path(self.dem_path)
                        norm_image = rvt.blend_func_dask.dask_normalize_image(image = rvt.default.get_raster_dask_arr(image_path)["array"],
                                                                         visualization =  visualization, min_norm = min_norm, max_norm = max_norm, normalization = normalization)
                    else:
                        image = default.get_dask_sky_view_factor(dem_arr=self.dem_arr, resolution=self.dem_resolution,
                                                            compute_svf=True, compute_asvf=False,
                                                            compute_opns=False, no_data=no_data)
                        norm_image = rvt.blend_func_dask.dask_normalize_image(image= image, visualization = visualization, min_norm = 
                                                                            min_norm, max_norm = max_norm, normalization = normalization)
                elif self.layers[i_img].vis.lower() == "anisotropic sky-view factor":
                    if save_visualizations:
                        default.save_dask_sky_view_factor(dem_path=self.dem_path, save_svf=False, save_asvf=True,
                                                     save_opns=False, custom_dir=save_render_directory, save_float=True,
                                                     save_8bit=False)
                        image_path = default.get_asvf_path(self.dem_path)
                        norm_image = rvt.blend_func_dask.dask_normalize_image(visualization, rvt.default.get_raster_arr(image_path)["array"], min_norm, max_norm, normalization)
                    else:
                        image = default.get_dask_sky_view_factor(dem_arr=self.dem_arr, resolution=self.dem_resolution,
                                                            compute_svf=False, compute_asvf=True,
                                                            compute_opns=False, no_data=no_data)["asvf"]
                        norm_image = rvt.blend_func_dask.dask_normalize_image(visualization, image, min_norm, max_norm, normalization)
                elif self.layers[i_img].vis.lower() == "openness - positive":
                    if save_visualizations:
                        default.save_dask_sky_view_factor(dem_path=self.dem_path, save_svf=False, save_asvf=False,
                                                     save_opns=True, custom_dir=save_render_directory, save_float=True,
                                                     save_8bit=False)
                        image_path = default.get_opns_path(self.dem_path)
                        norm_image = rvt.blend_func_dask.dask_normalize_image(image = rvt.default.get_raster_dask_arr(image_path)["array"],
                                                                         visualization =  visualization, min_norm = min_norm, max_norm = max_norm, normalization = normalization)
                    else:
                        image = default.get_dask_sky_view_factor(dem_arr=self.dem_arr, resolution=self.dem_resolution,
                                                            compute_svf=False, compute_asvf=False,
                                                            compute_opns=True, no_data=no_data)
                        norm_image = rvt.blend_func_dask.dask_normalize_image(image= image, visualization = visualization, min_norm = 
                                                                            min_norm, max_norm = max_norm, normalization = normalization)
                elif self.layers[i_img].vis.lower() == "openness - negative":
                    if save_visualizations:
                        default.save_dask_neg_opns(dem_path=self.dem_path, custom_dir=save_render_directory, save_float=True,
                                              save_8bit=False)
                        image_path = default.get_neg_opns_path(self.dem_path)
                        norm_image = rvt.blend_func_dask.dask_normalize_image(visualization, rvt.default.get_raster_dask_arr(image_path)["array"],
                                                     min_norm, max_norm, normalization)
                    else:
                        image = default.get_dask_neg_opns(dem_arr=self.dem_arr, resolution=self.dem_resolution,
                                                     no_data=no_data)
                        norm_image = rvt.blend_func_dask.dask_normalize_image(visualization, image, min_norm, max_norm, normalization)
                elif self.layers[i_img].vis.lower() == "sky illumination":
                    raise Exception("Not implemented with dask yet!")
                    if save_visualizations:
                        default.save_sky_illumination(dem_path=self.dem_path, custom_dir=save_render_directory,
                                                      save_float=True, save_8bit=False)
                        image_path = default.get_sky_illumination_path(self.dem_path)
                        norm_image = normalize_image(visualization, rvt.default.get_raster_arr(image_path)["array"],
                                                     min_norm, max_norm, normalization)
                    else:
                        image = default.get_sky_illumination(dem_arr=self.dem_arr, resolution=self.dem_resolution,
                                                             no_data=no_data)
                        norm_image = normalize_image(visualization, image, min_norm, max_norm, normalization)
                elif self.layers[i_img].vis.lower() == "local dominance":
                    if save_visualizations:
                        default.save_dask_local_dominance(dem_path=self.dem_path, custom_dir=save_render_directory,
                                                     save_float=True, save_8bit=False)
                        image_path = default.get_local_dominance_path(self.dem_path)
                        norm_image = rvt.blend_func_dask.dask_normalize_image(visualization, rvt.default.get_raster_dask_arr(image_path)["array"], min_norm, max_norm, normalization)
                    else:
                        image = default.get_dask_local_dominance(dem_arr=self.dem_arr, no_data=no_data)
                        norm_image = rvt.blend_func_dask.dask_normalize_image(visualization, image, min_norm, max_norm, normalization)
                elif self.layers[i_img].vis.lower() == "multi-scale relief model":
                    if save_visualizations:
                        default.save_dask_msrm(dem_path=self.dem_path, custom_dir=save_render_directory,
                                          save_float=True, save_8bit=False)
                        image_path = default.get_msrm_path(self.dem_path)
                        norm_image = rvt.blend_func_dask.dask_normalize_image(visualization, rvt.default.get_raster_dask_arr(image_path)["array"], min_norm, max_norm, normalization)
                    else:
                        image = default.get_dask_msrm(dem_arr=self.dem_arr, resolution=self.dem_resolution, no_data=no_data)
                        norm_image = rvt.blend_func_dask.dask_normalize_image(visualization, image, min_norm, max_norm, normalization)
                elif self.layers[i_img].vis.lower() == "multi-scale topographic position":
                    if save_visualizations:
                        default.save_dask_mstp(dem_path=self.dem_path, custom_dir=save_render_directory)
                        image_path = default.get_mstp_path(self.dem_path)
                        norm_image = rvt.blend_func_dask.dask_normalize_image(visualization, rvt.default.get_raster_dask_arr(image_path)["array"], min_norm, max_norm, normalization)
                    else:
                        image = default.get_mstp(dem_arr=self.dem_arr, no_data=no_data)
                        norm_image = rvt.blend_func_dask.dask_normalize_image(visualization, image, min_norm, max_norm, normalization)

            colormap = self.layers[i_img].colormap
            min_colormap_cut = self.layers[i_img].min_colormap_cut
            max_colormap_cut = self.layers[i_img].max_colormap_cut
            if colormap is not None and len(image.shape) < 3:
                norm_image = rvt.blend_func_dask.dask_gray_scale_to_color_ramp(gray_scale=norm_image, colormap=colormap, output_8bit=False,
                                                      min_colormap_cut=min_colormap_cut,
                                                      max_colormap_cut=max_colormap_cut, alpha=False)

            # if current layer has visualization applied, but there has been no rendering
            # of images yet, than current layer will be the initial value of rendered_image
            if type(rendered_image) != da.Array: ### TODO: better check for if rendered_image == [] ?
                rendered_image = norm_image
                continue
            else:
                active = norm_image
                background = rendered_image
                blend_mode = self.layers[i_img].blend_mode
                opacity = self.layers[i_img].opacity

                active = da.map_blocks(_check_min_max, active, dtype = np.float32)
                background = da.map_blocks(_check_min_max, background, dtype = np.float32)
                top = rvt.blend_func_dask.dask_blend_images(active = active, background = background, blend_mode = blend_mode)
                rendered_image = rvt.blend_func_dask.dask_render_images(active = top, background = background, opacity = opacity)
                rendered_image = da.map_blocks(_check_min_max_rendered, rendered_image, dtype = np.float32)

        if save_render_path is not None:  # if paths presented it saves image
            if save_float:
                rvt.default.dask_save_raster_tif(src_raster_path=self.dem_path, out_raster_path=save_render_path,
                                        out_raster_arr=rendered_image, dtype_to_save= "float32")
            if save_8bit:
                rendered_image_8bit = rvt.vis_dask.dask_byte_scale(data = rendered_image, c_min = 0, c_max = 1)
                rvt.default.dask_save_raster_tif(src_raster_path=self.dem_path, out_raster_arr = rendered_image_8bit, 
                                                out_raster_path = save_render_8bit_path, dtype_to_save="uint8")
        return rendered_image  # returns float

    def create_log_file(self, dem_path, combination_name, render_path, default: rvt.default.DefaultValues,
                        terrain_sett_name=None, custom_dir=None, computation_time=None):
        """Creates log file in custom_dir, if custom_dir=None it creates it in dem directory (dem_path)."""
        dict_arr_res = rvt.default.get_raster_dask_arr(raster_path=dem_path)
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
        log_name = "{}_blend_log_{}.txt".format(dem_name, log_file_time_str)
        log_path = os.path.join(log_dir, log_name)
        dat = open(log_path, "w")
        dat.write(
            "===============================================================================================\n"
            "Relief Visualization Toolbox (python), blender log\n"
            "Copyright:\n"
            "\tResearch Centre of the Slovenian Academy of Sciences and Arts\n"
            "\tUniversity of Ljubljana, Faculty of Civil and Geodetic Engineering\n"
            "===============================================================================================\n")
        dat.write("\n\n\n")

        dat.write("Processing info about visualizations\n"
                  "===============================================================================================\n\n")
        dat.write("# Metadata of the input file\n\n")
        dat.write("\tInput filename:\t\t{}\n".format(dem_path))
        dat.write("\tNumber of rows:\t\t{}\n".format(nr_rows))
        dat.write("\tNumber of columns:\t{}\n".format(nr_cols))
        dat.write("\tNumber of bands:\t{}\n".format(nr_bands))
        dat.write("\tResolution (x, y):\t{}, {}\n".format(resolution[0], resolution[1]))
        dat.write("\n")

        dat.write("# Selected visualization parameters\n")
        dat.write("\tOverwrite: {}\n".format(default.overwrite))
        dat.write("\tVertical exaggeration factor: {}\n".format(default.ve_factor))
        dat.write("\n")

        dat.write("# Combination:\n\n")
        dat.write("Combination name: {}\n".format(combination_name))
        if terrain_sett_name is not None:
            dat.write("Terrain settings: {}\n".format(terrain_sett_name))
        dat.write("\t>> Output render file:\n")
        dat.write("\t\t{}\n\n".format(render_path))
        dat.write("=== LAYERS ===\n\n")
        i_layer = 1
        for layer in self.layers:
            dat.write("Layer: {}\n".format(i_layer))
            dat.write("Visualization: {}\n".format(layer.vis))
            if layer.vis.lower() == "hillshade":
                dat.write("\ths_sun_el=\t\t{}\n".format(default.hs_sun_el))
                dat.write("\ths_sun_azi=\t\t{}\n".format(default.hs_sun_azi))
            elif layer.vis.lower() == "shadow":
                dat.write("\ths_sun_el=\t\t{}\n".format(default.hs_sun_el))
                dat.write("\ths_sun_azi=\t\t{}\n".format(default.hs_sun_azi))
            elif layer.vis.lower() == "multiple directions hillshade":
                dat.write("\tmhs_sun_el=\t\t{}\n".format(default.mhs_sun_el))
                dat.write("\tmhs_nr_dir=\t\t{}\n".format(default.mhs_nr_dir))
            elif layer.vis.lower() == "slope gradient":
                dat.write("\tslp_output_units=\t\t{}\n".format(default.slp_output_units))
            elif layer.vis.lower() == "simple local relief model":
                dat.write("\tslrm_rad_cell=\t\t{}\n".format(default.slrm_rad_cell))
            elif layer.vis.lower() == "sky-view factor":
                dat.write("\tsvf_n_dir=\t\t{}\n".format(default.svf_n_dir))
                dat.write("\tsvf_noise=\t\t{}\n".format(default.svf_noise))
                dat.write("\tsvf_r_max=\t\t{}\n".format(default.svf_r_max))
            elif layer.vis.lower() == "anisotropic sky-view factor":
                dat.write("\tsvf_n_dir=\t\t{}\n".format(default.svf_n_dir))
                dat.write("\tsvf_noise=\t\t{}\n".format(default.svf_noise))
                dat.write("\tsvf_r_max=\t\t{}\n".format(default.svf_r_max))
                dat.write("\tasvf_level=\t\t{}\n".format(default.asvf_level))
                dat.write("\tasvf_dir=\t\t{}\n".format(default.asvf_dir))
            elif layer.vis.lower() == "openness - positive":
                dat.write("\tsvf_n_dir=\t\t{}\n".format(default.svf_n_dir))
                dat.write("\tsvf_noise=\t\t{}\n".format(default.svf_noise))
                dat.write("\tsvf_r_max=\t\t{}\n".format(default.svf_r_max))
            elif layer.vis.lower() == "openness - negative":
                dat.write("\tsvf_n_dir=\t\t{}\n".format(default.svf_n_dir))
                dat.write("\tsvf_noise=\t\t{}\n".format(default.svf_noise))
                dat.write("\tsvf_r_max=\t\t{}\n".format(default.svf_r_max))
            elif layer.vis.lower() == "sky illumination":
                dat.write("\tsim_sky_mod=\t\t{}\n".format(default.sim_sky_mod))
                dat.write("\tsim_compute_shadow=\t\t{}\n".format(default.sim_compute_shadow))
                dat.write("\tsim_shadow_az=\t\t{}\n".format(default.sim_shadow_az))
                dat.write("\tsim_shadow_el=\t\t{}\n".format(default.sim_shadow_el))
                dat.write("\tsim_nr_dir=\t\t{}\n".format(default.sim_nr_dir))
                dat.write("\tsim_shadow_dist=\t\t{}\n".format(default.sim_shadow_dist))
            elif layer.vis.lower() == "local dominance":
                dat.write("\t\tld_rad_inc=\t\t{}\n".format(default.ld_rad_inc))
                dat.write("\t\tld_min_rad=\t\t{}\n".format(default.ld_min_rad))
                dat.write("\t\tld_max_rad=\t\t{}\n".format(default.ld_max_rad))
                dat.write("\t\tld_anglr_res=\t\t{}\n".format(default.ld_anglr_res))
                dat.write("\t\tld_observer_h=\t\t{}\n".format(default.ld_observer_h))
            elif layer.vis.lower() == "multi-scale relief model":
                dat.write("\t\tmsrm_feature_min=\t\t{}\n".format(default.msrm_feature_min))
                dat.write("\t\tmsrm_feature_max=\t\t{}\n".format(default.msrm_feature_max))
                dat.write("\t\tmsrm_scaling_factor=\t\t{}\n".format(default.msrm_scaling_factor))
            elif layer.vis.lower() == "multi-scale topographic position":
                dat.write("\t\tmstp_local_scale=\t\t({}, {}, {})\n".format(
                    default.mstp_local_scale[0], default.mstp_local_scale[1], default.mstp_local_scale[2]))
                dat.write("\t\tmstp_meso_scale=\t\t({}, {}, {})\n".format(
                    default.mstp_meso_scale[0], default.mstp_meso_scale[1], default.mstp_meso_scale[2]))
                dat.write("\t\tmstp_broad_scale=\t\t({}, {}, {})\n".format(
                    default.mstp_broad_scale[0], default.mstp_broad_scale[1], default.mstp_broad_scale[2]))
                dat.write("\t\tmstp_lightness=\t\t{}\n".format(default.mstp_lightness))
            dat.write("Norm: {}\n".format(layer.normalization))
            dat.write("Linear normalization, min: {}, max: {}\n".format(layer.min, layer.max))
            dat.write("Opacity: {}\n".format(layer.opacity))
            if layer.colormap is not None:
                dat.write("Colormap: {}\n".format(layer.colormap))
                if layer.min_colormap_cut is not None:
                    dat.write("Minimum Colormap cut: {}\n".format(layer.min_colormap_cut))
                if layer.max_colormap_cut is not None:
                    dat.write("Maximum Colormap cut: {}\n".format(layer.max_colormap_cut))
            dat.write("\n")

            i_layer += 1

        if computation_time is not None:
            dat.write("# Computation time: {}".format(computation_time))
        dat.close()


def compare_2_combinations(combination1: BlenderCombination, combination2: BlenderCombination):
    if len(combination1.layers) != len(combination2.layers):
        return False
    for i_layer in range(len(combination1.layers)):
        if combination1.layers[i_layer].vis.lower() != combination2.layers[i_layer].vis.lower():
            return False
        # if combination1.layers[i_layer].normalization.lower() != combination2.layers[i_layer].normalization.lower():
        #     return False
        # if combination1.layers[i_layer].min != combination2.layers[i_layer].min:
        #     return False
        # if combination1.layers[i_layer].max != combination2.layers[i_layer].max:
        #     return False
        if combination1.layers[i_layer].blend_mode.lower() != combination2.layers[i_layer].blend_mode.lower():
            return False
        if combination1.layers[i_layer].opacity != combination2.layers[i_layer].opacity:
            return False
    return True


class BlenderCombinations:
    """
    Class for storing combinations.

    Attributes
    ----------
    combinations : [BlenderCombination]
        List of BlenderCombination instances.
    """

    def __init__(self):
        self.combinations = []  # list of BlenderCombination

    def add_combination(self, combination: BlenderCombination, name=None):
        """Adds cobmination if parameter name not None it renames combination."""
        if name is not None:
            combination.name = name
        self.combinations.append(combination)

    def remove_all_combinations(self):
        """Removes all combinations from self.combinations."""
        self.combinations = []

    def select_combination_by_name(self, name):
        """Select first combination where self.combinations.BlenderCombination.name = name."""
        for combination in self.combinations:
            if combination.name == name:
                return combination

    def remove_combination_by_name(self, name):
        """Removes all combinations where self.combinations.BlenderCombination.name = name.
        If combinations list is empty function returns 0, else 1."""
        new_combinations = []
        for combination in self.combinations:
            if combination.name != name:
                new_combinations.append(combination)
        self.combinations = new_combinations
        if not new_combinations:
            return 0
        else:
            return 1

    def read_from_file(self, file_path):
        """Reads combinations from .json file."""
        self.combinations = []
        dat = open(file_path, "r")
        json_data = json.load(dat)
        combinations_data = json_data["combinations"]
        for combination_data in combinations_data:
            combination = BlenderCombination()
            combination.read_from_json(combination_data)
            self.combinations.append(combination)
        dat.close()

    def save_to_file(self, file_path):
        """Saves combination to .json file."""
        json_data = {"combinations": []}
        for combination in self.combinations:
            json_data["combinations"].append(combination.to_json())
        dat = open(file_path, "w")
        dat.write(json.dumps(json_data, indent=4))
        dat.close()

    def combination_in_combinations(self, input_combination: BlenderCombination):
        """If input_combination (BlenderCombination) has same attributes as one of the combinations (self), method
         returns name of the combination (from combinations). If there is no equal one it returns None."""
        for combination in self.combinations:
            if compare_2_combinations(input_combination, combination):
                return combination.name
        return None

    def combinations_names(self):
        """Returns list of combinations names."""
        names_list = []
        for combination in self.combinations:
            names_list.append(combination.name)
        return names_list


class TerrainSettings:
    """Terrain settings for GUI."""

    def __init__(self):
        self.name = None
        # slope gradient
        self.slp_output_units = None
        # hillshade
        self.hs_sun_azi = None
        self.hs_sun_el = None
        # multi hillshade
        self.mhs_nr_dir = None
        self.mhs_sun_el = None
        # simple local relief model
        self.slrm_rad_cell = None
        # sky view factor
        self.svf_n_dir = None
        self.svf_r_max = None
        self.svf_noise = None
        # anisotropic sky-view factor
        self.asvf_dir = None
        self.asvf_level = None
        # positive openness
        # negative openness
        # sky_illum
        self.sim_sky_mod = None
        self.sim_compute_shadow = None
        self.sim_nr_dir = None
        self.sim_shadow_dist = None
        self.sim_shadow_az = None
        self.sim_shadow_el = None
        # local dominance
        self.ld_min_rad = None
        self.ld_max_rad = None
        self.ld_rad_inc = None
        self.ld_anglr_res = None
        self.ld_observer_h = None
        # multi-scale relief model
        self.msrm_feature_min = None
        self.msrm_feature_max = None
        self.msrm_scaling_factor = None
        # multi-scale topographic position
        self.mstp_local_scale = None
        self.mstp_meso_scale = None
        self.mstp_broad_scale = None
        self.mstp_lightness = None

        # linear histogram stretches tuple(min, max)
        self.hs_stretch = None
        self.mhs_stretch = None
        self.slp_stretch = None
        self.slrm_stretch = None
        self.svf_stretch = None
        self.asvf_stretch = None
        self.pos_opns_stretch = None
        self.neg_opns_stretch = None
        self.sim_stretch = None
        self.ld_stretch = None
        self.msrm_stretch = None
        self.mstp_stretch = None

    def read_from_file(self, file_path):
        """Reads combinations from .json file."""
        dat = open(file_path, "r")
        json_data = json.load(dat)
        self.read_from_json(json_data)
        dat.close()

    def read_from_json(self, json_data):
        """Reads json dict and fills self attributes."""
        self.__init__()
        terrain_data = json_data["terrain_settings"]
        self.name = terrain_data["name"]
        try:
            self.slp_output_units = str(terrain_data["Slope gradient"]["slp_output_units"]["value"])
        except KeyError:
            pass
        try:
            self.hs_sun_azi = int(terrain_data["Hillshade"]["hs_sun_azi"]["value"])
        except KeyError:
            pass
        try:
            self.hs_sun_el = int(terrain_data["Hillshade"]["hs_sun_el"]["value"])
        except KeyError:
            pass
        try:
            self.mhs_nr_dir = int(terrain_data["Multiple directions hillshade"]["mhs_nr_dir"]["value"])
        except KeyError:
            pass
        try:
            self.mhs_sun_el = int(terrain_data["Multiple directions hillshade"]["mhs_sun_el"]["value"])
        except KeyError:
            pass
        try:
            self.slrm_rad_cell = int(terrain_data["Simple local relief model"]["slrm_rad_cell"]["value"])
        except KeyError:
            pass
        try:
            self.svf_n_dir = int(terrain_data["Sky-View Factor"]["svf_n_dir"]["value"])
        except KeyError:
            pass
        try:
            self.svf_r_max = int(terrain_data["Sky-View Factor"]["svf_r_max"]["value"])
        except KeyError:
            pass
        try:
            self.svf_noise = int(terrain_data["Sky-View Factor"]["svf_noise"]["value"])
        except KeyError:
            pass
        try:
            self.asvf_dir = int(terrain_data["Anisotropic Sky-View Factor"]["asvf_dir"]["value"])
        except KeyError:
            pass
        try:
            self.asvf_level = int(terrain_data["Anisotropic Sky-View Factor"]["asvf_level"]["value"])
        except KeyError:
            pass
        try:
            self.sim_sky_mod = str(terrain_data["Sky illumination"]["sim_sky_mod"]["value"])
        except KeyError:
            pass
        try:
            self.sim_compute_shadow = int(terrain_data["Sky illumination"]["sim_compute_shadow"]["value"])
        except KeyError:
            pass
        try:
            self.sim_nr_dir = int(terrain_data["Sky illumination"]["sim_nr_dir"]["value"])
        except KeyError:
            pass
        try:
            self.sim_shadow_dist = int(terrain_data["Sky illumination"]["sim_shadow_dist"]["value"])
        except KeyError:
            pass
        try:
            self.sim_shadow_az = int(terrain_data["Sky illumination"]["sim_shadow_az"]["value"])
        except KeyError:
            pass
        try:
            self.sim_shadow_el = int(terrain_data["Sky illumination"]["sim_shadow_el"]["value"])
        except KeyError:
            pass
        try:
            self.ld_min_rad = int(terrain_data["Local dominance"]["ld_min_rad"]["value"])
        except KeyError:
            pass
        try:
            self.ld_max_rad = int(terrain_data["Local dominance"]["ld_max_rad"]["value"])
        except KeyError:
            pass
        try:
            self.ld_rad_inc = int(terrain_data["Local dominance"]["ld_rad_inc"]["value"])
        except KeyError:
            pass
        try:
            self.ld_anglr_res = int(terrain_data["Local dominance"]["ld_anglr_res"]["value"])
        except KeyError:
            pass
        try:
            self.ld_observer_h = float(terrain_data["Local dominance"]["ld_observer_h"]["value"])
        except KeyError:
            pass
        try:
            self.msrm_feature_min = float(terrain_data["Multi-scale relief model"]["msrm_feature_min"]["value"])
        except KeyError:
            pass
        try:
            self.msrm_feature_max = float(terrain_data["Multi-scale relief model"]["msrm_feature_max"]["value"])
        except KeyError:
            pass
        try:
            self.msrm_scaling_factor = int(terrain_data["Multi-scale relief model"]["msrm_scaling_factor"]["value"])
        except KeyError:
            pass
        try:
            self.mstp_local_scale = (int(terrain_data["Multi-scale topographic position"]["mstp_local_scale"]["min"]),
                                     int(terrain_data["Multi-scale topographic position"]["mstp_local_scale"]["max"]),
                                     int(terrain_data["Multi-scale topographic position"]["mstp_local_scale"]["step"]))
        except KeyError:
            pass
        try:
            self.mstp_meso_scale = (int(terrain_data["Multi-scale topographic position"]["mstp_meso_scale"]["min"]),
                                    int(terrain_data["Multi-scale topographic position"]["mstp_meso_scale"]["max"]),
                                    int(terrain_data["Multi-scale topographic position"]["mstp_meso_scale"]["step"]))
        except KeyError:
            pass
        try:
            self.mstp_broad_scale = (int(terrain_data["Multi-scale topographic position"]["mstp_broad_scale"]["min"]),
                                     int(terrain_data["Multi-scale topographic position"]["mstp_broad_scale"]["max"]),
                                     int(terrain_data["Multi-scale topographic position"]["mstp_broad_scale"]["step"]))
        except KeyError:
            pass
        try:
            self.mstp_lightness = float(terrain_data["Multi-scale topographic position"]["mstp_lightness"]["value"])
        except KeyError:
            pass
        try:
            self.slp_stretch = (float(terrain_data["Slope gradient"]["stretch"]["min"]),
                                float(terrain_data["Slope gradient"]["stretch"]["max"]))
        except KeyError:
            pass
        try:
            self.hs_stretch = (float(terrain_data["Hillshade"]["stretch"]["min"]),
                               float(terrain_data["Hillshade"]["stretch"]["max"]))
        except KeyError:
            pass
        try:
            self.mhs_stretch = (float(terrain_data["Multiple directions hillshade"]["stretch"]["min"]),
                                float(terrain_data["Multiple directions hillshade"]["stretch"]["max"]))
        except KeyError:
            pass
        try:
            self.slrm_stretch = (float(terrain_data["Simple local relief model"]["stretch"]["min"]),
                                 float(terrain_data["Simple local relief model"]["stretch"]["max"]))
        except KeyError:
            pass
        try:
            self.svf_stretch = (float(terrain_data["Sky-View Factor"]["stretch"]["min"]),
                                float(terrain_data["Sky-View Factor"]["stretch"]["max"]))
        except KeyError:
            pass
        try:
            self.asvf_stretch = (float(terrain_data["Anisotropic Sky-View Factor"]["stretch"]["min"]),
                                 float(terrain_data["Anisotropic Sky-View Factor"]["stretch"]["max"]))
        except KeyError:
            pass
        try:
            self.pos_opns_stretch = (float(terrain_data["Openness - Positive"]["stretch"]["min"]),
                                     float(terrain_data["Openness - Positive"]["stretch"]["max"]))
        except KeyError:
            pass
        try:
            self.neg_opns_stretch = (float(terrain_data["Openness - Negative"]["stretch"]["min"]),
                                     float(terrain_data["Openness - Negative"]["stretch"]["max"]))
        except KeyError:
            pass
        try:
            self.neg_opns_stretch = (float(terrain_data["Sky illumination"]["stretch"]["min"]),
                                     float(terrain_data["Sky illumination"]["stretch"]["max"]))
        except KeyError:
            pass
        try:
            self.neg_opns_stretch = (float(terrain_data["Local dominance"]["stretch"]["min"]),
                                     float(terrain_data["Local dominance"]["stretch"]["max"]))
        except KeyError:
            pass
        try:
            self.msrm_stretch = (float(terrain_data["Multi-scale relief model"]["stretch"]["min"]),
                                 float(terrain_data["Multi-scale relief model"]["stretch"]["max"]))
        except KeyError:
            pass
        try:
            self.mstp_stretch = (float(terrain_data["Multi-scale topographic position"]["stretch"]["min"]),
                                 float(terrain_data["Multi-scale topographic position"]["stretch"]["max"]))
        except KeyError:
            pass

    def apply_terrain(self, default: rvt.default.DefaultValues, combination: BlenderCombination):
        """It overwrites default (DefaultValues) and combination (BlenderCombination),
         with self values that are not None."""
        if self.slp_output_units is not None:
            default.slp_output_units = self.slp_output_units
        if self.hs_sun_azi is not None:
            default.hs_sun_azi = self.hs_sun_azi
        if self.hs_sun_el is not None:
            default.hs_sun_el = self.hs_sun_el
        if self.mhs_nr_dir is not None:
            default.mhs_nr_dir = self.mhs_nr_dir
        if self.mhs_sun_el is not None:
            default.mhs_sun_el = self.mhs_sun_el
        if self.slrm_rad_cell is not None:
            default.slrm_rad_cell = self.slrm_rad_cell
        if self.svf_n_dir is not None:
            default.svf_n_dir = self.svf_n_dir
        if self.svf_r_max is not None:
            default.svf_r_max = self.svf_r_max
        if self.svf_noise is not None:
            default.svf_noise = self.svf_noise
        if self.asvf_dir is not None:
            default.asvf_dir = self.asvf_dir
        if self.asvf_level is not None:
            default.asvf_level = self.asvf_level
        if self.sim_sky_mod is not None:
            default.sim_sky_mod = self.sim_sky_mod
        if self.sim_compute_shadow is not None:
            default.sim_compute_shadow = self.sim_compute_shadow
        if self.sim_nr_dir is not None:
            default.sim_nr_dir = self.sim_nr_dir
        if self.sim_shadow_dist is not None:
            default.sim_shadow_dist = self.sim_shadow_dist
        if self.sim_shadow_az is not None:
            default.sim_shadow_az = self.sim_shadow_az
        if self.sim_shadow_el is not None:
            default.sim_shadow_el = self.sim_shadow_el
        if self.ld_min_rad is not None:
            default.ld_min_rad = self.ld_min_rad
        if self.ld_max_rad is not None:
            default.ld_max_rad = self.ld_max_rad
        if self.ld_rad_inc is not None:
            default.ld_rad_inc = self.ld_rad_inc
        if self.ld_anglr_res is not None:
            default.ld_anglr_res = self.ld_anglr_res
        if self.ld_observer_h is not None:
            default.ld_observer_h = self.ld_observer_h
        if self.msrm_feature_min is not None:
            default.msrm_feature_min = self.msrm_feature_min
        if self.msrm_feature_max is not None:
            default.msrm_feature_max = self.msrm_feature_max
        if self.msrm_scaling_factor is not None:
            default.msrm_scaling_factor = self.msrm_scaling_factor
        if self.mstp_local_scale is not None:
            default.mstp_local_scale = self.mstp_local_scale
        if self.mstp_meso_scale is not None:
            default.mstp_meso_scale = self.mstp_meso_scale
        if self.mstp_broad_scale is not None:
            default.mstp_broad_scale = self.mstp_broad_scale
        if self.mstp_lightness is not None:
            default.mstp_lightness = self.mstp_lightness

        # linear histogram stretches, combination values overwrite
        for layer in combination.layers:
            if layer.vis.lower() == "slope gradient" and self.slp_stretch is not None:
                layer.min = self.slp_stretch[0]
                layer.max = self.slp_stretch[1]
            elif layer.vis.lower() == "hillshade" and self.hs_stretch is not None:
                layer.min = self.hs_stretch[0]
                layer.max = self.hs_stretch[1]
            elif layer.vis.lower() == "multiple directions hillshade" and self.mhs_stretch is not None:
                layer.min = self.mhs_stretch[0]
                layer.max = self.mhs_stretch[1]
            elif layer.vis.lower() == "simple local relief model" and self.slrm_stretch is not None:
                layer.min = self.slrm_stretch[0]
                layer.max = self.slrm_stretch[1]
            elif layer.vis.lower() == "sky-view factor" and self.svf_stretch is not None:
                layer.min = self.svf_stretch[0]
                layer.max = self.svf_stretch[1]
            elif layer.vis.lower() == "anisotropic sky-view factor" and self.asvf_stretch is not None:
                layer.min = self.asvf_stretch[0]
                layer.max = self.asvf_stretch[1]
            elif layer.vis.lower() == "openness - positive" and self.pos_opns_stretch is not None:
                layer.min = self.pos_opns_stretch[0]
                layer.max = self.pos_opns_stretch[1]
            elif layer.vis.lower() == "openness - negative" and self.neg_opns_stretch is not None:
                layer.min = self.neg_opns_stretch[0]
                layer.max = self.neg_opns_stretch[1]
            elif layer.vis.lower() == "sky illumination" and self.sim_stretch is not None:
                layer.min = self.sim_stretch[0]
                layer.max = self.sim_stretch[1]
            elif layer.vis.lower() == "local dominance" and self.ld_stretch is not None:
                layer.min = self.ld_stretch[0]
                layer.max = self.ld_stretch[1]
            elif layer.vis.lower() == "multi-scale relief model" and self.msrm_stretch is not None:
                layer.min = self.msrm_stretch[0]
                layer.max = self.msrm_stretch[1]
            elif layer.vis.lower() == "multi-scale topographic position" and self.mstp_stretch is not None:
                layer.min = self.mstp_stretch[0]
                layer.max = self.mstp_stretch[1]


class TerrainsSettings:
    """Multiple Terrain settings."""

    def __init__(self):
        self.terrains_settings = []

    def read_from_file(self, file_path):
        """Reads combinations from .json file."""
        dat = open(file_path, "r")
        json_data = json.load(dat)
        terrains_settings_json = json_data["terrains_settings"]
        for terrain_json in terrains_settings_json:
            terrain_settings = TerrainSettings()
            terrain_settings.read_from_json(terrain_json)
            self.terrains_settings.append(terrain_settings)
        dat.close()

    def select_terrain_settings_by_name(self, name):
        """Select first combination where self.combinations.BlenderCombination.name = name."""
        for terrain_setting in self.terrains_settings:
            if terrain_setting.name == name:
                return terrain_setting


# Advance blending combinations
def color_relief_image_map(dem, resolution, default: rvt.default.DefaultValues = rvt.default.DefaultValues(),
                           colormap="OrRd", min_colormap_cut=0, max_colormap_cut=1, no_data=None):
    """
    RVT Color relief image map (CRIM)
    Blending combination where layers are:
    1st: Openness positive - Openness negative, overlay, 50% opacity
    2nd: Openness positive - Openness negative, luminosity, 50% opacity
    3rd: Slope gradient, colored with matplotlib colormap

    Parameters
    ----------
    dem : numpy.ndarray
        Input digital elevation model as 2D numpy array.
    resolution : float
        DEM pixel size.
    default : rvt.default.DefaultValues
        Default values for visualization functions.
    slope_norm : tuple(mode, min, max)
        Cutoff normalization for slope.
        Mode can be 'value' or 'percent' (cut-off units).
        Values min and max define stretch borders (in mode units).
    op_on_norm : tuple(mode, min, max)
        Cutoff normalization for (openness positive - openness negative).
        Mode can be 'value' or 'percent' (cut-off units).
        Values min and max define stretch borders (in mode units).
    colormap : str
        Colormap form matplotlib (https://matplotlib.org/3.3.2/tutorials/colors/colormaps.html).
    min_colormap_cut : float
        What lower part of colormap to cut. Between 0 and 1, if 0.2 it cuts off (deletes) 20% of lower colors in colormap.
        If None cut is not applied.
    max_colormap_cut : float
        What upper part of colormap to cut. Between 0 and 1, if 0.8 it cuts off (deletes) 20% of upper colors in colormap.
        If None cut is not applied.
    no_data : int or float
        Value that represents no_data, all pixels with this value are changed to np.nan .

    Returns
    -------
    crim_out : numpy.ndarray
        2D numpy result array of Color relief image map. 
    """
    slope_norm = ("value", 0, 0.8)  # ("value", 0, 50) deg
    op_on_norm = ("value", -28, 28)
    default.svf_r_max = 10

    if no_data is not None:
        dem[dem == no_data] = np.nan

    opns_pos_arr = default.get_sky_view_factor(dem_arr=dem, resolution=resolution,
                                               compute_svf=False, compute_asvf=False, compute_opns=True,
                                               no_data=None)["opns"]
    opns_neg_arr = default.get_sky_view_factor(dem_arr=-1 * dem, resolution=resolution,
                                               compute_svf=False, compute_asvf=False, compute_opns=True,
                                               no_data=None)["opns"]
    opns_pos_neg_arr = opns_pos_arr - opns_neg_arr
    slope_arr = rvt.vis.slope_aspect(dem=dem, resolution_x=resolution, resolution_y=resolution, output_units="radian")[
        "slope"]

    blend_combination = rvt.blend.BlenderCombination()
    blend_combination.create_layer(vis_method="Openness_Pos-Neg", normalization=op_on_norm[0], minimum=op_on_norm[1],
                                   maximum=op_on_norm[2], blend_mode="overlay", opacity=50,
                                   image=opns_pos_neg_arr)
    blend_combination.create_layer(vis_method="Openness_Pos-Neg", normalization=op_on_norm[0], minimum=op_on_norm[1],
                                   maximum=op_on_norm[2], blend_mode="luminosity", opacity=50,
                                   image=opns_pos_neg_arr)
    blend_combination.create_layer(vis_method="slope gradient rad", normalization=slope_norm[0], minimum=slope_norm[1],
                                   maximum=slope_norm[2], blend_mode="normal", opacity=100, colormap=colormap,
                                   min_colormap_cut=min_colormap_cut, max_colormap_cut=max_colormap_cut,
                                   image=slope_arr)
    crim_out = blend_combination.render_all_images()
    return crim_out


def e3mstp(dem, resolution, default: rvt.default.DefaultValues = rvt.default.DefaultValues(), no_data=None):
    """
    RVT enahanced version 3 Multi-scale topographic position (e3MSTP)
    Blending combination where layers are:
    1st: Simple local relief model (SLRM), sreen, 25% opacity
    2nd: Color relief image map where cmap=Reds_r(0.5-1) (CRIM_Reds_r), soft_light, 70% opacity
    3rd: Multi-scale topographic postion (MSTP)

    Parameters
    ----------
    dem : numpy.ndarray
        Input digital elevation model as 2D numpy array.
    resolution : float
        DEM pixel size.
    default : rvt.default.DefaultValues
        Default values for visualization functions.
    no_data : int or float
        Value that represents no_data, all pixels with this value are changed to np.nan .

    Returns
    -------
    crim_out : numpy.ndarray
        2D numpy result array of Color relief image map.
    """
    if no_data is not None:
        dem[dem == no_data] = np.nan
    slrm_arr = default.get_slrm(dem_arr=dem)
    crim_red_arr = color_relief_image_map(dem=dem, resolution=resolution, default=default,
                                          colormap="OrRd", min_colormap_cut=0, max_colormap_cut=1)
    mstp_arr = default.get_mstp(dem_arr=dem)

    blend_combination = rvt.blend.BlenderCombination()
    blend_combination.create_layer(vis_method="slrm", normalization="value",
                                   minimum=-0.5,
                                   maximum=0.5, blend_mode="screen", opacity=25,
                                   image=slrm_arr)
    blend_combination.create_layer(vis_method="crim_red", normalization="value",
                                   minimum=0,
                                   maximum=1, blend_mode="soft_light", opacity=70,
                                   image=crim_red_arr)
    blend_combination.create_layer(vis_method="mstp", normalization="value",
                                   minimum=0,
                                   maximum=255, blend_mode="normal", opacity=100,
                                   image=mstp_arr)
    e3mstp_out = blend_combination.render_all_images()
    return e3mstp_out

