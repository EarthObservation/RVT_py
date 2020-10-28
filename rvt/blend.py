"""
Relief Visualization Toolbox – Visualization Functions

Contains functions for blending.

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

# TODO: more testing, find and fix bugs if they exists

# python libraries
import numpy as np
import warnings
import rvt.default
import matplotlib as mpl
import matplotlib.cm
import os
import json


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


def gray_scale_to_color_ramp(gray_scale, colormap, alpha=False):
    """
    Turns normalized gray scale np.array to rgba (np.array of 4 np.arrays r, g, b, a).

    Parameters
    ----------
    gray_scale : np.array (2D)
        Normalized gray_scale img as np.array (0-1)
    colormap : str
        colormap form matplotlib (https://matplotlib.org/3.3.2/tutorials/colors/colormaps.html)
    alpha : bool
        If True outputs 4D array RGBA, if False outputs 3D array RGB
    Returns
    -------
    rgba_out : np.array (3D: red 0-255, green 0-255, blue 0-255)
            If alpha False: np.array (4D: red 0-255, green 0-255, blue 0-255, alpha 0-255)
    """
    cm = mpl.cm.get_cmap(colormap)
    rgba_out = cm(gray_scale)  # normalized rgb
    rgba_mtpl_out = np.uint8(rgba_out * 255)  # 0-1 scale to 0-255 and change type to uint8
    if alpha:
        rgba_out = np.array([rgba_mtpl_out[:, :, 0], rgba_mtpl_out[:, :, 1], rgba_mtpl_out[:, :, 2]])
    else:
        rgba_out = np.array([rgba_mtpl_out[:, :, 0], rgba_mtpl_out[:, :, 1], rgba_mtpl_out[:, :, 2],
                             rgba_mtpl_out[:, :, 3]])
    return rgba_out


def normalize_lin(image, minimum, maximum):
    # linear cut off
    idx_min = np.where(image < minimum)
    idx_max = np.where(image > maximum)
    if idx_min[0].size == 0 and idx_max[0].size == 0:
        return image
    if idx_min[0].size > 0:
        image[idx_min] = minimum
    if idx_max[0].size > 0:
        image[idx_max] = maximum

    # stretch to 0.0 - 1.0 interval
    image = (image - minimum) / (maximum - minimum)
    return image


def lin_cutoff_calc_from_perc(image, minimum, maximum):
    if 1 < minimum < 100:
        minimum = minimum / 100
    if 1 < maximum < 100:
        maximum = maximum / 100

    distribution = np.percentile(a=image, q=np.array([minimum, 1 - maximum]))
    min_lin = np.amin(distribution)
    max_lin = np.amax(distribution)
    return {"min_lin": min_lin, "max_lin": max_lin}


def normalize_perc(image, minimum, maximum):
    min_max_lin_dict = lin_cutoff_calc_from_perc(image, minimum, maximum)
    min_lin = min_max_lin_dict["min_lin"]
    max_lin = min_max_lin_dict["max_lin"]
    return normalize_lin(image, min_lin, max_lin)


def advanced_normalization(image, minimum, maximum, normalization):
    equ_image = image
    if normalization.lower() == "value":
        equ_image = normalize_lin(image, minimum, maximum)
    elif normalization.lower() == "perc":
        equ_image = normalize_perc(image, minimum, maximum)
    elif normalization is None:
        equ_image = image
    return equ_image


def image_join_channels(r, g, b):
    if r.shape != g.shape or r.shape != b.shape or g.shape != b.shape:
        raise Exception("rvt.blend.image_join_channels: r, g, b must me same dimensions!")
    return np.array([r, g, b])


def lum(img):
    if len(img.shape) == 3:
        r = img[0]
        g = img[1]
        b = img[2]
        return (0.3 * r) + (0.59 * g) + (0.11 * b)
    else:
        return img


def matrix_eq_min_lt_zero(r, idx_min_lt_zero, lum_c, min_c):
    r[idx_min_lt_zero] = lum_c[idx_min_lt_zero] + (((r[idx_min_lt_zero] - lum_c[idx_min_lt_zero]) *
                                                    lum_c[idx_min_lt_zero])
                                                   / (lum_c[idx_min_lt_zero] - min_c[idx_min_lt_zero]))
    return r


def matrix_eq_max_gt_one(r, idx_max_c_gt_one, lum_c, max_c):
    r[idx_max_c_gt_one] = lum_c[idx_max_c_gt_one] + (((r[idx_max_c_gt_one] - lum_c[idx_max_c_gt_one]) *
                                                      (1.0 - lum_c[idx_max_c_gt_one]))
                                                     / (max_c[idx_max_c_gt_one] - lum_c[idx_max_c_gt_one]))
    return r


def channel_min(r, g, b):
    min_c = r
    idx_min = np.where(g < min_c)
    min_c[idx_min] = g[idx_min]
    idx_min = np.where(b < min_c)
    min_c[idx_min] = b[idx_min]
    return min_c


def channel_max(r, g, b):
    max_c = r
    idx_min = np.where(g > max_c)
    max_c[idx_min] = g[idx_min]
    idx_min = np.where(b > max_c)
    max_c[idx_min] = b[idx_min]
    return max_c


def clip_color(c, min_c=None, max_c=None):
    lum_c = lum(c)
    r = c[0]
    g = c[1]
    b = c[2]

    if min_c is None and max_c is None:
        min_c = channel_min(r, g, b)
        max_c = channel_max(r, g, b)

    idx_min_lt_zero = np.where(min_c < 0)
    r = matrix_eq_min_lt_zero(r, idx_min_lt_zero, lum_c, min_c)
    g = matrix_eq_min_lt_zero(g, idx_min_lt_zero, lum_c, min_c)
    b = matrix_eq_min_lt_zero(b, idx_min_lt_zero, lum_c, min_c)

    idx_max_c_gt_one = np.where(max_c > 1)
    r = matrix_eq_max_gt_one(r, idx_max_c_gt_one, lum_c, max_c)
    g = matrix_eq_max_gt_one(g, idx_max_c_gt_one, lum_c, max_c)
    b = matrix_eq_max_gt_one(b, idx_max_c_gt_one, lum_c, max_c)

    c[0, :, :] = r
    c[1, :, :] = g
    c[2, :, :] = b

    return c


def blend_normal(active, background):
    return active


def blend_screen(active, background):
    return 1 - (1 - active) * (1 - background)


def blend_multiply(active, background):
    return active * background


def blend_overlay(active, background):
    idx1 = np.where(background > 0.5)
    idx2 = np.where(background <= 0.5)
    background[idx1[0], idx1[1]] = (
            1 - (1 - 2 * (background[idx1[0], idx1[1]] - 0.5)) * (1 - active[idx1[0], idx1[1]]))
    background[idx2[0], idx2[1]] = ((2 * background[idx2[0], idx2[1]]) * active[idx2[0], idx2[1]])
    return background


def blend_luminosity(active, background, min_c=None, max_c=None):
    lum_active = lum(active)
    lum_background = lum(background)
    luminosity = lum_active - lum_background
    if len(background.shape) < 3:
        return lum_active

    r = background[0] + luminosity
    g = background[1] + luminosity
    b = background[2] + luminosity

    c = np.zeros(background.shape)
    c[0, :, :] = r
    c[1, :, :] = g
    c[2, :, :] = b

    clipped_image = clip_color(c, min_c, max_c)

    return clipped_image


def equation_blend(blend_mode, active, background):
    if blend_mode.lower() == "screen":
        return blend_screen(active, background)
    elif blend_mode.lower() == "multiply":
        return blend_multiply(active, background)
    elif blend_mode.lower() == "overlay":
        return blend_overlay(active, background)


def blend_multi_dim_images(blend_mode, active, background):
    a_rgb = len(active.shape) == 3  # bool, is active rgb
    b_rgb = len(background.shape) == 3  # bool, is background rgb
    blended_image = None
    if a_rgb and b_rgb:
        blended_image = np.zeros(background.shape)
        for i in range(3):
            blended_image[i, :, :] = equation_blend(blend_mode, active[i, :, :], background[i, :, :])
    if a_rgb and not b_rgb:
        blended_image = np.zeros(active.shape)
        for i in range(3):
            blended_image[i, :, :] = equation_blend(blend_mode, active[i, :, :], background)
    if not a_rgb and b_rgb:
        blended_image = np.zeros(background.shape)
        for i in range(3):
            blended_image[i, :, :] = equation_blend(blend_mode, active, background[i, :, :])
    if not a_rgb and not b_rgb:
        blended_image = equation_blend(blend_mode, active, background)

    return blended_image


def blend_images(blend_mode, active, background, min_c=None, max_c=None):
    if blend_mode.lower() == "multiply" or blend_mode.lower() == "overlay" or blend_mode.lower() == "screen":
        return blend_multi_dim_images(blend_mode, active, background)
    elif blend_mode.lower() == "luminosity":
        return blend_luminosity(active, background, min_c, max_c)
    else:
        return blend_normal(active, background)


def render_images(active, background, opacity):
    if np.nanmin(active) < 0 or np.nanmax(active) > 1.1:
        active = scale_0_to_1(active)
    if np.nanmin(background) < 0 or np.nanmax(background) > 1.1:
        background = scale_0_to_1(background)

    a_rgb = len(active.shape) == 3
    b_rgb = len(background.shape) == 3
    render_image = 0
    if a_rgb and b_rgb:
        render_image = np.zeros(background.shape)
        for i in range(3):
            render_image[i, :, :] = apply_opacity(active[i, :, :], background[i, :, :], opacity)
    if a_rgb and not b_rgb:
        render_image = np.zeros(active.shape)
        for i in range(3):
            render_image[i, :, :] = apply_opacity(active[i, :, :], background, opacity)
    if not a_rgb and b_rgb:
        render_image = np.zeros(background.shape)
        for i in range(3):
            render_image[i, :, :] = apply_opacity(active, background[i, :, :], opacity)
    if not a_rgb and not b_rgb:
        render_image = apply_opacity(active, background, opacity)
    return render_image


def scale_within_0_and_1(numeric_value):
    if np.nanmin(numeric_value) >= 0 and np.nanmax(numeric_value) <= 1:
        return numeric_value

    numeric_value[np.isnan(numeric_value)] = np.nanmin(numeric_value)  # nan change to nanmin

    actual_min = np.nanmin(numeric_value)
    norm_min_value = np.nanmax(np.array(0, actual_min))

    actual_max = np.nanmax(numeric_value)
    norm_max_value = np.nanmin(np.array(1, actual_max))

    # do not scale values where max is between 1 and 255 if the max-min values diffrence is at least 30 and min >0
    # and numeric values are integer type
    if 255 >= actual_max > 1:
        if actual_max - actual_min > 30 and actual_min > 0:
            scaled = numeric_value / 255
            return scaled

    scaled = (numeric_value - norm_min_value) / (norm_max_value - norm_min_value)

    if np.nanmin(scaled) > -0.01:
        scaled[0 > scaled > -0.01] = 0

    return scaled


def scale_strict_0_to_1(numeric_value):
    if np.nanmin(numeric_value) == 0 and np.nanmax(numeric_value) == 1:
        return numeric_value

    numeric_value[np.isnan(numeric_value)] = 0  # nan change to 0

    min_value = np.nanmin(numeric_value)
    max_value = np.nanmax(numeric_value)

    scaled = (numeric_value - min_value) / (max_value - min_value)

    if np.nanmin(scaled) > -0.01:
        scaled[0 > scaled > -0.01] = 0

    return scaled


def scale_0_to_1(numeric_value):
    if 1 >= np.nanmax(numeric_value) > 0.9 and np.nanmin(numeric_value) == 0:
        return numeric_value
    elif np.nanmax(numeric_value) - np.nanmin(numeric_value) > 0.3:
        return scale_within_0_and_1(numeric_value)
    else:
        return scale_strict_0_to_1(numeric_value)


def apply_opacity(active, background, opacity):
    if opacity > 1:
        opacity = opacity / 100
    return active * opacity + background * (1 - opacity)


def normalize_image(visualization, image, min_norm, max_norm, normalization):
    if visualization is None:
        return None
    # workaround for RGB images because they are on scale [0, 255] not [0, 1],
    # we use multiplier to get proper values
    if normalization.lower() == "value" and visualization.lower() == "hillshade":
        if np.nanmax(image) > 100.0 and len(image.shape) == 3:
            # limit normalization 0 to 1
            # all numbers below are 0
            # numbers above are 1
            if min_norm < 0:
                min_norm = 0
            if max_norm > 1:
                max_norm = 1

            min_norm = round(min_norm * 255)
            max_norm = round(max_norm * 255)

    norm_image = advanced_normalization(image, min_norm, max_norm, normalization)
    # make sure it scales 0 to 1
    if np.nanmax(norm_image) > 1:
        if visualization.lower() == "multiple directions hillshade":
            norm_image = scale_0_to_1(norm_image)
        else:
            norm_image = scale_0_to_1(norm_image)
            warnings.warn("rvt.blend.normalize_images_on_layers: unexpected values! max > 1")
        if np.nanmin(norm_image) < 0:
            warnings.warn("rvt.blend.normalize_images_on_layers: unexpected values! min < 0")

    # for slope and neg openness, invert scale
    # meaning high slopes will be black
    if visualization.lower() == "openness - negative" or visualization.lower() == "slope gradient":
        norm_image = 1 - norm_image
    return norm_image


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
    image_path : str
        Path to DEM. Doesn't matter if image is not None. Leave None if you would like for blender to compute it.
    image : numpy.array (2D)
        Visualization raster. Leave None if you would like for blender to compute it.
    """

    def __init__(self, vis_method=None, normalization="value", minimum=None, maximum=None,
                 blend_mode="normal", opacity=100, image=None, image_path=None):
        self.vis = vis_method
        self.normalization = normalization
        self.min = minimum
        self.max = maximum
        self.blend_mode = blend_mode
        self.opacity = opacity
        self.image_path = image_path
        self.image = image

    def check_data(self):
        """ Check Attributes """
        if self.vis is None:  # if vis is None everything is None
            self.normalization = None
            self.min = None
            self.max = None
            self.blend_mode = None
            self.opacity = None
            self.image = None  # leave None if you wish for blender to compute visualization
            self.image_path = None  # leave None if you wish for blender to compute visualization
        else:
            if self.normalization.lower() != "value" and self.normalization.lower() != "perc":
                raise Exception("rvt.blend.BlenderLayer.check_data: normalization value incorrect!")
            if 0 > self.min > 100:
                raise Exception("rvt.blend.BlenderLayer.check_data: min value incorrect [0-100]!")
            if 0 > self.max > 100:
                raise Exception("rvt.blend.BlenderLayer.check_data: max value incorrect [0-100]!")
            if self.min > self.max:
                raise Exception("rvt.blend.BlenderLayer.check_data: min bigger than max!")
            if self.blend_mode.lower() != "normal" and self.blend_mode.lower() != "multiply" and \
                    self.blend_mode.lower() != "overlay" and self.blend_mode.lower() != "luminosity" and \
                    self.blend_mode.lower() != "screen":
                raise Exception("rvt.blend.BlenderLayer.check_data: blend_mode incorrect!")
            if 0 > self.opacity > 100:
                raise Exception("rvt.blend.BlenderLayer.check_data: opacity incorrect [0-100]!")
            if self.image is None and self.image_path is None:
                if self.vis.lower() != "slope gradient" and self.vis.lower() != "hillshade" and \
                        self.vis.lower() != "multiple directions hillshade" and \
                        self.vis.lower() != "simple local relief model" and self.vis.lower() != "sky-view factor" and \
                        self.vis.lower() != "anisotropic sky-view factor" and \
                        self.vis.lower() != "openness - positive" and self.vis.lower() != "openness - negative" and \
                        self.vis.lower() != "sky illumination" and self.vis.lower() != "local dominance":
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
                     blend_mode="normal", opacity=100, image=None, image_path=None):
        """Create BlenderLayer and adds it to layers attribute."""
        if vis_method is not None:
            layer = BlenderLayer(vis_method=vis_method, normalization=normalization, minimum=minimum, maximum=maximum,
                                 blend_mode=blend_mode, opacity=opacity, image=image, image_path=image_path)
            self.layers.append(layer)

    def add_layer(self, layer: BlenderLayer):
        """Add BlenderLayer instance to layers attribute."""
        if layer.vis is not None:
            self.layers.append(layer)

    def remove_all_layers(self):
        """Empties layers attribute."""
        self.layers = []

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
            self.add_layer(BlenderLayer(vis_method=vis_method, normalization=norm, minimum=norm_min, maximum=norm_max,
                                        blend_mode=blend_mode, opacity=opacity))

    def save_to_file(self, file_path):
        """Save layers (manually) to .json file. Parameters image and image_path in each layer have to be None,
         vis has to be correct!"""
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
            json_data["combination"]["layers"].append({"layer": str(i_layer), "visualization_method": layer.vis,
                                                       "norm": layer.normalization, "min": layer.min,
                                                       "max": layer.max, "blend_mode": layer.blend_mode,
                                                       "opacity": layer.opacity})
            i_layer += 1
        return json_data

    def check_data(self):
        for layer in self.layers:
            layer.check_data()

    def render_all_images(self, default=None, save_visualizations=False, save_render_path=None):
        """Render all layers and returns blended image. If specific layer (BlenderLayer) in layers has image
        (is not None), method uses this image, if image is None and layer has image_path method reads image from
        path. If both image and image_path are None method calculates visualization. If save_visualization is True
        method needs dem_path and saves each visualization (if it doesn't exists) in directory of dem_path,
        else (save_visualization=False) method needs dem_arr, dem_resolution and calculates each visualization
        simultaneously (in memory). Be careful save_visualisation applies only if specific BlenderLayer
        image and image_path are None"""

        # check data
        self.check_data()

        if save_render_path is not None and self.dem_path is None:
            raise Exception(
                "rvt.blend.BlenderCombination.render_all_images: If you would like to save rendered image (blender), "
                "you have to define dem_path (BlenderCombination.add_dem_path())!")

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
                norm_image = normalize_image(visualization, rvt.default.get_raster_arr(image_path)["array"], min_norm,
                                             max_norm, normalization)
            # if image is presented
            elif image is not None:
                norm_image = normalize_image(visualization, image, min_norm, max_norm, normalization)
            # they are both none
            else:
                if self.layers[i_img].vis.lower() == "slope gradient":
                    if save_visualizations:
                        default.save_slope(dem_path=self.dem_path)
                        image_path = default.get_slope_path(self.dem_path)
                        norm_image = normalize_image(visualization, rvt.default.get_raster_arr(image_path)["array"],
                                                     min_norm, max_norm, normalization)
                    else:
                        image = default.get_slope(dem_arr=self.dem_arr, resolution_x=self.dem_resolution,
                                                  resolution_y=self.dem_resolution)["slope"]
                        norm_image = normalize_image(visualization, image, min_norm, max_norm, normalization)
                elif self.layers[i_img].vis.lower() == "hillshade":
                    if save_visualizations:
                        default.save_hillshade(dem_path=self.dem_path)
                        image_path = default.get_hillshade_path(self.dem_path)
                        norm_image = normalize_image(visualization, rvt.default.get_raster_arr(image_path)["array"],
                                                     min_norm, max_norm, normalization)
                    else:
                        image = default.get_hillshade(dem_arr=self.dem_arr, resolution_x=self.dem_resolution,
                                                      resolution_y=self.dem_resolution)
                        norm_image = normalize_image(visualization, image, min_norm, max_norm, normalization)
                elif self.layers[i_img].vis.lower() == "multiple directions hillshade":
                    if save_visualizations:
                        default.save_multi_hillshade(dem_path=self.dem_path)
                        image_path = default.get_multi_hillshade_path(self.dem_path)
                        norm_image = normalize_image(visualization, rvt.default.get_raster_arr(image_path)["array"],
                                                     min_norm, max_norm, normalization)
                    else:
                        image = default.get_multi_hillshade(dem_arr=self.dem_arr, resolution_x=self.dem_resolution,
                                                            resolution_y=self.dem_resolution)
                        norm_image = normalize_image(visualization, image, min_norm, max_norm, normalization)
                elif self.layers[i_img].vis.lower() == "simple local relief model":
                    if save_visualizations:
                        default.save_slrm(dem_path=self.dem_path)
                        image_path = default.get_slrm_path(self.dem_path)
                        norm_image = normalize_image(visualization, rvt.default.get_raster_arr(image_path)["array"],
                                                     min_norm, max_norm, normalization)
                    else:
                        image = default.get_slrm(dem_arr=self.dem_arr)
                        norm_image = normalize_image(visualization, image, min_norm, max_norm, normalization)
                elif self.layers[i_img].vis.lower() == "sky-view factor":
                    if save_visualizations:
                        default.save_sky_view_factor(dem_path=self.dem_path, save_svf=True, save_asvf=False,
                                                     save_opns=False)
                        image_path = default.get_svf_path(self.dem_path)
                        norm_image = normalize_image(visualization, rvt.default.get_raster_arr(image_path)["array"],
                                                     min_norm, max_norm, normalization)
                    else:
                        image = default.get_sky_view_factor(dem_arr=self.dem_arr, resolution=self.dem_resolution,
                                                            compute_svf=True, compute_asvf=False,
                                                            compute_opns=False)["svf"]
                        norm_image = normalize_image(visualization, image, min_norm, max_norm, normalization)
                elif self.layers[i_img].vis.lower() == "anisotropic sky-view factor":
                    if save_visualizations:
                        default.save_sky_view_factor(dem_path=self.dem_path, save_svf=False, save_asvf=True,
                                                     save_opns=False)
                        image_path = default.get_asvf_path(self.dem_path)
                        norm_image = normalize_image(visualization, rvt.default.get_raster_arr(image_path)["array"],
                                                     min_norm, max_norm, normalization)
                    else:
                        image = default.get_sky_view_factor(dem_arr=self.dem_arr, resolution=self.dem_resolution,
                                                            compute_svf=False, compute_asvf=True,
                                                            compute_opns=False)["asvf"]
                        norm_image = normalize_image(visualization, image, min_norm, max_norm, normalization)
                elif self.layers[i_img].vis.lower() == "openness - positive":
                    if save_visualizations:
                        default.save_sky_view_factor(dem_path=self.dem_path, save_svf=False, save_asvf=False,
                                                     save_opns=True)
                        image_path = default.get_opns_path(self.dem_path)
                        norm_image = normalize_image(visualization, rvt.default.get_raster_arr(image_path)["array"],
                                                     min_norm, max_norm, normalization)
                    else:
                        image = default.get_sky_view_factor(dem_arr=self.dem_arr, resolution=self.dem_resolution,
                                                            compute_svf=False, compute_asvf=False,
                                                            compute_opns=True)["opns"]
                        norm_image = normalize_image(visualization, image, min_norm, max_norm, normalization)
                elif self.layers[i_img].vis.lower() == "openness - negative":
                    if save_visualizations:
                        default.save_neg_opns(dem_path=self.dem_path)
                        image_path = default.get_neg_opns_path(self.dem_path)
                        norm_image = normalize_image(visualization, rvt.default.get_raster_arr(image_path)["array"],
                                                     min_norm, max_norm, normalization)
                    else:
                        image = default.get_neg_opns(dem_arr=self.dem_arr, resolution=self.dem_resolution)
                        norm_image = normalize_image(visualization, image, min_norm, max_norm, normalization)
                elif self.layers[i_img].vis.lower() == "sky illumination":
                    if save_visualizations:
                        default.save_sky_illumination(dem_path=self.dem_path)
                        image_path = default.get_sky_illumination_path(self.dem_path)
                        norm_image = normalize_image(visualization, rvt.default.get_raster_arr(image_path)["array"],
                                                     min_norm, max_norm, normalization)
                    else:
                        image = default.get_sky_illumination(dem_arr=self.dem_arr, resolution=self.dem_resolution)
                        norm_image = normalize_image(visualization, image, min_norm, max_norm, normalization)
                elif self.layers[i_img].vis.lower() == "local dominance":
                    if save_visualizations:
                        default.save_local_dominance(dem_path=self.dem_path)
                        image_path = default.get_local_dominance_path(self.dem_path)
                        norm_image = normalize_image(visualization, rvt.default.get_raster_arr(image_path)["array"],
                                                     min_norm, max_norm, normalization)
                    else:
                        image = default.get_local_dominance(dem_arr=self.dem_arr)
                        norm_image = normalize_image(visualization, image, min_norm, max_norm, normalization)

            # if current layer has visualization applied, but there has been no rendering
            # of images yet, than current layer will be the initial value of rendered_image
            if rendered_image == []:
                rendered_image = norm_image
                continue
            else:
                active = norm_image
                background = rendered_image
                blend_mode = self.layers[i_img].blend_mode
                opacity = self.layers[i_img].opacity
                if np.nanmin(active) < 0 or np.nanmax(active) > 1:
                    active = scale_0_to_1(active)
                if np.nanmin(background) < 0 or np.nanmax(background) > 1:
                    background = scale_0_to_1(background)
                top = blend_images(blend_mode, active, background)
                rendered_image = render_images(top, background, opacity)

                if np.nanmin(background) < 0 or np.nanmax(background > 1):
                    warnings.warn("rvt.blend.BlenderCombination.render_all_images: Rendered image scale distorted")
        if save_render_path is not None:  # if paths presented it saves image
            rvt.default.save_raster(src_raster_path=self.dem_path, out_raster_path=save_render_path,
                                    out_raster_arr=rendered_image)
        return rendered_image


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

    def add_combination(self, combination: BlenderCombination):
        """Adds cobmination."""
        self.combinations.append(combination)

    def remove_all_combinations(self):
        """Removes all combinations from self.combinations."""
        self.combinations = []

    def select_combination_by_name(self, name):
        """Select first combination where self.combinations.BlenderCombination.name = name."""
        for combination in self.combinations:
            if combination.name == name:
                return combination

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
        self.sim_samp_pnts = None
        self.sim_shadow_dist = None
        self.sim_shadow_az = None
        self.sim_shadow_el = None
        # local dominance
        self.ld_min_rad = None
        self.ld_max_rad = None
        self.ld_rad_inc = None
        self.ld_anglr_res = None
        self.ld_observer_h = None
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
            self.sim_samp_pnts = int(terrain_data["Sky illumination"]["sim_samp_pnts"]["value"])
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
        if self.sim_samp_pnts is not None:
            default.sim_samp_pnts = self.sim_samp_pnts
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

    def terrain_sett_in_terrains_sett(self, input_terrain_sett: TerrainSettings):
        """Checks if terrain consist one of default terrains in terrains."""
        for terrain_sett in self.terrains_settings:
            if compare_2_terrains_settings(terrain_sett=input_terrain_sett, default_terrain_sett=terrain_sett):
                return terrain_sett.name
        return None


def compare_2_terrains_settings(terrain_sett: TerrainSettings, default_terrain_sett: TerrainSettings):
    """Checks if terrain consist elements of default terrain. Loop ignores default terrain None attributes."""
    class_attributes = list(vars(TerrainSettings()).keys())
    dict_terrain_sett = vars(terrain_sett)  # all class attributes to dict
    dict_default_terrain_sett = vars(default_terrain_sett)  # all class attributes to dict
    for attribute in class_attributes:
        if dict_default_terrain_sett[attribute] is not None and \
                dict_default_terrain_sett[attribute] != dict_terrain_sett[attribute]:
            return False
    return True
