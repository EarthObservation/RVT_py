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
import rasterio as rio
import matplotlib as mpl
import matplotlib.cm
import os
import json


def create_blender_file_example(file_path=None):
    """Create blender .json file example (can be changed and read). Example is VAT - Archaeological combination"""
    data = {"blender_layers": {"name": "VAT - Archaeological",
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


def get_raster_array(raster_path):
    raster_dataset = rio.open(raster_path)
    raster_arr = raster_dataset.read()[0]
    raster_dataset.close()
    return raster_arr


def save_rendered_image(rendered_image, dem_path, save_render_path):
    dem_dataset = rio.open(dem_path)
    profile = dem_dataset.profile
    dem_dataset.close()
    profile.update(dtype='float32')
    rendered_img_dataset = rio.open(save_render_path, "w", **profile)
    rendered_image = rendered_image.astype('float32')
    rendered_img_dataset.write(np.array([rendered_image]))
    rendered_img_dataset.close()


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
    Class which define layer for blending. BlenderLayer is basic element in BlenderLayers.layers list.

    Attributes
    ----------
    dem_arr : np.array (2D)
        Array with DEM data, needed for calculating visualization functions in memory.
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
        self.check_data()

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


class BlenderLayers:
    """
    Class for storing layers (rasters, parameters  for blending) and rendering(blending) into blended raster.

    Attributes
    ----------
    dem_arr : np.array (2D)
        Array with DEM data, needed for calculating visualization functions in memory.
    dem_path : str
        Path to DEM, needed for calculating visualization functions and saving them.
    layers : [BlenderLayer]
        List of BlenderLayer instances which will be blended together.
    """

    def __init__(self, dem_arr=None, dem_resolution=None, dem_path=None):
        self.dem_arr = dem_arr
        self.dem_resolution = dem_resolution
        self.dem_path = dem_path
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
        layer = BlenderLayer(vis_method=vis_method, normalization=normalization, minimum=minimum, maximum=maximum,
                             blend_mode=blend_mode, opacity=opacity, image=image, image_path=image_path)
        self.layers.append(layer)

    def add_layer(self, layer: BlenderLayer):
        """Add BlenderLayer instance to layers attribute."""
        self.layers.append(layer)

    def remove_all_layers(self):
        """Empties layers attribute."""
        self.layers = []

    def build_blender_layers_from_file(self, file_path):
        """Fill layers from .json file."""
        self.layers = []
        # Example file (for file_path) in dir settings: blender_file_example.txt
        dat = open(file_path, "r")
        data = json.load(dat)
        layers_data = data["blender_layers"]["layers"]
        for layer in layers_data:
            layer_name = layer["layer"]
            if layer["visualization_method"] == "None":
                vis_method = None
                continue
            else:
                vis_method = layer["visualization_method"]
            norm = layer["norm"]
            norm_min = layer["min"]
            norm_max = layer["max"]
            blend_mode = layer["blend_mode"]
            opacity = layer["opacity"]
            self.add_layer(BlenderLayer(vis_method=vis_method, normalization=norm, minimum=norm_min, maximum=norm_max,
                                        blend_mode=blend_mode, opacity=opacity))
        dat.close()

    def save_blender_layers_to_file(self, file_path):
        """Save layers (manually) to .json file. Param image and image_path have to be None, vis has to be correct!"""
        data = {"blender_layers": {"name": "VAT - Archaeological",
                                   "layers": []
                                   }}
        i = 1
        for layer in self.layers:
            data["blender_layers"]["layers"].append({"layer": str(i), "visualization_method": layer.vis,
                                                     "norm": layer.normalization, "min": layer.min, "max": layer.max,
                                                     "blend_mode": layer.blend_mode, "opacity": layer.opacity})
            i += 1
        dat = open(file_path, "w")
        dat.write(json.dumps(data, indent=4))
        dat.close()

    def render_all_images(self, default=None, save_visualizations=False, save_render_path=None):
        """Render all layers and returns blended image. If specific layer (BlenderLayer) in layers has image
            (is not None), method uses this image, if image is None and layer has image_path method reads image from path.
            If both image and image_path are None method calculates visualization. If save_visualization is True method
            needs dem_path and saves each visualization (if it doesn't exists) in directory of dem_path,
            else (save_visualization=False) method needs dem_arr, dem_resolution and calculates each visualization
            simultaneously (in memory). Be careful save_visualisation applies only if specific BlenderLayer
            image and image_path are None"""
        if save_render_path is not None and self.dem_path is None:
            raise Exception(
                "rvt.blend.BlenderLayers.render_all_images: If you would like to save rendered image (blender), "
                "you have to define dem_path (BlenderLayers.add_dem_path())!")

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
                raise Exception("rvt.blend.BlenderLayers.render_all_images: If you would like to save visualizations, "
                                "you have to define dem_path (BlenderLayers.add_dem_path())!")
            if not save_visualizations and self.dem_arr is None and self.dem_resolution is None and \
                    image_path is None and image is None:
                raise Exception(
                    "rvt.blend.BlenderLayers.render_all_images: If you would like to compute visualizations, "
                    "you have to define dem_arr and its resolution (BlenderLayers.add_dem_arr())!")

            # normalize images
            # if image is not presented and image_path is
            if image is None and image_path is not None:
                norm_image = normalize_image(visualization, get_raster_array(image_path), min_norm,
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
                        norm_image = normalize_image(visualization, get_raster_array(image_path), min_norm, max_norm,
                                                     normalization)
                    else:
                        image = default.get_slope(dem_arr=self.dem_arr, resolution_x=self.dem_resolution,
                                                  resolution_y=self.dem_resolution)["slope"]
                        norm_image = normalize_image(visualization, image, min_norm, max_norm, normalization)
                elif self.layers[i_img].vis.lower() == "hillshade":
                    if save_visualizations:
                        default.save_hillshade(dem_path=self.dem_path)
                        image_path = default.get_hillshade_path(self.dem_path)
                        norm_image = normalize_image(visualization, get_raster_array(image_path), min_norm, max_norm,
                                                     normalization)
                    else:
                        image = default.get_hillshade(dem_arr=self.dem_arr, resolution_x=self.dem_resolution,
                                                      resolution_y=self.dem_resolution)
                        norm_image = normalize_image(visualization, image, min_norm, max_norm, normalization)
                elif self.layers[i_img].vis.lower() == "multiple directions hillshade":
                    if save_visualizations:
                        default.save_multi_hillshade(dem_path=self.dem_path)
                        image_path = default.get_multi_hillshade_path(self.dem_path)
                        norm_image = normalize_image(visualization, get_raster_array(image_path), min_norm, max_norm,
                                                     normalization)
                    else:
                        image = default.get_multi_hillshade(dem_arr=self.dem_arr, resolution_x=self.dem_resolution,
                                                            resolution_y=self.dem_resolution)
                        norm_image = normalize_image(visualization, image, min_norm, max_norm, normalization)
                elif self.layers[i_img].vis.lower() == "simple local relief model":
                    if save_visualizations:
                        default.save_slrm(dem_path=self.dem_path)
                        image_path = default.get_slrm_path(self.dem_path)
                        norm_image = normalize_image(visualization, get_raster_array(image_path), min_norm, max_norm,
                                                     normalization)
                    else:
                        image = default.get_slrm(dem_arr=self.dem_arr)
                        norm_image = normalize_image(visualization, image, min_norm, max_norm, normalization)
                elif self.layers[i_img].vis.lower() == "sky-view factor":
                    if save_visualizations:
                        default.save_sky_view_factor(dem_path=self.dem_path, save_svf=True, save_asvf=False,
                                                     save_opns=False)
                        image_path = default.get_svf_path(self.dem_path)
                        norm_image = normalize_image(visualization, get_raster_array(image_path), min_norm, max_norm,
                                                     normalization)
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
                        norm_image = normalize_image(visualization, get_raster_array(image_path), min_norm, max_norm,
                                                     normalization)
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
                        norm_image = normalize_image(visualization, get_raster_array(image_path), min_norm, max_norm,
                                                     normalization)
                    else:
                        image = default.get_sky_view_factor(dem_arr=self.dem_arr, resolution=self.dem_resolution,
                                                            compute_svf=False, compute_asvf=False,
                                                            compute_opns=True)["opns"]
                        norm_image = normalize_image(visualization, image, min_norm, max_norm, normalization)
                elif self.layers[i_img].vis.lower() == "openness - negative":
                    if save_visualizations:
                        default.save_neg_opns(dem_path=self.dem_path)
                        image_path = default.get_neg_opns_path(self.dem_path)
                        norm_image = normalize_image(visualization, get_raster_array(image_path), min_norm, max_norm,
                                                     normalization)
                    else:
                        image = default.get_neg_opns(dem_arr=self.dem_arr, resolution=self.dem_resolution)
                        norm_image = normalize_image(visualization, image, min_norm, max_norm, normalization)
                elif self.layers[i_img].vis.lower() == "sky illumination":
                    if save_visualizations:
                        default.save_sky_illumination(dem_path=self.dem_path)
                        image_path = default.get_sky_illumination_path(self.dem_path)
                        norm_image = normalize_image(visualization, get_raster_array(image_path), min_norm, max_norm,
                                                     normalization)
                    else:
                        image = default.get_sky_illumination(dem_arr=self.dem_arr, resolution=self.dem_resolution)
                        norm_image = normalize_image(visualization, image, min_norm, max_norm, normalization)
                elif self.layers[i_img].vis.lower() == "local dominance":
                    if save_visualizations:
                        default.save_local_dominance(dem_path=self.dem_path)
                        image_path = default.get_local_dominance_path(self.dem_path)
                        norm_image = normalize_image(visualization, get_raster_array(image_path), min_norm, max_norm,
                                                     normalization)
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
                    warnings.warn("rvt.blend.BlenderLayers.render_all_images: Rendered image scale distorted")
        if save_render_path is not None:  # if paths presented it saves image
            save_rendered_image(rendered_image, dem_path=self.dem_path, save_render_path=save_render_path)
        return rendered_image
