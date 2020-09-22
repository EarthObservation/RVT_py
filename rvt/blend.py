"""
NAME:
    rvt_py, rvt.blend - rvt blend functions

DESCRIPTION:
    Contains all functions for blending. Functions are rewritten from RVT (IDL).
    https://github.com/EarthObservation/RVT

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
# TODO: more testing, find and fix bugs if they exists

# python libraries
import numpy as np
import warnings
import rvt.default


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
    if normalization.lower() == "value":
        equ_image = normalize_lin(image, minimum, maximum)
    elif normalization.lower() == "perc":
        equ_image = normalize_perc(image, minimum, maximum)
    elif normalization is None:
        equ_image = image
    return equ_image


def image_join_channels(r, g, b):
    if r.shape != g.shape or r.shape != b.shape or g.shape != b.shape:
        raise Exception("RVT image_join_channels: r, g, b must me same dimensions!")
    return np.array([r, g, b])


def lum(img):
    if img.shape == 3:
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


class BlenderLayer:
    def __init__(self, vis_method=None, normalization="value", minimum=None, maximum=None,
                 blend_mode="normal", opacity=100, image=None):
        self.vis = vis_method
        self.normalization = normalization
        self.min = minimum
        self.max = maximum
        self.blend_mode = blend_mode
        self.opacity = opacity
        self.image = image
        self.norm_image = None
        self.check_data()

    def check_data(self):
        if self.vis is None:
            self.normalization = None
            self.min = None
            self.max = None
            self.blend_mode = None
            self.opacity = None
            return
        if self.normalization.lower() != "value" and self.normalization.lower() != "perc" and \
                self.normalization is not None:
            raise Exception("RVT BlenderLayer check_data: normalization value incorrect!")
        if 0 > self.min > 100 and self.min is not None:
            raise Exception("RVT BlenderLayer check_data: min value incorrect [0-100]!")
        if 0 > self.max > 100 and self.max is not None:
            raise Exception("RVT BlenderLayer check_data: max value incorrect [0-100]!")
        if self.min > self.max:
            raise Exception("RVT BlenderLayer check_data: min bigger than max!")
        if self.blend_mode.lower() != "normal" and self.blend_mode.lower() != "multiply" and \
                self.blend_mode.lower() != "overlay" and self.blend_mode.lower() != "luminosity" and \
                self.blend_mode.lower() != "screen":
            raise Exception("RVT BlenderLayer check_data: blend_mode incorrect!")
        if 0 > self.opacity > 100:
            raise Exception("RVT BlenderLayer check_data: opacity incorrect [0-100]!")
        if self.vis is not None and self.image is None:
            raise Exception("RVT BlenderLayer check_data: Layer needs image!")


class BlenderLayers:
    layers = []
    default = rvt.default.DefaultValues()

    # create and add layer
    def create_layer(self, vis_method=None, normalization="value", minimum=None, maximum=None,
                     blend_mode="normal", opacity=100, image=None):
        layer = BlenderLayer(vis_method=vis_method, normalization=normalization, minimum=minimum, maximum=maximum,
                             blend_mode=blend_mode, opacity=opacity, image=image)
        self.layers.append(layer)

    def add_layer(self, layer : BlenderLayer):
        self.layers.append(layer)

    def normalize_images(self):
        nr_layers = sum(lyr.vis is not None for lyr in self.layers)
        nr_images = sum(lyr.image is not None for lyr in self.layers)

        if nr_layers != nr_images:
            raise Exception("RVT normalize_images_on_layers: layers and images number don't match")

        for i_img in range(nr_layers):
            visualization = self.layers[i_img].vis
            if visualization is None:
                continue
            image = self.layers[i_img].image
            min_norm = self.layers[i_img].min
            max_norm = self.layers[i_img].max
            normalization = self.layers[i_img].normalization

            # workaround for RGB images because they are on scale [0, 255] not [0, 1],
            # we use multiplier to get proper values
            if normalization.lower() == "value" and visualization.lower() == "hillshade":
                if np.nanmax(image) > 100.0 and image.shape[0] == 3:
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
                    warnings.warn("RVT normalize_images_on_layers: unexpected values! max > 1")
                if np.nanmin(norm_image) < 0:
                    warnings.warn("RVT normalize_images_on_layers: unexpected values! min < 0")

            # for slope and neg openness, invert scale
            # meaning high slopes will be black
            if visualization.lower() == "openness - negative" or visualization.lower() == "slope gradient":
                norm_image = 1 - norm_image

            self.layers[i_img].norm_image = norm_image

    def check_for_normalization(self):
        for lyr in self.layers:
            if lyr.norm_image is None and lyr.vis is not None:
                raise Exception("RVT BlenderLayers: normalize_images_on_layers before render!")

    def blend_normal(self, active, background):
        return active

    def blend_screen(self, active, background):
        return 1 - (1 - active) * (1 - background)

    def blend_multiply(self, active, background):
        return active * background

    def blend_overlay(self, active, background):
        idx1 = np.where(background > 0.5)
        idx2 = np.where(background <= 0.5)
        background[idx1[0], idx1[1]] = (
                    1 - (1 - 2 * (background[idx1[0], idx1[1]] - 0.5)) * (1 - active[idx1[0], idx1[1]]))
        background[idx2[0], idx2[1]] = ((2 * background[idx2[0], idx2[1]]) * active[idx2[0], idx2[1]])
        return background

    def equation_blend(self, blend_mode, active, background):
        if blend_mode.lower() == "screen":
            return self.blend_screen(active, background)
        elif blend_mode.lower() == "multiply":
            return self.blend_multiply(active, background)
        elif blend_mode.lower() == "overlay":
            return self.blend_overlay(active, background)

    def blend_luminosity(self, active, background, min_c=None, max_c=None):
        lum_active = lum(active)
        lum_background = lum(background)
        luminosity = lum_active - lum_background

        if background.shape[0] < 3:
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

    def blend_images(self, blend_mode, active, background, min_c=None, max_c=None):
        if blend_mode.lower() == "multiply" or blend_mode.lower() == "overlay" or blend_mode.lower() == "screen":
            return self.blend_multi_dim_images(blend_mode, active, background)
        elif blend_mode.lower() == "luminosity":
            return self.blend_luminosity(active, background, min_c, max_c)
        else:
            return self.blend_normal(active, background)

    def apply_opacity(self, active, background, opacity):
        if opacity > 1:
            opacity = opacity / 100
        return active * opacity + background * (1 - opacity)

    def blend_multi_dim_images(self, blend_mode, active, background):
        a_rgb = active.shape[0] == 3  # bool, is active rgb
        b_rgb = background.shape[0] == 3  # bool, is background rgb

        if a_rgb and b_rgb:
            blended_image = np.zeros(background.shape)
            for i in range(3):
                blended_image[i, :, :] = self.equation_blend(blend_mode, active[i, :, :], background[i, :, :])
        if a_rgb and not b_rgb:
            blended_image = np.zeros(active.shape)
            for i in range(3):
                blended_image[i, :, :] = self.equation_blend(blend_mode, active[i, :, :], background)
        if not a_rgb and b_rgb:
            blended_image = np.zeros(background.shape)
            for i in range(3):
                blended_image[i, :, :] = self.equation_blend(blend_mode, active, background[i, :, :])
        if not a_rgb and not b_rgb:
            blended_image = self.equation_blend(blend_mode, active, background)

        return blended_image

    def render_images(self, active, background, opacity):
        if np.nanmin(active) < 0 or np.nanmax(active) > 1.1:
            active = scale_0_to_1(active)
        if np.nanmin(background) < 0 or np.nanmax(background) > 1.1:
            background = scale_0_to_1(background)

        a_rgb = active.shape[0] == 3
        b_rgb = background.shape[0] == 3
        render_image = 0
        if a_rgb and b_rgb:
            render_image = np.zeros(background.shape)
            for i in range(3):
                render_image[i, :, :] = self.apply_opacity(active[i, :, :], background[i, :, :], opacity)
        if a_rgb and not b_rgb:
            render_image = np.zeros(active.shape)
            for i in range(3):
                render_image[i, :, :] = self.apply_opacity(active[i, :, :], background, opacity)
        if not a_rgb and b_rgb:
            render_image = np.zeros(background.shape)
            for i in range(3):
                render_image[i, :, :] = self.apply_opacity(active, background[i, :, :], opacity)
        if not a_rgb and not b_rgb:
            render_image = self.apply_opacity(active, background, opacity)
        return render_image

    # rendering across all layers - form last to first layer
    # each image is either: (1) image of visualization + blending or (2) original image (if vis == None for that layer)
    def render_all_images(self):
        self.check_for_normalization()

        rendered_image = []

        for i_img in range(len(self.layers) - 1, -1, -1):
            # if current layer visualization applied, skip
            visualization = self.layers[i_img].vis
            if visualization is None:
                continue

            # if current layer has visualization applied, but there has been no rendering
            # of images yet, than current layer will be the initial value of rendered_image
            if rendered_image == []:
                rendered_image = self.layers[i_img].norm_image
                continue
            else:
                active = self.layers[i_img].norm_image
                background = rendered_image
                blend_mode = self.layers[i_img].blend_mode
                opacity = self.layers[i_img].opacity
                if np.nanmin(active) < 0 or np.nanmax(active) > 1:
                    active = scale_0_to_1(active)
                if np.nanmin(background) < 0 or np.nanmax(background) > 1:
                    background = scale_0_to_1(background)
                top = self.blend_images(blend_mode, active, background)
                rendered_image = self.render_images(top, background, opacity)

                if np.nanmin(background) < 0 or np.nanmax(background > 1):
                    warnings.warn("RVT render_all_images: Rendered image scale distorted")

        return rendered_image

    def build_blender_layers_from_file(self, file_path):
        # Example file in dir settings: blender_file_example.txt
        dat = open(file_path, "r")
        for line in dat:
            line = line.strip()
            if line == "":
                continue
            if line[0] == "#":
                continue
            line_list = line.split(",")
            vis_method = None
            normalization = None
            minimum = None
            maximum = None
            blend_mode = None
            opacity = None
            image = None
            for element in line_list:
                element = element.strip()
                parameter_name = element.split(":")[0].strip()
                parameter_value = element.split(":")[1].strip()
                if parameter_name == "vis":
                    vis_method = parameter_value
                elif parameter_name == "norm":
                    normalization = parameter_value
                elif parameter_name == "min":
                    minimum = float(parameter_value)
                elif parameter_name == "max":
                    maximum = float(parameter_value)
                elif parameter_name == "blend_mode":
                    blend_mode = parameter_value
                elif parameter_name == "opacity":
                    opacity = float(parameter_value)
            # TODO: fill image (using rvt.default), first finish default
            #image
            layer = BlenderLayer(vis_method=vis_method, normalization=normalization, minimum=minimum, maximum=maximum,
                                 blend_mode=blend_mode, opacity=opacity, image=image)
            self.layers.append(layer)
        dat.close()


