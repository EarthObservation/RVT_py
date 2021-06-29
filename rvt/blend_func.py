"""
Relief Visualization Toolbox – Visualization Functions

Contains core functions for blending.

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
import matplotlib as mpl
import matplotlib.cm
import matplotlib.colors
import numpy as np
import warnings


def gray_scale_to_color_ramp(gray_scale, colormap, min_colormap_cut=None, max_colormap_cut=None, alpha=False,
                             output_8bit=True):
    """
    Turns normalized gray scale np.array to rgba (np.array of 4 np.arrays r, g, b, a).

    Parameters
    ----------
    gray_scale : np.array (2D)
        Normalized gray_scale img as np.array (0-1)
    colormap : str
        Colormap form matplotlib (https://matplotlib.org/3.3.2/tutorials/colors/colormaps.html)
    min_colormap_cut : float
        What lower part of colormap to cut to select part of colormap.
        Valid values are between 0 and 1, if 0.2 it cuts off (deletes) 20% of lower colors in colormap.
        If None cut is not applied.
    max_colormap_cut : float
        What upper part of colormap to cut to select part of colormap.
        Valid values are between 0 and 1, if 0.8 it cuts off (deletes) 20% of upper colors in colormap.
        If None cut is not applied.
    alpha : bool
        If True outputs 4D array RGBA, if False outputs 3D array RGB
    output_8bit : bool
        If true output values will be int 0-255 instead of normalized values.
    Returns
    -------
    rgba_out : np.array (3D: red 0-255, green 0-255, blue 0-255)
            If alpha False: np.array (4D: red 0-255, green 0-255, blue 0-255, alpha 0-255)
    """
    cm = mpl.cm.get_cmap(colormap)
    if min_colormap_cut is not None or max_colormap_cut is not None:
        if min_colormap_cut is None:
            min_colormap_cut = 0.0
        if max_colormap_cut is None:
            max_colormap_cut = 1.0
        if min_colormap_cut > 1 or min_colormap_cut < 0 or max_colormap_cut > 1 or max_colormap_cut < 0:
            raise Exception("rvt.blend_funct.gray_scale_to_color_ramp: min_colormap_cut and max_colormap_cut must be"
                            " between 0 and 1!")
        if min_colormap_cut >= max_colormap_cut:
            raise Exception("rvt.blend_funct.gray_scale_to_color_ramp: min_colormap_cut can't be smaller than"
                            " max_colormap_cut!")
        cm = truncate_colormap(cmap=cm, minval=min_colormap_cut, maxval=max_colormap_cut)
    rgba_mtpl_out = cm(gray_scale)  # normalized rgb
    if output_8bit:
        rgba_mtpl_out = np.uint8(rgba_mtpl_out * 255)  # 0-1 scale to 0-255 and change type to uint8
    if not alpha:
        rgba_out = np.array([rgba_mtpl_out[:, :, 0], rgba_mtpl_out[:, :, 1], rgba_mtpl_out[:, :, 2]])
    else:
        rgba_out = np.array([rgba_mtpl_out[:, :, 0], rgba_mtpl_out[:, :, 1], rgba_mtpl_out[:, :, 2],
                             rgba_mtpl_out[:, :, 3]])
    return rgba_out


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def normalize_lin(image, minimum, maximum):
    # linear cut off
    image[image > maximum] = maximum
    image[image < minimum] = minimum

    # stretch to 0.0 - 1.0 interval
    image = (image - minimum) / (maximum - minimum)
    image[image > 1] = 1
    image[image < 0] = 0
    return np.float32(image)


def lin_cutoff_calc_from_perc(image, minimum, maximum):
    """Minimum cutoff in percent, maximum cutoff in percent (0%-100%). Returns min and max values for linear
    stretch (cut-off)."""
    if minimum < 0 or maximum < 0 or minimum > 100 or maximum > 100:
        raise Exception("rvt.blend_funct.lin_cutoff_calc_from_perc: minimum, maximum are percent and have to be in "
                        "range 0-100!")
    if minimum + maximum > 100:
        raise Exception("rvt.blend_funct.lin_cutoff_calc_from_perc: if minimum + maximum > 100% then there are no"
                        " values left! You can't cutoff whole image!")
    distribution = np.nanpercentile(a=image, q=np.array([minimum, 100 - maximum]))
    min_lin = distribution[0]
    max_lin = distribution[1]
    if min_lin == max_lin:
        min_lin = np.nanmin(image)
        max_lin = np.nanmax(image)
    return {"min_lin": min_lin, "max_lin": max_lin}


def normalize_perc(image, minimum, maximum):
    min_max_lin_dict = lin_cutoff_calc_from_perc(image, minimum, maximum)
    min_lin = min_max_lin_dict["min_lin"]
    max_lin = min_max_lin_dict["max_lin"]
    return normalize_lin(image, min_lin, max_lin)


def advanced_normalization(image, minimum, maximum, normalization):
    equ_image = image
    if minimum == maximum and normalization == "value":
        raise Exception("rvt.blend_func.advanced_normalization: If normalization == value, min and max cannot be the"
                        " same!")
    if minimum > maximum and normalization == "value":
        raise Exception("rvt.blend_func.advanced_normalization: If normalization == value, max can't be smaller"
                        " than min!")
    if normalization.lower() == "value":
        equ_image = normalize_lin(image=image, minimum=minimum, maximum=maximum)
    elif normalization.lower() == "perc":
        equ_image = normalize_perc(image=image, minimum=minimum, maximum=maximum)
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
        lum_img = np.float32((0.3 * r) + (0.59 * g) + (0.11 * b))
    else:
        lum_img = img

    return lum_img


def matrix_eq_min_lt_zero(r, idx_min_lt_zero, lum_c, min_c):
    r[idx_min_lt_zero] = lum_c[idx_min_lt_zero] + (((r[idx_min_lt_zero] - lum_c[idx_min_lt_zero]) *
                                                    lum_c[idx_min_lt_zero]) / (lum_c[idx_min_lt_zero] -
                                                                               min_c[idx_min_lt_zero]))
    return r


def matrix_eq_max_gt_one(r, idx_max_c_gt_one, lum_c, max_c):
    r[idx_max_c_gt_one] = lum_c[idx_max_c_gt_one] + (((r[idx_max_c_gt_one] - lum_c[idx_max_c_gt_one]) *
                                                      (1.0 - lum_c[idx_max_c_gt_one])) / (max_c[idx_max_c_gt_one]
                                                                                          - lum_c[idx_max_c_gt_one]))
    return r


def channel_min(r, g, b):
    min_c = r * 1.0
    idx_min = np.where(g < min_c)
    min_c[idx_min] = g[idx_min]
    idx_min = np.where(b < min_c)
    min_c[idx_min] = b[idx_min]
    return min_c


def channel_max(r, g, b):
    max_c = r * 1.0
    idx_max = np.where(g > max_c)
    max_c[idx_max] = g[idx_max]
    idx_max = np.where(b > max_c)
    max_c[idx_max] = b[idx_max]
    return max_c


def clip_color(c, min_c=None, max_c=None):
    lum_c = lum(c)

    r = np.float32(c[0])
    g = np.float32(c[1])
    b = np.float32(c[2])

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

    c_out = np.zeros(c.shape)
    c_out[0, :, :] = r
    c_out[1, :, :] = g
    c_out[2, :, :] = b
    return c_out


def blend_normal(active, background):
    return active


def blend_screen(active, background):
    return 1 - (1 - active) * (1 - background)


def blend_multiply(active, background):
    return active * background


def blend_overlay(active, background):
    idx1 = np.where(background > 0.5)
    idx2 = np.where(background <= 0.5)
    background[idx1] = (1 - (1 - 2 * (background[idx1] - 0.5)) * (1 - active[idx1]))
    background[idx2] = ((2 * background[idx2]) * active[idx2])
    return background


def blend_soft_light(active, background):
    idx1 = np.where(background >= 0.5)
    idx2 = np.where(background < 0.5)
    background[idx1] = (active[idx1] + (background[idx1] - 0.5) * (0.5 - (0.5 - active[idx1]) ** 2))
    background[idx2] = (active[idx2] - (0.5 - background[idx2]) * (0.5 - (0.5 - active[idx2]) ** 2))
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
    elif blend_mode.lower() == "soft_light":
        return blend_soft_light(active, background)


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
    if blend_mode.lower() == "multiply" or blend_mode.lower() == "overlay" or blend_mode.lower() == "screen" \
            or blend_mode.lower() == "soft_light":
        return blend_multi_dim_images(blend_mode, active, background)
    elif blend_mode.lower() == "luminosity":
        return blend_luminosity(active, background, min_c, max_c)
    else:
        return blend_normal(active, background)


def render_images(active, background, opacity):
    if np.nanmin(active) < 0 or np.nanmax(active) > 1:
        active = scale_0_to_1(active)
    if np.nanmin(background) < 0 or np.nanmax(background) > 1:
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
        scaled[(0 > scaled) & (scaled > -0.01)] = 0

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
    if normalization == "percent":
        normalization = "perc"

    norm_image = advanced_normalization(image=image, minimum=min_norm, maximum=max_norm, normalization=normalization)

    # make sure it scales 0 to 1
    if np.nanmax(norm_image) > 1:
        if visualization.lower() == "multiple directions hillshade":
            norm_image = scale_0_to_1(norm_image)
        else:
            norm_image = scale_0_to_1(norm_image)
            warnings.warn("rvt.blend.normalize_images_on_layers: unexpected values! max > 1")
        if np.nanmin(norm_image) < 0:
            warnings.warn("rvt.blend.normalize_images_on_layers: unexpected values! min < 0")

    # for slope, invert scale
    # meaning high slopes will be black
    if visualization.lower() == "slope gradient" or visualization.lower() == "openness - negative":
        norm_image = 1 - norm_image
    return norm_image


def cut_off_normalize(image, mode, min=None, max=None, bool_norm=True):
    """
    One band image cut-off or normalization or both. Image is 2D np.ndarray of raster, mode is perc or value
    (min and max units), min and max are minimum value to cutoff and maximum value to cutoff.
    (e.x. percent min=2 and max=3 -> cutoff lower 2% values and higher 3% values;
     e.x. value min=10 and max=60 -> cutoff bellow 10 and above 60, image values will be 10-60)
    """
    if min is not None and max is not None:
        if min == max and mode == "value":
            raise Exception("rvt.blend_func.cut_off_normalize: If normalization == value, min and max cannot be the"
                            " same!")
        if min > max and mode == "value":
            raise Exception("rvt.blend_func.cut_off_normalize: If normalization == value, max can't be smaller"
                            " than min!")

    cut_off_arr = image
    if min is None and mode.lower() == "value":
        min = np.amin(image)
    if max is None and mode.lower() == "value":
        max = np.amax(image)
    if min is None and (mode.lower() == "perc" or mode.lower() == "percent"):
        min = 0
    if max is None and (mode.lower() == "perc" or mode.lower() == "percent"):
        max = 0
    if bool_norm:
        if mode.lower() == "value":
            cut_off_arr = normalize_lin(cut_off_arr, min, max)
        elif mode.lower() == "perc" or mode.lower() == "percent":
            cut_off_arr = normalize_perc(cut_off_arr, min, max)
    else:
        if mode.lower() == "value":
            cut_off_arr[cut_off_arr > max] = max
            cut_off_arr[cut_off_arr < min] = min
        elif mode.lower() == "perc" or mode.lower() == "percent":
            min_max_value_dict = lin_cutoff_calc_from_perc(cut_off_arr, min, max)
            min_value = min_max_value_dict["min_lin"]
            max_value = min_max_value_dict["max_lin"]
            cut_off_arr[cut_off_arr > max_value] = max_value
            cut_off_arr[cut_off_arr < min_value] = min_value
    return cut_off_arr
