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

# python libraries
import numpy as np


def image_join_channels(r, g, b):
    if r.shape != g.shape or r.shape != b.shape or g.shape != b.shape:
        raise Exception("RVT image_join_channels: r, g, b must me same dimensions!")
    return np.array([r, g, b])


def blend_normal(active, background):
    return active


def blend_screen(active, background):
    return 1 - (1 - active) * (1 - background)


def blend_multiply(active, background):
    return active * background


def blend_overlay(active, background):
    idx1 = np.where(background > 0.5)
    idx2 = np.where(background <= 0.5)
    background[idx1[0], idx1[1]] = (1 - (1 - 2 * (background[idx1[0], idx1[1]] - 0.5)) * (1 - active[idx1[0], idx1[1]]))
    background[idx2[0], idx2[1]] = ((2 * background[idx2[0], idx2[1]]) * active[idx2[0], idx2[1]])
    return background


def equation_blend(blend_mode, active, background):
    if blend_mode.lower() == "screen":
        return blend_screen(active, background)
    elif blend_mode.lower() == "multiply":
        return blend_multiply(active, background)
    elif blend_mode.lower() == "overlay":
        return blend_overlay(active, background)


def blend_multi_dim_images(blend_mode, active, background):
    a_rgb = active.shape[0] == 3  # bool, is active rgb
    b_rgb = background.shape[0] == 3  # bool, is background rgb

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


def blend_luminosity(active, background, min_c=None, max_c=None):
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


def blend_images(blend_mode, active, background, min_c=None, max_c=None):
    if blend_mode.lower() == "multiply" or blend_mode.lower() == "overlay" or blend_mode.lower() == "screen":
        return blend_multi_dim_images(blend_mode, active, background)
    elif blend_mode.lower() == "luminosity":
        return blend_luminosity(active, background, min_c, max_c)
    else:
        return blend_normal(active, background)


def apply_opacity(active, background, opacity):
    if opacity > 1:
        opacity = opacity/100
    return active * opacity + background * (1 - opacity)


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


def render_images(active, background, opacity):
    if np.nanmin(active) < 0 or np.nanmax(active) > 1.1:
        active = scale_0_to_1(active)
    if np.nanmin(background) < 0 or np.nanmax(background) > 1.1:
        background = scale_0_to_1(background)

    a_rgb = active.shape[0] == 3
    b_rgb = background.shpae[0] == 3

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
    if not a_rgb and b_rgb:
        render_image = apply_opacity(active, background, opacity)

    return render_image


# rendering across all layers - form last to first layer
# each image is either: (1) image of visualization + blending or (2) original image (if vis == None for that layer)
def render_all_images(layers, images):
    rendered_image = []

    for i_img in range(len(layers)-1, -1, -1):
        # if current layer visualization applied, skip
        visualization = layers[i_img].vis
        if visualization is None:
            continue

        # if current layer has visualization applied, but there has been no rendering of images yet, than current layer
        # will be the initial value of rendered_image
        if not rendered_image:
            rendered_image = images[i_img]
            continue
        else:
            active = images[i_img]
            background = rendered_image
            blend_mode = layers[i_img].blend_mode
            opacity = layers[i_img].opacity

            if np.nanmin(active) < 0 or np.nanmax(active) > 1:
                active = scale_0_to_1(active)
            if np.nanmin(background) < 0 or np.nanmax(background) > 1:
                background = scale_0_to_1(background)

            top = blend_images(blend_mode, active, background)
            rendered_image = render_images(top, background, opacity)

            if np.nanmin(background) < 0 or np.nanmax(background > 1):
                raise Warning("RVT render_all_images: Rendered omage scale distorted")

        return rendered_image
