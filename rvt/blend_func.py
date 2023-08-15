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
    Nejc Čož

Copyright:
    2010-2020 Research Centre of the Slovenian Academy of Sciences and Arts
    2016-2020 University of Ljubljana, Faculty of Civil and Geodetic Engineering
"""

# TODO: more testing, find and fix bugs if they exists

import warnings
from typing import Any, Optional, Dict, Tuple

import numpy as np
import numpy.typing as npt
from matplotlib.cm import get_cmap
from matplotlib.colors import LinearSegmentedColormap, Colormap


def gray_scale_to_color_ramp(
    gray_scale: npt.NDArray[Any],
    colormap: str,
    min_colormap_cut: Optional[float] = None,
    max_colormap_cut: Optional[float] = None,
    alpha: bool = False,
    output_8bit: bool = True,
) -> npt.NDArray[Any]:
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
    cm = get_cmap(colormap)

    # Truncate colormap if required
    if min_colormap_cut is not None or max_colormap_cut is not None:
        if min_colormap_cut is None:
            min_colormap_cut = 0.0
        if max_colormap_cut is None:
            max_colormap_cut = 1.0
        if (
            min_colormap_cut > 1
            or min_colormap_cut < 0
            or max_colormap_cut > 1
            or max_colormap_cut < 0
        ):
            raise Exception(
                "rvt.blend_func.gray_scale_to_color_ramp: min_colormap_cut and max_colormap_cut must be"
                " between 0 and 1!"
            )
        if min_colormap_cut >= max_colormap_cut:
            raise Exception(
                "rvt.blend_func.gray_scale_to_color_ramp: min_colormap_cut can't be smaller than"
                " max_colormap_cut!"
            )
        cm = truncate_colormap(
            cmap=cm, minval=min_colormap_cut, maxval=max_colormap_cut
        )

    # Compute normalized RGBA
    rgba_mtpl_out = cm(gray_scale)

    if output_8bit:
        nan_mask = np.isnan(gray_scale)
        rgba_mtpl_out[nan_mask] = 0  # Change nan to 0
        rgba_mtpl_out = np.uint8(
            rgba_mtpl_out * 255
        )  # 0-1 scale to 0-255 and change type to uint8

    # Move array axes to correct positions, i.e. (x, y, bands) to (bands, x, y)
    rgba_out = rgba_mtpl_out.transpose(2, 0, 1)

    # Discard 4th band if not using Alpha
    if not alpha:
        rgba_out = rgba_out[:3, ...]

    return rgba_out


def truncate_colormap(
    cmap: Colormap, minval: float = 0.0, maxval: float = 1.0, n: int = 100
) -> LinearSegmentedColormap:
    new_cmap = LinearSegmentedColormap.from_list(
        "trunc({n},{a:.2f},{b:.2f})".format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)),
    )
    return new_cmap


def normalize_lin(
    image: npt.NDArray[Any], minimum: float, maximum: float
) -> npt.NDArray[Any]:
    # linear cut off
    image[image > maximum] = maximum
    image[image < minimum] = minimum

    # stretch to 0.0 - 1.0 interval
    image = (image - minimum) / (maximum - minimum)
    image[image > 1] = 1
    image[image < 0] = 0
    return image.astype(np.float32)


def lin_cutoff_calc_from_perc(
    image: npt.NDArray[Any], minimum: float, maximum: float
) -> Dict[str, float]:
    """
    Minimum cutoff in percent, maximum cutoff in percent (0%-100%). Returns min and max values for linear
    stretch (cut-off).
    """
    if minimum < 0 or maximum < 0 or minimum > 100 or maximum > 100:
        raise Exception(
            "rvt.blend_func.lin_cutoff_calc_from_perc: minimum, maximum are percent and have to be in "
            "range 0-100!"
        )
    if minimum + maximum > 100:
        raise Exception(
            "rvt.blend_func.lin_cutoff_calc_from_perc: if minimum + maximum > 100% then there are no"
            " values left! You can't cutoff whole image!"
        )
    distribution = np.nanpercentile(a=image, q=np.array([minimum, 100 - maximum]))
    min_lin = distribution[0]
    max_lin = distribution[1]
    if min_lin == max_lin:
        min_lin = np.nanmin(image)
        max_lin = np.nanmax(image)
    return {"min_lin": min_lin, "max_lin": max_lin}


def normalize_perc(
    image: npt.NDArray[Any], minimum: float, maximum: float
) -> npt.NDArray[Any]:
    min_max_lin_dict = lin_cutoff_calc_from_perc(image, minimum, maximum)
    min_lin = min_max_lin_dict["min_lin"]
    max_lin = min_max_lin_dict["max_lin"]
    return normalize_lin(image, min_lin, max_lin)


def advanced_normalization(
    image: npt.NDArray[Any], minimum: float, maximum: float, normalization: str
) -> npt.NDArray[Any]:
    """Runs normalization based on the selected normalization type: value or percent."""

    # Preform checks if correct values were given
    if minimum == maximum and normalization == "value":
        raise Exception(
            "rvt.blend_func.advanced_normalization: If normalization == value, min and max cannot be the"
            " same!"
        )

    if minimum > maximum and normalization == "value":
        raise Exception(
            "rvt.blend_func.advanced_normalization: If normalization == value, max can't be smaller"
            " than min!"
        )

    # Select normalization type
    if normalization.lower() == "value":
        equ_image = normalize_lin(image=image, minimum=minimum, maximum=maximum)
    elif normalization.lower() == "perc":
        equ_image = normalize_perc(image=image, minimum=minimum, maximum=maximum)
    elif normalization is None:
        equ_image = image
    else:
        raise Exception(
            f"rvt.blend_func.advanced_normalization: Unknown normalization type: {normalization}"
        )

    return equ_image


def lum(img: npt.NDArray[Any]) -> npt.NDArray[Any]:
    if len(img.shape) == 3:
        r = img[0]
        g = img[1]
        b = img[2]
        lum_img = ((0.3 * r) + (0.59 * g) + (0.11 * b)).astype(np.float32)
    else:
        lum_img = img  # type: ignore

    return lum_img  # type: ignore


def matrix_eq_min_lt_zero(
    r: npt.NDArray[Any],
    idx_min_lt_zero: Tuple[npt.NDArray[np.signedinteger], ...],
    lum_c: npt.NDArray[Any],
    min_c: npt.NDArray[Any],
) -> npt.NDArray[Any]:
    r[idx_min_lt_zero] = lum_c[idx_min_lt_zero] + (
        ((r[idx_min_lt_zero] - lum_c[idx_min_lt_zero]) * lum_c[idx_min_lt_zero])
        / (lum_c[idx_min_lt_zero] - min_c[idx_min_lt_zero])
    )
    return r


def matrix_eq_max_gt_one(
    r: npt.NDArray[Any],
    idx_max_c_gt_one: Tuple[npt.NDArray[np.signedinteger], ...],
    lum_c: npt.NDArray[Any],
    max_c: npt.NDArray[Any],
) -> npt.NDArray[Any]:
    r[idx_max_c_gt_one] = lum_c[idx_max_c_gt_one] + (
        (
            (r[idx_max_c_gt_one] - lum_c[idx_max_c_gt_one])
            * (1.0 - lum_c[idx_max_c_gt_one])
        )
        / (max_c[idx_max_c_gt_one] - lum_c[idx_max_c_gt_one])
    )
    return r


def channel_min(
    r: npt.NDArray[Any], g: npt.NDArray[Any], b: npt.NDArray[Any]
) -> npt.NDArray[Any]:
    min_c = r * 1.0
    idx_min = np.where(g < min_c)
    min_c[idx_min] = g[idx_min]
    idx_min = np.where(b < min_c)
    min_c[idx_min] = b[idx_min]
    return min_c


def channel_max(
    r: npt.NDArray[Any], g: npt.NDArray[Any], b: npt.NDArray[Any]
) -> npt.NDArray[Any]:
    max_c = r * 1.0
    idx_max = np.where(g > max_c)
    max_c[idx_max] = g[idx_max]
    idx_max = np.where(b > max_c)
    max_c[idx_max] = b[idx_max]
    return max_c


def clip_color(
    c: npt.NDArray[Any],
    min_c: Optional[npt.NDArray[Any]] = None,
    max_c: Optional[npt.NDArray[Any]] = None,
) -> npt.NDArray[Any]:
    lum_c = lum(c)

    r = c[0].astype(np.float32)
    g = c[1].astype(np.float32)
    b = c[2].astype(np.float32)

    if min_c is None:
        min_c = channel_min(r, g, b)
    if max_c is None:
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


def blend_normal(
    active: npt.NDArray[Any], background: npt.NDArray[Any]
) -> npt.NDArray[Any]:
    return active


def blend_screen(
    active: npt.NDArray[Any], background: npt.NDArray[Any]
) -> npt.NDArray[Any]:
    return 1 - (1 - active) * (1 - background)


def blend_multiply(
    active: npt.NDArray[Any], background: npt.NDArray[Any]
) -> npt.NDArray[Any]:
    return active * background


def blend_overlay(
    active: npt.NDArray[Any], background: npt.NDArray[Any]
) -> npt.NDArray[Any]:
    idx1 = np.where(background > 0.5)
    idx2 = np.where(background <= 0.5)
    background[idx1] = 1 - (1 - 2 * (background[idx1] - 0.5)) * (1 - active[idx1])
    background[idx2] = (2 * background[idx2]) * active[idx2]
    return background


def blend_soft_light(
    active: npt.NDArray[Any], background: npt.NDArray[Any]
) -> npt.NDArray[Any]:
    # idx1 = np.where(active > 0.5)
    # idx2 = np.where(active <= 0.5)
    # background[idx1] = 1 - (1-background[idx1]) * (1-(active[idx1]-0.5))
    # background[idx2] = background[idx2] * (active[idx2]+0.5)
    idx1 = np.where(active < 0.5)
    idx2 = np.where(active >= 0.5)
    background[idx1] = 2 * background[idx1] * active[idx1] + background[idx1] ** 2 * (
        1.0 - 2 * active[idx1]
    )
    background[idx2] = 2 * background[idx2] * (1.0 - active[idx2]) + np.sqrt(
        background[idx2]
    ) * (2 * active[idx2] - 1.0)
    return background


def blend_luminosity(
    active: npt.NDArray[Any],
    background: npt.NDArray[Any],
    min_c: Optional[npt.NDArray[Any]] = None,
    max_c: Optional[npt.NDArray[Any]] = None,
) -> npt.NDArray[Any]:
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


def equation_blend(
    blend_mode: str, active: npt.NDArray[Any], background: npt.NDArray[Any]
) -> npt.NDArray[Any]:
    if blend_mode.lower() == "screen":
        return blend_screen(active, background)
    elif blend_mode.lower() == "multiply":
        return blend_multiply(active, background)
    elif blend_mode.lower() == "overlay":
        return blend_overlay(active, background)
    elif blend_mode.lower() == "soft_light":
        return blend_soft_light(active, background)
    else:
        raise ValueError(f"Invalid blend mode, {blend_mode=}!")


def blend_multi_dim_images(
    blend_mode: str, active: npt.NDArray[Any], background: npt.NDArray[Any]
) -> npt.NDArray[Any]:
    a_rgb = len(active.shape) == 3  # bool, is active rgb
    b_rgb = len(background.shape) == 3  # bool, is background rgb
    blended_image = None
    if a_rgb and b_rgb:
        blended_image = np.zeros(background.shape)
        for i in range(3):
            blended_image[i, :, :] = equation_blend(
                blend_mode, active[i, :, :], background[i, :, :]
            )
    if a_rgb and not b_rgb:
        blended_image = np.zeros(active.shape)
        for i in range(3):
            blended_image[i, :, :] = equation_blend(
                blend_mode, active[i, :, :], background
            )
    if not a_rgb and b_rgb:
        blended_image = np.zeros(background.shape)
        for i in range(3):
            blended_image[i, :, :] = equation_blend(
                blend_mode, active, background[i, :, :]
            )
    if not a_rgb and not b_rgb:
        blended_image = equation_blend(blend_mode, active, background)

    assert blended_image is not None

    return blended_image


def blend_images(
    blend_mode: str,
    active: npt.NDArray[Any],
    background: npt.NDArray[Any],
    min_c: Optional[npt.NDArray[Any]] = None,
    max_c: Optional[npt.NDArray[Any]] = None,
) -> npt.NDArray[Any]:
    if (
        blend_mode.lower() == "multiply"
        or blend_mode.lower() == "overlay"
        or blend_mode.lower() == "screen"
        or blend_mode.lower() == "soft_light"
    ):
        return blend_multi_dim_images(blend_mode, active, background)
    elif blend_mode.lower() == "luminosity":
        return blend_luminosity(active, background, min_c, max_c)
    else:
        return blend_normal(active, background)


def render_images(
    active: npt.NDArray[Any], background: npt.NDArray[Any], opacity: float
) -> npt.NDArray[Any]:

    # Both active and background image have to be between 0 and 1, scale if not
    if np.nanmin(active) < 0 or np.nanmax(active) > 1:
        active = scale_0_to_1(active)
    if np.nanmin(background) < 0 or np.nanmax(background) > 1:
        background = scale_0_to_1(background)

    # True if image has 3 bands (RGB), false if single band
    a_rgb = len(active.shape) == 3
    b_rgb = len(background.shape) == 3

    # Apply opacity
    if a_rgb and b_rgb:
        # Both images 3 bands
        render_image = np.zeros(background.shape)
        for i in range(3):
            render_image[i, :, :] = apply_opacity(
                active[i, :, :], background[i, :, :], opacity
            )
    elif a_rgb and not b_rgb:
        # Active image 3 bands
        render_image = np.zeros(active.shape)
        for i in range(3):
            render_image[i, :, :] = apply_opacity(active[i, :, :], background, opacity)
    elif not a_rgb and b_rgb:
        # Background image 3 bands
        render_image = np.zeros(background.shape)
        for i in range(3):
            render_image[i, :, :] = apply_opacity(active, background[i, :, :], opacity)
    else:
        render_image = apply_opacity(active, background, opacity)

    return render_image


def scale_within_0_and_1(numeric_value: npt.NDArray[Any]) -> npt.NDArray[Any]:
    if np.nanmin(numeric_value) >= 0 and np.nanmax(numeric_value) <= 1:
        return numeric_value

    # Create mask for NaN values
    # nan_mask = np.isnan(numeric_value)

    numeric_value[np.isnan(numeric_value)] = np.nanmin(
        numeric_value
    )  # nan change to nanmin

    actual_min = np.nanmin(numeric_value)
    norm_min_value = np.nanmax(np.array(0, actual_min))

    actual_max = np.nanmax(numeric_value)
    norm_max_value = np.nanmin(np.array(1, actual_max))

    # Do not scale values where max is between 1 and 255 if the max-min values difference is at least 30 and min >0
    # and numeric values are integer type
    if 255 >= actual_max > 1:
        if actual_max - actual_min > 30 and actual_min > 0:
            scaled = numeric_value / 255
            return scaled

    scaled = (numeric_value - norm_min_value) / (norm_max_value - norm_min_value)

    if np.nanmin(scaled) > -0.01:
        scaled[(0 > scaled) & (scaled > -0.01)] = 0

    # scaled[nan_mask] = np.nan

    return scaled


def scale_strict_0_to_1(numeric_value: npt.NDArray[Any]) -> npt.NDArray[Any]:
    if np.nanmin(numeric_value) == 0 and np.nanmax(numeric_value) == 1:
        return numeric_value

    numeric_value[np.isnan(numeric_value)] = 0  # nan change to 0

    min_value = np.nanmin(numeric_value)
    max_value = np.nanmax(numeric_value)

    scaled = (numeric_value - min_value) / (max_value - min_value)

    if np.nanmin(scaled) > -0.01:
        scaled[0 > scaled > -0.01] = 0

    return scaled


def scale_0_to_1(numeric_value: npt.NDArray[Any]) -> npt.NDArray[Any]:
    if 1 >= np.nanmax(numeric_value) > 0.9 and np.nanmin(numeric_value) == 0:
        return numeric_value
    elif np.nanmax(numeric_value) - np.nanmin(numeric_value) > 0.3:
        return scale_within_0_and_1(numeric_value)
    else:
        return scale_strict_0_to_1(numeric_value)


def apply_opacity(
    active: npt.NDArray[Any], background: npt.NDArray[Any], opacity: float
) -> npt.NDArray[Any]:
    if opacity > 1:
        opacity = opacity / 100
    return active * opacity + background * (1 - opacity)


def normalize_image(
    visualization: Optional[str],
    image: npt.NDArray[Any],
    min_norm: float,
    max_norm: float,
    normalization: str,
) -> Optional[npt.NDArray[Any]]:
    """Main function for normalization. Runs advanced normalization on the array and preforms special operations for
    some visualization types (e.g. invert scale for slope, scale for mhs, etc.).
    """
    if visualization is None:
        return None

    if normalization == "percent":
        normalization = "perc"

    norm_image = advanced_normalization(
        image=image, minimum=min_norm, maximum=max_norm, normalization=normalization
    )

    # Make sure it scales 0 to 1
    if np.nanmax(norm_image) > 1:
        if (
            visualization.lower() == "multiple directions hillshade"
            or visualization == "mhs"
        ):
            norm_image = scale_0_to_1(norm_image)
        else:
            norm_image = scale_0_to_1(norm_image)
            warnings.warn("rvt.blend_func.normalize_image: unexpected values! max > 1")

    if np.nanmin(norm_image) < 0:
        norm_image = scale_0_to_1(norm_image)
        warnings.warn("rvt.blend_func.normalize_image: unexpected values! min < 0")

    # For slope invert scale (high slopes will be black)
    if (
        visualization.lower() == "slope gradient"
        or visualization.lower() == "openness - negative"
        or visualization == "slp"
        or visualization == "neg_opns"
    ):
        norm_image = 1 - norm_image

    return norm_image


def cut_off_normalize(
    image: npt.NDArray[Any],
    mode: str,
    cutoff_min: Optional[float] = None,
    cutoff_max: Optional[float] = None,
    bool_norm: bool = True,
) -> npt.NDArray[Any]:
    """
    One band image cut-off or normalization or both. Image is 2D np.ndarray of raster, mode is perc or value
    (min and max units), min and max are minimum value to cutoff and maximum value to cutoff.
    (e.x. percent min=2 and max=3 -> cutoff lower 2% values and higher 3% values;
     e.x. value min=10 and max=60 -> cutoff bellow 10 and above 60, image values will be 10-60)
    """
    if cutoff_min is not None and cutoff_max is not None:
        if cutoff_min == cutoff_max and mode == "value":
            raise Exception(
                "rvt.blend_func.cut_off_normalize: If normalization == value, min and max cannot be the"
                " same!"
            )
        if cutoff_min > cutoff_max and mode == "value":
            raise Exception(
                "rvt.blend_func.cut_off_normalize: If normalization == value, max can't be smaller"
                " than min!"
            )

    cut_off_arr = image
    if cutoff_min is None and mode.lower() == "value":
        cutoff_min = np.amin(image)
    if cutoff_max is None and mode.lower() == "value":
        cutoff_max = np.amax(image)
    if cutoff_min is None and (mode.lower() == "perc" or mode.lower() == "percent"):
        cutoff_min = 0
    if cutoff_max is None and (mode.lower() == "perc" or mode.lower() == "percent"):
        cutoff_max = 0
    assert cutoff_min is not None
    assert cutoff_max is not None
    if bool_norm:
        if mode.lower() == "value":
            cut_off_arr = normalize_lin(cut_off_arr, cutoff_min, cutoff_max)
        elif mode.lower() == "perc" or mode.lower() == "percent":
            cut_off_arr = normalize_perc(cut_off_arr, cutoff_min, cutoff_max)
    else:
        if mode.lower() == "value":
            cut_off_arr[cut_off_arr > cutoff_max] = cutoff_max
            cut_off_arr[cut_off_arr < cutoff_min] = cutoff_min
        elif mode.lower() == "perc" or mode.lower() == "percent":
            min_max_value_dict = lin_cutoff_calc_from_perc(
                cut_off_arr, cutoff_min, cutoff_max
            )
            min_value = min_max_value_dict["min_lin"]
            max_value = min_max_value_dict["max_lin"]
            cut_off_arr[cut_off_arr > max_value] = max_value
            cut_off_arr[cut_off_arr < min_value] = min_value
    return cut_off_arr
