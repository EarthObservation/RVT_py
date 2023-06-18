"""
Relief Visualization Toolbox – Visualization Functions

Contains functions for computing the visualizations.

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
from typing import Dict, Optional, Tuple, Union, Iterable, List, Any

# python libraries
import numpy as np
import numpy.typing as npt
from scipy.interpolate import griddata, RectBivariateSpline
from scipy.ndimage.morphology import distance_transform_edt
from scipy.spatial import cKDTree


def byte_scale(
    data: npt.NDArray[Any],
    c_min: Optional[float] = None,
    c_max: Optional[float] = None,
    high: int = 255,
    low: int = 0,
    no_data: Optional[float] = None,
) -> npt.NDArray[Any]:
    """
    Remade old scipy function.
    Byte scales an array (image). Linear scale.

    Byte scaling means converting the input image to uint8 dtype and scaling
    the range to ``(low, high)`` (default 0-255).

    Parameters
    ----------
    data : numpy.ndarray
        Input data (visualization) as 2D or multi-D numpy array.
    c_min : int or float
        Scalar, Bias scaling of small values. Default is ``data.min()``.
    c_max : int or float
        Scalar, Bias scaling of large values. Default is ``data.max()``.
    high : int
        Scalar, Scale max value to `high`.  Default is 255.
    low : int
        Scalar, Scale min value to `low`.  Default is 0.
    no_data : int or float
        Value that represents no_data, it is changed to np.nan .

    Returns
    -------
    img_array : uint8 numpy.ndarray
        The byte-scaled array.
    """
    is_2d_arr = False
    data_bands = data
    if len(data.shape) == 2:
        is_2d_arr = True
        data_bands = np.array([data])

    c_min_orig = c_min
    c_max_orig = c_max

    byte_data_bands = []
    for i_band in data_bands:
        data = i_band
        c_min = c_min_orig
        c_max = c_max_orig
        if high < low:
            raise ValueError("`high` should be larger than `low`.")

        if no_data is not None:  # change no data to np.nan
            data[data == no_data] = np.nan

        if c_min is None:
            c_min = np.nanmin(data)
        if c_max is None:
            c_max = np.nanmax(data)

        c_scale = c_max - c_min
        if c_scale < 0:
            raise ValueError("`cmax` should be larger than `cmin`.")
        elif c_scale == 0:
            c_scale = 1

        if data.dtype == np.uint8:
            # TODO: the following line seems not good to  me - if cmin=0, then that pixel will get negative value
            byte_data = (
                (high + 1) * (data - c_min - 1) / (c_max - c_min)
            )  # copied from IDL BYTSCL
            byte_data[byte_data > high] = high
            byte_data[byte_data < 0] = 0
            byte_data[np.isnan(byte_data)] = 0  # change no_data to 0
            return np.cast[np.uint8](byte_data) + np.cast[np.uint8](low)

        # scale = float(high - low) / cscale  # old scipy fn
        # byte_data = (data * 1.0 - cmin) * scale + 0.4999  # old scipy fn

        byte_data = (
            (high + 0.9999) * (data - c_min) / (c_max - c_min)
        )  # copied from IDL BYTSCL
        byte_data[byte_data > high] = high
        byte_data[byte_data < 0] = 0
        byte_data[np.isnan(byte_data)] = 255  # change no_data to 255
        byte_data = np.cast[np.uint8](byte_data) + np.cast[np.uint8](low)
        byte_data_bands.append(byte_data)

    if is_2d_arr:  # if only one band
        return byte_data_bands[0]
    else:  # multiple bands
        return np.array(byte_data_bands)


def slope_aspect(
    dem: npt.NDArray[Any],
    resolution_x: float = 1,
    resolution_y: float = 1,
    output_units: str = "radian",
    ve_factor: float = 1,
    no_data: Optional[float] = None,
) -> Dict[str, npt.NDArray[Any]]:
    """
    Procedure can return terrain slope and aspect in radian units (default) or in alternative units (if specified).
    Available alternative units are 'degree' and 'percent'.
    Slope is defined as 0 for horizontal plane and pi/2 for vertical plane.
    Aspect is defined as geographic azimuth: clockwise increasing, 0 or 2pi for the North direction.
         0
     270    90
        180
    Currently applied finite difference method.

    Parameters
    ----------
    dem : numpy.ndarray
        Input digital elevation model as 2D numpy array.
    resolution_x : int
        DEM resolution in X direction.
    resolution_y : int
        DEM resolution in Y direction.
    output_units : str
        Output units, you can choose between: percent, degree, radian. Default value is radian.
    ve_factor : int or float
        Vertical exaggeration factor.
    no_data : int or float
        Value that represents no_data, all pixels with this value are changed to np.nan. Only has to be specified if
        a numerical value is used for nodata (e.g. -9999).

    Returns
    -------
    dict_out: dict
        Returns {"slope": slope_out, "aspect": aspect_out};
        slope_out, slope gradient : 2D numpy array (numpy.ndarray) of slope;
        aspect_out, aspect : 2D numpy array (numpy.ndarray) of aspect.
    """
    if dem.ndim != 2:
        raise Exception("rvt.visualization.slope_aspect: dem has to be 2D np.array!")
    if not (10000 >= ve_factor >= -10000):
        raise Exception(
            "rvt.visualization.slope_aspect: ve_factor must be between -10000 and 10000!"
        )
    if resolution_x < 0 or resolution_y < 0:
        raise Exception(
            "rvt.visualization.slope_aspect: resolution must be a positive number!"
        )

    # Make sure array has the correct dtype!
    dem = dem.astype(np.float32)

    # Change no_data to np.nan
    if no_data is not None:
        dem[dem == no_data] = np.nan

    # Save NaN mask
    nan_dem = np.isnan(dem)

    # Add 1 pixel edge padding
    dem = np.pad(array=dem, pad_width=1, mode="edge")

    # Vertical exaggeration
    dem = dem * ve_factor

    # Derivatives in X and Y direction
    dzdx = (
        (roll_fill_nans(dem, 1, axis=1) - roll_fill_nans(dem, -1, axis=1)) / 2
    ) / resolution_x
    dzdy = (
        (roll_fill_nans(dem, -1, axis=0) - roll_fill_nans(dem, 1, axis=0)) / 2
    ) / resolution_y
    tan_slope = np.sqrt(dzdx**2 + dzdy**2)

    # Compute slope
    if output_units == "percent":
        slope_out = tan_slope * 100
    elif output_units == "degree":
        slope_out = np.rad2deg(np.arctan(tan_slope))
    elif output_units == "radian":
        slope_out = np.arctan(tan_slope)
    else:
        raise Exception(
            "rvt.visualization.calculate_slope: Wrong function input 'output_units'!"
        )

    # Compute Aspect
    # aspect identifies the down slope direction of the maximum rate of change in value from each cell to its neighbors:
    #     0
    # 270    90
    #    180
    dzdy[
        dzdy == 0
    ] = 10e-9  # important for numeric stability - where dzdy is zero, make tangent to really high value
    aspect_out = np.arctan2(dzdx, dzdy)  # atan2 took care of the quadrants
    if output_units == "degree":
        aspect_out = np.rad2deg(aspect_out)

    # Remove padding
    aspect_out = aspect_out[1:-1, 1:-1]
    slope_out = slope_out[1:-1, 1:-1]

    # Apply NaN mask
    slope_out[nan_dem] = np.nan
    aspect_out[nan_dem] = np.nan

    return {"slope": slope_out, "aspect": aspect_out}


def roll_fill_nans(dem: npt.NDArray[Any], shift: int, axis: int) -> npt.NDArray[Any]:
    """
    Uses numpy.roll() function to roll array, then checks element-wise if new array has NaN value, but there was a
    numerical value in the source array, then use the original value instead of NaN. It is equivalent to edge padding.

    https://numpy.org/doc/stable/reference/generated/numpy.roll.html#numpy.roll
    """
    out = np.roll(dem, shift, axis=axis)

    # Fill NaNs with values from dem
    out[np.isnan(out) != np.isnan(dem)] = dem[np.isnan(out) != np.isnan(dem)]

    return out


def hillshade(
    dem: npt.NDArray[Any],
    resolution_x: float,
    resolution_y: float,
    sun_azimuth: float = 315,
    sun_elevation: float = 35,
    slope: Optional[npt.NDArray[Any]] = None,
    aspect: Optional[npt.NDArray[Any]] = None,
    ve_factor: float = 1,
    no_data: Optional[float] = None,
) -> npt.NDArray[Any]:
    """
    Compute hillshade.

    Parameters
    ----------
    dem : numpy.ndarray
        Input digital elevation model as 2D numpy array.
    resolution_x : int
        DEM resolution in X direction.
    resolution_y : int
        DEM resolution in Y direction.
    sun_azimuth : int or float
        Solar azimuth angle (clockwise from North) in degrees.
    sun_elevation : int or float
        Solar vertical angle (above the horizon) in degrees.
    slope : numpy.ndarray
        Slope arr in radians if you don't input it, it is calculated.
    aspect : numpy.ndarray
        Aspect arr in radians if you don't input it, it is calculated.
    ve_factor : int or float
        Vertical exaggeration factor.
    no_data : int or float
        Value that represents no_data, all pixels with this value are changed to np.nan.

    Returns
    -------
    hillshade_out : numpy.ndarray
        Result hillshade 2D numpy array.
    """
    if dem.ndim != 2:
        raise Exception("rvt.visualization.hillshade: dem has to be 2D np.array!")
    if not (10000 >= ve_factor >= -10000):
        raise Exception(
            "rvt.visualization.hillshade: ve_factor must be between -10000 and 10000!"
        )
    if sun_azimuth > 360 or sun_elevation > 90 or sun_azimuth < 0 or sun_elevation < 0:
        raise Exception(
            "rvt.visualization.hillshade: sun_azimuth must be [0-360] and sun_elevation [0-90]!"
        )
    if resolution_x < 0 or resolution_y < 0:
        raise Exception(
            "rvt.visualization.hillshade: resolution must be a positive number!"
        )

    # change no_data to np.nan
    if no_data is not None:
        dem[dem == no_data] = np.nan

    dem = dem.astype(np.float32)
    # add 1 pixel edge padding
    dem = np.pad(array=dem, pad_width=1, mode="edge")
    dem = dem * ve_factor

    # Convert solar position (degrees) to radians
    sun_azimuth_rad = np.deg2rad(sun_azimuth)
    sun_elevation_rad = np.deg2rad(sun_elevation)

    # Convert to solar zenith angle
    sun_zenith_rad = np.pi / 2 - sun_elevation_rad

    # are slope and aspect already calculated and presented
    if slope is None or aspect is None:
        # calculates slope and aspect
        dict_slp_asp = slope_aspect(
            dem=dem,
            resolution_x=resolution_x,
            resolution_y=resolution_y,
            output_units="radian",
        )
        slope = dict_slp_asp["slope"]
        aspect = dict_slp_asp["aspect"]

    # Compute solar incidence angle, hillshading
    hillshade_out = np.cos(sun_zenith_rad) * np.cos(slope) + np.sin(
        sun_zenith_rad
    ) * np.sin(slope) * np.cos(aspect - sun_azimuth_rad)

    hillshade_out[hillshade_out < 0] = 0  # set all negative to 0

    # remove padding
    hillshade_out = hillshade_out[1:-1, 1:-1]

    return hillshade_out


def multi_hillshade(
    dem: npt.NDArray[Any],
    resolution_x: float,
    resolution_y: float,
    nr_directions: int = 16,
    sun_elevation: float = 35,
    slope: Optional[npt.NDArray[Any]] = None,
    aspect: Optional[npt.NDArray[Any]] = None,
    ve_factor: float = 1,
    no_data: Optional[float] = None,
) -> npt.NDArray[Any]:
    """
    Calculates hillshades from multiple directions.

    Parameters
    ----------
    dem : numpy.ndarray
        Input digital elevation model as 2D numpy array.
    resolution_x : int
        DEM resolution in X direction.
    resolution_y : int
        DEM resolution in Y direction.
    nr_directions : int
        Number of solar azimuth angles (clockwise from North).
    sun_elevation : int or float
        Solar vertical angle (above the horizon) in degrees.
    slope : numpy.ndarray
        Slope in radians if you don't input it, it is calculated.
    aspect : numpy.ndarray
        Aspect in radians if you don't input it, it is calculated.
    ve_factor : int or float
        Vertical exaggeration factor.
    no_data : int or float
        Value that represents no_data, all pixels with this value are changed to np.nan .

    Returns
    -------
    multi_hillshade_out : numpy.ndarray
        Result multiple direction hillshade multidimensional (nr_directions=dimensions) numpy array.
    """
    if dem.ndim != 2:
        raise Exception("rvt.visualization.multi_hillshade: dem has to be 2D np.array!")
    if sun_elevation > 90 or sun_elevation < 0:
        raise Exception(
            "rvt.visualization.multi_hillshade: sun_elevation must be [0-90]!"
        )
    if resolution_x < 0 or resolution_y < 0:
        raise Exception(
            "rvt.visualization.multi_hillshade: resolution must be a positive number!"
        )
    if nr_directions < 1:
        raise Exception(
            "rvt.visualization.multi_hillshade: nr_directions must be a positive number!"
        )
    if not (10000 >= ve_factor >= -10000):
        raise Exception(
            "rvt.visualization.multi_hillshade: ve_factor must be between -10000 and 10000!"
        )

    # change no_data to np.nan
    if no_data is not None:
        dem[dem == no_data] = np.nan

    dem = dem.astype(np.float32)
    dem = dem * ve_factor

    # calculates slope and aspect if they are not added
    if (
        slope is None or aspect is None
    ):  # slope and aspect are the same, so we have to calculate it once
        dict_slp_asp = slope_aspect(
            dem=dem,
            resolution_x=resolution_x,
            resolution_y=resolution_y,
            output_units="radian",
        )
        slope = dict_slp_asp["slope"]
        aspect = dict_slp_asp["aspect"]

    hillshades_arr_list = []  # list of all hillshades in different directions
    for i_direction in range(nr_directions):
        sun_azimuth = (360 / nr_directions) * i_direction
        hillshading = hillshade(
            dem=dem,
            resolution_x=resolution_x,
            resolution_y=resolution_y,
            sun_elevation=sun_elevation,
            sun_azimuth=sun_azimuth,
            slope=slope,
            aspect=aspect,
        )
        hillshades_arr_list.append(hillshading)
    multi_hillshade_out = np.asarray(hillshades_arr_list)

    return multi_hillshade_out


def mean_filter(dem: npt.NDArray[Any], kernel_radius: int) -> npt.NDArray[Any]:
    """Applies mean filter (low pass filter) on DEM. Kernel radius is in pixels. Kernel size is 2 * kernel_radius + 1.
    It uses matrix shifting (roll) instead of convolutional approach (works faster).
    It returns mean filtered dem as numpy.ndarray (2D numpy array)."""
    radius_cell = int(kernel_radius)

    if kernel_radius == 0:
        return dem

    # store nans
    idx_nan_dem = np.isnan(dem)

    # mean filter
    dem_pad = np.pad(dem, (radius_cell + 1, radius_cell), mode="edge")
    # store nans
    idx_nan_dem_pad = np.isnan(dem_pad)
    # change nan to 0
    dem_pad[idx_nan_dem_pad] = 0

    # kernel nr pixel integral image
    dem_i_nr_pixels = np.ones(dem_pad.shape)
    dem_i_nr_pixels[idx_nan_dem_pad] = 0
    dem_i_nr_pixels = integral_image(dem_i_nr_pixels, np.int64)

    dem_i1 = integral_image(dem_pad)

    kernel_nr_pix_arr = (
        np.roll(dem_i_nr_pixels, (radius_cell, radius_cell), axis=(0, 1))
        + np.roll(dem_i_nr_pixels, (-radius_cell - 1, -radius_cell - 1), axis=(0, 1))
        - np.roll(dem_i_nr_pixels, (-radius_cell - 1, radius_cell), axis=(0, 1))
        - np.roll(dem_i_nr_pixels, (radius_cell, -radius_cell - 1), axis=(0, 1))
    )
    mean_out = (
        np.roll(dem_i1, (radius_cell, radius_cell), axis=(0, 1))
        + np.roll(dem_i1, (-radius_cell - 1, -radius_cell - 1), axis=(0, 1))
        - np.roll(dem_i1, (-radius_cell - 1, radius_cell), axis=(0, 1))
        - np.roll(dem_i1, (radius_cell, -radius_cell - 1), axis=(0, 1))
    )
    mean_out = mean_out / kernel_nr_pix_arr
    mean_out = mean_out.astype(np.float32)
    mean_out = mean_out[
        radius_cell : -(radius_cell + 1), radius_cell : -(radius_cell + 1)
    ]  # remove padding
    # nan back to nan
    mean_out[idx_nan_dem] = np.nan

    return mean_out


def slrm(
    dem: npt.NDArray[Any],
    radius_cell: int = 20,
    ve_factor: float = 1,
    no_data: Optional[float] = None,
) -> npt.NDArray[Any]:
    """
    Calculates Simple local relief model.

    Parameters
    ----------
    dem : numpy.ndarray
        Input digital elevation model as 2D numpy array.
    radius_cell : int
        Radius for trend assessment in pixels.
    ve_factor : int or float
        Vertical exaggeration factor.
    no_data : int or float
        Value that represents no_data, all pixels with this value are changed to np.nan .

    Returns
    -------
    slrm_out : numpy.ndarray
        Simple local relief model 2D numpy array.
    """
    if dem.ndim != 2:
        raise Exception("rvt.visualization.slrm: dem has to be 2D np.array!")
    if radius_cell < 10 or radius_cell > 50:
        raise Exception(
            "rvt.visualization.slrm: Radius for trend assessment needs to be in interval 10-50 pixels!"
        )
    if not (10000 >= ve_factor >= -10000):
        raise Exception(
            "rvt.visualization.slrm: ve_factor must be between -10000 and 10000!"
        )

    # change no_data to np.nan
    if no_data is not None:
        dem[dem == no_data] = np.nan

    dem = dem.astype(np.float32)
    dem = dem * ve_factor

    # mean filter
    dem_mean_filter = mean_filter(dem=dem, kernel_radius=radius_cell)
    slrm_out = dem - dem_mean_filter

    return slrm_out


def horizon_shift_vector(
    num_directions: int = 16, radius_pixels: int = 10, min_radius: int = 1
) -> Dict[float, Dict[str, Union[List[float]]]]:
    """
    Calculates Sky-View determination movements.

    Parameters
    ----------
    num_directions : int
        Number of directions as input.
    radius_pixels : int
        Radius to consider in pixels (not in meters).
    min_radius : int
        Radius to start searching for horizon in pixels (not in meters).

    Returns
    -------
    shift : dict
        Dict with keys corresponding to the directions of search azimuths rounded to 1 decimal number
            - for each key, a subdict contains a key "shift":
                values for this key is a list of tuples prepared for np.roll - shift along lines and columns
            - the second key is "distance":
                values for this key is a list of search radius used for the computation of the elevation angle
    """

    # Initialize the output dict
    shift = {}

    # Generate angles and corresponding normal shifts in X (columns)
    # and Y (lines) direction
    angles = (2 * np.pi / num_directions) * np.arange(num_directions)
    x = np.cos(angles)
    y = np.sin(angles)
    angles = np.round(np.degrees(angles), decimals=1)

    # Generate a range of radius values in pixels.
    # Make it finer for the selected scaling.
    # By adding the last constant we make sure that we do not start with
    # point (0,0).
    scale = 3.0
    radii = np.arange((radius_pixels - min_radius) * scale + 1) / scale + min_radius

    # For each direction compute all possible horizon point position
    # and round them to integers
    for i in range(num_directions):
        x_int = np.round(x[i] * radii, decimals=0)
        y_int = np.round(y[i] * radii, decimals=0)
        # consider only the minimal number of points
        # use the trick with set and complex number as the input
        coord_complex = set(x_int + 1j * y_int)
        # to sort proportional with increasing radius,
        # set has to be converted to numpy array
        shift_pairs = np.array([(k.real, k.imag) for k in coord_complex]).astype(int)
        distance = np.sqrt(np.sum(shift_pairs**2, axis=1))
        sort_index = np.argsort(distance)
        # write for each direction shifts and corresponding distances
        shift[angles[i]] = {
            "shift": [(k[0], k[1]) for k in shift_pairs[sort_index]],
            "distance": distance[sort_index],
        }

    return shift


def sky_view_factor_compute(
    height_arr: npt.NDArray[Any],
    radius_max: int = 10,
    radius_min: int = 1,
    num_directions: int = 16,
    compute_svf: bool = True,
    compute_opns: bool = False,
    compute_asvf: bool = False,
    a_main_direction: float = 315.0,
    a_poly_level: int = 4,
    a_min_weight: float = 0.4,
) -> Dict[str, npt.NDArray[Any]]:
    """
    Calculates horizon based visualizations: Sky-view factor, Anisotropic SVF and Openness.

    SVF processing is using search radius, that looks at values beyond the edge of an array. Consider using a buffered
    array as an input, with the buffer size equal to the radius_max.
    To prevent erosion of the edge, function applies mirrored padding in all four directions, however, this means that
    edge values are "averaged over half of the hemisphere". Similarly, the edges of the dataset (i.e. areas with NaN
    values), will be considered as fully open (SFV angle 0, Openness angle -90).

    Input array should use np.nan as nodata value.

    Parameters
    ----------
    height_arr : numpy.ndarray
        Elevation (DEM) as 2D numpy array.
    radius_max : int
        Maximal search radius in pixels/cells (not in meters).
    radius_min : int
        Minimal search radius in pixels/cells (not in meters), for noise reduction.
    num_directions : int
        Number of directions as input.
    compute_svf : bool
        If true it computes and outputs svf.
    compute_asvf : bool
        If true it computes and outputs asvf.
    compute_opns : bool
        If true it computes and outputs opns.
    a_main_direction : int or float
        Main direction of anisotropy.
    a_poly_level : int
        Level of polynomial that determines the anisotropy.
    a_min_weight : float
        Weight to consider anisotropy:
                 0 - low anisotropy,
                 1 - high  anisotropy (no illumination from the direction opposite the main direction)

    Returns
    -------
    dict_out : dictionary
        Return {"svf": svf_out, "asvf": asvf_out, "opns": opns_out};
        svf_out, skyview factor : 2D numpy array (numpy.ndarray) of skyview factor;
        asvf_out, anisotropic skyview factor : 2D numpy array (numpy.ndarray) of anisotropic skyview factor;
        opns_out, openness : 2D numpy array (numpy.ndarray) openness (elevation angle of horizon).
    """

    # Pad the array for the radius_max on all 4 sides
    height = np.pad(height_arr, pad_width=radius_max, mode="reflect")  # type: ignore

    # Compute the vector of movement and corresponding distances
    move = horizon_shift_vector(
        num_directions=num_directions, radius_pixels=radius_max, min_radius=radius_min
    )

    # Initiate the output for SVF
    if compute_svf:
        svf_out = (
            height * 0
        )  # Multiply with 0 instead of using np.zeros to preserve nodata
    else:
        svf_out = None

    # Initiate the output for azimuth dependent SVF
    if compute_asvf:
        asvf_out = (
            height * 0
        )  # Multiply with 0 instead of using np.zeros to preserve nodata
        w_m = a_min_weight
        w_a = np.deg2rad(a_main_direction)
        weight = np.arange(num_directions) * (2 * np.pi / num_directions)
        weight = (1 - w_m) * (np.cos((weight - w_a) / 2)) ** a_poly_level + w_m
    else:
        asvf_out = None
        weight = None

    # Initiate the output for Openness
    if compute_opns:
        opns_out = (
            height * 0
        )  # Multiply with 0 instead of using np.zeros to preserve nodata
    else:
        opns_out = None

    # Search for horizon in each direction...
    for i_dir, direction in enumerate(move):
        # Reset maximum at each iteration (i.e. at the start of new direction),
        # smallest possible elevation angle is -1000 rad (i.e. -90 deg)
        max_slope = np.zeros(height.shape, dtype=np.float32) - 1000

        # ... and for each search radius
        for i_rad, radius in enumerate(move[direction]["distance"]):
            # Get shift index from move dictionary
            shift_indx = move[direction]["shift"][i_rad]
            # Estimate the slope
            _ = (np.roll(height, shift=shift_indx, axis=(0, 1)) - height) / radius  # type: ignore
            # Compare to the previous max slope and keep the largest values (element wise). Use np.fmax to prevent NaN
            # values contaminating the edge of the image (if one of the elements is NaN, pick non-NaN element)
            max_slope = np.fmax(max_slope, _)

        # Convert to angle in radians and compute directional output
        max_slope = np.arctan(max_slope)

        # Sum max angle for all directions
        if compute_svf:
            assert svf_out is not None
            # For SVF minimum possible angle is 0 (hemisphere), use np.fmax() to change NaNs to 0
            svf_out = svf_out + (1 - np.sin(np.fmax(max_slope, 0)))
        if compute_asvf:
            assert asvf_out is not None
            assert weight is not None
            # For SVF minimum possible angle is 0 (hemisphere), use np.fmax() to change NaNs to 0
            asvf_out = asvf_out + (1 - np.sin(np.fmax(max_slope, 0))) * weight[i_dir]
        if compute_opns:
            assert opns_out is not None
            # For Openness taking the entire sphere
            opns_out = opns_out + max_slope

    # Cut to original extent and average the directional output over all directions
    if compute_svf:
        assert svf_out is not None
        svf_out = (
            svf_out[radius_max:-radius_max, radius_max:-radius_max] / num_directions  # type: ignore
        )
    if compute_asvf:
        assert asvf_out is not None
        asvf_out = asvf_out[radius_max:-radius_max, radius_max:-radius_max] / np.sum(
            weight
        )
    if compute_opns:
        assert opns_out is not None
        opns_out = np.rad2deg(
            0.5 * np.pi
            - (
                opns_out[radius_max:-radius_max, radius_max:-radius_max]
                / num_directions
            )
        )

    # Return results within dict
    dict_svf_asvf_opns = {"svf": svf_out, "asvf": asvf_out, "opns": opns_out}
    dict_svf_asvf_opns = {
        k: v for k, v in dict_svf_asvf_opns.items() if v is not None
    }  # filter out none

    return dict_svf_asvf_opns  # type: ignore


def sky_view_factor(
    dem: npt.NDArray[Any],
    resolution: float,
    compute_svf: bool = True,
    compute_opns: bool = False,
    compute_asvf: bool = False,
    svf_n_dir: int = 16,
    svf_r_max: int = 10,
    svf_noise: int = 0,
    asvf_dir: float = 315,
    asvf_level: int = 1,
    ve_factor: float = 1,
    no_data: Optional[float] = None,
) -> Dict[str, npt.NDArray[Any]]:
    """
    Prepare the data, call sky_view_factor_compute, reformat and return back 2D arrays.

    Parameters
    ----------
    dem : numpy.ndarray
        Input digital elevation model as 2D numpy array.
    compute_svf : bool
        Compute SVF (True) or not (False).
    compute_opns : bool
        Compute OPENNESS (True) or not (False).
    resolution : float
        Pixel resolution.
    svf_n_dir : int
        Number of directions.
    svf_r_max : int
        Maximal search radius in pixels.
    svf_noise : int
        The level of noise remove (0-don't remove, 1-low, 2-med, 3-high).
    compute_asvf : bool
        Compute anisotropic SVF (True) or not (False).
    asvf_level : int
        Level of anisotropy, 1-low, 2-high.
    asvf_dir : int or float
        Direction of anisotropy.
    ve_factor : int or float
        Vertical exaggeration factor.
    no_data : int or float
        Value that represents no_data, all pixels with this value are changed to np.nan. Use this parameter when nodata
        is not np.nan.

    Returns
    -------
    dict_out : dictionary
        Return {"svf": svf_out, "asvf": asvf_out, "opns": opns_out};
        svf_out, skyview factor : 2D numpy array (numpy.ndarray) of skyview factor;
        asvf_out, anisotropic skyview factor : 2D numpy array (numpy.ndarray) of anisotropic skyview factor;
        opns_out, openness : 2D numpy array (numpy.ndarray) openness (elevation angle of horizon).
    """

    # Checks for input parameters
    if dem.ndim != 2:
        raise Exception("rvt.visualization.sky_view_factor: dem has to be 2D np.array!")
    if not (10000 >= ve_factor >= -10000):
        raise Exception(
            "rvt.visualization.sky_view_factor: ve_factor must be between -10000 and 10000!"
        )
    if svf_noise != 0 and svf_noise != 1 and svf_noise != 2 and svf_noise != 3:
        raise Exception(
            "rvt.visualization.sky_view_factor: svf_noise must be one of the following"
            "values (0-don't remove, 1-low, 2-med, 3-high)!"
        )
    if asvf_level != 1 and asvf_level != 2:
        raise Exception(
            "rvt.visualization.sky_view_factor: asvf_leve must be one of the following"
            "values (1-low, 2-high)!"
        )
    if not compute_svf and not compute_asvf and not compute_opns:
        raise Exception("rvt.visualization.sky_view_factor: All computes are false!")
    if resolution < 0:
        raise Exception(
            "rvt.visualization.sky_view_factor: resolution must be a positive number!"
        )

    # Make sure array has the correct dtype!
    dem = dem.astype(np.float32)

    # CONSTANTS
    # Level of polynomial that determines the anisotropy, selected with asvf_level (1 - low, 2 - high)
    sc_asvf_pol = [4, 8]
    sc_asvf_min = [0.4, 0.1]
    # The portion (percent) of the maximal search radius to ignore in horizon estimation; for each noise level,
    # selected with svf_noise (0-3)
    sc_svf_r_min = [0.0, 10.0, 20.0, 40.0]

    # Before doing anything to the array, make sure all NODATA values are set to np.nan
    if no_data is not None:
        dem[dem == no_data] = np.nan
    # Save NaN mask (processing may change NaNs to arbitrary values)
    nan_mask = np.isnan(dem)

    # Vertical exaggeration
    dem = dem * ve_factor
    # Pixel size (adjust elevation to correctly calculate the vertical elevation angle, calculation thinks 1px == 1m)
    dem = dem / resolution

    # Minimal search radius depends on the noise level, it has to be an integer not smaller than 1
    svf_r_min = max(np.round(svf_r_max * sc_svf_r_min[svf_noise] * 0.01, decimals=0), 1)

    # Set anisotropy parameters
    poly_level = sc_asvf_pol[asvf_level - 1]
    min_weight = sc_asvf_min[asvf_level - 1]

    # Main routine for SVF processing
    dict_svf_asvf_opns = sky_view_factor_compute(
        height_arr=dem,
        radius_max=svf_r_max,
        radius_min=svf_r_min,
        num_directions=svf_n_dir,
        compute_svf=compute_svf,
        compute_opns=compute_opns,
        compute_asvf=compute_asvf,
        a_main_direction=asvf_dir,
        a_poly_level=poly_level,
        a_min_weight=min_weight,
    )

    # Apply NaN mask to outputs
    for item in dict_svf_asvf_opns.values():
        item[nan_mask] = np.nan

    return dict_svf_asvf_opns


def local_dominance(
    dem: npt.NDArray[Any],
    min_rad: int = 10,
    max_rad: int = 20,
    rad_inc: int = 1,
    angular_res: int = 15,
    observer_height: float = 1.7,
    ve_factor: float = 1,
    no_data: Optional[float] = None,
) -> npt.NDArray[Any]:
    """
    Compute Local Dominance dem visualization.
    Adapted from original version that is part of the Lidar Visualisation Toolbox LiVT developed by Ralf Hesse.

    Parameters
    ----------
    dem : numpy.ndarray
        Input digital elevation model as 2D numpy array.
    min_rad : int
        Minimum radial distance (in pixels) at which the algorithm starts with visualization computation.
    max_rad : int
        Maximum radial distance (in pixels) at which the algorithm ends with visualization computation.
    rad_inc : int
        Radial distance steps in pixels.
    angular_res : int
        Angular step for determination of number of angular directions.
    observer_height : int or float
        Height at which we observe the terrain.
    ve_factor : int or float
        Vertical exaggeration factor.
    no_data : int or float
        Value that represents no_data, all pixels with this value are changed to np.nan .

    Returns
    -------
    local_dom_out : numpy.ndarray
        2D numpy array of local dominance
    """
    if dem.ndim != 2:
        raise Exception("rvt.visualization.local_dominance: dem has to be 2D np.array!")
    if not (10000 >= ve_factor >= -10000):
        raise Exception(
            "rvt.visualization.local_dominance: ve_factor must be between -10000 and 10000!"
        )

    # change no_data to np.nan
    if no_data is not None:
        dem[dem == no_data] = np.nan

    dem = dem.astype(np.float32)
    # add max_rad pixel edge padding
    pad_width = max_rad
    dem = np.pad(array=dem, pad_width=pad_width, mode="edge")
    dem = dem * ve_factor

    # create a vector with possible distances
    n_dist = int((max_rad - min_rad) / rad_inc + 1)
    distances = np.arange(n_dist * rad_inc, step=rad_inc) + min_rad  # type: ignore
    # create vector with possible angles
    n_ang = int(359 / angular_res + 1)
    angles = np.arange(n_ang * angular_res, step=angular_res)  # type: ignore
    # determine total area within radius range
    norma = np.sum((observer_height / distances) * (2 * distances + rad_inc)) * n_ang

    # image shifts
    n_shifts = distances.size * angles.size
    x_t = (np.outer(np.cos(np.deg2rad(angles)), distances)).reshape(n_shifts)
    y_t = (np.outer(np.sin(np.deg2rad(angles)), distances)).reshape(n_shifts)
    distances = (np.outer(np.ones(n_ang), distances)).reshape(n_shifts)
    dist_factor = 2 * distances + rad_inc

    local_dom_out = dem * 0
    for i_s in range(n_shifts):
        dem_moved = np.roll(dem, int(round(y_t[i_s])), axis=0)
        dem_moved = np.roll(dem_moved, int(round(x_t[i_s])), axis=1)
        idx_lower = np.where((dem + observer_height) > dem_moved)
        if idx_lower[0].size > 0:
            local_dom_out[idx_lower[0], idx_lower[1]] = (
                local_dom_out[idx_lower[0], idx_lower[1]]
                + (
                    dem[idx_lower[0], idx_lower[1]]
                    + observer_height
                    - dem_moved[idx_lower[0], idx_lower[1]]
                )
                / distances[i_s]
                * dist_factor[i_s]
            )
    local_dom_out = local_dom_out / norma

    # Remove padding
    local_dom_out = local_dom_out[pad_width:-pad_width, pad_width:-pad_width]

    return local_dom_out


def horizon_generate_coarse_dem(
    dem_fine: npt.NDArray[Any],
    pyramid_scale: int,
    conv_from: int,
    conv_to: int,
    max_radius: int,
) -> npt.NDArray[Any]:
    # first reduce the size for the edge required for horizon search
    dem_fine = dem_fine[max_radius:-max_radius, max_radius:-max_radius]

    # get and adjust the array sizes
    in_shape = dem_fine.shape
    n_col_fine = in_shape[1]
    n_lin_fine = in_shape[0]
    n_lin_coarse = int(np.floor(n_lin_fine / pyramid_scale)) + 1
    n_col_coarse = int(np.floor(n_col_fine / pyramid_scale)) + 1
    # The corner points must fit in the new grid.
    # This is always the case with the left most column or the upper line.
    # But you have to adjust ne number of columns to the right and number of lines below.
    # The final number of columns/lines has to fulfill:
    #     n_coarse = pyramid_scale * n_fine + 1
    # columns
    mod_col = n_col_fine % pyramid_scale
    pad_col = 0
    if mod_col != 1:
        pad_col = np.abs(1 - mod_col)
    # lines
    mod_lin = n_lin_fine % pyramid_scale
    pad_lin = 0
    if mod_lin != 1:
        pad_lin = np.abs(1 - mod_lin)
    # Here we extend it to right and below, so padding with edge is OK
    # Edge-mode otherwise creates artefacts on left and above.
    dem_fine = np.pad(dem_fine, ((0, pad_lin), (0, pad_col)), mode="edge")

    # Once you have data in the shape appropriate for resizing,
    # pad the data to support np.move.
    dem_fine = np.pad(
        dem_fine, ((-conv_from, conv_to), (-conv_from, conv_to)), mode="symmetric"
    )

    # Convolution (keep maximum)
    dem_convolve = np.zeros(dem_fine.shape)
    for i in np.arange(pyramid_scale) + conv_from:
        for j in np.arange(pyramid_scale) + conv_from:
            dem_convolve = np.maximum(
                dem_convolve, np.roll(dem_fine, (i, j), axis=(0, 1))
            )
    # Divide by pyramid_scale to account for the change of resolution
    # (important for the angle computation later on)
    dem_convolve = dem_convolve / pyramid_scale

    # Consider only the selected convoluted points according to the scale change.
    # As we select slice's end point make sure to consider at least 1 point more
    # to the right / below to really include it (Python way of considering end index).
    dem_coarse = dem_convolve[
        -conv_from : (n_lin_coarse * pyramid_scale + 1) : pyramid_scale,
        -conv_from : (n_col_coarse * pyramid_scale + 1) : pyramid_scale,
    ]

    # Final padding to enable searching the horizon over the edge:
    # use constant-mode set to the minimal height, so it doesn't
    # affect the horizon estimation.
    dem_coarse = np.pad(
        dem_coarse,
        ((max_radius, max_radius), (max_radius, max_radius)),
        mode="constant",
        constant_values=dem_coarse.min(),
    )

    return dem_coarse


def horizon_generate_pyramids(
    dem: npt.NDArray[Any],
    num_directions: int = 4,
    max_fine_radius: int = 100,
    max_pyramid_radius: int = 7,
    pyramid_scale: int = 3,
) -> Dict[int, Dict[str, Any]]:
    # In the levels higher than 1, determine the minimal search distance
    # and number of search distances.
    # If you have for instance
    #     pyramid_scale = 3
    #     max_pyramid_radius = 10
    #     num_directions = 8
    # then you have original distances in level 0:
    # 1, 2, 3, ... 9, 10
    # In level 1, your resolution is 3-times coarser.
    # The first pixel that takes that this new resolution,
    # has in original distance value 12 (in coarse resolution 4):
    # 12->4, 15->5, 18->6 ... 27->9, 30->10
    # So you start in the level 1 with min_pyramid_radius=4
    # and you search from 4 to 10 distances (n_pyramid_radius=7)
    min_pyramid_radius = int(np.floor(max_pyramid_radius / pyramid_scale)) + 1
    n_pyramid_radius = max_pyramid_radius - min_pyramid_radius + 1

    # get the convolution window indices
    conv_to = int(np.floor(pyramid_scale / 2.0))
    if (pyramid_scale % 2) == 0:
        conv_from = 1 - conv_to
    else:
        conv_from = -conv_to

    # initializations
    pyramid_levels = 0
    work = True
    pyramid = {}

    last_radius = None
    dem_coarse = None

    # Determine the number of levels and
    # the last radius to be used in the highest level.
    while work:
        _ = max_fine_radius / pyramid_scale**pyramid_levels
        if _ > max_pyramid_radius:
            pyramid_levels = pyramid_levels + 1
        else:
            work = False
            last_radius = np.round(
                max_fine_radius / pyramid_scale**pyramid_levels, decimals=0
            )

    # fill out the pyramid dict with the metadata required for horizon searching.
    for level in range(pyramid_levels + 1):
        # the level 0 contains the other min_radius as the rest of levels
        if level == 0:
            min_radius = 1
            dem_fine = np.copy(
                np.pad(
                    dem, max_pyramid_radius, mode="constant", constant_values=dem.min()
                )
            )
        else:
            assert dem_coarse is not None
            min_radius = min_pyramid_radius - 1
            dem_fine = np.copy(dem_coarse)
        # the last level contains the other radius_pixels as the rest of levels
        if level == pyramid_levels:
            assert last_radius is not None
            max_radius = last_radius
        else:
            max_radius = max_pyramid_radius
        # determine the dict of shifts
        shift = horizon_shift_vector(num_directions, max_radius, min_radius)
        dem_coarse = horizon_generate_coarse_dem(
            dem_fine, pyramid_scale, conv_from, conv_to, max_pyramid_radius
        )
        i_lin = np.arange(dem_fine.shape[0])
        i_col = np.arange(dem_fine.shape[1])

        pyramid[level] = {
            "num_directions": num_directions,
            "radius_pixels": max_radius,
            "min_radius": min_radius,
            "shift": shift,
            "dem": dem_fine,
            "i_lin": i_lin,
            "i_col": i_col,
        }

    return pyramid


def sky_illumination(
    dem: npt.NDArray[Any],
    resolution: float,
    sky_model: str = "overcast",
    compute_shadow: bool = False,
    shadow_horizon_only: bool = False,
    max_fine_radius: int = 100,
    num_directions: int = 32,
    shadow_az: float = 315,
    shadow_el: float = 35,
    ve_factor: float = 1,
    no_data: Optional[float] = None,
) -> Union[npt.NDArray[Any], Dict[str, npt.NDArray[Any]]]:
    """
    Compute topographic corrections for sky illumination.

    Parameters
    ----------
    dem : numpy.ndarray
        Input digital elevation model as 2D numpy array.
    resolution : float
        DEM pixel size.
    sky_model : str
        Sky model, it can be 'overcast' or 'uniform'.
    compute_shadow : bool
        If True it computes and adds shadow.
    shadow_horizon_only : bool
        Returns dict {"shadow": shadow, "horizon": horizon}
    max_fine_radius : int
        Max shadow modeling distance in pixels.
    num_directions : int
        Number of directions to search for horizon.
    shadow_az : int or float
        Shadow azimuth.
    shadow_el : int or float
        Shadow elevation.
    ve_factor : int or float
        Vertical exaggeration factor.
    no_data : int or float
        Value that represents no_data, all pixels with this value are changed to np.nan .

    Returns
    -------
    sky_illumination : numpy.ndarray
        2D numpy result array of Sky illumination.
    """
    # standard pyramid settings
    pyramid_scale = 2
    max_pyramid_radius = 20

    if not (10000 >= ve_factor >= -10000):
        raise Exception(
            "rvt.visualization.sky_illumination: ve_factor must be between -10000 and 10000!"
        )
    if shadow_az > 360 or shadow_az < 0:
        raise Exception(
            "rvt.visualization.sky_illumination: shadow_az must be between 0 and 360!"
        )
    if shadow_el > 90 or shadow_el < 0:
        raise Exception(
            "rvt.visualization.sky_illumination: shadow_el must be between 0 and 90!"
        )
    if resolution < 0:
        raise Exception(
            "rvt.visualization.sky_illumination: resolution must be a positive number!"
        )

    # change no_data to np.nan
    if no_data is not None:
        dem[dem == no_data] = np.nan

    dem = dem.astype(np.float32)
    dem = dem * ve_factor

    if sky_model.lower() == "overcast":
        compute_overcast = True
        compute_uniform = False
    elif sky_model.lower() == "uniform":
        compute_overcast = False
        compute_uniform = True
    else:
        raise Exception(
            "rvt.visualization.sky_illumination: sky_model must be overcast or uniform!"
        )

    # generate slope and aspect
    _ = slope_aspect(
        np.pad(dem, max_pyramid_radius, mode="symmetric"), resolution, resolution
    )
    slope = _["slope"]
    aspect = _["aspect"]

    # build DEM pyramids
    pyramid = horizon_generate_pyramids(
        dem,
        num_directions=num_directions,
        max_fine_radius=max_fine_radius,
        max_pyramid_radius=max_pyramid_radius,
        pyramid_scale=pyramid_scale,
    )
    n_levels = np.max([i for i in pyramid])

    # get the convolution window indices
    conv_to = int(np.floor(pyramid_scale / 2.0))
    if (pyramid_scale % 2) == 0:
        conv_from = 1 - conv_to
    else:
        conv_from = -conv_to
    # directional halve-resolution for integration limits
    da = np.pi / num_directions

    # init the intermediate results for uniform SI
    uniform_a = np.zeros(
        (dem.shape[0] + 2 * max_pyramid_radius, dem.shape[1] + 2 * max_pyramid_radius),
        dtype=np.float32,
    )
    uniform_b = np.copy(uniform_a)
    # init the output for overcast SI
    if compute_overcast:
        overcast_out = np.zeros(dem.shape, dtype=np.float32)
        overcast_c = np.zeros(
            (
                dem.shape[0] + 2 * max_pyramid_radius,
                dem.shape[1] + 2 * max_pyramid_radius,
            ),
            dtype=np.float32,
        )
        overcast_d = np.copy(overcast_c)
    else:
        overcast_out = None

    shadow_out = None
    horizon_out = None
    overcast_sh_out = None
    uniform_sh_out = None
    # init the output for shadows
    if compute_shadow:
        # use closest direction from pyramids as proxy for shadow azimuth
        # (just in case it is not the same as standard directions)
        _ = np.array([d for d in pyramid[0]["shift"]])
        i = np.argmin(np.abs(_ - (360 - shadow_az)))
        shadow_az = _[i]
        # binary shadows
        shadow_out = np.zeros(dem.shape, dtype=np.float32)
        # height of horizon in degrees
        horizon_out = np.zeros(dem.shape, dtype=np.float32)
        # overcast model + binary shadow
        if compute_overcast:
            overcast_sh_out = np.zeros(dem.shape, dtype=np.float32)
        else:
            overcast_sh_out = None
        # uniform model + binary shadow
        uniform_sh_out = np.zeros(dem.shape, dtype=np.float32)

    # search for horizon in each direction...
    for i_dir, direction in enumerate(pyramid[0]["shift"]):
        dir_rad = np.radians(direction)
        # reset maximum at each iteration (direction)
        max_slope = np.zeros(pyramid[n_levels]["dem"].shape, dtype=np.float32) - 1000

        for i_level in reversed(range(n_levels + 1)):
            height = pyramid[i_level]["dem"]
            move = pyramid[i_level]["shift"]

            # ... and to the search radius
            for i_rad, radius in enumerate(move[direction]["distance"]):
                # get shift index from move dictionary
                shift_indx = move[direction]["shift"][i_rad]
                # estimate the slope
                _ = np.maximum(
                    (np.roll(height, shift_indx, axis=(0, 1)) - height) / radius, 0.0
                )
                # compare to the previous max slope and keep the larges
                max_slope = np.maximum(max_slope, _)

            # resample the max_slope to a lower pyramid level
            if i_level > 0:
                lin_fine = pyramid[i_level - 1]["i_lin"] + (
                    conv_from + max_pyramid_radius * pyramid_scale - max_pyramid_radius
                )
                col_fine = pyramid[i_level - 1]["i_col"] + (
                    conv_from + max_pyramid_radius * pyramid_scale - max_pyramid_radius
                )
                lin_coarse = pyramid[i_level]["i_lin"] * pyramid_scale
                col_coarse = pyramid[i_level]["i_col"] * pyramid_scale
                interp_spline = RectBivariateSpline(
                    lin_coarse, col_coarse, max_slope, kx=1, ky=1
                )
                max_slope = interp_spline(lin_fine, col_fine)

        # convert to angle in radians and compute directional output
        _ = np.arctan(max_slope)
        uniform_a = uniform_a + (np.cos(_)) ** 2
        _d_aspect = -2 * np.sin(da) * np.cos(dir_rad - aspect)
        uniform_b = uniform_b + np.maximum(
            _d_aspect * (np.pi / 4.0 - _ / 2.0 - np.sin(2.0 * _) / 4.0), 0
        )
        if compute_overcast:
            _cos3 = (np.cos(_)) ** 3
            overcast_c = overcast_c + np.maximum(_cos3, 0)
            overcast_d = overcast_d + np.maximum(
                _d_aspect * (2.0 / 3.0 - np.cos(_) + _cos3 / 3.0), 0
            )
        if compute_shadow and (direction == shadow_az):
            horizon_out = np.degrees(
                _[
                    max_pyramid_radius:-max_pyramid_radius,
                    max_pyramid_radius:-max_pyramid_radius,
                ]
            )
            shadow_out = (horizon_out < shadow_el) * 1  # type: ignore
            if shadow_horizon_only:
                assert shadow_out is not None
                assert horizon_out is not None
                return {"shadow": shadow_out, "horizon": horizon_out}

    # because of numeric stability check if the uniform_b is less then pi
    uniform_out = da * np.cos(slope) * uniform_a + np.sin(slope) * np.minimum(
        uniform_b, np.pi
    )
    uniform_out = uniform_out[
        max_pyramid_radius:-max_pyramid_radius, max_pyramid_radius:-max_pyramid_radius
    ]

    if compute_overcast:
        overcast_out = (2.0 * da / 3.0) * np.cos(slope) * overcast_c + np.sin(
            slope
        ) * overcast_d
        overcast_out = overcast_out[  # type: ignore
            max_pyramid_radius:-max_pyramid_radius,
            max_pyramid_radius:-max_pyramid_radius,
        ]
        overcast_out = 0.33 * uniform_out + 0.67 * overcast_out
        overcast_out = overcast_out / overcast_out.max()  # type: ignore
    if compute_shadow:
        uniform_sh_out = 0.8 * uniform_out + 0.2 * shadow_out  # type: ignore
        if compute_overcast:
            overcast_sh_out = 0.8 * overcast_out + 0.2 * shadow_out  # type: ignore

    # normalize
    uniform_out = uniform_out / np.pi

    # # return results within dict
    # dict_sky_illumination = {"uniform": uniform_out,
    #                          "overcast": overcast_out,
    #                          "shadow": shadow_out,
    #                          "horizon": horizon_out,
    #                          "uniform_shaded": uniform_sh_out,
    #                          "overcast_shaded": overcast_sh_out,
    #                          }
    # dict_sky_illumination = {k: v for k, v in dict_sky_illumination.items() if v is not None}  # filter out none
    # return dict_sky_illumination

    # output
    if compute_uniform and not compute_shadow:
        assert uniform_out is not None
        return uniform_out
    elif compute_uniform and compute_shadow:
        assert uniform_sh_out is not None
        return uniform_sh_out
    elif compute_overcast and not compute_shadow:
        assert overcast_out is not None
        return overcast_out
    elif compute_overcast and compute_shadow:
        assert overcast_sh_out is not None
        return overcast_sh_out
    else:
        raise RuntimeError


def shadow_horizon(
    dem: npt.NDArray[Any], resolution: float, shadow_az: float=315, shadow_el: float=35, ve_factor: float=1, no_data: Optional[float]=None
) -> Dict[str, npt.NDArray[Any]]:
    """
    Compute shadow and horizon.

    Parameters
    ----------
    dem : numpy.ndarray
        Input digital elevation model as 2D numpy array.
    resolution : float
        DEM pixel size.
    shadow_az : int or float
        Shadow azimuth.
    shadow_el : int or float
        Shadow elevation.
    ve_factor : int or float
        Vertical exaggeration factor.
    no_data : int or float
        Value that represents no_data, all pixels with this value are changed to np.nan .

    Returns
    -------
    dict_out : dict
        Returns {"shadow": shadow, "horizon": horizon};
        shadow : 2D binary numpy array (numpy.ndarray) of shadows;
        horizon; 2D numpy array (numpy.ndarray) of horizon.
    """
    if not (10000 >= ve_factor >= -10000):
        raise Exception(
            "rvt.visualization.shadow_horizon: ve_factor must be between -10000 and 10000!"
        )
    if shadow_az > 360 or shadow_az < 0:
        raise Exception(
            "rvt.visualization.shadow_horizon: shadow_az must be between 0 and 360!"
        )
    if shadow_el > 90 or shadow_el < 0:
        raise Exception(
            "rvt.visualization.shadow_horizon: shadow_el must be between 0 and 90!"
        )
    if resolution < 0:
        raise Exception(
            "rvt.visualization.shadow_horizon: resolution must be a positive number!"
        )

    return sky_illumination(
        dem=dem,
        resolution=resolution,
        compute_shadow=True,
        shadow_horizon_only=True,
        shadow_el=shadow_el,
        shadow_az=shadow_az,
        ve_factor=ve_factor,
        no_data=no_data,
    )


def msrm(
    dem, resolution, feature_min, feature_max, scaling_factor, ve_factor=1, no_data=None
):
    """
    Compute Multi-scale relief model (MSRM).

    Parameters
    ----------
    dem : numpy.ndarray
        Input digital elevation model as 2D numpy array.
    resolution : float
        DEM pixel size.
    feature_min: float
        Minimum size of the feature you want to detect in meters.
    feature_max: float
        Maximum size of the feature you want to detect in meters.
    scaling_factor: int
        Scaling factor, if larger than 1 it provides larger range of MSRM values (increase contrast and visibility),
        but could result in a loss of sensitivity for intermediate sized features.
    ve_factor : int or float
        Vertical exaggeration factor.
    no_data : int or float
        Value that represents no_data, all pixels with this value are changed to np.nan .

    Returns
    -------
    msrm_out : numpy.ndarray
        2D numpy result array of Multi-scale relief model.
    """
    if not (10000 >= ve_factor >= -10000):
        raise Exception(
            "rvt.visualization.msrm: ve_factor must be between -10000 and 10000!"
        )
    if resolution < 0:
        raise Exception("rvt.visualization.msrm: resolution must be a positive number!")

    # change no_data to np.nan
    if no_data is not None:
        dem[dem == no_data] = np.nan

    dem = dem.astype(np.float32)
    dem = dem * ve_factor

    if feature_min < resolution:  # feature_min can't be smaller than resolution
        feature_min = resolution

    scaling_factor = int(scaling_factor)  # has to be integer

    # calculation of i and n (from article)
    i = int(
        np.floor(
            ((feature_min - resolution) / (2 * resolution)) ** (1 / scaling_factor)
        )
    )
    n = int(
        np.ceil(((feature_max - resolution) / (2 * resolution)) ** (1 / scaling_factor))
    )

    # lpf = low pass filter
    relief_models_sum = np.zeros(dem.shape)  # sum of all substitution of 2 consecutive
    nr_relief_models = (
        0  # number of additions (substitutions of 2 consecutive surfaces)
    )
    last_lpf_surface = 0

    # generation of filtered surfaces (lpf_surface)
    for ndx in range(i, n + 1, 1):
        kernel_radius = ndx**scaling_factor
        # calculate mean filtered surface
        lpf_surface = mean_filter(dem=dem, kernel_radius=kernel_radius)
        if not ndx == i:  # if not first surface
            relief_models_sum += (
                last_lpf_surface - lpf_surface
            )  # substitution of 2 consecutive lpf_surface
            nr_relief_models += 1
        last_lpf_surface = lpf_surface

    msrm_out = relief_models_sum / nr_relief_models

    return msrm_out


def integral_image(dem, data_type=np.float64):
    """
    Calculates integral image (summed-area table), where origin is left upper corner.

    Parameters
    ----------
    dem : numpy.ndarray
        Input digital elevation model as 2D numpy array.
    data_type : np.__class__
        dtype as numpy data type class (np.float64, np.int8, etc.)

    Returns
    -------
    msrm_out : numpy.ndarray
        Cumulative sum of the elements along each axis of a 2D array.

    References
    ----------
    https://en.wikipedia.org/wiki/Summed-area_table

    Examples
    --------
    In: print(integral_image(np.array([[7, 4, 7, 2],
    ... [6, 9, 9, 5],
    ... [6, 6, 7, 6]])))

    Out: [[ 7. 11. 18. 20.]
          [13. 26. 42. 49.]
          [19. 38. 61. 74.]]
    """
    dem = dem.astype(data_type)
    return dem.cumsum(axis=0).cumsum(axis=1)


def topographic_dev(dem, dem_i_nr_pixels, dem_i1, dem_i2, kernel_radius):
    """
    Calculates topographic DEV - Deviation from mean elevation. DEV(D) = (z0 - zmD) / sD.
    Where D is radius of kernel, z0 is center pixel value, zmD is mean of all kernel values,
    sD is standard deviation of kernel.

    Parameters
    ----------
    dem : numpy.ndarray
        Input digital elevation model as 2D numpy array.
    dem_i_nr_pixels : numpy.ndarray
        Summed area table (integral image) of number of pixels.
    dem_i1 : numpy.ndarray
        Summed area table (integral image) of dem.
    dem_i2 : numpy.ndarray
        Summed area table (integral image) of dem squared (dem**2).
    kernel_radius : int
        Kernel radius (D).

    Returns
    -------
    dev_out : numpy.ndarray
        2D numpy result array of topographic DEV - Deviation from mean elevation.
    """
    radius_cell = int(kernel_radius)
    if radius_cell <= 0:
        return dem

    kernel_nr_pix_arr = (
        np.roll(dem_i_nr_pixels, (radius_cell, radius_cell), axis=(0, 1))
        + np.roll(dem_i_nr_pixels, (-radius_cell - 1, -radius_cell - 1), axis=(0, 1))
        - np.roll(dem_i_nr_pixels, (-radius_cell - 1, radius_cell), axis=(0, 1))
        - np.roll(dem_i_nr_pixels, (radius_cell, -radius_cell - 1), axis=(0, 1))
    )

    # sum
    dem_mean = (
        np.roll(dem_i1, (radius_cell, radius_cell), axis=(0, 1))
        + np.roll(dem_i1, (-radius_cell - 1, -radius_cell - 1), axis=(0, 1))
        - np.roll(dem_i1, (-radius_cell - 1, radius_cell), axis=(0, 1))
        - np.roll(dem_i1, (radius_cell, -radius_cell - 1), axis=(0, 1))
    )
    # divide with nr of pixels inside kernel
    with np.errstate(
        divide="ignore", invalid="ignore"
    ):  # Suppress warning for dividing by zero
        dem_mean = dem_mean / kernel_nr_pix_arr

    # std
    dem_std = (
        np.roll(dem_i2, (radius_cell, radius_cell), axis=(0, 1))
        + np.roll(dem_i2, (-radius_cell - 1, -radius_cell - 1), axis=(0, 1))
        - np.roll(dem_i2, (-radius_cell - 1, radius_cell), axis=(0, 1))
        - np.roll(dem_i2, (radius_cell, -radius_cell - 1), axis=(0, 1))
    )

    with np.errstate(
        divide="ignore", invalid="ignore"
    ):  # Suppress warning for dividing by zero
        dem_std = np.sqrt(np.abs(dem_std / kernel_nr_pix_arr - dem_mean**2))
        # returns nan values where division by zero happens

    dev_out = (np.roll(dem, (-1, -1), axis=(0, 1)) - dem_mean) / (
        dem_std + 1e-6
    )  # add 1e-6 to prevent division with 0

    return dev_out


def max_elevation_deviation(dem, minimum_radius, maximum_radius, step):
    """
    Calculates maximum deviation from mean elevation, dev_max (Maximum Deviation from mean elevation) for each
    grid cell in a digital elevation model (DEM) across a range specified spatial scales.

    Parameters
    ----------
    dem : numpy.ndarray
        Input digital elevation model as 2D numpy array.
    minimum_radius : int
        Minimum radius to calculate DEV (topographic_dev).
    maximum_radius : int
        Maximum radius to calculate DEV (topographic_dev).
    step : int
        Step from minimum to maximum radius to calc DEV (topographic_dev).

    Returns
    -------
    dev_out : numpy.ndarray
        2D numpy result array of maxDEV - Maximum Deviation from mean elevation.
    """
    minimum_radius = int(minimum_radius)
    maximum_radius = int(maximum_radius)
    step = int(step)

    # store positions of nan
    idx_nan_dem = np.isnan(dem)

    dem_pad = np.pad(dem, (maximum_radius + 1, maximum_radius), mode="symmetric")
    # store nans
    idx_nan_dem_pad = np.isnan(dem_pad)
    # change nan to 0
    dem_pad[idx_nan_dem_pad] = 0

    # number of pixels for summed area table
    dem_i_nr_pixels = np.ones(dem_pad.shape)
    dem_i_nr_pixels[idx_nan_dem_pad] = 0
    dem_i_nr_pixels = integral_image(dem_i_nr_pixels, np.int64)

    # This outputs float64, which is by design. Change final array to float32 at the end of the function (at return)
    dem_i1 = integral_image(dem_pad)
    dem_i2 = integral_image(dem_pad**2)

    for kernel_radius in range(minimum_radius, maximum_radius + 1, step):
        dev = topographic_dev(dem_pad, dem_i_nr_pixels, dem_i1, dem_i2, kernel_radius)[
            maximum_radius : -(maximum_radius + 1),
            maximum_radius : -(maximum_radius + 1),
        ]
        if kernel_radius == minimum_radius:
            dev_max_out = dev
            rad_max_out = np.zeros_like(dev, dtype=np.float32) + kernel_radius
        else:
            rad_max_out = np.where(
                np.abs(dev_max_out) >= np.abs(dev), rad_max_out, kernel_radius
            )
            dev_max_out = np.where(np.abs(dev_max_out) >= np.abs(dev), dev_max_out, dev)
    # rad_max_out, radius of DEV for maxDEV (for each pixel)

    # change where dem nan back to nan
    dev_max_out[idx_nan_dem] = np.nan
    rad_max_out[idx_nan_dem] = np.nan

    return dev_max_out.astype(np.float32)


def mstp(
    dem,
    local_scale=(3, 21, 2),
    meso_scale=(23, 203, 18),
    broad_scale=(223, 2023, 180),
    lightness=1.2,
    ve_factor=1,
    no_data=None,
):
    """
    Compute Multi-scale topographic position (MSTP).

    Parameters
    ----------
    dem : numpy.ndarray
        Input digital elevation model as 2D numpy array.
    local_scale : tuple(int, int, int)
        Input local scale minimum radius (local_scale[0]), maximum radius (local_scale[1]), step (local_scale[2]).
    meso_scale : tuple(int, int, int)
        Input meso scale minimum radius (meso_scale[0]), maximum radius (meso_scale[1]), step (meso_scale[2]).
    broad_scale : tuple(int, int, int)
        Input broad scale minimum radius (broad_scale[0]), maximum radius (broad_scale[1]), step (broad_scale[2]).
    lightness : float
        Lightness of image.
    ve_factor : int or float
        Vertical exaggeration factor.
    no_data : int or float
        Value that represents no_data, all pixels with this value are changed to np.nan .

    Returns
    -------
    msrm_out : numpy.ndarray
        3D numpy RGB result array of Multi-scale topographic position.
    """
    if (
        local_scale[0] > local_scale[1]
        or meso_scale[0] > meso_scale[1]
        or broad_scale[0] > broad_scale[1]
    ):
        raise Exception(
            "rvt.visualization.mstp: local_scale, meso_scale, broad_scale min has to be smaller than max!"
        )
    if (
        (local_scale[1] - local_scale[0] < local_scale[2])
        or (meso_scale[1] - meso_scale[0] < meso_scale[2])
        or (broad_scale[1] - broad_scale[0] < broad_scale[2])
    ):
        raise Exception(
            "rvt.visualization.mstp: local_scale, meso_scale, broad_scale step has"
            " to be within min and max!"
        )
    if not (10000 >= ve_factor >= -1000):
        raise Exception(
            "rvt.visualization.mstp: ve_factor must be between -10000 and 10000!"
        )

    # change no_data to np.nan
    if no_data is not None:
        dem[dem == no_data] = np.nan

    dem = dem.astype(np.float32)
    dem = dem * ve_factor

    local_dev = max_elevation_deviation(
        dem=dem,
        minimum_radius=local_scale[0],
        maximum_radius=local_scale[1],
        step=local_scale[2],
    )
    meso_dev = max_elevation_deviation(
        dem=dem,
        minimum_radius=meso_scale[0],
        maximum_radius=meso_scale[1],
        step=meso_scale[2],
    )
    broad_dev = max_elevation_deviation(
        dem=dem,
        minimum_radius=broad_scale[0],
        maximum_radius=broad_scale[1],
        step=broad_scale[2],
    )

    cutoff = lightness
    # RGB order - broad, meso, local
    red = 1 - np.exp(-cutoff * np.abs(broad_dev))
    green = 1 - np.exp(-cutoff * np.abs(meso_dev))
    blue = 1 - np.exp(-cutoff * np.abs(local_dev))

    red[red < 0] = 0
    green[green < 0] = 0
    blue[blue < 0] = 0

    red[red > 1] = 1
    green[green > 1] = 1
    blue[blue > 1] = 1

    return np.asarray([red, green, blue])  # RGB float32 (3 x 32bit)


def fill_where_nan(dem, method="idw"):
    """
    Replaces np.nan values, with interpolation (extrapolation).

    Parameters
    -------
    dem : numpy.ndarray
        Input digital elevation model as 2D numpy array.
    method : {'linear_row', 'idw_r_p', 'kd_tree', 'nearest_neighbour'}
        'linear_row', Linear row interpolation, array is flattened and then linear interpolation is performed.
        This method is fast but very inaccurate.
        'idw_r_p', Inverse Distance Weighting interpolation. If you only input idw it will take default parameters
        (r=20, p=2). You can also input interpolation radius (r) and power (p) for weights. (Example:
        idw_5_2 means radius = 5, power = 2.)
        'kd_tree', K-D Tree interpolation.
        'nearest_neighbour', Nearest neighbour interpolation.
    """
    if np.all(~np.isnan(dem)):  # if there is no nan return dem
        return dem

    dem_out = np.copy(dem)
    mask = np.isnan(dem_out)

    if method == "linear_row":
        # 1D row linear interpolation
        dem_out[mask] = np.interp(
            np.flatnonzero(mask), np.flatnonzero(~mask), dem_out[~mask]
        )

    elif method.split("_")[0] == "idw":
        radius = 20
        power = 2
        if len(method.split("_")) == 3:
            radius = int(method.split("_")[1])
            power = float(method.split("_")[2])
        nan_idx = zip(*np.where(mask))  # find nan positions
        for i_row, i_column in nan_idx:  # iterate through nans
            # nan surrounding array based on radius
            i_row_start = i_row - radius  # start row idx of nan surrounding
            i_column_start = i_column - radius  # start col idx of nan surrounding
            i_row_end = i_row + radius + 1  # end row idx of nan surrounding
            i_column_end = i_column + radius + 1  # end col idx of nan surrounding
            # row idx center pixel (nan pixel to idw) of nan surrounding array
            i_row_center = i_row - i_row_start
            # col idx center pixel (nan pixel to idw) of nan surrounding array
            i_column_center = i_column - i_column_start
            if i_row_start < 0:  # edge
                i_row_center = i_row
                i_row_start = 0
            if i_column_start < 0:  # edge
                i_column_center = i_column
                i_column_start = 0
            if i_row_end > dem.shape[0]:  # edge
                i_row_end = dem.shape[0]
            if i_column_end > dem.shape[1]:  # edge
                i_column_end = dem.shape[0]
            nan_surrounding_arr = dem[
                i_row_start:i_row_end, i_column_start:i_column_end
            ]
            if np.all(np.isnan(nan_surrounding_arr)):  # whole surrounding array is nan
                dem_out[i_row, i_column] = np.nan
            else:
                # calculate distance array (wight matrix)
                dist_arr = np.ones(nan_surrounding_arr.shape)  # all ones
                # center pixel is 0 to calc distance matrix around 0 pixel with distance_transform_edt
                dist_arr[i_row_center, i_column_center] = 0
                dist_arr = distance_transform_edt(dist_arr)
                dist_arr[dist_arr == 0] = np.nan  # can't divide with zero
                dist_arr = 1 / dist_arr**power
                nan_mask = np.isnan(nan_surrounding_arr)
                dist_arr[nan_mask] = 0  # where nan weight matrix is zero
                # calculate idw for one nan value
                dem_out[i_row, i_column] = np.nansum(
                    nan_surrounding_arr * dist_arr
                ) / np.nansum(dist_arr)

    elif (
        method == "kd_tree"
        or method == "nearest_neighbour"
        or method == "nearest_neighbor"
    ):
        x, y = np.mgrid[0 : dem_out.shape[0], 0 : dem_out.shape[1]]
        xy_good = np.array((x[~mask], y[~mask])).T
        xy_bad = np.array((x[mask], y[mask])).T

        # cKD-Tree (K-D Tree) interpolation
        # https://stackoverflow.com/questions/3662361/fill-in-missing-values-with-nearest-neighbour-in-python-numpy-masked-arrays
        if method == "kd_tree":
            leaf_size = 1000
            dem_out[mask] = dem_out[~mask][
                cKDTree(data=xy_good, leafsize=leaf_size).query(xy_bad)[1]
            ]

        # Nearest neighbour interpolation
        elif method == "nearest_neighbour" or method == "nearest_neighbor":
            dem_out[mask] = griddata(xy_good, dem_out[~mask], xy_bad, method="nearest")

    else:
        raise Exception("rvt.visualization.fill_where_nan: Wrong method!")

    return dem_out
