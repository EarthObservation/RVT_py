"""
Relief Visualization Toolbox – Visualization Functions

Contains functions for computing the visualizations.

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

# python libraries
import numpy as np
import scipy.interpolate
import warnings
import math


def byte_scale(data,
               c_min=None,
               c_max=None,
               high=255,
               low=0,
               no_data=None
               ):
    """
    Remade old scipy function.
    Byte scales an array (image). Linear scale.

    Byte scaling means converting the input image to uint8 dtype and scaling
    the range to ``(low, high)`` (default 0-255).

    Parameters
    ----------
    data : numpy image data array.
    c_min : scalar, Bias scaling of small values. Default is ``data.min()``.
    c_max : scalar, Bias scaling of large values. Default is ``data.max()``.
    high : scalar, Scale max value to `high`.  Default is 255.
    low : scalar, Scale min value to `low`.  Default is 0.
    no_data : value that represents no_data, it is changed to np.nan

    Returns
    -------
    img_array : uint8 ndarray
        The byte-scaled array.
    """

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
        # TODO: the following line seems not good to me - if cmin=0, then that pixel will get negative value
        byte_data = (high + 1) * (data - c_min - 1) / (c_max - c_min)  # copied from IDL BYTSCL
        byte_data[byte_data > high] = high
        byte_data[byte_data < 0] = 0
        return np.cast[np.uint8](byte_data) + np.cast[np.uint8](low)

    # scale = float(high - low) / cscale  # old scipy fn
    # byte_data = (data * 1.0 - cmin) * scale + 0.4999  # old scipy fn

    byte_data = (high + 0.9999) * (data - c_min) / (c_max - c_min)  # copied from IDL BYTSCL
    byte_data[byte_data > high] = high
    byte_data[byte_data < 0] = 0
    byte_data[np.isnan(byte_data)] = 0  # change no_data to 0
    return np.cast[np.uint8](byte_data) + np.cast[np.uint8](low)


def slope_aspect(dem,
                 resolution_x,
                 resolution_y,
                 output_units="radian",
                 ve_factor=1,
                 no_data=None,
                 fill_no_data=False,
                 keep_original_no_data=False
                 ):
    """
    Procedure can return terrain slope and aspect in radian units (default) or in alternative units (if specified).
    Slope is defined as 0 for Hz plane and pi/2 for vertical plane.
    Aspect iz defined as geographic azimuth: clockwise increasing, 0 or 2pi for the North direction.
    Currently applied finite difference method.

    Parameters
    ----------
    dem : input dem 2D numpy array
    resolution_x : dem resolution in X direction
    resolution_y : DEM resolution in Y direction
    output_units : percent, degree, radians
    ve_factor : vertical exaggeration factor
    no_data : value that represents no_data, all pixels with this value are changed to np.nan
    fill_no_data : if True it fills where np.nan (no_data) with mean of surrounding pixels (3x3)
    keep_original_no_data : if True it changes all output pixels to np.nan where dem has no_data

    Returns
    -------
    {"slope": slope_out, "aspect": aspect_out} : dictionaries with 2D numpy arrays
    """
    if dem.ndim != 2:
        raise Exception("rvt.vis.slope_aspect: dem has to be 2D np.array!")
    if not (1000 >= ve_factor >= -1000):
        raise Exception("rvt.vis.slope_aspect: ve_factor must be between -1000 and 1000!")
    if resolution_x < 0 or resolution_y < 0:
        raise Exception("rvt.vis.slope_aspect: resolution must be a positive number!")
    if no_data is None and fill_no_data:
        warnings.warn("rvt.vis.slope_aspect: In order to fill no data (fill_no_data = True) you have to input"
                      " no_data!")
    if (no_data is None or not fill_no_data) and keep_original_no_data:
        warnings.warn("rvt.vis.slope_aspect: In order to keep original no data (keep_original_no_data = True)"
                      " you have to input no_data and fill_no_data has to be True!")

    # change no_data to np.nan
    idx_no_data = None
    if no_data is not None:
        if keep_original_no_data:  # save indexes where is no_data
            if np.isnan(no_data):
                idx_no_data = np.where(np.isnan(dem))
            else:
                idx_no_data = np.where(dem == no_data)
        dem[dem == no_data] = np.nan

    dem = dem.astype(np.float32)
    dem = dem * ve_factor

    # fill no data with mean of surrounding pixels
    if fill_no_data:
        dem = fill_where_nan(dem)

    # derivatives in X and Y direction
    dzdx = ((np.roll(dem, 1, axis=1) - np.roll(dem, -1, axis=1)) / 2) / resolution_x
    dzdy = ((np.roll(dem, -1, axis=0) - np.roll(dem, 1, axis=0)) / 2) / resolution_y
    tan_slope = np.sqrt(dzdx ** 2 + dzdy ** 2)

    # Compute slope
    if output_units == "percent":
        slope_out = tan_slope * 100
    elif output_units == "degree":
        slope_out = np.rad2deg(np.arctan(tan_slope))
    elif output_units == "radian":
        slope_out = np.arctan(tan_slope)
    else:
        raise Exception("rvt.vis.calculate_slope: Wrong function input 'output_units'!")

    # compute Aspect
    # aspect identifies the down slope direction of the maximum rate of change in value from each cell to its neighbors:
    #     0
    # 270    90
    #    180
    dzdy[dzdy == 0] = 10e-9  # important for numeric stability - where dzdy is zero, make tangens to really high value

    aspect_out = np.arctan2(dzdx, dzdy)  # atan2 took care of the quadrants
    if output_units == "degree":
        aspect_out = np.rad2deg(aspect_out)

    # edges to -1
    slope_out[:, 0] = -1
    slope_out[0, :] = -1
    slope_out[:, -1] = -1
    slope_out[-1, :] = -1

    # change result to np.nan where dem is no_data
    if no_data is not None and keep_original_no_data:
        slope_out[idx_no_data] = np.nan
        aspect_out[idx_no_data] = np.nan

    return {"slope": slope_out, "aspect": aspect_out}


def hillshade(dem,
              resolution_x,
              resolution_y,
              sun_azimuth=315,
              sun_elevation=35,
              slope=None,
              aspect=None,
              ve_factor=1,
              no_data=None,
              fill_no_data=False,
              keep_original_no_data=False
              ):
    """
    Compute hillshade.

    Parameters
    ----------
    dem : input DEM 2D numpy array
    resolution_x : DEM resolution in X direction
    resolution_y : DEM resolution in Y direction
    sun_azimuth : solar azimuth angle (clockwise from North) in degrees
    sun_elevation : solar vertical angle (above the horizon) in degrees
    slope : slope arr in radians if you don't input it, it is calculated
    aspect : aspect arr in radians if you don't input it, it is calculated
    ve_factor : vertical exaggeration factor
    no_data : value that represents no_data, all pixels with this value are changed to np.nan
    fill_no_data : if True it fills where np.nan (no_data) with mean of surrounding pixels (3x3)
    keep_original_no_data : if True it changes all output pixels to np.nan where dem has no_data

    Returns
    -------
    hillshade_out : result numpy array
    """
    if dem.ndim != 2:
        raise Exception("rvt.vis.hillshade: dem has to be 2D np.array!")
    if not (1000 >= ve_factor >= -1000):
        raise Exception("rvt.vis.hillshade: ve_factor must be between -1000 and 1000!")
    if sun_azimuth > 360 or sun_elevation > 90 or sun_azimuth < 0 or sun_elevation < 0:
        raise Exception("rvt.vis.hillshade: sun_azimuth must be [0-360] and sun_elevation [0-90]!")
    if resolution_x < 0 or resolution_y < 0:
        raise Exception("rvt.vis.hillshade: resolution must be a positive number!")
    if no_data is None and fill_no_data:
        warnings.warn("rvt.vis.hillshade: In order to fill no data (fill_no_data = True) you have to input"
                      " no_data!")
    if (no_data is None or not fill_no_data) and keep_original_no_data:
        warnings.warn("rvt.vis.hillshade: In order to keep original no data (keep_original_no_data ="
                      " True) you have to input no_data and fill_no_data has to be True!")

    # change no_data to np.nan
    idx_no_data = None
    if no_data is not None:
        if keep_original_no_data:  # save indexes where is no_data
            if np.isnan(no_data):
                idx_no_data = np.where(np.isnan(dem))
            else:
                idx_no_data = np.where(dem == no_data)
        dem[dem == no_data] = np.nan

    dem = dem.astype(np.float32)
    dem = dem * ve_factor

    # fill no data with mean of surrounding pixels
    if fill_no_data:
        dem = fill_where_nan(dem)

    # Convert solar position (degrees) to radians
    sun_azimuth_rad = np.deg2rad(sun_azimuth)
    sun_elevation_rad = np.deg2rad(sun_elevation)

    # Convert to solar zenith angle
    sun_zenith_rad = np.pi / 2 - sun_elevation_rad

    # are slope and aspect already calculated and presented
    if slope is None or aspect is None:
        # calculates slope and aspect
        dict_slp_asp = slope_aspect(dem=dem, resolution_x=resolution_x, resolution_y=resolution_y,
                                    output_units="radian")
        slope = dict_slp_asp["slope"]
        aspect = dict_slp_asp["aspect"]

    # Compute solar incidence angle, hillshading
    hillshade_out = np.cos(sun_zenith_rad) * np.cos(slope) + np.sin(sun_zenith_rad) * np.sin(slope) * np.cos(
        aspect - sun_azimuth_rad)

    # edges to -1
    hillshade_out[:, 0] = -1
    hillshade_out[0, :] = -1
    hillshade_out[:, -1] = -1
    hillshade_out[-1, :] = -1

    # change result to np.nan where dem is no_data
    if no_data is not None and keep_original_no_data:
        hillshade_out[idx_no_data] = np.nan

    return hillshade_out


def multi_hillshade(dem,
                    resolution_x,
                    resolution_y,
                    nr_directions=16,
                    sun_elevation=35,
                    slope=None,
                    aspect=None,
                    ve_factor=1,
                    no_data=None,
                    fill_no_data=False,
                    keep_original_no_data=False
                    ):
    """
    Calculates hillshades from multiple directions.

    Parameters
    ----------
    dem : input DEM 2D numpy array
    resolution_x : DEM resolution in X direction
    resolution_y : DEM resolution in Y direction
    nr_directions : number of solar azimuth angles (clockwise from North)
    sun_elevation : solar vertical angle (above the horizon) in degrees
    slope : slope in radians if you don't input it, it is calculated
    aspect : aspect in radians if you don't input it, it is calculated
    ve_factor : vertical exaggeration factor
    no_data : value that represents no_data, all pixels with this value are changed to np.nan
    fill_no_data : if True it fills where np.nan (no_data) with mean of surrounding pixels (3x3)
    keep_original_no_data : if True it changes all output pixels to np.nan where dem has no_data

    Returns
    -------
    multi_hillshade_out : numpy array containing numpy_arrays of hillshades in different directions
    """
    if dem.ndim != 2:
        raise Exception("rvt.vis.multi_hillshade: dem has to be 2D np.array!")
    if sun_elevation > 90 or sun_elevation < 0:
        raise Exception("rvt.vis.multi_hillshade: sun_elevation must be [0-90]!")
    if resolution_x < 0 or resolution_y < 0:
        raise Exception("rvt.vis.multi_hillshade: resolution must be a positive number!")
    if nr_directions < 1:
        raise Exception("rvt.vis.multi_hillshade: nr_directions must be a positive number!")
    if not (1000 >= ve_factor >= -1000):
        raise Exception("rvt.vis.multi_hillshade: ve_factor must be between -1000 and 1000!")
    if no_data is None and fill_no_data:
        warnings.warn("rvt.vis.multi_hillshade: In order to fill no data (fill_no_data = True) you have to input"
                      " no_data!")
    if (no_data is None or not fill_no_data) and keep_original_no_data:
        warnings.warn("rvt.vis.multi_hillshade: In order to keep original no data (keep_original_no_data = True)"
                      " you have to input no_data and fill_no_data has to be True!")

    # change no_data to np.nan
    idx_no_data = None
    if no_data is not None:
        if keep_original_no_data:  # save indexes where is no_data
            if np.isnan(no_data):
                idx_no_data = np.where(np.isnan(dem))
            else:
                idx_no_data = np.where(dem == no_data)
        dem[dem == no_data] = np.nan

    dem = dem.astype(np.float32)
    dem = dem * ve_factor

    # fill no data with mean of surrounding pixels
    if fill_no_data:
        dem = fill_where_nan(dem)

    # calculates slope and aspect if they are not added
    if slope is None or aspect is None:  # slope and aspect are the same, so we have to calculate it once
        dict_slp_asp = slope_aspect(dem=dem, resolution_x=resolution_x, resolution_y=resolution_y,
                                    output_units="radian")
        slope = dict_slp_asp["slope"]
        aspect = dict_slp_asp["aspect"]

    hillshades_arr_list = []  # list of all hillshades in diffrent directions
    for i_direction in range(nr_directions):
        sun_azimuth = (360 / nr_directions) * i_direction
        hillshading = hillshade(dem=dem, resolution_x=resolution_x, resolution_y=resolution_y,
                                sun_elevation=sun_elevation, sun_azimuth=sun_azimuth, slope=slope, aspect=aspect)
        # change result to np.nan where dem is no_data
        if no_data is not None and keep_original_no_data:
            hillshading[idx_no_data] = np.nan
        hillshades_arr_list.append(hillshading)
    multi_hillshade_out = np.asarray(hillshades_arr_list)

    return multi_hillshade_out


def slrm(dem,
         radius_cell=20,
         ve_factor=1,
         no_data=None,
         fill_no_data=False,
         keep_original_no_data=False
         ):
    """
    Calculates Simple local relief model.

    Parameters
    ----------
    dem : input DEM 2D numpy array
    radius_cell : Radius for trend assessment [pixels]
    ve_factor : vertical exaggeration factor
    no_data : value that represents no_data, all pixels with this value are changed to np.nan
    fill_no_data : if True it fills where np.nan (no_data) with mean of surrounding pixels (3x3)
    keep_original_no_data : if True it changes all output pixels to np.nan where dem has no_data

    Returns
    -------
    slrm_out : slrm 2D numpy array
    """
    if dem.ndim != 2:
        raise Exception("rvt.vis.slrm: dem has to be 2D np.array!")
    if radius_cell < 10 or radius_cell > 50:
        raise Exception("rvt.vis.slrm: Radius for trend assessment needs to be in interval 10-50 pixels!")
    if not (1000 >= ve_factor >= -1000):
        raise Exception("rvt.vis.slrm: ve_factor must be between -1000 and 1000!")
    if no_data is None and fill_no_data:
        warnings.warn("rvt.vis.slrm: In order to fill no data (fill_no_data = True) you have to input"
                      " no_data!")
    if (no_data is None or not fill_no_data) and keep_original_no_data:
        warnings.warn("rvt.vis.slrm: In order to keep original no data (keep_original_no_data = True)"
                      " you have to input no_data and fill_no_data has to be True!")

    # change no_data to np.nan
    idx_no_data = None
    if no_data is not None:
        if keep_original_no_data:  # save indexes where is no_data
            if np.isnan(no_data):
                idx_no_data = np.where(np.isnan(dem))
            else:
                idx_no_data = np.where(dem == no_data)
        dem[dem == no_data] = np.nan

    dem = dem.astype(np.float32)
    dem = dem * ve_factor

    # fill no data with mean of surrounding pixels
    if fill_no_data:
        dem = fill_where_nan(dem)

    # mean filter
    radius_cell = int(radius_cell)  # nr_rolls in each direction
    dem_padded = np.pad(array=dem, pad_width=radius_cell, mode="edge")  # padding
    slrm_out = np.copy(dem_padded)
    for i_y_roll in range(radius_cell):
        roll = i_y_roll + 1  # y direction roll
        slrm_out += np.roll(np.copy(dem_padded), roll, axis=0)  # roll positive direction
        slrm_out += np.roll(np.copy(dem_padded), -roll, axis=0)  # roll negative direction
    y_rolls_sum = np.copy(slrm_out)  # sum of all rolls in y direction
    for i_x_roll in range(radius_cell):  # x direction roll
        roll = i_x_roll + 1
        slrm_out += np.roll(np.copy(y_rolls_sum), roll, axis=1)  # roll positive direction
        slrm_out += np.roll(np.copy(y_rolls_sum), -roll, axis=1)  # roll negative direction
    del y_rolls_sum
    slrm_out = slrm_out / ((2 * radius_cell + 1) ** 2)  # calculate mean, 1=current pixel
    slrm_out = slrm_out[radius_cell:-radius_cell, radius_cell:-radius_cell]  # remove padding
    slrm_out = dem - slrm_out

    # change result to np.nan where dem is no_data
    if no_data is not None and keep_original_no_data:
        slrm_out[idx_no_data] = np.nan

    return slrm_out


def horizon_shift_vector(num_directions=16,
                         radius_pixels=10,
                         min_radius=1
                         ):
    """
    Calculates Sky-View determination movements.

    Parameters
    ----------
    num_directions : number of directions as input
    radius_pixels : radius to consider in pixels (not in meters)
    min_radius : radius to start searching for horizon in pixels (not in meters)

    Returns
    -------
    shift : dict with keys corresponding to the directions of search azimuths rounded to 1 decimal number
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
    # Make it finer for the selcted scaling.
    # By adding the last constant we make sure that we do not start with
    # point (0,0).
    scale = 3.
    radii = np.arange((radius_pixels - min_radius) * scale + 1) / scale + min_radius

    # For each direction compute all possible horizont point position
    # and round them to integers
    for i in range(num_directions):
        x_int = np.round(x[i] * radii, decimals=0)
        y_int = np.round(y[i] * radii, decimals=0)
        # consider only the minimal number of points
        # use the trick with set and complex nuber as the input
        coord_complex = set(x_int + 1j * y_int)
        # to sort proportional with increasing radius, 
        # set has to be converted to numpy array
        shift_pairs = np.array([(k.real, k.imag) for k in coord_complex]).astype(int)
        distance = np.sqrt(np.sum(shift_pairs ** 2, axis=1))
        sort_index = np.argsort(distance)
        # write for each direction shifts and corresponding distances
        shift[angles[i]] = {
            "shift": [(k[0], k[1]) for k in shift_pairs[sort_index]],
            "distance": distance[sort_index],
        }

    return shift


def sky_view_factor_compute(height_arr,
                            radius_max=10,
                            radius_min=1,
                            num_directions=16,
                            compute_svf=True,
                            compute_opns=False,
                            compute_asvf=False,
                            a_main_direction=315.,
                            a_poly_level=4,
                            a_min_weight=0.4,
                            no_data=None,
                            fill_no_data=False,
                            keep_original_no_data=False
                            ):
    """
    Calculates horizon based visualizations: Sky-view factor, Anisotopic SVF and Openess.

    Parameters
    ----------
    height_arr : elevation (DEM) as 2D numpy array
    radius_max : maximal search radius in pixels/cells (not in meters)
    radius_min : minimal search radius in pixels/cells (not in meters), for noise reduction
    num_directions : number of directions as input
    compute_svf : if true it computes and outputs svf
    compute_asvf : if true it computes and outputs asvf
    compute_opns : if true it computes and outputs opns
    a_main_direction : main direction of anisotropy
    a_poly_level : level of polynomial that determines the anisotropy
    a_min_weight : weight to consider anisotropy:
                 0 - low anisotropy, 
                 1 - high  anisotropy (no illumination from the direction opposite the main direction)
    no_data : value that represents no_data, all pixels with this value are changed to np.nan
    fill_no_data : if True it fills where np.nan (no_data) with mean of surrounding pixels (3x3)
    keep_original_no_data : if True it changes all output pixels to np.nan where dem has no_data

    Returns
    -------
    {"svf": svf_out, "asvf": asvf_out, "opns": opns_out} : dictionary
        svf_out, skyview factor : 2D array of skyview factor.
        asvf_out, anisotropic skyview factor : 2D array of anisotropic skyview factor.
        opns_out, openness : 2D array openness (elevation angle of horizon)
    """
    # change no_data to np.nan
    idx_no_data = None
    if no_data is not None:
        if keep_original_no_data:  # save indexes where is no_data
            if np.isnan(no_data):
                idx_no_data = np.where(np.isnan(height_arr))
            else:
                idx_no_data = np.where(height_arr == no_data)
        height_arr[height_arr == no_data] = np.nan

    # fill no data with mean of surrounding pixels
    if fill_no_data:
        height_arr = fill_where_nan(height_arr)

    # pad the array for the radius_max on all 4 sides
    height = np.pad(height_arr, radius_max, mode='symmetric')

    # compute the vector of movement and corresponding distances
    move = horizon_shift_vector(num_directions=num_directions, radius_pixels=radius_max, min_radius=radius_min)

    # init the output for usual SVF
    if compute_svf:
        svf_out = np.zeros(height.shape, dtype=np.float32)
    else:
        svf_out = None
    # init the output for azimuth dependent SVF
    if compute_asvf:
        asvf_out = np.zeros(height.shape, dtype=np.float32)
        w_m = a_min_weight
        w_a = np.deg2rad(a_main_direction)
        weight = np.arange(num_directions) * (2 * np.pi / num_directions)
        weight = (1 - w_m) * (np.cos((weight - w_a) / 2)) ** a_poly_level + w_m
    else:
        asvf_out = None
    # init the output for Openess
    if compute_opns:
        opns_out = np.zeros(height.shape, dtype=np.float32)
    else:
        opns_out = None

        # search for horizon in each direction...
    for i_dir, direction in enumerate(move):
        # reset maximum at each iteration (direction)
        max_slope = np.zeros(height.shape, dtype=np.float32) - 1000

        # ... and to the search radius
        for i_rad, radius in enumerate(move[direction]["distance"]):
            # get shift index from move dictionary
            shift_indx = move[direction]["shift"][i_rad]
            # estimate the slope
            _ = (np.roll(height, shift_indx, axis=(0, 1)) - height) / radius
            # compare to the previus max slope and keep the larges
            max_slope = np.maximum(max_slope, _)

        # convert to angle in radians and compute directional output
        _ = np.arctan(max_slope)
        if compute_svf:
            svf_out = svf_out + (1 - np.sin(np.maximum(_, 0)))
        if compute_asvf:
            asvf_out = asvf_out + (1 - np.sin(np.maximum(_, 0))) * weight[i_dir]
        if compute_opns:
            opns_out = opns_out + _

    # cut to original extent and 
    # average the directional output over all directions
    if compute_svf:
        svf_out = svf_out[radius_max:-radius_max, radius_max:-radius_max] / num_directions
        # change result to np.nan where dem is no_data
        if no_data is not None and keep_original_no_data:
            svf_out[idx_no_data] = np.nan
    if compute_asvf:
        asvf_out = asvf_out[radius_max:-radius_max, radius_max:-radius_max] / np.sum(weight)
        # change result to np.nan where dem is no_data
        if no_data is not None and keep_original_no_data:
            asvf_out[idx_no_data] = np.nan
    if compute_opns:
        opns_out = np.rad2deg(0.5 * np.pi - (opns_out[radius_max:-radius_max, radius_max:-radius_max] / num_directions))
        # change result to np.nan where dem is no_data
        if no_data is not None and keep_original_no_data:
            opns_out[idx_no_data] = np.nan

    # return results within dict
    dict_svf_asvf_opns = {"svf": svf_out, "asvf": asvf_out, "opns": opns_out}
    dict_svf_asvf_opns = {k: v for k, v in dict_svf_asvf_opns.items() if v is not None}  # filter out none

    return dict_svf_asvf_opns


def sky_view_factor(dem,
                    resolution,
                    compute_svf=True,
                    compute_opns=False,
                    compute_asvf=False,
                    svf_n_dir=16,
                    svf_r_max=10,
                    svf_noise=0,
                    asvf_dir=315,
                    asvf_level=1,
                    ve_factor=1,
                    no_data=None,
                    fill_no_data=False,
                    keep_original_no_data=False
                    ):
    """
    Prepare the data, call sky_view_factor_compute, reformat and return back 2D arrays.

    Parameters
    ----------
    dem : input DEM 2D numpy array (Ve Exaggeration and pixel size already considered)
    compute_svf : compute SVF (True) or not (False)
    compute_opns : compute OPENNESS (True) or not (False)
    resolution : pixel resolution
    svf_n_dir : number of directions
    svf_r_max : maximal search radius in pixels
    svf_noise : the level of noise remove (0-don't remove, 1-low, 2-med, 3-high)
    compute_asvf : compute anisotropic SVF (True) or not (False)
    asvf_level : level of anisotropy, 1-low, 2-high,
    asvf_dir : dirction of anisotropy
    ve_factor : vertical exaggeration factor
    no_data : value that represents no_data, all pixels with this value are changed to np.nan
    fill_no_data : if True it fills where np.nan (no_data) with mean of surrounding pixels (3x3)
    keep_original_no_data : if True it changes all output pixels to np.nan where dem has no_data
        CONSTANTS:
            sc_asvf_min : level of polynomial that determines the anisotropy, selected with in_asvf_level
            sc_asvf_pol : level of polynomial that determines the anisotropy, selected with in_asvf_level
            sc_svf_r_min : the portion (percent) of the maximal search radius to ignore in horizon estimation;
                           for each noise level, selected with in_svf_noise

    Returns
    -------
    {"svf": svf_out, "asvf": asvf_out, "opns": opns_out} : dictionary
        svf_out, skyview factor : 2D numpy vector of skyview factor.
        asvf_out, anisotropic skyview factor : 2D numpy vector of anisotropic skyview factor.
        opns_out, openness : 2D numpy openness (elevation angle of horizon)
    """
    if dem.ndim != 2:
        raise Exception("rvt.vis.sky_view_factor: dem has to be 2D np.array!")
    if not (1000 >= ve_factor >= -1000):
        raise Exception("rvt.vis.sky_view_factor: ve_factor must be between -1000 and 1000!")
    if svf_noise != 0 and svf_noise != 1 and svf_noise != 2 and svf_noise != 3:
        raise Exception("rvt.vis.sky_view_factor: svf_noise must be one of the following values (0-don't remove, 1-low,"
                        " 2-med, 3-high)!")
    if asvf_level != 1 and asvf_level != 2:
        raise Exception("rvt.vis.sky_view_factor: asvf_leve must be one of the following values (1-low, 2-high)!")
    if not compute_svf and not compute_asvf and not compute_opns:
        raise Exception("rvt.vis.sky_view_factor: All computes are false!")
    if no_data is None and fill_no_data:
        warnings.warn("rvt.vis.sky_view_factor: In order to fill no data (fill_no_data = True) you have to input"
                      " no_data!")
    if (no_data is None or not fill_no_data) and keep_original_no_data:
        warnings.warn("rvt.vis.sky_view_factor: In order to keep original no data (keep_original_no_data = True)"
                      " you have to input no_data and fill_no_data has to be True!")

    # TODO: proper check of input data: DEM 2D nummeric array, resolution, max_radius....

    dem = dem.astype(np.float32)
    dem = dem * ve_factor

    # CONSTANTS
    # level of polynomial that determines the anisotropy, selected with in_asvf_level (1 - low, 2 - high)
    sc_asvf_pol = [4, 8]
    sc_asvf_min = [0.4, 0.1]
    # the portion (percent) of the maximal search radius to ignore in horizon estimation; for each noise level,
    # selected with in_svf_noise (0-3)
    sc_svf_r_min = [0., 10., 20., 40.]

    # pixel size
    dem = dem / resolution

    # minimal search radious depends on the noise level, it has to be an integer not smaller than 1
    svf_r_min = max(np.round(svf_r_max * sc_svf_r_min[svf_noise] * 0.01, decimals=0), 1)

    # set anisotropy parameters
    poly_level = sc_asvf_pol[asvf_level - 1]
    min_weight = sc_asvf_min[asvf_level - 1]

    dict_svf_asvf_opns = sky_view_factor_compute(height_arr=dem,
                                                 radius_max=svf_r_max,
                                                 radius_min=svf_r_min,
                                                 num_directions=svf_n_dir,
                                                 compute_svf=compute_svf,
                                                 compute_opns=compute_opns,
                                                 compute_asvf=compute_asvf,
                                                 a_main_direction=asvf_dir,
                                                 a_poly_level=poly_level,
                                                 a_min_weight=min_weight,
                                                 no_data=no_data,
                                                 fill_no_data=fill_no_data,
                                                 keep_original_no_data=keep_original_no_data
                                                 )

    return dict_svf_asvf_opns


def local_dominance(dem,
                    min_rad=10,
                    max_rad=20,
                    rad_inc=1,
                    angular_res=15,
                    observer_height=1.7,
                    ve_factor=1,
                    no_data=None,
                    fill_no_data=False,
                    keep_original_no_data=False
                    ):
    """
    Compute Local Dominance dem visualization.
    Adapted from original version that is part of the Lida Visualisation Toolbox LiVT developed by Ralf Hesse.

    Parameters
    ----------
    dem : input DEM 2D numpy array
    min_rad : minimum radial distance (in pixels) at which the algorithm starts with visualization computation
    max_rad : maximum radial distance (in pixels) at which the algorithm ends with visualization computation
    rad_inc : radial distance steps in pixels
    angular_res : angular step for determination of number of angular directions
    observer_height : height at which we observe the terrain
    ve_factor : vertical exaggeration factor
    no_data : value that represents no_data, all pixels with this value are changed to np.nan
    fill_no_data : if True it fills where np.nan (no_data) with mean of surrounding pixels (3x3)
    keep_original_no_data : if True it changes all output pixels to np.nan where dem has no_data

    Returns
    -------
    local_dom_out - 2D numpy array of local dominance
    """
    if dem.ndim != 2:
        raise Exception("rvt.vis.local_dominance: dem has to be 2D np.array!")
    if not (1000 >= ve_factor >= -1000):
        raise Exception("rvt.vis.local_dominance: ve_factor must be between -1000 and 1000!")
    if no_data is None and fill_no_data:
        warnings.warn("rvt.vis.local_dominance: In order to fill no data (fill_no_data = True) you have to input"
                      " no_data!")
    if (no_data is None or not fill_no_data) and keep_original_no_data:
        warnings.warn("rvt.vis.local_dominance: In order to keep original no data (keep_original_no_data = True)"
                      " you have to input no_data and fill_no_data has to be True!")

    # change no_data to np.nan
    idx_no_data = None
    if no_data is not None:
        if keep_original_no_data:  # save indexes where is no_data
            if np.isnan(no_data):
                idx_no_data = np.where(np.isnan(dem))
            else:
                idx_no_data = np.where(dem == no_data)
        dem[dem == no_data] = np.nan

    dem = dem.astype(np.float32)
    dem = dem * ve_factor

    # fill no data with mean of surrounding pixels
    if fill_no_data:
        dem = fill_where_nan(dem)

    # create a vector with possible distances
    n_dist = int((max_rad - min_rad) / rad_inc + 1)
    distances = np.arange(n_dist * rad_inc, step=rad_inc) + min_rad
    # create vector with possible angles
    n_ang = int(359 / angular_res + 1)
    angles = np.arange(n_ang * angular_res, step=angular_res)
    # determine total area within radius range
    norma = np.sum((observer_height / distances) * (2 * distances + rad_inc)) * n_ang

    # image shifts
    n_shifts = distances.size * angles.size
    x_t = (np.outer(np.cos(np.deg2rad(angles)), distances)).reshape(n_shifts)
    y_t = (np.outer(np.sin(np.deg2rad(angles)), distances)).reshape(n_shifts)
    distances = (np.outer(np.ones(n_ang), distances)).reshape(n_shifts)
    dist_factr = 2 * distances + rad_inc

    local_dom_out = dem * 0
    for i_s in range(n_shifts):
        dem_moved = np.roll(dem, int(round(y_t[i_s])), axis=0)
        dem_moved = np.roll(dem_moved, int(round(x_t[i_s])), axis=1)
        idx_lower = np.where((dem + observer_height) > dem_moved)
        if idx_lower[0].size > 0:
            local_dom_out[idx_lower[0], idx_lower[1]] = local_dom_out[idx_lower[0], idx_lower[1]] + \
                                                        (dem[idx_lower[0], idx_lower[1]] + observer_height -
                                                         dem_moved[idx_lower[0], idx_lower[1]]) / \
                                                        distances[i_s] * dist_factr[i_s]
    local_dom_out = local_dom_out / norma

    # change result to np.nan where dem is no_data
    if no_data is not None and keep_original_no_data:
        local_dom_out[idx_no_data] = np.nan

    return local_dom_out


def horizon_generate_coarse_dem(dem_fine,
                                pyramid_scale,
                                conv_from,
                                conv_to,
                                max_radius
                                ):
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
    # The final number of columns/lines has to fullfil:
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
    dem_fine = np.pad(dem_fine, ((-conv_from, conv_to), (-conv_from, conv_to)), mode="symmetric")

    # Convolution (keep maximum)
    dem_convolve = np.zeros(dem_fine.shape)
    for i in np.arange(pyramid_scale) + conv_from:
        for j in np.arange(pyramid_scale) + conv_from:
            dem_convolve = np.maximum(dem_convolve, np.roll(dem_fine, (i, j), axis=(0, 1)))
    # Divide by pyramid_scale to account for the chage of resolution
    # (important for the angle computation later on)
    dem_convolve = dem_convolve / pyramid_scale

    # Consider only the selceted convoluted points according to the scale change.
    # As we select slice's end point make sure to consider at least 1 point more 
    # to the right / below to really include it (Python way of considering end index).
    dem_coarse = dem_convolve[-conv_from:(n_lin_coarse * pyramid_scale + 1):pyramid_scale,
                 -conv_from:(n_col_coarse * pyramid_scale + 1):pyramid_scale]

    # Final padding to enable searching the horizon over the edge:
    # use constant-mode set to the minimal height, so it doesn't 
    # affect the horizon estimation.
    dem_coarse = np.pad(dem_coarse, ((max_radius, max_radius), (max_radius, max_radius)), mode="constant",
                        constant_values=dem_coarse.min())

    return dem_coarse


def horizon_generate_pyramids(dem,
                              num_directions=4,
                              max_fine_radius=100,
                              max_pyramid_radius=10,
                              pyramid_scale=5
                              ):
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
    # So you start in the level 1 with tmin_pyramid_radius=4
    # and you search from 4 to 10 distances (n_pyramid_radius=7)
    min_pyramid_radius = int(np.floor(max_pyramid_radius / pyramid_scale)) + 1
    n_pyramid_radius = max_pyramid_radius - min_pyramid_radius + 1

    # get the convolution window indices
    conv_to = int(np.floor(pyramid_scale / 2.))
    if (pyramid_scale % 2) == 0:
        conv_from = 1 - conv_to
    else:
        conv_from = -conv_to

    # initializations
    pyramid_levels = 0
    work = True
    pyramid = {}

    # Determine the number of levels and
    # the last radius to be used in the highest level.
    while work == True:
        _ = max_fine_radius / pyramid_scale ** pyramid_levels
        if _ > max_pyramid_radius:
            pyramid_levels = pyramid_levels + 1
        else:
            work = False
            last_radius = np.round(max_fine_radius / pyramid_scale ** pyramid_levels, decimals=0)

    # fill out the pyramid dict with the metadata required for horizont searching.
    for level in np.arange(pyramid_levels + 1):
        # the level 0 contains the other min_radius as the rest of levels
        if level == 0:
            min_radius = 1
            dem_fine = np.copy(np.pad(dem, max_pyramid_radius, mode="constant", constant_values=dem.min()))
        else:
            min_radius = min_pyramid_radius
            dem_fine = np.copy(dem_coarse)
        # the last level contains the other radius_pixels as the rest of levels
        if level == pyramid_levels:
            max_radius = last_radius
        else:
            max_radius = max_pyramid_radius
        # determine the dict of shifts
        shift = horizon_shift_vector(num_directions, max_radius, min_radius)
        dem_coarse = horizon_generate_coarse_dem(dem_fine, pyramid_scale, conv_from, conv_to, max_pyramid_radius)
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


def sky_illumination(dem,
                     resolution,
                     sky_model="overcast",
                     compute_shadow=True,
                     shadow_horizon_only=False,
                     max_fine_radius=100,
                     num_directions=32,
                     shadow_az=315,
                     shadow_el=35,
                     ve_factor=1,
                     no_data=None,
                     fill_no_data=False,
                     keep_original_no_data=False
                     ):
    """
    Compute topographic corrections for sky illumination.

    Parameters
    ----------
    shadow_horizon_only
    dem : numpy 2D array of elevation (DEM)
    resolution : dem pixel size
    sky_model : sky model, it can be 'overcast' or 'uniform'
    compute_shadow : bool compute shadow
    shadow_horizon_only : returns only dict {"shadow": shadow, "horizon": horizon}
    max_fine_radius : max shadow modeling distance [pixels]
    num_directions : number of directions to search for horizon
    shadow_az : shadow azimuth
    shadow_el : shadow elevation
    ve_factor : vertical exaggeration factor
    no_data : value that represents no_data, all pixels with this value are changed to np.nan
    fill_no_data : if True it fills where np.nan (no_data) with mean of surrounding pixels (3x3)
    keep_original_no_data : if True it changes all output pixels to np.nan where dem has no_data

    Returns
    -------
    sky_illum_out : 2D numpy result array
    """
    # standard pyramid settings
    pyramid_scale = 3
    max_pyramid_radius = 10

    if not (1000 >= ve_factor >= -1000):
        raise Exception("rvt.vis.sky_illumination: ve_factor must be between -1000 and 1000!")
    if shadow_az > 360 or shadow_az < 0:
        raise Exception("rvt.vis.sky_illumination: shadow_az must be between 0 and 360!")
    if shadow_el > 90 or shadow_el < 0:
        raise Exception("rvt.vis.sky_illumination: shadow_el must be between 0 and 90!")
    if no_data is None and fill_no_data:
        warnings.warn("rvt.vis.sky_illumination: In order to fill no data (fill_no_data = True) you have to input"
                      " no_data!")
    if (no_data is None or not fill_no_data) and keep_original_no_data:
        warnings.warn("rvt.vis.sky_illumination: In order to keep original no data (keep_original_no_data = True)"
                      " you have to input no_data and fill_no_data has to be True!")

    # change no_data to np.nan
    idx_no_data = None
    if no_data is not None:
        if keep_original_no_data:  # save indexes where is no_data
            if np.isnan(no_data):
                idx_no_data = np.where(np.isnan(dem))
            else:
                idx_no_data = np.where(dem == no_data)
        dem[dem == no_data] = np.nan

    dem = dem.astype(np.float32)
    dem = dem * ve_factor

    # fill no data with mean of surrounding pixels
    if fill_no_data:
        dem = fill_where_nan(dem)

    if sky_model.lower() == "overcast":
        compute_overcast = True
        compute_uniform = False
    elif sky_model.lower() == "uniform":
        compute_overcast = False
        compute_uniform = True
    else:
        raise Exception("rvt.vis.sky_illumination: sky_model must be overcast or uniform!")

    # generate slope and aspect
    _ = slope_aspect(np.pad(dem, max_pyramid_radius, mode="symmetric"), resolution, resolution)
    slope = _["slope"]
    aspect = _["aspect"]

    # build DEM pyramids
    pyramid = horizon_generate_pyramids(dem,
                                        num_directions=num_directions,
                                        max_fine_radius=max_fine_radius,
                                        max_pyramid_radius=max_pyramid_radius,
                                        pyramid_scale=pyramid_scale, )
    n_levels = np.max([i for i in pyramid])

    # get the convolution window indices
    conv_to = int(np.floor(pyramid_scale / 2.))
    if (pyramid_scale % 2) == 0:
        conv_from = 1 - conv_to
    else:
        conv_from = -conv_to
    # directional halve-resolution for integration limits
    da = np.pi / num_directions

    # init the intermediate results for uniform SI
    uniform_a = np.zeros((dem.shape[0] + 2 * max_pyramid_radius, dem.shape[1] + 2 * max_pyramid_radius),
                         dtype=np.float32)
    uniform_b = np.copy(uniform_a)
    # init the output for overcast SI
    if compute_overcast:
        overcast_out = np.zeros(dem.shape, dtype=np.float32)
        overcast_c = np.zeros((dem.shape[0] + 2 * max_pyramid_radius, dem.shape[1] + 2 * max_pyramid_radius),
                              dtype=np.float32)
        overcast_d = np.copy(overcast_c)
    else:
        overcast_out = None
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
    else:
        shadow_out = None
        horizon_out = None
        overcast_sh_out = None
        uniform_sh_out = None

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
                _ = (np.roll(height, shift_indx, axis=(0, 1)) - height) / radius
                # compare to the previus max slope and keep the larges
                max_slope = np.maximum(max_slope, _)

            # resample the max_slope to a lower pyramid level
            if i_level > 0:
                lin_fine = pyramid[i_level - 1]["i_lin"] + (
                        conv_from + max_pyramid_radius * pyramid_scale - max_pyramid_radius)
                col_fine = pyramid[i_level - 1]["i_col"] + (
                        conv_from + max_pyramid_radius * pyramid_scale - max_pyramid_radius)
                lin_coarse = pyramid[i_level]["i_lin"] * pyramid_scale
                col_coarse = pyramid[i_level]["i_col"] * pyramid_scale
                interp_spline = scipy.interpolate.RectBivariateSpline(lin_coarse, col_coarse, max_slope)
                max_slope = interp_spline(lin_fine, col_fine)

        # convert to angle in radians and compute directional output
        _ = np.arctan(max_slope)
        uniform_a = uniform_a + (np.cos(_)) ** 2
        _d_aspect = np.sin((dir_rad - da) - aspect) - np.sin((dir_rad + da) - aspect)
        uniform_b = uniform_b + _d_aspect * (np.pi / 4. - _ / 2. - np.sin(2. * _) / 4.)
        if compute_overcast:
            _cos3 = (np.cos(_)) ** 3
            overcast_c = overcast_c + _cos3
            overcast_d = overcast_d + _d_aspect * (2. / 3. - np.cos(_) + _cos3 / 3.)
        if compute_shadow and (direction == shadow_az):
            horizon_out = np.degrees(_[max_pyramid_radius:-max_pyramid_radius, max_pyramid_radius:-max_pyramid_radius])
            shadow_out = (horizon_out < shadow_el) * 1
            if shadow_horizon_only:
                # change result to np.nan where dem is no_data
                if no_data is not None and keep_original_no_data:
                    shadow_out[idx_no_data] = np.nan
                    horizon_out[idx_no_data] = np.nan
                return {"shadow": shadow_out, "horizon": horizon_out}

    # because of numeric stabilty check if the uniform_b is less then pi and greater than 0
    uniform_out = (da) * np.cos(slope) * uniform_a + np.sin(slope) * np.minimum(np.maximum(uniform_b, 0), np.pi)
    uniform_out = uniform_out[max_pyramid_radius:-max_pyramid_radius, max_pyramid_radius:-max_pyramid_radius]
    if compute_overcast:
        # because of numeric stabilty check if the uniform_b is less then pi and greater than 0
        overcast_out = (2. * da / 3.) * np.cos(slope) * overcast_c + np.sin(slope) * np.maximum(overcast_d, 0)
        overcast_out = overcast_out[max_pyramid_radius:-max_pyramid_radius, max_pyramid_radius:-max_pyramid_radius]
        overcast_out = 0.33 * uniform_out + 0.67 * overcast_out
        overcast_out = overcast_out / overcast_out.max()
    if compute_shadow:
        uniform_sh_out = (0.8 * uniform_out + 0.2 * shadow_out)
        if compute_overcast:
            overcast_sh_out = (0.8 * overcast_out + 0.2 * shadow_out)

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
        # change result to np.nan where dem is no_data
        if no_data is not None and keep_original_no_data:
            uniform_out[idx_no_data] = np.nan
        return uniform_out
    elif compute_uniform and compute_shadow:
        # change result to np.nan where dem is no_data
        if no_data is not None and keep_original_no_data:
            uniform_sh_out[idx_no_data] = np.nan
        return uniform_sh_out
    elif compute_overcast and not compute_shadow:
        # change result to np.nan where dem is no_data
        if no_data is not None and keep_original_no_data:
            overcast_out[idx_no_data] = np.nan
        return overcast_out
    elif compute_overcast and compute_shadow:
        # change result to np.nan where dem is no_data
        if no_data is not None and keep_original_no_data:
            overcast_sh_out[idx_no_data] = np.nan
        return overcast_sh_out


def shadow_horizon(dem,
                   resolution,
                   shadow_az=315,
                   shadow_el=35,
                   ve_factor=1,
                   no_data=None,
                   fill_no_data=False,
                   keep_original_no_data=False
                   ):
    """
    Compute shadow and horizon.

    Parameters
    ----------
    dem : numpy 2D array of elevation (DEM)
    resolution : raster resolution
    shadow_az : shadow azimuth
    shadow_el : shadow elevation
    ve_factor : vertical exaggeration factor
    no_data : value that represents no_data, all pixels with this value are changed to np.nan
    fill_no_data : if True it fills where np.nan (no_data) with mean of surrounding pixels (3x3)
    keep_original_no_data : if True it changes all output pixels to np.nan where dem has no_data

    Returns
    -------
    {"shadow": shadow 2D np.array, "horizon": horizon 2D np.array}
    """
    if not (1000 >= ve_factor >= -1000):
        raise Exception("rvt.vis.shadow_horizon: ve_factor must be between -1000 and 1000!")
    if shadow_az > 360 or shadow_az < 0:
        raise Exception("rvt.vis.shadow_horizon: shadow_az must be between 0 and 360!")
    if shadow_el > 90 or shadow_el < 0:
        raise Exception("rvt.vis.shadow_horizon: shadow_el must be between 0 and 90!")
    if no_data is None and fill_no_data:
        warnings.warn("rvt.vis.shadow_horizon: In order to fill no data (fill_no_data = True) you have to input"
                      " no_data!")
    if (no_data is None or not fill_no_data) and keep_original_no_data:
        warnings.warn("rvt.vis.shadow_horizon: In order to keep original no data (keep_original_no_data = True)"
                      " you have to input no_data and fill_no_data has to be True!")

    return sky_illumination(dem=dem, resolution=resolution, shadow_horizon_only=True, shadow_el=shadow_el,
                            shadow_az=shadow_az, ve_factor=ve_factor, no_data=no_data, fill_no_data=fill_no_data,
                            keep_original_no_data=keep_original_no_data)


def fill_where_nan(dem):
    """Replaces np.nan values. It takes nan surrounding (neighbor) cells (3x3), if number of nans <= max_nan_nr it
    calculates mean of not nans and writes mean to the center nan pixel.
    Function iterates and first calculates areas where is less nan values in the surrounding (neighbour) cells."""
    dem_out = np.copy(dem)
    rows = dem_out.shape[0]
    cols = dem_out.shape[1]
    max_nan_nr = 4  # maximum number of nans in surrounding array (3x3) to calculate mean
    for nr_limit_nan in range(2, max_nan_nr + 1, 1):  # iterate
        changes = True
        while changes:
            idx_nans = np.where(np.isnan(dem_out))  # index of all nans
            changes = False
            for i_nan in range(len(idx_nans[0])):
                row_idx_nan = idx_nans[0][i_nan]
                col_idx_nan = idx_nans[1][i_nan]
                if row_idx_nan == 0 or col_idx_nan == 0 or row_idx_nan == rows or col_idx_nan == cols:
                    continue  # skip if edge
                surr_arr = dem_out[row_idx_nan - 1:row_idx_nan + 2,
                           col_idx_nan - 1:col_idx_nan + 2]  # surrounding array 3x3
                nr_surr_nan = np.count_nonzero(np.isnan(surr_arr))  # number of nans in surrounding array
                if nr_surr_nan <= max_nan_nr:
                    dem_out[row_idx_nan, col_idx_nan] = np.nanmean(surr_arr)
                    changes = True
    return dem_out
