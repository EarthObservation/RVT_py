"""
NAME:
    rvt_py, rvt.vis - rvt visualization functions

DESCRIPTION:
    Contains all functions for visualization. Functions are rewritten from RVT (IDL).
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
import scipy.ndimage


# TODO: check speed sky_ilumination
# TODO: check IDL vectorisation remains:
#  (sky_view_factor, sky_view_det_move, sky_view_compute, morph_shade_move, morph_shade, sky_illumination)


def byte_scale(data, c_min=None, c_max=None, high=255, low=0):
    """
    Remade old scipy function.
    Byte scales an array (image).

    Byte scaling means converting the input image to uint8 dtype and scaling
    the range to ``(low, high)`` (default 0-255).

    Parameters
    ----------
    data : numpy image data array.
    c_min : scalar, Bias scaling of small values. Default is ``data.min()``.
    c_max : scalar, Bias scaling of large values. Default is ``data.max()``.
    high : scalar, Scale max value to `high`.  Default is 255.
    low : scalar, Scale min value to `low`.  Default is 0.

    Returns
    -------
    img_array : uint8 ndarray
        The byte-scaled array.
    """

    if high < low:
        raise ValueError("`high` should be larger than `low`.")

    if c_min is None:
        c_min = data.min()
    if c_max is None:
        c_max = data.max()

    c_scale = c_max - c_min
    if c_scale < 0:
        raise ValueError("`cmax` should be larger than `cmin`.")
    elif c_scale == 0:
        c_scale = 1

    if data.dtype == np.uint8:
        byte_data = (high + 1) * (data - c_min - 1) / (c_max - c_min)  # copied from IDL BYTSCL
        byte_data[byte_data > high] = high
        byte_data[byte_data < 0] = 0
        return np.cast[np.uint8](byte_data) + np.cast[np.uint8](low)

    # scale = float(high - low) / cscale  # old scipy fn
    # byte_data = (data * 1.0 - cmin) * scale + 0.4999  # old scipy fn

    byte_data = (high + 0.9999) * (data - c_min) / (c_max - c_min)  # copied from IDL BYTSCL
    byte_data[byte_data > high] = high
    byte_data[byte_data < 0] = 0
    return np.cast[np.uint8](byte_data) + np.cast[np.uint8](low)


def slope_aspect(dem, resolution_x, resolution_y, ve_factor=1, output_units="radian"):
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
    ve_factor : vertical exaggeration factor (must be greater than 0)
    output_units : percent, degree, radians

    Returns
    -------
    {"slope": slope_out, "aspect": aspect_out} : dictionaries with 2D numpy arrays
    """

    if ve_factor <= 0:
        raise Exception("rvt.vis.slope_aspect: ve_factor must be a positive number!")

    if resolution_x < 0 or resolution_y < 0:
        raise Exception("rvt.vis.slope_aspect: resolution must be a positive number!")

    dem = dem.astype(np.float32)
    if ve_factor != 1:
        dem = dem * ve_factor

    # add frame of 0 (additional row up bottom and column left right)
    dem = np.pad(dem, pad_width=1, mode="constant", constant_values=0)

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

    # remove the frame (padding)
    slope_out = slope_out[1:-1, 1:-1]
    aspect_out = aspect_out[1:-1, 1:-1]

    # edges to -1
    slope_out[:, 0] = -1
    slope_out[0, :] = -1
    slope_out[:, -1] = -1
    slope_out[-1, :] = -1
    aspect_out[:, 0] = -1
    aspect_out[0, :] = -1
    aspect_out[:, -1] = -1
    aspect_out[-1, :] = -1

    return {"slope": slope_out, "aspect": aspect_out}


def hillshade(dem, resolution_x, resolution_y, sun_azimuth=315, sun_elevation=35,
              slope=None, aspect=None, ve_factor=1):
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
    ve_factor : vertical exaggeration factor (must be greater than 0)

    Returns
    -------
    hillshade_out : result numpy array
    """
    if ve_factor <= 0:
        raise Exception("rvt.vis.slope_aspect: ve_factor must be a positive number!")

    if sun_azimuth > 360 or sun_elevation > 90 or sun_azimuth < 0 or sun_elevation < 0:
        raise Exception("rvt.vis.analytical_hillshading: sun_azimuth must be [0-360] and sun_elevation [0-90]!")

    if resolution_x < 0 or resolution_y < 0:
        raise Exception("rvt.vis.analytical_hillshading: resolution must be a positive number!")

    dem = dem.astype(np.float32)
    if ve_factor != 1:
        dem = dem * ve_factor

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

    return hillshade_out


def multi_hillshade(dem, resolution_x, resolution_y, nr_directions=16, sun_elevation=35,
                    slope=None, aspect=None, ve_factor=1):
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
    ve_factor : vertical exaggeration factor (must be greater than 0)

    Returns
    -------
    multi_hillshade_out : numpy array containing numpy_arrays of hillshades in different directions
    """

    if sun_elevation > 90 or sun_elevation < 0:
        raise Exception("rvt.vis.multiple_directions_hillshading: sun_elevation must be [0-90]!")

    if resolution_x < 0 or resolution_y < 0:
        raise Exception("rvt.vis.multiple_directions_hillshading: resolution must be a positive number!")

    if nr_directions < 1:
        raise Exception("rvt.vis.multiple_directions_hillshading: nr_directions must be a positive number!")

    if ve_factor <= 0:
        raise Exception("rvt.vis.slope_aspect: ve_factor must be a positive number!")

    dem = dem.astype(np.float32)
    if ve_factor != 1:
        dem = dem * ve_factor

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
        hillshades_arr_list.append(hillshading)
    multi_hillshade_out = np.asarray(hillshades_arr_list)

    return multi_hillshade_out


def slrm(dem, radius_cell=20, ve_factor=1):
    """
    Calculates Simple local relief model.

    Parameters
    ----------
    dem : input DEM 2D numpy array
    radius_cell : Radius for trend assessment [pixels]
    ve_factor : vertical exaggeration factor (must be greater than 0)

    Returns
    -------
    slrm_out : slrm 2D numpy array
    """
    if radius_cell < 10 or radius_cell > 50:
        raise Exception("rvt.vis.slrm: Radius for trend assessment needs to be in interval 10-50 pixels!")

    if ve_factor <= 0:
        raise Exception("rvt.vis.slope_aspect: ve_factor must be a positive number!")

    dem = dem.astype(np.float32)
    if ve_factor != 1:
        dem = dem * ve_factor

    dem[dem < -1200] = np.float64(np.NaN)
    dem[dem > 2000] = np.float64(np.NaN)

    # mean filter
    slrm_out = dem - scipy.ndimage.uniform_filter(dem, mode='nearest', size=radius_cell * 2)

    return slrm_out


def azimuth(xa, ya, xb, yb):
    """
    Determine the azimuth in the range of [0,2pi).

    Parameters
    ----------
    xa, ya, xb, yb (point coordinates)

    Returns
    -------
    a : outputs the azimuth in radians
    """
    north = ya - yb
    east = xb - xa
    if north == 0:
        if east > 0:
            a = np.pi / 2
        else:
            if east < 0:
                a = 3 * np.pi / 2
            else:
                a = np.float64(np.NaN)
    else:
        a0 = np.arctan(east / north)
        if north > 0 and east >= 0:
            a = a0
        elif north < 0:
            a = a0 + np.pi
        else:
            a = a0 + 2 * np.pi
    return a


def sky_view_det_move(num_directions, radius_cell, n_col):
    """
    Calculates Sky-View determination movement matrix.

    Parameters
    ----------
    num_directions : number of directions as input
    radius_cell : radius to consider in cells (not in meters)
    n_col : number of columns of the input DEM

    Returns
    -------
    move : 2D numpy array movement matrix
    """

    # parameters
    look_angle = 2 * np.pi / num_directions

    # matrix initialization
    move = np.zeros((3, radius_cell + 1, num_directions), dtype=np.float32)

    # for each direction
    for i_direction in range(num_directions):
        angle = i_direction * look_angle
        d = 0
        x0 = 0
        y0 = 0
        xt = x0
        yt = y0
        rad = 0

        if 0 <= angle < np.pi / 2:
            quad = 1
        elif np.pi / 2 <= angle < np.pi:
            quad = 2
        elif np.pi <= angle < 3 * np.pi / 2:
            quad = 3
        elif 3 * np.pi / 2 <= angle < 2 * np.pi:
            quad = 4

        # while within range
        while d <= radius_cell:
            # compute direction
            if quad == 1:
                # Right
                xa = xt + 1
                ya = yt
                # Up
                xb = xt
                yb = yt - 1
                # diagonal right up
                xc = xt + 1
                yc = yt - 1
            elif quad == 2:
                # Right
                xa = xt + 1
                ya = yt
                # Diagonal right down
                xb = xt + 1
                yb = yt + 1
                # Down
                xc = xt
                yc = yt + 1
            elif quad == 3:
                # Left
                xa = xt - 1
                ya = yt
                # Diagonal left down
                xb = xt - 1
                yb = yt + 1
                # Down
                xc = xt
                yc = yt + 1
            elif quad == 4:
                # Left
                xa = xt - 1
                ya = yt
                # Up
                xb = xt
                yb = yt - 1
                # Diagonal left up
                xc = xt - 1
                yc = yt - 1

            # azimuths of possible movements (nearest neighbor, no interpolation)
            k_a = azimuth(x0, y0, xa, ya)
            k_b = azimuth(x0, y0, xb, yb)
            k_c = azimuth(x0, y0, xc, yc)

            # Minimum difference in angle for new point
            if abs(k_a - angle) <= abs(k_b - angle):
                if abs(k_a - angle) <= abs(k_c - angle):
                    xt = xa
                    yt = ya
                else:
                    xt = xc
                    yt = yc
            else:
                if abs(k_b - angle) <= abs(k_c - angle):
                    xt = xb
                    yt = yb
                else:
                    xt = xc
                    yt = yc

            # Output
            move[0, rad, i_direction] = xt - x0
            move[1, rad, i_direction] = yt - y0
            d = np.sqrt((xt - x0) ** 2 + (yt - y0) ** 2)
            move[2, rad, i_direction] = d

            # next cell
            rad += 1

    # reformat the radius:
    # first row tells you, how many valid cells are available, bellow is actual radius in cells
    move[2, 1:radius_cell + 1, :] = move[2, 0:radius_cell, :]
    for i_direction in range(num_directions):
        tmp = move[2, 1:radius_cell + 1, i_direction]
        if tmp[tmp > radius_cell].size > 0:
            move[2, 0, i_direction] = np.min(tmp[tmp > radius_cell])
        else:
            move[2, 0, i_direction] = radius_cell

    # convert 2D index into 1D index
    move_t = np.zeros((2, radius_cell + 1, num_directions), dtype=np.float32)
    move_t[0, :, :] = move[1, :, :] * n_col + move[0, :, :]
    move_t[1, :, :] = move[2, :, :]
    move = move_t

    return move


def sky_view_factor_compute(height_arr, i_valid, radius_max, radius_min, num_directions,
                            a_main_direction, a_poly_level, a_min_weight, compute_svf=True, compute_asvf=False,
                            compute_opns=False):
    """
    Calculates Sky-View 1D vectors.

    Parameters
    ----------
    height_arr : elevation (DEM) as 2D numpy array
    i_valid : index of valid pixels to be processed
    radius_max : maximal search radius in pixels/cells (not in meters)
    radius_min : minimal search radius in pixels/cells (not in meters), for noise reduction
    num_directions : number of directions as input
    a_main_direction : main direction of anisotropy
    a_poly_level : level of polynomial that determines the anisotropy
    a_min_wight : weight to consider anisotropy (0 - isotropic, 1 - no illumination from the direction
                                                             opposite the main direction)
    compute_svf : if true it computes and outputs svf
    compute_asvf : if true it computes and outputs asvf
    compute_opns : if true it computes and outputs opns

    Returns
    -------
    {"svf": svf_out, "asvf": asvf_out, "opns": opns_out} : dictionary
        svf_out, skyview factor : 1D numpy vector of skyview factor.
        asvf_out, anisotropic skyview factor : 1D numpy vector of anisotropic skyview factor.
        opns_out, openness : 1D numpy openness (elevation angle of horizon)
    """

    height = height_arr

    svf_out = None  # if compute_svf is False then function doesn't compute svf and returns False
    asvf_out = None  # if compute_asvf is False then function doesn't compute asvf and returns False
    opns_out = None  # if compute_opns is False then function doesn't compute opns and returns False

    # directional step
    dir_step = 2 * np.pi / num_directions

    # vector of movement
    n_col = height.shape[1]  # number of columns
    count_height = i_valid[0].size  # number of all elements
    move = sky_view_det_move(num_directions=num_directions, radius_cell=radius_max, n_col=n_col)

    # init the outputs
    if compute_svf:
        svf_out = np.zeros(count_height, dtype=np.float32)
    if compute_opns:
        opns_out = np.zeros(count_height, dtype=np.float32)
    if compute_asvf:
        asvf_out = np.zeros(count_height, dtype=np.float32)
        w_m = a_min_weight  # compute weights
        w_a = np.deg2rad(a_main_direction)
        weight = np.arange(num_directions) * dir_step
        weight = (1 - w_m) * (np.cos((weight - w_a) / 2)) ** a_poly_level + w_m

    # look into each direction...
    for i_dir in range(num_directions):
        # reset maximum at each iteration - at each new direction
        max_slope = np.zeros(count_height, dtype=np.float32) - 1000

        # ... and to the search radius - this depends on the direction - radius is written in the first row
        for i_rad in range(1, int(move[1, 0, i_dir])):
            # TODO: Krištof
            # ignore radios if smaller than minimal defined radius
            if radius_min > move[1, i_rad, i_dir]:
                continue
            # search for max (sky)
            h_flt = height.flatten()
            np.seterr(divide='ignore', invalid='ignore')  # ignore warnings for dividing with zero
            m_slp = (h_flt[i_valid[0] + int(move[0, int(i_rad - 1), i_dir])] - h_flt[i_valid[0]]) / move[
                1, i_rad, i_dir]
            np.seterr(divide='warn', invalid='warn')  # reset warnings
            if np.isnan(m_slp).any():  # skip non existing m_slp
                continue
            np.clip(a=max_slope, a_min=m_slp, a_max=None, out=max_slope)

        max_slope = np.arctan(max_slope)

        if compute_opns:
            opns_out = opns_out + max_slope
        if compute_svf:
            svf_out = svf_out + (1 - np.sin(np.clip(a=max_slope, a_min=0, a_max=None)))
        if compute_asvf:
            asvf_out = asvf_out + (1 - np.sin(np.clip(a=max_slope, a_min=0, a_max=None))) * weight[i_dir]

    # Normalize to the number of directions / weights
    if compute_svf:
        svf_out = svf_out / num_directions
    if compute_asvf:
        asvf_out = asvf_out / np.sum(weight)
    if compute_opns:
        opns_out = np.pi / 2 - (opns_out / num_directions)

    dict_svf_asvf_opns = {"svf": svf_out, "asvf": asvf_out, "opns": opns_out}
    dict_svf_asvf_opns = {k: v for k, v in dict_svf_asvf_opns.items() if v is not None}  # filter out none

    return dict_svf_asvf_opns


def sky_view_factor(dem, resolution, compute_svf=True, compute_opns=False, compute_asvf=False,
                    svf_n_dir=16, svf_r_max=10, svf_noise=0, asvf_dir=315, asvf_level=1, ve_factor=1):
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
    a_min_weight : weight to consider anisotropy (0 - isotropic, 1 - no illumination from the direction
                   opposite the main direction)
    ve_factor : vertical exaggeration factor (must be greater than 0)
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

    if ve_factor <= 0:
        raise Exception("rvt.vis.slope_aspect: ve_factor must be a positive number!")

    dem = dem.astype(np.float32)
    if ve_factor != 1:
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
    # increase edge - visualization foes to edge, so mirror them, leave blank corners
    n_col = dem.shape[1]
    n_lin = dem.shape[0]
    tmp = np.zeros((n_lin + 2 * svf_r_max, n_col + 2 * svf_r_max), dtype=np.float32)

    # prepare indx for output
    tmp[svf_r_max:n_lin + svf_r_max, svf_r_max:n_col + svf_r_max] = 1

    indx_all = np.where(tmp.flatten() == 1)

    # fill center
    tmp[svf_r_max:n_lin + svf_r_max, svf_r_max:n_col + svf_r_max] = dem

    # final dem with increased edges
    dem = tmp
    del tmp

    # minimal search radious depends on the noise level
    svf_r_min = svf_r_max * sc_svf_r_min[svf_noise] * 0.01

    # set anisotropy parameters
    poly_level = sc_asvf_pol[asvf_level - 1]
    min_weight = sc_asvf_min[asvf_level - 1]

    svf_out = None
    asvf_out = None
    opns_out = None

    dict_svf_asvf_opns = sky_view_factor_compute(height_arr=dem, i_valid=indx_all, radius_max=svf_r_max,
                                                 radius_min=svf_r_min, num_directions=svf_n_dir,
                                                 a_main_direction=asvf_dir, a_poly_level=poly_level,
                                                 a_min_weight=min_weight, compute_svf=compute_svf,
                                                 compute_asvf=compute_asvf, compute_opns=compute_opns)
    # reshape 1D to 2D
    if compute_svf:
        svf_out = dict_svf_asvf_opns["svf"].reshape((n_lin, n_col))
    if compute_asvf:
        asvf_out = dict_svf_asvf_opns["asvf"].reshape((n_lin, n_col))
    if compute_opns:
        opns_out = dict_svf_asvf_opns["opns"].reshape((n_lin, n_col))
        opns_out = np.rad2deg(opns_out)

    dict_svf_asvf_opns = {"svf": svf_out, "asvf": asvf_out, "opns": opns_out}
    dict_svf_asvf_opns = {k: v for k, v in dict_svf_asvf_opns.items() if v is not None}  # filter out none

    return dict_svf_asvf_opns


def morph_shade_move(d_max, angle):
    """
    Calculates morph_shade movement matrix for sky illumination.

    Parameters
    ----------
    d_max : maximum search distance in pixels
    angle

    Returns
    -------
    move
    """

    # init
    move = np.zeros((d_max + 1, 3), dtype=np.float32)
    d = 0
    x0 = 0
    y0 = 0
    xt = x0
    yt = y0
    i_rad = 0
    quad = 0

    # determine quadrant number
    if 0 <= angle < np.pi / 2:
        quad = 1
    elif np.pi / 2 <= angle < np.pi:
        quad = 2
    elif np.pi <= angle < 3 * np.pi / 2:
        quad = 3
    elif 3 * np.pi / 2 <= angle < 2 * np.pi:
        quad = 4

    # while within range
    while d <= d_max:
        # compute direction
        if quad == 1:
            # Right
            xa = xt + 1
            ya = yt
            # Up
            xb = xt
            yb = yt - 1
            # diagonal right up
            xc = xt + 1
            yc = yt - 1
        elif quad == 2:
            # Right
            xa = xt + 1
            ya = yt
            # Diagonal right down
            xb = xt + 1
            yb = yt + 1
            # Down
            xc = xt
            yc = yt + 1
        elif quad == 3:
            # Left
            xa = xt - 1
            ya = yt
            # Diagonal left down
            xb = xt - 1
            yb = yt + 1
            # Down
            xc = xt
            yc = yt + 1
        elif quad == 4:
            # Left
            xa = xt - 1
            ya = yt
            # Up
            xb = xt
            yb = yt - 1
            # Diagonal left up
            xc = xt - 1
            yc = yt - 1

        # azimuths of possible movements (nearest neighbor, no interpolation)
        k_a = azimuth(x0, y0, xa, ya)
        k_b = azimuth(x0, y0, xb, yb)
        k_c = azimuth(x0, y0, xc, yc)

        # Minimum difference in angle for new point
        if abs(k_a - angle) <= abs(k_b - angle):
            if abs(k_a - angle) <= abs(k_c - angle):
                xt = xa
                yt = ya
            else:
                xt = xc
                yt = yc
        else:
            if abs(k_b - angle) <= abs(k_c - angle):
                xt = xb
                yt = yb
            else:
                xt = xc
                yt = yc

        # Output
        # TODO: Krištof
        # This is just a function to determine the direction of horizon search.
        # It is called only once, so there is no need to make any optimization here.
        move[i_rad, 0] = xt - x0
        move[i_rad, 1] = yt - y0
        d = np.sqrt((xt - x0) ** 2 + (yt - y0) ** 2)
        move[i_rad, 2] = d

        # next cell
        i_rad += 1

    move = move[0:i_rad, :]
    return move


def morph_shade(height, sol_z, sol_a, d_max, nrows, ncols, resolution):
    """
    Compute topographic corrections for sky illumination.

    Parameters
    ----------
    height : elevation 2D np array
    sol_z : solar zenith angle in rad (0 for vertical and pi/2 for horizontal surface)
    sol_a : solar azimuth angle
    d_max : maximum search distance in pixels
    resolution : pixel size

    Returns
    -------
    mask : determines those cells that are in its own (hillshade) or thrown (cast shade) shadow
    """

    # initialize the results
    mask = np.zeros((nrows + 2 * d_max, ncols + 2 * d_max), dtype=np.float32)
    mask[d_max:(nrows + d_max), d_max:(ncols + d_max)] = 1
    i_valid = np.where(mask.flatten() == 1)
    tmp = np.zeros((nrows + 2 * d_max, ncols + 2 * d_max), dtype=np.float32)
    tmp[d_max:(nrows + d_max), d_max:(ncols + d_max)] = height
    height = tmp
    del tmp
    # TODO: Krištof
    # determine the direction of movment
    move = morph_shade_move(d_max, sol_a)
    move_s = move.shape
    move1di = move[:, 1] * (ncols + 2 * d_max) + move[:, 0]
    move1dd = move[:, 2]
    # set the maximal allowed horizon angle (if it is greater, then the area in the shadow)
    max_slope = 0
    for i_rad in range(int(move_s[0])):
        height_flt = height.flatten()
        # m_slp = ((height_flt[i_valid[0]+int(move1di[i_rad])] - height_flt[i_valid[0]]) / move1dd[i_rad])
        sel = i_valid[0] + int(move1di[i_rad])
        nr_zeros = sel[sel > height_flt.size].size  # can't call indexes that are bigger than num of elements
        zeros = np.zeros(nr_zeros, dtype=np.float32)
        sel = sel[sel <= height_flt.size]
        m_slp = ((height_flt[sel] - height_flt[i_valid[0][:sel.size]]) / move1dd[i_rad])
        m_slp = np.append(m_slp, zeros)  # add zeros, for non existing indexes
        max_slope = np.clip(a=max_slope, a_min=m_slp, a_max=None)

    # update mask
    max_slope = np.arctan(max_slope / resolution)
    indx_mask = np.where(max_slope > (np.pi / 2 - sol_z))

    if indx_mask[0].size > 0:
        mask_size = mask.shape
        mask_flt = mask.flatten()
        mask_flt[i_valid[0][indx_mask[0]]] = 0
        mask = mask_flt.reshape((mask_size[0], mask_size[1]))
    mask = mask[d_max:(nrows + d_max), d_max:(ncols + d_max)]
    return mask


def sky_illumination(dem, resolution, sky_model="overcast", sampling_points=250, shadow_dist=100,
                     shadow_az=315, shadow_el=35, shadow_only=False, ve_factor=1):
    """
    Compute topographic corrections for sky illumination.

    Parameters
    ----------
    dem : numpy 2D array of elevation (DEM)
    resolution : dem pixel size
    sky_model : sky model [overcast, uniform]
    sampling_points : number of sampling points
    shadow_dist : max shadow modeling distance [pixels]
    shadow_az : shadow azimuth
    shadow_el : shadow elevation
    shadow_only : bool compute shadow only
    ve_factor : vertical exaggeration factor (must be greater than 0)

    Returns
    -------
    sky_illum_out : 2D numpy result array
    """

    if sampling_points != 250 and sampling_points != 500:
        raise Exception("rvt.vis.sky_illumination: sampling_points needs to be 250 or 500!")
    if sky_model != "overcast" and sky_model != "uniform":
        raise Exception("rvt.vis.sky_illumination: sky_model needs to be overcast or uniform!")
    if ve_factor <= 0:
        raise Exception("rvt.vis.slope_aspect: ve_factor must be a positive number!")

    dem = dem.astype(np.float32)
    if ve_factor != 1:
        dem = dem * ve_factor

    indx_no_values = np.where(dem < 0)
    dem[indx_no_values[0], indx_no_values[1]] = np.float64(np.NaN)

    # determine the max search distance in pixels
    h_min = np.nanmin(dem)
    h_max = np.nanmax(dem)
    dh = h_max - h_min
    dem_size = dem.shape
    dem[indx_no_values[0], indx_no_values[1]] = 0

    if shadow_az and shadow_el:
        sh_z = np.pi / 2 - np.deg2rad(shadow_el)
        sh_az = np.deg2rad(shadow_az)
        d_max = 1
        d_max = round(d_max * dh * np.tan(sh_z) / resolution)
        dem_tmp = dem
        out_shadow = morph_shade(dem_tmp, sh_z, sh_az, d_max, dem_size[0], dem_size[1], resolution)
        # out_shadow = adams_shadows(in_array=dem, az=shadow_az, alt=shadow_el, res=resolution, overlap=d_max)

    if shadow_only:
        dem[indx_no_values[0], indx_no_values[1]] = np.float64(np.NaN)
        return out_shadow

    else:
        dat_hillset = open(r'settings\{}_{}sp.txt'.format(sky_model, sampling_points), 'r')
        sky_illum_out = np.zeros((dem_size[0], dem_size[1]), dtype=np.float32)
        dict_slope_aspect = slope_aspect(dem=dem, resolution_x=resolution,
                                         resolution_y=resolution,
                                         ve_factor=1, output_units="radian")
        slope = dict_slope_aspect["slope"]
        aspect = dict_slope_aspect["aspect"]

        for line in dat_hillset:
            if line.strip() == "":  # empty line
                continue
            print(line)
            d_max = 1
            line = line.rstrip().split(",")
            azim = int(line[0])
            elev = int(line[1])
            weight = float(line[2])
            hillshade_tmp = hillshade(dem=dem, resolution_x=resolution, resolution_y=resolution,
                                      sun_azimuth=azim, sun_elevation=elev, slope=slope, aspect=aspect)
            sh_z = np.pi / 2 - np.deg2rad(elev)
            sh_az = np.deg2rad(azim)
            d_max = round(d_max * dh * np.tan(sh_z) / resolution)

            if d_max > 1:
                if shadow_dist != 'unlimited':
                    if d_max < int(shadow_dist):
                        d_max = int(shadow_dist)
                dem_tmp = dem
                out_shadow = morph_shade(dem_tmp, sh_z, sh_az, d_max, dem_size[0], dem_size[1], resolution)
                # out_shadow = adams_shadows(in_array=dem, az=azim, alt=elev, res=resolution, overlap=d_max)
                sky_illum_out += hillshade_tmp * out_shadow * weight
            else:
                sky_illum_out += hillshade_tmp * weight
        del hillshade_tmp

        if shadow_az and shadow_el:
            sky_illum_out = 0.8 * sky_illum_out + 0.2 * out_shadow

        dem[indx_no_values[0], indx_no_values[1]] = np.float64(np.NaN)

        return sky_illum_out


def local_dominance(dem, min_rad=10, max_rad=20, rad_inc=1, angular_res=15, observer_height=1.7, ve_factor=1):
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
    ve_factor : vertical exaggeration factor (must be greater than 0)

    Returns
    -------
    local_dom_out - 2D numpy array of local dominance
    """
    if ve_factor <= 0:
        raise Exception("rvt.vis.slope_aspect: ve_factor must be a positive number!")

    dem = dem.astype(np.float32)
    if ve_factor != 1:
        dem = dem * ve_factor
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

    return local_dom_out


# If we don't fix morp_shade we could use adams_shadows function (for that we would need to use numba)
# https://github.com/jacobdadams/general_scripts/blob/master/raster_chunk_processing.py
# import numba
#
#
# @numba.jit(nopython=True)
# def adams_shadows(in_array, az, alt, res, overlap, nodata=-1):
#     # Rows = i = y values, cols = j = x values
#     rows = in_array.shape[0]
#     cols = in_array.shape[1]
#     shadow_array = np.ones(in_array.shape)  # init to 1 (not shadowed), change to 0 if shadowed
#     max_elev = np.max(in_array)
#
#     az = 90. - az  # convert from 0 = north, cw to 0 = east, ccw
#
#     azrad = az * np.pi / 180.
#     altrad = alt * np.pi / 180.
#     delta_j = np.cos(azrad)
#     delta_i = -1. * np.sin(azrad)
#     tanaltrad = np.tan(altrad)
#
#     mult_size = 1
#     max_steps = 600
#
#     already_shadowed = 0
#
#     # precompute idx distances
#     distances = []
#     for d in range(1, max_steps):
#         distance = d * res
#         step_height = distance * tanaltrad
#         i_distance = delta_i * d
#         j_distance = delta_j * d
#         distances.append((step_height, i_distance, j_distance))
#
#     # Only compute shadows for the actual chunk area in a super_array
#     # We don't care about the overlap areas in the output array, they just get
#     # overwritten by the nodata value
#     if overlap > 0:
#         y_start = overlap - 1
#         y_end = rows - overlap
#         x_start = overlap - 1
#         x_end = cols - overlap
#     else:
#         y_start = 0
#         y_end = rows
#         x_start = 0
#         x_end = cols
#
#     for i in range(y_start, y_end):
#         for j in range(x_start, x_end):
#
#             point_elev = in_array[i, j]  # the point we want to determine if in shadow
#
#             for step in range(1, max_steps):  # start at a step of 1- a point cannot be shadowed by itself
#
#                 # No need to continue if it's already shadowed
#                 if shadow_array[i, j] == 0:
#                     already_shadowed += 1
#                     # print("shadow break")
#                     break
#
#                 critical_height = distances[step - 1][0] + point_elev
#
#                 # idx_i/j are indices of array corresponding to current position + y/x distances
#                 idx_i = int(round(i + distances[step - 1][1]))
#                 idx_j = int(round(j + distances[step - 1][2]))
#
#                 in_bounds = idx_i >= 0 and idx_i < rows and idx_j >= 0 and idx_j < cols
#                 in_height = critical_height < max_elev
#
#                 if in_bounds and in_height:
#                     next_elev = in_array[idx_i, idx_j]
#                     # Bail out if we hit a nodata area
#                     if next_elev == nodata or next_elev == np.nan:
#                         break
#
#                     if next_elev > point_elev and next_elev > critical_height:
#                         shadow_array[i, j] = 0
#
#                         # set all array indices in between our found shadowing index and the source index to shadowed
#                         for step2 in range(1, step):
#                             i2 = int(round(i + distances[step2 - 1][1]))
#                             j2 = int(round(j + distances[step2 - 1][2]))
#                             shadow_array[i2, j2] = 0
#
#                         break  # We're done with this point, move on to the next
#
#     return shadow_array
