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
import scipy.ndimage


# TODO: check speed sky_ilumination
# TODO: check IDL vectorisation remains:

def byte_scale(data,
               c_min=None,
               c_max=None,
               high=255,
               low=0,
               ):
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
    return np.cast[np.uint8](byte_data) + np.cast[np.uint8](low)


def slope_aspect(dem,
                 resolution_x,
                 resolution_y,
                 ve_factor=1,
                 output_units="radian",
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
    ve_factor : vertical exaggeration factor (must be greater than 0)
    output_units : percent, degree, radians

    Returns
    -------
    {"slope": slope_out, "aspect": aspect_out} : dictionaries with 2D numpy arrays
    """
    if dem.ndim != 2:
        raise Exception("rvt.vis.slope_aspect: dem has to be 2D np.array!")
    if ve_factor <= 0:
        raise Exception("rvt.vis.slope_aspect: ve_factor must be a positive number!")
    if resolution_x < 0 or resolution_y < 0:
        raise Exception("rvt.vis.slope_aspect: resolution must be a positive number!")

    dem = dem.astype(np.float32)
    dem = dem * ve_factor

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

    return {"slope": slope_out, "aspect": aspect_out}


def hillshade(dem,
              resolution_x,
              resolution_y,
              sun_azimuth=315,
              sun_elevation=35,
              slope=None,
              aspect=None,
              ve_factor=1
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
    ve_factor : vertical exaggeration factor (must be greater than 0)

    Returns
    -------
    hillshade_out : result numpy array
    """
    if dem.ndim != 2:
        raise Exception("rvt.vis.hillshade: dem has to be 2D np.array!")
    if ve_factor <= 0:
        raise Exception("rvt.vis.hillshade: ve_factor must be a positive number!")
    if sun_azimuth > 360 or sun_elevation > 90 or sun_azimuth < 0 or sun_elevation < 0:
        raise Exception("rvt.vis.analytical_hillshading: sun_azimuth must be [0-360] and sun_elevation [0-90]!")
    if resolution_x < 0 or resolution_y < 0:
        raise Exception("rvt.vis.analytical_hillshading: resolution must be a positive number!")

    dem = dem.astype(np.float32)
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

    # edges to -1
    hillshade_out[:, 0] = -1
    hillshade_out[0, :] = -1
    hillshade_out[:, -1] = -1
    hillshade_out[-1, :] = -1

    return hillshade_out


def multi_hillshade(dem,
                    resolution_x,
                    resolution_y,
                    nr_directions=16,
                    sun_elevation=35,
                    slope=None,
                    aspect=None,
                    ve_factor=1
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
    ve_factor : vertical exaggeration factor (must be greater than 0)

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
    if ve_factor <= 0:
        raise Exception("rvt.vis.multi_hillshade: ve_factor must be a positive number!")

    dem = dem.astype(np.float32)
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


def slrm(dem,
         radius_cell=20,
         ve_factor=1
         ):
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
    if dem.ndim != 2:
        raise Exception("rvt.vis.slrm: dem has to be 2D np.array!")

    if radius_cell < 10 or radius_cell > 50:
        raise Exception("rvt.vis.slrm: Radius for trend assessment needs to be in interval 10-50 pixels!")

    if ve_factor <= 0:
        raise Exception("rvt.vis.slrm: ve_factor must be a positive number!")

    dem = dem.astype(np.float32)
    dem = dem * ve_factor

    dem[dem < -1200] = np.float64(np.NaN)
    dem[dem > 2000] = np.float64(np.NaN)

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
    slrm_out = slrm_out / ((2*radius_cell+1)**2)  # calculate mean, 1=current pixel
    slrm_out = slrm_out[radius_cell:-radius_cell, radius_cell:-radius_cell]  # remove padding
    slrm_out = dem - slrm_out
    return slrm_out


def azimuth(xa,
            ya,
            xb,
            yb,
            ):
    """
    Determine the azimuth in the range of [0,2pi).

    Parameters
    ----------
    xa, ya, xb, yb (point coordinates)

    Returns
    -------
    a : outputs the azimuth in radians
    """
    # TODO, this is probably an obsolete function
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


def horizon_shift_vector(num_directions=16,
                         radius_pixels=10,
                         min_radius=1,
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

    Returns
    -------
    {"svf": svf_out, "asvf": asvf_out, "opns": opns_out} : dictionary
        svf_out, skyview factor : 2D array of skyview factor.
        asvf_out, anisotropic skyview factor : 2D array of anisotropic skyview factor.
        opns_out, openness : 2D array openness (elevation angle of horizon)
    """

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
    if compute_asvf:
        asvf_out = asvf_out[radius_max:-radius_max, radius_max:-radius_max] / np.sum(weight)
    if compute_opns:
        opns_out = np.rad2deg(0.5 * np.pi - (opns_out[radius_max:-radius_max, radius_max:-radius_max] / num_directions))

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
                    ve_factor=1
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
    if dem.ndim != 2:
        raise Exception("rvt.vis.sky_view_factor: dem has to be 2D np.array!")
    if ve_factor <= 0:
        raise Exception("rvt.vis.sky_view_factor: ve_factor must be a positive number!")
    if svf_noise != 0 and svf_noise != 1 and svf_noise != 2 and svf_noise != 3:
        raise Exception("rvt.vis.sky_view_factor: svf_noise must be one of the following values (0-don't remove, 1-low,"
                        " 2-med, 3-high)!")
    if asvf_level != 1 and asvf_level != 2:
        raise Exception("rvt.vis.sky_view_factor: asvf_leve must be one of the following values (1-low, 2-high)!")
    if not compute_svf and not compute_asvf and not compute_opns:
        raise Exception("rvt.vis.sky_view_factor: All computes are false!")

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
                                                 )

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
    ncols : number of columns
    nrows : number of rows
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
    if dem.ndim != 2:
        raise Exception("rvt.vis.sky_illumination: dem has to be 2D np.array!")
    if sampling_points != 250 and sampling_points != 500:
        raise Exception("rvt.vis.sky_illumination: sampling_points needs to be 250 or 500!")
    if sky_model != "overcast" and sky_model != "uniform":
        raise Exception("rvt.vis.sky_illumination: sky_model needs to be overcast or uniform!")
    if ve_factor <= 0:
        raise Exception("rvt.vis.sky_illumination: ve_factor must be a positive number!")

    dem = dem.astype(np.float32)
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


def local_dominance(dem,
                    min_rad=10,
                    max_rad=20,
                    rad_inc=1,
                    angular_res=15,
                    observer_height=1.7,
                    ve_factor=1):
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
    if dem.ndim != 2:
        raise Exception("rvt.vis.local_dominance: dem has to be 2D np.array!")
    if ve_factor <= 0:
        raise Exception("rvt.vis.local_dominance: ve_factor must be a positive number!")

    dem = dem.astype(np.float32)
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
