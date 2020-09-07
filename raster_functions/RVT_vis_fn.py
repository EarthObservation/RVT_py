"""
NAME:
    RVT_py visualization functions

PURPOSE:
    Contains all functions for visualization. All functions are rewritten from RVT (IDL)
    https://github.com/EarthObservation/RVT
"""

# python libraries
import numpy as np
import scipy.ndimage


def bytescale(data, cmin=None, cmax=None, high=255, low=0):
    """
    Byte scales an array (image).

    Byte scaling means converting the input image to uint8 dtype and scaling
    the range to ``(low, high)`` (default 0-255).

    Parameters
    ----------
    data : ndarray
        numpy image data array.
    cmin : scalar, optional
        Bias scaling of small values. Default is ``data.min()``.
    cmax : scalar, optional
        Bias scaling of large values. Default is ``data.max()``.
    high : scalar, optional
        Scale max value to `high`.  Default is 255.
    low : scalar, optional
        Scale min value to `low`.  Default is 0.

    Returns
    -------
    img_array : uint8 ndarray
        The byte-scaled array.

    """

    if high < low:
        raise ValueError("`high` should be larger than `low`.")

    if cmin is None:
        cmin = data.min()
    if cmax is None:
        cmax = data.max()

    cscale = cmax - cmin
    if cscale < 0:
        raise ValueError("`cmax` should be larger than `cmin`.")
    elif cscale == 0:
        cscale = 1

    if data.dtype == np.uint8:
        bytedata = (high + 1) * (data - cmin - 1) / (cmax - cmin)  # copied from IDL BYTSCL
        bytedata[bytedata > high] = high
        bytedata[bytedata < 0] = 0
        return np.cast[np.uint8](bytedata) + np.cast[np.uint8](low)

    # scale = float(high - low) / cscale  # old scipy fn
    # bytedata = (data * 1.0 - cmin) * scale + 0.4999  # old scipy fn

    bytedata = (high + 0.9999) * (data - cmin) / (cmax - cmin)  # copied from IDL BYTSCL
    bytedata[bytedata > high] = high
    bytedata[bytedata < 0] = 0
    return np.cast[np.uint8](bytedata) + np.cast[np.uint8](low)


"""
NAME:
    Slope, Aspect
    slope_aspect

DESCRIPTION:
    Procedure can return terrain slope and aspect in radian units (default) or in alternative units (if specified).
    Slope is defined as 0 for Hz plane and pi/2 for vertical plane.
    Aspect iz defined as geographic azimuth: clockwise increasing, 0 or 2pi for the North direction.
    Currently applied finite difference method.

INPUTS:
    input_DEM_arr       - input DEM 2D numpy array
    resolution_x      - DEM resolution in X direction
    resolution_y      - DEM resolution in Y direction
    ve_factor           - vertical exaggeration factor (must be greater than 0)
    output_units        - percent, degree, radians
    is_padding_applied  - is padding already applied on input array (needed for ArcGIS Pro which applies padding)

OUTPUTS:
    slope, aspect
        slope - result numpy array
        aspect - result numpy array
    
KEYWORDS:
    DEM_slope       - output terrain slope
    DEM_aspect      - output terrain aspect
    output_units    - percent, degree, radians
        percent         - output terrain slope in percent (0% for HZ surface, 100% for 45 degree tilted plane)
        degree          - output terrain slope and aspect in degrees
        radian          - output terrain slope and aspect in radians
    
DEPENDENCIES:
    /

AUTHOR:
    RVT:
        Klemen Zakšek (klemen.zaksek@zmaw.de)
    RVT_py:
        Žiga Maroh (ziga.maroh@icloud.com)

MODIFICATION HISTORY:
    RVT:
        Written by Klemen Zakšek, 2011.
    RVT_py:
        Written by Žiga Maroh, 2020.
"""


def slope_aspect(input_DEM_arr, resolution_x, resolution_y, ve_factor, is_padding_applied=False, output_units="radian"):
    if ve_factor <= 0:
        raise Exception("RVT slope_aspect: Vertical exagerration must be positive number!")
    dem = input_DEM_arr * ve_factor
    # add frame of 0 (additional row up bottom and column left right)
    if not is_padding_applied:
        dem = np.pad(dem, pad_width=1, mode="constant", constant_values=0)

    # Derivates in X and Y direction
    dzdx = ((np.roll(dem, 1, axis=1) - np.roll(dem, -1, axis=1)) / 2) / resolution_x
    dzdy = ((np.roll(dem, -1, axis=0) - np.roll(dem, 1, axis=0)) / 2) / resolution_y
    tan_slope = np.sqrt(dzdx ** 2 + dzdy ** 2)

    # Compute slope
    if output_units == "percent":
        slope = tan_slope * 100
    elif output_units == "degree":
        slope = np.rad2deg(np.arctan(tan_slope))
    elif output_units == "radian":
        slope = np.arctan(tan_slope)
    else:
        raise Exception("RVT calculate_slope: Wrong function input 'output_units'!")

    # Compute Aspect
    # Aspect identifies the downslope direction of the maximum rate of change in value from each cell to its neighbors:
    #     0
    # 270    90
    #    180
    dzdy[dzdy == 0] = 10e-9  # important for numeric stability - where dzdy is zero, make tangens to really high value

    aspect = np.arctan2(dzdx, dzdy)  # atan2 took care of the quadrants
    if output_units == "degree":
        aspect = np.rad2deg(aspect)

    # remove the frame (padding)
    slope = slope[1:-1, 1:-1]
    aspect = aspect[1:-1, 1:-1]

    # edges to -1
    slope[:, 0] = -1
    slope[0, :] = -1
    slope[:, -1] = -1
    slope[-1, :] = -1
    aspect[:, 0] = -1
    aspect[0, :] = -1
    aspect[:, -1] = -1
    aspect[-1, :] = -1

    return slope, aspect


"""
NAME:
    Analytical hillshading
    analytical_hillshading

DESCRIPTION:
    Compute hillshade.

INPUTS:
    input_DEM_arr       - input DEM 2D numpy array
    resolution_x        - DEM resolution in X direction
    resolution_y        - DEM resolution in Y direction
    sun_azimuth         - solar azimuth angle (clockwise from North) in degrees
    sun_elevation       - solar vertical angle (above the horizon) in degrees
    bytscl              - byte scale, if True scale values to 0-255 (u1, uint8)
    bytscl_min_max      - tuple(min, max) for bytscl (RVT: sc_hls_ev)
    is_padding_applied  - is padding already applied on input array (needed for ArcGIS Pro which applies padding)
    slope               - slope in radians if you don't input it, it is calculated
    aspect              - aspect in radians if you don't input it, it is calculated

OUTPUTS:
    hillshade - result numpy array

KEYWORDS:
    /

DEPENDENCIES:
    slope_aspect

AUTHOR:
    RVT:
        Klemen Zakšek (klemen.zaksek@zmaw.de)
    RVT_py:
        Žiga Maroh (ziga.maroh@icloud.com)

MODIFICATION HISTORY:
    RVT:
        1.0     Written by Klemen Zakšek, 2013.
        1.1     September 2014: Suppress_output and cosi keywords added to the procedure
    RVT_py:
        Written by Žiga Maroh, 2020.
"""


def analytical_hillshading(input_DEM_arr, resolution_x, resolution_y, sun_azimuth=315, sun_elevation=35, bytscl=True,
                           bytscl_min_max=(0, 1), is_padding_applied=False, slope=None, aspect=None):
    ve_factor = 1
    if sun_azimuth > 360 or sun_elevation > 90 or sun_azimuth < 0 or sun_elevation < 0:
        raise Exception("RVT analytical_hillshading: sun_azimuth must be [0-360] and sun_elevation [0-90]!")

    # Convert solar position (degrees) to radians
    sun_azimuth_rad = np.deg2rad(sun_azimuth)
    sun_elevation_rad = np.deg2rad(sun_elevation)

    # Convert to solar zenith angle
    sun_zenith_rad = np.pi / 2 - sun_elevation_rad

    # are slope and aspect already calculated and presented
    if slope is None or aspect is None:
        # calculates slope and aspect
        slope, aspect = slope_aspect(input_DEM_arr=input_DEM_arr, resolution_x=resolution_x, resolution_y=resolution_y,
                                     ve_factor=ve_factor, is_padding_applied=is_padding_applied, output_units="radian")

    # Compute solar incidence angle, hillshading
    hillshading = np.cos(sun_zenith_rad) * np.cos(slope) + np.sin(sun_zenith_rad) * np.sin(slope) * np.cos(
        aspect - sun_azimuth_rad)
    if bytscl:
        hillshading = bytescale(hillshading, cmin=bytscl_min_max[0], cmax=bytscl_min_max[1])

    return hillshading


"""
NAME:
    Multiple directions hillshading
    multiple_directions_hillshading

DESCRIPTION:
    Calculates hillshades from multiple directions.

INPUTS:
    input_DEM_arr           - input DEM 2D numpy array
    resolution_x            - DEM resolution in X direction
    resolution_y            - DEM resolution in Y direction
    nr_directions           - number of solar azimuth angles (clockwise from North)
    sun_elevation           - solar vertical angle (above the horizon) in degrees
    bytscl                  - byte scale, if True scale values to 0-255 (u1, uint8)
    bytscl_min_max          - tuple(min, max) for bytscl (RVT: sc_hls_ev)
    is_padding_applied      - is padding already applied on input array (needed for ArcGIS Pro which applies padding)
    slope               - slope in radians if you don't input it, it is calculated
    aspect              - aspect in radians if you don't input it, it is calculated

OUTPUTS:
    multi_hillshade - numpy array containing numpy_arrays of hillshades in different directions

KEYWORDS:
    /

DEPENDENCIES:
    slope_aspect
    analytical_hillshade

AUTHOR:
    RVT:
        Klemen Zakšek (klemen.zaksek@zmaw.de)
    RVT_py:
        Žiga Maroh (ziga.maroh@icloud.com)

MODIFICATION HISTORY:
    RVT:
        Written by Klemen Zakšek, 2013.
    RVT_py:
        Written by Žiga Maroh, 2020.
"""


def multiple_directions_hillshading(input_DEM_arr, resolution_x, resolution_y, nr_directions=16, sun_elevation=35,
                                    bytscl=True, bytscl_min_max=(0.00, 1.00),
                                    is_padding_applied=False, slope=None, aspect=None):
    if sun_elevation > 90 or sun_elevation < 0:
        raise Exception("RVT multiple_directions_hillshading: sun_elevation must be [0-90]!")

    ve_factor = 1

    # calculates slope and aspect if they are not added
    if slope is None or aspect is None:  # slope and aspect are the same, so we have to calculate it once
        slope, aspect = slope_aspect(input_DEM_arr=input_DEM_arr, resolution_x=resolution_x,
                                     resolution_y=resolution_y,
                                     ve_factor=ve_factor, is_padding_applied=is_padding_applied,
                                     output_units="radian")

    hillshades_arr_list = []  # list of all hillshades in diffrent directions
    for i_direction in range(nr_directions):
        sun_azimuth = (360 / nr_directions) * i_direction
        hillshade = analytical_hillshading(input_DEM_arr=input_DEM_arr, resolution_x=resolution_x,
                                           resolution_y=resolution_y, sun_elevation=sun_elevation,
                                           sun_azimuth=sun_azimuth, is_padding_applied=is_padding_applied, slope=slope,
                                           aspect=aspect, bytscl=bytscl, bytscl_min_max=bytscl_min_max)
        hillshades_arr_list.append(hillshade)

    multi_hillshade = np.asarray(hillshades_arr_list)

    return multi_hillshade


"""
NAME:
    Simple local relief model - SLRM
    SLRM

DESCRIPTION:
    Calculates simple local relief model.

INPUTS:
    input_DEM_arr           - input DEM 2D numpy array
    radius             - Radius for trend assessment [pixels]
    bytscl                  - byte scale, if True scale values to 0-255 (u1, uint8)
    bytscl_min_max          - tuple(min, max) for bytscl (RVT: sc_hls_ev)


OUTPUTS:
    slrm - SLRM 2D numpy array

KEYWORDS:
    /

DEPENDENCIES:
    /

AUTHOR:
    RVT:
        Klemen Zakšek (klemen.zaksek@zmaw.de)
    RVT_py:
        Žiga Maroh (ziga.maroh@icloud.com)

MODIFICATION HISTORY:
    RVT:
        Written by Klemen Zakšek, 2013.
    RVT_py:
        Written by Žiga Maroh, 2020.
"""


def SLRM(input_DEM_arr, radius_cell=20, bytscl=True, bytscl_min_max=(-2, 2)):
    if radius_cell < 10 or radius_cell > 50:
        import warnings
        raise Exception("RVT SLRM: Radius for trend assessment needs to be in interval 10-50 pixels!")

    dem = input_DEM_arr
    dem[dem < -1200] = np.nan
    dem[dem > 2000] = np.nan

    # mean filter
    slrm = dem - scipy.ndimage.uniform_filter(input_DEM_arr, mode='nearest', size=radius_cell * 2)

    if bytscl:
        slrm = bytescale(slrm, cmin=bytscl_min_max[0], cmax=bytscl_min_max[1])

    return slrm


"""
NAME:
    azimuth

DESCRIPTION:
     Determine the azimuth in the range of [0,2pi).

INPUTS:
    This function needs point coordinates xa,ya,xb,yb

OUTPUTS:
    a - outputs the azimuth in radians.

KEYWORDS:
    /

DEPENDENCIES:
    /

AUTHOR:
    RVT:
        Klemen Zakšek (klemen.zaksek@zmaw.de)
        Krištof Oštir (kristof.ostir@fgg.uni-lj.si)
    RVT_py:
        Žiga Maroh (ziga.maroh@icloud.com)

MODIFICATION HISTORY:
    RVT:
        Written by Klemen Zakšek, 2004.
        Implemented in IDL by Krištof Oštir, 2008.
    RVT_py:
        Written by Žiga Maroh, 2020.
"""


def azimuth(xa, ya, xb, yb):
    north = ya - yb
    east = xb - xa
    if north == 0:
        if east > 0:
            a = np.pi / 2
        else:
            if east < 0:
                a = 3 * np.pi / 2
            else:
                a = np.nan
    else:
        a0 = np.arctan(east / north)
        if north > 0 and east >= 0:
            a = a0
        elif north < 0:
            a = a0 + np.pi
        else:
            a = a0 + 2 * np.pi
    return a


"""
NAME:
    Sky-View determination movement matrix
    sky_view_det_move

DESCRIPTION:
    Determine the movement matrix for Sky-View computation.

INPUTS:
    num_directions - number of directions as input
    radius_cell - radius to consider in cells (not in meters)
    ncol - number of columns of the input DEM

OUTPUTS:
    move - 2D numpy array movement matrix.

KEYWORDS:
    /

DEPENDENCIES:
    azimuth

AUTHOR:
    RVT:
        Klemen Zakšek (klemen.zaksek@zmaw.de)
        Krištof Oštir (kristof.ostir@fgg.uni-lj.si)
    RVT_py:
        Žiga Maroh (ziga.maroh@icloud.com)

MODIFICATION HISTORY:
    RVT:
        Written by Klemen Zakšek, 2005.
        Implemented in IDL by Krištof Oštir, 2008.
        Optimized by Klemen Zakšek, 2013.
    RVT_py:
        Written by Žiga Maroh, 2020.
"""


def sky_view_det_move(num_directions, radius_cell, ncol):
    # Parameters
    look_angle = 2 * np.pi / num_directions

    # Matrix initialization
    move = np.zeros((3, radius_cell + 1, num_directions))

    # For each direction
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

    # Reformat the radius:
    # First row tells you, how many valid cells are available, bellow is actual radius in cells
    move[2, 1:radius_cell + 1, :] = move[2, 0:radius_cell, :]
    for i_direction in range(num_directions):
        tmp = move[2, 1:radius_cell + 1, i_direction]
        if tmp[tmp > radius_cell].size > 0:
            move[2, 0, i_direction] = np.min(tmp[tmp > radius_cell])
        else:
            move[2, 0, i_direction] = radius_cell

    # Convert 2D index into 1D index
    move_t = np.zeros((2, radius_cell + 1, num_directions))
    move_t[0, :, :] = move[1, :, :] * ncol + move[0, :, :]
    move_t[1, :, :] = move[2, :, :]
    move = move_t

    return move


"""
NAME:
    Sky-View factor compute
    sky_view_factor_compute

DESCRIPTION:
    Compute the Sky-View Factor.

INPUTS:
    height_arr              - elevation (DEM) as 2D numpy array
    i_valid                 - index of valid pixels to be processed
    radius_max              - maximal search radius in pixels/cells (not in meters)
    radius_min              - minimal search radius in pixels/cells (not in meters), for noise reduction
    num_directions          - number of directions as input
    a_main_direction        - main direction of anisotropy
    a_poly_level            - level of polynomial that determines the anisotropy
    a_min_wight             - weight to consider anisotropy (0 - isotropic, 1 - no illumination from the direction 
                                                             opposite the main direction)
    compute_svf             - if true it computes svf, if false output svf = False
    compute_asvf            - if true it computes asvf, if false output asvf = False
    compute_opns            - if true it computes opns, if false output opns = False
    
OUTPUTS:
    svf, asvf, opns
        svf, skyview_factor               - 2D numpy array of skyview factor.
        asvf, anisotropic_skyview_factor  - 2D numpy array of anisotropic skyview factor.
        opns, openness                    - openness (elevation angle of horizon)

KEYWORDS:
    /

DEPENDENCIES:
    azimuth
    sky_view_det_move

AUTHOR:
    RVT:
        Klemen Zakšek (klemen.zaksek@zmaw.de)
        Krištof Oštir (kristof.ostir@fgg.uni-lj.si)
    RVT_py:
        Žiga Maroh (ziga.maroh@icloud.com)

MODIFICATION HISTORY:
    RVT:
        Written by Klemen Zakšek, 2005.
        Implemented in IDL by Krištof Oštir, 2008.
        Optimized by Klemen Zakšek, 2009.
        Rewritten (cleaner code + option of anisometric SCF and openess) by Klemen Zakšek, 2013.
    RVT_py:
        Written by Žiga Maroh, 2020.
"""


def sky_view_factor_compute(height_arr, i_valid, radius_max, radius_min, num_directions,
                            a_main_direction, a_poly_level, a_min_weight, compute_svf=True, compute_asvf=True,
                            compute_opns=True):
    height = height_arr

    svf = False  # if compute_svf is False then function doesn't compute svf and returns False
    asvf = False  # if compute_asvf is False then function doesn't compute asvf and returns False
    opns = False  # if compute_opns is False then function doesn't compute opns and returns False

    # directional step
    dir_step = 2 * np.pi / num_directions

    # vector of movement
    ncol = height.shape[1]  # number of columns
    count_height = i_valid[0].size  # number of all elements
    move = sky_view_det_move(num_directions=num_directions, radius_cell=radius_max, ncol=ncol)
    # init the outputs
    if compute_svf:
        svf = np.zeros(count_height)
    if compute_opns:
        opns = np.zeros(count_height)
    if compute_asvf:
        asvf = np.zeros(count_height)
        w_m = a_min_weight  # compute weights
        w_a = np.deg2rad(a_main_direction)
        weight = np.arange(num_directions) * dir_step
        weight = (1 - w_m) * (np.cos((weight - w_a) / 2)) ** a_poly_level + w_m

    # look into each direction...
    for i_dir in range(num_directions):
        # reset maximum at each iteration - at each new direction
        max_slope = np.zeros(count_height) - 1000

        # ... and to the search radius - this depends on the direction - radius is written in the first row
        for i_rad in range(1, int(move[1, 0, i_dir])):
            # ignore radios if smaller than minimal defined radius
            if radius_min > move[1, i_rad, i_dir]:
                continue
            # search for max (sky)
            h_flt = height.flatten()
            np.seterr(divide='ignore', invalid='ignore')  # ignore warnings for dividing with zero
            m_slp = (h_flt[i_valid[0] + int(move[0, int(i_rad - 1), i_dir])] - h_flt[i_valid[0]]) / move[
                1, i_rad, i_dir]
            np.seterr(divide='warn', invalid='warn')  # reset warnings
            max_slope = (max_slope < m_slp).choose(max_slope, m_slp)

        max_slope = np.arctan(max_slope)

        if compute_opns:
            opns = opns + max_slope
        if compute_svf:
            svf = svf + (1 - np.sin((max_slope < 0).choose(max_slope, 0)))
        if compute_asvf:
            asvf = asvf + (1 - np.sin((max_slope < 0).choose(max_slope, 0))) * weight[i_dir]

    # Normalize to the number of directions / weights
    if compute_svf:
        svf = svf / num_directions
    if compute_asvf:
        asvf = asvf / np.sum(weight)
    if compute_opns:
        opns = np.pi / 2 - (opns / num_directions)

    return svf, asvf, opns


"""
NAME:
    Sky-View factor
    sky_view_factor

DESCRIPTION:
    Prepare the data and compute the Sky-View Factor.

INPUTS:
    input_DEM_arr           - input DEM 2D numpy array (Ve Exaggeration and pixel size already considered)
    compute_svf             - compute SVF (True) or not (False)
    compute_opns            - compute OPENNESS (True) or not (False)
    resolution              - pixel resolution
    svf_n_dir            - number of directions
    svf_r_max            - maximal search radius in pixels
    svf_noise            - the level of noise remove (0-3)
    compute_asvf            - compute anisotropic SVF (True) or not (False)
    asvf_level           - level of anisotropy, 1-low, 2-high,
    a_min_weight            - weight to consider anisotropy (0 - isotropic, 1 - no illumination from the direction opposite the main direction)
    bytscl                  - byte scale, if True scale values to 0-255 (u1, uint8)
    bytscl_min_max_svf      - tuple(min, max) for bytscl (RVT: sc_svf_ev)
    bytscl_min_max_opns     - tuple(min, max) for bytscl (RVT: sc_opns_ev)
    bytscl_min_max_asvf     - tuple(min, max) for bytscl (RVT: sc_asvf_ev)
        CONSTANTS:
            sc_asvf_min             - level of polynomial that determines the anisotropy, selected with in_asvf_level
            sc_asvf_pol             - level of polynomial that determines the anisotropy, selected with in_asvf_level
            sc_svf_r_min            - the portion (percent) of the maximal search radius to ignore in horizon estimation; for each noise level, selected with in_svf_noise
    
    
OUTPUTS:
    svf, asvf, opns
        svf, skyview_factor               - 2D numpy array of skyview factor.
        asvf, anisotropic_skyview_factor  - 2D numpy array of anisotropic skyview factor.
        opns, openness                    - openness (elevation angle of horizon)

KEYWORDS:
    /

DEPENDENCIES:
    azimuth
    sky_view_det_move
    sky_view_compute

AUTHOR:
    RVT:
        Klemen Zakšek (klemen.zaksek@zmaw.de)
        Krištof Oštir (kristof.ostir@fgg.uni-lj.si)
    RVT_py:
        Žiga Maroh (ziga.maroh@icloud.com)

MODIFICATION HISTORY:
    RVT:
        Written by Klemen Zakšek, 2013.
    RVT_py:
        Written by Žiga Maroh, 2020.
"""


def sky_view_factor(input_DEM_arr, resolution, compute_svf=True, compute_opns=True, compute_asvf=True,
                    svf_n_dir=16, svf_r_max=10, svf_noise=0, asvf_dir=315, asvf_level=1,
                    bytscl=True, bytscl_min_max_svf=(0.6375, 1.00),
                    bytscl_min_max_opns=(60, 95.), bytscl_min_max_asvf=(0.6375, 1.00)):
    # CONSTANTS
    # level of polynomial that determines the anisotropy, selected with in_asvf_level (1 - low, 2 - high)
    sc_asvf_pol = [4, 8]
    sc_asvf_min = [0.4, 0.1]
    # the portion (percent) of the maximal search radius to ignore in horizon estimation; for each noise level,
    # selected with in_svf_noise (0-3)
    sc_svf_r_min = [0., 10., 20., 40.]

    # pixel size
    dem = input_DEM_arr / resolution
    # increase edge - visualization foes to edge, so mirror them, leave blank corners
    ncol = dem.shape[1]
    nlin = dem.shape[0]
    tmp = np.zeros((nlin + 2 * svf_r_max, ncol + 2 * svf_r_max))

    # prepare indx for output
    tmp[svf_r_max:nlin + svf_r_max, svf_r_max:ncol + svf_r_max] = 1

    indx_all = np.where(tmp.flatten() == 1)

    # fill center
    tmp[svf_r_max:nlin + svf_r_max, svf_r_max:ncol + svf_r_max] = dem

    # final dem with increased edges
    dem = tmp
    del tmp

    # minimal search radious depends on the noise level
    svf_r_min = svf_r_max * sc_svf_r_min[svf_noise] * 0.01

    # set anisotropy parameters
    poly_level = sc_asvf_pol[asvf_level - 1]
    min_weight = sc_asvf_min[asvf_level - 1]

    svf, asvf, opns = sky_view_factor_compute(height_arr=dem, i_valid=indx_all, radius_max=svf_r_max,
                                              radius_min=svf_r_min, num_directions=svf_n_dir,
                                              a_main_direction=asvf_dir, a_poly_level=poly_level,
                                              a_min_weight=min_weight, compute_svf=compute_svf,
                                              compute_asvf=compute_asvf, compute_opns=compute_opns)
    # reshape 1D to 2D
    if compute_svf:
        svf = svf.reshape((nlin, ncol))
    if compute_asvf:
        asvf = asvf.reshape((nlin, ncol))
    if compute_opns:
        opns = opns.reshape((nlin, ncol))

    if compute_svf and bytscl:
        svf = bytescale(svf, cmin=bytscl_min_max_svf[0], cmax=bytscl_min_max_svf[1])
    if compute_asvf and bytscl:
        asvf = bytescale(asvf, cmin=bytscl_min_max_asvf[0], cmax=bytscl_min_max_asvf[1])
    if compute_opns and bytscl:
        opns = bytescale(opns, cmin=bytscl_min_max_opns[0], cmax=bytscl_min_max_opns[1])

    return svf, asvf, opns


"""
NAME:
    morph_shade_move

DESCRIPTION:
    Calculates azimuth for morph_shade

INPUTS:
    d_max  - maximum search distance in pixels
    angle

OUTPUTS:
    move

KEYWORDS:
    /

DEPENDENCIES:
    /

AUTHOR:
    RVT:
        Klemen Čotar
    RVT_py:
        Žiga Maroh (ziga.maroh@icloud.com)

MODIFICATION HISTORY:
    RVT:
        Initial version written by Klemen Čotar, 2014.
    RVT_py:
        Written by Žiga Maroh, 2020.
"""


def morph_shade_move(d_max, angle):
    # init
    move = np.zeros((d_max + 1, 3))
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


"""
NAME:
    morph_shade

DESCRIPTION:
    Compute topographic corrections for sky illumination.

INPUTS:
    height - elevation 2D np array
    sol_z  - solar zenith angle in rad (0 for vertical and pi/2 for horizontal surface)
    sol_a  - solar azimuth angle
    d_max  - maximum search distance in pixels
    resolution  - pixel size

OUTPUTS:
    This procedure determines those cells that are in its own (hillshade) or thrown (cast shade) shadow

KEYWORDS:
    /

DEPENDENCIES:
    morph_shade_move

AUTHOR:
    RVT:
        Klemen Čotar
    RVT_py:
        Žiga Maroh (ziga.maroh@icloud.com)

MODIFICATION HISTORY:
    RVT:
        Initial version written by Klemen Čotar, 2014.
    RVT_py:
        Written by Žiga Maroh, 2020.
"""


def morph_shade(height, sol_z, sol_a, d_max, nrows, ncols, resolution):
    # initialize the results
    # print(d_max)
    mask = np.zeros((nrows + 2 * d_max, ncols + 2 * d_max))
    mask[d_max:(nrows + d_max), d_max:(ncols + d_max)] = 1
    i_valid = np.where(mask.flatten() == 1)
    tmp = np.zeros((nrows + 2 * d_max, ncols + 2 * d_max))
    tmp[d_max:(nrows + d_max), d_max:(ncols + d_max)] = height
    height = tmp
    del tmp

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
        zeros = np.zeros(nr_zeros)
        sel = sel[sel <= height_flt.size]
        m_slp = ((height_flt[sel] - height_flt[i_valid[0][:sel.size]]) / move1dd[i_rad])
        print(m_slp)
        m_slp = np.append(m_slp, zeros)  # add zeros, for non existing indexes
        max_slope = (max_slope < m_slp).choose(max_slope, m_slp)

    # update mask
    max_slope = np.arctan(max_slope / resolution)
    indx_mask = np.where(max_slope > (np.pi / 2 - sol_z))

    # print(mask[i_valid[0][indx_mask[0]]])
    if indx_mask[0].size > 0:
        mask_size = mask.shape
        mask_flt = mask.flatten()
        mask_flt[i_valid[0][indx_mask[0]]] = 0
        mask = mask_flt.reshape((mask_size[0], mask_size[1]))
    mask = mask[d_max:(nrows + d_max), d_max:(ncols + d_max)]
    return mask


"""
NAME:
    Sky illumination
    sky_illumination

DESCRIPTION:
    Computes Sky ilumination.

INPUTS:
    input_DEM_arr   - numpy 2D array of elevation (DEM)
    resolution      - dem pixel size
    sky_model       - sky model [overcast, uniform]
    sampling_points - number of sampling points
    shadow_dist     - max shadow modeling distance [pixels]
    shadow_az       - shadow azimuth
    shadow_el       - shadow elevation
    shadow_only     - bool compute shadow only

OUTPUTS:
    skyilumination

KEYWORDS:
    /

DEPENDENCIES:
    morph_shade_move
    morph_shade

AUTHOR:
    RVT:
        Klemen Čotar
    RVT_py:
        Žiga Maroh (ziga.maroh@icloud.com)

MODIFICATION HISTORY:
    RVT:
        Initial version written by Klemen Čotar, 2014.
    RVT_py:
        Written by Žiga Maroh, 2020.
"""


def sky_illumination(input_DEM_arr, resolution, sky_model="overcast", sampling_points=250, shadow_dist=100,
                     shadow_az=315, shadow_el=35, shadow_only=False):
    if sampling_points != 250 and sampling_points != 500:
        raise Exception("RVT sky_illumination: sampling_points needs to be 250 or 500!")
    if sky_model != "overcast" and sky_model != "uniform":
        raise Exception("RVT sky_illumination: sky_model needs to be overcast or uniform!")

    dem = input_DEM_arr
    indx_novalues = np.where(dem < 0)
    dem[indx_novalues[0], indx_novalues[1]] = np.nan

    # determine the max search distance in pixels
    h_min = np.nanmin(dem)
    h_max = np.nanmax(dem)
    dh = h_max - h_min
    dem_size = dem.shape
    dem[indx_novalues[0], indx_novalues[1]] = 0

    if shadow_az and shadow_el:
        sh_z = np.pi / 2 - np.deg2rad(shadow_el)
        sh_az = np.deg2rad(shadow_az)
        d_max = 1
        d_max = round(d_max * dh * np.tan(sh_z) / resolution)
        dem_tmp = dem
        out_shadow = morph_shade(dem_tmp, sh_z, sh_az, d_max, dem_size[0], dem_size[1], resolution)

    if shadow_only:
        dem[indx_novalues[0], indx_novalues[1]] = np.nan
        return out_shadow

    else:
        sc_skyilu_ev = [0.25, 0.]  # percent
        scale_lower = sc_skyilu_ev[0]
        scale_upper = sc_skyilu_ev[1]

        dat_hillset = open(r'settings\{}_{}sp.txt'.format(sky_model, sampling_points), 'r')
        skyillumination = np.zeros((dem_size[0], dem_size[1]))
        slope, aspect = slope_aspect(input_DEM_arr=input_DEM_arr, resolution_x=resolution,
                                     resolution_y=resolution,
                                     ve_factor=1, is_padding_applied=False,
                                     output_units="radian")
        for line in dat_hillset:
            if line.strip() == "":  # empty line
                continue
            print(line)
            d_max = 1
            line = line.rstrip().split(",")
            azim = int(line[0])
            elev = int(line[1])
            weight = float(line[2])
            hillshade = analytical_hillshading(input_DEM_arr=dem, resolution_x=resolution, resolution_y=resolution,
                                               sun_azimuth=azim, sun_elevation=elev, slope=slope, aspect=aspect,
                                               is_padding_applied=False, bytscl=False)
            sh_z = np.pi / 2 - np.deg2rad(elev)
            sh_az = np.deg2rad(azim)
            d_max = round(d_max * dh * np.tan(sh_z) / resolution)

            if d_max > 1:
                if shadow_dist != 'unlimited':
                    if d_max < int(shadow_dist):
                        d_max = int(shadow_dist)
                dem_tmp = dem
                out_shadow = morph_shade(dem_tmp, sh_z, sh_az, d_max, dem_size[0], dem_size[1], resolution)
                skyillumination += hillshade * out_shadow * weight
            else:
                skyillumination += hillshade * weight
        del hillshade

        if shadow_az and shadow_el:
            skyillumination = 0.8 * skyillumination + 0.2 * out_shadow

        dem[indx_novalues[0], indx_novalues[1]] = np.nan

        return skyillumination


"""
NAME:
    Local dominance
    local_dominance

DESCRIPTION:
    Compute Local Dominance dem visualization.
    
    Adapted from original version that is part of the Lida Visualisation Toolbox LiVT developed by Ralf Hesse.

INPUTS:
    input_DEM_arr   - input DEM 2D numpy array
    min_rad         - minimum radial distance (in pixels) at which the algorithm starts with visualization computation
    max_rad         - maximum radial distance (in pixels) at which the algorithm ends with visualization computation
    rad_inc         - radial distance steps in pixels
    angular_res     - angular step for determination of number of angular directions
    observer_height - height at which we observe the terrain

OUTPUTS:
    localdominance - 2D numpy array of localdominance

KEYWORDS:
    /

DEPENDENCIES:
    /

AUTHOR:
    RVT:
        Klemen Čotar

MODIFICATION HISTORY:
    RVT:
        Initial version written by Klemen Čotar, 2016.
    RVT_py:
        Written by Žiga Maroh, 2020.
"""


def local_dominance(input_DEM_arr, min_rad=10, max_rad=20, rad_inc=1, angular_res=15, observer_height=1.7):
    dem = input_DEM_arr

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

    localdominance = dem * 0
    for i_s in range(n_shifts):
        dem_moved = np.roll(dem, int(round(y_t[i_s])), axis=0)
        dem_moved = np.roll(dem_moved, int(round(x_t[i_s])), axis=1)
        idx_lower = np.where((dem + observer_height) > dem_moved)
        if idx_lower[0].size > 0:
            localdominance[idx_lower[0], idx_lower[1]] += (dem[idx_lower[0], idx_lower[1]] + observer_height -
                                                           dem_moved[idx_lower[0], idx_lower[1]]) / \
                                                          distances[i_s] * dist_factr[i_s]
    localdominance /= norma

    return localdominance
