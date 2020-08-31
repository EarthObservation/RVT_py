"""
NAME:
    RVT_py visualization functions

PURPOSE:
    Contains all functions for visualization.
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
        Klemen Zaksek (klemen.zaksek@zmaw.de)
    RVT_py:
        Žiga Maroh (ziga.maroh@icloud.com)

MODIFICATION HISTORY:
    RVT:
        Written by Klemen Zaksek, 2011.
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
        Klemen Zaksek (klemen.zaksek@zmaw.de)
    RVT_py:
        Žiga Maroh (ziga.maroh@icloud.com)

MODIFICATION HISTORY:
    RVT:
        1.0     Written by Klemen Zaksek, 2013.
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
        Klemen Zaksek (klemen.zaksek@zmaw.de)
    RVT_py:
        Žiga Maroh (ziga.maroh@icloud.com)

MODIFICATION HISTORY:
    RVT:
        Written by Klemen Zaksek, 2013.
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
        Klemen Zaksek (klemen.zaksek@zmaw.de)
    RVT_py:
        Žiga Maroh (ziga.maroh@icloud.com)

MODIFICATION HISTORY:
    RVT:
        Written by Klemen Zaksek, 2013.
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
        Klemen Zaksek (klemen.zaksek@zmaw.de)
        Kristof Ostir (kristof.ostir@fgg.uni-lj.si)
    RVT_py:
        Žiga Maroh (ziga.maroh@icloud.com)

MODIFICATION HISTORY:
    RVT:
        Written by Klemen Zaksek, 2004.
        Implemented in IDL by Kristof Ostir, 2008.
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
        Klemen Zaksek (klemen.zaksek@zmaw.de)
        Kristof Ostir (kristof.ostir@fgg.uni-lj.si)
    RVT_py:
        Žiga Maroh (ziga.maroh@icloud.com)

MODIFICATION HISTORY:
    RVT:
        Written by Klemen Zaksek, 2005.
        Implemented in IDL by Kristof Ostir, 2008.
        Optimized by Klemen Zaksek, 2013.
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

