"""
NAME:
    RVT_py visualization functions

PURPOSE:
    Contains all functions for visualization.
"""

# python libraries
import numpy as np

"""
NAME:
    Analytical hillshading

DESCRIPTION:
    Compute hillshade.

INPUTS:
    input_DEM_arr   - input DEM numpy array
    resolution_x      - DEM resolution in X direction
    resolution_y      - DEM resolution in Y direction
    sun_azimuth     - solar azimuth angle (clockwise from North) in degrees
    sun_elevation   - solar vertical angle (above the horizon) in degrees
    is_padding_applied  - is padding already applied on input array (needed for ArcGIS Pro which applies padding)

OUTPUTS:
    hillshade - result numpy array
    
KEYWORDS:
    /

DEPENDENCIES:
    slope_aspect

AUTHOR:
    Žiga Maroh (ziga.maroh@icloud.com)

MODIFICATION HISTORY:
    1.0 Written by Žiga Maroh, 2020.
"""


def analytical_hillshading(input_DEM_arr, resolution_x, resolution_y, sun_azimuth=315, sun_elevation=35,
                           is_padding_applied=False):
    ve_factor = 1
    if sun_azimuth > 360 or sun_elevation > 90 or sun_azimuth < 0 or sun_elevation < 0:
        raise Exception("RVT analytical_hillshading: sun_azimuth must be [0-360] and sun_elevation [0-90]")

    # Convert solar position (degrees) to radians
    sun_azimuth_rad = np.deg2rad(sun_azimuth)
    sun_elevation_rad = np.deg2rad(sun_elevation)

    # Convert to solar zenith angle
    sun_zenith_rad = np.pi / 2 - sun_elevation_rad

    # Compute solar incidence angle, hillshading
    slope, aspect = slope_aspect(input_DEM_arr=input_DEM_arr, resolution_x=resolution_x, resolution_y=resolution_y,
                                 ve_factor=ve_factor, is_padding_applied=is_padding_applied, output_units="radian")

    hillshading = np.cos(sun_zenith_rad) * np.cos(slope) + np.sin(sun_zenith_rad) * np.sin(slope) * np.cos(aspect - sun_azimuth_rad)

    return hillshading


"""
NAME:
    Slope, Aspect

DESCRIPTION:
    Procedure can return terrain slope and aspect in radian units (default) or in alternative units (if specified).
    Slope is defined as 0 for Hz plane and pi/2 for vertical plane.
    Aspect iz defined as geographic azimuth: clockwise increasing, 0 or 2pi for the North direction.
    Currently applied finite difference method.

INPUTS:
    input_DEM_arr       - input DEM numpy array
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
    Žiga Maroh (ziga.maroh@icloud.com)

MODIFICATION HISTORY:
    1.0 Written by Žiga Maroh, 2020.
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

    # change no data to -1
    slope[slope == np.NaN] = -1
    aspect[aspect == np.NaN] = -1

    return slope, aspect
