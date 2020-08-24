"""
NAME:
    RVT Analytical hillshading, ESRI ArcGis Pro

DESCRIPTION:
    https://github.com/Esri/raster-functions/wiki/PythonRasterFunction#anatomy-of-a-python-raster-function
    Python raster function for ESRI ArcGis which computes hillshade.

INPUTS:
    input_DEM_arr   - input DEM numpy array
    resolution      - DEM resolution
    sun_azimuth     - solar azimuth angle (clockwise from North) in degrees
    sun_elevation   - solar vertical angle (above the horizon) in degrees

OUTPUTS:
    hillshade - result numpy array

KEYWORDS:
    /

DEPENDENCIES:
    RVT_vis_fn.slope_aspect
    RVT_vis_fn.analytical_hillshading

AUTHOR:
    Žiga Maroh (ziga.maroh@icloud.com)

MODIFICATION HISTORY:
    1.0 Written by Žiga Maroh, 2020.
"""

import numpy as np


class RVT_Hillshade():
    def __init__(self):
        self.name = 'RVT Hillshade Function'
        self.description = 'Computes hillshade.'
        self.zf = 'f'
        self.azimuth = 'f'
        self.elevation = 'f'

    def getParameterInfo(self):
        return [
            {
                'name': 'raster',
                'dataType': 'raster',
                'value': None,
                'required': True,
                'displayName': "Input Raster",
                'description': "The primary input raster where pixel values represent elevation.",
            },
            {
                'name': 'zf',
                'dataType': 'numeric',
                'value': 1.,
                'required': False,
                'displayName': "Z Factor",
                'description': (
                    "The multiplicative factor that converts elevation values to the units of the horizontal (xy-) coordinate system. "
                    "Or use larger values to add vertical exaggeration."),
            },
            {
                'name': 'sun_azimuth',
                'dataType': 'numeric',
                'value': 315.,
                'required': False,
                'displayName': "Sun azimuth",
                'description': (
                    "Solar azimuth angle (clockwise from North) in degrees"),
            },
            {
                'name': 'sun_elevation',
                'dataType': 'numeric',
                'value': 35.,
                'required': False,
                'displayName': "Sun elevation",
                'description': (
                    "Solar vertical angle (above the horizon) in degrees"),
            },
        ]

    def getConfiguration(self, **scalars):
        return {
            'inheritProperties': 4 | 8,  # inherit everything but the pixel type (1) and NoData (2)
            'invalidateProperties': 2 | 4 | 8,
            # invalidate these aspects because we are modifying pixel values and updating key properties.
            'padding': 1,  # one extra on each each of the input pixel block
            'inputMask': False  # Don't need input raster mask in .updatePixels(). Simply use the inherited NoData.
        }

    def updateRasterInfo(self, **kwargs):
        kwargs['output_info']['bandCount'] = 1
        kwargs['output_info']['pixelType'] = 'f4'
        kwargs['output_info']['statistics'] = ({'minimum': 0., 'maximum': 255.},)
        kwargs['output_info']['histogram'] = ()

        r = kwargs['raster_info']
        if r['bandCount'] > 1:
            raise Exception("Input raster has more than one band. Only single-band raster datasets are supported")

        self.azimuth = kwargs.get('sun_azimuth')
        self.zf = kwargs.get('zf')
        self.elevation = kwargs.get('sun_elevation')

        return kwargs

    def updatePixels(self, tlc, shape, props, **pixelBlocks):
        input_arr = np.array(pixelBlocks['raster_pixels'], dtype='f4', copy=False)[0]
        hillshade = self.analytical_hillshading(input_DEM_arr=input_arr, resolution=pixelBlocks['raster_cell_size'],
                                                sun_azimuth=self.azimuth, sun_elevation=self.elevation,
                                                ve_factor=self.zf)

        pixelBlocks['output_pixels'] = hillshade[1:-1, 1:-1].astype(props['pixelType'], copy=False)  # remove padding
        return pixelBlocks

    def updateKeyMetadata(self, names, bandIndex, **keyMetadata):
        if bandIndex == -1:  # dataset-level properties
            keyMetadata['datatype'] = 'Processed'  # outgoing dataset is now 'Processed'
        elif bandIndex == 0:  # properties for the first band
            keyMetadata['wavelengthmin'] = None  # reset inapplicable band-specific key metadata
            keyMetadata['wavelengthmax'] = None
            keyMetadata['bandname'] = 'RVT_Hillshade'
        return keyMetadata

    def analytical_hillshading(self, input_DEM_arr, resolution, sun_azimuth=315, sun_elevation=35, ve_factor=1):
        # Convert solar position (degrees) to radians
        sun_azimuth_rad = np.deg2rad(sun_azimuth)
        sun_elevation_rad = np.deg2rad(sun_elevation)

        # Convert to solar zenith angle
        sun_zenith_rad = np.pi / 2 - sun_elevation_rad

        # Compute solar incidence angle, hillshading
        slope, aspect = self.slope_aspect(input_DEM_arr, resolution, ve_factor)
        hillshading = np.cos(sun_zenith_rad) * np.cos(slope) + np.sin(sun_zenith_rad) * np.sin(slope) * \
                      np.cos(aspect - sun_azimuth_rad)

        return hillshading

    def slope_aspect(self, input_DEM_arr, resolution, ve_factor):
        dem = input_DEM_arr * ve_factor

        # Derivates in X and Y direction
        dzdx = ((np.roll(dem, 1, axis=1) - np.roll(dem, -1, axis=1)) / 2) / resolution
        dzdy = ((np.roll(dem, -1, axis=0) - np.roll(dem, 1, axis=0)) / 2) / resolution
        tan_slope = np.sqrt(dzdx ** 2 + dzdy ** 2)

        slope = np.arctan(tan_slope)

        # Compute Aspect
        # Aspect identifies the downslope direction of the maximum rate of change in value from each cell to its
        # neighbors:
        #     0
        # 270    90
        #    180

        # important for numeric stability - where dzdy zero is, make tangens to really high value
        dzdy[dzdy == 0] = 10e-9

        aspect = np.arctan2(dzdx, dzdy)  # atan2 took care of the quadrants

        return slope, aspect
