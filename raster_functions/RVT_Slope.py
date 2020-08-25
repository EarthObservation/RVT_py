"""
NAME:
    RVT Slope, ESRI ArcGis Pro python function

DESCRIPTION:
    https://github.com/Esri/raster-functions/wiki/PythonRasterFunction#anatomy-of-a-python-raster-function
    Python raster function for ESRI ArcGis Pro which computes Slope.
    Slope is defined as 0 for Hz plane and pi/2 for vertical plane.

INPUTS:
    input_DEM_arr       - input DEM numpy array
    resolution          - DEM resolution
    ve_factor           - vertical exaggeration factor (must be greater than 0)
    output_units        - percent, degree, radians
    is_padding_applied  - is padding already applied on input array (needed for ArcGIS Pro which applies padding)


OUTPUTS:
    slope

KEYWORDS:
    DEM_slope       - output terrain slope
    DEM_aspect      - output terrain aspect
    output_units    - percent, degree, radians
        percent         - output terrain slope in percent (0% for HZ surface, 100% for 45 degree tilted plane)
        degree          - output terrain slope and aspect in degrees
        radian          - output terrain slope and aspect in radians

DEPENDENCIES:
    RVT_vis_fn.slope_aspect

AUTHOR:
    Žiga Maroh (ziga.maroh@icloud.com)

MODIFICATION HISTORY:
    1.0 Written by Žiga Maroh, 2020.
"""


import numpy as np
import RVT_vis_fn


class RVT_Slope():

    def __init(self):
        self.name = "RVT Slope"
        self.description = "Calculates Slope in selected unit."
        self.prepare()

    def getParameterInfo(self):
        return [
            {
                'name': 'raster',
                'dataType': 'raster',
                'value': None,
                'required': True,
                'displayName': "Input Raster",
                'description': "Input raster for which to create the slope map."
            },
            {
                'name': 've_factor',
                'dataType': 'numeric',
                'value': 1.,
                'required': False,
                'displayName': "Ve-factor",
                'description': ("Vertical exaggeration factor (must be greater than 0).")
            },
            {
                'name': 'output_unit',
                'dataType': 'string',
                'value': "degree",
                'required': False,
                'displayName': "Output unit",
                'domain': ('degree', 'radian', 'percent'),
                'description': ("Unit of the output raster.")
            }
        ]

    def getConfiguration(self, **scalars):
        return {
            'compositeRasters': False,
            'inheritProperties': 2 | 4 | 8,         # Inherit all but the pixel type and NoData from the input raster dataset.
            'invalidateProperties': 2 | 4 | 8,          # Invalidate the statistics and histogram on the parent dataset because the pixel values are modified.
            'inputMask': True,
            'resampling': False,
            'padding': 1                                # An input raster mask is not needed in .updatePixels() because the inherited area of NoData are used instead.
        }

    def updateRasterInfo(self, **kwargs):
        kwargs['output_info']['bandCount'] = 1
        r = kwargs['raster_info']
        kwargs['output_info']['noData'] = self.assignNoData(r['pixelType']) if not(r['noData']) else r['noData']
        kwargs['output_info']['pixelType'] = 'f4'
        kwargs['output_info']['histogram'] = ()
        self.prepare(ve_factor=kwargs.get('ve_factor'), output_unit=kwargs.get("output_unit"))
        return kwargs

    def updatePixels(self, tlc, shape, props, **pixelBlocks):
        dem = np.array(pixelBlocks['raster_pixels'], dtype='f4', copy=False)[0]                     # Input pixel array.
        m = np.array(pixelBlocks['raster_mask'], dtype='u1', copy=False)[0]                         # Input raster mask.
        self.noData = self.assignNoData(props['pixelType']) if not(props['noData']) else props['noData']
        dem = np.where(np.not_equal(dem, self.noData), dem, dem)

        pixel_size = props['cellSize']
        if (pixel_size[0] <= 0) | (pixel_size[1] <= 0):
            raise Exception("Input raster cell size is invalid.")

        slope, aspect = RVT_vis_fn.slope_aspect(input_DEM_arr=dem, resolution_x=pixel_size[0],
                                                resolution_y=pixel_size[1], ve_factor=self.ve_factor,
                                                is_padding_applied=True, output_units=self.output_unit)

        pixelBlocks['output_pixels'] = slope.astype(props['pixelType'])
        pixelBlocks['output_mask'] = \
            m[:-2, :-2]  & m[1:-1, :-2]  & m[2:, :-2]  \
            & m[:-2, 1:-1] & m[1:-1, 1:-1] & m[2:, 1:-1] \
            & m[:-2, 2:] & m[1:-1, 2:] & m[2:, 2:]

        return pixelBlocks

    def assignNoData(self, pixelType):
                                                        # assign noData depending on input pixelType
        if pixelType == 'f4':
            return np.array([-3.4028235e+038, ])        # float 32 bit
        elif pixelType == 'i4':
            return np.array([-65536, ])                 # signed integer 32 bit
        elif pixelType == 'i2':
            return np.array([-32768, ])                 # signed integer 16 bit
        elif pixelType == 'i1':
            return np.array([-256, ])                   # signed integer 8 bit
        elif pixelType == 'u4':
            return np.array([65535, ])                  # unsigned integer 32 bit
        elif pixelType == 'u2':
            return np.array([32767, ])                  # unsigned integer 16 bit
        elif pixelType == 'u1':
            return np.array([255, ])                    # unsigned integer 8 bit

    def prepare(self, ve_factor=1, output_unit="degree"):
        self.ve_factor = ve_factor
        self.output_unit = output_unit
