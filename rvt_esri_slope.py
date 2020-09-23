"""
NAME:
    RVT slope esri raster function
    rvt_py, rvt.vis.slope

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

import numpy as np
import rvt.vis


class RVTSlope():

    def __init__(self):
        self.name = "RVT slope"
        self.description = "Calculates slope(gradient)."
        # default values
        self.ve_factor = 1.
        self.output_unit = "degree"
        self.noData = None

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
                'value': self.ve_factor,
                'required': False,
                'displayName': "Ve-factor",
                'description': ("Vertical exaggeration factor (must be greater than 0).")
            },
            {
                'name': 'output_unit',
                'dataType': 'string',
                'value': self.output_unit,
                'required': False,
                'displayName': "Output unit",
                'domain': ('degree', 'radian', 'percent'),
                'description': ("Unit of the output raster.")
            }
        ]

    def getConfiguration(self, **scalars):
        return {
            'compositeRasters': False,
            'inheritProperties': 2 | 4 | 8,
            'invalidateProperties': 2 | 4 | 8,
            'inputMask': False,
            'resampling': False,
            'padding': 1
        }

    def updateRasterInfo(self, **kwargs):
        kwargs['output_info']['bandCount'] = 1
        r = kwargs['raster_info']
        kwargs['output_info']['noData'] = self.assignNoData(r['pixelType']) if not (r['noData']) else r['noData']
        kwargs['output_info']['pixelType'] = 'f4'
        kwargs['output_info']['histogram'] = ()
        kwargs['output_info']['statistics'] = ()
        self.prepare(ve_factor=kwargs.get('ve_factor'), output_unit=kwargs.get("output_unit"))
        return kwargs

    def updatePixels(self, tlc, shape, props, **pixelBlocks):
        dem = np.array(pixelBlocks['raster_pixels'], dtype='f4', copy=False)[0]  # Input pixel array.
        self.noData = self.assignNoData(props['pixelType']) if not (props['noData']) else props['noData']
        dem = np.where(np.not_equal(dem, self.noData), dem, dem)
        pixel_size = props['cellSize']
        if (pixel_size[0] <= 0) | (pixel_size[1] <= 0):
            raise Exception("Input raster cell size is invalid.")

        dict_slp_asp = rvt.vis.slope_aspect(dem=dem, resolution_x=pixel_size[0], resolution_y=pixel_size[1],
                                            ve_factor=self.ve_factor, is_padding_applied=True,
                                            output_units=self.output_unit)
        slope = dict_slp_asp["slope"]

        pixelBlocks['output_pixels'] = slope.astype(props['pixelType'], copy=False)
        return pixelBlocks

    def assignNoData(self, pixelType):
        # assign noData depending on input pixelType
        if pixelType == 'f4':
            return np.array([-3.4028235e+038, ])  # float 32 bit
        elif pixelType == 'i4':
            return np.array([-65536, ])  # signed integer 32 bit
        elif pixelType == 'i2':
            return np.array([-32768, ])  # signed integer 16 bit
        elif pixelType == 'i1':
            return np.array([-256, ])  # signed integer 8 bit
        elif pixelType == 'u4':
            return np.array([65535, ])  # unsigned integer 32 bit
        elif pixelType == 'u2':
            return np.array([32767, ])  # unsigned integer 16 bit
        elif pixelType == 'u1':
            return np.array([255, ])  # unsigned integer 8 bit

    def prepare(self, ve_factor=315, output_unit=35):
        self.ve_factor = ve_factor
        self.output_unit = output_unit
