"""
Relief Visualization Toolbox – Visualization Functions

RVT slope esri raster function
rvt_py, rvt.vis.slope

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

import numpy as np
import rvt.vis


class RVTSlope():
    def __init__(self):
        self.name = "RVT slope"
        self.description = "Calculates slope(gradient)."
        # default values
        self.ve_factor = 1.
        self.output_unit = "degree"

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
            'padding': 0
        }

    def updateRasterInfo(self, **kwargs):
        kwargs['output_info']['bandCount'] = 1
        r = kwargs['raster_info']
        kwargs['output_info']['noData'] = -3.4028235e+038
        kwargs['output_info']['pixelType'] = 'f4'
        kwargs['output_info']['histogram'] = ()
        kwargs['output_info']['statistics'] = ()
        self.prepare(ve_factor=kwargs.get('ve_factor'), output_unit=kwargs.get("output_unit"))
        return kwargs

    def updatePixels(self, tlc, shape, props, **pixelBlocks):
        dem = np.array(pixelBlocks['raster_pixels'], dtype='f4', copy=False)[0]  # Input pixel array.
        pixel_size = props['cellSize']
        if (pixel_size[0] <= 0) | (pixel_size[1] <= 0):
            raise Exception("Input raster cell size is invalid.")

        dict_slp_asp = rvt.vis.slope_aspect(dem=dem, resolution_x=pixel_size[0], resolution_y=pixel_size[1],
                                            ve_factor=self.ve_factor, output_units=self.output_unit)
        slope = dict_slp_asp["slope"]

        pixelBlocks['output_pixels'] = slope.astype(props['pixelType'], copy=False)
        return pixelBlocks

    def prepare(self, ve_factor=315, output_unit=35):
        self.ve_factor = ve_factor
        self.output_unit = output_unit
