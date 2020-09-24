"""
NAME:
    RVT simple local relief model esri raster function
    rvt_py, rvt.vis.slrm

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


class RVTSlrm():

    def __init__(self):
        self.name = "RVT simple local relief model."
        self.description = "Calculates simple local relief model."
        # default values
        self.radius_cell = 20.

    def getParameterInfo(self):
        return [
            {
                'name': 'raster',
                'dataType': 'raster',
                'value': None,
                'required': True,
                'displayName': "Input Raster",
                'description': "Input raster for which to create the simple local relief model."
            },
            {
                'name': 'radius_cell',
                'dataType': 'numeric',
                'value': self.radius_cell,
                'required': False,
                'displayName': "Radius cell",
                'description': "Radius for trend assessment [pixels]."
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
        self.prepare(radius_cell=kwargs.get('radius_cell'))
        return kwargs

    def updatePixels(self, tlc, shape, props, **pixelBlocks):
        dem = np.array(pixelBlocks['raster_pixels'], dtype='f4', copy=False)[0]  # Input pixel array.
        slrm = rvt.vis.slrm(dem=dem, radius_cell=self.radius_cell)
        pixelBlocks['output_pixels'] = slrm.astype(props['pixelType'], copy=False)
        return pixelBlocks

    def prepare(self, radius_cell=20.):
        self.radius_cell = int(radius_cell)

