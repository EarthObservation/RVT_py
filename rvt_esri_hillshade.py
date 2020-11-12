"""
Relief Visualization Toolbox – Visualization Functions

RVT hillshade esri raster function
rvt_py, rvt.vis.hillshade

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


class RVTHillshade:
    def __init__(self):
        self.name = "RVT hillshade"
        self.description = "Calculates hillshade."
        # default values
        self.azimuth = 315.
        self.elevation = 35.
        self.padding = 1

    def getParameterInfo(self):
        return [
            {
                'name': 'raster',
                'dataType': 'raster',
                'value': None,
                'required': True,
                'displayName': "Input Raster",
                'description': "Input raster for which to create the hillshade map."
            },
            {
                'name': 'sun_azimuth',
                'dataType': 'numeric',
                'value': self.azimuth,
                'required': False,
                'displayName': "Sun azimuth",
                'description': "Solar azimuth angle (clockwise from North) in degrees."
            },
            {
                'name': 'sun_elevation',
                'dataType': 'numeric',
                'value': self.elevation,
                'required': False,
                'displayName': "Sun elevation",
                'description': "Solar vertical angle (above the horizon) in degrees."
            }
        ]

    def getConfiguration(self, **scalars):
        return {
            'compositeRasters': False,
            'inheritProperties': 2 | 4 | 8,
            'invalidateProperties': 2 | 4 | 8,
            'inputMask': False,
            'resampling': False,
            'padding': self.padding
        }

    def updateRasterInfo(self, **kwargs):
        kwargs['output_info']['bandCount'] = 1
        r = kwargs['raster_info']
        kwargs['output_info']['noData'] = -3.4028235e+038
        kwargs['output_info']['pixelType'] = 'f4'
        kwargs['output_info']['histogram'] = ()
        kwargs['output_info']['statistics'] = ()
        self.prepare(azimuth=kwargs.get('sun_azimuth'), elevation=kwargs.get("sun_elevation"))
        return kwargs

    def updatePixels(self, tlc, shape, props, **pixelBlocks):
        dem = np.array(pixelBlocks['raster_pixels'], dtype='f4', copy=False)[0]  # Input pixel array.
        pixel_size = props['cellSize']
        if (pixel_size[0] <= 0) | (pixel_size[1] <= 0):
            raise Exception("Input raster cell size is invalid.")
        hillshade = rvt.vis.hillshade(dem=dem, resolution_x=pixel_size[0],
                                      resolution_y=pixel_size[1], sun_azimuth=self.azimuth,
                                      sun_elevation=self.elevation)
        hillshade = hillshade[self.padding:-self.padding, self.padding:-self.padding]  # remove padding
        pixelBlocks['output_pixels'] = hillshade.astype(props['pixelType'], copy=False)
        return pixelBlocks

    def prepare(self, azimuth=315, elevation=35):
        self.azimuth = azimuth
        self.elevation = elevation
