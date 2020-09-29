"""
NAME:
    RVT multiple directions hillshade esri raster function
    rvt_py, rvt.vis.multi_hillshade

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

# TODO: If nr_band > 10 function doesn't work, find out from Esri developers why
# Currently working 8bit, 3 bands (azimuth 315 for red, 22.5 for green, and 90 for the blue band)
# TODO: when solved remove temporarily working 8bit (azimuth 315 for red, 22.5 for green, and 90 for the blue band)

class RVTMultiHillshade():
    def __init__(self):
        self.name = "RVT multi hillshade"
        self.description = "Calculates multiple directions hillshade."
        # default values
        #self.nr_directions = 16.
        self.elevation = 35.

    def getParameterInfo(self):
        return [
            {
                'name': 'raster',
                'dataType': 'raster',
                'value': None,
                'required': True,
                'displayName': "Input Raster",
                'description': "Input raster for which to create the multiple directions hillshade map."
            },
            # {
            #     'name': 'nr_directions',
            #     'dataType': 'numeric',
            #     'value': self.nr_directions,
            #     'required': False,
            #     'displayName': "Number of directions",
            #     'description': "Number of directions for sun azimuth angle (clockwise from North)."
            # },
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
            'padding': 0
        }

    def updateRasterInfo(self, **kwargs):
        kwargs['output_info']['bandCount'] = 3 #int(kwargs.get('nr_directions'))
        r = kwargs['raster_info']
        kwargs['output_info']['noData'] = -3.4028235e+038
        kwargs['output_info']['pixelType'] = 'f4'
        kwargs['output_info']['histogram'] = ()
        kwargs['output_info']['statistics'] = ()
        #self.prepare(nr_directions=kwargs.get('nr_directions'), elevation=kwargs.get("sun_elevation"))
        self.prepare(elevation=kwargs.get("sun_elevation"))
        return kwargs

    def updatePixels(self, tlc, shape, props, **pixelBlocks):
        dem = np.array(pixelBlocks['raster_pixels'], dtype='f4', copy=False)[0]  # Input pixel array.
        pixel_size = props['cellSize']
        if (pixel_size[0] <= 0) | (pixel_size[1] <= 0):
            raise Exception("Input raster cell size is invalid.")
        dict_slp_asp = rvt.vis.slope_aspect(dem=dem, resolution_x=pixel_size[0], resolution_y=pixel_size[1],
                                            ve_factor=1)
        hillshade_r = rvt.vis.hillshade(dem=dem, resolution_x=pixel_size[0], resolution_y=pixel_size[1],
                                        sun_azimuth=315, sun_elevation=self.elevation,
                                        slope=dict_slp_asp["slope"], aspect=dict_slp_asp["aspect"])
        hillshade_g = rvt.vis.hillshade(dem=dem, resolution_x=pixel_size[0], resolution_y=pixel_size[1],
                                        sun_azimuth=22.5, sun_elevation=self.elevation,
                                        slope=dict_slp_asp["slope"], aspect=dict_slp_asp["aspect"])
        hillshade_b = rvt.vis.hillshade(dem=dem, resolution_x=pixel_size[0], resolution_y=pixel_size[1],
                                        sun_azimuth=90, sun_elevation=self.elevation,
                                        slope=dict_slp_asp["slope"], aspect=dict_slp_asp["aspect"])
        hillshade_rgb = np.array([hillshade_r, hillshade_g, hillshade_b])
        pixelBlocks['output_pixels'] = hillshade_rgb.astype(props['pixelType'], copy=False)
        # multihillshade = rvt.vis.multi_hillshade(dem=dem, resolution_x=pixel_size[0], resolution_y=pixel_size[1],
        #                                          nr_directions=self.nr_directions, sun_elevation=self.elevation,
        #                                          slope=dict_slp_asp["slope"], aspect=dict_slp_asp["aspect"])
        # pixelBlocks['output_pixels'] = multihillshade.astype(props['pixelType'], copy=False)
        return pixelBlocks

    # def prepare(self, nr_directions=16, elevation=35):
    #     self.nr_directions = int(nr_directions)
    #     self.elevation = elevation
    def prepare(self, elevation=35):
        self.elevation = elevation
