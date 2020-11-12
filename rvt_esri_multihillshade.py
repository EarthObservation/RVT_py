"""
Relief Visualization Toolbox – Visualization Functions

RVT multiple directions hillshade esri raster function
rvt_py, rvt.vis.multi_hillshade

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


class RVTMultiHillshade:
    def __init__(self):
        self.name = "RVT multi hillshade"
        self.description = "Calculates multiple directions hillshade."
        # default values
        self.nr_directions = 16.
        self.elevation = 35.
        self.calc_8_bit = False
        self.padding = 1

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
            {
                'name': 'calc_8_bit',
                'dataType': 'boolean',
                'value': self.calc_8_bit,
                'required': False,
                'displayName': "Calculate 8-bit",
                'description': "If True it only calculates 8-bit (nr_directions doesn't matter),"
                               " in 3 directions (sun_azimuith = 315, 22.5, 90)"
            },
            {
                'name': 'nr_directions',
                'dataType': 'numeric',
                'value': self.nr_directions,
                'required': False,
                'displayName': "Number of directions",
                'description': "Number of directions for sun azimuth angle (clockwise from North)."
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
        kwargs['output_info']['noData'] = -3.4028235e+038
        kwargs['output_info']['pixelType'] = 'f4'
        kwargs['output_info']['histogram'] = ()
        kwargs['output_info']['statistics'] = int(self.nr_directions) * ({'minimum': -1, 'maximum': 1},)
        self.prepare(nr_directions=kwargs.get('nr_directions'), elevation=kwargs.get("sun_elevation"),
                     calc_8_bit=kwargs.get("calc_8_bit"))
        if self.calc_8_bit:
            kwargs['output_info']['bandCount'] = 3
        else:
            kwargs['output_info']['bandCount'] = int(self.nr_directions)
        return kwargs

    def updatePixels(self, tlc, shape, props, **pixelBlocks):
        dem = np.array(pixelBlocks['raster_pixels'], dtype='f4', copy=False)[0]  # Input pixel array.
        pixel_size = props['cellSize']
        if (pixel_size[0] <= 0) | (pixel_size[1] <= 0):
            raise Exception("Input raster cell size is invalid.")
        dict_slp_asp = rvt.vis.slope_aspect(dem=dem, resolution_x=pixel_size[0], resolution_y=pixel_size[1],
                                            ve_factor=1)
        if self.calc_8_bit:  # calc 8 bit
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
            hillshade_rgb = hillshade_rgb[:, self.padding:-self.padding, self.padding:-self.padding]  # remove padding
            pixelBlocks['output_pixels'] = hillshade_rgb.astype(props['pixelType'], copy=False)
        else:  # calc nr_directions
            multihillshade = rvt.vis.multi_hillshade(dem=dem, resolution_x=pixel_size[0], resolution_y=pixel_size[1],
                                                     nr_directions=self.nr_directions, sun_elevation=self.elevation,
                                                     slope=dict_slp_asp["slope"], aspect=dict_slp_asp["aspect"])
            multihillshade = multihillshade[:, self.padding:-self.padding, self.padding:-self.padding]  # remove padding
            pixelBlocks['output_pixels'] = multihillshade.astype(props['pixelType'], copy=False)
        return pixelBlocks

    def prepare(self, nr_directions=16, elevation=35, calc_8_bit=False):
        self.nr_directions = int(nr_directions)
        self.elevation = elevation
        self.calc_8_bit = calc_8_bit
