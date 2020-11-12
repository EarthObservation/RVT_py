"""
Relief Visualization Toolbox – Visualization Functions

RVT Sky-view factor esri raster function
rvt_py, rvt.vis.sky_view_factor

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


class RVTSvf:
    def __init__(self):
        self.name = "RVT Sky-view factor"
        self.description = "Calculates Sky-view factor."
        # default values
        self.nr_directions = 16.
        self.max_rad = 10.
        self.noise = "0-don't remove"
        self.padding = int(self.max_rad/2)

    def getParameterInfo(self):
        return [
            {
                'name': 'raster',
                'dataType': 'raster',
                'value': None,
                'required': True,
                'displayName': "Input Raster",
                'description': "Input raster for which to create the sky-view factor map."
            },
            {
                'name': 'nr_directions',
                'dataType': 'numeric',
                'value': self.nr_directions,
                'required': False,
                'displayName': "Number of directions",
                'description': "Number of directions."
            },
            {
                'name': 'max_rad',
                'dataType': 'numeric',
                'value': self.max_rad,
                'required': False,
                'displayName': "Max radius",
                'description': "Maximal search radius in pixels."
            },
            {
                'name': 'noise_remove',
                'dataType': 'string',
                'value': self.noise,
                'required': False,
                'displayName': "Noise removal",
                'domain': ("0-don't remove", "1-low", "2-med", "3-high"),
                'description': ("The level of noise remove (0-don't remove, 1-low, 2-med, 3-high).")
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
        self.prepare(nr_directions=kwargs.get('nr_directions'), max_rad=kwargs.get("max_rad"),
                     noise=kwargs.get("noise_remove"))
        return kwargs

    def updatePixels(self, tlc, shape, props, **pixelBlocks):
        dem = np.array(pixelBlocks['raster_pixels'], dtype='f4', copy=False)[0]  # Input pixel array.
        pixel_size = props['cellSize']

        if (pixel_size[0] <= 0) | (pixel_size[1] <= 0):
            raise Exception("Input raster cell size is invalid.")

        dict_svf = rvt.vis.sky_view_factor(dem=dem, resolution=pixel_size[0], compute_svf=True, compute_asvf=False,
                                           compute_opns=False, svf_n_dir=self.nr_directions, svf_r_max=self.max_rad,
                                           svf_noise=self.noise)
        svf = dict_svf["svf"][self.padding:-self.padding, self.padding:-self.padding]  # remove padding

        pixelBlocks['output_pixels'] = svf.astype(props['pixelType'], copy=False)
        return pixelBlocks

    def prepare(self, nr_directions=16, max_rad=10, noise="0"):
        self.nr_directions = int(nr_directions)
        self.max_rad = int(max_rad)
        self.noise = int(noise[0])
        self.padding = int(max_rad/2)
