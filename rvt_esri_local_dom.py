"""
NAME:
    RVT local dominance esri raster function
    rvt_py, rvt.vis.local_dominance

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

# NOT WORKING!!! Needs to be debugged

class RVTLocalDominance:
    def __init__(self):
        self.name = "RVT Local dominance"
        self.description = "Calculates Local dominance."
        # default values
        self.min_rad = 10.
        self.max_rad = 20.
        self.rad_inc = 1.
        self.anglr_res = 15.
        self.observer_h = 1.7

    def getParameterInfo(self):
        return [
            {
                'name': 'raster',
                'dataType': 'raster',
                'value': None,
                'required': True,
                'displayName': "Input Raster",
                'description': "Input raster for which to create the local dominance map."
            },
            {
                'name': 'min_rad',
                'dataType': 'numeric',
                'value': self.min_rad,
                'required': False,
                'displayName': "Minimum radial distance",
                'description': "Minimum radial distance (in pixels) at which the algorithm starts with visualization"
                               " computation."
            },
            {
                'name': 'max_rad',
                'dataType': 'numeric',
                'value': self.max_rad,
                'required': False,
                'displayName': "Maximum radial distance",
                'description': "Maximum radial distance (in pixels) at which the algorithm ends with visualization"
                               " computation."
            },
            {
                'name': 'rad_inc',
                'dataType': 'numeric',
                'value': self.rad_inc,
                'required': False,
                'displayName': "Radial distance steps",
                'description': "Radial distance steps in pixels"
            },
            {
                'name': 'anglr_res',
                'dataType': 'numeric',
                'value': self.anglr_res,
                'required': False,
                'displayName': "Angular resolution",
                'description': "Angular step for determination of number of angular directions."
            },
            {
                'name': 'observer_h',
                'dataType': 'numeric',
                'value': self.observer_h,
                'required': False,
                'displayName': "Observer height",
                'description': "Height at which we observe the terrain in meters."
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
        self.prepare(min_rad=kwargs.get('min_rad'), max_rad=kwargs.get("max_rad"), rad_inc=kwargs.get("rad_inc"),
                     anglr_res=kwargs.get("anglr_res"), observer_h=kwargs.get("observer_h"))
        return kwargs

    def updatePixels(self, tlc, shape, props, **pixelBlocks):
        dem = np.array(pixelBlocks['raster_pixels'], dtype='f4', copy=False)[0]  # Input pixel array.
        pixel_size = props['cellSize']
        if (pixel_size[0] <= 0) | (pixel_size[1] <= 0):
            raise Exception("Input raster cell size is invalid.")
        local_dominance = rvt.vis.local_dominance(dem=dem, min_rad=self.min_rad, max_rad=self.max_rad,
                                                  rad_inc=self.rad_inc, angular_res=self.anglr_res,
                                                  observer_height=self.observer_h)
        pixelBlocks['output_pixels'] = local_dominance.astype(props['pixelType'], copy=False)
        return pixelBlocks

    def prepare(self, min_rad=10, max_rad=20, rad_inc=1, anglr_res=15, observer_h=1.7):
        self.min_rad = int(min_rad)
        self.max_rad = int(max_rad)
        self.rad_inc = int(rad_inc)
        self.anglr_res = int(anglr_res)
        self.observer_h = float(observer_h)
