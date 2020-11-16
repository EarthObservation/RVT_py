"""
Relief Visualization Toolbox – Visualization Functions

RVT blend esri raster function
rvt_py, rvt.blend_func.blend_images, rvt.blend_func.render_images

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
import rvt.blend_func


class RVTBlend:
    def __init__(self):
        self.name = "RVT blend"
        self.description = "Blend and render two images together."
        # default values
        self.blend_mode = "normal"
        self.opacity = 100.


    def getParameterInfo(self):
        return [
            {
                'name': 'top_raster',
                'dataType': 'raster',
                'value': None,
                'required': True,
                'displayName': "Input top Raster",
                'description': "Input top raster on which we apply opacity, blend and render with background raster."
            },
            {
                'name': 'background_raster',
                'dataType': 'raster',
                'value': None,
                'required': True,
                'displayName': "Input background Raster",
                'description': "Input background raster which we blend and render with top raster."
            },
            {
                'name': 'blend_mode',
                'dataType': 'string',
                'value': self.blend_mode,
                'required': True,
                'displayName': "Blend mode",
                'domain': ('normal', 'multiply', 'overlay', 'luminosity', 'screen'),
                'description': "Blending mode for blending top and background raster together."
            },
            {
                'name': 'opacity',
                'dataType': 'numeric',
                'value': self.opacity,
                'required': True,
                'displayName': "Opacity",
                'description': "Opacity in percent to apply on top raster (0-100)."
            }
        ]

    def getConfiguration(self, **scalars):
        return {
            'compositeRasters': False,
            'inheritProperties': 4 | 8,
            'invalidateProperties': 2 | 4 | 8,
            'inputMask': False,
            'resampling': False,
            'padding': 0
        }

    def updateRasterInfo(self, **kwargs):
        kwargs['output_info']['bandCount'] = 1
        # r = kwargs['raster_info']
        kwargs['output_info']['noData'] = -3.4028235e+038
        kwargs['output_info']['pixelType'] = 'f4'
        kwargs['output_info']['histogram'] = ()
        kwargs['output_info']['statistics'] = ({'minimum': 0.0, 'maximum': 1.0}, )
        self.prepare(blend_mode=kwargs.get('blend_mode'), opacity=kwargs.get("opacity"))
        return kwargs

    def updatePixels(self, tlc, shape, props, **pixelBlocks):
        top_raster = np.array(pixelBlocks['top_raster_pixels'], dtype='f4', copy=False)[0]
        background_raster = np.array(pixelBlocks['background_raster_pixels'], dtype='f4', copy=False)[0]
        pixel_size = props['cellSize']
        if (pixel_size[0] <= 0) | (pixel_size[1] <= 0):
            raise Exception("Input raster cell size is invalid.")

        top_raster = rvt.blend_func.blend_images(blend_mode=self.blend_mode, active=top_raster,
                                                 background=background_raster)
        rendered_image = rvt.blend_func.render_images(active=top_raster, background=background_raster,
                                                      opacity=self.opacity)

        pixelBlocks['output_pixels'] = rendered_image.astype(props['pixelType'], copy=False)
        return pixelBlocks

    def prepare(self, blend_mode="normal", opacity=100):
        self.blend_mode = blend_mode
        if opacity > 100:
            self.opacity = 100
        elif opacity < 100:
            self.opacity = 0
        else:
            self.opacity = int(opacity)


