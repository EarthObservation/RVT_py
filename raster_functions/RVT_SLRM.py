"""
NAME:
    RVT Simple local relief model - SLRM, ESRI ArcGis Pro

DESCRIPTION:
    https://github.com/Esri/raster-functions/wiki/PythonRasterFunction#anatomy-of-a-python-raster-function
    Python raster function for ESRI ArcGis which calculates simple local relief model.

INPUTS:
    raster      - input raster
    radius_cell - Radius for trend assessment [pixels]

OUTPUTS:
    slrm

KEYWORDS:
    /

DEPENDENCIES:
    RVT_vin_fn.SLRM

AUTHOR:
    RVT:
        Klemen Zaksek (klemen.zaksek@zmaw.de)
    RVT_py:
        Žiga Maroh (ziga.maroh@icloud.com)

MODIFICATION HISTORY:
    RVT:
        1.0 Written by Klemen Zaksek, 2013.
    RVT_py:
        1.0 Written by Žiga Maroh, 2020.
"""

import numpy as np
import RVT_vis_fn


class RVT_SLRM():

    def __init(self):
        self.name = "RVT SLRM"
        self.description = "Calculates simple local relief model."
        self.prepare()

    def getParameterInfo(self):
        return [
            {
                'name': 'raster',
                'dataType': 'raster',
                'value': None,
                'required': True,
                'displayName': "Input Raster",
                'description': "Input raster for which to create the SLRM."
            },
            {
                'name': 'radius_cell',
                'dataType': 'numeric',
                'value': 20.,
                'required': False,
                'displayName': "Radius [pixels]",
                'description': "Radius for trend assessment [pixels], allowed values 10-50. "
                               "If radius less than 10 program changes it to 10, "
                               "if more than 50 program changes it to 50"
            }
        ]

    def getConfiguration(self, **scalars):
        return {
            'compositeRasters': False,
            'inheritProperties': 2 | 4 | 8,  # Inherit all but the pixel type and NoData from the input raster dataset.
            'invalidateProperties': 2 | 4 | 8,
            # Invalidate the statistics and histogram on the parent dataset because the pixel values are modified.
            'inputMask': True,
            'resampling': False
            # An input raster mask is not needed in .updatePixels() because the inherited area of NoData are used instead.
        }

    def updateRasterInfo(self, **kwargs):
        kwargs['output_info']['bandCount'] = 1
        r = kwargs['raster_info']
        kwargs['output_info']['noData'] = self.assignNoData(r['pixelType']) if not (r['noData']) else r['noData']
        kwargs['output_info']['pixelType'] = 'f4'
        kwargs['output_info']['histogram'] = ()
        kwargs['output_info']['statistics'] = ()
        self.prepare(radius_cell=kwargs.get('radius_cell'))
        return kwargs

    def updatePixels(self, tlc, shape, props, **pixelBlocks):
        dem = np.array(pixelBlocks['raster_pixels'], dtype='f4', copy=False)[0]  # Input pixel array.
        m = np.array(pixelBlocks['raster_mask'], dtype='u1', copy=False)[0]  # Input raster mask.
        self.noData = self.assignNoData(props['pixelType']) if not (props['noData']) else props['noData']
        dem = np.where(np.not_equal(dem, self.noData), dem, dem)

        slrm = RVT_vis_fn.SLRM(den=dem, radius_cell=self.radius_cell)

        pixelBlocks['output_pixels'] = slrm.astype(props['pixelType'], copy=False)
        pixelBlocks['output_mask'] = m

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

    def prepare(self, radius_cell=20):
        if radius_cell < 10:
            self.radius_cell = 10
        elif radius_cell > 50:
            self.radius_cell = 50
        else:
            self.radius_cell = int(radius_cell)

