"""
NAME:
    RVT blender esri raster function
    rvt_py, rvt.blend

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
    University of Ljubljana, Faculty of Civil and Geodetic Engineering
"""

import numpy as np
import rvt.blend
import rvt.default
import rvt.vis
import os


class RVTBlender:
    def __init__(self):
        self.name = "RVT Blender"
        self.description = "Blend visualisations together using norm, min, max, blending mode, opacity."
        # default values
        self.dem_path = None
        self.default_val_file_path = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                                                   "settings/default_settings.json"))
        self.blend_from_file = False
        self.blend_file_path = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                                             "settings/blender_file_example.json"))

        self.layer1_vis = "Sky-View Factor"
        self.layer1_norm = "Value"
        self.layer1_min = 0.7
        self.layer1_max = 1.
        self.layer1_blending_mode = "Multiply"
        self.layer1_opacity = 25.

        self.layer2_vis = "Openness - Positive"
        self.layer2_norm = "Value"
        self.layer2_min = 68.
        self.layer2_max = 93.
        self.layer2_blending_mode = "Overlay"
        self.layer2_opacity = 50.

        self.layer3_vis = "Slope gradient"
        self.layer3_norm = "Value"
        self.layer3_min = 0.
        self.layer3_max = 50.
        self.layer3_blending_mode = "Luminosity"
        self.layer3_opacity = 50.

        self.layer4_vis = "Hillshade"
        self.layer4_norm = "Value"
        self.layer4_min = 0.
        self.layer4_max = 1.
        self.layer4_blending_mode = "Normal"
        self.layer4_opacity = 100.

        self.layer5_vis = "None"
        self.layer5_norm = "Value"
        self.layer5_min = 0.
        self.layer5_max = 0.
        self.layer5_blending_mode = "Normal"
        self.layer5_opacity = 0.

    def getParameterInfo(self):
        return [
            {
                'name': 'raster',
                'dataType': 'raster',
                'value': None,
                'required': True,
                'displayName': "Input Raster",
                'description': "Input raster for which to create the blender."
            },
            {
                'name': 'default_val_file_path',
                'dataType': 'string',
                'value': self.default_val_file_path,
                'required': False,
                'displayName': "Visualisation functions parameter values file path",
                'description': "File with parameter values for visualisation functions. If file doesn't exist"
                               " it takes default values and creates file"
                               " python_file_location/settings/default_settings.json."
            },
            {
                'name': 'blend_from_file',
                'dataType': 'boolean',
                'value': self.blend_from_file,
                'required': False,
                'displayName': "Blend from file",
                'description': "If True it uses blending defined in blending file (blend file path),"
                               " it ignores parameters defined below"
                               " and it saves individual vis function. If False it uses parameters defined down below."
            },
            {
                'name': 'blend_file_path',
                'dataType': 'string',
                'value': self.blend_file_path,
                'required': False,
                'displayName': "Blend file path",
                'description': "Useful only when blend from file is True. File where blending is defined. "
                               "If file doesn't exists and blend from file is true it creates default file in "
                               "python_file_location/settings/blender_file_example.json."
            },
            {
                'name': 'layer1_vis',
                'dataType': 'string',
                'value': self.layer1_vis,
                'required': False,
                'displayName': "Layer 1: Visualisation method.",
                'domain': ("None", "Hillshade", "Multiple directions hillshade", "Slope gradient",
                           "Simple local relief model", "Sky-View Factor", "Anisotropic Sky-View Factor",
                           "Openness - Positive", "Openness - Negative", "Sky illumination", "Local dominance"),
                'description': "Visualisation method (visualisation function)."
            },
            {
                'name': 'layer1_norm',
                'dataType': 'string',
                'value': self.layer1_norm,
                'required': False,
                'displayName': "Layer 1: Norm",
                'domain': ("Value", "Perc"),
                'description': "Normalization"
            },
            {
                'name': 'layer1_min',
                'dataType': 'numeric',
                'value': self.layer1_min,
                'required': False,
                'displayName': "Layer 1: Min",
                'description': "Minimum"
            },
            {
                'name': 'layer1_max',
                'dataType': 'numeric',
                'value': self.layer1_max,
                'required': False,
                'displayName': "Layer 1: Max",
                'description': "Maximum"
            },
            {
                'name': 'layer1_blending_mode',
                'dataType': 'string',
                'value': self.layer1_blending_mode,
                'required': False,
                'displayName': "Layer 1: Blending mode",
                'domain': ("Normal", "Multiply", "Overlay", "Luminosity", "Screen"),
                'description': "Blending mode"
            },
            {
                'name': 'layer1_opacity',
                'dataType': 'numeric',
                'value': self.layer1_opacity,
                'required': False,
                'displayName': "Layer 1: Opacity",
                'description': "Opacity"
            },
            {
                'name': 'layer2_vis',
                'dataType': 'string',
                'value': self.layer2_vis,
                'required': False,
                'displayName': "Layer 2: Visualisation method.",
                'domain': ("None", "Hillshade", "Multiple directions hillshade", "Slope gradient",
                           "Simple local relief model", "Sky-View Factor", "Anisotropic Sky-View Factor",
                           "Openness - Positive", "Openness - Negative", "Sky illumination", "Local dominance"),
                'description': "Visualisation method (visualisation function)."
            },
            {
                'name': 'layer2_norm',
                'dataType': 'string',
                'value': self.layer2_norm,
                'required': False,
                'displayName': "Layer 2: Norm",
                'domain': ("Value", "Perc"),
                'description': "Normalization"
            },
            {
                'name': 'layer2_min',
                'dataType': 'numeric',
                'value': self.layer2_min,
                'required': False,
                'displayName': "Layer 2: Min",
                'description': "Minimum"
            },
            {
                'name': 'layer2_max',
                'dataType': 'numeric',
                'value': self.layer2_max,
                'required': False,
                'displayName': "Layer 2: Max",
                'description': "Maximum"
            },
            {
                'name': 'layer2_blending_mode',
                'dataType': 'string',
                'value': self.layer2_blending_mode,
                'required': False,
                'displayName': "Layer 2: Blending mode",
                'domain': ("Normal", "Multiply", "Overlay", "Luminosity", "Screen"),
                'description': "Blending mode"
            },
            {
                'name': 'layer2_opacity',
                'dataType': 'numeric',
                'value': self.layer2_opacity,
                'required': False,
                'displayName': "Layer 2: Opacity",
                'description': "Opacity"
            },
            {
                'name': 'layer3_vis',
                'dataType': 'string',
                'value': self.layer3_vis,
                'required': False,
                'displayName': "Layer 3: Visualisation method.",
                'domain': ("None", "Hillshade", "Multiple directions hillshade", "Slope gradient",
                           "Simple local relief model", "Sky-View Factor", "Anisotropic Sky-View Factor",
                           "Openness - Positive", "Openness - Negative", "Sky illumination", "Local dominance"),
                'description': "Visualisation method (visualisation function)."
            },
            {
                'name': 'layer3_norm',
                'dataType': 'string',
                'value': self.layer3_norm,
                'required': False,
                'displayName': "Layer 3: Norm",
                'domain': ("Value", "Perc"),
                'description': "Normalization"
            },
            {
                'name': 'layer3_min',
                'dataType': 'numeric',
                'value': self.layer3_min,
                'required': False,
                'displayName': "Layer 3: Min",
                'description': "Minimum"
            },
            {
                'name': 'layer3_max',
                'dataType': 'numeric',
                'value': self.layer3_max,
                'required': False,
                'displayName': "Layer 3: Max",
                'description': "Maximum"
            },
            {
                'name': 'layer3_blending_mode',
                'dataType': 'string',
                'value': self.layer3_blending_mode,
                'required': False,
                'displayName': "Layer 3: Blending mode",
                'domain': ("Normal", "Multiply", "Overlay", "Luminosity", "Screen"),
                'description': "Blending mode"
            },
            {
                'name': 'layer3_opacity',
                'dataType': 'numeric',
                'value': self.layer3_opacity,
                'required': False,
                'displayName': "Layer 3: Opacity",
                'description': "Opacity"
            },
            {
                'name': 'layer4_vis',
                'dataType': 'string',
                'value': self.layer4_vis,
                'required': False,
                'displayName': "Layer 4: Visualisation method.",
                'domain': ("None", "Hillshade", "Multiple directions hillshade", "Slope gradient",
                           "Simple local relief model", "Sky-View Factor", "Anisotropic Sky-View Factor",
                           "Openness - Positive", "Openness - Negative", "Sky illumination", "Local dominance"),
                'description': "Visualisation method (visualisation function)."
            },
            {
                'name': 'layer4_norm',
                'dataType': 'string',
                'value': self.layer4_norm,
                'required': False,
                'displayName': "Layer 4: Norm",
                'domain': ("Value", "Perc"),
                'description': "Normalization"
            },
            {
                'name': 'layer4_min',
                'dataType': 'numeric',
                'value': self.layer4_min,
                'required': False,
                'displayName': "Layer4: Min",
                'description': "Minimum"
            },
            {
                'name': 'layer4_max',
                'dataType': 'numeric',
                'value': self.layer4_max,
                'required': False,
                'displayName': "Layer4: Max",
                'description': "Maximum"
            },
            {
                'name': 'layer4_blending_mode',
                'dataType': 'string',
                'value': self.layer4_blending_mode,
                'required': False,
                'displayName': "Layer 4: Blending mode",
                'domain': ("Normal", "Multiply", "Overlay", "Luminosity", "Screen"),
                'description': "Blending mode"
            },
            {
                'name': 'layer4_opacity',
                'dataType': 'numeric',
                'value': self.layer4_opacity,
                'required': False,
                'displayName': "Layer 4: Opacity",
                'description': "Opacity"
            },
            {
                'name': 'layer5_vis',
                'dataType': 'string',
                'value': self.layer5_vis,
                'required': False,
                'displayName': "Layer 5: Visualisation method.",
                'domain': ("None", "Hillshade", "Multiple directions hillshade", "Slope gradient",
                           "Simple local relief model", "Sky-View Factor", "Anisotropic Sky-View Factor",
                           "Openness - Positive", "Openness - Negative", "Sky illumination", "Local dominance"),
                'description': "Visualisation method (visualisation function)."
            },
            {
                'name': 'layer5_norm',
                'dataType': 'string',
                'value': self.layer5_norm,
                'required': False,
                'displayName': "Layer 5: Norm",
                'domain': ("Value", "Perc"),
                'description': "Normalization"
            },
            {
                'name': 'layer5_min',
                'dataType': 'numeric',
                'value': self.layer5_min,
                'required': False,
                'displayName': "Layer 5: Min",
                'description': "Minimum"
            },
            {
                'name': 'layer5_max',
                'dataType': 'numeric',
                'value': self.layer5_max,
                'required': False,
                'displayName': "Layer 5: Max",
                'description': "Maximum"
            },
            {
                'name': 'layer5_blending_mode',
                'dataType': 'string',
                'value': self.layer5_blending_mode,
                'required': False,
                'displayName': "Layer 5: Blending mode",
                'domain': ("Normal", "Multiply", "Overlay", "Luminosity", "Screen"),
                'description': "Blending mode"
            },
            {
                'name': 'layer5_opacity',
                'dataType': 'numeric',
                'value': self.layer5_opacity,
                'required': False,
                'displayName': "Layer 5: Opacity",
                'description': "Opacity"
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
        self.prepare(default_val_file_path=kwargs.get("default_val_file_path"),
                     blend_from_file=kwargs.get("blend_from_file"),
                     blend_file_path=kwargs.get("blend_file_path"),
                     layer1_vis=kwargs.get("layer1_vis"), layer1_norm=kwargs.get("layer1_norm"),
                     layer1_min=kwargs.get("layer1_min"), layer1_max=kwargs.get("layer1_max"),
                     layer1_blending_mode=kwargs.get("layer1_blending_mode"),
                     layer1_opacity=kwargs.get("layer1_opacity"), layer2_vis=kwargs.get("layer2_vis"),
                     layer2_norm=kwargs.get("layer2_norm"), layer2_min=kwargs.get("layer2_min"),
                     layer2_max=kwargs.get("layer2_max"), layer2_blending_mode=kwargs.get("layer2_blending_mode"),
                     layer2_opacity=kwargs.get("layer2_opacity"), layer3_vis=kwargs.get("layer3_vis"),
                     layer3_norm=kwargs.get("layer3_norm"), layer3_min=kwargs.get("layer3_min"),
                     layer3_max=kwargs.get("layer3_max"), layer3_blending_mode=kwargs.get("layer3_blending_mode"),
                     layer3_opacity=kwargs.get("layer3_opacity"), layer4_vis=kwargs.get("layer4_vis"),
                     layer4_norm=kwargs.get("layer4_norm"), layer4_min=kwargs.get("layer4_min"),
                     layer4_max=kwargs.get("layer4_max"), layer4_blending_mode=kwargs.get("layer4_blending_mode"),
                     layer4_opacity=kwargs.get("layer4_opacity"), layer5_vis=kwargs.get("layer5_vis"),
                     layer5_norm=kwargs.get("layer5_norm"), layer5_min=kwargs.get("layer5_min"),
                     layer5_max=kwargs.get("layer5_max"), layer5_blending_mode=kwargs.get("layer5_blending_mode"),
                     layer5_opacity=kwargs.get("layer5_opacity"))
        return kwargs

    def updatePixels(self, tlc, shape, props, **pixelBlocks):
        dem = np.array(pixelBlocks['raster_pixels'], dtype='f4', copy=False)[0]  # Input pixel array.
        pixel_size = props['cellSize']
        if (pixel_size[0] <= 0) | (pixel_size[1] <= 0):
            raise Exception("Input raster cell size is invalid.")

        # read visualisation parameters values in DefaultValues class
        default = rvt.default.DefaultValues()  # already contains default values
        layers = rvt.blend.BlenderLayers()  # containing layers
        layers.add_dem_arr(dem, pixel_size[0])

        if os.path.isfile(self.default_val_file_path):  # if settings file exists change values from file
            default.read_default_from_file(self.default_val_file_path)

        else:  # if doesn't exist create one with default values
            default_file_out_path = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                                                  "settings/default_settings.json"))
            if not os.path.exists(os.path.dirname(default_file_out_path)):
                os.makedirs(os.path.dirname(default_file_out_path))
            default.save_default_to_file(file_path=default_file_out_path)

        if self.blend_from_file:  # blending from file (blend form file checkbox checked)
            if os.path.isfile(self.blend_file_path):
                layers.build_blender_layers_from_file(file_path=self.blend_file_path)

            else:  # if blender file doesn't exist it creates example, which can be changed
                blend_file_path = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                                                "settings/blender_file_example.json"))
                if os.path.isfile(blend_file_path):
                    pass
                else:
                    if not os.path.exists(os.path.dirname(blend_file_path)):
                        os.makedirs(os.path.dirname(blend_file_path))
                    rvt.blend.create_blender_file_example(file_path=blend_file_path)

        else:  # blending from parameters defined in ArcGIS Pro
            layers.create_layer(self.layer1_vis, self.layer1_norm, self.layer1_min, self.layer1_max,
                                self.layer1_blending_mode, self.layer1_opacity)
            layers.create_layer(self.layer2_vis, self.layer2_norm, self.layer2_min, self.layer2_max,
                                self.layer2_blending_mode, self.layer2_opacity)
            layers.create_layer(self.layer3_vis, self.layer3_norm, self.layer3_min, self.layer3_max,
                                self.layer3_blending_mode, self.layer3_opacity)
            layers.create_layer(self.layer4_vis, self.layer4_norm, self.layer4_min, self.layer4_max,
                                self.layer4_blending_mode, self.layer4_opacity)
            layers.create_layer(self.layer5_vis, self.layer5_norm, self.layer5_min, self.layer5_max,
                                self.layer5_blending_mode, self.layer5_opacity)

        blender = layers.render_all_images(default=default)  # calculate blender
        pixelBlocks['output_pixels'] = blender.astype(props['pixelType'], copy=False)
        return pixelBlocks

    def prepare(self, default_val_file_path="settings/default_settings.json", save_individual_vis_fn=True,
                blend_from_file=False, blend_file_path="settings/blender_file_example.json",
                layer1_vis="Sky-View Factor", layer1_norm="Value", layer1_min=0.7, layer1_max=1.,
                layer1_blending_mode="Multiply", layer1_opacity=25., layer2_vis="Openness - Positive",
                layer2_norm="Value", layer2_min=68., layer2_max=93., layer2_blending_mode="Overlay",
                layer2_opacity=50., layer3_vis="Slope gradient", layer3_norm="Value", layer3_min=0.,
                layer3_max=50., layer3_blending_mode="Luminosity", layer3_opacity=50., layer4_vis="Hillshade",
                layer4_norm="Value", layer4_min=0., layer4_max=1., layer4_blending_mode="Normal",
                layer4_opacity=100., layer5_vis="None", layer5_norm="Value", layer5_min=0, layer5_max=0,
                layer5_blending_mode="Normal", layer5_opacity=0):

        self.default_val_file_path = default_val_file_path
        self.blend_from_file = blend_from_file
        self.blend_file_path = blend_file_path
        if layer1_vis == "None":
            self.layer1_vis = None
        else:
            self.layer1_vis = layer1_vis
        self.layer1_norm = layer1_norm
        self.layer1_min = float(layer1_min)
        self.layer1_max = float(layer1_max)
        self.layer1_blending_mode = layer1_blending_mode
        self.layer1_opacity = int(layer1_opacity)
        if layer2_vis == "None":
            self.layer2_vis = None
        else:
            self.layer2_vis = layer2_vis
        self.layer2_norm = layer2_norm
        self.layer2_min = float(layer2_min)
        self.layer2_max = float(layer2_max)
        self.layer2_blending_mode = layer2_blending_mode
        self.layer2_opacity = int(layer2_opacity)
        if layer3_vis == "None":
            self.layer3_vis = None
        else:
            self.layer3_vis = layer3_vis
        self.layer3_norm = layer3_norm
        self.layer3_min = float(layer3_min)
        self.layer3_max = float(layer3_max)
        self.layer3_blending_mode = layer3_blending_mode
        self.layer3_opacity = int(layer3_opacity)
        if layer4_vis == "None":
            self.layer4_vis = None
        else:
            self.layer4_vis = layer4_vis
        self.layer4_norm = layer4_norm
        self.layer4_min = float(layer4_min)
        self.layer4_max = float(layer4_max)
        self.layer4_blending_mode = layer4_blending_mode
        self.layer4_opacity = int(layer4_opacity)
        if layer5_vis == "None":
            self.layer5_vis = None
        else:
            self.layer5_vis = layer5_vis
        self.layer5_norm = layer5_norm
        self.layer5_min = float(layer5_min)
        self.layer5_max = float(layer5_max)
        self.layer5_blending_mode = layer5_blending_mode
        self.layer5_opacity = int(layer5_opacity)
