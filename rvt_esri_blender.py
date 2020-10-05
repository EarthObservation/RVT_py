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
    Space-SI
    University of Ljubljana, Faculty of Civil and Geodetic Engineering
"""

import numpy as np
import rvt.blend
import rvt.default
import rvt.vis
import os


# Working only when save_individual_vis_fn=No
# Needs to be debugged for when save_individual_vis_fn=Yes
# Find out problems with for loop why layers contains more than 5 layers (layers.add_layer(layer))

class RVTBlender:
    def __init__(self):
        self.name = "RVT Blender"
        self.description = "Blend visualisations together using norm, min, max, blending mode, opacity."
        # default values
        self.dem_path = None
        self.default_val_file_path = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                                                   "settings/default_settings.txt"))
        self.save_individual_vis_fn = "Yes"  # in prepare changed to True False

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
                'name': 'save_individual_vis_fn',
                'dataType': 'string',
                'value': self.save_individual_vis_fn,
                'required': False,
                'displayName': "Save individual visualisation function.",
                'domain': ("No", "Yes"),
                'description': "If Yes (works faster) it saves individual visualization function (layer) in raster"
                               " directory. If No (works very slow or crashes when"
                               " dealing with big raster), it holds specific raster function result (raster) in memory."
            },
            {
                'name': 'default_val_file_path',
                'dataType': 'string',
                'value': self.default_val_file_path,
                'required': False,
                'displayName': "Visualisation functions parameter values file path",
                'description': "File with parameter values for visualisation functions (if file doesn't exist"
                               " it takes default values and creates file /settings/default_settings.txt)."
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
        # self.dem_path = str(kwargs.get("_cachedproperties").get("tifftag_documentname"))  # ASK ESRI how to get raster path?
        self.dem_path = "D:/RVT_py/RVT_py/test_data/TM1_564_146.tif"
        self.prepare(default_val_file_path=kwargs.get("default_val_file_path"),
                     save_individual_vis_fn=kwargs.get("save_individual_vis_fn"),
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

        if os.path.isfile(self.default_val_file_path):  # if settings file exists change values from file
            default.read_default_from_file(self.default_val_file_path)

        else:  # if doesn't exists create one with default values
            default_file_out_path = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                                                  "settings/default_settings.txt"))
            if not os.path.exists(os.path.dirname(default_file_out_path)):
                os.makedirs(os.path.dirname(default_file_out_path))
            default.save_default_to_file(file_path=default_file_out_path)

        vis_list = [self.layer1_vis, self.layer2_vis, self.layer3_vis, self.layer4_vis, self.layer5_vis]
        norm_list = [self.layer1_norm, self.layer2_norm, self.layer3_norm, self.layer4_norm, self.layer5_norm]
        min_list = [self.layer1_min, self.layer2_min, self.layer3_min, self.layer4_min, self.layer5_min]
        max_list = [self.layer1_max, self.layer2_max, self.layer3_max, self.layer4_max, self.layer5_max]
        blending_mode_list = [self.layer1_blending_mode, self.layer2_blending_mode, self.layer3_blending_mode,
                              self.layer4_blending_mode, self.layer5_blending_mode]
        opacity_list = [self.layer1_opacity, self.layer2_opacity, self.layer3_opacity, self.layer4_opacity,
                        self.layer5_opacity]
        layers = rvt.blend.BlenderLayers()  # containing layers

        for i_lyr in range(5):  # iterate trough all layers and saves them in layers (rvt.blend.BlenderLayers())
            vis_method = vis_list[i_lyr]
            normalization = norm_list[i_lyr]
            minimum = min_list[i_lyr]
            maximum = max_list[i_lyr]
            blend_mode = blending_mode_list[i_lyr]
            opacity = opacity_list[i_lyr]
            image = None
            image_path = None
            if vis_method is not None:
                if vis_method.lower() == "slope gradient":
                    if self.save_individual_vis_fn:  # save individually visualisation function output as raster (tif)
                        default.save_slope(self.dem_path)
                        image_path = default.get_slope_path(self.dem_path)
                    else:
                        image = default.get_slope(dem_arr=dem, resolution_x=pixel_size[0],
                                                  resolution_y=pixel_size[1])["slope"]
                elif vis_method.lower() == "hillshade":
                    if self.save_individual_vis_fn:  # save individually visualisation function output as raster (tif)
                        default.save_hillshade(self.dem_path)
                        image_path = default.get_hillshade_path(self.dem_path)
                    else:
                        image = default.get_hillshade(dem_arr=dem, resolution_x=pixel_size[0],
                                                      resolution_y=pixel_size[1])
                elif vis_method.lower() == "multiple directions hillshade":
                    if self.save_individual_vis_fn:  # save individually visualisation function output as raster (tif)
                        default.save_multi_hillshade(self.dem_path)
                        image_path = default.get_multi_hillshade_path(self.dem_path)
                    else:
                        image = default.get_multi_hillshade(dem_arr=dem, resolution_x=pixel_size[0],
                                                            resolution_y=pixel_size[1])
                elif vis_method.lower() == "simple local relief model":
                    if self.save_individual_vis_fn:  # save individually visualisation function output as raster (tif)
                        default.save_slrm(self.dem_path)
                        image_path = default.get_slrm_path(self.dem_path)
                    else:
                        image = default.get_slrm(dem_arr=dem)
                elif vis_method.lower() == "sky-view factor":
                    if self.save_individual_vis_fn:  # save individually visualisation function output as raster (tif)
                        default.save_sky_view_factor(self.dem_path, save_svf=True, save_asvf=False, save_opns=False)
                        image_path = default.get_svf_path(self.dem_path)
                    else:
                        image = default.get_sky_view_factor(dem_arr=dem, resolution=pixel_size[0],
                                                            compute_svf=True,
                                                            compute_asvf=False, compute_opns=False)["svf"]
                elif vis_method.lower() == "anisotropic sky-view factor":
                    if self.save_individual_vis_fn:  # save individually visualisation function output as raster (tif)
                        default.save_sky_view_factor(self.dem_path, save_svf=False, save_asvf=True, save_opns=False)
                        image_path = default.get_asvf_path(self.dem_path)
                    else:
                        image = default.get_sky_view_factor(dem_arr=dem, resolution=pixel_size[0],
                                                            compute_svf=False,
                                                            compute_asvf=True, compute_opns=False)["asvf"]
                elif vis_method.lower() == "openness - positive":
                    if self.save_individual_vis_fn:  # save individually visualisation function output as raster (tif)
                        default.save_sky_view_factor(self.dem_path, save_svf=False, save_asvf=False, save_opns=True)
                        image_path = default.get_opns_path(self.dem_path)
                    else:
                        image = default.get_sky_view_factor(dem_arr=dem, resolution=pixel_size[0],
                                                            compute_svf=False,
                                                            compute_asvf=False, compute_opns=True)["opns"]
                elif vis_method.lower() == "openness - negative":
                    if self.save_individual_vis_fn:  # save individually visualisation function output as raster (tif)
                        default.save_neg_opns(self.dem_path)
                        image_path = default.get_neg_opns_path(self.dem_path)
                    else:
                        image = default.get_neg_opns(dem_arr=dem, resolution=pixel_size[0])
                elif vis_method.lower() == "sky illumination":
                    if self.save_individual_vis_fn:  # save individually visualisation function output as raster (tif)
                        default.save_sky_illumination(self.dem_path)
                        image_path = default.get_sky_illumination_path(self.dem_path)
                    else:
                        image = default.get_sky_illumination(dem_arr=dem, resolution=pixel_size[0])
                elif vis_method.lower() == "local dominance":
                    if self.save_individual_vis_fn:  # save individually visualisation function output as raster (tif)
                        default.save_local_dominance(self.dem_path)
                        image_path = default.get_local_dominance_path(self.dem_path)
                    else:
                        image = default.get_local_dominance(dem_arr=dem)

                layer = rvt.blend.BlenderLayer(vis_method=vis_method, normalization=normalization, minimum=minimum,
                                               maximum=maximum, blend_mode=blend_mode, opacity=opacity, image=image,
                                               image_path=image_path)
            else:
                layer = rvt.blend.BlenderLayer(vis_method=None)

            layers.add_layer(layer)

        layers.layers = layers.layers[-5:]  # ASK ESRI why is there more elements than 5
        blender = layers.render_all_images()  # calculate blender
        pixelBlocks['output_pixels'] = blender.astype(props['pixelType'], copy=False)
        return pixelBlocks


    def prepare(self, default_val_file_path="settings/default_settings.txt", save_individual_vis_fn="Yes",
                layer1_vis="Sky-View Factor", layer1_norm="Value", layer1_min=0.7, layer1_max=1.,
                layer1_blending_mode="Multiply", layer1_opacity=25., layer2_vis="Openness - Positive",
                layer2_norm="Value", layer2_min=68., layer2_max=93., layer2_blending_mode="Overlay",
                layer2_opacity=50., layer3_vis="Slope gradient", layer3_norm="Value", layer3_min=0.,
                layer3_max=50., layer3_blending_mode="Luminosity", layer3_opacity=50., layer4_vis="Hillshade",
                layer4_norm="Value", layer4_min=0., layer4_max=1., layer4_blending_mode="Normal",
                layer4_opacity=100., layer5_vis="None", layer5_norm="Value", layer5_min=0, layer5_max=0,
                layer5_blending_mode="Normal", layer5_opacity=0):

        self.default_val_file_path = default_val_file_path

        if save_individual_vis_fn == "Yes":
            self.save_individual_vis_fn = True
        elif save_individual_vis_fn == "No":
            self.save_individual_vis_fn = False

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
