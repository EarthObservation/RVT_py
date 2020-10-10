# Relief Visualization Toolbox in Python

Relief Visualization Toolbox was produced to help scientist visualize raster elevation model datasets. We have narrowed down the selection to include techniques that have proven to be effective for identification of small scale features. Default settings therefore assume working with high resolution digital elevation models, derived from airborne laser scanning missions (lidar).

Despite this, techniques are also used for different other purposes. Sky-view factor, for example, can be efficiently used in numerous studies where digital elevation model visualizations and automatic feature extraction techniques are indispensable, e.g. in geography, archeology,  geomorphology, cartography, hydrology, glaciology, forestry and disaster management. It can be used even in engineering applications, such as, predicting the availability of the GPS signal in urban areas.

Methods currently implemented are:

*   hillshading,
*   hillshading from multiple directions,
*   slope gradient,
*   simple local relief model,
*   sky illumination,
*   sky-view factor (as developed by our team),
*   anisotropic sky-view factor,
*   positive and negative openness,
*   local dominance.

For a more detailed description see references given at each method in the manual and a comparative paper describing them (e.g. Kokalj and Hesse 2017, see below).

RVT python library called rvt contains 3 modules: vis (rvt.vis), blend (rvt.blend) and default (rvt.default). Modules contains:
* vis       -   visualization functions (mentioned above), for computing visualizations;
* blend     -   blender (mixer), for blending visualizations;
* default   -   default values, class for defining default parameters with methods to compute and save visualization functions using set parameters.

For every visualization function directory also contains Python Esri raster functions for ArcGIS Pro (rvt_esri_*.py).

## References

When using the tools, please cite:

*   Kokalj, Ž., Somrak, M. 2019. Why Not a Single Image? Combining Visualizations to Facilitate Fieldwork and On-Screen Mapping. Remote Sensing 11(7): 747.
*   Zakšek, K., Oštir, K., Kokalj, Ž. 2011. Sky-View Factor as a Relief Visualization Technique. Remote Sensing 3: 398-415.
*   Kokalj, Ž., Zakšek, K., Oštir, K. 2011. Application of Sky-View Factor for the Visualization of Historic Landscape Features in Lidar-Derived Relief Models. Antiquity 85, 327: 263-273.

## Installation

Copy or clone the files to your environment.

### Python libraries
We suggest using an anaconda enviroment and python version 3.5 or higher.
You'll need libraries (could also work with other versions):
*   numpy 1.19.2
*   scipy 1.5.2
*   rasterio 1.1.7

### ArcGIS Pro
If you would like to use visualization functions in ArcGIS Pro you have to click Analysis->Raster functions->"three parallel lines"->Open Python Raster Function. Then select Python Module (rvt_esri_*.py) and Class Name (it is only one).

For rvt_esri_blender.py you will need to install rasterio into Python ArcGIS Pro conda environment. To do that open ArcGIS Pro click Python then click Manage Environments and Clone the default environment. After you clone default environment set new env to active and click ok. Then click Add Packages search for rasterio and install it. In case you are having problems try older version of rasterio.

## Usage

### Module vis
For reading raster data (dem) from files (tif) we suggest to use rasterio.
```python
import rasterio as rio

input_dem_path = "test_data/TM1_564_146.tif"  # path to your DEM
input_dem_dataset = rio.open(input_dem_path)  # open with rasterio
t = input_dem_dataset.transform
x_res = t[0]  # to get x resolution
y_res = -t[4] # to get y resolution
input_dem_arr = input_dem_dataset.read()[0]  # to get 2D numpy array of DEM
input_dem_dataset.close()
```

Import rvt.vis module:
```python
import rvt.vis
```

*   **Slope Gradient** <br>
Parameters <br>
---------- <br>
ve_factor : vertical exaggeration factor (must be greater than 0) <br>
output_units : percent, degree, radians <br> <br>
Returns <br>
------- <br>
{"slope": slope_out, "aspect": aspect_out} : dictionaries with 2D numpy arrays <br>
```python
dict_slp_asp = rvt.vis.slope_aspect(dem=input_dem_arr, resolution_x=x_res, resolution_y=y_res, ve_factor=ve_factor,
                                    output_units=output_units)
slope_arr = dict_slp_asp["slope"]
aspect_arr = dict_slo_asp["aspect"]
```

*   **Hillshade** <br>
Parameters <br>
---------- <br>
sun_azimuth : solar azimuth angle (clockwise from North) in degrees <br>
sun_elevation : solar vertical angle (above the horizon) in degrees <br> <br>
Returns <br>
------- <br>
hillshade_out : result 2D numpy array <br>
```python
hillshade_arr = rvt.vis.hillshade(dem=input_dem_arr, resolution_x=x_res, resolution_y=y_res, sun_azimuth=sun_azimuth,
                                  sun_elevation=sun_elevation)
```

*   **Hillshade from multiple directions** <br>
Parameters <br>
---------- <br>
nr_directions : number of solar azimuth angles (clockwise from North) <br>
sun_elevation : solar vertical angle (above the horizon) in degrees <br> <br>
Returns <br>
------- <br>
multi_hillshade_out : numpy array containing numpy_arrays of hillshades in different directions <br>
```python
multi_hillshade_arr = rvt.vis.multi_hillshade(dem=input_dem_arr, resolution_x=x_res, resolution_y=y_res,
                                              nr_directions=nr_directions, sun_elevation=sun_elevation)
```

*   **Simple Local Relief Model** <br>
Parameters <br>
---------- <br>
radius_cell : Radius for trend assessment [pixels] <br> <br>
Returns <br>
------- <br>
slrm_out : slrm 2D numpy array <br>
```python
slrm_arr = rvt.vis.slrm(dem=input_dem_arr, radious_cell=radious_cell)
```

*   **Sky Illumination** <br>
Parameters <br>
---------- <br>
sky_model : sky model [overcast, uniform] <br>
sampling_points : number of sampling points <br>
shadow_dist : max shadow modeling distance [pixels] <br>
shadow_az : shadow azimuth <br>
shadow_el : shadow elevation <br> <br>
Returns <br>
------- <br>
sky_illum_out : 2D numpy result array <br>
```python
sky_illumination_arr = rvt.vis.sky_illumination(dem=input_dem_arr, resolution=x_res, sky_model=sky_model,
                                                sampling_points=sampling_points, shadow_dist=shadow_dist,
                                                shadow_az=shadow_az, shadow_el=shadow_el)
```

*   **Sky-View Factor / Anisotropic Sky-View Factor / Positive and Negative Openness** <br>
Sky-view factor, Antisotropic Sky-view factor and Openness are all calculated in the same function (less computing if you have to compute all of them). Negative openness is openness where (-1)*dem. <br><br>
Parameters <br>
compute_svf : compute SVF (True) or not (False) <br>
compute_asvf : compute anisotropic SVF (True) or not (False) <br>
compute_opns : compute OPENNESS (True) or not (False) <br>
resolution : pixel resolution <br>
svf_n_dir : number of directions <br>
svf_r_max : maximal search radius in pixels <br>
svf_noise : the level of noise remove (0-don't remove, 1-low, 2-med, 3-high) <br>
asvf_level : level of anisotropy, 1-low, 2-high, <br>
a_min_weight : weight to consider anisotropy (0 - isotropic, 1 - no illumination from the direction
           opposite the main direction) <br> <br>
Returns <br>
------- <br>
{"svf": svf_out, "asvf": asvf_out, "opns": opns_out} : dictionary <br>
svf_out, skyview factor : 2D numpy vector of skyview factor. <br>
asvf_out, anisotropic skyview factor : 2D numpy vector of anisotropic skyview factor. <br>
opns_out, openness : 2D numpy openness (elevation angle of horizon) <br>
```python
# Compute Svf, Asvf, Pos Opns simultaneously
dict_svf_asvf_opns = rvt.vis.sky_view_factor(dem=input_dem_arr, resolution=x_res, compute_svf=True,
                                             compute_asvf=True, compute_opns=True, svf_n_dir=svf_n_dir,
                                             svf_r_max=svf_r_max, svf_noise=svf_noise, asvf_level=asvf_level,
                                             a_min_weight=a_min_weight)
svf_arr = dict_svf_asvf_opns["svf"]  # sky-view factor
asvf_arr = dict_svf_asvf_opns["asvf"]  # antiostropic sky-view factor
opns_arr = dict_svf_asvf_opns["opns"]  # positive openness

# Compute only Svf
svf_arr = rvt.vis.sky_view_factor(dem=input_dem_arr, resolution=x_res, compute_svf=True,
                                  compute_asvf=False, compute_opns=False, svf_n_dir=svf_n_dir,
                                  svf_r_max=svf_r_max, svf_noise=svf_noise)["svf"]  # be aware ["svf"], because result is dict

# Compute only Asvf
asvf_arr = rvt.vis.sky_view_factor(dem=input_dem_arr, resolution=x_res, compute_svf=False,
                                   compute_asvf=True, compute_opns=False, svf_n_dir=svf_n_dir,
                                   svf_r_max=svf_r_max, svf_noise=svf_noise, asvf_level=asvf_level,
                                   a_min_weight=a_min_weight)["asvf"]  # be aware ["asvf"], because result is dict

# Compute only Positive Openness
opns_arr = rvt.vis.sky_view_factor(dem=input_dem_arr, resolution=x_res, compute_svf=False,
                                   compute_asvf=False, compute_opns=True, svf_n_dir=svf_n_dir,
                                   svf_r_max=svf_r_max, svf_noise=svf_noise)["opns"]  # be aware ["opns"], because result is dict

# Compute Negative Openness, be aware (-1)*dem
neg_opns_arr = rvt.vis.sky_view_factor(dem=(-1)*input_dem_arr, resolution=x_res, compute_svf=False,
                                   compute_asvf=False, compute_opns=True, svf_n_dir=svf_n_dir,
                                   svf_r_max=svf_r_max, svf_noise=svf_noise)["opns"]  # be aware ["opns"], because result is dict
```

*   **Local Dominance** <br>
Parameters <br>
---------- <br>
min_rad : minimum radial distance (in pixels) at which the algorithm starts with visualization computation <br>
max_rad : maximum radial distance (in pixels) at which the algorithm ends with visualization computation <br>
rad_inc : radial distance steps in pixels <br> 
angular_res : angular step for determination of number of angular directions <br>
observer_height : height at which we observe the terrain <br><br>
Returns <br>
------- <br>
local_dom_out - 2D numpy array of local dominance <br>
```python
local_dominance_arr = rvt.vis.local_dominance(dem=input_dem_arr, min_rad=min_rad, max_rad=max_rad, rad_inc=rad_inc,
                                              angular_res=angular_res, observer_height=observer_height)
```

To save computed array we also recommend to use rasterio.
```python
output_path = "test_data/TM1_564_146_computed.tif"  # output tif path (new file)
profile = input_dem_dataset.profile  # use profile from input DEM dataset
profile.update(dtype='float32')  # if multihilshade add count=nr_bands
output_computed_dataset = rio.open(output_path, "w", **profile)
# if multihillshade it is already 3D array .write(multihillshade)
output_computed_dataset.write(np.array([computed_arr]))  # computed_arr is array with computed visualization, ex: slope_arr
output_computed_dataset.close()
```

### Module default
Default module contains class DefaultValues() where we can store our visualization functions parameters, this class also has methods for saving and computing visualization functions with that parameters. To import it we use:
```python
import rvt.default
```
We have to create class instance first.
```python
default = rvt.default.DefaultValues()
```
When we create instance parameter values are already populated with default values. We can store them in file, change them in file and read them back from file. We can also change them in program.
```python
default.save_default_to_file(file_path=r"settings\default_settings.txt")  # save parameters in file
# we can make a copy of our own parameters values and load them to default
default.read_default_from_file(file_path="settings\default_settings.txt")  # read parameters from file
# example to change specific parameter programmatically
default.hs_sun_el = 45  # change hillshade sun elevation
```
After we have set our parameters, we can compute visualization functions and get results in array or we can save them in file. If we save them result will be saved where input raster location is. You don't need to use get methods if you only need saved results.
```python
# Slope
dict_slp_asp = default.get_slope(dem_arr=dem_arr, resolution_x=x_res, resolution_y=y_res)  # you have to read dem and input dem_arr, it returns dict, same as rvt.vis.slope_aspect()
default.save_slope(dem_path=dem_path)  # reads dem, computes slope and saves slope in directory where input raster is (dem_path)

# Hillshade
hillshade_arr = default.get_hillshade(dem_arr=dem_arr, resolution_x=x_res, resolution_y=y_res)
default.save_hillshade(dem_path=dem_path)

# Multi hillshade
multi_hillshade_arr = default.get_multi_hillshade(dem_arr=dem_arr, resolution_x=x_res, resolution_y=y_res)
default.save_multi_hillshade(dem_path=dem_path)

# Simple local relief model
slrm_arr = get_slrm(dem_arr=dem_arr)
default.save_slrm(dem_path=dem_path)

# Svf, ASvf, Pos Opns, You have to define with compute and save (True/False) which ones you need
dict_svf_asvf_opns = default.get_sky_view_factor(dem_arr=dem_arr, resolution=res, compute_svf=True, compute_asvf=True, compute_opns=True)  # it returns dict same as rvt.vis.sky_view_factor()
default.save_sky_view_factor(dem_path=dem_path, save_svf=True, save_asvf=True, save_opns=True)  # define (True/False) which ones to save

# Negative Openness
neg_opns_arr = default.get_neg_opns(dem_arr=dem_arr, resolution=res)
default.save_neg_opns(dem_path=dem_path)

# Sky illumination
sky_illumination_arr = default.get_sky_illumination(dem_arr=dem_arr, resolution=res)
default.save_sky_illumination(dem_path=dem_path)

# Local dominance
local_dominance_arr = default.get_local_dominance(dem_arr=dem_arr)
default.save_local_dominance(dem_path=dem_path)

# To save visualizations with one line of code
default.save_visualizations(dem_path=dem_path, sav_slope=True, sav_hillshade=True, sav_mulit_hillshade=True, sav_slrm=True,
                            sav_svf=True, sav_asvf=True, sav_opns=True, sav_neg_opns=True, sav_sky_illumination=True,
                            sav_local_dominance=True)  # if you want to save specific vis set it to True else to False
```

### Module blend
Blend is module for blending different visualizations together. To import it we use:
```python
import rvt.blend
```
We could use manual blending (we compute visualisations) or we could use automatic blending from file (automatically computed visualisations with rvt.vis and stored in location where dem is).
Manual blending depends on module rvt.default (look section Module default) and rvt.vis . To start blending we first need to create instance of class BlenderLayers() which contains list of layers (BlenderLayer()). Single layer is defined in BlenderLayer() class.
#### Manual blending
```python
layers_manual = rvt.blend.BlenderLayers()  # create class which will hold layers
# you have two options to add layer:
# option 1, create with method
layers_manual.create_layer(vis_method="Sky-View Factor", normalization="value", minimum=0.7, maximum=1,
                           blend_mode="multiply", opacity=25,
                           image=svf_arr)  # automatically creates BlenderLayer() and appends it to BlenderLayers()
# option 2, create class BlenderLayer instance and then add with method
layer1 = rvt.blend.BlenderLayer(vis_method="Sky-View Factor", normalization="value", minimum=0.7, maximum=1,
                                blend_mode="multiply", opacity=25,
                                image=svf_arr)
layers_manual.add_layer(layer1)
```
You can add as many layers as you need. When adding / creating layers you can define image or image_path parameter. If you define image_path (you have to save image first) and not image then blending will work faster because it will not hold all images (from all layers) in memory. It will read them simultaneously.
```python
layers_manual.create_layer(vis_method="Sky-View Factor", normalization="value", minimum=0.7, maximum=1,
                           blend_mode="multiply", opacity=25,
                           image_path=svf_path)  # image_path instead of image
```
After you added all the layers you would like to blend. You call method render_all(). If you define method parameters input_dem_path (needed for profile info) and output_blend_path. Result will be saved in output_blend_path else it will only return result raster array.
```python
render_arr = layers_manual.render_all_images(input_dem_path, output_blend_path)  # to save rendered array in output_blend_path
render_arr = layers_manual.render_all_images()  # to only get result render array (render_arr)
```

#### Automatic blending
Automatic blending depends on rvt.default, so you have to import rvt.default.
```python
import rvt.blend
import rvt.default
```
Automatic blending is blending from file. All the visualisations we need for blending will be computed and stored in the directory where input raster is. To create example file where we can later change parameters we call function create_blender_file_example.
```python
blender_file = rvt.blend.create_blender_file_example(file_path=r"settings\blender_file_example.txt")
```
To blend from file we also need function parameters values which we define in class rvt.default.DefaultValues() (look Module default).
```python
default = rvt.default.DefaultValues()
```
To blend from file we create BlenderLayers() class, call method build_blender_layers_from_file and then render_all_images.
```python
layers_auto = rvt.blend.BlenderLayers()
layers_auto.build_blender_layers_from_file(dem_path=input_dem_path, file_path=blender_file, default=default)  # we can make our own blender_file (change example), we can change default
render_arr = layers_auto.render_all_images(dem_path=input_dem_path, save_render_path=output_blend_path)  # if we dont input dem_path and save_rander_path result won't be stored as raster (tif)
```

## Examples

A sample dataset for trying RVT python is available here ("TM1_564_146.tif"):  
https://www.dropbox.com/sh/p7ia8fk6mywa8y3/AABWuw4wFUvULU7SeNXyWhjka?dl=0

Download it, save it in test_data directory and try the visualisations.

Examples how to use are in files:
```python
test_vis.py

test_blend.py

test_default.py

test_custom_color_scheme.py
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please report any bugs and suggestions for improvements.

## Acknowledgment

Development of RVT was part financed by the European Commission's Culture Programme through the ArchaeoLandscapes Europe project and by the Slovenian Research Agency core funding No. P2-0406, and by research projects No. J6-7085 and No. J6-9395.

## License
This project is licensed under the terms of the [Apache License](LICENSE).

