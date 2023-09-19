# rvt.blend.gray_scale_to_color_ramp quick TEST

import rasterio as rio
from rasterio.enums import ColorInterp
import rvt.visualizations
import rvt.blender
import numpy as np


#  Computed slope for example and saved as custom color scheme
input_dem_path = r"test_data\TM1_564_146.tif"
output_path = r"test_data\TM1_564_146_test_slp_clr.tif"

input_dem_dataset = rio.open(input_dem_path)
t = input_dem_dataset.transform
x_res = t[0]
y_res = -t[4]
input_dem_arr = input_dem_dataset.read()[0]
# compute slope for example
dict_slp_asp = rvt.vis.slope_aspect(
    dem=input_dem_arr,
    resolution_x=x_res,
    resolution_y=y_res,
    vertical_exaggeration_factor=1,
    output_units="degree",
)
slope_arr = dict_slp_asp["slope"]
slope_arr_norm = np.interp(
    slope_arr, (np.nanmin(slope_arr), np.nanmax(slope_arr)), (0, 1)
)  # normalize (0-1)
# more colors: https://matplotlib.org/3.3.2/tutorials/colors/colormaps.html
slope_arr_rgba = rvt.blend.gray_scale_to_color_ramp(
    gray_scale=slope_arr_norm, colormap="BuGn"
)  # compute color scheme

profile = input_dem_dataset.profile
profile.update(dtype="uint8", count=4, photometric="RGBA")
output_slope_dataset = rio.open(output_path, "w", **profile)
output_slope_dataset.colorinterp = [
    ColorInterp.red,
    ColorInterp.green,
    ColorInterp.blue,
    ColorInterp.alpha,
]
output_slope_dataset.write(slope_arr_rgba)
input_dem_dataset.close()
output_slope_dataset.close()
