import rvt.default
import rvt.vis
import rvt.blend_func
import numpy as np


# Prepare
dem_path = r"test_data\TM1_564_146.tif"

dem_dict = rvt.default.get_raster_arr(dem_path)
dem_arr = dem_dict["array"]
dem_res = dem_dict["resolution"][0]
dem_nodata = dem_dict["no_data"]

default = rvt.default.DefaultValues()

default.mstp_local_scale = (3, 21, 2)
default.mstp_meso_scale = (23, 203, 18)
default.mstp_broad_scale = (223, 2023, 180)
default.mstp_lightness = 1.2

# Run visualization function
rendered_image = rvt.vis.mstp(dem=dem_arr, no_data=dem_nodata)

# # Save result (float32)
# output_mstp_path = dem_path.replace(".tif", "_mstp_float.tif")
# rvt.default.save_raster(src_raster_path=dem_path, out_raster_arr=rendered_image, out_raster_path=output_mstp_path,
#                         e_type=6)

# Save 8-bit
output_mstp_path_8bit = dem_path.replace(".tif", "_mstp_8bit2.tif")
# # Can skip this as data is already between 0 and 1
# norm_arr = rvt.blend_func.normalize_image(
#                 visualization="mstp",
#                 image=rendered_image,
#                 min_norm=0,
#                 max_norm=1,
#                 normalization="value"
#             )
# Convert to 8 bit
image_8bit = rvt.vis.byte_scale(data=rendered_image, no_data=np.nan, c_min=0, c_max=1)
rvt.default.save_raster(src_raster_path=dem_path, out_raster_path=output_mstp_path_8bit, out_raster_arr=image_8bit,
                        no_data=np.nan, e_type=1)  # e_type has to be 1 because visualization is 8-bit (0-255)

# # Try saving as 8bit using the default method
# default.mstp_save_float = False
# default.mstp_save_8bit = True
# default.save_mstp(dem_path)
