import rvt.default
import rvt.blend_func
import rvt.blend
import rvt.vis
import numpy as np


# Test enhanced version 3 multi-scale topographic position (e3MSTP)
dem_path = r"D:\RVT_py\test\E3MSTP_TEST\Aquileia_DBM_05m.tif" #r"D:\RVT_py\test\E3MSTP_TEST\Skolj_dem_05m.tif"
dem_dict = rvt.default.get_raster_arr(dem_path)
dem_arr = dem_dict["array"]
dem_res = dem_dict["resolution"][0]

output_e3mstp_path = "{}_{}.tif".format(dem_path.rstrip(".tif"), "e3mstp")
rendered_image = rvt.blend.e3mstp(dem=dem_arr, resolution=dem_res)
rvt.default.save_raster(src_raster_path=dem_path, out_raster_arr=rendered_image, out_raster_path=output_e3mstp_path,
                        e_type=6)