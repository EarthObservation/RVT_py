import rvt.default
import rvt.vis
import numpy as np


dem_path = r"test_data\TM1_564_146.tif"
raster_dict = rvt.default.get_raster_arr(dem_path)
dem_arr = raster_dict["array"]
dem_res = raster_dict["resolution"][0]
dem_no_data = raster_dict["no_data"]

msrm_arr = rvt.vis.msrm(dem=dem_arr, resolution=dem_res, feature_min=5, feature_max=20, scale_factor=2,
                        no_data=dem_no_data)

rvt.default.save_raster(src_raster_path=dem_path, out_raster_arr=msrm_arr,
                        out_raster_path=r"test_data\TM1_564_146_msrm_test.tif", no_data=np.nan, e_type=6)