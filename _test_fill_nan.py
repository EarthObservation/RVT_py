import rvt.visualizations
import rvt.default
import numpy as np
import time

rast_path = r"D:\RVT_py\test\manhattan\Manhattan_DSM_1m_max.tif"
dict_rast = rvt.default.get_raster_arr(rast_path)
arr = dict_rast["array"]
nodata = dict_rast["no_data"]
arr[arr == nodata] = np.nan

# idw
start_time = time.time()
arr_out = rvt.vis.fill_where_nan(arr, "idw_20_2.0")
print("idw = {}s".format(time.time() - start_time))
rast_path_out = r"D:\RVT_py\test\manhattan\Manhattan_DSM_1m_max_tst_idw.tif"
rvt.default.save_raster(src_raster_path=rast_path, out_raster_path=rast_path_out,
                        out_raster_arr=arr_out, no_data=np.nan)

# kd_tree
start_time = time.time()
arr_out = rvt.vis.fill_where_nan(arr, "kd_tree")
print("kd_tree = {}s".format(time.time() - start_time))
rast_path_out = r"D:\RVT_py\test\manhattan\Manhattan_DSM_1m_max_tst_kd_tree.tif"
rvt.default.save_raster(src_raster_path=rast_path, out_raster_path=rast_path_out,
                        out_raster_arr=arr_out, no_data=np.nan)

# nearest_neighbour
start_time = time.time()
arr_out = rvt.vis.fill_where_nan(arr, "nearest_neighbour")
print("nearest_neighbour = {}s".format(time.time() - start_time))
rast_path_out = r"D:\RVT_py\test\manhattan\Manhattan_DSM_1m_max_tst_nearest_neighbour.tif"
rvt.default.save_raster(src_raster_path=rast_path, out_raster_path=rast_path_out,
                        out_raster_arr=arr_out, no_data=np.nan)
