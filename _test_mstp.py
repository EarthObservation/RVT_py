import rvt.default
import rvt.vis


dem_path = r"test_data\TM1_564_146.tif"

dem_dict = rvt.default.get_raster_arr(dem_path)
dem_arr = dem_dict["array"]
dem_res = dem_dict["resolution"][0]
dem_nodata = dem_dict["no_data"]

output_e3mstp_path = "{}_{}.tif".format(dem_path.rstrip(".tif"), "mstp")
default = rvt.default.DefaultValues()
default.mstp_local_scale = (1, 5, 1)
default.mstp_meso_scale = (5, 50, 5)
default.mstp_broad_scale = (50, 500, 50)
default.mstp_lightness = 1.2
default.slrm_rad_cell = 10

rendered_image = rvt.vis.mstp(dem=dem_arr, no_data=dem_nodata)

rvt.default.save_raster(src_raster_path=dem_path, out_raster_arr=rendered_image, out_raster_path=output_e3mstp_path,
                        e_type=6)
