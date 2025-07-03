import rvt.default
import rvt.blend
import rvt.vis


# Test enhanced version 4 multi-scale topographic position (e3MSTP)
dem_path = r"c:\test_data\test_small\test_small.tif"
dem_dict = rvt.default.get_raster_arr(dem_path)
dem_arr = dem_dict["array"]
dem_res = dem_dict["resolution"][0]

# Output path
output_e4mstp_path = "{}_{}.tif".format(dem_path.rstrip(".tif"), "e4mstp_TEST")


default = rvt.default.DefaultValues()
default.read_default_from_file(r"examples\tiled_processing\default_1.json")
# default.mstp_local_scale = (1, 5, 1)
# default.mstp_meso_scale = (5, 50, 5)
# default.mstp_broad_scale = (50, 500, 50)
# default.mstp_lightness = 0.9
# default.slrm_rad_cell = 10

rendered_image = rvt.blend.e3mstp(dem=dem_arr, resolution=dem_res, default=default)

rvt.default.save_raster(
    src_raster_path=dem_path,
    out_raster_arr=rendered_image,
    out_raster_path=output_e4mstp_path,
    e_type=6
)