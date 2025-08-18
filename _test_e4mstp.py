"""
Test enhanced version 4 multi-scale topographic position (e4MSTP)

Copyright:
    2025 Research Centre of the Slovenian Academy of Sciences and Arts
    2025 University of Ljubljana, Faculty of Civil and Geodetic Engineering
"""
import rvt.default
import rvt.blend
import rvt.vis


# Prepare input data for test
dem_path = r"c:\test_data\test_small\test_small.tif"
dem_dict = rvt.default.get_raster_arr(dem_path)
dem_arr = dem_dict["array"]
dem_res = dem_dict["resolution"][0]

output_e4mstp_path = "{}_{}.tif".format(dem_path.rstrip(".tif"), "e4mstp_TEST")

default = rvt.default.DefaultValues()
default.read_default_from_file(r"examples\tiled_processing\default_1.json")


rendered_image = rvt.blend.e4mstp(dem=dem_arr, resolution=dem_res, default=default)

# Run test
rvt.default.save_raster(
    src_raster_path=dem_path,
    out_raster_arr=rendered_image,
    out_raster_path=output_e4mstp_path,
    e_type=6
)
