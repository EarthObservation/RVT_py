# rvt.default quick TEST

import rvt.default

# class with default values for all visualization functions
default = rvt.default.DefaultValues()
default.read_default_from_file("D:\RVT_py\RVT_py\settings\default_settings.txt")

# save specific visualization function
default.save_hillshade(r"D:\RVT_py\test\test.tif")
default.save_slope(r"D:\RVT_py\test\test.tif")
default.save_multi_hillshade(r"D:\RVT_py\test\test.tif")
default.save_slrm(r"D:\RVT_py\test\test.tif")
default.save_sky_view_factor(r"D:\RVT_py\test\test.tif", save_svf=True, save_asvf=True, save_opns=True)
default.save_neg_opns(r"D:\RVT_py\test\test.tif")
#default.save_sky_illumination(r"D:\RVT_py\test\test.tif")
default.save_local_dominance(r"D:\RVT_py\test\test.tif")

# save all visualization functions
default.save_visualizations(r"D:\RVT_py\test\test.tif", sav_sky_illumination=False)
