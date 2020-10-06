# rvt.default quick TEST

import rvt.default

# class with default values for all visualization functions
default = rvt.default.DefaultValues()
default.read_default_from_file("settings\default_settings.txt")
default.save_default_to_file()
# save specific visualization function
default.save_hillshade(r"test_data\TM1_564_146.tif")
default.save_slope(r"test_data\TM1_564_146.tif")
default.save_multi_hillshade(r"test_data\TM1_564_146.tif")
default.save_slrm(r"test_data\TM1_564_146.tif")
default.save_sky_view_factor(r"test_data\TM1_564_146.tif", save_svf=True, save_asvf=True, save_opns=True)
default.save_neg_opns(r"test_data\TM1_564_146.tif")
#default.save_sky_illumination(r"test_data\test.tif")
default.save_local_dominance(r"test_data\TM1_564_146.tif")

# save all visualization functions
default.save_visualizations(r"test_data\TM1_564_146.tif", sav_sky_illumination=False)