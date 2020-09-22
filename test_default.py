# rvt.default quick TEST

import rvt.default


default = rvt.default.DefaultValues()
default.read_default_from_file("D:\RVT_py\RVT_py\settings\default_settings.txt")
default.save_hillshade(r"D:\RVT_py\test\test.tif")
default.save_slope(r"D:\RVT_py\test\test.tif")