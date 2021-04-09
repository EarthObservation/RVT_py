# rvt.default quick TEST

import rvt.default

if __name__ == "__main__":
    # class with default values for all visualization functions
    default = rvt.default.DefaultValues()
    default.save_default_to_file()
    default.read_default_from_file(r"settings\default_settings.json")
    dem_path = r"test_data\TM1_564_146.tif"
    default.overwrite = 1
    default.fill_no_data = 0
    default.fill_method = "kd_tree"
    default.keep_original_no_data = 0

    # save specific visualization function
    default.save_hillshade(dem_path, save_float=True, save_8bit=True)
    default.save_slope(dem_path, save_float=True, save_8bit=True)
    default.save_multi_hillshade(dem_path, save_float=True, save_8bit=True)
    default.save_slrm(dem_path, save_float=True, save_8bit=True)
    default.save_sky_view_factor(dem_path, save_svf=True, save_asvf=True, save_opns=True,
                                 save_float=True, save_8bit=True)
    default.save_neg_opns(dem_path, save_float=True, save_8bit=True)
    default.save_sky_illumination(dem_path, save_float=True, save_8bit=True)
    default.save_local_dominance(dem_path, save_float=True, save_8bit=True)
    default.save_msrm(dem_path, save_float=True, save_8bit=True)

    # save visualization functions, where default.compute = 1
    # for example if default.hs_compute = 1, it will calculate and save hillshade, it also generates log file
    default.save_visualizations(dem_path)
