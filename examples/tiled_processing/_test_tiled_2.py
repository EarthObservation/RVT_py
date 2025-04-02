from tiled_multiprocess import run_main

if __name__ == "__main__":
    # List of INPUT DATASETS (have to be in this format)
    list_tifs = [
        r"c:\test_data\SLRM\BAR_dem_05m_SLRM_R20.tif",
        # r"c:\test_data\testmx\2022MxCP_C_dem_05m.tif",
        # r"c:\test_data\testmx_s\2022MxCP_C_dem_05m_s.tif"
    ]

    # Select visualizations
    # 'slp', 'hs', 'slrm', 'svf', 'opns', 'neg_opns', 'ld', 'mstp'
    vis_types = [
        'slrm'
    ]
    # "vat_combined" / "e2MSTP" / "e4MSTP" / "crim" / "rrim"
    blend_types = [
        "new_blend"
    ]

    run_main(list_tifs, vis_types, blend_types, save_float=False)
