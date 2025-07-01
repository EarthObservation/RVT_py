from tiled_multiprocess import run_main

if __name__ == "__main__":
    # List of INPUT DATASETS (have to be in this format)
    list_tifs = [
        #r"c:\test_data\SLRM\BAR_dem_05m_SLRM_R20.tif",
        # r"c:\test_data\testmx\2022MxCP_C_dem_05m.tif",
        # r"c:\test_data\testmx_s\2022MxCP_C_dem_05m_s.tif",
        # r"c:\test_data\2022MxCP_1GB\2022MxCP_E_dem_05m.tif",  # Large 1GB sample (66 tiles)
        # r"c:\test_data_tiled\4000\Chactun_izrez_dem_05m.tif", # 4000 px sample
        r"c:\test_data_tiled\6000\ISA-88_Gorteen_dem_05m_clipped.tif" # 6000 px sample
    ]

    # Select visualizations
    # 'slp', 'hs', 'slrm', 'svf', 'opns', 'neg_opns', 'ld', 'mstp'
    vis_types = [
        "slrm"
    ]
    # "vat_combined" / "e2MSTP" / "e4MSTP" / "crim" / "rrim"
    blend_types = [
        # "vat_combined"
    ]

    # save_vrt (if True saves VRT, if False saves as TIF mosaic), only used when tiling
    run_main(list_tifs, vis_types, blend_types, save_float=False, save_vrt=True)
