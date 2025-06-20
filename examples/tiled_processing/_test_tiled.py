from tiled_multiprocess import run_main

if __name__ == "__main__":
    # List of INPUT DATASETS (have to be in this format)
    list_tifs = [
        r"c:\path-to-file\test-file.tif",
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
    run_main(list_tifs, vis_types, blend_types, save_float=False, save_vrt=False)
