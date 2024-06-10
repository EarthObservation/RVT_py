from tiled_multiprocess import run_main

if __name__ == "__main__":
    # List of INPUT DATASETS (have to be in this format)
    list_tifs = [
        r"c:\Users\ncoz\GitHub\erc_potencial\test_data\test_small.tif"
    ]

    # Select visualizations
    # 'slope_aspect', 'hillshade', 'multi_hillshade', 'slrm', 'sky_view_factor', 'openness', 'neg_openness',
    # 'local_dominance', 'sky_illumination', 'shadow_horizon', 'msrm', 'mstp'
    vis_types = [
        "slrm"
    ]
    # "vat_combined_8bit" / "e2MSTP" / "e4MSTP" / "CRIM" / "RRIM"
    blend_types = [
        "vat_combined_8bit"
    ]

    run_main(list_tifs, vis_types, blend_types)
