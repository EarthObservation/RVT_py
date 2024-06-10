from tiled_multiprocess import run_main

if __name__ == "__main__":
    # List of INPUT DATASETS (have to be in this format)
    list_tifs = [
        r"c:\Users\ncoz\GitHub\erc_potencial\test_data\test_small.tif"
    ]

    # Select visualizations
    # 'slp', 'hs', 'slrm', 'svf', 'opns', 'neg_opns', 'ld', 'mstp'
    vis_types = [
        'slrm'
    ]
    # "vat_combined_8bit" / "e2MSTP" / "e4MSTP" / "crim" / "rrim"
    blend_types = [
        "vat_combined_8bit"
    ]

    run_main(list_tifs, vis_types, blend_types)
