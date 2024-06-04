from pathlib import Path

import rasterio

import grid_tools as gt
from tiled_multiprocess import tiled_blending


# List of INPUT DATASETS (have to be in this format)
list_tifs = [
    Path(r"c:\Users\ncoz\GitHub\erc_potencial\test_data\test_small.tif")
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


for in_file in list_tifs:

    # (1) Instead of folder, the tif/vrt is given
    input_vrt = in_file.as_posix()
    ds = in_file.parent

    print("Start --- " + ds.name)

    # (2) Check if larger than 4000 px
    with rasterio.open(input_vrt) as src:
        if src.shape > (6000, 6000):
            # If at least 1 dim is larger than 6000 pixels
            tile_size = 4000
        elif src.shape > (4000, 4000):
            tile_size = 2000
        else:
            tile_size = None

    if tile_size:
        # (3) To filter we need polygon covering valid data
        vdm_file = list(ds.glob("*_validDataMask*"))

        if vdm_file:
            print("* validDataMask exists; Removing and creating new...")
            Path(vdm_file[0]).unlink()
        else:
            # If it doesn't exist, try creating from VRT
            print("* validDataMask doesn't exists; Creating...")

        valid_data_outline = gt.poly_from_valid(input_vrt, save_gpkg=False)

        # (4) Create reference grid, filter it and save it to disk
        print("Create REF GRID")

        refgrid_name = input_vrt[:-4] + "_refgrid.gpkg"
        if Path(refgrid_name).exists():
            Path(refgrid_name).unlink()

        tiles_extents = gt.bounding_grid(input_vrt, tile_size, tag=True)
        tiles_extents = gt.filter_by_outline(
            tiles_extents, valid_data_outline.as_posix(),
            save_gpkg=False, save_path=refgrid_name
        )

        tiles_list = tiles_extents[["minx", "miny", "maxx", "maxy"]].values.tolist()
    else:
        tiles_list = None

    # (4) Run tiled blending (visualizations calculated on the go, stored in memory)
    tiled_blending(
        vis_types=vis_types,
        blend_types=blend_types,
        input_vrt_path=in_file,
        tiles_list=tiles_list
    )
