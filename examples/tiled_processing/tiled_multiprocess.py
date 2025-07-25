"""
Tiled processing with RVT_py
Created on 6 May 2024
@author: Nejc Čož, ZRC SAZU, Novi trg 2, 1000 Ljubljana, Slovenia
"""

import multiprocessing as mp
import os
import shutil
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import rasterio
from geopandas import GeoDataFrame
from osgeo import gdal
from rasterio.merge import merge
from rasterio.windows import from_bounds
from shapely.geometry import box

import grid_tools as gt
import rvt.blend
import rvt.default
import rvt.vis
from rvt.blend_func import normalize_image

gdal.UseExceptions()


def run_main(list_tifs, vis_types, blend_types, save_float=False, save_vrt=False):
    for in_file in list_tifs:
        in_file = Path(in_file)

        # (1) Instead of folder, the tif/vrt is given
        input_vrt = in_file.as_posix()

        print("Start --- " + in_file.as_posix())

        # (2) Check if larger than 4000 px
        with rasterio.open(input_vrt) as src:
            if src.shape > (6000, 6000):
                # If at least 1 dim is larger than 6000 pixels
                tile_size = 4000
            elif src.shape > (4000, 4000):
                tile_size = 2000
            else:
                tile_size = None
            raster_bounds = src.bounds

        if tile_size:
            # # (3) To filter we need polygon covering valid data
            # valid_data_outline = gt.poly_from_valid(input_vrt, save_gpkg=False)

            # (4) Create reference grid
            refgrid_name = input_vrt[:-4] + "_refgrid.gpkg"
            if Path(refgrid_name).exists():
                Path(refgrid_name).unlink()
            tiles_extents = gt.bounding_grid(input_vrt, tile_size, tag=False)

            # # (5) Filter the grid
            # tiles_extents = gt.filter_by_outline(
            #     tiles_extents,
            #     valid_data_outline,
            #     save_gpkg=False
            # )

            # Crop outer tiles to the extents of raster (IMPORTANT so output size is same as input when tiling!)
            raster_bbox = box(*raster_bounds)
            tiles_extents = GeoDataFrame(geometry=tiles_extents.geometry.intersection(raster_bbox), crs=tiles_extents.crs)

            # Add extents column: Extents are (L, B, R, T).
            tiles_extents["extents"] = tiles_extents.bounds.apply(lambda x: (x.minx, x.miny, x.maxx, x.maxy), axis=1)

            # Extract list from GeoDataFrame
            tiles_list = tiles_extents["extents"].values.tolist()
        else:
            tiles_list = None

        # (4) Run tiled blending (visualizations calculated on the go, stored in memory)
        tiled_blending(
            vis_types=vis_types,
            blend_types=blend_types,
            input_vrt_path=in_file,
            tiles_list=tiles_list,
            save_float=save_float,
            save_vrt=save_vrt
        )


def tiled_blending(vis_types, blend_types, input_vrt_path, tiles_list, save_float, save_vrt):
    t0 = time.time()

    # Prepare paths
    src_tif_path = Path(input_vrt_path)
    ll_path = src_tif_path.parent

    # Run multiprocessing if tiles_list is provided (else single image processing)
    if tiles_list:
        # Determine nr_processes from available CPUs (leave two free)
        nr_processes = os.cpu_count() - 2
        if nr_processes < 1:
            nr_processes = 1
        print('Threads running on:', nr_processes)

        results = []
        # # HERE MULTIPROCESSING STARTS
        #
        # - src_tif_path  (const)
        # - ll_path  (const)
        # - vis_types  (const)
        # - blend_types  (const)
        # - save_float  (const)
        # - one_tile in tiles_list  (variable)
        #
        # ----------------------------------
        input_process_list = [(src_tif_path, ll_path, vis_types, blend_types, i, save_float) for i in tiles_list]
        with mp.Pool(nr_processes) as p:
            realist = [p.apply_async(compute_save_blends, r) for r in input_process_list]
            for i, result in enumerate(realist):
                pool_out = result.get()
                results.append(pool_out)
                print(f"Finished tile {i + 1} of {len(tiles_list)}")

        # # SINGLE-PROCESS FOR DEBUG
        # for i, one_tile in enumerate(tiles_list):
        #     result = compute_save_blends(src_tif_path, ll_path, vis_types, blend_types, one_tile, save_float)
        #     results.append(result)
        #     print(f"Finished tile {i + 1} of {len(tiles_list)}")

        # Collect all paths into a list for each visualisation to be merged
        result_dict = defaultdict(list)
        # Empty dicts (of skipped all-nan tiles) are ignored
        for d in results:
            for key, value in d.items():
                result_dict[key].append(value)
        result_dict = dict(result_dict)

        # Tiles to mosaic and crop to original size (nodata edges)
        for result, result_paths in result_dict.items():
            # Old: Combine tiles into mosaic (old way saved mosaic to file)
            # create_mosaic(result_dict["slrm"], ll_path / "test.tif")

            # Remove file suffix
            result = result[:-4]
            # Now: Use VRT instead of saving another file to disk
            result_parent = result_paths[0].parent.as_posix()
            # Create path for VRT file
            vrt_path = ll_path / f"{result}.vrt"
            # Create VRT from list of files
            if result[-6:] == "_float":
                sf = True
            else:
                sf = False
            mosaic_path = build_vrt(result_paths, vrt_path, save_float=sf)

            # Build a mosaic from VRT if selected #todo: before app starts give warning about size if VRT is not ON
            if not save_vrt:
                # Create mosaic and crop to the original extents
                vrt_to_mosaic(
                    vrt_path=mosaic_path,  # ll_path / "test.tif",
                    output_mosaic_path=ll_path / f"{result}.tif",
                )

                # Delete VRT files
                vrt_path.unlink()
                shutil.rmtree(result_parent)

            # print(results)

    else:
        print("Processing as single tile...")
        # Run the whole image as single tile
        one_tile = None
        _ = compute_save_blends(src_tif_path, ll_path, vis_types, blend_types, one_tile, save_float)

    t1 = time.time() - t0
    print(f"Done with computing blends in {round(t1/60, ndigits=None)} min.")


def compute_save_blends(src_path, low_levels_path, vis_types, blend_types, one_extent, save_float=False):

    # Prepare filenames for saving
    if one_extent:
        # Determine name of the tile (coordinates)
        one_tile_name = f"{one_extent[0]:.0f}_{one_extent[1]:.0f}"
    else:
        one_tile_name = None

    # This dictionary will be used to store save paths
    out_path_dict = {}

    # Filename that will be use RVT built-in function for naming
    filename_rvt = src_path.name

    # ********** COMPUTE LOW-LEVELS *****************************************************

    # Determine req. low-level vis (vis_types from input + req. by blends)
    req_arrays = get_required_arrays(vis_types, blend_types)

    # Default 1 is for GENERAL
    default_1 = rvt.default.DefaultValues()
    # Read from file, path is relative to the Current script directory
    def1_pth = Path(__file__).resolve().parent / "default_1.json"
    default_1.read_default_from_file(def1_pth)

    default_1.fill_no_data = 1
    default_1.keep_original_no_data = 0

    # Default 2 is for FLAT
    default_2 = rvt.default.DefaultValues()
    # Read from file, path is relative to the Current script directory
    def2_pth = Path(__file__).resolve().parent / "default_2.json"
    default_2.read_default_from_file(def2_pth)

    default_2.fill_no_data = 1
    default_2.keep_original_no_data = 0

    # Only compute required visualizations
    in_arrays = compute_low_levels(
        default_1,
        default_2,
        src_path,
        one_extent,
        req_arrays
    )

    # SKIP IF ALL-NANs ARE RETURNED
    if in_arrays["all_nan"]:
        # If all-nan tile encountered, this function will return an empty dict
        return out_path_dict

    # ********** SAVE SELECTED VISUALIZATIONS *****************************************************
    for vis in vis_types:
        # Visualization keyword and name used for get filename can be different
        if vis == "slp":
            vis_name = "slope"
        elif vis == "hs":
            vis_name = "hillshade"
        elif vis == "ld":
            vis_name = "local_dominance"
        else:
            vis_name = vis

        # Determine save path
        rvt_save_name = getattr(default_1, "get_" + vis_name + "_path")(filename_rvt)
        if not one_tile_name:
            # Use RVT naming if this is a single image
            save_path = low_levels_path / rvt_save_name
        else:
            # Use tile naming if this is only one tile
            save_path = low_levels_path / vis / f"{one_tile_name}_rvt_{vis}.tif"
            save_path.parent.mkdir(exist_ok=True)

        # Add path to output dictionary
        out_path_dict[rvt_save_name] = save_path

        # Convert to byte scale and save to disk
        if save_float:
            # Determine file path and create parent folder
            spf = save_path_float(save_path)
            spf.parent.mkdir(exist_ok=True)

            # Save path to the list for creating VRTs
            out_path_dict[rvt_save_name[:-4] + "_float.tif"] = spf

            # Adapt to visualization keywords used in in_arrays
            vis_1 = vis + "_1"

            rasterio_save(
                in_arrays[vis_1],
                in_arrays["profile"],
                save_path=spf,
                nodata=np.nan
            )

        vis_bytscl_save(in_arrays, vis, default_1, save_path)

    # ********** COMPUTE & SAVE SELECTED BLENDS *****************************************************

    # Calculate selected BLENDS
    if "vat_combined" in blend_types:
        # TODO: VAT flat and general may already be calculated (if combined, flat and general are all selected)
        # Determine save path
        save_path, rvt_save_name = save_path_for_blend(
            save_filename="VAT_combined",
            save_dir=low_levels_path,
            source_filename=filename_rvt,
            save_tile_name=one_tile_name
        )
        # Add path to output dictionary
        out_path_dict[rvt_save_name] = save_path
        if save_float:
            rvt_save_name_float = rvt_save_name[:-4] + "_float.tif"
            out_path_dict[rvt_save_name_float] = save_path_float(save_path)
        # Run VAT Combined blend
        in_arrays["vat_combined"] = vat_combined(in_arrays, save_path, save_float)

    if "vat_general" in blend_types:
        # Determine save path
        save_path, rvt_save_name = save_path_for_blend(
            save_filename="VAT_general",
            save_dir=low_levels_path,
            source_filename=filename_rvt,
            save_tile_name=one_tile_name
        )
        # Add path to output dictionary
        out_path_dict[rvt_save_name] = save_path
        if save_float:
            rvt_save_name_float = rvt_save_name[:-4] + "_float.tif"
            out_path_dict[rvt_save_name_float] = save_path_float(save_path)
        # Run VAT general blend
        in_arrays["vat_general"] = vat_general(in_arrays, save_path, save_float)

    if "vat_flat" in blend_types:
        # Determine save path
        save_path, rvt_save_name = save_path_for_blend(
            save_filename="VAT_flat",
            save_dir=low_levels_path,
            source_filename=filename_rvt,
            save_tile_name=one_tile_name
        )
        # Add path to output dictionary
        out_path_dict[rvt_save_name] = save_path
        if save_float:
            rvt_save_name_float = rvt_save_name[:-4] + "_float.tif"
            out_path_dict[rvt_save_name_float] = save_path_float(save_path)
        # Run VAT general blend
        in_arrays["vat_flat"] = vat_flat(in_arrays, save_path, save_float)

    # if "VAT_flat_3B" in blend_types:
    #     # Determine save path
    #     save_path = low_levels_path / "VAT_flat_3B" / f"{one_tile}_rvt_VAT_flat_3B.tif"
    #     save_path.parent.mkdir(exist_ok=True)
    #     in_arrays["vat_flat_3bands"] = vat_flat_3bands(in_arrays, save_path)

    # if "VAT_3B" in blend_types:
    #     # Determine save path
    #     save_path = low_levels_path / "VAT_3B" / f"{one_tile}_rvt_VAT_3B.tif"
    #     save_path.parent.mkdir(exist_ok=True)
    #     in_arrays["vat_3bands"] = vat_3bands(in_arrays, save_path)

    # if "VAT_combined_3B" in blend_types:
    #     # Determine save path
    #     save_path = low_levels_path / "VAT_combined_3B" / f"{one_tile}_rvt_VAT_combined_3B.tif"
    #     save_path.parent.mkdir(exist_ok=True)
    #     vat_combined_3bands(in_arrays, save_path)

    if "rrim" in blend_types:
        # Determine save path
        save_path, rvt_save_name = save_path_for_blend(
            save_filename="RRIM",
            save_dir=low_levels_path,
            source_filename=filename_rvt,
            save_tile_name=one_tile_name
        )
        # Add path to output dictionary
        out_path_dict[rvt_save_name] = save_path
        if save_float:
            rvt_save_name_float = rvt_save_name[:-4] + "_float.tif"
            out_path_dict[rvt_save_name_float] = save_path_float(save_path)
        # Run RRIM blend
        in_arrays["rrim"] = blend_rrim(in_arrays, save_path, save_float)

    if "e2MSTP" in blend_types:
        if "rrim" not in in_arrays.keys():
            # We need RRIM as component of e2MSTP, but don't need to save it
            in_arrays["rrim"] = blend_rrim(in_arrays)
        # Determine save path
        save_path, rvt_save_name = save_path_for_blend(
            save_filename="e2MSTP",
            save_dir=low_levels_path,
            source_filename=filename_rvt,
            save_tile_name=one_tile_name
        )
        # Add path to output dictionary
        out_path_dict[rvt_save_name] = save_path
        if save_float:
            rvt_save_name_float = rvt_save_name[:-4] + "_float.tif"
            out_path_dict[rvt_save_name_float] = save_path_float(save_path)
        # Run e2MSTP blend
        blend_e2mstp(in_arrays, save_path, save_float)

    if "crim" in blend_types:
        # Determine save path
        save_path, rvt_save_name = save_path_for_blend(
            save_filename="CRIM",
            save_dir=low_levels_path,
            source_filename=filename_rvt,
            save_tile_name=one_tile_name
        )
        # Add path to output dictionary
        out_path_dict[rvt_save_name] = save_path
        if save_float:
            rvt_save_name_float = rvt_save_name[:-4] + "_float.tif"
            out_path_dict[rvt_save_name_float] = save_path_float(save_path)
        # Run CRIM blend
        in_arrays["crim"] = blend_crim(in_arrays, save_path, save_float)

    if "e3MSTP" in blend_types:
        # Check if CRIM was already calculated
        if "crim" not in in_arrays.keys():
            in_arrays["crim"] = blend_crim(in_arrays)
        # Determine save path
        save_path, rvt_save_name = save_path_for_blend(
            save_filename="e3MSTP",
            save_dir=low_levels_path,
            source_filename=filename_rvt,
            save_tile_name=one_tile_name
        )
        # Add path to output dictionary
        out_path_dict[rvt_save_name] = save_path
        if save_float:
            rvt_save_name_float = rvt_save_name[:-4] + "_float.tif"
            out_path_dict[rvt_save_name_float] = save_path_float(save_path)
        # Run e3MSTP blend
        blend_e3mstp(in_arrays, save_path, save_float)

    if "e4MSTP" in blend_types:
        # Determine save path
        save_path, rvt_save_name = save_path_for_blend(
            save_filename="e4MSTP",
            save_dir=low_levels_path,
            source_filename=filename_rvt,
            save_tile_name=one_tile_name
        )
        # Add path to output dictionary
        out_path_dict[rvt_save_name] = save_path
        if save_float:
            rvt_save_name_float = rvt_save_name[:-4] + "_float.tif"
            out_path_dict[rvt_save_name_float] = save_path_float(save_path)
        # Run e4MSTP blend
        blend_e4mstp(in_arrays, save_path, save_float)

    return out_path_dict


def vat_general(dict_arrays, save_path=None, save_float=False):
    # BLEND VAT GENERAL
    vat_combination_general = rvt.blend.BlenderCombination()
    vat_combination_general.create_layer(
        vis_method="Sky-View Factor",
        normalization="Value",
        minimum=0.7,
        maximum=1.0,
        blend_mode="Multiply",
        opacity=25,
        image=dict_arrays['svf_1'].squeeze()
    )
    vat_combination_general.create_layer(
        vis_method="Openness - Positive",
        normalization="Value",
        minimum=68,
        maximum=93,
        blend_mode="Overlay",
        opacity=50,
        image=dict_arrays['opns_1'].squeeze()
    )
    vat_combination_general.create_layer(
        vis_method="Slope gradient",
        normalization="Value",
        minimum=0,
        maximum=50,
        blend_mode="Luminosity",
        opacity=50,
        image=dict_arrays['slp_1'].squeeze()
    )
    vat_combination_general.create_layer(
        vis_method="Hillshade",
        normalization="Value",
        minimum=0,
        maximum=1,
        blend_mode="Normal",
        opacity=100,
        image=dict_arrays['hs_1'].squeeze()
    )
    vat_1 = vat_combination_general.render_all_images(
        save_visualizations=False,
        save_render_path=None,
        no_data=np.nan
    )
    out_vat_general = vat_1.astype("float32")

    # Save GeoTIF
    out_profile = dict_arrays['profile'].copy()

    if save_float:
        # out_profile.update(dtype='uint8')
        rasterio_save(
            out_vat_general,
            out_profile,
            save_path=save_path_float(save_path),
            nodata=np.nan
        )

    if save_path:
        # Convert to 8bit image
        out_vat_general_8bit = rvt.vis.byte_scale(
            out_vat_general,
            c_min=0,
            c_max=1
        )

        out_profile.update(dtype='uint8')
        rasterio_save(
            out_vat_general_8bit,
            out_profile,
            save_path=save_path,
            nodata=None
        )

    return out_vat_general


def vat_flat(dict_arrays, save_path=None, save_float=False):
    # BLEND VAT FLAT
    vat_combination_flat = rvt.blend.BlenderCombination()
    vat_combination_flat.create_layer(
        vis_method="Sky-View Factor",
        normalization="Value",
        minimum=0.9,
        maximum=1.0,
        blend_mode="Multiply",
        opacity=25,
        image=dict_arrays['svf_2'].squeeze()
    )
    vat_combination_flat.create_layer(
        vis_method="Openness - Positive",
        normalization="Value",
        minimum=85,
        maximum=93,
        blend_mode="Overlay",
        opacity=50,
        image=dict_arrays['opns_2'].squeeze()
    )
    vat_combination_flat.create_layer(
        vis_method="Slope gradient",
        normalization="Value",
        minimum=0,
        maximum=15,
        blend_mode="Luminosity",
        opacity=50,
        image=dict_arrays['slp_1'].squeeze()
    )
    vat_combination_flat.create_layer(
        vis_method="Hillshade",
        normalization="Value",
        minimum=0,
        maximum=1,
        blend_mode="Normal",
        opacity=100,
        image=dict_arrays['hs_2'].squeeze()
    )
    vat_2 = vat_combination_flat.render_all_images(
        save_visualizations=False,
        save_render_path=None,
        no_data=np.nan
    )
    out_vat_flat = vat_2.astype("float32")

    # Save GeoTIF
    out_profile = dict_arrays['profile'].copy()

    if save_float:
        # out_profile.update(dtype='uint8')
        rasterio_save(
            out_vat_flat,
            out_profile,
            save_path=save_path_float(save_path),
            nodata=np.nan
        )

    if save_path:
        # Convert to 8bit image
        out_vat_flat_8bit = rvt.vis.byte_scale(
            out_vat_flat,
            c_min=0,
            c_max=1
        )

        out_profile.update(dtype='uint8')
        rasterio_save(
            out_vat_flat_8bit,
            out_profile,
            save_path=save_path,
            nodata=None
        )

    return out_vat_flat


def vat_combined(dict_arrays, save_path, save_float=False):
    """
    VAT Combined 8bit
    """
    # BLEND VAT GENERAL
    vat_1 = vat_general(dict_arrays)

    # BLEND VAT FLAT
    vat_2 = vat_flat(dict_arrays)

    # BLEND VAT COMBINED
    comb_vat_combined = rvt.blend.BlenderCombination()
    comb_vat_combined.create_layer(
        vis_method="vat_general", normalization="value",
        minimum=0, maximum=1,
        blend_mode="normal", opacity=50,
        image=vat_1
    )
    comb_vat_combined.create_layer(
        vis_method="vat_flat", normalization="value",
        minimum=0, maximum=1,
        blend_mode="normal", opacity=100,
        image=vat_2
    )
    out_vat_combined = comb_vat_combined.render_all_images(
        save_visualizations=False,
        save_render_path=None,
        no_data=np.nan
    )

    # Save GeoTIF
    out_profile = dict_arrays['profile'].copy()

    if save_float:
        # out_profile.update(dtype='uint8')
        rasterio_save(
            out_vat_combined,
            out_profile,
            save_path=save_path_float(save_path),
            nodata=np.nan
        )

    # Convert to 8bit image
    out_vat_combined_8bit = rvt.vis.byte_scale(
        out_vat_combined,
        c_min=0,
        c_max=1
    )

    out_profile.update(dtype='uint8')
    rasterio_save(
        out_vat_combined_8bit,
        out_profile,
        save_path=save_path,
        nodata=None
    )

    return out_vat_combined


def vat_flat_3bands(dict_arrays, save_path):
    """
    SVF (normalised 0.9 – 1)
    Openness positive (normalised 85 – 93)
    Slope gradient (normalised 0° – 15°)
    """

    svf = normalize_image(
        visualization="sky-view factor",
        image=dict_arrays["svf_2"].squeeze(),
        min_norm=0.9,
        max_norm=1,
        normalization="value"
    )
    opns = normalize_image(
        visualization="openness - positive",
        image=dict_arrays["opns_2"].squeeze(),
        min_norm=85,
        max_norm=93,
        normalization="value"
    )
    slope = normalize_image(
        visualization="slope gradient",
        image=dict_arrays["slp_1"].squeeze(),
        min_norm=0,
        max_norm=15,
        normalization="value"
    )
    out_flat_vat3 = np.stack(
        [slope, svf, opns],
        axis=0, out=None
    )

    # Save GeoTIF
    rasterio_save(
        out_flat_vat3,
        dict_arrays['profile'],
        save_path=save_path,
        nodata=np.nan
    )

    return out_flat_vat3


def vat_3bands(dict_arrays, save_path):
    """
    SVF (normalised 0.7 – 1)
    Openness positive (normalised 68 – 93)
    Slope gradient (normalised 0° – 50°)
    """

    svf = normalize_image(
        visualization="sky-view factor",
        image=dict_arrays["svf_1"].squeeze(),
        min_norm=0.7,
        max_norm=1,
        normalization="value"
    )
    opns = normalize_image(
        visualization="openness - positive",
        image=dict_arrays["opns_1"].squeeze(),
        min_norm=68,
        max_norm=93,
        normalization="value"
    )
    slope = normalize_image(
        visualization="slope gradient",
        image=dict_arrays["slp_1"].squeeze(),
        min_norm=0,
        max_norm=50,
        normalization="value"
    )
    out_vat3 = np.stack(
        [slope, svf, opns],
        axis=0, out=None
    )

    # Save GeoTIF
    rasterio_save(
        out_vat3,
        dict_arrays['profile'],
        save_path=save_path,
        nodata=np.nan
    )

    return out_vat3


def vat_combined_3bands(dict_arrays, save_path):
    out_vat_combined_3bands = np.zeros(dict_arrays["vat_3bands"].shape)
    for i in range(3):
        comb_vat_combined_3bands = rvt.blend.BlenderCombination()
        comb_vat_combined_3bands.create_layer(
            vis_method="band1_vat_3bands", normalization="value",
            minimum=0, maximum=1,
            blend_mode="normal", opacity=50,
            image=dict_arrays["vat_3bands"][i, :, :].squeeze()
        )
        comb_vat_combined_3bands.create_layer(
            vis_method="band1_vat_flat_3bands", normalization="value",
            minimum=0, maximum=1,
            blend_mode="normal", opacity=100,
            image=dict_arrays["vat_flat_3bands"][i, :, :].squeeze()
        )
        out_vat_combined_3bands[i, :, :] = comb_vat_combined_3bands.render_all_images(
            save_visualizations=False,
            save_render_path=None,
            no_data=np.nan)

    out_vat_combined_3bands = out_vat_combined_3bands.astype("float32")
    # Save GeoTIF
    rasterio_save(
        out_vat_combined_3bands,
        dict_arrays['profile'],
        save_path=save_path,
        nodata=np.nan
    )

    return out_vat_combined_3bands


def vis_bytscl_save(image_arrays, visualization, defaults, save_path):
    # Adapt to visualization keywords used in in_arrays
    vis_1 = visualization + "_1"

    # Determine min/max values for norm from default.json "bytscl"
    if visualization == "opns":
        bytscl_value = getattr(defaults, "pos_opns_bytscl")
    else:
        bytscl_value = getattr(defaults, visualization + "_bytscl")

    # Use RVT function for normalization
    norm_image = normalize_image(
        visualization=visualization,
        image=image_arrays[vis_1].squeeze(),
        min_norm=bytscl_value[1],
        max_norm=bytscl_value[2],
        normalization=bytscl_value[0]
    )

    # Use RVT function for bytscale
    out_image = rvt.vis.byte_scale(
        norm_image,
        c_min=0,
        c_max=1
    )

    # Save GeoTIF
    out_profile = image_arrays['profile'].copy()
    out_profile.update(dtype='uint8')
    rasterio_save(
        out_image,
        out_profile,
        save_path=save_path,
        nodata=None
    )
    return out_image


def blend_rrim(dict_arrays, save_path=None, save_float=False):
    comb_rrim = rvt.blend.BlenderCombination()
    comb_rrim.create_layer(vis_method="Slope gradient", normalization="Value",
                           minimum=0, maximum=45,
                           blend_mode="Normal", opacity=50,
                           colormap="Reds_r", min_colormap_cut=0, max_colormap_cut=1,
                           image=dict_arrays['slp_1'].squeeze()
                           )
    comb_rrim.create_layer(vis_method="Opns_Pos_Neg/2", normalization="Value",
                           minimum=-25, maximum=25,
                           blend_mode="Normal", opacity=100,
                           colormap="Greys_r", min_colormap_cut=0, max_colormap_cut=1,
                           image=((dict_arrays['opns_1'] - dict_arrays['neg_opns_1'])/2).squeeze()
                           )
    out_rrim = comb_rrim.render_all_images(save_visualizations=False,
                                           save_render_path=None,
                                           no_data=np.nan)
    out_rrim = out_rrim.astype("float32")

    # Save GeoTIF
    if save_path:
        out_profile = dict_arrays['profile'].copy()

        if save_float:
            # out_profile.update(dtype='uint8')
            rasterio_save(
                out_rrim,
                out_profile,
                save_path=save_path_float(save_path),
                nodata=np.nan
            )

        # Convert to 8bit image
        out_rrim_8bit = rvt.vis.byte_scale(
            out_rrim,
            c_min=0,
            c_max=1
        )

        out_profile.update(dtype='uint8')
        rasterio_save(
            out_rrim_8bit,
            out_profile,
            save_path=save_path,
            nodata=None
        )

    return out_rrim


def blend_e2mstp(dict_arrays, save_path, save_float=False):
    comb_e2mstp = rvt.blend.BlenderCombination()
    comb_e2mstp.create_layer(vis_method="slrm", normalization="value",
                             minimum=-0.5, maximum=0.5,
                             blend_mode="screen", opacity=25,
                             image=dict_arrays["slrm_1"].squeeze()
                             )
    comb_e2mstp.create_layer(vis_method="rrim", normalization="value",
                             minimum=0, maximum=1,
                             blend_mode="soft_light", opacity=70,
                             image=dict_arrays["rrim"]
                             )
    comb_e2mstp.create_layer(vis_method="mstp", normalization="value",
                             minimum=0, maximum=1,
                             blend_mode="normal", opacity=100,
                             image=dict_arrays["mstp_1"]
                             )
    out_e2mstp = comb_e2mstp.render_all_images(save_visualizations=False,
                                               save_render_path=None,
                                               no_data=np.nan)
    out_e2mstp = out_e2mstp.astype("float32")
    out_e2mstp[np.isnan(dict_arrays["rrim"])] = np.nan
    out_e2mstp[out_e2mstp > 1] = 1

    # Save GeoTIF
    out_profile = dict_arrays['profile'].copy()

    if save_float:
        # out_profile.update(dtype='uint8')
        rasterio_save(
            out_e2mstp,
            out_profile,
            save_path=save_path_float(save_path),
            nodata=np.nan
        )

    # Convert to 8bit image
    out_8bit = rvt.vis.byte_scale(
        out_e2mstp,
        c_min=0,
        c_max=1
    )

    out_profile.update(dtype='uint8')
    rasterio_save(
        out_8bit,
        out_profile,
        save_path=save_path,
        nodata=None
    )


def blend_crim(dict_arrays, save_path=None, save_float=False):
    comb_crim = rvt.blend.BlenderCombination()
    comb_crim.create_layer(vis_method="Openness_Pos-Neg", normalization="Value",
                           minimum=-28, maximum=28,
                           blend_mode="overlay", opacity=50,
                           image=(dict_arrays['opns_1'] - dict_arrays['neg_opns_1']).squeeze()
                           )
    comb_crim.create_layer(vis_method="Openness_Pos-Neg", normalization="Value",
                           minimum=-28, maximum=28,
                           blend_mode="luminosity", opacity=50,
                           image=(dict_arrays['opns_1'] - dict_arrays['neg_opns_1']).squeeze()
                           )
    comb_crim.create_layer(vis_method="slope gradient red", normalization="Value",
                           minimum=0, maximum=45,
                           blend_mode="normal", opacity=100,
                           colormap="OrRd", min_colormap_cut=0, max_colormap_cut=1,
                           image=dict_arrays['slp_1'].squeeze()
                           )
    out_crim = comb_crim.render_all_images(save_visualizations=False,
                                           save_render_path=None,
                                           no_data=np.nan)
    out_crim = out_crim.astype("float32")

    # Save GeoTIF
    if save_path:
        out_profile = dict_arrays['profile'].copy()

        if save_float:
            # out_profile.update(dtype='uint8')
            rasterio_save(
                out_crim,
                out_profile,
                save_path=save_path_float(save_path),
                nodata=np.nan
            )

        # Convert to 8bit image
        out_crim_8bit = rvt.vis.byte_scale(
            out_crim,
            c_min=0,
            c_max=1
        )

        # Save GeoTIF
        out_profile.update(dtype='uint8')
        rasterio_save(
            out_crim_8bit,
            out_profile,
            save_path=save_path,
            nodata=None
        )

    return out_crim


def blend_e3mstp(dict_arrays, save_path, save_float=False):
    comb_e3mstp = rvt.blend.BlenderCombination()
    comb_e3mstp.create_layer(vis_method="slrm", normalization="value",
                             minimum=-0.5, maximum=0.5,
                             blend_mode="screen", opacity=25,
                             image=dict_arrays["slrm_1"].squeeze()
                             )
    comb_e3mstp.create_layer(vis_method="crim", normalization="value",
                             minimum=0, maximum=1,
                             blend_mode="soft_light", opacity=70,
                             image=dict_arrays["crim"]
                             )
    comb_e3mstp.create_layer(vis_method="mstp", normalization="value",
                             minimum=0, maximum=1,
                             blend_mode="normal", opacity=100,
                             image=dict_arrays["mstp_1"]
                             )
    out_e3mstp = comb_e3mstp.render_all_images(save_visualizations=False,
                                               save_render_path=None,
                                               no_data=np.nan)
    out_e3mstp = out_e3mstp.astype("float32")
    out_e3mstp[np.isnan(dict_arrays["crim"])] = np.nan
    out_e3mstp[out_e3mstp > 1] = 1

    # Save GeoTIF
    out_profile = dict_arrays['profile'].copy()

    if save_float:
        # out_profile.update(dtype='uint8')
        rasterio_save(
            out_e3mstp,
            out_profile,
            save_path=save_path_float(save_path),
            nodata=np.nan
        )

    # Convert to 8bit image
    out_e3mstp = rvt.vis.byte_scale(
        out_e3mstp,
        c_min=0,
        c_max=1
    )

    out_profile.update(dtype='uint8')
    rasterio_save(
        out_e3mstp,
        out_profile,
        save_path=save_path,
        nodata=None
    )


def blend_e4mstp(dict_arrays, save_path, save_float=False):
    # Get SVF combined
    dict_arrays['svf_combined'] = blend_svf_combined(dict_arrays)
    # Get Openness + LD
    dict_arrays['opns_ld'] = blend_opns_ld(dict_arrays)

    comb_nv = rvt.blend.BlenderCombination()
    comb_nv.create_layer(
        vis_method="mstp",
        normalization="value", minimum=0, maximum=1,
        blend_mode="overlay", opacity=90,
        image=dict_arrays['mstp_1']
    )
    comb_nv.create_layer(
        vis_method="Comb svf",
        normalization="value", minimum=-0.5, maximum=0.5,
        blend_mode="multiply", opacity=25,
        image=dict_arrays['svf_combined']
    )
    comb_nv.create_layer(
        vis_method="Comb openness LD",
        normalization="value", minimum=0, maximum=1,
        blend_mode="multiply", opacity=100,
        image=dict_arrays['opns_ld']
    )
    comb_nv.create_layer(
        vis_method="Slope gradient",
        normalization="value", minimum=0, maximum=55,
        blend_mode="normal", opacity=100,
        colormap="Reds_r", min_colormap_cut=0, max_colormap_cut=1,
        image=dict_arrays['slp_1'].squeeze()
    )
    out_e4mstp = comb_nv.render_all_images(
        save_visualizations=False,
        save_render_path=None,
        no_data=np.nan
    )
    out_e4mstp = out_e4mstp.astype("float32")
    out_e4mstp[np.isnan(dict_arrays['mstp_1'])] = np.nan
    out_e4mstp[out_e4mstp > 1] = 1

    # Save GeoTIF
    out_profile = dict_arrays['profile'].copy()

    if save_float:
        # out_profile.update(dtype='uint8')
        rasterio_save(
            out_e4mstp,
            out_profile,
            save_path=save_path_float(save_path),
            nodata=np.nan
        )

    # Convert to 8bit image
    out_e4mstp = rvt.vis.byte_scale(
        out_e4mstp,
        c_min=0,
        c_max=1
    )

    out_profile.update(dtype='uint8')
    rasterio_save(
        out_e4mstp,
        out_profile,
        save_path=save_path,
        nodata=None
    )

    return out_e4mstp


def blend_opns_ld(dict_arrays, save_path=None):
    """Required dict_arrays:
    - opns
    - neg_opns
    - ld
    """
    comb = rvt.blend.BlenderCombination()
    comb.create_layer(
        vis_method="Openness difference", normalization="Value",
        minimum=-15, maximum=15,
        blend_mode="normal", opacity=50,
        image=(dict_arrays['opns_1'] - dict_arrays['neg_opns_1']).squeeze()
    )
    comb.create_layer(
        vis_method="Local dominance", normalization="Value",
        minimum=0.5, maximum=1.8,
        blend_mode="normal", opacity=100,
        image=dict_arrays['ld_1'].squeeze()
    )
    opns_ld = comb.render_all_images(
        save_visualizations=False,
        save_render_path=None,
        no_data=np.nan
    )
    opns_ld = opns_ld.astype("float32")

    if save_path:
        # Save GeoTIF
        rasterio_save(
            opns_ld,
            dict_arrays['profile'],
            save_path=save_path,
            nodata=np.nan
        )

    return opns_ld


def blend_svf_combined(dict_arrays, save_path=None):
    """Required dict_arrays:
    - svf_1
    - svf_2
    """
    comb_svf = rvt.blend.BlenderCombination()
    comb_svf.create_layer(vis_method="Sky-view factor", normalization="Value",
                          minimum=0.7, maximum=1,
                          blend_mode="normal", opacity=50,
                          image=dict_arrays['svf_1'].squeeze()
                          )
    comb_svf.create_layer(vis_method="Sky-view factor", normalization="Value",
                          minimum=0.9, maximum=1,
                          blend_mode="normal", opacity=100,
                          image=dict_arrays['svf_2'].squeeze()
                          )
    cs_svf = comb_svf.render_all_images(save_visualizations=False,
                                        save_render_path=None,
                                        no_data=np.nan)
    cs_svf = cs_svf.astype("float32")

    if save_path:
        # Save GeoTIF
        rasterio_save(
            cs_svf,
            dict_arrays['profile'],
            save_path=save_path,
            nodata=np.nan
        )

    return cs_svf


def rasterio_save(array, profile, save_path, nodata=None):
    if len(array.shape) == 2:
        array = np.expand_dims(array, axis=0)
    profile.update(dtype=array.dtype,
                   count=array.shape[0],
                   nodata=nodata,
                   compress="LZW",
                   predictor=2)
    save_path.parent.mkdir(exist_ok=True)
    with rasterio.open(save_path, "w", **profile) as dst:
        dst.write(array)


def get_required_arrays(vis_types, blend_types):
    # Initialize dict with all possible visualizations
    # NOTE: The keys with "_1" have to match the input values of visualizations in GUI!!!
    req_arrays = {
        # These are visualizations (also GENERAL for VAT):
        "slp_1": False,
        "hs_1": False,
        # "mhs_1": False,
        "slrm_1": False,
        "svf_1": False,  # large = FLAT = 5m
        "opns_1": False,
        "neg_opns_1": False,
        "ld_1": False,
        # "sky_illumination_1": False,
        # "shadow_horizon_1": False,
        # "msrm_1": False,
        "mstp_1": False,

        # Flat terrain:
        "hs_2": False,
        "svf_2": False,  # large = FLAT = 10m
        "opns_2": False,
        "neg_opns_2": False
    }

    # Update dictionary based on given visualizations:
    for key in vis_types:
        key = key + "_1"
        if key in req_arrays:
            req_arrays[key] = True

    # Update dictionary based on given blends:
    # if "VAT_3B" in blend_types:
    #     req_arrays["svf_1"] = True
    #     req_arrays["opns_1"] = True
    #     req_arrays["slp_1"] = True

    # if "VAT_flat_3B" in blend_types:
    #     req_arrays["svf_2"] = True
    #     req_arrays["opns_2"] = True
    #     req_arrays["slp_1"] = True

    # if "VAT_combined_3B" in blend_types:
    #     req_arrays["svf_2"] = True
    #     req_arrays["opns_2"] = True
    #     req_arrays["svf_1"] = True
    #     req_arrays["opns_1"] = True
    #     req_arrays["slp_1"] = True

    if ("e3MSTP" in blend_types) or ("e2MSTP" in blend_types):
        req_arrays["slrm_1"] = True
        req_arrays["mstp_1"] = True
        req_arrays["slp_1"] = True
        req_arrays["opns_1"] = True
        req_arrays["neg_opns_1"] = True

    if "e4MSTP" in blend_types:
        req_arrays["ld_1"] = True
        req_arrays["svf_1"] = True,
        req_arrays["mstp_1"] = True
        req_arrays["svf_2"] = True
        req_arrays["opns_2"] = True
        req_arrays["opns_1"] = True
        req_arrays["neg_opns_1"] = True
        req_arrays["slp_1"] = True

    if "vat_combined" in blend_types:
        req_arrays["svf_2"] = True
        req_arrays["opns_2"] = True
        req_arrays["svf_1"] = True
        req_arrays["opns_1"] = True
        req_arrays["slp_1"] = True
        req_arrays["hs_1"] = True
        req_arrays["hs_2"] = True

    if "vat_general" in blend_types:
        req_arrays["svf_1"] = True
        req_arrays["opns_1"] = True
        req_arrays["slp_1"] = True
        req_arrays["hs_1"] = True

    if "vat_flat" in blend_types:
        req_arrays["svf_2"] = True
        req_arrays["opns_2"] = True
        req_arrays["slp_1"] = True
        req_arrays["hs_2"] = True

    if "rrim" in blend_types:
        req_arrays["slp_1"] = True
        req_arrays["opns_1"] = True
        req_arrays["neg_opns_1"] = True

    if "new_blend" in blend_types:
        req_arrays["slrm_1"] = True

    return req_arrays


def compute_low_levels(
        default_1,
        default_2,
        vrt_path,
        input_dem_extents,
        vis_types
):
    # Read buffer values from defaults (already changed from meters to pixels!!!)!
    all_buffers = {
        "slp_1": 0,
        "slrm_1": default_1.slrm_rad_cell,
        "ld_1": default_1.ld_max_rad,
        "mstp_1": default_1.mstp_broad_scale[1],

        "hs_1": 1,
        "hs_2": 1,

        "svf_GEN": default_1.svf_r_max,  # SVF, OPNS+ and OPNS- for (GENERAL = SMALL = 5m)

        "svf_FLAT": default_2.svf_r_max,  # SVF, OPNS+ and OPNS- for (FLAT = LARGE = 10m)
    }

    # Get required visualizations (account for SVF, opns and neg.opns, they are calculated in the same function)
    req_visualizations = [a for (a, v) in vis_types.items() if v]
    if any(value in req_visualizations for value in ["svf_1", "opns_1", "neg_opns_1"]):
        req_visualizations.append("svf_GEN")
    if any(value in req_visualizations for value in ["svf_2", "opns_2", "neg_opns_2"]):
        req_visualizations.append("svf_FLAT")

    # Filter buffer_dict based on required visualizations
    buffer_dict = {key: all_buffers[key] for key in req_visualizations if key in all_buffers}

    # Select the largest required buffer
    max_buff = max(buffer_dict, key=buffer_dict.get)
    buffer = buffer_dict[max_buff]

    # Read array into RVT dictionary format
    dict_arrays = get_raster_vrt(vrt_path, input_dem_extents, buffer)

    # Change nodata value to np.nan, to avoid problems later
    dict_arrays["array"][dict_arrays["array"] == dict_arrays["no_data"]] = np.nan
    dict_arrays["no_data"] = np.nan

    # Skip if all pixels are nodata (remove buffer when checking)
    all_nan_check = np.all(np.isnan(dict_arrays["array"][buffer: -buffer, buffer: -buffer]))
    # True means we have all-nan array and have to skip this tile
    if all_nan_check:
        # SKIP THIS TILE
        dict_arrays["all_nan"] = True
        return dict_arrays
    else:
        dict_arrays["all_nan"] = False

    # --- START VISUALIZATION WITH RVT ---
    vis_out = dict()
    for vis_type in buffer_dict:
        # Obtain buffer for current visualization type
        arr_buff = buffer_dict[vis_type]
        # Slice raster to minimum required size
        arr_slice = buffer_dict[max_buff] - arr_buff
        if arr_slice == 0:
            sliced_arr = dict_arrays["array"]
        else:
            sliced_arr = dict_arrays["array"][arr_slice:-arr_slice, arr_slice:-arr_slice]

        # Run visualization
        if vis_type == "slp_1":
            vis_out = {
                vis_type: default_1.get_slope(
                    sliced_arr,
                    resolution_x=dict_arrays["resolution"][0],
                    resolution_y=dict_arrays["resolution"][1]
                )
            }
        elif vis_type == "slrm_1":
            vis_out = {
                vis_type: default_1.get_slrm(sliced_arr)
            }
        elif vis_type == "ld_1":
            vis_out = {
                vis_type: default_1.get_local_dominance(sliced_arr)
            }
        elif vis_type == "mstp_1":
            # test = default_1.get_slrm(sliced_arr)
            vis_out = {
                # vis_type: np.stack((test, test, test), 0)
                vis_type: default_1.get_mstp(sliced_arr)
            }
        elif vis_type == "hs_1":
            vis_out = {
                vis_type: default_1.get_hillshade(
                    sliced_arr,
                    resolution_x=dict_arrays["resolution"][0],
                    resolution_y=dict_arrays["resolution"][1]
                )
            }
        elif vis_type == "hs_2":
            vis_out = {
                vis_type: default_2.get_hillshade(
                    sliced_arr,
                    resolution_x=dict_arrays["resolution"][0],
                    resolution_y=dict_arrays["resolution"][1]
                )
            }
        elif vis_type == "svf_GEN":  # small, GENERAL, 5m
            # Check which of the 3 to be computed
            compute_svf = True if "svf_1" in req_visualizations else False
            compute_opns = True if "opns_1" in req_visualizations else False
            compute_neg_opns = True if "neg_opns_1" in req_visualizations else False

            vis_out = {}

            if compute_svf or compute_opns:
                vis_out = default_1.get_sky_view_factor(
                    sliced_arr,
                    dict_arrays["resolution"][0],
                    compute_svf=compute_svf,
                    compute_opns=compute_opns
                )
                # Rename to correct vis_type (which svf is this?)
                for k in list(vis_out.keys()):
                    vis_out[f"{k}_1"] = vis_out.pop(k)

            if compute_neg_opns:
                vis_out["neg_opns_1"] = default_1.get_neg_opns(
                    sliced_arr,
                    dict_arrays["resolution"][0]
                )
        elif vis_type == "svf_FLAT":  # large, FLAT, 10m
            # Check which of the 3 to be computed
            compute_svf = True if "svf_2" in req_visualizations else False
            compute_opns = True if "opns_2" in req_visualizations else False
            compute_neg_opns = True if "neg_opns_2" in req_visualizations else False

            vis_out = {}

            if compute_svf or compute_opns:
                vis_out = default_2.get_sky_view_factor(
                    sliced_arr,
                    dict_arrays["resolution"][0],
                    compute_svf=compute_svf,
                    compute_opns=compute_opns
                )
                # Rename to correct vis_type (which svf is this?)
                for k in list(vis_out.keys()):
                    vis_out[f"{k}_2"] = vis_out.pop(k)

            if compute_neg_opns:
                vis_out["neg_opns_2"] = default_2.get_neg_opns(
                    sliced_arr,
                    dict_arrays["resolution"][0]
                )
        else:
            raise ValueError("Wrong vis_type in the visualization for loop")

        # Remove buffer and Store visualization in dictionary
        for i, array in vis_out.items():
            # Slice away buffer
            if arr_buff == 0:
                arr_out = array
            else:
                arr_out = array[..., arr_buff:-arr_buff, arr_buff:-arr_buff]

            # Make sure the dimensions of array are correct
            if arr_out.ndim == 2:
                arr_out = np.expand_dims(arr_out, axis=0)

            # Add to results dictionary
            dict_arrays[i] = arr_out

    return dict_arrays


def get_raster_vrt(vrt_path, extents, buffer):
    """
    Extents have to be transformed into rasterio Window object, it is passed into the function as tuple.
    (left, bottom, right, top)


    Parameters
    ----------
    vrt_path : str
        Path to raster file. Can be any rasterio readable format.
    extents : tuple
        Extents to be read (left, bottom, right, top).
    buffer : int
        Buffer in pixels.

    Returns
    -------
        A dictionary containing the raster array and all the required metadata.

    """
    with rasterio.open(vrt_path) as vrt:
        # If extents are not given, use source extents
        if not extents:
            extents = list(vrt.bounds)

        # Read VRT metadata
        vrt_res = vrt.res
        vrt_nodata = vrt.nodata
        vrt_transform = vrt.transform
        vrt_crs = vrt.crs
        nodata_val = vrt.nodata

        # ADD BUFFER TO EXTENTS (LBRT) - transform pixels to meters!
        buffer_m = buffer * vrt_res[0]
        buff_extents = (
            extents[0] - buffer_m,
            extents[1] - buffer_m,
            extents[2] + buffer_m,
            extents[3] + buffer_m
        )

        # Pack extents into rasterio's Window object
        buff_window = from_bounds(*buff_extents, vrt_transform)
        orig_window = from_bounds(*extents, vrt_transform)

        # Read windowed array (with added buffer)
        # boundless - if window falls out of bounds, read it and fill with NaNs
        win_array = vrt.read(window=buff_window, boundless=True)

        # Deal with nodata
        if nodata_val is not None and not np.isnan(nodata_val):
            # Ensure array is float to support np.nan
            win_array = win_array.astype("float32", copy=False)
            win_array[win_array == nodata_val] = np.nan

        # Save transform object of both extents (original and buffered)
        orig_transform = vrt.window_transform(orig_window)

    # For raster with only one band, remove first axis from the array (RVT requirement)
    if win_array.shape[0] == 1:
        win_array = np.squeeze(win_array, axis=0)

    # Prepare output metadata profile
    out_profile = {
        'driver': 'GTiff',
        'nodata': np.nan,
        'width':  win_array.shape[1] - 2 * buffer,
        'height':  win_array.shape[0] - 2 * buffer,
        'count':  1,
        'crs': vrt_crs,
        'transform': orig_transform,
        "compress": "lzw"
    }

    output = {
        "array": win_array,
        "resolution": vrt_res,
        "no_data": vrt_nodata,
        "orig_transform": orig_transform,
        "crs": vrt_crs,
        "profile": out_profile
    }

    return output


def build_vrt(tif_list, vrt_path, save_float=False):
    """Creates a vrt file from a list of tif files. The user has to provide the path of the final VRT file.

    By default, it saves byte arrays setting all nodata pixels to 255 and metadata as nodata=None. If user wants to save
    float array, this needs to be specified by setting save_float flat to True, in this case np.nan is used for nodata
    and

    """
    # Make sure all paths are in pathlib format
    vrt_path = Path(vrt_path)
    tif_list = [Path(a) for a in tif_list]

    # Different nodata for different dtypes
    if save_float:
        vrt_nodata = "nan"  # TODO: Check if this is correct way to set Nan
        hide_nodata = False
    else:
        vrt_nodata = 255
        hide_nodata = True

    # Change paths in tif list to strings
    tif_list = [a.as_posix() for a in tif_list]

    # Create VRT
    vrt_options = gdal.BuildVRTOptions(VRTNodata=vrt_nodata, hideNodata=hide_nodata)
    my_vrt = gdal.BuildVRT(vrt_path.as_posix(), tif_list, options=vrt_options)
    my_vrt = None

    return vrt_path


def save_path_for_blend(save_filename: str, save_dir, source_filename, save_tile_name=None):
    """
    Two options for creating save path:

    - if single TIF, save the raster using the path to the final save file (located in the same folder as source file
        and using RVT naming conventions). In this case the save_tile_name variable is given as None

    - in the case where we want to save only the one tile, use the constructed save_tile_name and save into a child
        directory of the same name as the final name of the save file
    """
    save_dir = Path(save_dir)
    source_filename = Path(source_filename)
    rvt_save_name = f"{source_filename.stem}_{save_filename}.tif"
    # Determine save path
    if not save_tile_name:
        # Use RVT naming if this is a single image
        save_path = save_dir / rvt_save_name
    else:
        # Use tile naming if this is only one tile
        save_path = save_dir / save_filename / f"{save_tile_name}_rvt_{save_filename}.tif"
        save_path.parent.mkdir(exist_ok=True)
    return save_path, rvt_save_name


def save_path_float(save_path):
    """Adds suffix _float to the existing name"""
    save_path = Path(save_path)

    parent_folder = Path(save_path.parent.as_posix() + "_float")
    file_name = save_path.stem + "_float" + save_path.suffix

    spf = parent_folder / file_name

    return spf


def create_mosaic(input_files_list, output_file):
    # Get a list of all GeoTIFF files in the directory
    # input_files_list = glob.glob(os.path.join(input_dir, '*.tif'))

    # Open the GeoTIFF files
    src_files_to_mosaic = [rasterio.open(fp) for fp in input_files_list]
    # for fp in input_files_list:
    #     src = rasterio.open(fp)
    #     src_files_to_mosaic.append(src)

    # Merge the GeoTIFF files
    mosaic, out_trans = merge(src_files_to_mosaic)

    # Copy the metadata from one of the input files
    out_meta = src_files_to_mosaic[0].meta.copy()
    out_crs = src_files_to_mosaic[0].crs

    # Update the metadata with the new dimensions, transform, and CRS
    out_meta.update({
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_trans,
        "crs": out_crs,
        "compress": "LZW",
        "predictor": 2,
        'BIGTIFF': 'IF_SAFER'
    })

    # Write the mosaic to a new GeoTIFF file
    with rasterio.open(output_file, "w", **out_meta) as dest:
        dest.write(mosaic)

    # Close the source files
    for src in src_files_to_mosaic:
        src.close()

    # # Delete tiles (temp files)
    # shutil.rmtree(input_files_list[0].parent)

    # print(f"Mosaic created and saved as {output_file}")

    return output_file


def vrt_to_mosaic(vrt_path, output_mosaic_path):

    with rasterio.open(vrt_path) as vrt:
        mosaic = vrt.read()
        out_meta = vrt.profile.copy()

    # Update metadata and store to file
    out_meta.update({
        'driver': 'GTiff',
        "compress": "LZW",
        "predictor": 2,
        'tiled': True,
        'blockxsize': 256,
        'blockysize': 256,
        'BIGTIFF': 'IF_SAFER'
    })

    with rasterio.open(output_mosaic_path, "w", **out_meta) as dest:
        dest.write(mosaic)

    # print(f"Mosaic saved as {output_mosaic_path}")
