"""
Tiled processing with RVT_py
Created on 6 May 2024
@author: Nejc Čož, ZRC SAZU, Novi trg 2, 1000 Ljubljana, Slovenia
"""

import glob
import multiprocessing as mp
import time
from math import ceil
from pathlib import Path

import numpy as np
import rasterio
import rvt.blend
import rvt.default
import rvt.vis
from osgeo import gdal
from rasterio.windows import from_bounds
from rvt.blend_func import normalize_image


def tiled_blending(blend_types, input_vrt_path, tiles_list, nr_processes):
    t0 = time.time()

    # ds_path = Path(input_vrt_path).parent
    # # ds_name = ds_path.name
    # # We are saving low-levels to a different disk, but we are keeping the folder structure the same
    # ll_path = ds_path  # Path(low_levels_main_dir) / ds_name

    src_tif_path = Path(input_vrt_path)
    ll_path = src_tif_path.parent

    # We need to know the resolution, to select correct radius R# for svf/opns
    with rasterio.open(src_tif_path) as src:
        res = src.res[0]

    # # HERE MULTIPROCESSING STARTS
    #
    # - ds_path  (const) ... for creating output file name
    # - ll_path  (const)
    # - blend_types  (const)
    # - res  (const)
    # - one_tile  (variable)
    #
    # ----------------------------------
    input_process_list = [(src_tif_path, ll_path, blend_types, res, i) for i in tiles_list]
    with mp.Pool(nr_processes) as p:
        realist = [p.apply_async(compute_save_blends, r) for r in input_process_list]
        for result in realist:
            pool_out = result.get()
            print("Finished tile:", pool_out[1])

    # # SINGLE-PROCESS FOR DEBUG
    # for one_tile in [tiles_list[3]]:
    #     result = compute_save_blends(src_tif_path, ll_path, blend_types, res, one_tile)
    #     print("Finished tile:", result[1])

    # Build VRTs
    for ds_dir in [ll_path / i for i in blend_types]:
        # Name = <original DEM name> + <visualization (subdir name)> + .vrt
        if ds_dir == "e2MSTP" or ds_dir == "e3MSTP":
            vrt_name = Path(input_vrt_path).stem + "_" + Path(ds_dir).name + ".vrt"
            ds_dir = ds_dir
        else:
            vrt_name = Path(input_vrt_path).stem + "_" + Path(ds_dir).name + ".vrt"
        out_path = build_vrt(ds_dir, vrt_name)
        print("  - Created:", out_path)
    if "e2MSTP" in blend_types:
        ds_dir = ll_path / "rrim"
        vrt_name = Path(input_vrt_path).stem + "_" + Path(ds_dir).name + ".vrt"
        vrt_path = ds_dir.parents[0] / vrt_name
        if not vrt_path.exists():
            out_path = build_vrt(ds_dir, vrt_name)
            print("  - Created:", out_path)
    if "e3MSTP" in blend_types:
        ds_dir = ll_path / "crim"
        vrt_name = Path(input_vrt_path).stem + "_" + Path(ds_dir).name + ".vrt"
        vrt_path = ds_dir.parents[0] / vrt_name
        if not vrt_path.exists():
            out_path = build_vrt(ds_dir, vrt_name)
            print("  - Created:", out_path)

    t1 = time.time() - t0
    print(f"Done with computing blends in {round(t1/60, ndigits=None)} min.")


def compute_save_blends(src_path, low_levels_path, blend_types, res, one_extent):

    # Determine name of the tile (coordinates)
    one_tile = f"{one_extent[0]:.0f}_{one_extent[1]:.0f}"

    # Determine required datasets paths!
    req_arrays, req_arr_dirs = get_required_arrays(blend_types, low_levels_path, res, one_tile, find_paths=False)

    # ***************************************************************
    # Default 1 is for GENERAL
    default_1 = rvt.default.DefaultValues()
    # Read from file, path is relative to the Current script directory
    def1_pth = Path(__file__).resolve().parent / "default_1.json"
    default_1.read_default_from_file(def1_pth)

    default_1.fill_no_data = 1
    default_1.keep_original_no_data = 0

    # default_1.hs_sun_el = 35  # HILLSHADE FOR VAT General
    # default_1.svf_r_max = ceil(5 / res)

    # default_1.slrm_rad_cell = ceil(10 / res) if res < 1 else 10

    # default_1.ld_min_rad = ceil(10 / res)
    # default_1.ld_max_rad = ceil(20 / res)

    # MSTP (default values divided by resolution to get meters)
    # default_1.mstp_local_scale = tuple(ceil(ti / res) for ti in default_1.mstp_local_scale)
    # default_1.mstp_meso_scale = tuple(ceil(ti / res) for ti in default_1.mstp_meso_scale)
    # default_1.mstp_broad_scale = tuple(ceil(ti / res) for ti in default_1.mstp_broad_scale)

    # Default 2 is for FLAT
    default_2 = rvt.default.DefaultValues()
    # Read from file, path is relative to the Current script directory
    def2_pth = Path(__file__).resolve().parent / "default_2.json"
    default_2.read_default_from_file(def2_pth)

    default_2.fill_no_data = 1
    default_2.keep_original_no_data = 0

    # default_2.hs_sun_el = 15
    # default_2.svf_r_max = ceil(10 / res)  # 10 m (divide by pixel size)
    # ***************************************************************

    # Only compute required visualizations
    in_arrays = compute_low_levels(
        default_1,
        default_2,
        src_path,
        one_extent,
        req_arrays
    )

    if "vat_combined_8bit" in blend_types:
        # Determine save path
        save_path = low_levels_path / "VAT_combined_8bit" / f"{one_tile}_rvt_VAT_comb_8bit.tif"
        save_path.parent.mkdir(exist_ok=True)
        in_arrays["vat_combined_8bit"] = vat_combined_8bit(in_arrays, save_path)

    if "VAT_flat_3B" in blend_types:
        # Determine save path
        save_path = low_levels_path / "VAT_flat_3B" / f"{one_tile}_rvt_VAT_flat_3B.tif"
        save_path.parent.mkdir(exist_ok=True)
        in_arrays["vat_flat_3bands"] = vat_flat_3bands(in_arrays, save_path)

    if "VAT_3B" in blend_types:
        # Determine save path
        save_path = low_levels_path / "VAT_3B" / f"{one_tile}_rvt_VAT_3B.tif"
        save_path.parent.mkdir(exist_ok=True)
        in_arrays["vat_3bands"] = vat_3bands(in_arrays, save_path)

    if "VAT_combined_3B" in blend_types:
        # Determine save path
        save_path = low_levels_path / "VAT_combined_3B" / f"{one_tile}_rvt_VAT_combined_3B.tif"
        save_path.parent.mkdir(exist_ok=True)
        vat_combined_3bands(in_arrays, save_path)

    if "SLRM" in blend_types:
        # Determine save path
        save_path = low_levels_path / "SLRM" / f"{one_tile}_rvt_SLRM.tif"
        save_path.parent.mkdir(exist_ok=True)
        slrm_normalize(in_arrays, save_path)

    if "e2MSTP" in blend_types:
        if not req_arrays["rrim"]:
            save_path = low_levels_path / "rrim" / f"{one_tile}_rvt_RRIM.tif"
            save_path.parent.mkdir(exist_ok=True)
            in_arrays["rrim"] = blend_rrim(in_arrays, save_path)

        # Determine save path
        save_path = low_levels_path / "e2MSTP" / f"{one_tile}_rvt_e2MSTP.tif"
        save_path.parent.mkdir(exist_ok=True)
        blend_e2mstp(in_arrays, save_path)

    if "e3MSTP" in blend_types:
        if not req_arrays["crim"]:
            save_path = low_levels_path / "crim" / f"{one_tile}_rvt_CRIM.tif"
            save_path.parent.mkdir(exist_ok=True)
            in_arrays["crim"] = blend_crim(in_arrays, save_path)

        # Determine save path
        save_path = low_levels_path / "e3MSTP" / f"{one_tile}_rvt_e3MSTP.tif"
        save_path.parent.mkdir(exist_ok=True)
        blend_e3mstp(in_arrays, save_path)

    if "e4MSTP" in blend_types:
        # Determine save path
        save_path = low_levels_path / "e4MSTP" / f"{one_tile}_rvt_e4MSTP.tif"
        save_path.parent.mkdir(exist_ok=True)
        blend_e4mstp(in_arrays, save_path)

    return 0, one_tile


def vat_combined_8bit(dict_arrays, save_path):
    """
    VAT Combined 8bit
    """

    # BLEND VAT GENERAL
    vat_combination_general = rvt.blend.BlenderCombination()
    vat_combination_general.create_layer(
        vis_method="Sky-View Factor",
        normalization="Value",
        minimum=0.7,
        maximum=1.0,
        blend_mode="Multiply",
        opacity=25,
        image=dict_arrays['svf_R_small'].squeeze()
    )
    vat_combination_general.create_layer(
        vis_method="Openness - Positive",
        normalization="Value",
        minimum=68,
        maximum=93,
        blend_mode="Overlay",
        opacity=50,
        image=dict_arrays['opns_R_small'].squeeze()
    )
    vat_combination_general.create_layer(
        vis_method="Slope gradient",
        normalization="Value",
        minimum=0,
        maximum=50,
        blend_mode="Luminosity",
        opacity=50,
        image=dict_arrays['slope'].squeeze()
    )
    vat_combination_general.create_layer(
        vis_method="Hillshade",
        normalization="Value",
        minimum=0,
        maximum=1,
        blend_mode="Normal",
        opacity=100,
        image=dict_arrays['hs_general'].squeeze()
    )
    vat_1 = vat_combination_general.render_all_images(
        save_visualizations=False,
        save_render_path=None,
        no_data=np.nan
    )
    vat_1 = vat_1.astype("float32")

    # BLEND VAT FLAT
    vat_combination_flat = rvt.blend.BlenderCombination()
    vat_combination_flat.create_layer(
        vis_method="Sky-View Factor",
        normalization="Value",
        minimum=0.9,
        maximum=1.0,
        blend_mode="Multiply",
        opacity=25,
        image=dict_arrays['svf_R_large'].squeeze()
    )
    vat_combination_flat.create_layer(
        vis_method="Openness - Positive",
        normalization="Value",
        minimum=85,
        maximum=93,
        blend_mode="Overlay",
        opacity=50,
        image=dict_arrays['opns_R_large'].squeeze()
    )
    vat_combination_flat.create_layer(
        vis_method="Slope gradient",
        normalization="Value",
        minimum=0,
        maximum=15,
        blend_mode="Luminosity",
        opacity=50,
        image=dict_arrays['slope'].squeeze()
    )
    vat_combination_flat.create_layer(
        vis_method="Hillshade",
        normalization="Value",
        minimum=0,
        maximum=1,
        blend_mode="Normal",
        opacity=100,
        image=dict_arrays['hs_flat'].squeeze()
    )
    vat_2 = vat_combination_flat.render_all_images(
        save_visualizations=False,
        save_render_path=None,
        no_data=np.nan
    )
    vat_2 = vat_2.astype("float32")

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

    out_vat_combined = rvt.vis.byte_scale(out_vat_combined, c_min=0, c_max=1)

    # Save GeoTIF
    out_profile = dict_arrays['profile'].copy()
    out_profile.update(dtype='uint8', nodata=None)
    rasterio_save(
        out_vat_combined,
        out_profile,
        save_path=save_path
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
        image=dict_arrays["svf_R_large"].squeeze(),
        min_norm=0.9,
        max_norm=1,
        normalization="value"
    )
    opns = normalize_image(
        visualization="openness - positive",
        image=dict_arrays["opns_R_large"].squeeze(),
        min_norm=85,
        max_norm=93,
        normalization="value"
    )
    slope = normalize_image(
        visualization="slope gradient",
        image=dict_arrays["slope"].squeeze(),
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
        image=dict_arrays["svf_R_small"].squeeze(),
        min_norm=0.7,
        max_norm=1,
        normalization="value"
    )
    opns = normalize_image(
        visualization="openness - positive",
        image=dict_arrays["opns_R_small"].squeeze(),
        min_norm=68,
        max_norm=93,
        normalization="value"
    )
    slope = normalize_image(
        visualization="slope gradient",
        image=dict_arrays["slope"].squeeze(),
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


def slrm_normalize(dict_arrays, save_path):
    out_slrm = normalize_image(
        visualization="slrm",
        image=dict_arrays["slrm_vis"].squeeze(),
        min_norm=-0.5,
        max_norm=0.5,
        normalization="value"
    )
    # Save GeoTIF
    rasterio_save(
        out_slrm,
        dict_arrays['profile'],
        save_path=save_path,
        nodata=np.nan
    )
    return out_slrm


def blend_rrim(dict_arrays, save_path):
    comb_rrim = rvt.blend.BlenderCombination()
    comb_rrim.create_layer(vis_method="Slope gradient", normalization="Value",
                           minimum=0, maximum=45,
                           blend_mode="Normal", opacity=50,
                           colormap="Reds_r", min_colormap_cut=0, max_colormap_cut=1,
                           image=dict_arrays['slope'].squeeze()
                           )
    comb_rrim.create_layer(vis_method="Opns_Pos_Neg/2", normalization="Value",
                           minimum=-25, maximum=25,
                           blend_mode="Normal", opacity=100,
                           colormap="Greys_r", min_colormap_cut=0, max_colormap_cut=1,
                           image=((dict_arrays['opns_R_small'] - dict_arrays['neg_opns_R_small'])/2).squeeze()
                           )
    out_rrim = comb_rrim.render_all_images(save_visualizations=False,
                                           save_render_path=None,
                                           no_data=np.nan)
    out_rrim = out_rrim.astype("float32")
    # Save GeoTIF
    rasterio_save(
        out_rrim,
        dict_arrays['profile'],
        save_path=save_path,
        nodata=np.nan
    )

    return out_rrim


def blend_e2mstp(dict_arrays, save_path):
    comb_e2mstp = rvt.blend.BlenderCombination()
    comb_e2mstp.create_layer(vis_method="slrm", normalization="value",
                             minimum=-0.5, maximum=0.5,
                             blend_mode="screen", opacity=25,
                             image=dict_arrays["slrm_vis"].squeeze()
                             )
    comb_e2mstp.create_layer(vis_method="rrim", normalization="value",
                             minimum=0, maximum=1,
                             blend_mode="soft_light", opacity=70,
                             image=dict_arrays["rrim"]
                             )
    comb_e2mstp.create_layer(vis_method="mstp", normalization="value",
                             minimum=0, maximum=1,
                             blend_mode="normal", opacity=100,
                             image=dict_arrays["mstp"]
                             )
    out_e2mstp = comb_e2mstp.render_all_images(save_visualizations=False,
                                               save_render_path=None,
                                               no_data=np.nan)
    out_e2mstp = out_e2mstp.astype("float32")
    out_e2mstp[np.isnan(dict_arrays["mstp"])] = np.nan
    out_e2mstp[out_e2mstp > 1] = 1
    # Save GeoTIF
    rasterio_save(
        out_e2mstp,
        dict_arrays['profile'],
        save_path=save_path,
        nodata=np.nan
    )

    return out_e2mstp


def blend_crim(dict_arrays, save_path):
    comb_crim = rvt.blend.BlenderCombination()
    comb_crim.create_layer(vis_method="Openness_Pos-Neg", normalization="Value",
                           minimum=-28, maximum=28,
                           blend_mode="overlay", opacity=50,
                           image=(dict_arrays['opns_R_small'] - dict_arrays['neg_opns_R_small']).squeeze()
                           )
    comb_crim.create_layer(vis_method="Openness_Pos-Neg", normalization="Value",
                           minimum=-28, maximum=28,
                           blend_mode="luminosity", opacity=50,
                           image=(dict_arrays['opns_R_small'] - dict_arrays['neg_opns_R_small']).squeeze()
                           )
    comb_crim.create_layer(vis_method="slope gradient red", normalization="Value",
                           minimum=0, maximum=45,
                           blend_mode="normal", opacity=100,
                           colormap="OrRd", min_colormap_cut=0, max_colormap_cut=1,
                           image=dict_arrays['slope'].squeeze()
                           )
    out_crim = comb_crim.render_all_images(save_visualizations=False,
                                           save_render_path=None,
                                           no_data=np.nan)
    out_crim = out_crim.astype("float32")
    # Save GeoTIF
    rasterio_save(
        out_crim,
        dict_arrays['profile'],
        save_path=save_path,
        nodata=np.nan
    )

    return out_crim


def blend_e3mstp(dict_arrays, save_path):
    comb_e3mstp = rvt.blend.BlenderCombination()
    comb_e3mstp.create_layer(vis_method="slrm", normalization="value",
                             minimum=-0.5, maximum=0.5,
                             blend_mode="screen", opacity=25,
                             image=dict_arrays["slrm_vis"].squeeze()
                             )
    comb_e3mstp.create_layer(vis_method="crim", normalization="value",
                             minimum=0, maximum=1,
                             blend_mode="soft_light", opacity=70,
                             image=dict_arrays["crim"]
                             )
    comb_e3mstp.create_layer(vis_method="mstp", normalization="value",
                             minimum=0, maximum=1,
                             blend_mode="normal", opacity=100,
                             image=dict_arrays["mstp"]
                             )
    out_e3mstp = comb_e3mstp.render_all_images(save_visualizations=False,
                                               save_render_path=None,
                                               no_data=np.nan)
    out_e3mstp = out_e3mstp.astype("float32")
    out_e3mstp[np.isnan(dict_arrays["mstp"])] = np.nan
    out_e3mstp[out_e3mstp > 1] = 1
    # Save GeoTIF
    rasterio_save(
        out_e3mstp,
        dict_arrays['profile'],
        save_path=save_path,
        nodata=np.nan
    )


def blend_e4mstp(dict_arrays, save_path):
    # Get Coloured Slope
    dict_arrays['cs'] = blend_coloured_slope(dict_arrays)
    # Get SVF combined
    dict_arrays['svf_combined'] = blend_svf_combined(dict_arrays)
    # Get Openness + LD
    dict_arrays['opns_ld'] = blend_opns_ld(dict_arrays)

    comb_nv = rvt.blend.BlenderCombination()
    comb_nv.create_layer(
        vis_method="mstp",
        normalization="value", minimum=0, maximum=1,
        blend_mode="overlay", opacity=90,
        image=dict_arrays['mstp']
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
        vis_method="coloured slope",
        normalization="value", minimum=0, maximum=1,
        blend_mode="normal", opacity=100,
        image=dict_arrays['cs']
    )
    out_e4mstp = comb_nv.render_all_images(
        save_visualizations=False,
        save_render_path=None,
        no_data=np.nan
    )
    out_e4mstp = out_e4mstp.astype("float32")
    out_e4mstp[np.isnan(dict_arrays['mstp'])] = np.nan
    out_e4mstp[out_e4mstp > 1] = 1

    # Save GeoTIF
    rasterio_save(
        out_e4mstp,
        dict_arrays['profile'],
        save_path=save_path,
        nodata=np.nan
    )

    return out_e4mstp


def blend_coloured_slope(dict_arrays, save_path=None):
    # Coloured slope
    comb_cs = rvt.blend.BlenderCombination()
    comb_cs.create_layer(
        vis_method="Slope gradient", normalization="Value",
        minimum=0, maximum=55,
        blend_mode="normal", opacity=100,
        colormap="Reds_r", min_colormap_cut=0, max_colormap_cut=1,
        image=dict_arrays['slope'].squeeze()
    )
    coloured_slope = comb_cs.render_all_images(
        save_visualizations=False,
        save_render_path=None,
        no_data=np.nan
    )
    coloured_slope = coloured_slope.astype("float32")

    if save_path:
        # Save GeoTIF
        rasterio_save(
            coloured_slope,
            dict_arrays['profile'],
            save_path=save_path,
            nodata=np.nan
        )

    return coloured_slope


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
        image=(dict_arrays['opns_R_small'] - dict_arrays['neg_opns_R_small']).squeeze()
    )
    comb.create_layer(
        vis_method="Local dominance", normalization="Value",
        minimum=0.5, maximum=1.8,
        blend_mode="normal", opacity=100,
        image=dict_arrays['ld'].squeeze()
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
                          image=dict_arrays['svf_R_small'].squeeze()
                          )
    comb_svf.create_layer(vis_method="Sky-view factor", normalization="Value",
                          minimum=0.9, maximum=1,
                          blend_mode="normal", opacity=100,
                          image=dict_arrays['svf_R_large'].squeeze()
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
    with rasterio.open(save_path, "w", **profile) as dst:
        dst.write(array)


def get_required_arrays(blend_types, low_levels_path, res, one_tile, find_paths=True):
    r_large = ceil(10 / res) if res < 1 else 10
    r_small = ceil(5 / res)

    req_arrays = {
        # These are visualizations:
        "slope": False,
        "slrm_vis": False,
        "ld": False,
        "mstp": False,

        # Hillshade general/flat
        "hs_general": False,
        "hs_flat": False,

        # small = GENERAL = 5m
        "svf_R_small": False,
        "opns_R_small": False,
        "neg_opns_R_small": False,

        # large = FLAT = 10m
        "svf_R_large": False,
        "opns_R_large": False,
        "neg_opns_R_large": False,

        # These are combinations:
        "rrim": False,
        "crim": False
    }

    if "VAT_3B" in blend_types:
        req_arrays["svf_R_small"] = True
        req_arrays["opns_R_small"] = True
        req_arrays["slope"] = True

    if "VAT_flat_3B" in blend_types:
        req_arrays["svf_R_large"] = True
        req_arrays["opns_R_large"] = True
        req_arrays["slope"] = True

    if "VAT_combined_3B" in blend_types:
        req_arrays["svf_R_large"] = True
        req_arrays["opns_R_large"] = True
        req_arrays["svf_R_small"] = True
        req_arrays["opns_R_small"] = True
        req_arrays["slope"] = True

    if "SLRM" in blend_types:
        req_arrays["slrm_vis"] = True

    if ("e3MSTP" in blend_types) or ("e2MSTP" in blend_types):
        req_arrays["slrm_vis"] = True
        req_arrays["mstp"] = True
        req_arrays["slope"] = True
        req_arrays["opns_R_small"] = True
        req_arrays["neg_opns_R_small"] = True

    if "e4MSTP" in blend_types:
        req_arrays["ld"] = True
        req_arrays["svf_R_small"] = True,
        req_arrays["mstp"] = True
        req_arrays["svf_R_large"] = True
        req_arrays["opns_R_large"] = True
        req_arrays["opns_R_small"] = True
        req_arrays["neg_opns_R_small"] = True
        req_arrays["slope"] = True

    if "vat_combined_8bit" in blend_types:
        req_arrays["svf_R_large"] = True
        req_arrays["opns_R_large"] = True
        req_arrays["svf_R_small"] = True
        req_arrays["opns_R_small"] = True
        req_arrays["slope"] = True
        req_arrays["hs_general"] = True
        req_arrays["hs_flat"] = True

    if find_paths:
        req_arr_dirs = {}
        for (k, v) in req_arrays.items():
            vis_type = k
            if k[-6:] == "_small":
                k = k.replace("_small", f"{r_small}")
            elif k[-6:] == "_large":
                k = k.replace("_large", f"{r_large}")

            if v:
                add_path = next(Path(low_levels_path / k).glob(f"{one_tile}*.tif"))
                req_arr_dirs[vis_type] = add_path
    else:
        req_arr_dirs = None

    return req_arrays, req_arr_dirs


def compute_low_levels(
        default_1,
        default_2,
        vrt_path,
        input_dem_extents,
        vis_types
):
    # Read buffer values from defaults (already changed from meters to pixels!!!)!
    all_buffers = {
        "slope": 0,
        "slrm_vis": default_1.slrm_rad_cell,
        "ld": default_1.ld_max_rad,
        # "mstp": default_1.mstp_broad_scale[1],

        "hs_general": 1,
        "hs_flat": 1,

        "svf_GEN": default_1.svf_r_max,  # SVF, OPNS+ and OPNS- for (GENERAL = SMALL = 5m)

        "svf_FLAT": default_2.svf_r_max,  # SVF, OPNS+ and OPNS- for (FLAT = LARGE = 10m)
    }

    # Filter buffer_dict based on required visualizations
    req_visualizations = [a for (a, v) in vis_types.items() if v]
    if any(value in req_visualizations for value in ["svf_R_small", "opns_R_small", "neg_opns_R_small"]):
        req_visualizations.append("svf_GEN")
    if any(value in req_visualizations for value in ["svf_R_large", "opns_R_large", "neg_opns_R_large"]):
        req_visualizations.append("svf_FLAT")
    buffer_dict = {key: all_buffers[key] for key in req_visualizations if key in all_buffers}

    # Select the largest required buffer
    max_buff = max(buffer_dict, key=buffer_dict.get)
    buffer = buffer_dict[max_buff]

    # Read array into RVT dictionary format
    dict_arrays = get_raster_vrt(vrt_path, input_dem_extents, buffer)

    # Change nodata value to np.nan, to avoid problems later
    dict_arrays["array"][dict_arrays["array"] == dict_arrays["no_data"]] = np.nan
    dict_arrays["no_data"] = np.nan

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
        if vis_type == "slope":
            vis_out = {
                vis_type: default_1.get_slope(
                    sliced_arr,
                    resolution_x=dict_arrays["resolution"][0],
                    resolution_y=dict_arrays["resolution"][1]
                )
            }
        elif vis_type == "slrm_vis":
            vis_out = {
                vis_type: default_1.get_slrm(sliced_arr)
            }
        elif vis_type == "ld":
            vis_out = {
                vis_type: default_1.get_local_dominance(sliced_arr)
            }
        elif vis_type == "mstp":
            vis_out = {
                vis_type: default_1.get_mstp(sliced_arr)
            }
        elif vis_type == "hs_general":
            vis_out = {
                vis_type: default_1.get_hillshade(
                    sliced_arr,
                    resolution_x=dict_arrays["resolution"][0],
                    resolution_y=dict_arrays["resolution"][1]
                )
            }
        elif vis_type == "hs_flat":
            vis_out = {
                vis_type: default_2.get_hillshade(
                    sliced_arr,
                    resolution_x=dict_arrays["resolution"][0],
                    resolution_y=dict_arrays["resolution"][1]
                )
            }
        elif vis_type == "svf_GEN":  # small, GENERAL, 5m
            # Check which of the 3 to be computed
            compute_svf = True if "svf_R_small" in req_visualizations else False
            compute_opns = True if "opns_R_small" in req_visualizations else False
            compute_neg_opns = True if "neg_opns_R_small" in req_visualizations else False

            if compute_svf or compute_opns:
                vis_out = default_1.get_sky_view_factor(
                    sliced_arr,
                    dict_arrays["resolution"][0],
                    compute_svf=compute_svf,
                    compute_opns=compute_opns
                )
                # Rename to correct vis_type (which svf is this?)
                for k in list(vis_out.keys()):
                    vis_out[f"{k}_R_small"] = vis_out.pop(k)

            if compute_neg_opns:
                vis_out = {
                    "neg_opns_R_small": default_1.get_neg_opns(
                        sliced_arr,
                        dict_arrays["resolution"][0]
                    )
                }
        elif vis_type == "svf_FLAT":  # large, FLAT, 10m
            # Check which of the 3 to be computed
            compute_svf = True if "svf_R_large" in req_visualizations else False
            compute_opns = True if "opns_R_large" in req_visualizations else False
            compute_neg_opns = True if "neg_opns_R_large" in req_visualizations else False

            if compute_svf or compute_opns:
                vis_out = default_1.get_sky_view_factor(
                    sliced_arr,
                    dict_arrays["resolution"][0],
                    compute_svf=compute_svf,
                    compute_opns=compute_opns
                )
                # Rename to correct vis_type (which svf is this?)
                for k in list(vis_out.keys()):
                    vis_out[f"{k}_R_large"] = vis_out.pop(k)

            if compute_neg_opns:
                vis_out = {
                    "neg_opns_R_large": default_1.get_neg_opns(
                        sliced_arr,
                        dict_arrays["resolution"][0]
                    )
                }
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
        # Read VRT metadata
        vrt_res = vrt.res
        vrt_nodata = vrt.nodata
        vrt_transform = vrt.transform
        vrt_crs = vrt.crs

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

        # Save transform object of both extents (original and buffered)
        buff_transform = vrt.window_transform(buff_window)
        orig_transform = vrt.window_transform(orig_window)

    # For raster with only one band, remove first axis from the array (RVT requirement)
    if win_array.shape[0] == 1:
        win_array = np.squeeze(win_array, axis=0)

    # Prepare output metadata profile
    out_profile = {
        'driver': 'GTiff',
        'nodata': None,
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
        "buff_transform": buff_transform,
        "orig_transform": orig_transform,
        "crs": vrt_crs,
        "profile": out_profile
    }

    return output


def build_vrt(ds_dir, vrt_name):
    ds_dir = Path(ds_dir)
    vrt_path = ds_dir.parents[0] / vrt_name
    tif_list = glob.glob(Path(ds_dir / "*.tif").as_posix())

    vrt_options = gdal.BuildVRTOptions()
    my_vrt = gdal.BuildVRT(vrt_path.as_posix(), tif_list, options=vrt_options)
    my_vrt = None

    return vrt_path
