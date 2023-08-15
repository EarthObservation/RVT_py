"""
Relief Visualization Toolbox – Visualization Functions

Contains functions to save visualizations tile by tile.

Credits:
    Žiga Kokalj (ziga.kokalj@zrc-sazu.si)
    Krištof Oštir (kristof.ostir@fgg.uni-lj.si)
    Klemen Zakšek
    Peter Pehani
    Klemen Čotar
    Maja Somrak
    Žiga Maroh
    Nejc Čož

Copyright:
    2010-2022 Research Centre of the Slovenian Academy of Sciences and Arts
    2016-2022 University of Ljubljana, Faculty of Civil and Geodetic Engineering
"""
from pathlib import Path
from typing import Callable, Dict, Any, Optional
import numpy as np
import numpy.typing as npt
from osgeo import gdal
import rvt.default


def _create_blank_raster(
    in_data_set: gdal.Dataset,
    out_raster_path: Path,
    nr_bands: int = 1,
    no_data: float = np.nan,
    e_type: int = 6,
) -> None:
    """Takes input data set and creates new raster. It copies input data set size, projection and geo info."""
    gtiff_driver = gdal.GetDriverByName("GTiff")
    band = in_data_set.GetRasterBand(1)
    x_size = band.XSize  # number of columns
    y_size = band.YSize  # number of rows
    out_ds = gtiff_driver.Create(
        out_raster_path.as_posix(),
        xsize=x_size,
        ysize=y_size,
        bands=nr_bands,
        eType=e_type,
        options=["BIGTIFF=IF_NEEDED"],
    )
    out_ds.SetProjection(in_data_set.GetProjection())
    out_ds.SetGeoTransform(in_data_set.GetGeoTransform())
    out_ds.GetRasterBand(1).SetNoDataValue(no_data)
    out_ds.FlushCache()
    out_ds = None


def save_visualization_tile_by_tile(
    visualization_function: Callable,
    function_parameters: Optional[Dict[str, Optional[Any]]],
    dem_path: Path,
    overlap: int,
    tile_size_x: int,
    tile_size_y: int,
    out_raster_path: Path,
    out_raster_nr_of_bands: int = 1,
    out_raster_e_type: int = 6,
    out_visualization_dict_key: Optional[str] = None,
) -> None:
    """
    Some DEMs are too large to load them into memory. This function reads dem raster tile by tile,
    calculates visualization on it tile by tile and than saves calculated visualization tile by tile in out raster.
    Note that visualization_function needs dem parameter but it shouldn't be inputted in function_parameters because it
    is read tile_by_tile from dem_path.

    Parameters
    ----------
    visualization_function : Callable
        Python function which represents visualization function. Function needs to have parameter called dem!
    function_parameters: Optional[Dict[str, Optional[Any]]]
        Visualization function parameters in form of dict where key represents parameter and key value parameter value.
        Parameter dem needs to be excluded because it is read from dem_path!
        If function_parameters contains no_data key and its value is None it will be taken from dem (dem_path),
        same goes for resolutions (resolution_x, resolution_y, resolution).
    dem_path : Path
        Path to a Digital elevation model.
    overlap : int
        When calculating visualization on tile we need some information from neighbouring tiles.
        This parameter defines number of pixels we need from neighbouring tiles (overlap, offset).
    tile_size_x : int
        Tile size in pixels in x direction (number of columns).
    tile_size_y : int
        Tile size in pixels in y direction (number of rows).
    out_raster_path : Path
        Path to output visualization.
    out_raster_nr_of_bands : int
        Output visualization number of bands.
    out_raster_e_type : GDALDataType
        Output visualization data type.
        https://gdal.org/api/raster_c_api.html#_CPPv412GDALDataType, (GDT_Float32 = 6, GDT_UInt8 = 1, ...)
    out_visualization_dict_key : Optional[str]
        Set to None if output of visualization is 2D numpy array.
        If output of visualization function is dictionary then this parameter is key,
        to define result 2D numpy array in dictionary.
        For example rvt.visualization.slope_aspect outputs dictionary with keys "slope" and "aspect".
        To select slope set this parameter to "slope".

    Returns
    -------
    out : None
    """
    if not dem_path.exists():
        Exception(
            "rvt.tile.save_visualization_tile_by_tile: Input dem path does not exist!"
        )
    if tile_size_x < 50 or tile_size_y < 50:
        Exception(
            "rvt.tile.save_visualization_tile_by_tile: Tile size too small (tile_size_x, tile_size_y),"
            " it needs to be bigger than 50 pixels!"
        )

    dem_ds = gdal.Open(dem_path.as_posix())
    gt = dem_ds.GetGeoTransform()
    x_res = gt[1]  # x_resolution
    y_res = -gt[5]  # y_resolution
    no_data = dem_ds.GetRasterBand(1).GetNoDataValue()
    band = dem_ds.GetRasterBand(1)
    x_size = band.XSize  # number of columns
    y_size = band.YSize  # number of rows

    # set resolution and no_data function_parameters if needed (are set to None) from dem
    if function_parameters is not None:
        if "resolution" in function_parameters:
            if function_parameters["resolution"] is None:
                function_parameters["resolution"] = x_res
        if "resolution_x" in function_parameters:
            if function_parameters["resolution_x"] is None:
                function_parameters["resolution_x"] = x_res
        if "resolution_y" in function_parameters:
            if function_parameters["resolution_y"] is None:
                function_parameters["resolution_y"] = y_res
        if "no_data" in function_parameters:
            if function_parameters["no_data"] is None:
                function_parameters["no_data"] = no_data

    _create_blank_raster(
        in_data_set=dem_ds,
        out_raster_path=out_raster_path,
        nr_bands=out_raster_nr_of_bands,
        e_type=out_raster_e_type,
    )

    for y in range(0, y_size, tile_size_y):
        if y + tile_size_y < y_size:  # if rows overlap
            rows = tile_size_y
        else:
            rows = y_size - y
        for x in range(0, x_size, tile_size_x):
            if x + tile_size_x < x_size:  # if cols overlap
                cols = tile_size_x
            else:
                cols = x_size - x

            # get offset for each tile, check edges
            left_offset = 0
            right_offset = 0
            top_offset = 0
            bottom_offset = 0
            if x != 0:  # left offset
                if x - overlap < 0:  # left overlap
                    left_offset = x
                else:
                    left_offset = overlap
            if x + cols != x_size:  # right offset
                if x + cols + overlap > x_size:  # right overlap
                    right_offset = x_size - x - cols
                else:
                    right_offset = overlap
            if y != 0:  # apply top offset
                if y - overlap < 0:  # top overlap
                    top_offset = y
                else:
                    top_offset = overlap
            if y + rows != y_size:  # bottom offset
                if y + rows + overlap > y_size:  # bottom overlap
                    bottom_offset = y_size - y - rows
                else:
                    bottom_offset = overlap

            # read tile
            x_off = x - left_offset
            y_off = y - top_offset
            cols_off = cols + left_offset + right_offset
            rows_off = rows + top_offset + bottom_offset

            tile_array = np.array(
                dem_ds.GetRasterBand(1).ReadAsArray(x_off, y_off, cols_off, rows_off)
            )
            if function_parameters is not None:
                visualization_array = visualization_function(
                    dem=tile_array, **function_parameters
                )
            else:
                visualization_array = visualization_function(dem=tile_array)

            if out_visualization_dict_key is not None:
                visualization_array = visualization_array[out_visualization_dict_key]

            # remove offset from visualization block
            if out_raster_nr_of_bands == 1:
                if right_offset == 0 and bottom_offset == 0:
                    visualization_array = visualization_array[top_offset:, left_offset:]
                elif right_offset == 0:
                    visualization_array = visualization_array[
                        top_offset:-bottom_offset, left_offset:
                    ]
                elif bottom_offset == 0:
                    visualization_array = visualization_array[
                        top_offset:, left_offset:-right_offset
                    ]
                else:
                    visualization_array = visualization_array[
                        top_offset:-bottom_offset, left_offset:-right_offset
                    ]
            else:
                if right_offset == 0 and bottom_offset == 0:
                    visualization_array = visualization_array[
                        :, top_offset:, left_offset:
                    ]
                elif right_offset == 0:
                    visualization_array = visualization_array[
                        :, top_offset:-bottom_offset, left_offset:
                    ]
                elif bottom_offset == 0:
                    visualization_array = visualization_array[
                        :, top_offset:, left_offset:-right_offset
                    ]
                else:
                    visualization_array = visualization_array[
                        :, top_offset:-bottom_offset, left_offset:-right_offset
                    ]

            # write tile
            out_ds = gdal.Open(out_raster_path.as_posix(), gdal.GA_Update)
            if out_raster_nr_of_bands == 1:  # one band_number
                out_ds.GetRasterBand(1).WriteArray(visualization_array, x, y)
                out_ds.FlushCache()
            else:  # multiple bands
                for i_band in range(out_raster_nr_of_bands):
                    band = i_band + 1
                    out_ds.GetRasterBand(band).WriteArray(
                        visualization_array[i_band], x, y
                    )
                    out_ds.FlushCache()
            out_ds = None
    dem_ds = None


def _create_rvt_visualization_blank_raster(
    rvt_visualization: "rvt.default.RVTVisualization",
    rvt_default: "rvt.default.DefaultValues",
    dem_path: Path,
    output_dir_path: Path,
    dem_ds: gdal.Dataset,
    save_float: bool,
    save_8bit: bool,
) -> None:
    """ "Create blank raster or rasters for rvt_visualization to later store visualization in it tile by tile."""
    if save_float:
        out_float_path = rvt_default.get_visualization_path(
            rvt_visualization=rvt_visualization,
            dem_path=dem_path,
            output_dir_path=output_dir_path,
            path_8bit=False,
        )
        nr_bands = 1
        if rvt_visualization == rvt.default.RVTVisualization.MULTI_HILLSHADE:
            nr_bands = rvt_default.mhs_nr_dir
        _create_blank_raster(
            in_data_set=dem_ds,
            out_raster_path=out_float_path,
            nr_bands=nr_bands,
            e_type=6,
        )
    if save_8bit:
        out_8bit_path = rvt_default.get_visualization_path(
            rvt_visualization=rvt_visualization,
            dem_path=dem_path,
            output_dir_path=output_dir_path,
            path_8bit=True,
        )
        nr_bands = 1
        if (
            rvt_visualization == rvt.default.RVTVisualization.MULTI_HILLSHADE
            or rvt_visualization
            == rvt.default.RVTVisualization.MULTI_SCALE_TOPOGRAPHIC_POSITION
        ):
            nr_bands = 3
        _create_blank_raster(
            in_data_set=dem_ds,
            out_raster_path=out_8bit_path,
            nr_bands=nr_bands,
            e_type=1,
        )


def _get_rvt_visualization_overlap(
    rvt_visualization: "rvt.default.RVTVisualization",
    rvt_default: "rvt.default.DefaultValues",
) -> int:
    if rvt_visualization == rvt.default.RVTVisualization.SLOPE:
        return 1
    elif rvt_visualization == rvt.default.RVTVisualization.HILLSHADE:
        return 1
    elif rvt_visualization == rvt.default.RVTVisualization.SHADOW:
        return 1
    elif rvt_visualization == rvt.default.RVTVisualization.MULTI_HILLSHADE:
        return 1
    elif rvt_visualization == rvt.default.RVTVisualization.SIMPLE_LOCAL_RELIEF_MODEL:
        return int(rvt_default.slrm_rad_cell)
    elif rvt_visualization == rvt.default.RVTVisualization.SKY_VIEW_FACTOR:
        return int(rvt_default.svf_r_max)
    elif rvt_visualization == rvt.default.RVTVisualization.ANISOTROPIC_SKY_VIEW_FACTOR:
        return int(rvt_default.svf_r_max)
    elif rvt_visualization == rvt.default.RVTVisualization.POSITIVE_OPENNESS:
        return int(rvt_default.svf_r_max)
    elif rvt_visualization == rvt.default.RVTVisualization.NEGATIVE_OPENNESS:
        return int(rvt_default.svf_r_max)
    elif rvt_visualization == rvt.default.RVTVisualization.SKY_ILLUMINATION:
        return int(rvt_default.sim_shadow_dist)
    elif rvt_visualization == rvt.default.RVTVisualization.LOCAL_DOMINANCE:
        return int(rvt_default.ld_max_rad)
    elif rvt_visualization == rvt.default.RVTVisualization.MULTI_SCALE_RELIEF_MODEL:
        return int(rvt_default.msrm_feature_max)
    elif (
        rvt_visualization
        == rvt.default.RVTVisualization.MULTI_SCALE_TOPOGRAPHIC_POSITION
    ):
        return int(rvt_default.mstp_broad_scale[1])
    else:
        raise ValueError(f"Overlap is not defined for {rvt_visualization.name}.")


def save_rvt_visualization_tile_by_tile(
    rvt_visualization: "rvt.default.RVTVisualization",
    rvt_default: "rvt.default.DefaultValues",
    dem_path: Path,
    output_dir_path: Optional[Path] = None,
    save_float: bool = True,
    save_8bit: bool = False,
) -> None:
    """
    Some DEMs are too large to load them into memory. This function reads dem raster tile by tile,
    calculates RVT visualization on it tile by tile and than saves calculated visualization tile by tile in out raster.
    This function can silmultaniously store float and 8bit version of visualization (where possible).

    Parameters
    ----------
    rvt_visualization : RVTVisualization
        RVT visualization.
    rvt_default : Default
        Class where RVT parameters are stored.
    dem_path : Path
        Path to a Digital elevation model.
    output_dir_path : Path
        Out directory to save visualizations. If None it uses dem_dir from dem_path.
    save_float : bool
        If save float.
    save_8bit : bool
        If save 8bit.

    Returns
    -------
    out : None
    """
    if not save_float and not save_8bit:
        Exception(
            "rvt.tile.save_visualization_tile_by_tile: At least one of save_float or save_8bit must be true!"
        )
    if not dem_path.exists():
        Exception(
            "rvt.tile.save_visualization_tile_by_tile: Input dem path does not exist!"
        )

    tile_size_x = rvt_default.tile_size[0]
    tile_size_y = rvt_default.tile_size[1]

    if tile_size_x < 50 or tile_size_y < 50:
        Exception(
            "rvt.tile.save_visualization_tile_by_tile: Tile size too small (tile_size_x, tile_size_y),"
            " it needs to be bigger than 50 pixels!"
        )

    if output_dir_path is None:
        output_dir_path = dem_path.parent

    dem_ds = gdal.Open(dem_path.as_posix())
    gt = dem_ds.GetGeoTransform()
    x_res = gt[1]  # x_resolution
    y_res = -gt[5]  # y_resolution
    no_data = dem_ds.GetRasterBand(1).GetNoDataValue()
    band = dem_ds.GetRasterBand(1)
    x_size = band.XSize  # number of columns
    y_size = band.YSize  # number of rows

    _create_rvt_visualization_blank_raster(
        rvt_visualization=rvt_visualization,
        rvt_default=rvt_default,
        dem_path=dem_path,
        output_dir_path=output_dir_path,
        dem_ds=dem_ds,
        save_float=save_float,
        save_8bit=save_8bit,
    )

    overlap = _get_rvt_visualization_overlap(
        rvt_visualization=rvt_visualization, rvt_default=rvt_default
    )

    for y in range(0, y_size, tile_size_y):
        if y + tile_size_y < y_size:  # if rows overlap
            rows = tile_size_y
        else:
            rows = y_size - y
        for x in range(0, x_size, tile_size_x):
            if x + tile_size_x < x_size:  # if cols overlap
                cols = tile_size_x
            else:
                cols = x_size - x

            # get offset for each tile, check edges
            left_offset = 0
            right_offset = 0
            top_offset = 0
            bottom_offset = 0
            if x != 0:  # left offset
                if x - overlap < 0:  # left overlap
                    left_offset = x
                else:
                    left_offset = overlap
            if x + cols != x_size:  # right offset
                if x + cols + overlap > x_size:  # right overlap
                    right_offset = x_size - x - cols
                else:
                    right_offset = overlap
            if y != 0:  # apply top offset
                if y - overlap < 0:  # top overlap
                    top_offset = y
                else:
                    top_offset = overlap
            if y + rows != y_size:  # bottom offset
                if y + rows + overlap > y_size:  # bottom overlap
                    bottom_offset = y_size - y - rows
                else:
                    bottom_offset = overlap

            # read tile
            x_off = x - left_offset
            y_off = y - top_offset
            cols_off = cols + left_offset + right_offset
            rows_off = rows + top_offset + bottom_offset

            tile_array = np.array(
                dem_ds.GetRasterBand(1).ReadAsArray(x_off, y_off, cols_off, rows_off)
            )

            (
                visualization_float_arr,
                visualization_8bit_arr,
            ) = rvt_default.calculate_visualization(
                visualization=rvt_visualization,
                dem=tile_array,
                resolution_x=x_res,
                resolution_y=y_res,
                no_data=no_data,
                save_float=save_float,
                save_8bit=save_8bit,
            )

            # remove offset from visualization block
            if visualization_float_arr is not None and save_float:
                if visualization_float_arr.ndim == 2:
                    if right_offset == 0 and bottom_offset == 0:
                        visualization_float_arr = visualization_float_arr[
                            top_offset:, left_offset:
                        ]
                    elif right_offset == 0:
                        visualization_float_arr = visualization_float_arr[
                            top_offset:-bottom_offset, left_offset:
                        ]
                    elif bottom_offset == 0:
                        visualization_float_arr = visualization_float_arr[
                            top_offset:, left_offset:-right_offset
                        ]
                    else:
                        visualization_float_arr = visualization_float_arr[
                            top_offset:-bottom_offset, left_offset:-right_offset
                        ]
                else:
                    if right_offset == 0 and bottom_offset == 0:
                        visualization_float_arr = visualization_float_arr[
                            :, top_offset:, left_offset:
                        ]
                    elif right_offset == 0:
                        visualization_float_arr = visualization_float_arr[
                            :, top_offset:-bottom_offset, left_offset:
                        ]
                    elif bottom_offset == 0:
                        visualization_float_arr = visualization_float_arr[
                            :, top_offset:, left_offset:-right_offset
                        ]
                    else:
                        visualization_float_arr = visualization_float_arr[
                            :, top_offset:-bottom_offset, left_offset:-right_offset
                        ]
            if visualization_8bit_arr is not None and save_8bit:
                if visualization_8bit_arr.ndim == 2:
                    if right_offset == 0 and bottom_offset == 0:
                        visualization_8bit_arr = visualization_8bit_arr[
                            top_offset:, left_offset:
                        ]
                    elif right_offset == 0:
                        visualization_8bit_arr = visualization_8bit_arr[
                            top_offset:-bottom_offset, left_offset:
                        ]
                    elif bottom_offset == 0:
                        visualization_8bit_arr = visualization_8bit_arr[
                            top_offset:, left_offset:-right_offset
                        ]
                    else:
                        visualization_8bit_arr = visualization_8bit_arr[
                            top_offset:-bottom_offset, left_offset:-right_offset
                        ]
                else:
                    if right_offset == 0 and bottom_offset == 0:
                        visualization_8bit_arr = visualization_8bit_arr[
                            :, top_offset:, left_offset:
                        ]
                    elif right_offset == 0:
                        visualization_8bit_arr = visualization_8bit_arr[
                            :, top_offset:-bottom_offset, left_offset:
                        ]
                    elif bottom_offset == 0:
                        visualization_8bit_arr = visualization_8bit_arr[
                            :, top_offset:, left_offset:-right_offset
                        ]
                    else:
                        visualization_8bit_arr = visualization_8bit_arr[
                            :, top_offset:-bottom_offset, left_offset:-right_offset
                        ]

            # write tile
            if visualization_float_arr is not None and save_float:
                out_visualization_float_path = rvt_default.get_visualization_path(
                    rvt_visualization=rvt_visualization,
                    dem_path=dem_path,
                    output_dir_path=output_dir_path,
                    path_8bit=False,
                )
                out_ds_float = gdal.Open(
                    out_visualization_float_path.as_posix(), gdal.GA_Update
                )
                if visualization_float_arr.ndim == 2:
                    out_ds_float.GetRasterBand(1).WriteArray(
                        visualization_float_arr, x, y
                    )
                    out_ds_float.FlushCache()
                else:
                    for i_band in range(visualization_float_arr.shape[0]):
                        band = i_band + 1
                        out_ds_float.GetRasterBand(band).WriteArray(
                            visualization_float_arr[i_band], x, y
                        )
                        out_ds_float.FlushCache()
                out_ds_float = None
            if visualization_8bit_arr is not None and save_8bit:  # multiple bands
                out_visualization_8bit_path = rvt_default.get_visualization_path(
                    rvt_visualization=rvt_visualization,
                    dem_path=dem_path,
                    output_dir_path=output_dir_path,
                    path_8bit=True,
                )
                out_ds_8bit = gdal.Open(
                    out_visualization_8bit_path.as_posix(), gdal.GA_Update
                )
                if visualization_8bit_arr.ndim == 2:
                    out_ds_8bit.GetRasterBand(1).WriteArray(
                        visualization_8bit_arr, x, y
                    )
                    out_ds_8bit.FlushCache()
                else:
                    for i_band in range(visualization_8bit_arr.shape[0]):
                        band = i_band + 1
                        out_ds_8bit.GetRasterBand(band).WriteArray(
                            visualization_8bit_arr[i_band], x, y
                        )
                        out_ds_8bit.FlushCache()
                out_ds_8bit = None

    dem_ds = None
