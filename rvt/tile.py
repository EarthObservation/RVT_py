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

Copyright:
    2010-2022 Research Centre of the Slovenian Academy of Sciences and Arts
    2016-2022 University of Ljubljana, Faculty of Civil and Geodetic Engineering
"""
from pathlib import Path
from typing import Callable, Dict, Any, Optional, Union
import numpy as np
from osgeo import gdal


def create_blank_raster(
        in_data_set: gdal.Dataset,
        out_raster_path: Path,
        nr_bands: int = 1,
        no_data: float = np.nan,
        e_type: int = 6,
):
    """Takes input data set and creates new raster. It copies input data set size, projection and geo info."""
    gtiff_driver = gdal.GetDriverByName("GTiff")
    band = in_data_set.GetRasterBand(1)
    x_size = band.XSize  # number of columns
    y_size = band.YSize  # number of rows
    out_ds = gtiff_driver.Create(out_raster_path.as_posix(), xsize=x_size, ysize=y_size, bands=nr_bands, eType=e_type,
                                 options=["BIGTIFF=IF_NEEDED"])
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
        For example rvt.vis.slope_aspect outputs dictionary with keys "slope" and "aspect".
        To select slope set this parameter to "slope".

    Returns
    -------
    out : None
    """
    if not dem_path.exists():
        Exception("rvt.tile.save_visualization_tile_by_tile: Input dem path does not exist!")
    if tile_size_x < 50 or tile_size_y < 50:
        Exception("rvt.tile.save_visualization_tile_by_tile: Tile size too small (tile_size_x, tile_size_y),"
                  " it needs to be bigger than 50 pixels!")

    dem_ds = gdal.Open(dem_path.as_posix())
    gt = dem_ds.GetGeoTransform()
    x_res = gt[1]  # x_resolution
    y_res = -gt[5]  # y_resolution
    no_data = dem_ds.GetRasterBand(1).GetNoDataValue()
    band = dem_ds.GetRasterBand(1)
    x_size = band.XSize  # number of columns
    y_size = band.YSize  # number of rows
    del band

    # set resolution and no_data function_parameters if needed (are set to None) from dem
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

    create_blank_raster(in_data_set=dem_ds, out_raster_path=out_raster_path, nr_bands=out_raster_nr_of_bands,
                        e_type=out_raster_e_type)

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

            tile_array = np.array(dem_ds.GetRasterBand(1).ReadAsArray(x_off, y_off, cols_off, rows_off))
            if function_parameters is not None:
                visualization_array = visualization_function(dem=tile_array, **function_parameters)
            else:
                visualization_array = visualization_function(dem=tile_array)

            if out_visualization_dict_key is not None:
                visualization_array = visualization_array[out_visualization_dict_key]

            # remove offset from visualization block
            if out_raster_nr_of_bands == 1:
                if right_offset == 0 and bottom_offset == 0:
                    visualization_array = visualization_array[top_offset:, left_offset:]
                elif right_offset == 0:
                    visualization_array = visualization_array[top_offset:-bottom_offset, left_offset:]
                elif bottom_offset == 0:
                    visualization_array = visualization_array[top_offset:, left_offset:-right_offset]
                else:
                    visualization_array = visualization_array[top_offset:-bottom_offset, left_offset:-right_offset]
            else:
                if right_offset == 0 and bottom_offset == 0:
                    visualization_array = visualization_array[:, top_offset:, left_offset:]
                elif right_offset == 0:
                    visualization_array = visualization_array[:, top_offset:-bottom_offset, left_offset:]
                elif bottom_offset == 0:
                    visualization_array = visualization_array[:, top_offset:, left_offset:-right_offset]
                else:
                    visualization_array = visualization_array[:, top_offset:-bottom_offset, left_offset:-right_offset]

            # write tile
            out_ds = gdal.Open(out_raster_path.as_posix(), gdal.GA_Update)
            if out_raster_nr_of_bands == 1:  # one band
                out_ds.GetRasterBand(1).WriteArray(visualization_array, x, y)
                out_ds.FlushCache()
            else:  # multiple bands
                for i_band in range(out_raster_nr_of_bands):
                    band = i_band + 1
                    out_ds.GetRasterBand(band).WriteArray(visualization_array[i_band], x, y)
                    out_ds.FlushCache()
            out_ds = None
    dem_ds = None
