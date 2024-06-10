"""
Script containing tools for working with reference grid.
Created on 26 May 2023
@author: Nejc Čož, ZRC SAZU, Novi trg 2, 1000 Ljubljana, Slovenia
"""

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import shapes
from shapely.geometry import box
from shapely.geometry import shape


def bounding_grid(raster_file, tile_size_pix, tag=False, grid_type="GDF", save_gdf=None):
    """Creates bounding grid based on the extents of VRT file.

    Parameters
    ----------
    raster_file : str
        Path to VRT file (Virtual Raster Mosaic)
    tile_size_pix : int
        Tile size in pixels
    tag : bool
        Target Aligned Grid - rounds up coordinates of grid cells to nearest multiple of the tile size.
    grid_type : str
        Type of output. Can be: "GDF" - GeoDataFrame; "extents" - list of extents;
    save_gdf : str
        Use to export GeoDataFrame to disk, same path as vrt_file. Use "SHP" or "GPKG" or "GeoJSON".

    Returns
    -------
    gpd.geodataframe.GeoDataFrame
        List of tuples, each tuple represents one grid cell with coordinates (left, bottom, right, top).

    """
    # Get information about VRT raster
    with rasterio.open(raster_file) as src:
        extents = src.bounds
        res = src.res
        crs = src.crs

    # transform pixels to meters
    tile_w = tile_size_pix * res[0]
    tile_h = tile_w  # square

    # TAG TO GRID (ROUND TO NEAREST x)
    if tag:
        left = np.floor(extents.left / tile_w) * tile_w
        # bottom = np.floor(extents.bottom / tile_h) * tile_h
        # right = np.ceil(extents.right / tile_w) * tile_w
        top = np.ceil(extents.top / tile_h) * tile_h
        _, bottom, right, _ = extents  # ONLY TOP-LEFT NEEDS TO BE ROUNDED
    else:
        left, bottom, right, top = extents

    # List all individual grid cells
    grid_cells = []
    for x0 in np.arange(left, right, tile_w):
        for y1 in np.arange(top, bottom, -tile_h):
            # bounds
            x1 = x0 + tile_w
            y0 = y1 - tile_h
            grid_cells.append((x0, y0, x1, y1))

    # Output in the correct format
    if grid_type == "GDF":
        out_grid = gpd.GeoDataFrame([box(*a) for a in grid_cells], columns=['geometry'], crs=crs)
        if save_gdf:
            output_name = raster_file.rstrip(".vrt") + f"_{tile_size_pix}pix." + save_gdf.lower()
            if save_gdf == "SHP":
                save_gdf = "ESRI Shapefile"
            out_grid.to_file(output_name, driver=save_gdf)
    elif grid_type == "extents":
        out_grid = grid_cells
    else:
        print("Error: select either 'GDF' or 'extents'!")
        out_grid = None

    return out_grid


def filter_by_outline(in_grid, outline_file, save_gpkg=False, save_path=None):
    """Filters the grid in GDF format to retain only the tiles that intersect with the outline.

    WARNING: this function only works for simple continuous areas/shapes.
    
    Parameters
    ----------
    in_grid : gpd.geodataframe.GeoDataFrame, str
        Grid made of polygons in GeoDataFormat
    outline_file : gpd.geodataframe.GeoDataFrame, str
        Path to the outline file (any format readable by GeoPandas).
    save_gpkg : bool, default False
        If true, grid is saved to disk.
    save_path : str or pathlib.Path, default None
        Path to the directory where the results (vector file) will be saved.

    Returns
    -------
    gpd.geodataframe.GeoDataFrame
        Filtered grid, with tile_ID and extents.
    """
    # If it is GDF then pass, else read from file
    if type(outline_file) is gpd.geodataframe.GeoDataFrame:
        outline = outline_file
    elif type(outline_file) is str:
        outline = gpd.read_file(outline_file)
    else:
        raise ValueError('Wrong "outline_file" format - has to be GDF or str (path to file)!')

    # Outline has to be a single feature (can be multipolygon).
    if outline.shape[0] != 1:
        outline = outline.dissolve()

    # Two options
    # a) without sindex
    # outline_filter = in_grid.intersects(outline.geometry[0])
    # out_grid = in_grid[outline_filter]

    # b) Using sindex
    if not in_grid.has_sindex:
        _ = in_grid.sindex
    outline_filter = in_grid.sindex.query(outline.geometry[0], predicate="intersects")
    out_grid = in_grid.iloc[outline_filter].reset_index(drop=True)

    # Add cell_ID and extents columns. Extents are (L, B, R, T).
    out_grid = out_grid.reset_index()
    out_grid = out_grid.rename(columns={'index': 'tile_ID'})
    out_grid["extents"] = out_grid.bounds.apply(lambda x: (x.minx, x.miny, x.maxx, x.maxy), axis=1)

    if save_gpkg:
        if save_path:
            # Because tuple can't be saved into file, split extents into separate columns
            out_grid[["minx", "miny", "maxx", "maxy"]] = pd.DataFrame(
                out_grid['extents'].tolist(),
                index=out_grid.index
            )
            out_grid = out_grid.drop(columns=['extents'])
            out_grid.to_file(save_path, driver="GPKG")
        else:
            raise ValueError("Specify save_path for the output file!")

    return out_grid


def poly_from_valid(tif_pth, save_gpkg=False):
    # Read raster data
    with rasterio.open(tif_pth) as src:
        raster = src.read()
        nodata = src.nodata
        crs = src.crs
        transform = src.transform

    # Set nodaata to 0 and valid data to 1
    if nodata is None:
        raster[:] = 1
    elif np.isnan(nodata):
        raster[~np.isnan(raster)] = 1
        raster[np.isnan(raster)] = 0
    else:
        raster[raster != nodata] = 1
        raster[raster == nodata] = 0

    # Outputs a list of (polygon, value) tuples
    output = list(shapes(raster, transform=transform))

    # Find polygon covering valid data (value = 1) and transform to GDF friendly format
    poly = []
    for polygon, value in output:
        if value == 1:
            poly.append(shape(polygon))

    # Make Geodataframe
    grid = gpd.GeoDataFrame(poly, columns=['geometry'], crs=crs)

    if save_gpkg:
        new_name = tif_pth[:-4] + "_validDataMask.gpkg"
        grid.to_file(new_name, driver="GPKG")
        return Path(new_name)
    else:
        return grid
