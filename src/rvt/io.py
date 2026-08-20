"""
Rasterio-backed raster input/output.

The compute layer works exclusively with float32 arrays where nodata is NaN;
this module converts to and from that canonical form at the file boundary.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import rasterio
from affine import Affine
from numpy.typing import NDArray
from rasterio.crs import CRS

from rvt.core import Resolution


@dataclass(frozen=True)
class Raster:
    """
    An in-memory georeferenced raster.

    :param array: 2D single-band array, or ``(bands, rows, columns)`` for multi-band.
        Float arrays use NaN as nodata; uint8 arrays have no nodata concept.
    :param transform: Affine geotransform mapping pixel to map coordinates.
    :param crs: Coordinate reference system, or ``None`` when the source has none.
    :param source_nodata: Nodata value of the file this raster was read from, kept
        for provenance only — the in-memory array already uses NaN.
    """

    array: NDArray[np.float32] | NDArray[np.uint8]
    transform: Affine
    crs: CRS | None
    source_nodata: float | None = None

    @property
    def resolution(self) -> Resolution:
        """
        Cell size derived from the geotransform, always positive.
        """
        return Resolution(abs(self.transform.a), abs(self.transform.e))

    @property
    def band_count(self) -> int:
        """
        Number of bands (1 for a 2D array).
        """
        if self.array.ndim == 2:
            return 1
        return self.array.shape[0]


def read(path: str | Path) -> Raster:
    """
    Read a raster into the canonical in-memory form: float32 with NaN nodata.

    The pixel values are read raw (no masking, no scale/offset application) and
    then converted, matching the legacy GDAL behaviour bit for bit.

    :param path: Path of the raster file to read.
    :returns: The raster as float32 with nodata converted to NaN.
    """
    with rasterio.open(path) as source:
        array = source.read().astype(np.float32)
        nodata = source.nodata
        if nodata is not None and not np.isnan(nodata):
            array[array == nodata] = np.nan
        if array.shape[0] == 1:
            array = array[0]
        return Raster(
            array=array,
            transform=source.transform,
            crs=source.crs,
            source_nodata=nodata,
        )


def write(raster: Raster, path: str | Path) -> None:
    """
    Write a raster as an internally tiled, LZW-compressed GeoTIFF.

    Float rasters are tagged with NaN nodata; uint8 rasters carry no nodata tag
    (NaN pixels are already encoded by the stretch step, legacy semantics).

    :param raster: The raster to write.
    :param path: Destination file path.
    """
    array = raster.array if raster.array.ndim == 3 else raster.array[np.newaxis, :, :]
    is_float = np.issubdtype(array.dtype, np.floating)
    with rasterio.open(
        path,
        mode="w",
        driver="GTiff",
        width=array.shape[2],
        height=array.shape[1],
        count=array.shape[0],
        dtype=array.dtype.name,
        transform=raster.transform,
        crs=raster.crs,
        nodata=np.nan if is_float else None,
        tiled=True,
        blockxsize=512,
        blockysize=512,
        compress="lzw",
        bigtiff="IF_SAFER",
    ) as destination:
        destination.write(array)
