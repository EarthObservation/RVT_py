"""
Tiled processing of DEMs too large to visualize in memory.

:func:`compute_tiled` reads a DEM window by window, computes a visualization
per tile, and writes each result into a fixed window of a pre-created output
raster. Every window is read expanded by the config's
``required_overlap(resolution)`` pixels, clamped to the raster bounds (never
boundless), and the overlap is sliced off before writing — so tile borders see
real data and raster borders rely on each visualization's own edge padding,
exactly as an untiled compute does. Tiled output therefore equals untiled
output.

Only single-image configs are supported: multi-output configs
(:class:`~rvt.visualization.SlopeAspect`,
:class:`~rvt.visualization.SkyViewFactor`,
:class:`~rvt.visualization.ShadowHorizon`) raise :class:`TypeError` — compute
them untiled via their own ``compute``. Configs whose result depends on data
outside any fixed neighbourhood (or whose output shape differs from the
input) declare ``TILEABLE = False`` and raise :class:`ValueError`.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import NamedTuple

import numpy as np
import rasterio
from affine import Affine
from numpy.typing import NDArray
from rasterio.io import DatasetWriter
from rasterio.windows import Window

from rvt.core import Resolution, Visualization
from rvt.visualization import ShadowHorizon, SkyViewFactor, SlopeAspect

DEFAULT_TILE_SIZE = (4000, 4000)

_MULTI_OUTPUT_CONFIGS = (ShadowHorizon, SkyViewFactor, SlopeAspect)


#: A visualization config :func:`compute_tiled` can run: the shared
#: :class:`rvt.core.Visualization` contract, specialized to configs whose
#: ``compute`` yields a single image (2D, or a fixed band stack of the DEM's
#: shape) rather than a multi-output result object.
TileableVisualization = Visualization[NDArray[np.float32]]


class _TileTask(NamedTuple):
    """One core (overlap-free) tile window, in pixel coordinates."""

    row_offset: int
    column_offset: int
    height: int
    width: int


def _read_resolution(transform: Affine) -> Resolution:
    """Return the absolute cell size of a geotransform."""
    return Resolution(abs(transform.a), abs(transform.e))


def _compute_tile(
    config: TileableVisualization,
    input_path: str,
    task: _TileTask,
    overlap: int,
) -> NDArray[np.float32]:
    """
    Read one tile expanded by ``overlap`` (clamped to the raster bounds),
    compute the visualization, and slice the result back to the core window.

    Module-level so a :class:`~concurrent.futures.ProcessPoolExecutor` worker
    can pickle it; the input raster is opened per call for the same reason.
    """
    with rasterio.open(input_path) as source:
        expanded = Window(
            task.column_offset - overlap,
            task.row_offset - overlap,
            task.width + 2 * overlap,
            task.height + 2 * overlap,
        ).intersection(Window(0, 0, source.width, source.height))
        dem = source.read(1, window=expanded).astype(np.float32)
        nodata = source.nodata
        resolution = _read_resolution(source.transform)
    if nodata is not None and not np.isnan(nodata):
        dem[dem == nodata] = np.nan

    result = np.asarray(config.compute(dem, resolution), dtype=np.float32)
    top = task.row_offset - int(expanded.row_off)
    left = task.column_offset - int(expanded.col_off)
    return result[..., top : top + task.height, left : left + task.width]


def _write_tile(
    destination: DatasetWriter, task: _TileTask, tile: NDArray[np.float32]
) -> None:
    """Write one computed tile into its core window of the output raster."""
    window = Window(task.column_offset, task.row_offset, task.width, task.height)
    if tile.ndim == 2:
        destination.write(tile, indexes=1, window=window)
    else:
        destination.write(tile, window=window)


def compute_tiled(
    config: TileableVisualization,
    input_path: str | Path,
    output_path: str | Path,
    tile_size: tuple[int, int] = DEFAULT_TILE_SIZE,
    workers: int = 1,
) -> Path:
    """
    Compute a visualization of a DEM raster tile by tile and write the result
    as a float32 GeoTIFF (NaN nodata, LZW-compressed, internally tiled),
    equal to the untiled compute of the whole DEM.

    Each tile is read expanded by ``config.required_overlap(resolution)``
    pixels, clamped to the raster bounds; the overlap is sliced off before the
    core window is written. With ``workers`` greater than 1 the tiles are
    computed by a process pool; the output is deterministic and identical to a
    serial run either way, because every tile is written to a fixed window of
    the pre-created output raster.

    :param config: The visualization config to compute. Must be a tileable,
        single-image config.
    :param input_path: Path of the DEM raster to read.
    :param output_path: Path of the GeoTIFF to create (overwritten if present).
    :param tile_size: Core tile size as ``(rows, columns)`` in pixels.
    :param workers: Number of worker processes; 1 computes serially in-process.
    :returns: ``output_path`` as a :class:`~pathlib.Path`.
    :raises TypeError: If the config produces multiple outputs (slope/aspect,
        the sky-view factor family, shadow/horizon) — compute those untiled
        via their own ``compute``.
    :raises ValueError: If the config is not tileable (``TILEABLE`` is false:
        its result depends on data outside any fixed neighbourhood, or its
        output shape differs from the input), or if ``tile_size`` or
        ``workers`` is not positive.
    """
    if isinstance(config, _MULTI_OUTPUT_CONFIGS):
        raise TypeError(
            f"{type(config).__name__} produces multiple outputs; compute_tiled"
            " supports single-image configs only — compute it untiled via its"
            " own compute method."
        )
    if not config.TILEABLE:
        raise ValueError(
            f"{type(config).__name__} cannot be computed tiled: its result is"
            " not independent of the tiling (data-dependent normalisation, or"
            " an output shape differing from the input); compute it untiled"
            " instead."
        )
    tile_height, tile_width = tile_size
    if tile_height < 1 or tile_width < 1:
        raise ValueError(f"tile_size must be positive, got '{tile_size}'!")
    if workers < 1:
        raise ValueError(f"workers must be positive, got '{workers}'!")

    input_path = Path(input_path)
    output_path = Path(output_path)
    with rasterio.open(input_path) as source:
        width, height = source.width, source.height
        transform, crs = source.transform, source.crs
    overlap = config.required_overlap(_read_resolution(transform))

    tasks = [
        _TileTask(
            row_offset=row_offset,
            column_offset=column_offset,
            height=min(tile_height, height - row_offset),
            width=min(tile_width, width - column_offset),
        )
        for row_offset in range(0, height, tile_height)
        for column_offset in range(0, width, tile_width)
    ]

    # The first tile is computed up front to learn the output band count
    # before the output raster is created.
    first_tile = _compute_tile(config, str(input_path), tasks[0], overlap)
    band_count = 1 if first_tile.ndim == 2 else first_tile.shape[0]

    with rasterio.open(
        output_path,
        mode="w",
        driver="GTiff",
        width=width,
        height=height,
        count=band_count,
        dtype="float32",
        transform=transform,
        crs=crs,
        nodata=np.nan,
        tiled=True,
        blockxsize=512,
        blockysize=512,
        compress="lzw",
        bigtiff="IF_SAFER",
    ) as destination:
        _write_tile(destination, tasks[0], first_tile)
        if workers == 1:
            for task in tasks[1:]:
                _write_tile(
                    destination,
                    task,
                    _compute_tile(config, str(input_path), task, overlap),
                )
        else:
            with ProcessPoolExecutor(max_workers=workers) as pool:
                futures = [
                    pool.submit(_compute_tile, config, str(input_path), task, overlap)
                    for task in tasks[1:]
                ]
                for task, future in zip(tasks[1:], futures, strict=True):
                    _write_tile(destination, task, future.result())
    return output_path
