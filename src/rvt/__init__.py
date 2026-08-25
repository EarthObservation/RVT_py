"""
Relief Visualization Toolbox: DEM visualizations and blends.
"""

from __future__ import annotations

from typing import Protocol

import numpy as np
from numpy.typing import NDArray

from rvt.core import Resolution
from rvt.core import Visualization as Visualization
from rvt.core import VisualizationOutput as VisualizationOutput
from rvt.io import Raster
from rvt.session import Session as Session
from rvt.session import SessionTask as SessionTask

__version__ = "3.0.0.dev0"


class SingleImageVisualization(Protocol):
    """
    A visualization config whose ``compute`` yields one array (2D, or a band
    stack).

    Deliberately narrower than :class:`rvt.core.Visualization`: :func:`compute`
    only ever calls ``compute``, so anything array-producing qualifies, without
    having to name files or declare tileability.
    """

    def compute(
        self, dem: NDArray[np.float32], resolution: Resolution
    ) -> NDArray[np.float32]: ...


def compute(config: SingleImageVisualization, raster: Raster) -> Raster:
    """
    Compute a visualization of a raster, keeping its georeferencing.

    Multi-output configs (e.g. slope/aspect, the sky-view factor family) return
    typed result objects from their own ``compute`` — call them directly and wrap
    the members you need.

    :param config: A visualization config producing a single image.
    :param raster: The DEM as read by :func:`rvt.io.read`.
    :returns: The visualization as a raster with the DEM's transform and CRS.
    :raises ValueError: If the raster is not a float DEM.
    """
    dem = raster.array
    if not np.issubdtype(dem.dtype, np.floating):
        raise ValueError(f"expected a float DEM raster, got dtype '{dem.dtype}'!")
    array = config.compute(dem.astype(np.float32, copy=False), raster.resolution)
    return Raster(array=array, transform=raster.transform, crs=raster.crs)
