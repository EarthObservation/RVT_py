"""
Relief Visualization Toolbox: DEM visualizations and blends.
"""

from __future__ import annotations

from typing import Protocol

import numpy as np
from numpy.typing import NDArray

from rvt.core import Resolution
from rvt.io import Raster

__version__ = "3.0.0.dev0"


class SingleImageVisualization(Protocol):
    """
    A visualization config whose ``compute`` yields one array (2D, or a band stack).
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
    """
    array = config.compute(raster.array, raster.resolution)
    return Raster(array=array, transform=raster.transform, crs=raster.crs)
