"""
Simple local relief model.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

import numpy as np
from numpy.typing import NDArray

from rvt.core import (
    Resolution,
    Visualization,
    VisualizationOutput,
    validate_dem_and_resolution,
    validate_vertical_exaggeration,
)
from rvt.stretch import Stretch, StretchMode

MINIMUM_RADIUS_CELLS = 10
MAXIMUM_RADIUS_CELLS = 50


def _integral_image(
    array: NDArray[np.floating], data_type: type = np.float64
) -> NDArray[np.floating]:
    """Summed-area table with the origin in the upper-left corner."""
    return array.astype(data_type).cumsum(axis=0).cumsum(axis=1)


def _mean_filter(dem: NDArray[np.float32], kernel_radius: int) -> NDArray[np.float32]:
    """
    Legacy mean (low-pass) filter via integral images, ported verbatim for
    golden-master parity. Kernel size is ``2 * kernel_radius + 1``.
    """
    radius_cells = int(kernel_radius)

    if kernel_radius == 0:
        return dem

    nan_mask = np.isnan(dem)

    dem_padded = np.pad(dem, (radius_cells + 1, radius_cells), mode="edge")
    padded_nan_mask = np.isnan(dem_padded)
    dem_padded[padded_nan_mask] = 0

    pixel_counts = np.ones(dem_padded.shape)
    pixel_counts[padded_nan_mask] = 0
    pixel_counts_integral = _integral_image(pixel_counts, np.int64)

    dem_integral = _integral_image(dem_padded)

    kernel_pixel_counts = (
        np.roll(pixel_counts_integral, (radius_cells, radius_cells), axis=(0, 1))
        + np.roll(
            pixel_counts_integral,
            (-radius_cells - 1, -radius_cells - 1),
            axis=(0, 1),
        )
        - np.roll(pixel_counts_integral, (-radius_cells - 1, radius_cells), axis=(0, 1))
        - np.roll(pixel_counts_integral, (radius_cells, -radius_cells - 1), axis=(0, 1))
    )
    mean = (
        np.roll(dem_integral, (radius_cells, radius_cells), axis=(0, 1))
        + np.roll(dem_integral, (-radius_cells - 1, -radius_cells - 1), axis=(0, 1))
        - np.roll(dem_integral, (-radius_cells - 1, radius_cells), axis=(0, 1))
        - np.roll(dem_integral, (radius_cells, -radius_cells - 1), axis=(0, 1))
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        mean = mean / kernel_pixel_counts
    mean = mean.astype(np.float32)
    mean = mean[radius_cells : -(radius_cells + 1), radius_cells : -(radius_cells + 1)]
    mean[nan_mask] = np.nan

    return mean


@dataclass(frozen=True)
class SimpleLocalReliefModel(Visualization[NDArray[np.float32]]):
    """
    Simple local relief model: the DEM minus its mean-filtered trend surface.

    :param radius_cells: Radius for trend assessment in pixels, 10-50.
    :param vertical_exaggeration: Factor the DEM heights are multiplied by.
    """

    EIGHT_BIT_STRETCH: ClassVar[Stretch] = Stretch(StretchMode.VALUE, -2.0, 2.0)
    TILEABLE: ClassVar[bool] = True

    radius_cells: int = 20
    vertical_exaggeration: float = 1.0

    def __post_init__(self) -> None:
        validate_vertical_exaggeration(self.vertical_exaggeration)
        if not (MINIMUM_RADIUS_CELLS <= self.radius_cells <= MAXIMUM_RADIUS_CELLS):
            raise ValueError(
                f"radius_cells must be between {MINIMUM_RADIUS_CELLS} and"
                f" {MAXIMUM_RADIUS_CELLS}, got '{self.radius_cells}'!"
            )

    def output_stem(self, dem_name: str) -> str:
        """
        Return the parameter-encoded output file stem, e.g. ``dem_SLRM_R20``.

        :param dem_name: Stem of the DEM file name.
        :returns: The output file stem, without extension.
        """
        return f"{dem_name}_SLRM_R{self.radius_cells}"

    def outputs(self, dem_name: str) -> tuple[VisualizationOutput, ...]:
        """
        Return the single output this config produces.

        :param dem_name: Stem of the DEM file name.
        :returns: One output: the computed array, under
            :meth:`output_stem`, rendered with ``EIGHT_BIT_STRETCH``.
        """
        return (
            VisualizationOutput(
                None, self.output_stem(dem_name), self.EIGHT_BIT_STRETCH
            ),
        )

    def required_overlap(self, resolution: Resolution) -> int:
        """
        Return the tiling overlap in pixels: ``radius_cells + 1`` — the legacy
        integral-image window is shifted one pixel up-left of the kernel
        centre, so the neighbourhood extends one pixel past the radius.

        :param resolution: Cell size of the DEM (unused; the radius is in cells).
        :returns: The neighbourhood radius in pixels.
        """
        return self.radius_cells + 1

    def compute(
        self, dem: NDArray[np.float32], resolution: Resolution
    ) -> NDArray[np.float32]:
        """
        Compute the simple local relief model of a DEM.

        NaN propagates from the DEM's nodata pixels.

        :param dem: 2D float32 DEM with NaN as nodata. Not modified.
        :param resolution: Cell size of the DEM.
        :returns: Simple local relief model array of the DEM's shape.
        :raises ValueError: If the DEM is not 2D or the resolution is not positive.
        """
        validate_dem_and_resolution(dem.ndim, resolution)

        dem = dem.astype(np.float32)
        dem = dem * self.vertical_exaggeration

        dem_trend = _mean_filter(dem=dem, kernel_radius=self.radius_cells)
        return dem - dem_trend
