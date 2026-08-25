"""
Multi-scale topographic position.
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
    format_parameter,
    validate_dem_and_resolution,
    validate_vertical_exaggeration,
)
from rvt.stretch import Stretch, StretchMode


def _integral_image(
    array: NDArray[np.floating], data_type: type = np.float64
) -> NDArray[np.floating]:
    """Summed-area table with the origin in the upper-left corner."""
    return array.astype(data_type).cumsum(axis=0).cumsum(axis=1)


def _topographic_deviation(
    dem: NDArray[np.floating],
    pixel_counts_integral: NDArray[np.floating],
    dem_integral: NDArray[np.floating],
    dem_squared_integral: NDArray[np.floating],
    kernel_radius: int,
) -> NDArray[np.floating]:
    """
    Legacy deviation from mean elevation, DEV(D) = (z0 - mean) / std over a kernel
    of radius D, ported verbatim for golden-master parity.
    """
    radius_cells = int(kernel_radius)
    if radius_cells <= 0:
        return dem

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

    dem_mean = (
        np.roll(dem_integral, (radius_cells, radius_cells), axis=(0, 1))
        + np.roll(dem_integral, (-radius_cells - 1, -radius_cells - 1), axis=(0, 1))
        - np.roll(dem_integral, (-radius_cells - 1, radius_cells), axis=(0, 1))
        - np.roll(dem_integral, (radius_cells, -radius_cells - 1), axis=(0, 1))
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        dem_mean = dem_mean / kernel_pixel_counts

    dem_standard_deviation = (
        np.roll(dem_squared_integral, (radius_cells, radius_cells), axis=(0, 1))
        + np.roll(
            dem_squared_integral,
            (-radius_cells - 1, -radius_cells - 1),
            axis=(0, 1),
        )
        - np.roll(dem_squared_integral, (-radius_cells - 1, radius_cells), axis=(0, 1))
        - np.roll(dem_squared_integral, (radius_cells, -radius_cells - 1), axis=(0, 1))
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        dem_standard_deviation = np.sqrt(
            np.abs(dem_standard_deviation / kernel_pixel_counts - dem_mean**2)
        )

    # The 1e-6 prevents division by zero (legacy behaviour).
    return (np.roll(dem, (-1, -1), axis=(0, 1)) - dem_mean) / (
        dem_standard_deviation + 1e-6
    )


def _max_elevation_deviation(
    dem: NDArray[np.float32], minimum_radius: int, maximum_radius: int, step: int
) -> NDArray[np.float32]:
    """
    Legacy maximum deviation from mean elevation across a range of kernel radii,
    ported verbatim for golden-master parity.
    """
    minimum_radius = int(minimum_radius)
    maximum_radius = int(maximum_radius)
    step = int(step)

    dem_nan_mask = np.isnan(dem)

    dem_padded = np.pad(dem, (maximum_radius + 1, maximum_radius), mode="symmetric")
    padded_nan_mask = np.isnan(dem_padded)
    dem_padded[padded_nan_mask] = 0

    pixel_counts = np.ones(dem_padded.shape)
    pixel_counts[padded_nan_mask] = 0
    pixel_counts_integral = _integral_image(pixel_counts, np.int64)

    # Float64 integral images by design; the result is cast to float32 at return.
    dem_integral = _integral_image(dem_padded)
    dem_squared_integral = _integral_image(dem_padded**2)

    maximum_deviation = np.empty(0)
    for kernel_radius in range(minimum_radius, maximum_radius + 1, step):
        deviation = _topographic_deviation(
            dem_padded,
            pixel_counts_integral,
            dem_integral,
            dem_squared_integral,
            kernel_radius,
        )[
            maximum_radius : -(maximum_radius + 1),
            maximum_radius : -(maximum_radius + 1),
        ]
        if kernel_radius == minimum_radius:
            maximum_deviation = deviation
        else:
            maximum_deviation = np.where(
                np.abs(maximum_deviation) >= np.abs(deviation),
                maximum_deviation,
                deviation,
            )

    maximum_deviation[dem_nan_mask] = np.nan

    return maximum_deviation.astype(np.float32)


def _validate_scale(name: str, scale: tuple[int, int, int]) -> None:
    minimum_radius, maximum_radius, step = scale
    if minimum_radius > maximum_radius:
        raise ValueError(
            f"{name} minimum radius has to be smaller than maximum radius,"
            f" got '{scale}'!"
        )
    if maximum_radius - minimum_radius < step:
        raise ValueError(
            f"{name} step has to be within minimum and maximum radius, got '{scale}'!"
        )


@dataclass(frozen=True)
class MultiScaleTopographicPosition(Visualization[NDArray[np.float32]]):
    """
    Multi-scale topographic position: an RGB composite of the maximum deviation
    from mean elevation at broad (red), meso (green) and local (blue) scales.

    :param local_scale: Local scale (minimum radius, maximum radius, step) in pixels.
    :param meso_scale: Meso scale (minimum radius, maximum radius, step) in pixels.
    :param broad_scale: Broad scale (minimum radius, maximum radius, step) in pixels.
    :param lightness: Lightness of the image.
    :param vertical_exaggeration: Factor the DEM heights are multiplied by.
    """

    EIGHT_BIT_STRETCH: ClassVar[Stretch] = Stretch(StretchMode.VALUE, 0.0, 1.0)
    TILEABLE: ClassVar[bool] = True

    local_scale: tuple[int, int, int] = (1, 5, 1)
    meso_scale: tuple[int, int, int] = (5, 50, 5)
    broad_scale: tuple[int, int, int] = (50, 500, 50)
    lightness: float = 1.2
    vertical_exaggeration: float = 1.0

    def __post_init__(self) -> None:
        validate_vertical_exaggeration(self.vertical_exaggeration)
        _validate_scale("local_scale", self.local_scale)
        _validate_scale("meso_scale", self.meso_scale)
        _validate_scale("broad_scale", self.broad_scale)

    def output_stem(self, dem_name: str) -> str:
        """
        Return the parameter-encoded output file stem — the maximum radius of
        each scale plus the lightness, e.g. ``dem_MSTP_5_50_500_L1.2``.

        :param dem_name: Stem of the DEM file name.
        :returns: The output file stem, without extension.
        """
        return (
            f"{dem_name}_MSTP_{self.local_scale[1]}_{self.meso_scale[1]}"
            f"_{self.broad_scale[1]}_L{format_parameter(self.lightness)}"
        )

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
        Return the tiling overlap in pixels: the broad scale's maximum radius
        plus one (the integral-image window sits one pixel up-left of the
        kernel centre).

        With the default broad scale this is 501 pixels, so tiles smaller than
        that read most of the raster anyway — tiling pays off only for rasters
        much larger than the overlap.

        :param resolution: Cell size of the DEM (unused; the scales are in pixels).
        :returns: The neighbourhood radius in pixels.
        """
        return self.broad_scale[1] + 1

    def compute(
        self, dem: NDArray[np.float32], resolution: Resolution
    ) -> NDArray[np.float32]:
        """
        Compute the multi-scale topographic position of a DEM.

        NaN propagates from the DEM's nodata pixels.

        :param dem: 2D float32 DEM with NaN as nodata. Not modified.
        :param resolution: Cell size of the DEM.
        :returns: RGB float array of shape ``(3, rows, columns)`` with values 0-1.
        :raises ValueError: If the DEM is not 2D or the resolution is not positive.
        """
        validate_dem_and_resolution(dem.ndim, resolution)

        dem = dem.astype(np.float32)
        dem = dem * self.vertical_exaggeration

        local_deviation = _max_elevation_deviation(
            dem,
            minimum_radius=self.local_scale[0],
            maximum_radius=self.local_scale[1],
            step=self.local_scale[2],
        )
        meso_deviation = _max_elevation_deviation(
            dem,
            minimum_radius=self.meso_scale[0],
            maximum_radius=self.meso_scale[1],
            step=self.meso_scale[2],
        )
        broad_deviation = _max_elevation_deviation(
            dem,
            minimum_radius=self.broad_scale[0],
            maximum_radius=self.broad_scale[1],
            step=self.broad_scale[2],
        )

        # RGB order - broad, meso, local (legacy behaviour).
        red = 1 - np.exp(-self.lightness * np.abs(broad_deviation))
        green = 1 - np.exp(-self.lightness * np.abs(meso_deviation))
        blue = 1 - np.exp(-self.lightness * np.abs(local_deviation))

        red[red < 0] = 0
        green[green < 0] = 0
        blue[blue < 0] = 0

        red[red > 1] = 1
        green[green > 1] = 1
        blue[blue > 1] = 1

        return np.asarray([red, green, blue])
