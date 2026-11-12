"""
Multi-direction hillshade.
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
from rvt.visualization.hillshade import Hillshade
from rvt.visualization.slope import SlopeUnits, _slope_aspect

FULL_CIRCLE_DEGREES = 360


@dataclass(frozen=True)
class MultiHillshade(Visualization[NDArray[np.float32]]):
    """
    Hillshading from multiple evenly spaced sun azimuths.

    Azimuths start at 0 (north) and step by ``360 / directions`` clockwise.

    :param directions: Number of solar azimuth angles, at least 1.
    :param sun_elevation: Solar elevation angle above the horizon in degrees, 0-90.
    :param vertical_exaggeration: Factor the DEM heights are multiplied by.

    The 8-bit rendering is not a stretch of the float stack: like version 2,
    it is a 3-band RGB of single hillshades at ``RGB_AZIMUTHS``, each band
    stretched with ``EIGHT_BIT_STRETCH``.
    """

    RGB_AZIMUTHS: ClassVar[tuple[float, float, float]] = (315.0, 22.5, 90.0)
    EIGHT_BIT_STRETCH: ClassVar[Stretch] = Stretch(StretchMode.VALUE, 0.0, 1.0)
    # Not tileable: the legacy pipeline trims a 1-pixel border, so the output
    # shape is (directions, rows - 2, cols - 2) — it cannot be written back
    # into fixed windows of an input-shaped raster.
    TILEABLE: ClassVar[bool] = False

    directions: int = 16
    sun_elevation: float = 35.0
    vertical_exaggeration: float = 1.0

    def __post_init__(self) -> None:
        validate_vertical_exaggeration(self.vertical_exaggeration)
        if self.directions < 1:
            raise ValueError(
                f"directions must be a positive number, got '{self.directions}'!"
            )
        if not (0 <= self.sun_elevation <= 90):
            raise ValueError(
                f"sun_elevation must be between 0 and 90, got '{self.sun_elevation}'!"
            )

    def output_stem(self, dem_name: str) -> str:
        """
        Return the parameter-encoded output file stem, e.g. ``dem_MULTI-HS_D16_H35``.

        :param dem_name: Stem of the DEM file name.
        :returns: The output file stem, without extension.
        """
        return (
            f"{dem_name}_MULTI-HS_D{self.directions}"
            f"_H{format_parameter(self.sun_elevation)}"
        )

    def outputs(self, dem_name: str) -> tuple[VisualizationOutput, ...]:
        """
        Return the single output this config produces.

        :param dem_name: Stem of the DEM file name.
        :returns: One output: the band stack under :meth:`output_stem`, whose
            8-bit form is the RGB rendering of :meth:`eight_bit_image`.
        """
        return (
            VisualizationOutput(
                None, self.output_stem(dem_name), self.EIGHT_BIT_STRETCH
            ),
        )

    def rgb_image(
        self, dem: NDArray[np.float32], resolution: Resolution
    ) -> NDArray[np.float32]:
        """
        Compute the bands-first RGB stack of hillshades at ``RGB_AZIMUTHS``.

        :param dem: 2D float32 DEM with NaN as nodata. Not modified.
        :param resolution: Cell size of the DEM.
        :returns: A ``(3, rows, cols)`` stack of hillshades, values 0-1.
        """
        return np.array(
            [
                Hillshade(
                    sun_azimuth=azimuth,
                    sun_elevation=self.sun_elevation,
                    vertical_exaggeration=self.vertical_exaggeration,
                ).compute(dem, resolution)
                for azimuth in self.RGB_AZIMUTHS
            ]
        )

    def eight_bit_image(
        self,
        image: NDArray[np.float32],
        dem: NDArray[np.float32],
        resolution: Resolution,
    ) -> NDArray[np.float32]:
        """
        Return the RGB hillshade stack from :meth:`rgb_image`, which the 8-bit
        rendering stretches instead of the float stack.

        :param image: The computed float stack (unused).
        :param dem: The DEM it was computed from.
        :param resolution: Cell size of the DEM.
        :returns: The ``(3, rows, cols)`` RGB hillshade stack.
        """
        return self.rgb_image(dem, resolution)

    def required_overlap(self, resolution: Resolution) -> int:
        """
        Return the tiling overlap in pixels: 1, the 3x3 derivative kernel's
        radius (informational — the config is not tileable, see ``TILEABLE``).

        :param resolution: Cell size of the DEM (unused; the kernel is fixed).
        :returns: The neighbourhood radius in pixels.
        """
        return 1

    def compute(
        self, dem: NDArray[np.float32], resolution: Resolution
    ) -> NDArray[np.float32]:
        """
        Compute the multi-direction hillshade stack of a DEM.

        Values range 0-1 (negative incidence clamped to 0); NaN propagates from
        the DEM's nodata pixels.

        :param dem: 2D float32 DEM with NaN as nodata. Not modified.
        :param resolution: Cell size of the DEM.
        :returns: Stack of shape ``(directions, rows - 2, cols - 2)`` — the legacy
            pipeline trims a 1-pixel border, kept for golden-master parity.
        :raises ValueError: If the DEM is not 2D or the resolution is not positive.
        """
        validate_dem_and_resolution(dem.ndim, resolution)

        # Legacy pipeline, ported verbatim for golden-master parity: slope/aspect
        # are computed once on the exaggerated (unpadded) DEM, then each azimuth's
        # hillshade is shaded on them and trimmed by one pixel on every edge.
        dem = dem.astype(np.float32)
        dem = dem * self.vertical_exaggeration

        slope, aspect = _slope_aspect(
            dem, resolution, SlopeUnits.RADIAN, vertical_exaggeration=1.0
        )

        sun_zenith_radians = np.pi / 2 - np.deg2rad(self.sun_elevation)

        hillshades = []
        for direction_index in range(self.directions):
            sun_azimuth = (FULL_CIRCLE_DEGREES / self.directions) * direction_index
            sun_azimuth_radians = np.deg2rad(sun_azimuth)
            hillshade = np.cos(sun_zenith_radians) * np.cos(slope) + np.sin(
                sun_zenith_radians
            ) * np.sin(slope) * np.cos(aspect - sun_azimuth_radians)
            hillshade[hillshade < 0] = 0
            hillshades.append(hillshade[1:-1, 1:-1])

        # The solar-position scalars promote the arithmetic to float64; the cast
        # back matches the float32 the legacy pipeline wrote to disk, and keeps
        # the compute layer's float32 contract.
        return np.asarray(hillshades, dtype=np.float32)
