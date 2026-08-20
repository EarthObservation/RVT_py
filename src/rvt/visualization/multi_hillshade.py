"""
Multi-direction hillshade.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from rvt.core import (
    Resolution,
    validate_dem_and_resolution,
    validate_vertical_exaggeration,
)
from rvt.visualization.slope import SlopeUnits, _slope_aspect

FULL_CIRCLE_DEGREES = 360


@dataclass(frozen=True)
class MultiHillshade:
    """
    Hillshading from multiple evenly spaced sun azimuths.

    Azimuths start at 0 (north) and step by ``360 / directions`` clockwise.

    :param directions: Number of solar azimuth angles, at least 1.
    :param sun_elevation: Solar elevation angle above the horizon in degrees, 0-90.
    :param vertical_exaggeration: Factor the DEM heights are multiplied by.
    """

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

        return np.asarray(hillshades)
