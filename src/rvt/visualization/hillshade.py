"""
Hillshade.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

import numpy as np
from numpy.typing import NDArray

from rvt.core import (
    Resolution,
    format_parameter,
    validate_dem_and_resolution,
    validate_vertical_exaggeration,
)
from rvt.stretch import Stretch, StretchMode
from rvt.visualization.slope import SlopeUnits, _slope_aspect


@dataclass(frozen=True)
class Hillshade:
    """
    Analytical hillshading from a single sun position.

    :param sun_azimuth: Solar azimuth angle in degrees, clockwise from north, 0–360.
    :param sun_elevation: Solar elevation angle above the horizon in degrees, 0–90.
    :param vertical_exaggeration: Factor the DEM heights are multiplied by.
    """

    EIGHT_BIT_STRETCH: ClassVar[Stretch] = Stretch(StretchMode.VALUE, 0.0, 1.0)
    TILEABLE: ClassVar[bool] = True

    sun_azimuth: float = 315.0
    sun_elevation: float = 35.0
    vertical_exaggeration: float = 1.0

    def __post_init__(self) -> None:
        validate_vertical_exaggeration(self.vertical_exaggeration)
        if not (0 <= self.sun_azimuth <= 360):
            raise ValueError(
                f"sun_azimuth must be between 0 and 360, got '{self.sun_azimuth}'!"
            )
        if not (0 <= self.sun_elevation <= 90):
            raise ValueError(
                f"sun_elevation must be between 0 and 90, got '{self.sun_elevation}'!"
            )

    def output_stem(self, dem_name: str) -> str:
        """
        Return the parameter-encoded output file stem, e.g. ``dem_HS_A315_H35``.

        :param dem_name: Stem of the DEM file name.
        :returns: The output file stem, without extension.
        """
        return (
            f"{dem_name}_HS_A{format_parameter(self.sun_azimuth)}"
            f"_H{format_parameter(self.sun_elevation)}"
        )

    def required_overlap(self, resolution: Resolution) -> int:
        """
        Return the tiling overlap in pixels: 1, the 3x3 derivative kernel's
        radius.

        :param resolution: Cell size of the DEM (unused; the kernel is fixed).
        :returns: The neighbourhood radius in pixels.
        """
        return 1

    def compute(
        self, dem: NDArray[np.float32], resolution: Resolution
    ) -> NDArray[np.float32]:
        """
        Compute the hillshade of a DEM.

        Values range 0–1 (negative incidence clamped to 0); NaN propagates from
        the DEM's nodata pixels.

        :param dem: 2D float32 DEM with NaN as nodata. Not modified.
        :param resolution: Cell size of the DEM.
        :returns: Hillshade array of the DEM's shape.
        :raises ValueError: If the DEM is not 2D or the resolution is not positive.
        """
        validate_dem_and_resolution(dem.ndim, resolution)

        # Legacy pipeline, ported verbatim for golden-master parity: the padded DEM
        # is exaggerated and slope/aspect are computed on it (padding again inside).
        dem = dem.astype(np.float32)
        dem = np.pad(array=dem, pad_width=1, mode="edge")
        dem = dem * self.vertical_exaggeration

        sun_azimuth_radians = np.deg2rad(self.sun_azimuth)
        sun_zenith_radians = np.pi / 2 - np.deg2rad(self.sun_elevation)

        slope, aspect = _slope_aspect(
            dem, resolution, SlopeUnits.RADIAN, vertical_exaggeration=1.0
        )

        hillshade = np.cos(sun_zenith_radians) * np.cos(slope) + np.sin(
            sun_zenith_radians
        ) * np.sin(slope) * np.cos(aspect - sun_azimuth_radians)
        hillshade[hillshade < 0] = 0
        return hillshade[1:-1, 1:-1]
