"""
Cast shadow and horizon elevation angle from long-distance horizon search.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar, NamedTuple

import numpy as np
from numpy.typing import NDArray

from rvt.core import (
    Resolution,
    validate_dem_and_resolution,
    validate_vertical_exaggeration,
)
from rvt.stretch import Stretch, StretchMode
from rvt.visualization.sky_illumination import _compute_sky_illumination


class ShadowHorizonResult(NamedTuple):
    """
    Result of a shadow/horizon computation.

    :param shadow: Binary cast shadow, 1 where the terrain is lit (the horizon is
        below the sun elevation) and 0 where it is in shadow.
    :param horizon: Horizon elevation angle towards the sun azimuth in degrees.
    """

    shadow: NDArray[np.float32]
    horizon: NDArray[np.float32]


@dataclass(frozen=True)
class ShadowHorizon:
    """
    Binary cast shadow and horizon elevation angle for a single sun position.

    :param azimuth: Sun azimuth in degrees, clockwise from north, 0-360.
    :param elevation: Sun elevation above the horizon in degrees, 0-90.
    :param vertical_exaggeration: Factor the DEM heights are multiplied by.
    """

    EIGHT_BIT_STRETCH: ClassVar[Stretch] = Stretch(StretchMode.VALUE, 0.0, 1.0)

    azimuth: float = 315.0
    elevation: float = 35.0
    vertical_exaggeration: float = 1.0

    def __post_init__(self) -> None:
        validate_vertical_exaggeration(self.vertical_exaggeration)
        if not (0 <= self.azimuth <= 360):
            raise ValueError(
                f"azimuth must be between 0 and 360, got '{self.azimuth}'!"
            )
        if not (0 <= self.elevation <= 90):
            raise ValueError(
                f"elevation must be between 0 and 90, got '{self.elevation}'!"
            )

    def compute(
        self, dem: NDArray[np.float32], resolution: Resolution
    ) -> ShadowHorizonResult:
        """
        Compute the cast shadow and horizon angle of a DEM.

        A DEM containing any NaN yields an all-NaN horizon (the NaN poisons the
        horizon pyramid padding, a legacy behaviour kept for parity) and therefore
        an all-zero shadow.

        :param dem: 2D float32 DEM with NaN as nodata. Not modified.
        :param resolution: Cell size of the DEM.
        :returns: Shadow and horizon arrays of the DEM's shape.
        :raises ValueError: If the DEM is not 2D or the resolution is not positive.
        """
        validate_dem_and_resolution(dem.ndim, resolution)
        # The legacy shadow_horizon wrapper forwards its own defaults for the
        # remaining sky-illumination parameters (overcast model, 100 px, 32
        # directions); the early return makes the sky model irrelevant here.
        outputs = _compute_sky_illumination(
            dem,
            resolution,
            compute_overcast=True,
            compute_shadow=True,
            shadow_and_horizon_only=True,
            maximum_fine_radius=100,
            directions=32,
            shadow_azimuth=self.azimuth,
            shadow_elevation=self.elevation,
            vertical_exaggeration=self.vertical_exaggeration,
        )
        assert outputs.shadow is not None
        assert outputs.horizon is not None
        return ShadowHorizonResult(shadow=outputs.shadow, horizon=outputs.horizon)
