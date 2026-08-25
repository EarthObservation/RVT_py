"""
Slope gradient and aspect.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum, auto
from typing import ClassVar, NamedTuple

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


class SlopeUnits(StrEnum):
    """
    Output units for slope: radians (default), degrees, or percent.
    """

    RADIAN = auto()
    DEGREE = auto()
    PERCENT = auto()


class SlopeAspectResult(NamedTuple):
    """
    Result of a slope/aspect computation.

    :param slope: Slope gradient, 0 for a horizontal plane, in the requested units.
    :param aspect: Downslope direction as geographic azimuth (0 north, clockwise),
        in radians unless degrees were requested.
    """

    slope: NDArray[np.float32]
    aspect: NDArray[np.float32]


def _roll_fill_nans(
    dem: NDArray[np.float32], shift: int, axis: int
) -> NDArray[np.float32]:
    """Roll the array, replacing NaNs the roll introduced with the source values (edge padding)."""
    rolled = np.roll(dem, shift, axis=axis)
    introduced_nan = np.isnan(rolled) != np.isnan(dem)
    rolled[introduced_nan] = dem[introduced_nan]
    return rolled


def _slope_aspect(
    dem: NDArray[np.float32],
    resolution: Resolution,
    units: SlopeUnits,
    vertical_exaggeration: float,
) -> SlopeAspectResult:
    """
    Legacy slope/aspect math, ported verbatim for golden-master parity.
    """
    dem = dem.astype(np.float32)
    nan_mask = np.isnan(dem)
    dem = np.pad(array=dem, pad_width=1, mode="edge")
    dem = dem * vertical_exaggeration

    dzdx = (
        (_roll_fill_nans(dem, 1, axis=1) - _roll_fill_nans(dem, -1, axis=1)) / 2
    ) / resolution.x
    dzdy = (
        (_roll_fill_nans(dem, -1, axis=0) - _roll_fill_nans(dem, 1, axis=0)) / 2
    ) / resolution.y
    tan_slope = np.sqrt(dzdx**2 + dzdy**2)

    match units:
        case SlopeUnits.PERCENT:
            slope = tan_slope * 100
        case SlopeUnits.DEGREE:
            slope = np.rad2deg(np.arctan(tan_slope))
        case SlopeUnits.RADIAN:
            slope = np.arctan(tan_slope)

    # Where dzdy is zero nudge it, for numeric stability of the quadrant computation.
    dzdy[dzdy == 0] = 10e-9
    aspect = np.arctan2(dzdx, dzdy)
    if units == SlopeUnits.DEGREE:
        aspect = np.rad2deg(aspect)

    slope = slope[1:-1, 1:-1]
    aspect = aspect[1:-1, 1:-1]
    slope[nan_mask] = np.nan
    aspect[nan_mask] = np.nan
    return SlopeAspectResult(slope=slope, aspect=aspect)


@dataclass(frozen=True)
class SlopeAspect(Visualization[SlopeAspectResult]):
    """
    Terrain slope gradient and aspect via finite differences.

    :param units: Output units for slope (aspect converts only for degrees).
    :param vertical_exaggeration: Factor the DEM heights are multiplied by.
    """

    # Legacy 8-bit slope is stretched over degrees 0-51 and inverted (steep is dark).
    EIGHT_BIT_STRETCH: ClassVar[Stretch] = Stretch(
        StretchMode.VALUE, 0.0, 51.0, invert=True
    )
    TILEABLE: ClassVar[bool] = True

    units: SlopeUnits = SlopeUnits.RADIAN
    vertical_exaggeration: float = 1.0

    def __post_init__(self) -> None:
        validate_vertical_exaggeration(self.vertical_exaggeration)

    def outputs(self, dem_name: str) -> tuple[VisualizationOutput, ...]:
        """
        Return the two outputs this config produces: ``dem_SLOPE`` and
        ``dem_ASPECT``. Aspect has no legacy 8-bit form, so it is float-only.

        :param dem_name: Stem of the DEM file name.
        :returns: The slope and aspect outputs, in that order.
        """
        return (
            VisualizationOutput("slope", f"{dem_name}_SLOPE", self.EIGHT_BIT_STRETCH),
            VisualizationOutput("aspect", f"{dem_name}_ASPECT", None),
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
    ) -> SlopeAspectResult:
        """
        Compute slope and aspect of a DEM.

        :param dem: 2D float32 DEM with NaN as nodata. Not modified.
        :param resolution: Cell size of the DEM.
        :returns: Slope and aspect arrays with NaN where the DEM had NaN.
        :raises ValueError: If the DEM is not 2D or the resolution is not positive.
        """
        validate_dem_and_resolution(dem.ndim, resolution)
        return _slope_aspect(dem, resolution, self.units, self.vertical_exaggeration)
