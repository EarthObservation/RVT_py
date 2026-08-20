"""
Shared core types for the RVT library.
"""

from __future__ import annotations

from typing import NamedTuple

MAXIMUM_VERTICAL_EXAGGERATION = 10000.0


class Resolution(NamedTuple):
    """
    Raster cell size in map units.

    :param x: Cell size in the X direction (always positive).
    :param y: Cell size in the Y direction (always positive).
    """

    x: float
    y: float


def format_parameter(value: float) -> str:
    """
    Format a numeric parameter for a file name: integral values render without
    a decimal part (``315.0`` becomes ``315``), others as-is (``22.5``).

    :param value: The parameter value to encode.
    :returns: The compact string form of the value.
    """
    if value == int(value):
        return str(int(value))
    return str(value)


def validate_vertical_exaggeration(vertical_exaggeration: float) -> None:
    """
    Raise ``ValueError`` when the exaggeration factor is outside ±10000.

    :param vertical_exaggeration: Vertical exaggeration factor applied to the DEM.
    :raises ValueError: If the value is outside the supported range.
    """
    if not (
        -MAXIMUM_VERTICAL_EXAGGERATION
        <= vertical_exaggeration
        <= MAXIMUM_VERTICAL_EXAGGERATION
    ):
        raise ValueError(
            f"vertical_exaggeration must be between -{MAXIMUM_VERTICAL_EXAGGERATION:g}"
            f" and {MAXIMUM_VERTICAL_EXAGGERATION:g}, got '{vertical_exaggeration}'!"
        )


def validate_dem_and_resolution(dem_ndim: int, resolution: Resolution) -> None:
    """
    Raise ``ValueError`` for a non-2D DEM or a non-positive resolution.

    :param dem_ndim: Number of dimensions of the DEM array.
    :param resolution: Cell size of the DEM.
    :raises ValueError: If the DEM is not 2D or the resolution is not positive.
    """
    if dem_ndim != 2:
        raise ValueError(f"dem has to be a 2D array, got {dem_ndim} dimensions!")
    if resolution.x <= 0 or resolution.y <= 0:
        raise ValueError(f"resolution must be positive, got '{resolution}'!")
