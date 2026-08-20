"""
Conversion of float visualizations to 8-bit via typed stretch specifications.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum, auto

import numpy as np
from numpy.typing import NDArray


class StretchMode(StrEnum):
    """
    How the stretch bounds are interpreted: absolute values, or percentile cutoffs.
    """

    VALUE = auto()
    PERCENT = auto()


@dataclass(frozen=True)
class Stretch:
    """
    Mapping of a float visualization onto the 0–255 range.

    :param mode: Whether ``minimum``/``maximum`` are absolute values or percentile
        cutoffs (percent of pixels cut from each end).
    :param minimum: Lower bound (value) or lower cutoff percentage.
    :param maximum: Upper bound (value) or upper cutoff percentage.
    :param invert: Invert the scale after normalization (legacy behaviour for slope
        gradient and negative openness, where high values render dark).
    """

    mode: StretchMode
    minimum: float
    maximum: float
    invert: bool = False

    def __post_init__(self) -> None:
        if self.mode == StretchMode.VALUE:
            if self.minimum >= self.maximum:
                raise ValueError(
                    f"minimum must be smaller than maximum, got '{self.minimum}' and '{self.maximum}'!"
                )
        else:
            if (
                self.minimum < 0
                or self.maximum < 0
                or self.minimum > 100
                or self.maximum > 100
            ):
                raise ValueError(
                    f"percent cutoffs must be in range 0-100, got '{self.minimum}' and '{self.maximum}'!"
                )
            if self.minimum + self.maximum > 100:
                raise ValueError(
                    f"percent cutoffs must sum below 100, got '{self.minimum}' and '{self.maximum}'!"
                )


def _normalize_linear(
    image: NDArray[np.float32],
    minimum: float | np.floating,
    maximum: float | np.floating,
) -> NDArray[np.float32]:
    """Legacy linear cut-off and stretch to 0-1, ported verbatim (made pure)."""
    image = image.copy()
    image[image > maximum] = maximum
    image[image < minimum] = minimum
    image = (image - minimum) / (maximum - minimum)
    image[image > 1] = 1
    image[image < 0] = 0
    return np.float32(image)


def _normalize_percent(
    image: NDArray[np.float32], minimum: float, maximum: float
) -> NDArray[np.float32]:
    """Legacy percentile cut-off, ported verbatim: falls back to min/max when the percentiles coincide."""
    distribution = np.nanpercentile(a=image, q=np.array([minimum, 100 - maximum]))
    minimum_linear = distribution[0]
    maximum_linear = distribution[1]
    if minimum_linear == maximum_linear:
        minimum_linear = np.nanmin(image)
        maximum_linear = np.nanmax(image)
    # The numpy float64 scalars are passed through unconverted: legacy promoted the
    # stretch arithmetic to float64 here, and parity requires the same rounding.
    return _normalize_linear(image, minimum_linear, maximum_linear)


def _byte_scale(
    data: NDArray[np.float32], c_min: float, c_max: float
) -> NDArray[np.uint8]:
    """Legacy IDL BYTSCL port: linear scale to 0-255, NaN encoded as 255."""
    high = 255
    byte_data = (high + 0.9999) * (data - c_min) / (c_max - c_min)
    byte_data[byte_data > high] = high
    byte_data[byte_data < 0] = 0
    byte_data[np.isnan(byte_data)] = 255
    return np.asarray(byte_data, dtype=np.uint8)


def to_8bit(array: NDArray[np.float32], stretch: Stretch) -> NDArray[np.uint8]:
    """
    Convert a float visualization to 8-bit: normalize with the stretch, then byte
    scale 0–255 with NaN encoded as 255 (legacy semantics).

    :param array: 2D float visualization with NaN as nodata. Not modified.
    :param stretch: The stretch specification; each visualization config carries a
        suitable default as ``EIGHT_BIT_STRETCH``.
    :returns: The uint8 array.
    """
    if stretch.mode == StretchMode.VALUE:
        normalized = _normalize_linear(array, stretch.minimum, stretch.maximum)
    else:
        normalized = _normalize_percent(array, stretch.minimum, stretch.maximum)
    if stretch.invert:
        normalized = 1 - normalized
    return _byte_scale(normalized, c_min=0, c_max=1)
