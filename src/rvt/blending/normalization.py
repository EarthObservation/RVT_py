"""
Normalization of visualization images onto the 0-1 range before blending.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from enum import StrEnum, auto

import numpy as np
from numpy.typing import NDArray


class NormalizationMode(StrEnum):
    """
    How the normalization bounds are interpreted: absolute values, or
    percentile cutoffs (percent of pixels cut from each end).
    """

    VALUE = auto()
    PERCENT = auto()


@dataclass(frozen=True)
class Normalization:
    """
    Linear normalization of a visualization onto the 0-1 range.

    :param mode: Whether ``minimum``/``maximum`` are absolute values or
        percentile cutoffs.
    :param minimum: Lower bound (value) or lower cutoff percentage.
    :param maximum: Upper bound (value) or upper cutoff percentage.
    """

    mode: NormalizationMode
    minimum: float
    maximum: float

    def __post_init__(self) -> None:
        if self.mode == NormalizationMode.VALUE:
            if self.minimum >= self.maximum:
                raise ValueError(
                    f"minimum must be smaller than maximum, got '{self.minimum}'"
                    f" and '{self.maximum}'!"
                )
        else:
            if (
                self.minimum < 0
                or self.maximum < 0
                or self.minimum > 100
                or self.maximum > 100
            ):
                raise ValueError(
                    f"percent cutoffs must be in range 0-100, got '{self.minimum}'"
                    f" and '{self.maximum}'!"
                )
            if self.minimum + self.maximum > 100:
                raise ValueError(
                    f"percent cutoffs must sum below 100, got '{self.minimum}'"
                    f" and '{self.maximum}'!"
                )

    def apply(self, image: NDArray[np.floating]) -> NDArray[np.float32]:
        """
        Normalize ``image`` to 0-1. Pure: the input array is not modified.

        :param image: Visualization image (single band, or bands-first RGB).
        :returns: The normalized float32 array.
        """
        if self.mode == NormalizationMode.VALUE:
            return _normalize_linear(image, self.minimum, self.maximum)
        return _normalize_percent(image, self.minimum, self.maximum)


def _normalize_linear(
    image: NDArray[np.floating],
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
    return image.astype(np.float32)


def _normalize_percent(
    image: NDArray[np.floating], minimum: float, maximum: float
) -> NDArray[np.float32]:
    """Legacy percentile cut-off, ported verbatim: falls back to min/max when the percentiles coincide."""
    distribution = np.nanpercentile(a=image, q=np.array([minimum, 100 - maximum]))
    minimum_linear = distribution[0]
    maximum_linear = distribution[1]
    if minimum_linear == maximum_linear:
        minimum_linear = np.nanmin(image)
        maximum_linear = np.nanmax(image)
    # The numpy float64 scalars are passed through unconverted: legacy promoted
    # the stretch arithmetic to float64 here, and parity requires the same rounding.
    return _normalize_linear(image, minimum_linear, maximum_linear)


def normalize_image(
    image: NDArray[np.floating],
    normalization: Normalization,
    invert: bool = False,
) -> NDArray[np.float32]:
    """
    Normalize a layer image for blending: apply the normalization, rescue the
    scale when the result leaves 0-1 (legacy heuristic, with a warning), and
    optionally invert (legacy identity for slope gradient and negative
    openness, where high values render dark).

    :param image: Visualization image (single band, or bands-first RGB). Not
        modified.
    :param normalization: The normalization specification.
    :param invert: Invert the scale after normalization.
    :returns: The normalized float32 array.
    """
    normalized = normalization.apply(image)

    if np.nanmax(normalized) > 1:
        normalized = scale_0_to_1(normalized)
        warnings.warn("rvt.blending.normalize_image: unexpected values! max > 1")
    if np.nanmin(normalized) < 0:
        normalized = scale_0_to_1(normalized)
        warnings.warn("rvt.blending.normalize_image: unexpected values! min < 0")

    if invert:
        normalized = 1 - normalized
    return normalized


def scale_0_to_1(image: NDArray[np.floating]) -> NDArray[np.floating]:
    """
    Legacy heuristic rescue of an image onto the 0-1 range, ported verbatim
    (made pure): already-scaled images pass through, wide-range images are
    scaled within 0-1, narrow-range images are stretched strictly to 0-1.

    :param image: Image to rescale. Not modified.
    :returns: The rescaled array.
    """
    if 1 >= np.nanmax(image) > 0.9 and np.nanmin(image) == 0:
        return image
    if np.nanmax(image) - np.nanmin(image) > 0.3:
        return _scale_within_0_and_1(image)
    return _scale_strict_0_to_1(image)


def _scale_within_0_and_1(image: NDArray[np.floating]) -> NDArray[np.floating]:
    """Legacy ``scale_within_0_and_1``, ported verbatim (made pure)."""
    if np.nanmin(image) >= 0 and np.nanmax(image) <= 1:
        return image

    image = image.copy()
    image[np.isnan(image)] = np.nanmin(image)

    actual_minimum = np.nanmin(image)
    actual_maximum = np.nanmax(image)

    # 8-bit-looking images are scaled by 255 instead of min/max.
    if (
        255 >= actual_maximum > 1
        and actual_maximum - actual_minimum > 30
        and actual_minimum > 0
    ):
        return image / 255

    # The legacy bounds computation `np.nanmax(np.array(0, actual_min))` treats
    # the second argument as a dtype, so the bounds are the constants 0 and 1.
    scaled = (image - 0.0) / (1.0 - 0.0)

    if np.nanmin(scaled) > -0.01:
        scaled[(0 > scaled) & (scaled > -0.01)] = 0
    return scaled


def _scale_strict_0_to_1(image: NDArray[np.floating]) -> NDArray[np.floating]:
    """Legacy ``scale_strict_0_to_1``, ported verbatim (made pure)."""
    if np.nanmin(image) == 0 and np.nanmax(image) == 1:
        return image

    image = image.copy()
    image[np.isnan(image)] = 0

    minimum = np.nanmin(image)
    maximum = np.nanmax(image)
    scaled = (image - minimum) / (maximum - minimum)

    # Legacy wrote `scaled[0 > scaled > -0.01]`, which raises on arrays; the
    # intended elementwise snap of tiny negatives to 0 is a no-op here anyway.
    if np.nanmin(scaled) > -0.01:
        scaled[(0 > scaled) & (scaled > -0.01)] = 0
    return scaled
