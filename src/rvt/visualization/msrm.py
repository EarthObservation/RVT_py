"""
Multi-scale relief model.
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
from rvt.visualization.slrm import _mean_filter


@dataclass(frozen=True)
class MultiScaleReliefModel:
    """
    Multi-scale relief model (MSRM): the mean of differences between
    consecutive mean-filtered surfaces over a range of feature sizes.

    :param feature_minimum: Minimum size of the features to detect, in map
        units. Values below the cell size are raised to the cell size.
    :param feature_maximum: Maximum size of the features to detect, in map
        units.
    :param scaling_factor: Scaling factor; values larger than 1 provide a
        larger range of output values (more contrast and visibility) but can
        lose sensitivity for intermediate-sized features.
    :param vertical_exaggeration: Factor the DEM heights are multiplied by.
    """

    EIGHT_BIT_STRETCH: ClassVar[Stretch] = Stretch(StretchMode.VALUE, -2.5, 2.5)

    feature_minimum: float
    feature_maximum: float
    scaling_factor: int
    vertical_exaggeration: float = 1.0

    def __post_init__(self) -> None:
        validate_vertical_exaggeration(self.vertical_exaggeration)

    def output_stem(self, dem_name: str) -> str:
        """
        Return the parameter-encoded output file stem, e.g. ``dem_MSRM_F_M0-20_S2``.

        :param dem_name: Stem of the DEM file name.
        :returns: The output file stem, without extension.
        """
        return (
            f"{dem_name}_MSRM_F_M{format_parameter(self.feature_minimum)}"
            f"-{format_parameter(self.feature_maximum)}_S{self.scaling_factor}"
        )

    def compute(
        self, dem: NDArray[np.float32], resolution: Resolution
    ) -> NDArray[np.float32]:
        """
        Compute the multi-scale relief model of a DEM.

        NaN propagates from the DEM's nodata pixels.

        :param dem: 2D float32 DEM with NaN as nodata. Not modified.
        :param resolution: Cell size of the DEM.
        :returns: Multi-scale relief model array of the DEM's shape.
        :raises ValueError: If the DEM is not 2D or the resolution is not positive.
        """
        validate_dem_and_resolution(dem.ndim, resolution)

        # Legacy pipeline, ported verbatim for golden-master parity. The legacy
        # function takes a single scalar resolution; the X cell size is used.
        cell_size = resolution.x

        dem = dem.astype(np.float32)
        dem = dem * self.vertical_exaggeration

        feature_minimum = self.feature_minimum
        if feature_minimum < cell_size:
            feature_minimum = cell_size

        scaling_factor = int(self.scaling_factor)

        # i and n from the MSRM article.
        minimum_step = int(
            np.floor(
                ((feature_minimum - cell_size) / (2 * cell_size))
                ** (1 / scaling_factor)
            )
        )
        maximum_step = int(
            np.ceil(
                ((self.feature_maximum - cell_size) / (2 * cell_size))
                ** (1 / scaling_factor)
            )
        )

        relief_models_sum = np.zeros(dem.shape)
        number_of_relief_models = 0
        last_filtered_surface = 0

        for step in range(minimum_step, maximum_step + 1, 1):
            kernel_radius = step**scaling_factor
            filtered_surface = _mean_filter(dem=dem, kernel_radius=kernel_radius)
            if not step == minimum_step:
                # Difference of two consecutive low-pass-filtered surfaces.
                relief_models_sum += last_filtered_surface - filtered_surface
                number_of_relief_models += 1
            last_filtered_surface = filtered_surface

        # The legacy accumulator is float64; the final cast back to float32
        # matches the float32 golden rasters bit-for-bit.
        return (relief_models_sum / number_of_relief_models).astype(np.float32)
