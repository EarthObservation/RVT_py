"""
Local dominance.
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


@dataclass(frozen=True)
class LocalDominance:
    """
    Local dominance: how much a pixel dominates the terrain around it.

    Adapted from the original version that is part of the Lidar Visualisation
    Toolbox LiVT developed by Ralf Hesse.

    :param minimum_radius: Minimum radial distance in pixels at which the
        computation starts.
    :param maximum_radius: Maximum radial distance in pixels at which the
        computation ends.
    :param radius_increment: Radial distance step in pixels.
    :param angular_resolution: Angular step in degrees determining the number
        of angular directions.
    :param observer_height: Height in map units at which the terrain is
        observed.
    :param vertical_exaggeration: Factor the DEM heights are multiplied by.
    """

    minimum_radius: int = 10
    maximum_radius: int = 20
    radius_increment: int = 1
    angular_resolution: int = 15
    observer_height: float = 1.7
    vertical_exaggeration: float = 1.0

    def __post_init__(self) -> None:
        validate_vertical_exaggeration(self.vertical_exaggeration)

    def compute(
        self, dem: NDArray[np.float32], resolution: Resolution
    ) -> NDArray[np.float32]:
        """
        Compute the local dominance of a DEM.

        NaN propagates from the DEM's nodata pixels.

        :param dem: 2D float32 DEM with NaN as nodata. Not modified.
        :param resolution: Cell size of the DEM.
        :returns: Local dominance array of the DEM's shape.
        :raises ValueError: If the DEM is not 2D or the resolution is not positive.
        """
        validate_dem_and_resolution(dem.ndim, resolution)

        # Legacy pipeline, ported verbatim for golden-master parity.
        dem = dem.astype(np.float32)
        pad_width = self.maximum_radius
        dem = np.pad(array=dem, pad_width=pad_width, mode="edge")
        dem = dem * self.vertical_exaggeration

        number_of_distances = int(
            (self.maximum_radius - self.minimum_radius) / self.radius_increment + 1
        )
        distances = (
            np.arange(
                number_of_distances * self.radius_increment,
                step=self.radius_increment,
            )
            + self.minimum_radius
        )
        number_of_angles = int(359 / self.angular_resolution + 1)
        angles = np.arange(
            number_of_angles * self.angular_resolution, step=self.angular_resolution
        )
        # Total area within the radius range, used for normalisation.
        normalisation = (
            np.sum(
                (self.observer_height / distances)
                * (2 * distances + self.radius_increment)
            )
            * number_of_angles
        )

        number_of_shifts = distances.size * angles.size
        shifts_x = (np.outer(np.cos(np.deg2rad(angles)), distances)).reshape(
            number_of_shifts
        )
        shifts_y = (np.outer(np.sin(np.deg2rad(angles)), distances)).reshape(
            number_of_shifts
        )
        distances = (np.outer(np.ones(number_of_angles), distances)).reshape(
            number_of_shifts
        )
        distance_factors = 2 * distances + self.radius_increment

        local_dominance = dem * 0
        for shift_index in range(number_of_shifts):
            dem_moved = np.roll(dem, int(round(shifts_y[shift_index])), axis=0)
            dem_moved = np.roll(dem_moved, int(round(shifts_x[shift_index])), axis=1)
            lower_indices = np.where((dem + self.observer_height) > dem_moved)
            if lower_indices[0].size > 0:
                local_dominance[lower_indices[0], lower_indices[1]] = (
                    local_dominance[lower_indices[0], lower_indices[1]]
                    + (
                        dem[lower_indices[0], lower_indices[1]]
                        + self.observer_height
                        - dem_moved[lower_indices[0], lower_indices[1]]
                    )
                    / distances[shift_index]
                    * distance_factors[shift_index]
                )
        # The legacy expression order promotes to float64; the final cast back to
        # float32 matches the float32 golden rasters bit-for-bit.
        local_dominance = local_dominance / normalisation

        return local_dominance[pad_width:-pad_width, pad_width:-pad_width].astype(
            np.float32
        )
