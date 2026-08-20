"""
Sky illumination: topographic correction under a uniform or overcast sky model,
with an optional cast-shadow component, via long-distance horizon search on a
max-pooling DEM pyramid.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum, auto
from typing import Any, ClassVar, NamedTuple

import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import RectBivariateSpline

from rvt.core import (
    Resolution,
    validate_dem_and_resolution,
    validate_vertical_exaggeration,
)
from rvt.stretch import Stretch, StretchMode
from rvt.visualization.slope import SlopeUnits, _slope_aspect

_PYRAMID_SCALE = 2
_MAXIMUM_PYRAMID_RADIUS = 20


class SkyModel(StrEnum):
    """
    Sky luminance distribution model.
    """

    OVERCAST = auto()
    UNIFORM = auto()


class _DirectionSearch(NamedTuple):
    """Horizon search offsets for one azimuth: np.roll shifts and their distances."""

    azimuth: float
    shifts: list[tuple[int, int]]
    distances: NDArray[np.float64]


class _PyramidLevel(NamedTuple):
    """One DEM pyramid level with its per-direction horizon search offsets."""

    dem: NDArray[Any]
    searches: list[_DirectionSearch]
    line_indices: NDArray[Any]
    column_indices: NDArray[Any]


class _SkyIlluminationOutputs(NamedTuple):
    """Internal outputs of the shared sky-illumination core."""

    illumination: NDArray[np.float32] | None
    shadow: NDArray[np.float32] | None
    horizon: NDArray[np.float32] | None


def _horizon_shift_vector(
    num_directions: int, radius_pixels: float, min_radius: float
) -> list[_DirectionSearch]:
    """Legacy horizon search-movement construction, ported verbatim for parity."""
    angles = (2 * np.pi / num_directions) * np.arange(num_directions)
    x = np.cos(angles)
    y = np.sin(angles)
    # The legacy code keys the shifts by the azimuth rounded to one decimal and
    # later derives the direction angle from that rounded value; keep both quirks.
    rounded_azimuths = np.round(np.degrees(angles), decimals=1)

    # A finer range of radii; the added constant avoids starting at (0, 0).
    scale = 3.0
    radii = np.arange((radius_pixels - min_radius) * scale + 1) / scale + min_radius

    searches = []
    for i in range(num_directions):
        x_int = np.round(x[i] * radii, decimals=0)
        y_int = np.round(y[i] * radii, decimals=0)
        # Deduplicate the rounded points via a set of complex numbers.
        coord_complex = set(x_int + 1j * y_int)
        shift_pairs = np.array([(k.real, k.imag) for k in coord_complex]).astype(int)
        distance = np.sqrt(np.sum(shift_pairs**2, axis=1))
        sort_index = np.argsort(distance)
        searches.append(
            _DirectionSearch(
                azimuth=float(rounded_azimuths[i]),
                shifts=[(k[0], k[1]) for k in shift_pairs[sort_index]],
                distances=distance[sort_index],
            )
        )
    return searches


def _horizon_generate_coarse_dem(
    dem_fine: NDArray[Any],
    pyramid_scale: int,
    conv_from: int,
    conv_to: int,
    max_radius: int,
) -> NDArray[Any]:
    """Legacy max-pooling downsampling of one pyramid level, ported verbatim."""
    # First reduce the size for the edge required for horizon search.
    dem_fine = dem_fine[max_radius:-max_radius, max_radius:-max_radius]

    number_of_columns_fine = dem_fine.shape[1]
    number_of_lines_fine = dem_fine.shape[0]
    number_of_lines_coarse = int(np.floor(number_of_lines_fine / pyramid_scale)) + 1
    number_of_columns_coarse = int(np.floor(number_of_columns_fine / pyramid_scale)) + 1
    # The corner points must fit in the new grid, so pad columns to the right and
    # lines below until n_coarse = pyramid_scale * n_fine + 1 holds.
    modulo_columns = number_of_columns_fine % pyramid_scale
    pad_columns = 0
    if modulo_columns != 1:
        pad_columns = np.abs(1 - modulo_columns)
    modulo_lines = number_of_lines_fine % pyramid_scale
    pad_lines = 0
    if modulo_lines != 1:
        pad_lines = np.abs(1 - modulo_lines)
    # Extending only right and below, so edge padding is artefact-free here.
    dem_fine = np.pad(dem_fine, ((0, pad_lines), (0, pad_columns)), mode="edge")

    # Pad the data to support the rolling window of the convolution.
    dem_fine = np.pad(
        dem_fine, ((-conv_from, conv_to), (-conv_from, conv_to)), mode="symmetric"
    )

    # Convolution (keep maximum).
    dem_convolve = np.zeros(dem_fine.shape)
    for i in np.arange(pyramid_scale) + conv_from:
        for j in np.arange(pyramid_scale) + conv_from:
            dem_convolve = np.maximum(
                dem_convolve, np.roll(dem_fine, (i, j), axis=(0, 1))
            )
    # Divide by pyramid_scale to account for the change of resolution
    # (important for the angle computation later on).
    dem_convolve = dem_convolve / pyramid_scale

    # Keep only the convoluted points that survive the scale change; the slice end
    # is one past the last index so the last point is really included.
    dem_coarse = dem_convolve[
        -conv_from : (number_of_lines_coarse * pyramid_scale + 1) : pyramid_scale,
        -conv_from : (number_of_columns_coarse * pyramid_scale + 1) : pyramid_scale,
    ]

    # Final padding to enable searching the horizon over the edge: constant at the
    # minimal height so it never becomes the horizon (data-dependent legacy quirk;
    # an all-NaN-adjacent DEM minimum poisons the whole level).
    dem_coarse = np.pad(
        dem_coarse,
        ((max_radius, max_radius), (max_radius, max_radius)),
        mode="constant",
        constant_values=dem_coarse.min(),
    )

    return dem_coarse


def _horizon_generate_pyramids(
    dem: NDArray[np.float32],
    num_directions: int,
    max_fine_radius: int,
    max_pyramid_radius: int,
    pyramid_scale: int,
) -> list[_PyramidLevel]:
    """Legacy DEM pyramid construction for horizon search, ported verbatim."""
    # In levels above 0 the search starts at the first radius whose fine-level
    # distance is only representable at the coarser resolution.
    min_pyramid_radius = int(np.floor(max_pyramid_radius / pyramid_scale)) + 1

    conv_to = int(np.floor(pyramid_scale / 2.0))
    if (pyramid_scale % 2) == 0:
        conv_from = 1 - conv_to
    else:
        conv_from = -conv_to

    # Determine the number of levels and the last radius used in the highest level.
    pyramid_levels = 0
    last_radius = float(max_fine_radius)
    work = True
    while work:
        remaining_radius = max_fine_radius / pyramid_scale**pyramid_levels
        if remaining_radius > max_pyramid_radius:
            pyramid_levels = pyramid_levels + 1
        else:
            work = False
            last_radius = np.round(
                max_fine_radius / pyramid_scale**pyramid_levels, decimals=0
            )

    pyramid = []
    dem_coarse = dem
    for level in range(pyramid_levels + 1):
        # Level 0 uses a different minimal radius than the rest of the levels.
        if level == 0:
            min_radius = 1
            dem_fine = np.copy(
                np.pad(
                    dem, max_pyramid_radius, mode="constant", constant_values=dem.min()
                )
            )
        else:
            min_radius = min_pyramid_radius - 1
            dem_fine = np.copy(dem_coarse)
        # The last level uses a different maximal radius than the rest of the levels.
        if level == pyramid_levels:
            max_radius = last_radius
        else:
            max_radius = max_pyramid_radius
        searches = _horizon_shift_vector(num_directions, max_radius, min_radius)
        dem_coarse = _horizon_generate_coarse_dem(
            dem_fine, pyramid_scale, conv_from, conv_to, max_pyramid_radius
        )
        pyramid.append(
            _PyramidLevel(
                dem=dem_fine,
                searches=searches,
                line_indices=np.arange(dem_fine.shape[0]),
                column_indices=np.arange(dem_fine.shape[1]),
            )
        )

    return pyramid


def _compute_sky_illumination(
    dem: NDArray[np.float32],
    resolution: Resolution,
    compute_overcast: bool,
    compute_shadow: bool,
    shadow_and_horizon_only: bool,
    maximum_fine_radius: int,
    directions: int,
    shadow_azimuth: float,
    shadow_elevation: float,
    vertical_exaggeration: float,
) -> _SkyIlluminationOutputs:
    """
    Legacy ``sky_illumination`` math, ported verbatim for golden-master parity.

    Reproduced quirks: the overcast result is divided by its own maximum
    (data-dependent), the uniform+shadow blend skips the /pi normalisation the
    plain uniform result gets, direction angles are the azimuths rounded to one
    decimal, the horizon slopes ignore the map resolution (pixel distances), and a
    DEM containing NaN poisons the pyramid padding (``dem.min()`` is NaN) so the
    illumination and horizon come out all-NaN while the binary shadow is all zeros.
    """
    pyramid_scale = _PYRAMID_SCALE
    max_pyramid_radius = _MAXIMUM_PYRAMID_RADIUS

    dem = dem.astype(np.float32)
    dem = dem * vertical_exaggeration

    slope, aspect = _slope_aspect(
        np.pad(dem, max_pyramid_radius, mode="symmetric"),
        resolution,
        SlopeUnits.RADIAN,
        vertical_exaggeration=1.0,
    )

    pyramid = _horizon_generate_pyramids(
        dem,
        num_directions=directions,
        max_fine_radius=maximum_fine_radius,
        max_pyramid_radius=max_pyramid_radius,
        pyramid_scale=pyramid_scale,
    )
    number_of_levels = len(pyramid) - 1

    conv_to = int(np.floor(pyramid_scale / 2.0))
    if (pyramid_scale % 2) == 0:
        conv_from = 1 - conv_to
    else:
        conv_from = -conv_to
    # Directional halve-resolution for the integration limits.
    direction_half_angle = np.pi / directions

    padded_shape = (
        dem.shape[0] + 2 * max_pyramid_radius,
        dem.shape[1] + 2 * max_pyramid_radius,
    )
    uniform_a: NDArray[Any] = np.zeros(padded_shape, dtype=np.float32)
    uniform_b: NDArray[Any] = np.copy(uniform_a)
    overcast_c: NDArray[Any] = np.zeros(padded_shape, dtype=np.float32)
    overcast_d: NDArray[Any] = np.copy(overcast_c)

    # Use the closest pyramid direction as a proxy for the shadow azimuth.
    azimuths = np.array([search.azimuth for search in pyramid[0].searches])
    shadow_direction_index = int(np.argmin(np.abs(azimuths - (360 - shadow_azimuth))))
    shadow_out: NDArray[Any] | None = None
    horizon_out: NDArray[Any] | None = None

    for direction_index, direction in enumerate(pyramid[0].searches):
        direction_radians = np.radians(direction.azimuth)
        # Reset the maximum at each direction.
        max_slope: NDArray[Any] = (
            np.zeros(pyramid[number_of_levels].dem.shape, dtype=np.float32) - 1000
        )

        for level in reversed(range(number_of_levels + 1)):
            height = pyramid[level].dem
            search = pyramid[level].searches[direction_index]

            for radius_index, radius in enumerate(search.distances):
                shift_index = search.shifts[radius_index]
                elevation_slope = np.maximum(
                    (np.roll(height, shift_index, axis=(0, 1)) - height) / radius, 0.0
                )
                max_slope = np.maximum(max_slope, elevation_slope)

            # Resample the maximal slope down to the next finer pyramid level.
            if level > 0:
                lines_fine = pyramid[level - 1].line_indices + (
                    conv_from + max_pyramid_radius * pyramid_scale - max_pyramid_radius
                )
                columns_fine = pyramid[level - 1].column_indices + (
                    conv_from + max_pyramid_radius * pyramid_scale - max_pyramid_radius
                )
                lines_coarse = pyramid[level].line_indices * pyramid_scale
                columns_coarse = pyramid[level].column_indices * pyramid_scale
                interpolation_spline = RectBivariateSpline(
                    lines_coarse, columns_coarse, max_slope, kx=1, ky=1
                )
                max_slope = interpolation_spline(lines_fine, columns_fine)

        # Convert to an angle in radians and accumulate the directional output.
        horizon_angle = np.arctan(max_slope)
        uniform_a = uniform_a + (np.cos(horizon_angle)) ** 2
        aspect_term = (
            -2 * np.sin(direction_half_angle) * np.cos(direction_radians - aspect)
        )
        uniform_b = uniform_b + np.maximum(
            aspect_term
            * (np.pi / 4.0 - horizon_angle / 2.0 - np.sin(2.0 * horizon_angle) / 4.0),
            0,
        )
        if compute_overcast:
            cosine_cubed = (np.cos(horizon_angle)) ** 3
            overcast_c = overcast_c + np.maximum(cosine_cubed, 0)
            overcast_d = overcast_d + np.maximum(
                aspect_term * (2.0 / 3.0 - np.cos(horizon_angle) + cosine_cubed / 3.0),
                0,
            )
        if compute_shadow and direction_index == shadow_direction_index:
            horizon_out = np.degrees(
                horizon_angle[
                    max_pyramid_radius:-max_pyramid_radius,
                    max_pyramid_radius:-max_pyramid_radius,
                ]
            )
            shadow_out = (horizon_out < shadow_elevation) * 1
            if shadow_and_horizon_only:
                return _SkyIlluminationOutputs(
                    illumination=None,
                    shadow=shadow_out.astype(np.float32),
                    horizon=horizon_out.astype(np.float32),
                )

    # For numeric stability cap uniform_b at pi.
    uniform_out = direction_half_angle * np.cos(slope) * uniform_a + np.sin(
        slope
    ) * np.minimum(uniform_b, np.pi)
    uniform_out = uniform_out[
        max_pyramid_radius:-max_pyramid_radius, max_pyramid_radius:-max_pyramid_radius
    ]

    if compute_overcast:
        overcast_out = (2.0 * direction_half_angle / 3.0) * np.cos(
            slope
        ) * overcast_c + np.sin(slope) * overcast_d
        overcast_out = overcast_out[
            max_pyramid_radius:-max_pyramid_radius,
            max_pyramid_radius:-max_pyramid_radius,
        ]
        overcast_out = 0.33 * uniform_out + 0.67 * overcast_out
        overcast_out = overcast_out / overcast_out.max()
        if compute_shadow:
            assert shadow_out is not None
            illumination = 0.8 * overcast_out + 0.2 * shadow_out
        else:
            illumination = overcast_out
    else:
        if compute_shadow:
            assert shadow_out is not None
            # Legacy quirk: this blend uses the un-normalised uniform result.
            illumination = 0.8 * uniform_out + 0.2 * shadow_out
        else:
            illumination = uniform_out / np.pi

    return _SkyIlluminationOutputs(
        illumination=illumination.astype(np.float32),
        shadow=None if shadow_out is None else shadow_out.astype(np.float32),
        horizon=None if horizon_out is None else horizon_out.astype(np.float32),
    )


@dataclass(frozen=True)
class SkyIllumination:
    """
    Topographic correction for illumination under a uniform or overcast sky.

    :param sky_model: Sky luminance distribution to integrate over.
    :param compute_shadow: Blend a binary cast-shadow component into the result.
    :param maximum_fine_radius: Maximal shadow-modelling distance in pixels.
    :param directions: Number of azimuths searched for the horizon.
    :param shadow_azimuth: Cast-shadow sun azimuth in degrees, clockwise from
        north, 0-360.
    :param shadow_elevation: Cast-shadow sun elevation above the horizon in
        degrees, 0-90.
    :param vertical_exaggeration: Factor the DEM heights are multiplied by.
    """

    EIGHT_BIT_STRETCH: ClassVar[Stretch] = Stretch(StretchMode.PERCENT, 0.25, 0.0)

    sky_model: SkyModel = SkyModel.OVERCAST
    compute_shadow: bool = False
    maximum_fine_radius: int = 100
    directions: int = 32
    shadow_azimuth: float = 315.0
    shadow_elevation: float = 35.0
    vertical_exaggeration: float = 1.0

    def __post_init__(self) -> None:
        validate_vertical_exaggeration(self.vertical_exaggeration)
        if not (0 <= self.shadow_azimuth <= 360):
            raise ValueError(
                f"shadow_azimuth must be between 0 and 360, got '{self.shadow_azimuth}'!"
            )
        if not (0 <= self.shadow_elevation <= 90):
            raise ValueError(
                f"shadow_elevation must be between 0 and 90, got '{self.shadow_elevation}'!"
            )
        if self.directions < 1:
            raise ValueError(f"directions must be positive, got '{self.directions}'!")
        if self.maximum_fine_radius < 1:
            raise ValueError(
                f"maximum_fine_radius must be positive, got '{self.maximum_fine_radius}'!"
            )

    def output_stem(self, dem_name: str) -> str:
        """
        Return the parameter-encoded output file stem, e.g.
        ``dem_SIM_overcast_D32_100px``.

        :param dem_name: Stem of the DEM file name.
        :returns: The output file stem, without extension.
        """
        return (
            f"{dem_name}_SIM_{self.sky_model}_D{self.directions}"
            f"_{self.maximum_fine_radius}px"
        )

    def compute(
        self, dem: NDArray[np.float32], resolution: Resolution
    ) -> NDArray[np.float32]:
        """
        Compute the sky illumination of a DEM.

        The overcast result is normalised by its own maximum (a data-dependent
        legacy behaviour, kept for parity); a DEM containing any NaN yields an
        all-NaN result because the NaN poisons the horizon pyramid padding.

        :param dem: 2D float32 DEM with NaN as nodata. Not modified.
        :param resolution: Cell size of the DEM.
        :returns: Sky illumination array of the DEM's shape.
        :raises ValueError: If the DEM is not 2D or the resolution is not positive.
        """
        validate_dem_and_resolution(dem.ndim, resolution)
        outputs = _compute_sky_illumination(
            dem,
            resolution,
            compute_overcast=self.sky_model == SkyModel.OVERCAST,
            compute_shadow=self.compute_shadow,
            shadow_and_horizon_only=False,
            maximum_fine_radius=self.maximum_fine_radius,
            directions=self.directions,
            shadow_azimuth=self.shadow_azimuth,
            shadow_elevation=self.shadow_elevation,
            vertical_exaggeration=self.vertical_exaggeration,
        )
        assert outputs.illumination is not None
        return outputs.illumination
