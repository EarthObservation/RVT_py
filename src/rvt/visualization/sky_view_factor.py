"""
Sky-view factor, anisotropic sky-view factor, and openness.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum, auto
from typing import ClassVar, NamedTuple

import numpy as np
from numpy.typing import NDArray

from rvt.core import (
    Resolution,
    format_parameter,
    validate_dem_and_resolution,
    validate_vertical_exaggeration,
)
from rvt.stretch import Stretch, StretchMode


class NoiseRemoval(StrEnum):
    """
    Noise removal level for horizon searches; higher levels skip a larger
    share of the search radius closest to the pixel.
    """

    NONE = auto()
    LOW = auto()
    MEDIUM = auto()
    HIGH = auto()

    def _file_name_suffix(self) -> str:
        """Return the file name suffix encoding this level, empty for ``NONE``."""
        if self == NoiseRemoval.NONE:
            return ""
        return f"_NR{self.value}"

    def _minimum_radius_percent(self) -> float:
        """Return the percent of the maximal search radius ignored near the pixel."""
        match self:
            case NoiseRemoval.NONE:
                return 0.0
            case NoiseRemoval.LOW:
                return 10.0
            case NoiseRemoval.MEDIUM:
                return 20.0
            case NoiseRemoval.HIGH:
                return 40.0


class SkyViewFactorResult(NamedTuple):
    """
    Result of a sky-view factor computation.

    :param sky_view_factor: Sky-view factor, the visible portion of the sky
        hemisphere, 0-1.
    :param anisotropic: Anisotropic sky-view factor, or ``None`` when it was
        not requested.
    :param positive_openness: Positive openness (mean horizon elevation angle
        in degrees), or ``None`` when it was not requested.
    """

    sky_view_factor: NDArray[np.float32]
    anisotropic: NDArray[np.float32] | None
    positive_openness: NDArray[np.float32] | None


class _DirectionSearch(NamedTuple):
    """Horizon search offsets for one azimuth: np.roll shifts and their distances."""

    shifts: list[tuple[int, int]]
    distances: NDArray[np.float64]


def _horizon_shift_vector(
    num_directions: int, radius_pixels: int, min_radius: float
) -> list[_DirectionSearch]:
    """Legacy sky-view search-movement construction, ported verbatim for parity."""
    angles = (2 * np.pi / num_directions) * np.arange(num_directions)
    x = np.cos(angles)
    y = np.sin(angles)

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
                shifts=[(k[0], k[1]) for k in shift_pairs[sort_index]],
                distances=distance[sort_index],
            )
        )
    return searches


def _sky_view_factor_compute(
    height_array: NDArray[np.float32],
    radius_max: int,
    radius_min: float,
    num_directions: int,
    compute_svf: bool,
    compute_asvf: bool,
    compute_opns: bool,
    a_main_direction: float,
    a_poly_level: int,
    a_min_weight: float,
) -> tuple[
    NDArray[np.floating] | None,
    NDArray[np.floating] | None,
    NDArray[np.floating] | None,
]:
    """Legacy horizon-based SVF/ASVF/openness math, ported verbatim for parity."""
    height = np.pad(height_array, radius_max, mode="reflect")

    move = _horizon_shift_vector(
        num_directions=num_directions,
        radius_pixels=radius_max,
        min_radius=radius_min,
    )

    # Multiply with 0 instead of using np.zeros to preserve nodata.
    svf_out = height * 0 if compute_svf else None
    asvf_out = height * 0 if compute_asvf else None
    opns_out = height * 0 if compute_opns else None

    w_m = a_min_weight
    w_a = np.deg2rad(a_main_direction)
    weight = np.arange(num_directions) * (2 * np.pi / num_directions)
    weight = (1 - w_m) * (np.cos((weight - w_a) / 2)) ** a_poly_level + w_m

    for i_dir, direction in enumerate(move):
        # Smallest possible elevation angle is -1000 rad (i.e. -90 deg).
        max_slope = np.zeros(height.shape, dtype=np.float32) - 1000

        for i_rad, radius in enumerate(direction.distances):
            shift_index = direction.shifts[i_rad]
            slope_estimate = (np.roll(height, shift_index, axis=(0, 1)) - height) / (
                radius
            )
            # np.fmax prevents NaN values contaminating the edge of the image.
            max_slope = np.fmax(max_slope, slope_estimate)

        max_slope = np.arctan(max_slope)

        if svf_out is not None:
            # Minimum possible angle is 0 (hemisphere); np.fmax changes NaNs to 0.
            svf_out = svf_out + (1 - np.sin(np.fmax(max_slope, 0)))
        if asvf_out is not None:
            asvf_out = asvf_out + (1 - np.sin(np.fmax(max_slope, 0))) * weight[i_dir]
        if opns_out is not None:
            # For openness taking the entire sphere.
            opns_out = opns_out + max_slope

    # Cut to original extent and average the directional output over all directions.
    if svf_out is not None:
        svf_out = svf_out[radius_max:-radius_max, radius_max:-radius_max] / (
            num_directions
        )
    if asvf_out is not None:
        asvf_out = asvf_out[radius_max:-radius_max, radius_max:-radius_max] / np.sum(
            weight
        )
    if opns_out is not None:
        opns_out = np.rad2deg(
            0.5 * np.pi
            - (opns_out[radius_max:-radius_max, radius_max:-radius_max])
            / num_directions
        )

    return svf_out, asvf_out, opns_out


def _anisotropy_parameters(anisotropy_level: int) -> tuple[int, float]:
    """Return (polynomial level, minimum weight) for the anisotropy level."""
    match anisotropy_level:
        case 1:
            return 4, 0.4
        case 2:
            return 8, 0.1
        case _:
            raise ValueError(
                f"anisotropy_level must be 1 (low) or 2 (high),"
                f" got '{anisotropy_level}'!"
            )


def _compute_sky_view_family(
    dem: NDArray[np.float32],
    resolution: Resolution,
    directions: int,
    max_radius: int,
    noise_removal: NoiseRemoval,
    compute_svf: bool,
    compute_asvf: bool,
    compute_opns: bool,
    anisotropy_direction: float,
    anisotropy_level: int,
    vertical_exaggeration: float,
) -> tuple[
    NDArray[np.float32] | None,
    NDArray[np.float32] | None,
    NDArray[np.float32] | None,
]:
    """Legacy preparation around :func:`_sky_view_factor_compute`, ported for parity."""
    validate_dem_and_resolution(dem.ndim, resolution)
    if resolution.x != resolution.y:
        raise ValueError(
            f"sky-view computations require square pixels, got '{resolution}'!"
        )

    dem = dem.astype(np.float32)
    # Save the NaN mask (processing may change NaNs to arbitrary values).
    nan_mask = np.isnan(dem)

    dem = dem * vertical_exaggeration
    # Adjust elevation to the pixel size, the calculation thinks 1 px == 1 m.
    dem = dem / resolution.x

    # Minimal search radius has to be an integer not smaller than 1.
    radius_min = max(
        np.round(
            max_radius * noise_removal._minimum_radius_percent() * 0.01, decimals=0
        ),
        1,
    )
    poly_level, min_weight = _anisotropy_parameters(anisotropy_level)

    svf_out, asvf_out, opns_out = _sky_view_factor_compute(
        height_array=dem,
        radius_max=max_radius,
        radius_min=radius_min,
        num_directions=directions,
        compute_svf=compute_svf,
        compute_asvf=compute_asvf,
        compute_opns=compute_opns,
        a_main_direction=anisotropy_direction,
        a_poly_level=poly_level,
        a_min_weight=min_weight,
    )

    outputs = []
    for output in (svf_out, asvf_out, opns_out):
        if output is None:
            outputs.append(None)
            continue
        output = output.astype(np.float32)
        output[nan_mask] = np.nan
        outputs.append(output)
    return outputs[0], outputs[1], outputs[2]


def _validate_search_parameters(directions: int, max_radius: int) -> None:
    """Raise ``ValueError`` for a non-positive direction count or search radius."""
    if directions < 1:
        raise ValueError(f"directions must be a positive integer, got '{directions}'!")
    if max_radius < 1:
        raise ValueError(f"max_radius must be a positive integer, got '{max_radius}'!")


@dataclass(frozen=True)
class SkyViewFactor:
    """
    Sky-view factor with optional anisotropic variant and positive openness.

    All three visualizations share one horizon search, so requesting the
    optional outputs alongside the sky-view factor is much cheaper than
    computing them separately.

    :param directions: Number of horizon search directions.
    :param max_radius: Maximal search radius in pixels.
    :param noise_removal: Noise removal level; higher levels ignore a larger
        share of the search radius closest to the pixel.
    :param include_anisotropic: Also compute the anisotropic sky-view factor.
    :param include_positive_openness: Also compute positive openness.
    :param anisotropy_direction: Main direction of anisotropy in degrees.
    :param anisotropy_level: Level of anisotropy, 1 (low) or 2 (high).
    :param vertical_exaggeration: Factor the DEM heights are multiplied by.
    """

    EIGHT_BIT_STRETCH: ClassVar[Stretch] = Stretch(StretchMode.VALUE, 0.6375, 1.0)
    ANISOTROPIC_EIGHT_BIT_STRETCH: ClassVar[Stretch] = Stretch(
        StretchMode.VALUE, 0.7, 0.9
    )
    POSITIVE_OPENNESS_EIGHT_BIT_STRETCH: ClassVar[Stretch] = Stretch(
        StretchMode.VALUE, 60.0, 95.0
    )

    directions: int = 16
    max_radius: int = 10
    noise_removal: NoiseRemoval = NoiseRemoval.NONE
    include_anisotropic: bool = False
    include_positive_openness: bool = False
    anisotropy_direction: float = 315.0
    anisotropy_level: int = 1
    vertical_exaggeration: float = 1.0

    def __post_init__(self) -> None:
        validate_vertical_exaggeration(self.vertical_exaggeration)
        _validate_search_parameters(self.directions, self.max_radius)
        _anisotropy_parameters(self.anisotropy_level)

    def sky_view_factor_output_stem(self, dem_name: str) -> str:
        """
        Return the output file stem for the sky-view factor member, e.g.
        ``dem_SVF_R10_D16`` (a noise-removal suffix such as ``_NRhigh`` is
        appended when noise removal is enabled).

        :param dem_name: Stem of the DEM file name.
        :returns: The output file stem, without extension.
        """
        return (
            f"{dem_name}_SVF_R{self.max_radius}_D{self.directions}"
            f"{self.noise_removal._file_name_suffix()}"
        )

    def anisotropic_output_stem(self, dem_name: str) -> str:
        """
        Return the output file stem for the anisotropic member, e.g.
        ``dem_SVF-A_R10_D16_A315_ALlow``.

        :param dem_name: Stem of the DEM file name.
        :returns: The output file stem, without extension.
        """
        anisotropy_level = "low" if self.anisotropy_level == 1 else "high"
        return (
            f"{dem_name}_SVF-A_R{self.max_radius}_D{self.directions}"
            f"_A{format_parameter(self.anisotropy_direction)}_AL{anisotropy_level}"
            f"{self.noise_removal._file_name_suffix()}"
        )

    def positive_openness_output_stem(self, dem_name: str) -> str:
        """
        Return the output file stem for the positive-openness member, e.g.
        ``dem_OPEN-POS_R10_D16``.

        :param dem_name: Stem of the DEM file name.
        :returns: The output file stem, without extension.
        """
        return (
            f"{dem_name}_OPEN-POS_R{self.max_radius}_D{self.directions}"
            f"{self.noise_removal._file_name_suffix()}"
        )

    def compute(
        self, dem: NDArray[np.float32], resolution: Resolution
    ) -> SkyViewFactorResult:
        """
        Compute the sky-view factor (and requested companions) of a DEM.

        :param dem: 2D float32 DEM with NaN as nodata. Not modified.
        :param resolution: Cell size of the DEM; pixels must be square.
        :returns: Result with the optional members ``None`` when not requested.
        :raises ValueError: If the DEM is not 2D, the resolution is not
            positive, or the pixels are not square.
        """
        svf_array, anisotropic_array, positive_openness_array = (
            _compute_sky_view_family(
                dem,
                resolution,
                directions=self.directions,
                max_radius=self.max_radius,
                noise_removal=self.noise_removal,
                compute_svf=True,
                compute_asvf=self.include_anisotropic,
                compute_opns=self.include_positive_openness,
                anisotropy_direction=self.anisotropy_direction,
                anisotropy_level=self.anisotropy_level,
                vertical_exaggeration=self.vertical_exaggeration,
            )
        )
        assert svf_array is not None  # compute_svf is always requested
        return SkyViewFactorResult(
            sky_view_factor=svf_array,
            anisotropic=anisotropic_array,
            positive_openness=positive_openness_array,
        )


@dataclass(frozen=True)
class NegativeOpenness:
    """
    Negative openness: positive openness of the negated DEM (mean horizon
    depression angle in degrees).

    :param directions: Number of horizon search directions.
    :param max_radius: Maximal search radius in pixels.
    :param noise_removal: Noise removal level; higher levels ignore a larger
        share of the search radius closest to the pixel.
    :param vertical_exaggeration: Factor the DEM heights are multiplied by.
    """

    # Legacy 8-bit negative openness is inverted (deep is dark).
    EIGHT_BIT_STRETCH: ClassVar[Stretch] = Stretch(
        StretchMode.VALUE, 60.0, 95.0, invert=True
    )

    directions: int = 16
    max_radius: int = 10
    noise_removal: NoiseRemoval = NoiseRemoval.NONE
    vertical_exaggeration: float = 1.0

    def __post_init__(self) -> None:
        validate_vertical_exaggeration(self.vertical_exaggeration)
        _validate_search_parameters(self.directions, self.max_radius)

    def output_stem(self, dem_name: str) -> str:
        """
        Return the parameter-encoded output file stem, e.g. ``dem_OPEN-NEG_R10_D16``
        (the noise-removal suffix for level 3 is ``_NRhigh``, consistently with the
        other sky-view visualizations, where legacy used ``_NRstrong`` here).

        :param dem_name: Stem of the DEM file name.
        :returns: The output file stem, without extension.
        """
        return (
            f"{dem_name}_OPEN-NEG_R{self.max_radius}_D{self.directions}"
            f"{self.noise_removal._file_name_suffix()}"
        )

    def compute(
        self, dem: NDArray[np.float32], resolution: Resolution
    ) -> NDArray[np.float32]:
        """
        Compute the negative openness of a DEM.

        :param dem: 2D float32 DEM with NaN as nodata. Not modified.
        :param resolution: Cell size of the DEM; pixels must be square.
        :returns: Negative openness array of the DEM's shape, in degrees.
        :raises ValueError: If the DEM is not 2D, the resolution is not
            positive, or the pixels are not square.
        """
        # The legacy pipeline negates the DEM before vertical exaggeration.
        dem = -1 * dem
        _, _, positive_openness_array = _compute_sky_view_family(
            dem,
            resolution,
            directions=self.directions,
            max_radius=self.max_radius,
            noise_removal=self.noise_removal,
            compute_svf=False,
            compute_asvf=False,
            compute_opns=True,
            anisotropy_direction=315.0,
            anisotropy_level=1,
            vertical_exaggeration=self.vertical_exaggeration,
        )
        assert positive_openness_array is not None  # compute_opns is requested
        return positive_openness_array
