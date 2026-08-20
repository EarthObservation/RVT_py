"""
Blend layers and combinations: the typed shell around the legacy blending
engine (bottom-up rendering, last layer is the base).
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Protocol

import numpy as np
from numpy.typing import NDArray

from rvt.blending.blend_modes import BlendMode, blend_images_inplace
from rvt.blending.normalization import Normalization, normalize_image, scale_0_to_1
from rvt.core import Resolution
from rvt.visualization import NegativeOpenness, SkyViewFactor, SlopeAspect, SlopeUnits


class _ImageVisualization(Protocol):
    """A visualization config whose ``compute`` returns a bare image array."""

    def compute(
        self, dem: NDArray[np.float32], resolution: Resolution
    ) -> NDArray[np.float32]: ...


@dataclass(frozen=True)
class ImageSource:
    """
    A precomputed image used as a blend layer.

    :param image: The image array (single band, or bands-first RGB).
    :param inverts_normalization: Invert the normalized image (set for images
        that carry the slope-gradient or negative-openness identity, where
        legacy renders high values dark).
    """

    image: NDArray[np.floating]
    inverts_normalization: bool = False

    def compute_image(
        self, dem: NDArray[np.float32], resolution: Resolution
    ) -> NDArray[np.floating]:
        """Return the precomputed image; the DEM is ignored."""
        return self.image


@dataclass(frozen=True)
class VisualizationSource:
    """
    A blend layer computed from any visualization config that returns a bare
    image array (e.g. :class:`rvt.visualization.Hillshade`,
    :class:`rvt.visualization.SimpleLocalReliefModel`,
    :class:`rvt.visualization.LocalDominance`).

    :param config: The visualization config.
    :param inverts_normalization: Invert the normalized image.
    """

    config: _ImageVisualization
    inverts_normalization: bool = False

    def compute_image(
        self, dem: NDArray[np.float32], resolution: Resolution
    ) -> NDArray[np.floating]:
        """Compute the visualization image from the DEM."""
        return self.config.compute(dem, resolution)


@dataclass(frozen=True)
class SlopeGradientSource:
    """
    The slope band of :class:`rvt.visualization.SlopeAspect` as a blend
    layer. Defaults to degrees (the legacy blender's slope unit) and to the
    legacy inverted normalization (steep is dark).

    :param config: The slope/aspect config.
    :param inverts_normalization: Invert the normalized image.
    """

    config: SlopeAspect = field(default=SlopeAspect(units=SlopeUnits.DEGREE))
    inverts_normalization: bool = True

    def compute_image(
        self, dem: NDArray[np.float32], resolution: Resolution
    ) -> NDArray[np.floating]:
        """Compute the slope image from the DEM."""
        return self.config.compute(dem, resolution).slope


@dataclass(frozen=True)
class SkyViewFactorSource:
    """
    The sky-view factor band of :class:`rvt.visualization.SkyViewFactor` as a
    blend layer.

    :param config: The sky-view factor config.
    :param inverts_normalization: Invert the normalized image.
    """

    config: SkyViewFactor = field(default=SkyViewFactor())
    inverts_normalization: bool = False

    def compute_image(
        self, dem: NDArray[np.float32], resolution: Resolution
    ) -> NDArray[np.floating]:
        """Compute the sky-view factor image from the DEM."""
        return self.config.compute(dem, resolution).sky_view_factor


@dataclass(frozen=True)
class PositiveOpennessSource:
    """
    The positive-openness band of :class:`rvt.visualization.SkyViewFactor` as
    a blend layer.

    :param config: The sky-view factor config; must request positive openness.
    :param inverts_normalization: Invert the normalized image.
    """

    config: SkyViewFactor = field(default=SkyViewFactor(include_positive_openness=True))
    inverts_normalization: bool = False

    def __post_init__(self) -> None:
        if not self.config.include_positive_openness:
            raise ValueError(
                "PositiveOpennessSource requires a config with"
                " include_positive_openness=True!"
            )

    def compute_image(
        self, dem: NDArray[np.float32], resolution: Resolution
    ) -> NDArray[np.floating]:
        """Compute the positive-openness image from the DEM."""
        positive_openness = self.config.compute(dem, resolution).positive_openness
        assert positive_openness is not None  # guaranteed by __post_init__
        return positive_openness


@dataclass(frozen=True)
class NegativeOpennessSource:
    """
    :class:`rvt.visualization.NegativeOpenness` as a blend layer, with the
    legacy inverted normalization (deep is dark).

    :param config: The negative-openness config.
    :param inverts_normalization: Invert the normalized image.
    """

    config: NegativeOpenness = field(default=NegativeOpenness())
    inverts_normalization: bool = True

    def compute_image(
        self, dem: NDArray[np.float32], resolution: Resolution
    ) -> NDArray[np.floating]:
        """Compute the negative-openness image from the DEM."""
        return self.config.compute(dem, resolution)


type LayerSource = (
    ImageSource
    | VisualizationSource
    | SlopeGradientSource
    | SkyViewFactorSource
    | PositiveOpennessSource
    | NegativeOpennessSource
)


@dataclass(frozen=True)
class BlendLayer:
    """
    One layer of a blend combination.

    :param source: What image the layer contributes (a visualization computed
        from the DEM, or a precomputed image).
    :param normalization: How the image is normalized to 0-1 before blending.
    :param blend_mode: How the layer combines with the layers below it.
    :param opacity: Layer opacity in percent, 0-100.
    """

    source: LayerSource
    normalization: Normalization
    blend_mode: BlendMode = BlendMode.NORMAL
    opacity: float = 100.0

    def __post_init__(self) -> None:
        if not (0 <= self.opacity <= 100):
            raise ValueError(
                f"opacity must be between 0 and 100, got '{self.opacity}'!"
            )


@dataclass(frozen=True)
class BlendCombination:
    """
    An ordered stack of blend layers, rendered bottom-up: the last layer is
    the base, and each earlier layer is blended over the result.

    :param name: Human-readable name of the combination.
    :param layers: The layer stack, top first (legacy order).
    """

    name: str
    layers: tuple[BlendLayer, ...]

    def __post_init__(self) -> None:
        if not self.layers:
            raise ValueError("A blend combination needs at least one layer!")

    def render(
        self, dem: NDArray[np.float32], resolution: Resolution
    ) -> NDArray[np.float32]:
        """
        Render the combination over a DEM: per layer compute, normalize, then
        blend and composite with the layer's opacity (legacy engine, ported
        verbatim but pure).

        :param dem: 2D float32 DEM with NaN as nodata. Not modified.
        :param resolution: Cell size of the DEM.
        :returns: The rendered float image on the 0-1 range (single band, or
            bands-first RGB when a layer contributes RGB).
        """
        rendered: NDArray[np.floating] | None = None
        for layer in reversed(self.layers):
            image = layer.source.compute_image(dem, resolution)
            normalized = normalize_image(
                image, layer.normalization, invert=layer.source.inverts_normalization
            )

            if rendered is None:
                rendered = normalized
                continue

            active = normalized
            background = rendered
            if np.nanmin(active) < 0 or np.nanmax(active) > 1:
                active = scale_0_to_1(active)
            if np.nanmin(background) < 0 or np.nanmax(background) > 1:
                background = scale_0_to_1(background)

            # The engine owns `background` (rendered state), so the legacy
            # in-place blend is safe here; overlay and soft-light overwrite it
            # with the blended values, which makes the opacity mix below blend
            # the result with itself — the legacy quirk parity requires.
            top = blend_images_inplace(layer.blend_mode, active, background)
            rendered = _composite(top, background, layer.opacity)

            if np.nanmin(rendered) < 0 or np.nanmax(rendered) > 1:
                warnings.warn(
                    "rvt.blending.BlendCombination.render: rendered image"
                    " scale distorted"
                )
        assert rendered is not None  # at least one layer is guaranteed
        return rendered

    def render_8bit(
        self, dem: NDArray[np.float32], resolution: Resolution
    ) -> NDArray[np.uint8]:
        """
        Render the combination and byte-scale the 0-1 result to 0-255 with
        NaN encoded as 255 (legacy semantics).

        :param dem: 2D float32 DEM with NaN as nodata. Not modified.
        :param resolution: Cell size of the DEM.
        :returns: The uint8 array.
        """
        return _byte_scale_0_to_1(self.render(dem, resolution))


def _apply_opacity(
    active: NDArray[np.floating], background: NDArray[np.floating], opacity: float
) -> NDArray[np.floating]:
    # Legacy treated opacity values of at most 1 as fractions; BlendLayer
    # validates a real 0-100 percentage, so the division is unconditional.
    fraction = opacity / 100
    return active * fraction + background * (1 - fraction)


def _composite(
    active: NDArray[np.floating], background: NDArray[np.floating], opacity: float
) -> NDArray[np.floating]:
    """Legacy ``render_images``: rescale out-of-range inputs, mix per band by opacity."""
    if np.nanmin(active) < 0 or np.nanmax(active) > 1:
        active = scale_0_to_1(active)
    if np.nanmin(background) < 0 or np.nanmax(background) > 1:
        background = scale_0_to_1(background)

    active_is_rgb = active.ndim == 3
    background_is_rgb = background.ndim == 3

    if not active_is_rgb and not background_is_rgb:
        return _apply_opacity(active, background, opacity)

    composited = np.zeros(background.shape if background_is_rgb else active.shape)
    for band in range(3):
        active_band = active[band] if active_is_rgb else active
        background_band = background[band] if background_is_rgb else background
        composited[band] = _apply_opacity(active_band, background_band, opacity)
    return composited


def _byte_scale_0_to_1(data: NDArray[np.floating]) -> NDArray[np.uint8]:
    """Legacy IDL BYTSCL port over c_min=0, c_max=1: linear scale to 0-255, NaN encoded as 255."""
    high = 255
    byte_data = (high + 0.9999) * data
    byte_data[byte_data > high] = high
    byte_data[byte_data < 0] = 0
    byte_data[np.isnan(byte_data)] = 255
    return np.asarray(byte_data, dtype=np.uint8)
