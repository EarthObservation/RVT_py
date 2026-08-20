"""
Blend modes: legacy blend formulas behind a typed enum.

Legacy quirk preserved for parity: ``blend_overlay`` and ``blend_soft_light``
overwrote the ``background`` array in place, so the opacity mix that follows
in the render loop blended the result with itself (the layer opacity is
effectively ignored for those modes). The engine replicates that dataflow on
arrays it owns; the public :meth:`BlendMode.blend` is pure.
"""

from __future__ import annotations

from enum import StrEnum, auto

import numpy as np
from numpy.typing import NDArray


class BlendMode(StrEnum):
    """
    How an active layer combines with the background under it.
    """

    NORMAL = auto()
    MULTIPLY = auto()
    OVERLAY = auto()
    LUMINOSITY = auto()
    SCREEN = auto()
    SOFT_LIGHT = auto()

    def blend(
        self, active: NDArray[np.floating], background: NDArray[np.floating]
    ) -> NDArray[np.floating]:
        """
        Blend the active image over the background image.

        Both images are expected on the 0-1 range, each single band or
        bands-first RGB. Pure: neither input array is modified.

        :param active: The layer being blended in.
        :param background: The image rendered so far.
        :returns: The blended array.
        """
        if self._mutates_background():
            background = background.copy()
        return blend_images_inplace(self, active, background)

    def _mutates_background(self) -> bool:
        """Return True for the modes whose legacy formula overwrites the background."""
        return self in (BlendMode.OVERLAY, BlendMode.SOFT_LIGHT)


def blend_images_inplace(
    blend_mode: BlendMode,
    active: NDArray[np.floating],
    background: NDArray[np.floating],
) -> NDArray[np.floating]:
    """
    Legacy ``blend_images``, ported verbatim including its mutation
    semantics: overlay and soft-light overwrite ``background`` (or its bands)
    with blended values, which the legacy render loop relied on. The caller
    must own ``background``.

    :param blend_mode: The blend mode.
    :param active: The layer being blended in. Not modified.
    :param background: The image rendered so far; overwritten by overlay and
        soft-light.
    :returns: The blended array.
    """
    match blend_mode:
        case BlendMode.NORMAL:
            return active
        case BlendMode.LUMINOSITY:
            return _blend_luminosity(active, background)
        case _:
            return _blend_per_band(blend_mode, active, background)


def _equation_blend(
    blend_mode: BlendMode,
    active: NDArray[np.floating],
    background: NDArray[np.floating],
) -> NDArray[np.floating]:
    """Legacy ``equation_blend`` over single bands; overlay and soft-light mutate ``background``."""
    match blend_mode:
        case BlendMode.SCREEN:
            return 1 - (1 - active) * (1 - background)
        case BlendMode.MULTIPLY:
            return active * background
        case BlendMode.OVERLAY:
            return _blend_overlay_inplace(active, background)
        case BlendMode.SOFT_LIGHT:
            return _blend_soft_light_inplace(active, background)
        case _:
            raise ValueError(f"Blend mode '{blend_mode}' is not separable per band!")


def _blend_overlay_inplace(
    active: NDArray[np.floating], background: NDArray[np.floating]
) -> NDArray[np.floating]:
    above = np.where(background > 0.5)
    below = np.where(background <= 0.5)
    background[above] = 1 - (1 - 2 * (background[above] - 0.5)) * (1 - active[above])
    background[below] = (2 * background[below]) * active[below]
    return background


def _blend_soft_light_inplace(
    active: NDArray[np.floating], background: NDArray[np.floating]
) -> NDArray[np.floating]:
    below = np.where(active < 0.5)
    above = np.where(active >= 0.5)
    background[below] = 2 * background[below] * active[below] + background[
        below
    ] ** 2 * (1.0 - 2 * active[below])
    background[above] = 2 * background[above] * (1.0 - active[above]) + np.sqrt(
        background[above]
    ) * (2 * active[above] - 1.0)
    return background


def _blend_per_band(
    blend_mode: BlendMode,
    active: NDArray[np.floating],
    background: NDArray[np.floating],
) -> NDArray[np.floating]:
    """Legacy ``blend_multi_dim_images``: broadcast the formula over RGB bands."""
    active_is_rgb = active.ndim == 3
    background_is_rgb = background.ndim == 3

    if not active_is_rgb and not background_is_rgb:
        return _equation_blend(blend_mode, active, background)

    blended = np.zeros(background.shape if background_is_rgb else active.shape)
    for band in range(3):
        active_band = active[band] if active_is_rgb else active
        # For a 2D background under an RGB active image, legacy mutated the one
        # background across the band loop, cascading each band's overlay or
        # soft-light result into the next band's input.
        background_band = background[band] if background_is_rgb else background
        blended[band] = _equation_blend(blend_mode, active_band, background_band)
    return blended


def luminance(image: NDArray[np.floating]) -> NDArray[np.floating]:
    """
    Legacy ``lum``: the 0.3R + 0.59G + 0.11B luminance of an RGB image;
    single-band images are their own luminance.

    :param image: Single band or bands-first RGB image.
    :returns: The single-band luminance array.
    """
    if image.ndim == 3:
        # Legacy wrapped this in np.float32(...); astype is the same cast.
        combined = (0.3 * image[0]) + (0.59 * image[1]) + (0.11 * image[2])
        return combined.astype(np.float32)
    return image


def _blend_luminosity(
    active: NDArray[np.floating], background: NDArray[np.floating]
) -> NDArray[np.floating]:
    lum_active = luminance(active)
    lum_background = luminance(background)
    luminosity = lum_active - lum_background

    if background.ndim < 3:
        return lum_active

    shifted = np.zeros(background.shape)
    for band in range(3):
        shifted[band] = background[band] + luminosity
    return _clip_color(shifted)


def _channel_extremum(
    red: NDArray[np.floating],
    green: NDArray[np.floating],
    blue: NDArray[np.floating],
    take_maximum: bool,
) -> NDArray[np.floating]:
    extremum = red * 1.0
    for channel in (green, blue):
        if take_maximum:
            index = np.where(channel > extremum)
        else:
            index = np.where(channel < extremum)
        extremum[index] = channel[index]
    return extremum


def _clip_color(color: NDArray[np.floating]) -> NDArray[np.floating]:
    """Legacy ``clip_color`` gamut clipping after a luminance shift."""
    lum_color = luminance(color)

    red = color[0].astype(np.float32)
    green = color[1].astype(np.float32)
    blue = color[2].astype(np.float32)

    minimum = _channel_extremum(red, green, blue, take_maximum=False)
    maximum = _channel_extremum(red, green, blue, take_maximum=True)

    below = np.where(minimum < 0)
    above = np.where(maximum > 1)
    for channel in (red, green, blue):
        channel[below] = lum_color[below] + (
            (channel[below] - lum_color[below]) * lum_color[below]
        ) / (lum_color[below] - minimum[below])
        channel[above] = lum_color[above] + (
            (channel[above] - lum_color[above]) * (1.0 - lum_color[above])
        ) / (maximum[above] - lum_color[above])

    clipped = np.zeros(color.shape)
    clipped[0] = red
    clipped[1] = green
    clipped[2] = blue
    return clipped
