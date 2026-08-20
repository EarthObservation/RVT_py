"""
Coloring of normalized grayscale images with matplotlib color ramps.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from matplotlib import colormaps
from matplotlib.colors import Colormap, LinearSegmentedColormap
from numpy.typing import NDArray

# The legacy truncation resampled the colormap at 100 points; parity requires
# the same quantization even for a full 0-1 cut.
_TRUNCATION_SAMPLES = 100


@dataclass(frozen=True)
class ColorRamp:
    """
    A matplotlib colormap, optionally cut to a sub-range, applied to a
    normalized (0-1) grayscale image to produce a bands-first RGB image.

    :param name: Matplotlib colormap name
        (https://matplotlib.org/stable/users/explain/colors/colormaps.html).
    :param minimum_cut: Lower part of the colormap to cut away, 0-1, or
        ``None`` for no lower cut. Note that setting any cut (even the full
        0-1 range) resamples the colormap at 100 points, exactly like the
        legacy blender.
    :param maximum_cut: Upper part of the colormap to cut away, 0-1, or
        ``None`` for no upper cut.
    """

    name: str
    minimum_cut: float | None = None
    maximum_cut: float | None = None

    def __post_init__(self) -> None:
        if self.name not in colormaps:
            raise ValueError(f"'{self.name}' is not a matplotlib colormap!")
        if self.minimum_cut is None and self.maximum_cut is None:
            return
        minimum, maximum = self._resolved_cuts()
        if not (0 <= minimum <= 1 and 0 <= maximum <= 1):
            raise ValueError(
                f"colormap cuts must be between 0 and 1, got '{minimum}'"
                f" and '{maximum}'!"
            )
        if minimum >= maximum:
            raise ValueError(
                f"minimum_cut must be smaller than maximum_cut, got '{minimum}'"
                f" and '{maximum}'!"
            )

    def apply(self, gray_scale: NDArray[np.floating]) -> NDArray[np.float64]:
        """
        Color a normalized grayscale image (legacy ``gray_scale_to_color_ramp``
        with float RGB output, ported verbatim).

        :param gray_scale: Normalized 2D grayscale image on the 0-1 range.
            NaN pixels map to the colormap's "bad" color (black). Not
            modified.
        :returns: Bands-first RGB float array on the 0-1 range.
        """
        color_map: Colormap = colormaps[self.name]
        if self.minimum_cut is not None or self.maximum_cut is not None:
            color_map = _truncate_colormap(color_map, *self._resolved_cuts())
        rgba = np.asarray(color_map(gray_scale))
        return rgba.transpose(2, 0, 1)[:3, ...]

    def _resolved_cuts(self) -> tuple[float, float]:
        """Return the cuts with ``None`` resolved to the full 0-1 range."""
        minimum = 0.0 if self.minimum_cut is None else self.minimum_cut
        maximum = 1.0 if self.maximum_cut is None else self.maximum_cut
        return minimum, maximum


def _truncate_colormap(color_map: Colormap, minimum: float, maximum: float) -> Colormap:
    """Legacy ``truncate_colormap``, ported verbatim (100-sample resampling)."""
    return LinearSegmentedColormap.from_list(
        f"trunc({color_map.name},{minimum:.2f},{maximum:.2f})",
        color_map(np.linspace(minimum, maximum, _TRUNCATION_SAMPLES)),
    )
