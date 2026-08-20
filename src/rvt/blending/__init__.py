"""
Blending of visualizations: typed layers, blend modes, and combinations with
golden-master parity to the legacy blender.
"""

from rvt.blending import presets
from rvt.blending.blend_modes import BlendMode
from rvt.blending.colormap import ColorRamp
from rvt.blending.combination import (
    BlendCombination,
    BlendedSource,
    BlendLayer,
    ColormappedSource,
    ImageSource,
    LayerSource,
    MultiDirectionalHillshadeRgbSource,
    NegativeOpennessSource,
    OpennessDifferenceSource,
    PositiveOpennessSource,
    SkyViewFactorSource,
    SlopeGradientSource,
    VisualizationSource,
)
from rvt.blending.normalization import Normalization, NormalizationMode
from rvt.blending.serialization import combination_from_json, combination_to_json

__all__ = [
    "BlendCombination",
    "BlendLayer",
    "BlendMode",
    "BlendedSource",
    "ColorRamp",
    "ColormappedSource",
    "ImageSource",
    "LayerSource",
    "MultiDirectionalHillshadeRgbSource",
    "NegativeOpennessSource",
    "Normalization",
    "NormalizationMode",
    "OpennessDifferenceSource",
    "PositiveOpennessSource",
    "SkyViewFactorSource",
    "SlopeGradientSource",
    "VisualizationSource",
    "combination_from_json",
    "combination_to_json",
    "presets",
]
