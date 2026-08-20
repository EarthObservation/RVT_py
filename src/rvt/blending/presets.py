"""
Shipped blend combinations as inspectable factory functions.
"""

from __future__ import annotations

from rvt.blending.blend_modes import BlendMode
from rvt.blending.combination import (
    BlendCombination,
    BlendLayer,
    PositiveOpennessSource,
    SkyViewFactorSource,
    SlopeGradientSource,
    VisualizationSource,
)
from rvt.blending.normalization import Normalization, NormalizationMode
from rvt.visualization import Hillshade


def vat() -> BlendCombination:
    """
    The VAT (Visualization for Archaeological Topography) combination:
    sky-view factor multiplied over positive openness, slope gradient
    luminosity, and a hillshade base (legacy ``blender_VAT.json``).

    :returns: The combination, ready to render.
    """
    return BlendCombination(
        name="VAT - Archaeological",
        layers=(
            BlendLayer(
                source=SkyViewFactorSource(),
                normalization=Normalization(NormalizationMode.VALUE, 0.7, 1.0),
                blend_mode=BlendMode.MULTIPLY,
                opacity=25.0,
            ),
            BlendLayer(
                source=PositiveOpennessSource(),
                normalization=Normalization(NormalizationMode.VALUE, 68.0, 93.0),
                blend_mode=BlendMode.OVERLAY,
                opacity=50.0,
            ),
            BlendLayer(
                source=SlopeGradientSource(),
                normalization=Normalization(NormalizationMode.VALUE, 0.0, 50.0),
                blend_mode=BlendMode.LUMINOSITY,
                opacity=50.0,
            ),
            BlendLayer(
                source=VisualizationSource(Hillshade()),
                normalization=Normalization(NormalizationMode.VALUE, 0.0, 1.0),
                blend_mode=BlendMode.NORMAL,
                opacity=100.0,
            ),
        ),
    )
