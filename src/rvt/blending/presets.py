"""
Shipped blend combinations as inspectable factory functions.
"""

from __future__ import annotations

from rvt.blending.blend_modes import BlendMode
from rvt.blending.colormap import ColorRamp
from rvt.blending.combination import (
    BlendCombination,
    BlendedSource,
    BlendLayer,
    ColormappedSource,
    NegativeOpennessSource,
    OpennessDifferenceSource,
    PositiveOpennessSource,
    ShadowSource,
    SkyViewFactorSource,
    SlopeGradientSource,
    VisualizationSource,
)
from rvt.blending.normalization import Normalization, NormalizationMode
from rvt.core import Resolution
from rvt.visualization import (
    Hillshade,
    LocalDominance,
    MultiScaleTopographicPosition,
    NoiseRemoval,
    SimpleLocalReliefModel,
    SkyIllumination,
    SkyViewFactor,
    SlopeAspect,
    SlopeUnits,
)


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


def prismatic_openness() -> BlendCombination:
    """
    The prismatic openness combination: positive openness over negative
    openness over a slope gradient base, each at 50% opacity (legacy
    ``default_blender_combinations.json``, "Prismatic openness").

    :returns: The combination, ready to render.
    """
    return BlendCombination(
        name="Prismatic openness",
        layers=(
            BlendLayer(
                source=PositiveOpennessSource(),
                normalization=Normalization(NormalizationMode.VALUE, 68.0, 93.0),
                blend_mode=BlendMode.NORMAL,
                opacity=50.0,
            ),
            BlendLayer(
                source=NegativeOpennessSource(),
                normalization=Normalization(NormalizationMode.VALUE, 68.0, 95.0),
                blend_mode=BlendMode.NORMAL,
                opacity=50.0,
            ),
            BlendLayer(
                source=SlopeGradientSource(),
                normalization=Normalization(NormalizationMode.VALUE, 0.0, 50.0),
                blend_mode=BlendMode.LUMINOSITY,
                opacity=50.0,
            ),
        ),
    )


def city() -> BlendCombination:
    """
    The city combination: cast shadow over sky-view factor over a sky
    illumination base (legacy ``default_blender_combinations.json``, "City").

    :returns: The combination, ready to render.
    """
    return BlendCombination(
        name="City",
        layers=(
            BlendLayer(
                source=ShadowSource(),
                normalization=Normalization(NormalizationMode.VALUE, 0.0, 1.0),
                blend_mode=BlendMode.NORMAL,
                opacity=70.0,
            ),
            BlendLayer(
                source=SkyViewFactorSource(),
                normalization=Normalization(NormalizationMode.VALUE, 0.1, 1.0),
                blend_mode=BlendMode.NORMAL,
                opacity=10.0,
            ),
            BlendLayer(
                source=VisualizationSource(SkyIllumination()),
                normalization=Normalization(NormalizationMode.PERCENT, 1.0, 0.0),
                blend_mode=BlendMode.NORMAL,
                opacity=100.0,
            ),
        ),
    )


_IDENTITY_NORMALIZATION = Normalization(NormalizationMode.VALUE, 0.0, 1.0)

# The legacy blends went through DefaultValues, whose mstp_lightness (0.9)
# differs from the rvt.vis.mstp function default (1.2).
_MSTP_LEGACY_DEFAULTS = MultiScaleTopographicPosition(lightness=0.9)


def crim(*, color_ramp: ColorRamp = ColorRamp("OrRd", 0.0, 1.0)) -> BlendCombination:
    """
    The CRIM (color relief image map) combination: two positive-minus-negative
    openness layers (overlay and luminosity) over a colormapped slope base
    (legacy ``rvt.blend.color_relief_image_map``).

    :param color_ramp: The color ramp for the slope base layer.
    :returns: The combination, ready to render.
    """
    openness_normalization = Normalization(NormalizationMode.VALUE, -28.0, 28.0)
    return BlendCombination(
        name="Color relief image map (CRIM)",
        layers=(
            BlendLayer(
                source=OpennessDifferenceSource(),
                normalization=openness_normalization,
                blend_mode=BlendMode.OVERLAY,
                opacity=50.0,
            ),
            BlendLayer(
                source=OpennessDifferenceSource(),
                normalization=openness_normalization,
                blend_mode=BlendMode.LUMINOSITY,
                opacity=50.0,
            ),
            BlendLayer(
                source=ColormappedSource(
                    # Legacy tagged this layer "slope gradient rad", which its
                    # normalizer does not invert (unlike plain slope layers).
                    source=SlopeGradientSource(
                        config=SlopeAspect(units=SlopeUnits.RADIAN),
                        inverts_normalization=False,
                    ),
                    normalization=Normalization(NormalizationMode.VALUE, 0.0, 0.8),
                    color_ramp=color_ramp,
                ),
                normalization=_IDENTITY_NORMALIZATION,
                blend_mode=BlendMode.NORMAL,
                opacity=100.0,
            ),
        ),
    )


def e3mstp() -> BlendCombination:
    """
    The e3MSTP combination: simple local relief model screened over a CRIM
    soft-light layer and a multi-scale topographic position base (legacy
    ``rvt.blend.e3mstp``).

    :returns: The combination, ready to render.
    """
    return BlendCombination(
        name="e3MSTP",
        layers=(
            BlendLayer(
                source=VisualizationSource(SimpleLocalReliefModel()),
                normalization=Normalization(NormalizationMode.VALUE, -0.5, 0.5),
                blend_mode=BlendMode.SCREEN,
                opacity=25.0,
            ),
            BlendLayer(
                source=BlendedSource(crim()),
                normalization=_IDENTITY_NORMALIZATION,
                blend_mode=BlendMode.SOFT_LIGHT,
                opacity=70.0,
            ),
            BlendLayer(
                source=VisualizationSource(_MSTP_LEGACY_DEFAULTS),
                normalization=_IDENTITY_NORMALIZATION,
                blend_mode=BlendMode.NORMAL,
                opacity=100.0,
            ),
        ),
    )


def e4mstp(*, resolution: Resolution = Resolution(1.0, 1.0)) -> BlendCombination:
    """
    The e4MSTP combination: multi-scale topographic position over a combined
    sky-view factor sub-blend, a combined openness/local-dominance sub-blend,
    and a colormapped slope base (legacy ``rvt.blend.e4mstp``).

    :param resolution: Cell size of the DEM the combination will render.
        Only used to derive the hardcoded flat-terrain sky-view factor search
        radius (legacy ``int(10 / resolution)``, i.e. 10 map units); pass the
        same resolution given to ``render``.
    :returns: The combination, ready to render.
    """
    flat_terrain_svf = SkyViewFactor(
        max_radius=int(10 / resolution.x), noise_removal=NoiseRemoval.HIGH
    )
    combined_svf = BlendCombination(
        name="e4MSTP combined sky-view factor",
        layers=(
            BlendLayer(
                source=SkyViewFactorSource(),
                normalization=Normalization(NormalizationMode.VALUE, 0.7, 1.0),
                blend_mode=BlendMode.NORMAL,
                opacity=50.0,
            ),
            BlendLayer(
                source=SkyViewFactorSource(flat_terrain_svf),
                normalization=Normalization(NormalizationMode.VALUE, 0.9, 1.0),
                blend_mode=BlendMode.NORMAL,
                opacity=100.0,
            ),
        ),
    )
    combined_openness_local_dominance = BlendCombination(
        name="e4MSTP combined openness and local dominance",
        layers=(
            BlendLayer(
                source=OpennessDifferenceSource(),
                normalization=Normalization(NormalizationMode.VALUE, -15.0, 15.0),
                blend_mode=BlendMode.NORMAL,
                opacity=50.0,
            ),
            BlendLayer(
                source=VisualizationSource(LocalDominance()),
                normalization=Normalization(NormalizationMode.VALUE, 0.5, 1.8),
                blend_mode=BlendMode.NORMAL,
                opacity=100.0,
            ),
        ),
    )
    return BlendCombination(
        name="e4MSTP",
        layers=(
            BlendLayer(
                source=VisualizationSource(_MSTP_LEGACY_DEFAULTS),
                normalization=_IDENTITY_NORMALIZATION,
                blend_mode=BlendMode.OVERLAY,
                opacity=90.0,
            ),
            BlendLayer(
                source=BlendedSource(combined_svf),
                normalization=Normalization(NormalizationMode.VALUE, -0.5, 0.5),
                blend_mode=BlendMode.MULTIPLY,
                opacity=25.0,
            ),
            BlendLayer(
                source=BlendedSource(combined_openness_local_dominance),
                normalization=_IDENTITY_NORMALIZATION,
                blend_mode=BlendMode.MULTIPLY,
                opacity=100.0,
            ),
            BlendLayer(
                source=ColormappedSource(
                    # Legacy computed slope with rvt.vis.slope_aspect's default
                    # radian units, yet labelled the layer "Slope gradient",
                    # which its normalizer inverts; parity keeps both.
                    source=SlopeGradientSource(
                        config=SlopeAspect(units=SlopeUnits.RADIAN)
                    ),
                    normalization=Normalization(NormalizationMode.VALUE, 0.0, 55.0),
                    color_ramp=ColorRamp("Reds_r", 0.0, 1.0),
                ),
                normalization=_IDENTITY_NORMALIZATION,
                blend_mode=BlendMode.NORMAL,
                opacity=100.0,
            ),
        ),
    )
