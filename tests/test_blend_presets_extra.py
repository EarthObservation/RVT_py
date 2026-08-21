from __future__ import annotations

import numpy as np
import pytest
from numpy.typing import NDArray

from rvt.blending import (
    BlendCombination,
    BlendMode,
    NegativeOpennessSource,
    Normalization,
    NormalizationMode,
    PositiveOpennessSource,
    ShadowSource,
    SkyViewFactorSource,
    SlopeGradientSource,
    VisualizationSource,
    combination_from_json,
    combination_to_json,
    presets,
)
from rvt.core import Resolution
from rvt.visualization import (
    NegativeOpenness,
    ShadowHorizon,
    SkyIllumination,
    SkyModel,
    SkyViewFactor,
    SlopeAspect,
    SlopeUnits,
)
from tests.conftest import read_input_dem

RESOLUTION = Resolution(1.0, 1.0)
SMOKE_TEST_SIZE = 120


def small_dem() -> NDArray[np.float32]:
    """Read a small nodata-free slice of the golden input DEM."""
    dem = read_input_dem("TM1_564_146_crop.tif")
    return dem[:SMOKE_TEST_SIZE, :SMOKE_TEST_SIZE].copy()


def test_prismatic_openness_layer_stack_matches_legacy_definition() -> None:
    combination = presets.prismatic_openness()

    assert combination.name == "Prismatic openness"
    assert len(combination.layers) == 3

    positive_openness, negative_openness, slope = combination.layers

    assert positive_openness.source == PositiveOpennessSource(
        config=SkyViewFactor(include_positive_openness=True)
    )
    assert positive_openness.source.inverts_normalization is False
    assert positive_openness.normalization == Normalization(
        NormalizationMode.VALUE, 68.0, 93.0
    )
    assert positive_openness.blend_mode == BlendMode.NORMAL
    assert positive_openness.opacity == 50.0

    assert negative_openness.source == NegativeOpennessSource(config=NegativeOpenness())
    # Legacy inverts the normalization of "openness - negative".
    assert negative_openness.source.inverts_normalization is True
    assert negative_openness.normalization == Normalization(
        NormalizationMode.VALUE, 68.0, 95.0
    )
    assert negative_openness.blend_mode == BlendMode.NORMAL
    assert negative_openness.opacity == 50.0

    # The base layer: legacy computes "slope gradient" in degrees and inverts it.
    assert slope.source == SlopeGradientSource(
        config=SlopeAspect(units=SlopeUnits.DEGREE)
    )
    assert slope.source.inverts_normalization is True
    assert slope.normalization == Normalization(NormalizationMode.VALUE, 0.0, 50.0)
    assert slope.blend_mode == BlendMode.LUMINOSITY
    assert slope.opacity == 50.0


def test_city_layer_stack_matches_legacy_definition() -> None:
    combination = presets.city()

    assert combination.name == "City"
    assert len(combination.layers) == 3

    shadow, sky_view_factor, sky_illumination = combination.layers

    # Legacy computes the "Shadow" layer with the hillshade sun position.
    assert shadow.source == ShadowSource(
        config=ShadowHorizon(azimuth=315.0, elevation=35.0)
    )
    assert shadow.source.inverts_normalization is False
    assert shadow.normalization == Normalization(NormalizationMode.VALUE, 0.0, 1.0)
    assert shadow.blend_mode == BlendMode.NORMAL
    assert shadow.opacity == 70.0

    assert sky_view_factor.source == SkyViewFactorSource(config=SkyViewFactor())
    assert sky_view_factor.source.inverts_normalization is False
    assert sky_view_factor.normalization == Normalization(
        NormalizationMode.VALUE, 0.1, 1.0
    )
    assert sky_view_factor.blend_mode == BlendMode.NORMAL
    assert sky_view_factor.opacity == 10.0

    # The base layer: legacy computes "sky illumination" with the overcast sky
    # model, no cast shadow, 32 directions and a 100 px shadow distance.
    assert sky_illumination.source == VisualizationSource(
        SkyIllumination(
            sky_model=SkyModel.OVERCAST,
            compute_shadow=False,
            maximum_fine_radius=100,
            directions=32,
        )
    )
    assert sky_illumination.source.inverts_normalization is False
    assert sky_illumination.normalization == Normalization(
        NormalizationMode.PERCENT, 1.0, 0.0
    )
    assert sky_illumination.blend_mode == BlendMode.NORMAL
    assert sky_illumination.opacity == 100.0


@pytest.mark.parametrize(
    "combination",
    [presets.prismatic_openness(), presets.city()],
    ids=["prismatic_openness", "city"],
)
def test_preset_renders_finite_single_band_image(combination: BlendCombination) -> None:
    dem = small_dem()

    rendered = combination.render(dem, RESOLUTION)

    assert rendered.shape == (SMOKE_TEST_SIZE, SMOKE_TEST_SIZE)
    assert np.all(np.isfinite(rendered))
    assert rendered.min() >= 0.0
    assert rendered.max() <= 1.0


@pytest.mark.parametrize(
    "combination",
    [presets.prismatic_openness(), presets.city()],
    ids=["prismatic_openness", "city"],
)
def test_preset_round_trips_through_json(combination: BlendCombination) -> None:
    assert combination_from_json(combination_to_json(combination)) == combination
