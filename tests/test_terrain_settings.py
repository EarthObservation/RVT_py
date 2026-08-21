from __future__ import annotations

import numpy as np
import pytest
from numpy.typing import NDArray

from rvt.blending import (
    BlendCombination,
    BlendedSource,
    BlendLayer,
    BlendMode,
    ColormappedSource,
    ColorRamp,
    MultiDirectionalHillshadeRgbSource,
    NegativeOpennessSource,
    Normalization,
    NormalizationMode,
    PositiveOpennessSource,
    ShadowSource,
    SkyViewFactorSource,
    SlopeGradientSource,
    Terrain,
    TerrainSettings,
    VisualizationSource,
    combination_from_json,
    combination_to_json,
    presets,
)
from rvt.core import Resolution
from rvt.visualization import (
    Hillshade,
    LocalDominance,
    MultiHillshade,
    NoiseRemoval,
    SimpleLocalReliefModel,
    SkyIllumination,
    SkyViewFactor,
)
from tests.conftest import read_input_dem

RESOLUTION = Resolution(1.0, 1.0)
SMOKE_TEST_SIZE = 120


def small_dem() -> NDArray[np.float32]:
    """Read a small nodata-free slice of the golden input DEM."""
    dem = read_input_dem("TM1_564_146_crop.tif")
    return dem[:SMOKE_TEST_SIZE, :SMOKE_TEST_SIZE].copy()


def value_bounds(minimum: float, maximum: float) -> Normalization:
    return Normalization(NormalizationMode.VALUE, minimum, maximum)


# Shipped terrain values


def test_general_terrain_carries_legacy_values() -> None:
    settings = Terrain.GENERAL.settings()

    assert settings.name == "general"
    assert settings.hillshade_sun_elevation == 35.0
    assert settings.multi_hillshade_sun_elevation == 35.0
    assert settings.simple_local_relief_model_radius_cells == 20
    assert settings.sky_view_factor_max_radius == 10
    assert settings.sky_view_factor_noise_removal == NoiseRemoval.NONE
    assert settings.local_dominance_minimum_radius == 10
    assert settings.local_dominance_maximum_radius == 20
    assert settings.local_dominance_observer_height == 1.7

    assert settings.slope_gradient_normalization == value_bounds(0.0, 50.0)
    assert settings.hillshade_normalization == value_bounds(0.0, 1.0)
    assert settings.multi_hillshade_normalization == value_bounds(0.0, 1.0)
    assert settings.simple_local_relief_model_normalization == value_bounds(-1.0, 1.0)
    assert settings.sky_view_factor_normalization == value_bounds(0.7, 1.0)
    assert settings.anisotropic_sky_view_factor_normalization == value_bounds(0.65, 1.0)
    assert settings.positive_openness_normalization == value_bounds(68.0, 93.0)
    assert settings.negative_openness_normalization == value_bounds(60.0, 95.0)
    # Version 2 parsed these two stretches into the negative openness field and
    # never applied them; here each one has its own field.
    assert settings.sky_illumination_normalization == value_bounds(0.25, 0.3)
    assert settings.local_dominance_normalization == value_bounds(0.5, 1.8)


def test_flat_terrain_carries_legacy_values() -> None:
    settings = Terrain.FLAT.settings()

    assert settings.name == "flat"
    assert settings.hillshade_sun_elevation == 15.0
    assert settings.multi_hillshade_sun_elevation == 15.0
    assert settings.simple_local_relief_model_radius_cells == 10
    assert settings.sky_view_factor_max_radius == 20
    assert settings.sky_view_factor_noise_removal == NoiseRemoval.HIGH
    assert settings.local_dominance_minimum_radius == 10
    assert settings.local_dominance_maximum_radius == 20
    assert settings.local_dominance_observer_height == 1.7

    assert settings.slope_gradient_normalization == value_bounds(0.0, 15.0)
    assert settings.sky_view_factor_normalization == value_bounds(0.9, 1.0)
    assert settings.positive_openness_normalization == value_bounds(85.0, 93.0)
    assert settings.negative_openness_normalization == value_bounds(75.0, 95.0)
    assert settings.local_dominance_normalization == value_bounds(0.5, 3.0)

    # The flat terrain has no stretch for these visualizations.
    assert settings.hillshade_normalization is None
    assert settings.multi_hillshade_normalization is None
    assert settings.simple_local_relief_model_normalization is None
    assert settings.anisotropic_sky_view_factor_normalization is None
    assert settings.sky_illumination_normalization is None


def test_steep_terrain_carries_legacy_values() -> None:
    settings = Terrain.STEEP.settings()

    assert settings.name == "steep"
    assert settings.hillshade_sun_elevation == 55.0
    assert settings.multi_hillshade_sun_elevation == 55.0
    assert settings.simple_local_relief_model_radius_cells == 10
    assert settings.sky_view_factor_max_radius == 10
    assert settings.sky_view_factor_noise_removal == NoiseRemoval.NONE
    assert settings.local_dominance_minimum_radius == 10
    assert settings.local_dominance_maximum_radius == 10
    assert settings.local_dominance_observer_height == 16.0

    assert settings.slope_gradient_normalization == value_bounds(0.0, 60.0)
    assert settings.sky_view_factor_normalization == value_bounds(0.55, 1.0)
    assert settings.positive_openness_normalization == value_bounds(55.0, 95.0)
    assert settings.negative_openness_normalization == value_bounds(45.0, 95.0)
    assert settings.local_dominance_normalization == value_bounds(0.55, 0.95)


def test_terrain_enum_values() -> None:
    assert [terrain.value for terrain in Terrain] == ["general", "flat", "steep"]


def test_empty_name_is_rejected() -> None:
    with pytest.raises(ValueError, match="must not be empty"):
        TerrainSettings(name="")


def test_inverted_local_dominance_radii_are_rejected() -> None:
    with pytest.raises(ValueError, match="local_dominance_minimum_radius"):
        TerrainSettings(
            name="broken",
            local_dominance_minimum_radius=20,
            local_dominance_maximum_radius=10,
        )


# Applying a terrain


def test_flat_terrain_applied_to_vat_overrides_parameters_and_bounds() -> None:
    adjusted = Terrain.FLAT.settings().apply(presets.vat())

    assert adjusted.name == "VAT - Archaeological"
    sky_view_factor, positive_openness, slope, hillshade = adjusted.layers

    assert sky_view_factor.source == SkyViewFactorSource(
        config=SkyViewFactor(max_radius=20, noise_removal=NoiseRemoval.HIGH)
    )
    assert sky_view_factor.normalization == value_bounds(0.9, 1.0)
    assert sky_view_factor.blend_mode == BlendMode.MULTIPLY
    assert sky_view_factor.opacity == 25.0

    assert positive_openness.source == PositiveOpennessSource(
        config=SkyViewFactor(
            max_radius=20,
            noise_removal=NoiseRemoval.HIGH,
            include_positive_openness=True,
        )
    )
    assert positive_openness.normalization == value_bounds(85.0, 93.0)

    # Slope gradient has no parameter override, only bounds.
    assert slope.source == SlopeGradientSource()
    assert slope.normalization == value_bounds(0.0, 15.0)

    assert hillshade.source == VisualizationSource(Hillshade(sun_elevation=15.0))
    assert hillshade.normalization == value_bounds(0.0, 1.0)


def test_general_terrain_applied_to_vat_keeps_the_preset_values() -> None:
    assert Terrain.GENERAL.settings().apply(presets.vat()) == presets.vat()


def test_applying_a_terrain_does_not_mutate_the_input() -> None:
    combination = presets.vat()
    before = combination_to_json(combination)

    Terrain.STEEP.settings().apply(combination)

    assert combination_to_json(combination) == before
    assert combination == presets.vat()


def test_empty_terrain_settings_are_a_no_op() -> None:
    combination = presets.e4mstp()

    assert TerrainSettings(name="empty").apply(combination) == combination


def test_terrain_reaches_every_source_kind() -> None:
    combination = BlendCombination(
        name="everything",
        layers=(
            BlendLayer(
                source=MultiDirectionalHillshadeRgbSource(),
                normalization=value_bounds(0.0, 2.0),
            ),
            BlendLayer(
                source=VisualizationSource(SimpleLocalReliefModel()),
                normalization=value_bounds(-2.0, 2.0),
            ),
            BlendLayer(
                source=VisualizationSource(MultiHillshade()),
                normalization=value_bounds(0.0, 2.0),
            ),
            BlendLayer(
                source=VisualizationSource(LocalDominance()),
                normalization=value_bounds(0.0, 2.0),
            ),
            BlendLayer(
                source=VisualizationSource(SkyIllumination()),
                normalization=value_bounds(0.0, 2.0),
            ),
            BlendLayer(
                source=SkyViewFactorSource(SkyViewFactor(include_anisotropic=True)),
                normalization=value_bounds(0.0, 2.0),
            ),
            BlendLayer(
                source=NegativeOpennessSource(),
                normalization=value_bounds(0.0, 2.0),
            ),
            BlendLayer(
                source=ColormappedSource(
                    source=SlopeGradientSource(),
                    normalization=value_bounds(0.0, 0.8),
                    color_ramp=ColorRamp("OrRd", 0.0, 1.0),
                ),
                normalization=value_bounds(0.0, 1.0),
            ),
            BlendLayer(
                source=ShadowSource(),
                normalization=value_bounds(0.0, 1.0),
            ),
        ),
    )

    adjusted = Terrain.GENERAL.settings().apply(combination)
    (
        multi_hillshade_rgb,
        simple_local_relief_model,
        multi_hillshade,
        local_dominance,
        sky_illumination,
        anisotropic,
        negative_openness,
        colormapped,
        shadow,
    ) = adjusted.layers

    assert multi_hillshade_rgb.source == MultiDirectionalHillshadeRgbSource(
        sun_elevation=35.0
    )
    assert multi_hillshade_rgb.normalization == value_bounds(0.0, 1.0)

    assert simple_local_relief_model.source == VisualizationSource(
        SimpleLocalReliefModel(radius_cells=20)
    )
    assert simple_local_relief_model.normalization == value_bounds(-1.0, 1.0)

    assert multi_hillshade.source == VisualizationSource(
        MultiHillshade(sun_elevation=35.0)
    )
    assert multi_hillshade.normalization == value_bounds(0.0, 1.0)

    assert local_dominance.source == VisualizationSource(
        LocalDominance(minimum_radius=10, maximum_radius=20, observer_height=1.7)
    )
    assert local_dominance.normalization == value_bounds(0.5, 1.8)

    assert sky_illumination.normalization == value_bounds(0.25, 0.3)
    assert anisotropic.normalization == value_bounds(0.65, 1.0)
    assert negative_openness.normalization == value_bounds(60.0, 95.0)

    # A colormapped source keeps its own normalization, the layer bounds and the
    # shadow layer are not terrain-dependent at all.
    assert colormapped.source == ColormappedSource(
        source=SlopeGradientSource(),
        normalization=value_bounds(0.0, 0.8),
        color_ramp=ColorRamp("OrRd", 0.0, 1.0),
    )
    assert colormapped.normalization == value_bounds(0.0, 1.0)
    assert shadow == BlendLayer(
        source=ShadowSource(), normalization=value_bounds(0.0, 1.0)
    )


def test_terrain_recurses_into_nested_combinations() -> None:
    nested = BlendCombination(
        name="nested",
        layers=(
            BlendLayer(
                source=VisualizationSource(Hillshade()),
                normalization=value_bounds(0.0, 2.0),
            ),
        ),
    )
    combination = BlendCombination(
        name="outer",
        layers=(
            BlendLayer(
                source=BlendedSource(nested), normalization=value_bounds(0.0, 1.0)
            ),
        ),
    )

    adjusted = Terrain.STEEP.settings().apply(combination)
    blended_source = adjusted.layers[0].source
    assert isinstance(blended_source, BlendedSource)
    inner_layer = blended_source.combination.layers[0]

    assert inner_layer.source == VisualizationSource(Hillshade(sun_elevation=55.0))
    # The steep terrain prescribes no hillshade bounds, so they stay as they were.
    assert inner_layer.normalization == value_bounds(0.0, 2.0)
    # The outer layer's own bounds are untouched.
    assert adjusted.layers[0].normalization == value_bounds(0.0, 1.0)


# VAT combined


def test_vat_combined_has_a_two_stage_structure() -> None:
    combination = presets.vat_combined()

    assert combination.name == "VAT combined - Archaeological"
    assert len(combination.layers) == 2
    general, flat = combination.layers

    assert general.blend_mode == BlendMode.NORMAL
    assert general.opacity == 50.0
    assert general.normalization == value_bounds(0.0, 1.0)
    assert flat.blend_mode == BlendMode.NORMAL
    assert flat.opacity == 100.0
    assert flat.normalization == value_bounds(0.0, 1.0)

    assert isinstance(general.source, BlendedSource)
    assert isinstance(flat.source, BlendedSource)
    assert general.source.combination == Terrain.GENERAL.settings().apply(presets.vat())
    assert flat.source.combination == Terrain.FLAT.settings().apply(presets.vat())


def test_vat_combined_opacity_and_terrains_are_configurable() -> None:
    combination = presets.vat_combined(
        general_opacity=70.0,
        general_terrain=Terrain.STEEP,
        flat_terrain=Terrain.GENERAL,
    )
    steep, general = combination.layers

    assert steep.opacity == 70.0
    assert isinstance(steep.source, BlendedSource)
    assert steep.source.combination == Terrain.STEEP.settings().apply(presets.vat())
    assert isinstance(general.source, BlendedSource)
    assert general.source.combination == Terrain.GENERAL.settings().apply(presets.vat())


def test_vat_combined_renders() -> None:
    dem = small_dem()

    image = presets.vat_combined().render(dem, RESOLUTION)

    assert image.shape == (SMOKE_TEST_SIZE, SMOKE_TEST_SIZE)
    assert np.all(np.isfinite(image))
    assert 0.0 <= float(np.min(image)) and float(np.max(image)) <= 1.0


# JSON round-trips


def test_vat_combined_round_trips_through_json() -> None:
    combination = presets.vat_combined()

    assert combination_from_json(combination_to_json(combination)) == combination


def test_terrain_adjusted_combination_round_trips_through_json() -> None:
    combination = Terrain.STEEP.settings().apply(presets.e3mstp())

    assert combination_from_json(combination_to_json(combination)) == combination
