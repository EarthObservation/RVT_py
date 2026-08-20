from __future__ import annotations

import numpy as np
import pytest

from rvt.blending import (
    BlendedSource,
    BlendMode,
    ColormappedSource,
    ColorRamp,
    ImageSource,
    MultiDirectionalHillshadeRgbSource,
    Normalization,
    NormalizationMode,
    OpennessDifferenceSource,
    SkyViewFactorSource,
    SlopeGradientSource,
    VisualizationSource,
    presets,
)
from rvt.core import Resolution
from rvt.visualization import (
    Hillshade,
    LocalDominance,
    MultiScaleTopographicPosition,
    NoiseRemoval,
    SimpleLocalReliefModel,
    SkyViewFactor,
    SlopeUnits,
)
from tests.conftest import read_golden_array, read_input_dem

RESOLUTION = Resolution(1.0, 1.0)

# Same rationale as tests/test_blending.py: goldens come from legacy float64
# SVF-family arrays, the new compute layer is float32-only, and the tight VALUE
# normalizations amplify the last-bit difference. Absolute error stays within
# one float32 ulp of the 0-1 render scale.
FLOAT_ABSOLUTE_TOLERANCE = 1e-6


def test_crim_preset_matches_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    rendered = presets.crim().render(dem, RESOLUTION)

    np.testing.assert_allclose(
        rendered.astype(np.float32),
        read_golden_array("blend/crim.tif"),
        rtol=1e-7,
        atol=FLOAT_ABSOLUTE_TOLERANCE,
    )


def test_e3mstp_preset_matches_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    rendered = presets.e3mstp().render(dem, RESOLUTION)

    np.testing.assert_allclose(
        rendered.astype(np.float32),
        read_golden_array("blend/e3mstp.tif"),
        rtol=1e-7,
        atol=FLOAT_ABSOLUTE_TOLERANCE,
    )


def test_e4mstp_preset_matches_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    rendered = presets.e4mstp(resolution=RESOLUTION).render(dem, RESOLUTION)

    np.testing.assert_allclose(
        rendered.astype(np.float32),
        read_golden_array("blend/e4mstp.tif"),
        rtol=1e-7,
        atol=FLOAT_ABSOLUTE_TOLERANCE,
    )


def test_crim_preset_structure() -> None:
    combination = presets.crim()

    assert combination.name == "Color relief image map (CRIM)"
    assert [layer.blend_mode for layer in combination.layers] == [
        BlendMode.OVERLAY,
        BlendMode.LUMINOSITY,
        BlendMode.NORMAL,
    ]
    assert [layer.opacity for layer in combination.layers] == [50.0, 50.0, 100.0]

    for openness_layer in combination.layers[:2]:
        assert isinstance(openness_layer.source, OpennessDifferenceSource)
        assert openness_layer.normalization == Normalization(
            NormalizationMode.VALUE, -28.0, 28.0
        )

    base = combination.layers[-1].source
    assert isinstance(base, ColormappedSource)
    assert base.color_ramp == ColorRamp("OrRd", 0.0, 1.0)
    assert base.normalization == Normalization(NormalizationMode.VALUE, 0.0, 0.8)
    assert isinstance(base.source, SlopeGradientSource)
    assert base.source.config.units == SlopeUnits.RADIAN
    # The legacy CRIM slope layer is tagged "slope gradient rad", which the
    # legacy normalizer does not invert (unlike plain "slope gradient").
    assert not base.source.inverts_normalization


def test_crim_preset_accepts_custom_color_ramp() -> None:
    color_ramp = ColorRamp("Reds_r", 0.0, 1.0)

    combination = presets.crim(color_ramp=color_ramp)

    base = combination.layers[-1].source
    assert isinstance(base, ColormappedSource)
    assert base.color_ramp == color_ramp


def test_e3mstp_preset_structure() -> None:
    combination = presets.e3mstp()

    assert combination.name == "e3MSTP"
    assert [layer.blend_mode for layer in combination.layers] == [
        BlendMode.SCREEN,
        BlendMode.SOFT_LIGHT,
        BlendMode.NORMAL,
    ]
    assert [layer.opacity for layer in combination.layers] == [25.0, 70.0, 100.0]

    slrm_layer, crim_layer, mstp_layer = combination.layers
    assert isinstance(slrm_layer.source, VisualizationSource)
    assert isinstance(slrm_layer.source.config, SimpleLocalReliefModel)
    assert slrm_layer.normalization == Normalization(NormalizationMode.VALUE, -0.5, 0.5)
    assert isinstance(crim_layer.source, BlendedSource)
    assert crim_layer.source.combination == presets.crim()
    assert isinstance(mstp_layer.source, VisualizationSource)
    # Legacy blends went through DefaultValues, whose mstp_lightness is 0.9.
    assert mstp_layer.source.config == MultiScaleTopographicPosition(lightness=0.9)


def test_e4mstp_preset_structure() -> None:
    combination = presets.e4mstp()

    assert combination.name == "e4MSTP"
    assert [layer.blend_mode for layer in combination.layers] == [
        BlendMode.OVERLAY,
        BlendMode.MULTIPLY,
        BlendMode.MULTIPLY,
        BlendMode.NORMAL,
    ]
    assert [layer.opacity for layer in combination.layers] == [
        90.0,
        25.0,
        100.0,
        100.0,
    ]

    mstp_layer, svf_layer, openness_layer, slope_layer = combination.layers
    assert isinstance(mstp_layer.source, VisualizationSource)
    assert mstp_layer.source.config == MultiScaleTopographicPosition(lightness=0.9)

    assert isinstance(svf_layer.source, BlendedSource)
    svf_top, svf_base = svf_layer.source.combination.layers
    assert isinstance(svf_top.source, SkyViewFactorSource)
    assert svf_top.normalization == Normalization(NormalizationMode.VALUE, 0.7, 1.0)
    assert isinstance(svf_base.source, SkyViewFactorSource)
    assert svf_base.source.config == SkyViewFactor(
        max_radius=10, noise_removal=NoiseRemoval.HIGH
    )
    assert svf_base.normalization == Normalization(NormalizationMode.VALUE, 0.9, 1.0)

    assert isinstance(openness_layer.source, BlendedSource)
    openness_top, openness_base = openness_layer.source.combination.layers
    assert isinstance(openness_top.source, OpennessDifferenceSource)
    assert openness_top.normalization == Normalization(
        NormalizationMode.VALUE, -15.0, 15.0
    )
    assert isinstance(openness_base.source, VisualizationSource)
    assert isinstance(openness_base.source.config, LocalDominance)
    assert openness_base.normalization == Normalization(
        NormalizationMode.VALUE, 0.5, 1.8
    )

    assert isinstance(slope_layer.source, ColormappedSource)
    assert slope_layer.source.color_ramp == ColorRamp("Reds_r", 0.0, 1.0)
    assert slope_layer.source.normalization == Normalization(
        NormalizationMode.VALUE, 0.0, 55.0
    )
    assert isinstance(slope_layer.source.source, SlopeGradientSource)
    # Legacy e4MSTP called rvt.vis.slope_aspect with its default radian units
    # while labelling the layer "Slope gradient" (which inverts).
    assert slope_layer.source.source.config.units == SlopeUnits.RADIAN
    assert slope_layer.source.source.inverts_normalization


def test_e4mstp_flat_terrain_radius_follows_resolution() -> None:
    combination = presets.e4mstp(resolution=Resolution(0.5, 0.5))

    svf_layer = combination.layers[1]
    assert isinstance(svf_layer.source, BlendedSource)
    svf_base = svf_layer.source.combination.layers[-1]
    assert isinstance(svf_base.source, SkyViewFactorSource)
    assert svf_base.source.config.max_radius == 20


def test_color_ramp_rejects_unknown_colormap() -> None:
    with pytest.raises(ValueError):
        ColorRamp("not_a_real_colormap")


def test_color_ramp_rejects_invalid_cuts() -> None:
    with pytest.raises(ValueError):
        ColorRamp("OrRd", -0.1, 1.0)
    with pytest.raises(ValueError):
        ColorRamp("OrRd", 0.0, 1.1)
    with pytest.raises(ValueError):
        ColorRamp("OrRd", 0.5, 0.5)
    with pytest.raises(ValueError):
        ColorRamp("OrRd", 0.8, 0.2)


def test_color_ramp_apply_matches_legacy_truncated_colormap() -> None:
    gray = np.array([[0.0, 0.5], [1.0, np.nan]], dtype=np.float32)

    rgb = ColorRamp("OrRd", 0.0, 1.0).apply(gray)

    assert rgb.shape == (3, 2, 2)
    assert rgb.min() >= 0 and rgb.max() <= 1
    # Matplotlib maps NaN to the transparent "bad" color; RGB becomes 0.
    np.testing.assert_array_equal(rgb[:, 1, 1], np.zeros(3))


def test_color_ramp_without_cuts_uses_untruncated_colormap() -> None:
    from matplotlib import colormaps

    gray = np.linspace(0, 1, 5, dtype=np.float32).reshape(1, 5)

    rgb = ColorRamp("OrRd").apply(gray)

    expected = np.asarray(colormaps["OrRd"](gray)).transpose(2, 0, 1)[:3]
    np.testing.assert_array_equal(rgb, expected)


def test_colormapped_source_rejects_rgb_source_image() -> None:
    source = ColormappedSource(
        source=ImageSource(np.zeros((3, 2, 2), dtype=np.float32)),
        normalization=Normalization(NormalizationMode.VALUE, 0.0, 1.0),
        color_ramp=ColorRamp("OrRd"),
    )

    with pytest.raises(ValueError):
        source.compute_image(np.zeros((2, 2), dtype=np.float32), RESOLUTION)


def test_colormapped_source_normalizes_and_inverts_before_coloring() -> None:
    image = np.array([[0.0, 25.0, 50.0]], dtype=np.float32)
    source = ColormappedSource(
        source=ImageSource(image, inverts_normalization=True),
        normalization=Normalization(NormalizationMode.VALUE, 0.0, 50.0),
        color_ramp=ColorRamp("OrRd", 0.0, 1.0),
    )

    rgb = source.compute_image(np.zeros((1, 3), dtype=np.float32), RESOLUTION)
    expected = ColorRamp("OrRd", 0.0, 1.0).apply(
        np.array([[1.0, 0.5, 0.0]], dtype=np.float32)
    )

    np.testing.assert_array_equal(rgb, expected)


def test_openness_difference_source_requires_positive_openness() -> None:
    with pytest.raises(ValueError):
        OpennessDifferenceSource(
            positive_config=SkyViewFactor(include_positive_openness=False)
        )


def test_blended_source_renders_nested_combination_as_float32() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")
    combination = presets.vat()

    image = BlendedSource(combination).compute_image(dem, RESOLUTION)

    assert image.dtype == np.float32
    np.testing.assert_array_equal(
        image, combination.render(dem, RESOLUTION).astype(np.float32)
    )


def test_multi_directional_hillshade_rgb_source_stacks_three_azimuths() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    image = MultiDirectionalHillshadeRgbSource().compute_image(dem, RESOLUTION)

    assert image.shape == (3, *dem.shape)
    for band, azimuth in enumerate((315.0, 22.5, 90.0)):
        np.testing.assert_array_equal(
            image[band],
            Hillshade(sun_azimuth=azimuth, sun_elevation=35.0).compute(dem, RESOLUTION),
        )


def test_multi_directional_hillshade_rgb_source_does_not_invert() -> None:
    assert not MultiDirectionalHillshadeRgbSource().inverts_normalization
