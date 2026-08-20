from __future__ import annotations

import numpy as np
import pytest

from rvt.blending import (
    BlendCombination,
    BlendLayer,
    BlendMode,
    ImageSource,
    Normalization,
    NormalizationMode,
    PositiveOpennessSource,
    SkyViewFactorSource,
    SlopeGradientSource,
    VisualizationSource,
    presets,
)
from rvt.core import Resolution
from rvt.visualization import Hillshade, SkyViewFactor
from tests.conftest import read_golden_array, read_input_dem

RESOLUTION = Resolution(1.0, 1.0)

# The goldens were rendered from legacy float64 sky-view factor arrays; the new
# compute layer is float32-only (ticket 06), and the VALUE 0.7-1.0 normalization
# amplifies that last-bit difference near its lower bound. Absolute error stays
# within one float32 ulp of the 0-1 render scale.
FLOAT_ABSOLUTE_TOLERANCE = 1e-6


def _assert_8bit_close(
    computed: np.ndarray, golden: np.ndarray, max_mismatched_fraction: float = 1e-4
) -> None:
    """8-bit goldens match except off-by-one bins on isolated pixels (float32 SVF, see above)."""
    difference = computed.astype(np.int16) - golden.astype(np.int16)
    assert np.abs(difference).max() <= 1
    assert (difference != 0).mean() <= max_mismatched_fraction


def _custom_combination() -> BlendCombination:
    return BlendCombination(
        name="custom",
        layers=(
            BlendLayer(
                source=VisualizationSource(Hillshade()),
                normalization=Normalization(NormalizationMode.PERCENT, 2.0, 3.0),
                blend_mode=BlendMode.MULTIPLY,
                opacity=60.0,
            ),
            BlendLayer(
                source=SlopeGradientSource(),
                normalization=Normalization(NormalizationMode.VALUE, 0.0, 50.0),
                blend_mode=BlendMode.OVERLAY,
                opacity=50.0,
            ),
            BlendLayer(
                source=SkyViewFactorSource(),
                normalization=Normalization(NormalizationMode.VALUE, 0.7, 1.0),
                blend_mode=BlendMode.NORMAL,
                opacity=100.0,
            ),
        ),
    )


def test_vat_preset_matches_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    rendered = presets.vat().render(dem, RESOLUTION)

    np.testing.assert_allclose(
        rendered.astype(np.float32),
        read_golden_array("blend/vat.tif"),
        rtol=1e-7,
        atol=FLOAT_ABSOLUTE_TOLERANCE,
    )


def test_vat_preset_8bit_matches_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    rendered_8bit = presets.vat().render_8bit(dem, RESOLUTION)

    _assert_8bit_close(rendered_8bit, read_golden_array("blend/vat_8bit.tif"))


def test_vat_preset_nodata_matches_golden() -> None:
    dem = read_input_dem("TM1_564_146_nodata_crop.tif")

    rendered = presets.vat().render(dem, RESOLUTION)

    np.testing.assert_allclose(
        rendered.astype(np.float32),
        read_golden_array("blend_nodata/vat.tif"),
        rtol=1e-7,
        atol=FLOAT_ABSOLUTE_TOLERANCE,
    )


def test_vat_preset_nodata_8bit_matches_golden() -> None:
    dem = read_input_dem("TM1_564_146_nodata_crop.tif")

    rendered_8bit = presets.vat().render_8bit(dem, RESOLUTION)

    _assert_8bit_close(rendered_8bit, read_golden_array("blend_nodata/vat_8bit.tif"))


def test_custom_combination_matches_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    rendered = _custom_combination().render(dem, RESOLUTION)

    np.testing.assert_allclose(
        rendered.astype(np.float32),
        read_golden_array("blend/custom.tif"),
        rtol=1e-7,
        atol=FLOAT_ABSOLUTE_TOLERANCE,
    )


def test_custom_combination_8bit_matches_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    rendered_8bit = _custom_combination().render_8bit(dem, RESOLUTION)

    _assert_8bit_close(rendered_8bit, read_golden_array("blend/custom_8bit.tif"))


def test_vat_preset_structure() -> None:
    combination = presets.vat()

    assert combination.name == "VAT - Archaeological"
    assert [layer.blend_mode for layer in combination.layers] == [
        BlendMode.MULTIPLY,
        BlendMode.OVERLAY,
        BlendMode.LUMINOSITY,
        BlendMode.NORMAL,
    ]
    assert [layer.opacity for layer in combination.layers] == [25.0, 50.0, 50.0, 100.0]


def test_does_not_mutate_input_dem() -> None:
    dem = read_input_dem("TM1_564_146_nodata_crop.tif")
    dem_before = dem.copy()

    presets.vat().render(dem, RESOLUTION)

    assert np.array_equal(dem, dem_before, equal_nan=True)


def test_does_not_mutate_layer_input_images() -> None:
    rng = np.random.default_rng(42)
    slope_like = (rng.random((8, 8)) * 60).astype(np.float32)
    slope_like[0, 0] = np.nan
    base_like = rng.random((8, 8)).astype(np.float32)
    slope_before = slope_like.copy()
    base_before = base_like.copy()
    combination = BlendCombination(
        name="images",
        layers=(
            BlendLayer(
                source=ImageSource(slope_like, inverts_normalization=True),
                normalization=Normalization(NormalizationMode.VALUE, 0.0, 50.0),
                blend_mode=BlendMode.OVERLAY,
                opacity=50.0,
            ),
            BlendLayer(
                source=ImageSource(base_like),
                normalization=Normalization(NormalizationMode.VALUE, 0.0, 1.0),
            ),
        ),
    )

    combination.render(np.zeros((8, 8), dtype=np.float32), RESOLUTION)

    assert np.array_equal(slope_like, slope_before, equal_nan=True)
    assert np.array_equal(base_like, base_before, equal_nan=True)


def test_rejects_invalid_opacity() -> None:
    source = ImageSource(np.zeros((2, 2), dtype=np.float32))
    normalization = Normalization(NormalizationMode.VALUE, 0.0, 1.0)

    with pytest.raises(ValueError):
        BlendLayer(source=source, normalization=normalization, opacity=-1.0)
    with pytest.raises(ValueError):
        BlendLayer(source=source, normalization=normalization, opacity=100.5)


def test_rejects_invalid_normalization() -> None:
    with pytest.raises(ValueError):
        Normalization(NormalizationMode.VALUE, 1.0, 1.0)
    with pytest.raises(ValueError):
        Normalization(NormalizationMode.VALUE, 2.0, 1.0)
    with pytest.raises(ValueError):
        Normalization(NormalizationMode.PERCENT, -1.0, 3.0)
    with pytest.raises(ValueError):
        Normalization(NormalizationMode.PERCENT, 60.0, 50.0)


def test_rejects_empty_combination() -> None:
    with pytest.raises(ValueError):
        BlendCombination(name="empty", layers=())


def test_positive_openness_source_requires_openness_output() -> None:
    with pytest.raises(ValueError):
        PositiveOpennessSource(SkyViewFactor(include_positive_openness=False))


def test_blend_mode_normal() -> None:
    active = np.array([[0.2, 0.8]], dtype=np.float32)
    background = np.array([[0.5, 0.5]], dtype=np.float32)

    np.testing.assert_array_equal(BlendMode.NORMAL.blend(active, background), active)


def test_blend_mode_screen() -> None:
    active = np.array([[0.0, 0.5, 1.0]], dtype=np.float32)
    background = np.array([[0.4, 0.4, 0.4]], dtype=np.float32)

    np.testing.assert_allclose(
        BlendMode.SCREEN.blend(active, background),
        np.array([[0.4, 0.7, 1.0]], dtype=np.float32),
        rtol=1e-7,
    )


def test_blend_mode_multiply() -> None:
    active = np.array([[0.0, 0.5, 1.0]], dtype=np.float32)
    background = np.array([[0.4, 0.4, 0.4]], dtype=np.float32)

    np.testing.assert_allclose(
        BlendMode.MULTIPLY.blend(active, background),
        np.array([[0.0, 0.2, 0.4]], dtype=np.float32),
        rtol=1e-7,
    )


def test_blend_mode_overlay() -> None:
    # background <= 0.5: 2 * background * active; background > 0.5:
    # 1 - (1 - 2 * (background - 0.5)) * (1 - active).
    active = np.array([[0.5, 0.5]], dtype=np.float32)
    background = np.array([[0.25, 0.75]], dtype=np.float32)

    np.testing.assert_allclose(
        BlendMode.OVERLAY.blend(active, background),
        np.array([[0.25, 0.75]], dtype=np.float32),
        rtol=1e-6,
    )


def test_blend_mode_soft_light() -> None:
    # active < 0.5: 2 * background * active + background^2 * (1 - 2 * active);
    # active >= 0.5: 2 * background * (1 - active) + sqrt(background) * (2 * active - 1).
    active = np.array([[0.25, 0.75]], dtype=np.float32)
    background = np.array([[0.64, 0.64]], dtype=np.float32)

    np.testing.assert_allclose(
        BlendMode.SOFT_LIGHT.blend(active, background),
        np.array([[2 * 0.64 * 0.25 + 0.64**2 * 0.5, 2 * 0.64 * 0.25 + 0.8 * 0.5]]),
        rtol=1e-6,
    )


def test_blend_mode_luminosity_gray_background_returns_active_luminance() -> None:
    active = np.array([[0.2, 0.8]], dtype=np.float32)
    background = np.array([[0.5, 0.5]], dtype=np.float32)

    np.testing.assert_array_equal(
        BlendMode.LUMINOSITY.blend(active, background), active
    )


def test_blend_mode_luminosity_rgb_background_shifts_channels() -> None:
    active = np.full((1, 1), 0.5, dtype=np.float32)
    background = np.array(
        [[[0.2]], [[0.4]], [[0.6]]], dtype=np.float32
    )  # lum = 0.3*0.2 + 0.59*0.4 + 0.11*0.6 = 0.362

    blended = BlendMode.LUMINOSITY.blend(active, background)

    np.testing.assert_allclose(
        blended,
        np.array([[[0.338]], [[0.538]], [[0.738]]]),
        rtol=1e-6,
    )


def test_blend_modes_do_not_mutate_inputs() -> None:
    rng = np.random.default_rng(7)
    active = rng.random((6, 6)).astype(np.float32)
    background = rng.random((6, 6)).astype(np.float32)
    active_before = active.copy()
    background_before = background.copy()

    for blend_mode in BlendMode:
        blend_mode.blend(active, background)

    np.testing.assert_array_equal(active, active_before)
    np.testing.assert_array_equal(background, background_before)


def test_normalization_inverts_for_slope_and_negative_openness_sources() -> None:
    from rvt.blending import NegativeOpennessSource

    assert SlopeGradientSource().inverts_normalization
    assert NegativeOpennessSource().inverts_normalization
    assert not SkyViewFactorSource().inverts_normalization
    assert not PositiveOpennessSource().inverts_normalization
