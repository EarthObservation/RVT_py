from __future__ import annotations

import numpy as np
import pytest
from numpy.typing import NDArray

from rvt.blending import Normalization, NormalizationMode
from rvt.blending.normalization import (
    _scale_strict_0_to_1,
    _scale_within_0_and_1,
    normalize_image,
    scale_0_to_1,
)


def test_apply_value_mode_clips_and_stretches() -> None:
    image = np.array([-1.0, 0.0, 5.0, 10.0, 11.0], dtype=np.float64)

    normalized = Normalization(NormalizationMode.VALUE, 0.0, 10.0).apply(image)

    np.testing.assert_allclose(normalized, [0.0, 0.0, 0.5, 1.0, 1.0])
    assert normalized.dtype == np.float32


def test_apply_percent_mode_cuts_at_percentiles() -> None:
    # 11 evenly spaced values: the 10th percentile is 1.0 and the 80th
    # (= 100 - 20 upper cutoff) is 8.0, so the stretch is (v - 1) / 7.
    image = np.arange(11, dtype=np.float64)

    normalized = Normalization(NormalizationMode.PERCENT, 10.0, 20.0).apply(image)

    np.testing.assert_allclose(
        normalized,
        [0.0, 0.0, 1 / 7, 2 / 7, 3 / 7, 4 / 7, 5 / 7, 6 / 7, 1.0, 1.0, 1.0],
        rtol=1e-6,
    )


def test_apply_percent_mode_falls_back_to_min_max_when_percentiles_coincide() -> None:
    # A plateau puts both percentiles on the same value (5.0); the legacy code
    # then falls back to nanmin/nanmax, i.e. 0.0 and 10.0.
    image = np.array([0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 10], dtype=np.float64)

    normalized = Normalization(NormalizationMode.PERCENT, 10.0, 20.0).apply(image)

    np.testing.assert_allclose(normalized, [0.0] + [0.5] * 9 + [1.0])


def test_apply_percent_mode_on_constant_image_yields_nan() -> None:
    # Preserved legacy quirk: for a constant image even the nanmin/nanmax
    # fallback coincides, so the stretch divides by zero and everything is NaN.
    image = np.full((2, 2), 5.0, dtype=np.float64)

    with np.errstate(invalid="ignore"):
        normalized = Normalization(NormalizationMode.PERCENT, 10.0, 20.0).apply(image)

    assert np.all(np.isnan(normalized))


def test_apply_propagates_nan() -> None:
    image = np.array([0.0, np.nan, 10.0], dtype=np.float64)

    value = Normalization(NormalizationMode.VALUE, 0.0, 10.0).apply(image)
    percent = Normalization(NormalizationMode.PERCENT, 0.0, 0.0).apply(image)

    np.testing.assert_array_equal(np.isnan(value), [False, True, False])
    np.testing.assert_array_equal(np.isnan(percent), [False, True, False])


def test_apply_does_not_mutate_input() -> None:
    image = np.array([-1.0, 0.0, 5.0, 11.0], dtype=np.float64)
    original = image.copy()

    Normalization(NormalizationMode.VALUE, 0.0, 10.0).apply(image)
    Normalization(NormalizationMode.PERCENT, 10.0, 20.0).apply(image)

    np.testing.assert_array_equal(image, original)


def test_rejects_invalid_bounds() -> None:
    with pytest.raises(ValueError):
        Normalization(NormalizationMode.VALUE, 1.0, 1.0)
    with pytest.raises(ValueError):
        Normalization(NormalizationMode.VALUE, 2.0, 1.0)
    with pytest.raises(ValueError):
        Normalization(NormalizationMode.PERCENT, -1.0, 3.0)
    with pytest.raises(ValueError):
        Normalization(NormalizationMode.PERCENT, 3.0, -1.0)
    with pytest.raises(ValueError):
        Normalization(NormalizationMode.PERCENT, 101.0, 0.0)
    with pytest.raises(ValueError):
        Normalization(NormalizationMode.PERCENT, 0.0, 101.0)
    with pytest.raises(ValueError):
        Normalization(NormalizationMode.PERCENT, 60.0, 50.0)


def test_normalize_image_inverts() -> None:
    image = np.array([0.0, 2.5, 10.0], dtype=np.float64)
    normalization = Normalization(NormalizationMode.VALUE, 0.0, 10.0)

    np.testing.assert_allclose(
        normalize_image(image, normalization, invert=True), [1.0, 0.75, 0.0]
    )


class _PassThroughNormalization(Normalization):
    """
    Test double that skips the 0-1 clip so ``normalize_image``'s rescue
    branches can be reached at all (see the note on the tests below).
    """

    def apply(self, image: NDArray[np.floating]) -> NDArray[np.float32]:
        return image.astype(np.float32)


# NOTE: ``Normalization.apply`` always ends with an unconditional clip to 0-1,
# so the two rescue branches in ``normalize_image`` cannot trigger for any real
# normalization -- they are a preserved legacy guard. The test double above
# feeds them the out-of-range input the guard was written for, and the assertions
# pin down the fact that the "rescue" leaves the values untouched: the legacy
# ``_scale_within_0_and_1`` bounds bug makes its rescaling the identity
# ``(image - 0) / (1 - 0)``.
def test_normalize_image_warns_and_rescues_when_result_exceeds_one() -> None:
    normalization = _PassThroughNormalization(NormalizationMode.VALUE, 0.0, 1.0)
    image = np.array([0.0, 0.5, 2.0], dtype=np.float32)

    with pytest.warns(UserWarning, match="max > 1"):
        normalized = normalize_image(image, normalization)

    np.testing.assert_allclose(normalized, [0.0, 0.5, 2.0])


def test_normalize_image_warns_and_rescues_when_result_below_zero() -> None:
    normalization = _PassThroughNormalization(NormalizationMode.VALUE, 0.0, 1.0)
    image = np.array([-0.5, 0.5], dtype=np.float32)

    with pytest.warns(UserWarning, match="min < 0"):
        normalized = normalize_image(image, normalization)

    np.testing.assert_allclose(normalized, [-0.5, 0.5])


def test_normalize_image_warns_twice_when_result_leaves_both_bounds() -> None:
    normalization = _PassThroughNormalization(NormalizationMode.VALUE, 0.0, 1.0)
    image = np.array([-0.5, 2.0], dtype=np.float32)

    with pytest.warns(UserWarning) as records:
        normalized = normalize_image(image, normalization)

    messages = [str(record.message) for record in records]
    assert any("max > 1" in message for message in messages)
    assert any("min < 0" in message for message in messages)
    np.testing.assert_allclose(normalized, [-0.5, 2.0])


def test_scale_0_to_1_passes_through_already_scaled_image() -> None:
    # First guard: max in (0.9, 1] and min exactly 0 counts as already scaled.
    image = np.array([0.0, 0.5, 0.95])

    assert scale_0_to_1(image) is image


def test_scale_0_to_1_uses_within_branch_for_wide_range() -> None:
    # Range > 0.3 routes to ``_scale_within_0_and_1``, whose own early return
    # passes an already-in-0-1 image straight through.
    image = np.array([0.0, 0.8])

    assert scale_0_to_1(image) is image


def test_scale_0_to_1_divides_8bit_looking_image_by_255() -> None:
    # Preserved legacy heuristic: max in (1, 255], span > 30 and min > 0 is
    # taken to mean "this is really an 8-bit image".
    image = np.array([10.0, 200.0])

    np.testing.assert_allclose(scale_0_to_1(image), [10 / 255, 200 / 255])


def test_scale_0_to_1_within_branch_replaces_nan_with_minimum() -> None:
    image = np.array([np.nan, -1.0, 1.0])

    # Preserved legacy bug: the bounds are the constants 0 and 1 rather than the
    # actual min/max, so the rescaling is the identity and values stay outside 0-1.
    np.testing.assert_allclose(scale_0_to_1(image), [-1.0, -1.0, 1.0])


def test_scale_0_to_1_within_branch_snaps_tiny_negatives_to_zero() -> None:
    image = np.array([-0.005, 2.0])

    np.testing.assert_allclose(scale_0_to_1(image), [0.0, 2.0])


def test_scale_0_to_1_uses_strict_branch_for_narrow_range() -> None:
    # Range <= 0.3 routes to ``_scale_strict_0_to_1``, which stretches to
    # exactly 0-1 using the real min/max.
    image = np.array([0.5, 0.55, 0.6])

    np.testing.assert_allclose(scale_0_to_1(image), [0.0, 0.5, 1.0])


def test_scale_strict_0_to_1_passes_through_exact_0_1_image() -> None:
    image = np.array([0.0, 1.0, 0.5])

    assert _scale_strict_0_to_1(image) is image


def test_scale_strict_0_to_1_maps_nan_to_zero() -> None:
    # NaN is set to 0 *before* the min/max are taken, so it becomes the new
    # minimum and ends up at 0.0 in the output.
    image = np.array([np.nan, 0.4, 0.5])

    np.testing.assert_allclose(_scale_strict_0_to_1(image), [0.0, 0.8, 1.0])


def test_scale_helpers_do_not_mutate_input() -> None:
    image = np.array([np.nan, 0.4, 0.5])
    original = image.copy()

    _scale_strict_0_to_1(image)
    _scale_within_0_and_1(np.array([np.nan, -1.0, 1.0]))
    scale_0_to_1(image)

    np.testing.assert_array_equal(image, original)
