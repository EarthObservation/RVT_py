from __future__ import annotations

import numpy as np
import pytest

from rvt.core import Resolution
from rvt.stretch import Stretch, StretchMode, to_8bit
from rvt.visualization import Hillshade, SimpleLocalReliefModel, SlopeAspect, SlopeUnits
from tests.conftest import read_golden_array, read_input_dem

RESOLUTION = Resolution(1.0, 1.0)


@pytest.mark.parametrize(
    "input_name, golden_directory",
    [
        ("TM1_564_146_crop.tif", "vis"),
        ("TM1_564_146_nodata_crop.tif", "vis_nodata"),
    ],
)
def test_hillshade_8bit_matches_golden(input_name: str, golden_directory: str) -> None:
    dem = read_input_dem(input_name)

    hillshade = Hillshade().compute(dem, RESOLUTION).astype(np.float32)
    eight_bit = to_8bit(hillshade, Hillshade.EIGHT_BIT_STRETCH)

    np.testing.assert_array_equal(
        eight_bit,
        read_golden_array(f"{golden_directory}/hillshade/defaults__8bit.tif"),
    )


@pytest.mark.parametrize(
    "input_name, golden_directory",
    [
        ("TM1_564_146_crop.tif", "vis"),
        ("TM1_564_146_nodata_crop.tif", "vis_nodata"),
    ],
)
def test_slope_8bit_matches_golden(input_name: str, golden_directory: str) -> None:
    dem = read_input_dem(input_name)

    slope = (
        SlopeAspect(units=SlopeUnits.DEGREE)
        .compute(dem, RESOLUTION)
        .slope.astype(np.float32)
    )
    eight_bit = to_8bit(slope, SlopeAspect.EIGHT_BIT_STRETCH)

    np.testing.assert_array_equal(
        eight_bit,
        read_golden_array(f"{golden_directory}/slope_aspect/defaults_slope__8bit.tif"),
    )


@pytest.mark.parametrize(
    "input_name, golden_directory",
    [
        ("TM1_564_146_crop.tif", "vis"),
        ("TM1_564_146_nodata_crop.tif", "vis_nodata"),
    ],
)
def test_slrm_8bit_matches_golden(input_name: str, golden_directory: str) -> None:
    dem = read_input_dem(input_name)

    relief = SimpleLocalReliefModel().compute(dem, RESOLUTION).astype(np.float32)
    eight_bit = to_8bit(relief, SimpleLocalReliefModel.EIGHT_BIT_STRETCH)

    np.testing.assert_array_equal(
        eight_bit,
        read_golden_array(f"{golden_directory}/slrm/defaults__8bit.tif"),
    )


def test_to_8bit_does_not_mutate_input() -> None:
    array = np.array([[0.2, np.nan], [0.8, 1.5]], dtype=np.float32)
    original = array.copy()

    to_8bit(array, Stretch(StretchMode.VALUE, 0.0, 1.0))

    np.testing.assert_array_equal(array, original)


def test_nan_encodes_as_255() -> None:
    array = np.array([[0.0, np.nan]], dtype=np.float32)

    eight_bit = to_8bit(array, Stretch(StretchMode.VALUE, 0.0, 1.0))

    assert eight_bit[0, 1] == 255


def test_invalid_stretch_bounds_raise() -> None:
    with pytest.raises(ValueError):
        Stretch(StretchMode.VALUE, 1.0, 1.0)
    with pytest.raises(ValueError):
        Stretch(StretchMode.PERCENT, 60.0, 60.0)
    with pytest.raises(ValueError):
        Stretch(StretchMode.PERCENT, -1.0, 5.0)
    with pytest.raises(ValueError):
        Stretch(StretchMode.PERCENT, 101.0, 0.0)
    with pytest.raises(ValueError):
        Stretch(StretchMode.PERCENT, 0.0, 101.0)


def test_percent_stretch_cuts_at_percentiles() -> None:
    # 11 evenly spaced values: the 10th percentile is 1.0 and the 80th
    # (= 100 - 20 upper cutoff) is 8.0, so the normalized value is (v - 1) / 7
    # and the byte value is floor(255.9999 * that).
    array = np.arange(11, dtype=np.float32)

    eight_bit = to_8bit(array, Stretch(StretchMode.PERCENT, 10.0, 20.0))

    np.testing.assert_array_equal(
        eight_bit, [0, 0, 36, 73, 109, 146, 182, 219, 255, 255, 255]
    )


def test_percent_stretch_inverts() -> None:
    array = np.arange(11, dtype=np.float32)

    eight_bit = to_8bit(array, Stretch(StretchMode.PERCENT, 10.0, 20.0, invert=True))

    np.testing.assert_array_equal(
        eight_bit, [255, 255, 219, 182, 146, 109, 73, 36, 0, 0, 0]
    )


def test_percent_stretch_falls_back_to_min_max_when_percentiles_coincide() -> None:
    # A plateau puts both percentiles on the same value (5.0); the legacy code
    # then falls back to nanmin/nanmax, i.e. 0.0 and 10.0, so the plateau lands
    # at 0.5 -> floor(255.9999 * 0.5) = 127.
    array = np.array([0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 10], dtype=np.float32)

    eight_bit = to_8bit(array, Stretch(StretchMode.PERCENT, 10.0, 20.0))

    np.testing.assert_array_equal(eight_bit, [0] + [127] * 9 + [255])


def test_percent_stretch_on_constant_array_is_all_nodata() -> None:
    # Preserved legacy quirk: for a constant array even the nanmin/nanmax
    # fallback coincides, so the stretch divides by zero, every pixel becomes
    # NaN, and byte scaling encodes NaN as 255.
    array = np.full((2, 2), 5.0, dtype=np.float32)

    with np.errstate(invalid="ignore"):
        eight_bit = to_8bit(array, Stretch(StretchMode.PERCENT, 10.0, 20.0))

    np.testing.assert_array_equal(eight_bit, np.full((2, 2), 255, dtype=np.uint8))


def test_percent_stretch_ignores_nan_but_keeps_it_as_nodata() -> None:
    # The percentiles are computed over the 11 valid values only, so they match
    # test_percent_stretch_cuts_at_percentiles; the NaN comes out as 255.
    array = np.array([*range(11), np.nan], dtype=np.float32)

    eight_bit = to_8bit(array, Stretch(StretchMode.PERCENT, 10.0, 20.0))

    np.testing.assert_array_equal(
        eight_bit, [0, 0, 36, 73, 109, 146, 182, 219, 255, 255, 255, 255]
    )


def test_percent_stretch_does_not_mutate_input() -> None:
    array = np.array([0.0, 1.0, np.nan, 9.0, 10.0], dtype=np.float32)
    original = array.copy()

    to_8bit(array, Stretch(StretchMode.PERCENT, 10.0, 20.0))

    np.testing.assert_array_equal(array, original)
