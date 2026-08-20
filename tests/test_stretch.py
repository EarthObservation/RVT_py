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
