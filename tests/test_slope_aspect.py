from __future__ import annotations

import numpy as np
import pytest

from rvt.core import Resolution
from rvt.visualization import SlopeAspect, SlopeUnits
from tests.conftest import read_golden_array, read_input_dem

RESOLUTION = Resolution(1.0, 1.0)


def test_defaults_match_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    result = SlopeAspect().compute(dem, RESOLUTION)

    np.testing.assert_allclose(
        result.slope.astype(np.float32),
        read_golden_array("vis/slope_aspect/defaults_slope__float.tif"),
        rtol=1e-7,
    )
    np.testing.assert_allclose(
        result.aspect.astype(np.float32),
        read_golden_array("vis/slope_aspect/defaults_aspect__float.tif"),
        rtol=1e-7,
    )


def test_degree_with_exaggeration_matches_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    result = SlopeAspect(units=SlopeUnits.DEGREE, vertical_exaggeration=2.0).compute(
        dem, RESOLUTION
    )

    np.testing.assert_allclose(
        result.slope.astype(np.float32),
        read_golden_array("vis/slope_aspect/degree_ve2_slope__float.tif"),
        rtol=1e-7,
    )
    np.testing.assert_allclose(
        result.aspect.astype(np.float32),
        read_golden_array("vis/slope_aspect/degree_ve2_aspect__float.tif"),
        rtol=1e-7,
    )


def test_nodata_dem_matches_golden() -> None:
    dem = read_input_dem("TM1_564_146_nodata_crop.tif")

    result = SlopeAspect().compute(dem, RESOLUTION)

    np.testing.assert_allclose(
        result.slope.astype(np.float32),
        read_golden_array("vis_nodata/slope_aspect/defaults_slope__float.tif"),
        rtol=1e-7,
    )
    np.testing.assert_allclose(
        result.aspect.astype(np.float32),
        read_golden_array("vis_nodata/slope_aspect/defaults_aspect__float.tif"),
        rtol=1e-7,
    )


def test_does_not_mutate_input() -> None:
    dem = read_input_dem("TM1_564_146_nodata_crop.tif")
    original = dem.copy()

    SlopeAspect().compute(dem, RESOLUTION)

    np.testing.assert_array_equal(dem, original)


def test_rejects_out_of_range_exaggeration() -> None:
    with pytest.raises(ValueError):
        SlopeAspect(vertical_exaggeration=10001.0)


def test_rejects_non_2d_dem() -> None:
    with pytest.raises(ValueError):
        SlopeAspect().compute(np.zeros((2, 3, 4), dtype=np.float32), RESOLUTION)
