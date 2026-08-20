from __future__ import annotations

import numpy as np
import pytest

from rvt.core import Resolution
from rvt.visualization import SimpleLocalReliefModel
from tests.conftest import read_golden_array, read_input_dem

RESOLUTION = Resolution(1.0, 1.0)


def test_defaults_match_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    simple_local_relief_model = SimpleLocalReliefModel().compute(dem, RESOLUTION)

    np.testing.assert_allclose(
        simple_local_relief_model.astype(np.float32),
        read_golden_array("vis/slrm/defaults__float.tif"),
        rtol=1e-7,
    )


def test_radius_50_matches_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    simple_local_relief_model = SimpleLocalReliefModel(radius_cells=50).compute(
        dem, RESOLUTION
    )

    np.testing.assert_allclose(
        simple_local_relief_model.astype(np.float32),
        read_golden_array("vis/slrm/rad50__float.tif"),
        rtol=1e-7,
    )


def test_nodata_dem_matches_golden() -> None:
    dem = read_input_dem("TM1_564_146_nodata_crop.tif")

    simple_local_relief_model = SimpleLocalReliefModel().compute(dem, RESOLUTION)

    np.testing.assert_allclose(
        simple_local_relief_model.astype(np.float32),
        read_golden_array("vis_nodata/slrm/defaults__float.tif"),
        rtol=1e-7,
    )


def test_does_not_mutate_input() -> None:
    dem = read_input_dem("TM1_564_146_nodata_crop.tif")
    original = dem.copy()

    SimpleLocalReliefModel().compute(dem, RESOLUTION)

    np.testing.assert_array_equal(dem, original)


def test_rejects_out_of_range_radius() -> None:
    with pytest.raises(ValueError):
        SimpleLocalReliefModel(radius_cells=9)
    with pytest.raises(ValueError):
        SimpleLocalReliefModel(radius_cells=51)


def test_rejects_out_of_range_exaggeration() -> None:
    with pytest.raises(ValueError):
        SimpleLocalReliefModel(vertical_exaggeration=10001.0)


def test_rejects_non_2d_dem() -> None:
    with pytest.raises(ValueError):
        SimpleLocalReliefModel().compute(
            np.zeros((2, 3, 4), dtype=np.float32), RESOLUTION
        )
