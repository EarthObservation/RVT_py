from __future__ import annotations

import numpy as np
import pytest

from rvt.core import Resolution
from rvt.visualization.local_dominance import LocalDominance
from tests.conftest import read_golden_array, read_input_dem

RESOLUTION = Resolution(1.0, 1.0)


def test_defaults_match_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    local_dominance = LocalDominance().compute(dem, RESOLUTION)

    np.testing.assert_allclose(
        local_dominance.astype(np.float32),
        read_golden_array("vis/local_dominance/defaults__float.tif"),
        rtol=1e-7,
    )


def test_custom_radii_and_observer_match_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    local_dominance = LocalDominance(
        minimum_radius=5,
        maximum_radius=50,
        radius_increment=2,
        angular_resolution=30,
        observer_height=3,
    ).compute(dem, RESOLUTION)

    np.testing.assert_allclose(
        local_dominance.astype(np.float32),
        read_golden_array(
            "vis/local_dominance/minrad5_maxrad50_inc2_ang30_oh3__float.tif"
        ),
        rtol=1e-7,
    )


def test_nodata_dem_matches_golden() -> None:
    dem = read_input_dem("TM1_564_146_nodata_crop.tif")

    local_dominance = LocalDominance().compute(dem, RESOLUTION)

    np.testing.assert_allclose(
        local_dominance.astype(np.float32),
        read_golden_array("vis_nodata/local_dominance/defaults__float.tif"),
        rtol=1e-7,
    )


def test_does_not_mutate_input() -> None:
    dem = read_input_dem("TM1_564_146_nodata_crop.tif")
    original = dem.copy()

    LocalDominance().compute(dem, RESOLUTION)

    np.testing.assert_array_equal(dem, original)


def test_rejects_out_of_range_exaggeration() -> None:
    with pytest.raises(ValueError):
        LocalDominance(vertical_exaggeration=10001.0)
    with pytest.raises(ValueError):
        LocalDominance(vertical_exaggeration=-10001.0)


def test_rejects_non_2d_dem() -> None:
    with pytest.raises(ValueError):
        LocalDominance().compute(np.zeros((2, 3, 4), dtype=np.float32), RESOLUTION)


def test_rejects_non_positive_resolution() -> None:
    with pytest.raises(ValueError):
        LocalDominance().compute(
            np.zeros((5, 5), dtype=np.float32), Resolution(0.0, 1.0)
        )
