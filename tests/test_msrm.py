from __future__ import annotations

import numpy as np
import pytest

from rvt.core import Resolution
from rvt.visualization.msrm import MultiScaleReliefModel
from tests.conftest import read_golden_array, read_input_dem

RESOLUTION = Resolution(1.0, 1.0)


def test_defaults_match_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    multi_scale_relief_model = MultiScaleReliefModel(
        feature_minimum=0, feature_maximum=20, scaling_factor=2
    ).compute(dem, RESOLUTION)

    np.testing.assert_allclose(
        multi_scale_relief_model.astype(np.float32),
        read_golden_array("vis/msrm/defaults__float.tif"),
        rtol=1e-7,
    )


def test_custom_features_and_scaling_match_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    multi_scale_relief_model = MultiScaleReliefModel(
        feature_minimum=1, feature_maximum=50, scaling_factor=3
    ).compute(dem, RESOLUTION)

    np.testing.assert_allclose(
        multi_scale_relief_model.astype(np.float32),
        read_golden_array("vis/msrm/fmin1_fmax50_sf3__float.tif"),
        rtol=1e-7,
    )


def test_nodata_dem_matches_golden() -> None:
    dem = read_input_dem("TM1_564_146_nodata_crop.tif")

    multi_scale_relief_model = MultiScaleReliefModel(
        feature_minimum=0, feature_maximum=20, scaling_factor=2
    ).compute(dem, RESOLUTION)

    np.testing.assert_allclose(
        multi_scale_relief_model.astype(np.float32),
        read_golden_array("vis_nodata/msrm/defaults__float.tif"),
        rtol=1e-7,
    )


def test_does_not_mutate_input() -> None:
    dem = read_input_dem("TM1_564_146_nodata_crop.tif")
    original = dem.copy()

    MultiScaleReliefModel(
        feature_minimum=0, feature_maximum=20, scaling_factor=2
    ).compute(dem, RESOLUTION)

    np.testing.assert_array_equal(dem, original)


def test_rejects_out_of_range_exaggeration() -> None:
    with pytest.raises(ValueError):
        MultiScaleReliefModel(
            feature_minimum=0,
            feature_maximum=20,
            scaling_factor=2,
            vertical_exaggeration=10001.0,
        )
    with pytest.raises(ValueError):
        MultiScaleReliefModel(
            feature_minimum=0,
            feature_maximum=20,
            scaling_factor=2,
            vertical_exaggeration=-10001.0,
        )


def test_rejects_non_2d_dem() -> None:
    with pytest.raises(ValueError):
        MultiScaleReliefModel(
            feature_minimum=0, feature_maximum=20, scaling_factor=2
        ).compute(np.zeros((2, 3, 4), dtype=np.float32), RESOLUTION)


def test_rejects_non_positive_resolution() -> None:
    with pytest.raises(ValueError):
        MultiScaleReliefModel(
            feature_minimum=0, feature_maximum=20, scaling_factor=2
        ).compute(np.zeros((5, 5), dtype=np.float32), Resolution(1.0, -1.0))
