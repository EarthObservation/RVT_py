from __future__ import annotations

import numpy as np
import pytest

from rvt.core import Resolution
from rvt.visualization.mstp import MultiScaleTopographicPosition
from tests.conftest import read_golden_array, read_input_dem

RESOLUTION = Resolution(1.0, 1.0)


def test_defaults_match_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    mstp = MultiScaleTopographicPosition().compute(dem, RESOLUTION)

    np.testing.assert_allclose(
        mstp.astype(np.float32),
        read_golden_array("vis/mstp/defaults__float.tif"),
        rtol=1e-7,
    )


def test_variant_scales_and_lightness_match_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    mstp = MultiScaleTopographicPosition(
        local_scale=(1, 10, 2),
        meso_scale=(10, 100, 10),
        broad_scale=(100, 500, 100),
        lightness=0.9,
    ).compute(dem, RESOLUTION)

    np.testing.assert_allclose(
        mstp.astype(np.float32),
        read_golden_array("vis/mstp/scales_var_light09__float.tif"),
        rtol=1e-7,
    )


def test_nodata_dem_matches_golden() -> None:
    dem = read_input_dem("TM1_564_146_nodata_crop.tif")

    mstp = MultiScaleTopographicPosition().compute(dem, RESOLUTION)

    np.testing.assert_allclose(
        mstp.astype(np.float32),
        read_golden_array("vis_nodata/mstp/defaults__float.tif"),
        rtol=1e-7,
    )


def test_does_not_mutate_input() -> None:
    dem = read_input_dem("TM1_564_146_nodata_crop.tif")
    original = dem.copy()

    MultiScaleTopographicPosition().compute(dem, RESOLUTION)

    np.testing.assert_array_equal(dem, original)


def test_rejects_scale_minimum_above_maximum() -> None:
    with pytest.raises(ValueError):
        MultiScaleTopographicPosition(local_scale=(5, 1, 1))
    with pytest.raises(ValueError):
        MultiScaleTopographicPosition(meso_scale=(50, 5, 5))
    with pytest.raises(ValueError):
        MultiScaleTopographicPosition(broad_scale=(500, 50, 50))


def test_rejects_scale_step_wider_than_range() -> None:
    with pytest.raises(ValueError):
        MultiScaleTopographicPosition(local_scale=(1, 5, 5))
    with pytest.raises(ValueError):
        MultiScaleTopographicPosition(meso_scale=(5, 50, 46))
    with pytest.raises(ValueError):
        MultiScaleTopographicPosition(broad_scale=(50, 500, 451))


def test_rejects_out_of_range_exaggeration() -> None:
    with pytest.raises(ValueError):
        MultiScaleTopographicPosition(vertical_exaggeration=-10001.0)


def test_rejects_non_2d_dem() -> None:
    with pytest.raises(ValueError):
        MultiScaleTopographicPosition().compute(
            np.zeros((2, 3, 4), dtype=np.float32), RESOLUTION
        )
