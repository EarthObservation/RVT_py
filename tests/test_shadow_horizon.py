from __future__ import annotations

import numpy as np
import pytest

from rvt.core import Resolution
from rvt.visualization.shadow import ShadowHorizon
from tests.conftest import read_golden_array, read_input_dem

RESOLUTION = Resolution(1.0, 1.0)


def test_defaults_match_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    result = ShadowHorizon().compute(dem, RESOLUTION)

    np.testing.assert_allclose(
        result.horizon.astype(np.float32),
        read_golden_array("vis/shadow_horizon/defaults_horizon__float.tif"),
        rtol=1e-7,
    )
    np.testing.assert_allclose(
        result.shadow.astype(np.float32),
        read_golden_array("vis/shadow_horizon/defaults_shadow__float.tif"),
        rtol=1e-7,
    )


def test_east_low_sun_matches_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    result = ShadowHorizon(azimuth=90.0, elevation=10.0).compute(dem, RESOLUTION)

    np.testing.assert_allclose(
        result.horizon.astype(np.float32),
        read_golden_array("vis/shadow_horizon/az90_el10_horizon__float.tif"),
        rtol=1e-7,
    )
    np.testing.assert_allclose(
        result.shadow.astype(np.float32),
        read_golden_array("vis/shadow_horizon/az90_el10_shadow__float.tif"),
        rtol=1e-7,
    )


def test_nodata_dem_matches_golden() -> None:
    dem = read_input_dem("TM1_564_146_nodata_crop.tif")

    result = ShadowHorizon().compute(dem, RESOLUTION)

    golden_horizon = read_golden_array(
        "vis_nodata/shadow_horizon/defaults_horizon__float.tif"
    )
    assert np.isnan(golden_horizon).all()
    assert np.isnan(result.horizon).all()
    # The legacy binary shadow of an all-NaN horizon is all zeros, not NaN.
    np.testing.assert_allclose(
        result.shadow.astype(np.float32),
        read_golden_array("vis_nodata/shadow_horizon/defaults_shadow__float.tif"),
        rtol=1e-7,
    )


def test_does_not_mutate_input_dem() -> None:
    dem = read_input_dem("TM1_564_146_nodata_crop.tif")
    dem_before = dem.copy()

    ShadowHorizon().compute(dem, RESOLUTION)

    np.testing.assert_array_equal(dem, dem_before)


def test_rejects_invalid_parameters() -> None:
    with pytest.raises(ValueError):
        ShadowHorizon(azimuth=-1.0)
    with pytest.raises(ValueError):
        ShadowHorizon(elevation=91.0)
    with pytest.raises(ValueError):
        ShadowHorizon(vertical_exaggeration=-10001.0)


def test_rejects_invalid_dem_and_resolution() -> None:
    with pytest.raises(ValueError):
        ShadowHorizon().compute(np.zeros(10, dtype=np.float32), RESOLUTION)
    with pytest.raises(ValueError):
        ShadowHorizon().compute(
            np.zeros((10, 10), dtype=np.float32), Resolution(1.0, 0.0)
        )
