from __future__ import annotations

import numpy as np
import pytest

from rvt.core import Resolution
from rvt.visualization.sky_illumination import SkyIllumination, SkyModel
from tests.conftest import read_golden_array, read_input_dem

RESOLUTION = Resolution(1.0, 1.0)


def test_defaults_match_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    illumination = SkyIllumination().compute(dem, RESOLUTION)

    np.testing.assert_allclose(
        illumination.astype(np.float32),
        read_golden_array("vis/sky_illumination/defaults__float.tif"),
        rtol=1e-7,
    )


def test_uniform_with_shadow_matches_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    illumination = SkyIllumination(
        sky_model=SkyModel.UNIFORM, compute_shadow=True
    ).compute(dem, RESOLUTION)

    np.testing.assert_allclose(
        illumination.astype(np.float32),
        read_golden_array("vis/sky_illumination/uniform_shadow__float.tif"),
        rtol=1e-7,
    )


def test_nodata_dem_matches_golden() -> None:
    dem = read_input_dem("TM1_564_146_nodata_crop.tif")

    illumination = SkyIllumination().compute(dem, RESOLUTION)

    golden = read_golden_array("vis_nodata/sky_illumination/defaults__float.tif")
    assert np.isnan(golden).all()
    assert np.isnan(illumination).all()


def test_does_not_mutate_input_dem() -> None:
    dem = read_input_dem("TM1_564_146_nodata_crop.tif")
    dem_before = dem.copy()

    SkyIllumination().compute(dem, RESOLUTION)

    np.testing.assert_array_equal(dem, dem_before)


def test_rejects_invalid_parameters() -> None:
    with pytest.raises(ValueError):
        SkyIllumination(shadow_azimuth=361.0)
    with pytest.raises(ValueError):
        SkyIllumination(shadow_elevation=-1.0)
    with pytest.raises(ValueError):
        SkyIllumination(vertical_exaggeration=10001.0)
    with pytest.raises(ValueError):
        SkyIllumination(directions=0)
    with pytest.raises(ValueError):
        SkyIllumination(maximum_fine_radius=0)


def test_rejects_invalid_dem_and_resolution() -> None:
    with pytest.raises(ValueError):
        SkyIllumination().compute(np.zeros(10, dtype=np.float32), RESOLUTION)
    with pytest.raises(ValueError):
        SkyIllumination().compute(
            np.zeros((10, 10), dtype=np.float32), Resolution(-1.0, 1.0)
        )
