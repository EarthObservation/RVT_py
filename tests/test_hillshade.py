from __future__ import annotations

import numpy as np
import pytest

from rvt.core import Resolution
from rvt.visualization import Hillshade
from tests.conftest import read_golden_array, read_input_dem

RESOLUTION = Resolution(1.0, 1.0)


def test_defaults_match_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    hillshade = Hillshade().compute(dem, RESOLUTION)

    np.testing.assert_allclose(
        hillshade.astype(np.float32),
        read_golden_array("vis/hillshade/defaults__float.tif"),
        rtol=1e-7,
    )


def test_low_sun_with_exaggeration_matches_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    hillshade = Hillshade(
        sun_azimuth=90.0, sun_elevation=10.0, vertical_exaggeration=2.0
    ).compute(dem, RESOLUTION)

    np.testing.assert_allclose(
        hillshade.astype(np.float32),
        read_golden_array("vis/hillshade/az90_el10_ve2__float.tif"),
        rtol=1e-7,
    )


def test_nodata_dem_matches_golden() -> None:
    dem = read_input_dem("TM1_564_146_nodata_crop.tif")

    hillshade = Hillshade().compute(dem, RESOLUTION)

    np.testing.assert_allclose(
        hillshade.astype(np.float32),
        read_golden_array("vis_nodata/hillshade/defaults__float.tif"),
        rtol=1e-7,
    )


def test_rejects_invalid_sun_position() -> None:
    with pytest.raises(ValueError):
        Hillshade(sun_azimuth=361.0)
    with pytest.raises(ValueError):
        Hillshade(sun_elevation=-1.0)
