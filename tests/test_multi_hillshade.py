from __future__ import annotations

import numpy as np
import pytest

from rvt.core import Resolution
from rvt.visualization import MultiHillshade
from tests.conftest import read_golden_array, read_input_dem

RESOLUTION = Resolution(1.0, 1.0)


def test_defaults_match_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    multi_hillshade = MultiHillshade().compute(dem, RESOLUTION)

    np.testing.assert_allclose(
        multi_hillshade.astype(np.float32),
        read_golden_array("vis/multi_hillshade/defaults__float.tif"),
        rtol=1e-7,
    )


def test_8_directions_low_sun_matches_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    multi_hillshade = MultiHillshade(directions=8, sun_elevation=10.0).compute(
        dem, RESOLUTION
    )

    np.testing.assert_allclose(
        multi_hillshade.astype(np.float32),
        read_golden_array("vis/multi_hillshade/dir8_el10__float.tif"),
        rtol=1e-7,
    )


def test_nodata_dem_matches_golden() -> None:
    dem = read_input_dem("TM1_564_146_nodata_crop.tif")

    multi_hillshade = MultiHillshade().compute(dem, RESOLUTION)

    np.testing.assert_allclose(
        multi_hillshade.astype(np.float32),
        read_golden_array("vis_nodata/multi_hillshade/defaults__float.tif"),
        rtol=1e-7,
    )


def test_does_not_mutate_input() -> None:
    dem = read_input_dem("TM1_564_146_nodata_crop.tif")
    original = dem.copy()

    MultiHillshade().compute(dem, RESOLUTION)

    np.testing.assert_array_equal(dem, original)


def test_rejects_non_positive_directions() -> None:
    with pytest.raises(ValueError):
        MultiHillshade(directions=0)


def test_rejects_out_of_range_sun_elevation() -> None:
    with pytest.raises(ValueError):
        MultiHillshade(sun_elevation=91.0)
    with pytest.raises(ValueError):
        MultiHillshade(sun_elevation=-1.0)


def test_rejects_out_of_range_exaggeration() -> None:
    with pytest.raises(ValueError):
        MultiHillshade(vertical_exaggeration=-10001.0)


def test_rejects_non_2d_dem() -> None:
    with pytest.raises(ValueError):
        MultiHillshade().compute(np.zeros((2, 3, 4), dtype=np.float32), RESOLUTION)
