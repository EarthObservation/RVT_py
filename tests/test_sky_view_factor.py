from __future__ import annotations

import numpy as np
import pytest

from rvt.core import Resolution
from rvt.visualization.sky_view_factor import (
    NegativeOpenness,
    NoiseRemoval,
    SkyViewFactor,
)
from tests.conftest import read_golden_array, read_input_dem

RESOLUTION = Resolution(1.0, 1.0)


def test_defaults_match_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    result = SkyViewFactor(
        include_anisotropic=True, include_positive_openness=True
    ).compute(dem, RESOLUTION)

    assert result.anisotropic is not None
    assert result.positive_openness is not None
    np.testing.assert_allclose(
        result.sky_view_factor.astype(np.float32),
        read_golden_array("vis/sky_view_factor/defaults_svf__float.tif"),
        rtol=1e-7,
    )
    np.testing.assert_allclose(
        result.anisotropic.astype(np.float32),
        read_golden_array("vis/sky_view_factor/defaults_asvf__float.tif"),
        rtol=1e-7,
    )
    np.testing.assert_allclose(
        result.positive_openness.astype(np.float32),
        read_golden_array("vis/sky_view_factor/defaults_opns__float.tif"),
        rtol=1e-7,
    )


def test_custom_parameters_match_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    result = SkyViewFactor(
        directions=8,
        max_radius=20,
        noise_removal=NoiseRemoval.MEDIUM,
        include_anisotropic=True,
        include_positive_openness=True,
        anisotropy_direction=90.0,
        anisotropy_level=2,
    ).compute(dem, RESOLUTION)

    assert result.anisotropic is not None
    assert result.positive_openness is not None
    np.testing.assert_allclose(
        result.sky_view_factor.astype(np.float32),
        read_golden_array(
            "vis/sky_view_factor/ndir8_rmax20_noise2_dir90_lvl2_svf__float.tif"
        ),
        rtol=1e-7,
    )
    np.testing.assert_allclose(
        result.anisotropic.astype(np.float32),
        read_golden_array(
            "vis/sky_view_factor/ndir8_rmax20_noise2_dir90_lvl2_asvf__float.tif"
        ),
        rtol=1e-7,
    )
    np.testing.assert_allclose(
        result.positive_openness.astype(np.float32),
        read_golden_array(
            "vis/sky_view_factor/ndir8_rmax20_noise2_dir90_lvl2_opns__float.tif"
        ),
        rtol=1e-7,
    )


def test_nodata_dem_matches_golden() -> None:
    dem = read_input_dem("TM1_564_146_nodata_crop.tif")

    result = SkyViewFactor(
        include_anisotropic=True, include_positive_openness=True
    ).compute(dem, RESOLUTION)

    assert result.anisotropic is not None
    assert result.positive_openness is not None
    np.testing.assert_allclose(
        result.sky_view_factor.astype(np.float32),
        read_golden_array("vis_nodata/sky_view_factor/defaults_svf__float.tif"),
        rtol=1e-7,
    )
    np.testing.assert_allclose(
        result.anisotropic.astype(np.float32),
        read_golden_array("vis_nodata/sky_view_factor/defaults_asvf__float.tif"),
        rtol=1e-7,
    )
    np.testing.assert_allclose(
        result.positive_openness.astype(np.float32),
        read_golden_array("vis_nodata/sky_view_factor/defaults_opns__float.tif"),
        rtol=1e-7,
    )


def test_negative_openness_defaults_match_golden() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    negative_openness = NegativeOpenness().compute(dem, RESOLUTION)

    np.testing.assert_allclose(
        negative_openness.astype(np.float32),
        read_golden_array("vis/sky_view_factor/neg_opns_defaults__float.tif"),
        rtol=1e-7,
    )


def test_negative_openness_nodata_dem_matches_golden() -> None:
    dem = read_input_dem("TM1_564_146_nodata_crop.tif")

    negative_openness = NegativeOpenness().compute(dem, RESOLUTION)

    np.testing.assert_allclose(
        negative_openness.astype(np.float32),
        read_golden_array("vis_nodata/sky_view_factor/neg_opns_defaults__float.tif"),
        rtol=1e-7,
    )


def test_optional_outputs_are_none_when_not_requested() -> None:
    dem = read_input_dem("TM1_564_146_crop.tif")

    result = SkyViewFactor().compute(dem, RESOLUTION)

    assert result.anisotropic is None
    assert result.positive_openness is None
    np.testing.assert_allclose(
        result.sky_view_factor.astype(np.float32),
        read_golden_array("vis/sky_view_factor/defaults_svf__float.tif"),
        rtol=1e-7,
    )


def test_does_not_mutate_input_dem() -> None:
    dem = read_input_dem("TM1_564_146_nodata_crop.tif")
    dem_before = dem.copy()

    SkyViewFactor(include_anisotropic=True, include_positive_openness=True).compute(
        dem, RESOLUTION
    )
    NegativeOpenness().compute(dem, RESOLUTION)

    assert np.array_equal(dem, dem_before, equal_nan=True)


def test_rejects_invalid_configuration() -> None:
    with pytest.raises(ValueError):
        SkyViewFactor(directions=0)
    with pytest.raises(ValueError):
        SkyViewFactor(max_radius=0)
    with pytest.raises(ValueError):
        SkyViewFactor(anisotropy_level=3)
    with pytest.raises(ValueError):
        SkyViewFactor(vertical_exaggeration=10001.0)
    with pytest.raises(ValueError):
        NegativeOpenness(directions=0)
    with pytest.raises(ValueError):
        NegativeOpenness(max_radius=0)
    with pytest.raises(ValueError):
        NegativeOpenness(vertical_exaggeration=-10001.0)


def test_rejects_invalid_dem_and_resolution() -> None:
    dem = np.zeros((5, 5), dtype=np.float32)

    with pytest.raises(ValueError):
        SkyViewFactor().compute(dem[0], RESOLUTION)
    with pytest.raises(ValueError):
        SkyViewFactor().compute(dem, Resolution(0.0, 1.0))
    with pytest.raises(ValueError):
        SkyViewFactor().compute(dem, Resolution(1.0, 2.0))
    with pytest.raises(ValueError):
        NegativeOpenness().compute(dem[0], RESOLUTION)
    with pytest.raises(ValueError):
        NegativeOpenness().compute(dem, Resolution(1.0, 2.0))
